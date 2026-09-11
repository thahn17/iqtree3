/***************************************************************************
 *   Reusable SPR (subtree pruning and regrafting) hill-climbing search.  *
 *                                                                        *
 *   This is the search machinery that used to live entirely inside       *
 *   search_experiments/spr_topology_test.cpp's anonymous namespace. It was lifted out  *
 *   verbatim so a SECOND caller could use it: IQTree::doSPRSearch(),     *
 *   which runs it as the per-iteration refinement stage of IQ-TREE's own *
 *   stochastic search, in place of doNNISearch -- see                    *
 *   IQTree::doTreeSearch and Params::refine_mode.                        *
 *                                                                        *
 *   Nothing here changed in the move: spr_topology_test --hillclimb runs *
 *   exactly the same code it always did, now by calling into this module *
 *   rather than by owning it. The one genuinely new thing is the         *
 *   SPRSearchOptions/SPRSearchState pair plus runSPRSteps/runSPRSweep at *
 *   the bottom, which package what used to be ~40 loose locals of        *
 *   runHillClimb so the step loop can be entered repeatedly (once per    *
 *   IQ-TREE perturbation) instead of exactly once per process.           *
 *                                                                        *
 *   Full flag reference: search_experiments/spr_topology_test_usage.txt.               *
 ***************************************************************************/

#ifndef SPRSEARCH_H
#define SPRSEARCH_H

#include "tree/phylotree.h"
#include "alignment/alignment.h"
#include <deque>
#include <string>
#include <unordered_map>
#include <vector>

namespace sprsearch {

using namespace std;

/**
    IQ-TREE's own short model label for this alignment's sequence type:
    JC or GTR+FO for DNA, LG or LG+FO for protein. `richerModel` picks the
    +FO ("gtr") variant. Only DNA and protein are supported.
 */
string modelNameFor(SeqType seqType, bool richerModel);

/**
    whether this sequence type is one the SPR search can handle at all
    (DNA or protein; binary/morphological/codon/PoMo are not).
 */
bool isSupportedSeqType(SeqType seqType);

/**
    build `t` as a fresh parse of `newickStr` against `aln`, with its own
    ModelFactory for `modelName`. Used for the throwaway scratch clones
    findopt scores against, and by runBranchLengthCompare's copies.
 */
void initClonedTree(PhyloTree &t, const string &newickStr, Alignment *aln, Params &params, string modelName);

/*==========================================================================
    Everything below is moved verbatim from spr_topology_test.cpp. The full
    explanatory comment for each function stays with its definition in
    search_experiments/sprsearch.cpp rather than being duplicated here.
 *========================================================================*/

void collectLeafNames(PhyloNode *node, PhyloNode *dad, vector<string> &names);

string describeEdge(PhyloNode *node, PhyloNode *dad);

string describeEdgeCompact(PhyloNode *node, PhyloNode *dad, size_t maxNames = 4);

/**
    a single legal SPR regraft target: the edge (node,dad), and how many
    hops away from the prune point it lies
 */
struct GraftCandidate {
    PhyloNode *node;
    PhyloNode *dad;
    int radius;
};

vector<GraftCandidate> findGraftPositions(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, int maxRadius);

/**
    every SPR-eligible edge in the tree, indexed by a stable slot number so
    choosePrune() can pick a uniformly random one in O(1) -- no traversal.
    "Eligible" excludes the one edge incident to the tree's arbitrary root
    leaf (tree.root), which isLegalSPR always rejects regardless of how it's
    addressed.

    Kept in sync by applySPRTracked/rollbackSPRTracked below instead of
    being rebuilt from scratch: an SPR move always destroys exactly 3 edges
    and creates exactly 3 new ones, reusing the same 3 freed slots (see the
    comment on applySPRTracked), so an update costs a small constant number
    of hash-map operations, never a tree traversal.
 */
struct EdgeRegistry {
    vector<pair<PhyloNode*, PhyloNode*> > slots;
    unordered_map<uint64_t, int> slotOfKey;

    // order-independent key for the edge {a,b}, from their (stable, unique)
    // node ids -- so looking up an edge doesn't care which side is passed
    // as a vs b
    static uint64_t key(PhyloNode *a, PhyloNode *b) {
        uint32_t lo = (uint32_t) min(a->id, b->id);
        uint32_t hi = (uint32_t) max(a->id, b->id);
        return (((uint64_t) hi) << 32) | lo;
    }

    void addEdge(PhyloNode *a, PhyloNode *b) {
        slotOfKey[key(a, b)] = (int) slots.size();
        slots.push_back(make_pair(a, b));
    }

    int slotOf(PhyloNode *a, PhyloNode *b) const {
        unordered_map<uint64_t, int>::const_iterator it = slotOfKey.find(key(a, b));
        ASSERT(it != slotOfKey.end());
        return it->second;
    }

    // overwrite 3 slots at once with 3 new edges. Must be batched like this
    // (erase all 3 stale keys first, only then insert all 3 new ones)
    // rather than done as 3 independent erase-then-insert calls: when a
    // regraft target lands next to one of the prune point's own siblings
    // (a common, legal case -- any radius-2 candidate), one of the 3 new
    // edges can have the exact same key as a DIFFERENT one of the 3 old
    // edges being replaced in the same batch. Inserting that new key
    // one-slot-at-a-time would silently overwrite the old mapping before
    // its own slot's turn to erase it comes up, so that erase call would
    // then wrongly delete the just-inserted (correct, still-needed) key
    // instead of the stale one -- permanently losing a key for an edge
    // that's still very much in the registry (this is exactly what caused
    // an intermittent "not found" assertion in slotOf() during testing).
    void replaceEdges3(int slotA, PhyloNode *a1, PhyloNode *a2,
                        int slotB, PhyloNode *b1, PhyloNode *b2,
                        int slotC, PhyloNode *c1, PhyloNode *c2) {
        slotOfKey.erase(key(slots[slotA].first, slots[slotA].second));
        slotOfKey.erase(key(slots[slotB].first, slots[slotB].second));
        slotOfKey.erase(key(slots[slotC].first, slots[slotC].second));

        slots[slotA] = make_pair(a1, a2);
        slots[slotB] = make_pair(b1, b2);
        slots[slotC] = make_pair(c1, c2);

        slotOfKey[key(a1, a2)] = slotA;
        slotOfKey[key(b1, b2)] = slotB;
        slotOfKey[key(c1, c2)] = slotC;
    }

    // same idea as replaceEdges3, batched the same way, but for the
    // 2-slot case used when one of dad1's siblings is the tree's root leaf
    // (see applySPRTracked) -- its edge to dad1, and the edge this move
    // creates in its place, are never tracked at all (buildEdgeRegistry
    // excludes root's one edge from the registry by design), so only 2 of
    // the usual 3 edges actually have real slots to update.
    void replaceEdges2(int slotA, PhyloNode *a1, PhyloNode *a2,
                        int slotB, PhyloNode *b1, PhyloNode *b2) {
        slotOfKey.erase(key(slots[slotA].first, slots[slotA].second));
        slotOfKey.erase(key(slots[slotB].first, slots[slotB].second));

        slots[slotA] = make_pair(a1, a2);
        slots[slotB] = make_pair(b1, b2);

        slotOfKey[key(a1, a2)] = slotA;
        slotOfKey[key(b1, b2)] = slotB;
    }
};

void collectEdgesExceptRoot(PhyloNode *node, PhyloNode *dad, EdgeRegistry &reg);

void buildEdgeRegistry(PhyloTree &tree, EdgeRegistry &reg);

bool subtreeContainsRoot(PhyloTree &tree, PhyloNode *node, PhyloNode *awayFrom);

bool isCherryAwayFrom(PhyloNode *node, PhyloNode *awayFrom);

/**
    how choosePrune picks which edge to prune.

    PRUNE_UNIFORM  every edge equally likely, regardless of length.
    PRUNE_LONG     weight proportional to the edge's own branch length --
                   "weightprune", or "weightprune long".
    PRUNE_SHORT    DEPRECATED -- accepted, warned about, and not worth
                   using. Weight proportional to the RECIPROCAL length,
                   softened by a floor at a tenth of the median so a
                   zero-length branch cannot swallow the whole
                   distribution ("weightprune short"). The bet was that
                   short internal branches, being the ones the data barely
                   resolves, are where the search is most likely to still
                   be wrong. Measured on sim.fa (100 taxa, 10k sites,
                   GTR+FO) it was the worst configuration tried: a mean
                   per-iteration gap of 18,264 logL to the optimum against
                   5,300 for an unweighted prune and 4,889 for PRUNE_LONG.
                   Relocating branches the data cannot place mostly
                   reshuffles pieces that were never the problem. Kept so
                   the comparison stays reproducible, not because it
                   works.
 */
enum PruneWeighting { PRUNE_UNIFORM, PRUNE_LONG, PRUNE_SHORT };

bool choosePrune(PhyloTree &tree, EdgeRegistry &reg, PhyloNode* &outNode, PhyloNode* &outDad,
        PruneWeighting weightPrune);

bool resolvePruneOrientationForSlot(PhyloTree &tree, PhyloNode *a, PhyloNode *b, bool preferA,
        PhyloNode* &outNode, PhyloNode* &outDad);

double graftDistanceWeight(int d);

int chooseGraftDistance(int r);

PhyloNode* collapsePruneDad(PhyloNode *from, PhyloNode *raw, PhyloNode *pruneDad, PhyloNode *pruneNode);

bool virtuallyAdjacent(PhyloNode *a, PhyloNode *b, PhyloNode *pruneDad);

bool chooseGraft(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, int r,
        PhyloNode* &outNode, PhyloNode* &outDad, int *outDistance = nullptr);

double distanceRadiusBudget(PhyloTree &tree, double radiusPercent);

bool chooseGraftByDistance(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, double radiusPercent,
        PhyloNode* &outNode, PhyloNode* &outDad, int *outHops = nullptr, double *outPercentUsed = nullptr);

/**
    bookkeeping stashed by applySPRTracked so rollbackSPRTracked can restore
    the registry (not just the tree) to its exact prior state.
 */
struct TrackedSPR {
    SPRMove move;
    SPRRollback rollback;
    PhyloNode *sibling1, *sibling2;
    int slotDad1Sib1, slotDad1Sib2, slotNode2Dad2; // -1 where not tracked (see below)
    bool sib1IsRoot, sib2IsRoot;
};

void applySPRTracked(PhyloTree &tree, EdgeRegistry &reg, const SPRMove &move, TrackedSPR &t);

void rollbackSPRTracked(PhyloTree &tree, EdgeRegistry &reg, const TrackedSPR &t);

void resetLikelihoodBuffers(PhyloTree &tree);

/**
    Local-recompute path for SPR trial scoring: instead of a full-tree
    delete/reinitialize/recompute for every candidate, install dedicated
    scratch buffers only on the directions whose represented subtree changes.
    Those are the rewired edges, prune_node->prune_dad, and one direction of
    every persistent edge on the removal-to-insertion path. The path entries
    are the essential correction to the six-direction attempt documented in
    PARTIAL_LIKELIHOOD_ATTEMPT.md: invalidating only the physically rewired
    edges left plausible but stale partials between the two locations.
 */
struct SPRLocalLhCache {
    vector<double*> lh;
    vector<UBYTE*> scaleNum;
    size_t lhCount, scaleCount;
};

void allocateSPRLocalLhCache(PhyloTree &tree, SPRLocalLhCache &cache);

void growSPRLocalLhCache(SPRLocalLhCache &cache, size_t count);

void freeSPRLocalLhCache(SPRLocalLhCache &cache);

/**
    per-call bookkeeping for beginLocalSPRInvalidation/end...Discard below:
    the affected PhyloNeighbor directions, and each
    one's saved (partial_lh, scale_num, lh_scale_factor, partial_lh_computed)
    from immediately before the swap -- exactly what it takes to put a
    direction back exactly as found, whether or not it happened to be null/
    uncomputed already (a real, common state under the default LM_PER_NODE
    memory mode: only "away from root" directions get a persistent buffer
    at startup, everything else starts null and is filled in on demand by
    PhyloTree::reorientPartialLh's buffer-stealing).
 */
struct SPRLocalInvalidationState {
    vector<PhyloNeighbor*> nei;
    vector<double*> savedLh;
    vector<UBYTE*> savedScaleNum;
    vector<double> savedScaleFactor;
    vector<int> savedComputed;
};

void beginLocalSPRInvalidation(SPRLocalLhCache &cache, SPRLocalInvalidationState &state,
        const vector<PhyloNeighbor*> &affected);

double computeLocalSPRLikelihood(PhyloTree &tree, PhyloNode *dad1, PhyloNode *dad2);

void endLocalSPRInvalidationDiscard(SPRLocalInvalidationState &state);

void findSPRSiblings(PhyloNode *node1, PhyloNode *dad1, PhyloNode *&sibling1, PhyloNode *&sibling2);

void clampBranchLengthForOptimization(PhyloNode *a, PhyloNode *b, double minLen);

void clampAllBranchLengthsForOptimization(PhyloTree &tree, double minLen);

void reoptimizeSPREdges(PhyloTree &tree, PhyloNode *dad1, PhyloNode *dad2, PhyloNode *node2, PhyloNode *node1,
        PhyloNode *sibling1, PhyloNode *sibling2);

double scoreTrialSPRMoveFullReset(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths);

bool findSPRNodePath(PhyloNode *node, PhyloNode *dad, PhyloNode *goal, PhyloNode *blocked,
        vector<PhyloNode*> &path);

void addUniqueSPRAffected(vector<PhyloNeighbor*> &affected, PhyloNeighbor *nei);

vector<PhyloNeighbor*> collectSPRAffectedBeforeApply(const SPRMove &move);

void addSPRChangedEdgeDirections(vector<PhyloNeighbor*> &affected,
        PhyloNode *dad1, PhyloNode *dad2, PhyloNode *node2,
        PhyloNode *sibling1, PhyloNode *sibling2);

double scoreTrialSPRMove(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths,
        SPRLocalLhCache *localCache = nullptr);

string buildRunId(int radius, int maxSteps, bool randomStart,
        bool useFastSelection, int numCandidates, bool reoptimizeBranchLengths,
        int fullReoptEveryNSteps, int fullReoptRounds, bool fullReoptInitialFit, bool useGtrModel,
        bool investigateFlag, int investigateRadius, bool alternateFlag, bool shrinkFlag,
        int shrinkStallThreshold, bool learnradiusFlag, int learnradiusN, bool sweepFlag, int sweepCount,
        int findoptEveryNSteps, bool iqtreeStart, int iqtreeStartPoolSize, bool weightpruneFlag);

string buildRecordTag(bool useFastSelection, bool useDistanceRadius, bool reoptimizeBranchLengths,
        int fullReoptEveryNSteps, bool investigateFlag, int investigateRadius, bool alternateFlag,
        bool shrinkFlag, bool learnradiusFlag, bool sweepFlag, int sweepCount, int findoptEveryNSteps,
        bool noTrueTree, bool weightpruneFlag, bool tunnelFlag = false, double tunnelTolerance = 0.0,
        bool slackFlag = false, double slackDelta = 0.0, bool slackAnneal = false);

string recordSpreadsheetPath(const string &modelName, const string &recordTag);

string topologySpreadsheetPath(const string &modelName, const string &recordTag);

void appendRecordRow(const string &modelName, const string &recordTag, const string &runId,
        long candidatesEvaluated, double timeElapsedSec, double logL, double trueTreeLogl,
        bool recordTopology, PhyloTree &treeForTopology);

string trajectoryTopologyPath(const string &runId);

void appendTrajectoryTopology(const string &runId, PhyloTree &tree);

void maybeRunPeriodicFullReopt(PhyloTree &tree, long successfulSteps, int fullReoptEveryNSteps, int fullReoptRounds,
        bool useGtrModel, bool quiet, bool recordProgress, bool recordTopology, const string &modelName,
        const string &recordTag, const string &runId, long candidatesEvaluated, double cpuClockStart,
        double trueTreeLogl, double &curScore);

void maybeRunFindopt(PhyloTree &tree, int step, int findoptEveryNSteps, bool quiet,
        bool recordProgress, bool recordTopology, const string &modelName, const string &recordTag,
        const string &runId, long candidatesEvaluated, double &cpuClockStart, double trueTreeLogl, double curScore,
        Alignment *aln, Params &params, const string &refitModelName = "");

void maybeShrinkRadius(bool improved, int shrinkStallThreshold, bool quiet,
        int &shrinkStallCount, int &shrinkCurrentRadius);

double sampleStandardNormal();

double sampleStandardGamma(double shape);

void recordLearnRadiusSample(deque<double> &window, int windowSize, double achievedRadius);

double sampleStartingRadius(double center, double maxPath);

double learnRadiusContinuous(const deque<double> &window, int windowSize, double center, double maxPath);

int edgeHopDistance(PhyloNode *seedA, PhyloNode *seedB, PhyloNode *targetA, PhyloNode *targetB,
        double *outLength = nullptr);

void maybeFinalizeLearnRadiusExcursion(PhyloTree &tree, bool useDistanceRadius, bool &excursionOpen,
        PhyloNode *anchorA, PhyloNode *anchorB, PhyloNode *excursionFinalNode, PhyloNode *currentDad,
        deque<double> &window, int windowSize);

bool computeSiblingCompatibilityScore(PhyloTree &tree, PhyloNode *p, PhyloNode *q, double curScore, double &outScore,
        PhyloNode *&outPruneNode, PhyloNode *&outPruneDad);

/*==========================================================================
    Re-entrant driver
    -----------------
    runHillClimb() used to hold all of the following as plain locals and
    run the step loop exactly once. IQ-TREE's stochastic search needs to
    re-enter that same loop once per perturbation ("kick"), on a tree
    object rebuilt from a Newick string each time, so the loop's
    configuration and its cross-iteration state are split into these two
    structs: SPRSearchOptions is fixed for a whole run, SPRSearchState
    carries what must survive from one pass to the next.
 *========================================================================*/

/**
    every --hillclimb flag that shapes an individual step. Deliberately
    NOT here: the flags that set up a whole standalone run rather than a
    refinement pass (random/iqtreestart/starttree/notree/model, and the
    choice of alignment) -- those belong to whoever built the tree, and
    the IQ-TREE-side parser rejects them outright. See parseRefineSpec
    at the bottom of this header.
 */
/**
    Stochastic acceptance rule: the probability of KEEPING a move, given the
    log-likelihood change it produces.

    Shape of the rule, for delta = newScore - curScore:

        delta >= 0            P = 1          (never refuse an improvement)
        delta <  0            P = exp(-(|delta| / T)^shape)

    so the whole (-inf, 0) side is a smooth, tunable curve and the [0, inf)
    side is flat at 1. T ("temperature") sets the scale of a tolerable loss
    and `shape` sets how sharply the curve falls off past it:

        shape = 1    Boltzmann / classic Metropolis -- P = exp(-|delta|/T)
        shape = 2    Gaussian-like: near-1 well inside T, then a fast cliff
        shape -> inf a hard threshold at |delta| = T, i.e. exactly what the
                     "slack" perturbation rule does, which makes slack the
                     limiting case of this same family rather than a
                     separate idea

    ANNEALING. With `anneal`, T is scaled by (1 - progress) as the run
    advances, floored at `tempFloor`. Progress is supplied by the caller
    (steps/budget for a continuous SPR stage; for the iterated search,
    IQTree::updateSearchProgress's clock, 0 on the first main-loop
    iteration and 1 on the last), so this struct stays agnostic about
    what is being counted. At T = 0 the rule degenerates to strict hill-climbing,
    which is why annealing to a floor of 0 is the sensible default: the run
    ends as a pure hill-climb no matter how exploratory it started.

    WHY THE DEFAULT T IS SMALL. The default (0.5) is chosen for final
    convergence rather than exploration: a move 0.5 logL worse is kept with
    P = 0.37, one 2 logL worse with P = 0.018, and one 5 logL worse with
    P = 4.5e-5. That admits the small backward steps that let a search
    round a barrier while making a genuinely destructive move essentially
    impossible. Raise T for a more exploratory run.
 */
struct AcceptDist {
    bool enabled;
    double temperature;            // T, the scale of a tolerable loss
    double shape;                  // exponent; 1 = Boltzmann
    bool anneal;                   // decay T toward tempFloor over the run
    double tempFloor;              // lowest T annealing may reach

    AcceptDist();

    /** T after annealing, for progress in [0, 1]. */
    double effectiveTemperature(double progress) const;

    /** P(keep) for this delta. Always 1 for delta >= 0. */
    double probability(double delta, double progress) const;

    /**
        Draw against probability(). Returns true for every improvement
        without consuming a random number, so an improve-only search under
        this rule stays bit-identical to one without it.
     */
    bool accept(double delta, double progress) const;
};

/**
 *  Per-kick state for --perturb-slack on IQ-TREE's NNI kick
 *  (IQTree::doRandomNNIs and the judgeNNIKickMove hook). active is false
 *  unless the flag was given, and then the kick loop is the stock one.
 */
struct KickSlack {
    enum Verdict { KEEP, REDRAW, GIVE_UP };
    static const int MAX_DRAWS = 8;     // redraws per slot before giving it up
    bool active = false;
    double delta = 0.0;                 // this kick's bound (after annealing)
    double curScore = 0.0;              // score after the last kept move
    long kept = 0, refused = 0;
    int draws = 0;
};

/**
 *  Per-call state for --accept-dist inside the NNI refiner
 *  (IQTree::optimizeNNI): the bounded stall counter, and the best tree the
 *  walk has passed through, restored on the way out.
 */
struct NNIAcceptState {
    static const int STALL_LIMIT = 5;   // consecutive non-improving steps
    int stall = 0;
    double bestScore = 0.0;
    std::string bestTree;
};

struct SPRSearchOptions {
    // --- candidate generation ---
    int radius;                    // <radius>: hop count, or (distradius) percent of tree length
    bool useFastSelection;         // "fast": random-walk draws instead of exhaustive enumeration
    int numCandidates;             // "fast N"
    bool useDistanceRadius;        // "distradius"
    PruneWeighting weightpruneFlag;  // "weightprune [long|short]"
    bool alternateFlag;            // "alternate"
    bool investigateFlag;          // "investigate"
    int investigateRadius;         // "investigate N"
    bool shrinkFlag;               // "shrink"
    int shrinkStallThreshold;      // "shrink N"
    bool learnradiusFlag;          // "learnradius"
    int learnradiusN;              // "learnradius N"

    // --- scoring / branch lengths ---
    bool reoptimizeBranchLengths;  // "reopt"
    int fullReoptEveryNSteps;      // "fullreopt M N" -- N
    int fullReoptRounds;           // "fullreopt M N" -- M
    bool useGtrModel;              // "gtr"

    // --- escaping a stall (continuous stage only) ---
    // "escape <span> <tries>": once `span` consecutive steps pass without
    // an improvement, save the current tree, apply a deliberately negative
    // random kick, and spend `tries` steps trying to beat the saved score.
    // Beat it and the excursion is kept; fail and the tree is rolled back
    // exactly. This is IQ-TREE's own perturb-then-refine idea moved INSIDE
    // one continuous SPR run -- with the difference that a failed
    // excursion reverts in place rather than restarting from the
    // candidate set.
    bool escapeFlag;
    int escapeSpan;                // non-improving steps that trigger a kick
    int escapeTries;               // steps allowed to beat the saved score

    /**
     *  "tunnel TOL [K]" -- do not roll back a MARGINALLY worse move
     *  immediately. When the step's best candidate loses no more than TOL
     *  log-likelihood, keep it applied and probe up to K random follow-up
     *  moves, ONE AT A TIME, testing whether the pair together beats the
     *  score the tree had before the first move. Each failed probe is
     *  rolled back on its own (a partial rollback, leaving the tolerated
     *  first move in place) so the next probe starts from the same
     *  tolerated position; only if every probe fails does the first move
     *  get rolled back too.
     *
     *  The point is that plain best-improvement SPR can only cross a
     *  barrier if some single move improves. A pair of moves whose first
     *  half is slightly downhill is invisible to it. TOL bounds how far
     *  downhill the search will step on spec, so the walk stays anchored
     *  near the incumbent instead of drifting: nothing is ever KEPT unless
     *  it beats the pre-move score outright.
     *
     *  Distinct from "escape", which kicks only after a stall and rolls
     *  the whole excursion back as a unit. This fires on an ordinary
     *  rejected step, costs at most K extra evaluations, and commits only
     *  a strict improvement.
     */
    bool tunnelFlag;
    double tunnelTolerance;        // max permissible logL decrease, > 0
    int tunnelTries;               // random follow-up probes, default 3

    /**
     *  "slack D [anneal]" -- threshold accepting. Lower the bar for what
     *  counts as an acceptable move: keep any candidate scoring above
     *  (curScore - D) instead of requiring a strict improvement. This is
     *  the general form of what "tunnel" does in a special case -- tunnel
     *  tolerates one downhill move only while a paired follow-up is being
     *  tested and never keeps it alone, whereas slack simply walks
     *  downhill whenever the step is small enough.
     *
     *  With "anneal", D decays linearly to 0 across the run's whole step
     *  budget (stepsPerPass), so the search starts exploratory and ends as
     *  a strict hill-climb. Without it, D is constant for the whole run.
     *
     *  This applies to the PERTURBATION only, never to refinement. A kick
     *  is supposed to move the tree somewhere worse; the point of slack is
     *  to bound HOW MUCH worse, so the perturbation lands in a nearby
     *  basin the refinement can still work with instead of anywhere at
     *  all. The refinement stage keeps its strict improve-only rule, so
     *  nothing downhill can survive into the reported tree.
     *
     *  Contrast --spr-perturb without it: moves are applied blind, with no
     *  score consulted, so a kick's damage is unbounded. Contrast "tunnel":
     *  that tolerates a downhill move inside the hill-climb, and only
     *  while a paired follow-up is being tested.
     */
    bool slackFlag;
    double slackDelta;             // permitted logL decrease, > 0
    bool slackAnneal;              // decay slackDelta to 0 over the budget

    /**
     *  --accept-dist: the stochastic acceptance rule applied to REFINEMENT
     *  moves (see AcceptDist). Disabled by default, in which case every
     *  refiner keeps its strict improve-only behaviour unchanged.
     */
    AcceptDist acceptDist;

    /**
     *  Annealing progress supplied by the CALLER, in [0, 1]. A continuous
     *  stage measures its own progress in steps against "steps N", but the
     *  iterated search has no step budget to measure against -- its clock
     *  is the iteration count, which only IQTree knows. It sets this each
     *  iteration so the cooling schedule advances in both modes; without
     *  it an in-loop SPR refinement would sit at the starting temperature
     *  forever and never converge.
     */
    double acceptProgressBase;

    // --- phases ---
    bool sweepFlag;                // "sweep"
    int sweepCount;                // "sweep N"
    bool findoptFlag;              // "findopt"
    int findoptEveryNSteps;        // "findopt N"

    // --- output ---
    bool quiet;                    // "quiet"
    bool recordProgress;           // "record"
    bool recordTopology;           // record's topology sidecar
    bool trajectoryFlag;           // "trajectory"
    bool noTrueTree;               // no ground-truth tree to measure a gap against

    // How many SPR steps one refinement pass gets. Only IQ-TREE reads it
    // ("steps N"): --hillclimb takes its step budget positionally and
    // passes it to runSPRSteps directly. 0 means "scale with the tree",
    // resolved by IQTree::doSPRSearch.
    int stepsPerPass;

    SPRSearchOptions();
};

/**
    what has to survive from one runSPRSteps() call to the next.

    Two lifetimes are mixed here on purpose. The numeric learning
    (learnRadiusWindow, the shrink counters, the running candidate/step
    totals, the record identity) is meant to accumulate across an entire
    run, IQ-TREE perturbations included. The PhyloNode* members are only
    ever valid within ONE pass, because IQ-TREE rebuilds its tree object
    from a Newick string on every kick, invalidating every node pointer;
    beginSPRSearchPass() clears exactly those and leaves the rest alone.
 */
struct SPRSearchState {
    // record/trajectory identity -- fixed once per run by initSPRSearchState
    string runId;
    string recordTag;
    string modelName;
    double cpuClockStart;
    double trueTreeLogl;            // NaN when there is no ground-truth tree
    double trueTreeLoglForFindopt;

    // running totals across every pass
    long candidatesEvaluated;

    /**
     *  "tunnel" accounting, reported at the end of a run: how many rejected
     *  steps were marginal enough to probe from, how many probe evaluations
     *  that cost, and how many of those excursions actually cleared the
     *  barrier. A large tunnelEntered with tunnelCommitted near zero means
     *  the flag is pure overhead on this data.
     */
    long tunnelEntered;
    long tunnelProbes;
    long tunnelCommitted;

    /**
     *  "slack" bookkeeping, accumulated across every kick in the run:
     *  slackAccepted counts perturbation moves actually applied,
     *  slackRejected counts proposals refused for costing more than the
     *  current threshold, and slackKicks counts kicks performed. The ratio
     *  of the first two says whether the threshold is doing any filtering
     *  at all -- if nothing is ever rejected, slack is just a slower way
     *  to write --spr-perturb.
     */
    long slackAccepted;
    long slackRejected;
    long slackKicks;
    double slackLastDelta;

    /**
     *  AcceptDist bookkeeping: how many moves were kept ONLY because the
     *  stochastic rule admitted them (a strict hill-climb would have
     *  refused), how many downhill moves were offered to it in total, and
     *  the last effective temperature used. acceptedDownhill near zero
     *  against a large offeredDownhill means the temperature is too cold
     *  to be doing anything.
     */
    long acceptedDownhill;
    long offeredDownhill;
    double lastTemperature;

    /**
     *  High-water mark for a stochastic run. Accepting downhill moves lets
     *  the CURRENT tree drift below the best one seen, so the best
     *  topology is remembered here and restored before the result is
     *  reported -- without this the run would report wherever the walk
     *  happened to stop.
     */
    double bestSeenScore;
    string bestSeenTree;
    bool bestSeenRestored;
    /** how many refinements ended in a best-seen restore, either refiner */
    long bestSeenRestores;
    long successfulSteps;
    int stepsRun;                   // total steps attempted, for findopt's cadence

    // "learnradius": the sliding window survives passes; the excursion
    // anchors do not (see the struct comment above)
    deque<double> learnRadiusWindow;
    double learnRadiusMaxPath;
    bool learnRadiusExcursionOpen;
    PhyloNode *learnRadiusAnchorA;
    PhyloNode *learnRadiusAnchorB;
    PhyloNode *learnRadiusFinalNode;

    // "investigate": pass-local, like the excursion anchors above
    bool investigateNext;
    PhyloNode *investigatePruneNode;
    PhyloNode *investigatePruneDad;

    // "shrink": survives passes
    int shrinkStallCount;
    int shrinkCurrentRadius;

    SPRSearchState();
};

/**
    fill in the run-scoped half of `st`: the record/trajectory identity
    strings, the CPU-clock baseline, the true-tree reference logL (pass
    NaN when there is none), and learnradius' own structural ceiling.
    Call once per run, before the first runSPRSteps.
 */
void initSPRSearchState(SPRSearchState &st, const SPRSearchOptions &opt, const EdgeRegistry &reg,
        const string &runId, const string &recordTag, const string &modelName,
        double cpuClockStart, double trueTreeLogl, double trueTreeLoglForFindopt);

/**
    reset the pass-local half of `st` -- every PhyloNode* it holds. Call
    at the start of each pass whose tree object is not the one the
    previous pass left behind (i.e. always, from IQ-TREE; a harmless
    no-op for --hillclimb's single pass).
 */
void beginSPRSearchPass(SPRSearchState &st);

/**
    the hill-climbing step loop itself, lifted verbatim out of
    runHillClimb. Runs up to `maxSteps` steps against `tree`/`reg`,
    starting from `curScore`, and returns the resulting logL. `reg` must
    already describe `tree` (see buildEdgeRegistry). `localCache` is the
    reusable path-partial scratch space -- pass nullptr under "reopt",
    which deliberately keeps the full-reset scoring path.

    `stepLabelPrefix` is prepended to each step's printed line so IQ-TREE
    can say which perturbation a step belongs to; empty for --hillclimb,
    which then prints exactly what it always did.
 */
double runSPRSteps(PhyloTree &tree, EdgeRegistry &reg, const SPRSearchOptions &opt, SPRSearchState &st,
        double curScore, int maxSteps, Alignment *aln, Params &params,
        SPRLocalLhCache *localCache, const string &stepLabelPrefix = "");

/**
    "sweep": the post-processing phase, also lifted verbatim. Returns the
    resulting logL.
 */
double runSPRSweep(PhyloTree &tree, EdgeRegistry &reg, const SPRSearchOptions &opt, SPRSearchState &st,
        double curScore);

/**
    the record CSV's unconditional final row -- see the comment at its
    definition for why one is always written regardless of whether the
    last event was accepted.
 */
void appendFinalRecordRow(const SPRSearchState &st, const SPRSearchOptions &opt, PhyloTree &tree, double curScore);

/**
    parse --spr-refine's / --nni-refine's flag string into `opt`.

    `spec` is passed through verbatim from the command line in exactly the
    vocabulary spr_topology_test --hillclimb takes its own trailing flags
    in, so the two tools share one flag language rather than two that can
    drift ("radius 10 fast quiet learnradius 30"). An empty spec is legal
    and means "this mode, all defaults".

    `sprMode` false is --nni-refine: IQ-TREE's own NNI search still does
    the refining, and only the refiner-independent flags (quiet, record,
    trajectory, findopt) are accepted -- anything that shapes an SPR step
    is an error naming --spr-refine.

    Flags that configure a whole standalone run rather than a refinement
    pass are rejected in BOTH modes, with an error pointing at the IQ-TREE
    option that already owns that job: random/iqtreestart/starttree (-t),
    model/gtr (-m), notree (implied -- there is never a ground-truth tree
    here). Returns false and fills `err` on any rejection; leaves `opt`
    unspecified in that case.
 */
bool parseRefineSpec(const string &spec, bool sprMode, SPRSearchOptions &opt, string &err);

/**
    parse --spr-perturb's flag string. A perturbation is a few random moves
    rather than a search, so only the flags describing a MOVE are accepted
    -- radius, distradius, weightprune [long|short]; everything else is an
    error. The kick's STRENGTH stays --perturb/initPS, exactly as for the
    NNI kick. Returns false and fills `err` on a rejection.
 */
bool parsePerturbSpec(const string &spec, SPRSearchOptions &opt, string &err);

/**
    apply `numMoves` random legal SPR moves to `tree` -- the SPR
    counterpart of IQTree::doRandomNNIs. Nothing is scored or rolled back;
    a kick is meant to make the tree worse. Returns how many moves actually
    landed (fewer than asked only if the topology runs out of legal prune
    or graft positions).
 */
int doRandomSPRs(PhyloTree &tree, const SPRSearchOptions &opt, int numMoves,
        SPRSearchState *st = nullptr, double annealScale = 1.0);

/**
    the record CSV's own tag for one refinement configuration, e.g.
    "_spr_fast_learnradius" or "_nni_findopt" -- buildRecordTag's output
    with the refiner's own name in front, so an SPR-refined run and an
    NNI-refined run of the same model never share a spreadsheet.
 */
string buildRefineRecordTag(const SPRSearchOptions &opt, bool sprMode);

} // namespace sprsearch

#endif // SPRSEARCH_H
