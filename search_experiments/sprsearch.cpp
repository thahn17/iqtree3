/***************************************************************************
 *   Reusable SPR hill-climbing search -- implementation.                 *
 *                                                                        *
 *   Moved verbatim out of search_experiments/spr_topology_test.cpp's anonymous         *
 *   namespace so both callers can share it (see search_experiments/sprsearch.h for     *
 *   why). The only genuinely new code in this file is the re-entrant     *
 *   driver at the bottom: runSPRSteps/runSPRSweep are runHillClimb's own *
 *   step loop and sweep phase, lifted unchanged, with the ~40 loose      *
 *   locals they used to read rebound onto SPRSearchOptions /             *
 *   SPRSearchState so the loop can be entered once per IQ-TREE           *
 *   perturbation instead of exactly once per process.                    *
 ***************************************************************************/

#include "sprsearch.h"
#include "model/modelfactory.h"
#include "utils/timeutil.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>

namespace sprsearch {

/**
    picks this tool's model name for an alignment's auto-detected sequence
    type (see Alignment::detectSequenceType, alignment/alignment.cpp --
    triggered by passing a null/empty sequence_type string into the
    Alignment constructor instead of a hardcoded one, which every alignment
    load in this file now does). Mirrors real IQ-TREE's own per-type
    "usual" default substitution model (getUsualModelSubst,
    main/phylotesting.cpp: "GTR" for DNA, "LG" for protein) rather than
    running a full ModelFinder search, since this is a lightweight test
    tool -- LG is IQ-TREE's own quick/default protein model, used the same
    way e.g. to build a fast guide tree without a full model-testing phase.
    `richerModel` mirrors this tool's existing "gtr" flag: false keeps the
    long-standing fixed-parameter default for that sequence type (JC for
    DNA, plain LG for protein); true asks for the ML-estimated-frequency
    variant ("+FO") instead -- GTR+FO already meant this for DNA; LG has no
    free rate-matrix parameters of its own the way GTR does, so "+FO" is
    the closest protein equivalent of "let the model fit more than just
    branch lengths". Only DNA and protein are supported; callers are
    expected to have already rejected any other detected SeqType.
 */
string modelNameFor(SeqType seqType, bool richerModel) {
    if (seqType == SEQ_PROTEIN)
        return richerModel ? "LG+FO" : "LG";
    return richerModel ? "GTR+FO" : "JC";
}

/**
    true iff `seqType` is one this tool's model handling
    (modelNameFor/ModelFactory setup) actually supports. Alignment content
    can auto-detect to other SeqType values (SEQ_BINARY, SEQ_MORPH,
    SEQ_CODON, ...) that this tool has never had any model logic for; call
    this right after loading an alignment and bail out with a clear error
    instead of silently mishandling one of those.
 */
bool isSupportedSeqType(SeqType seqType) {
    return seqType == SEQ_DNA || seqType == SEQ_PROTEIN;
}


/**
    build `t` as a fresh parse of `newickStr` against `aln`, with its own
    JC ModelFactory -- the common setup shared by trees B/C/D in
    runBranchLengthCompare (each an independent clone of tree A's own
    starting Newick text, see that function's own "starting point"
    comment for why a shared Newick parse, not a shared PhyloTree object,
    is what keeps them all bit-for-bit identical to start with).
 */
void initClonedTree(PhyloTree &t, const string &newickStr, Alignment *aln, Params &params, string modelName) {
    t.setParams(&params);
    t.read_TreeString(newickStr, false);
    t.setAlignment(aln);
    t.setNumThreads(1);
    t.setLikelihoodKernel(LK_SSE2);
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    t.setModelFactory(new ModelFactory(params, modelName, &t, modelsBlock));
    delete modelsBlock;
    t.setModel(t.getModelFactory()->model);
    t.setRate(t.getModelFactory()->site_rate);
    t.initializeAllPartialLh();
}

void collectLeafNames(PhyloNode *node, PhyloNode *dad, vector<string> &names) {
    if (node->isLeaf()) {
        names.push_back(node->name);
        return;
    }
    FOR_NEIGHBOR_IT(node, dad, it)
        collectLeafNames((PhyloNode*) (*it)->node, node, names);
}

string describeEdge(PhyloNode *node, PhyloNode *dad) {
    vector<string> names;
    collectLeafNames(node, dad, names);
    string result;
    for (size_t i = 0; i < names.size(); i++) {
        if (i)
            result += ",";
        result += names[i];
    }
    return result;
}

/**
    same leaf set as describeEdge, but capped at maxNames entries followed
    by "+K more" -- for compact one-line-per-step progress output where the
    full comma-joined list (potentially every leaf in a large clade) would
    be unreadable. Not meant to be fed back in as an edge spec argument;
    use describeEdge for that.
 */
string describeEdgeCompact(PhyloNode *node, PhyloNode *dad, size_t maxNames) {
    vector<string> names;
    collectLeafNames(node, dad, names);
    string result;
    size_t shown = min(names.size(), maxNames);
    for (size_t i = 0; i < shown; i++) {
        if (i)
            result += ",";
        result += names[i];
    }
    if (names.size() > shown)
        result += ",+" + to_string(names.size() - shown) + " more";
    return result;
}

/**
    enumerate every distinct, legal SPR regraft target within `maxRadius`
    edges of (pruneNode,pruneDad), where radius 1 is the nearest possible
    legal regraft target -- the edges incident to pruneDad's own two
    OTHER neighbors ("sibling1"/"sibling2", the ones that become directly
    connected once pruneDad is suppressed) -- equivalent to the smallest
    possible SPR move, an NNI; radius 2 is one hop further, and so on.
    pruneDad's own two immediate edges are never given a radius at all
    and never appear as a candidate, since isLegalSPR always rejects them
    (they're exactly the edges removed by suppressing pruneDad) -- this
    matches chooseGraft's own "d" convention (see its comment: "Step 1...
    pick uniformly among every real edge incident to sibling1 or
    sibling2"), so the SAME <radius> value passed to either "fast" or the
    exhaustive scan now means the same maximum real hop-distance in both.

    Implementation: a breadth-first search seeded at pruneDad's two OTHER
    neighbors directly (as if arriving FROM pruneDad, at distance 0),
    never crossing into pruneNode's subtree (so no candidate can ever lie
    inside the subtree being pruned), recording one candidate per edge the
    first time it is reached. A tree traversal can only ever visit each
    edge once, so duplicates are structurally impossible here -- but every
    edge is additionally checked against a canonical (min id, max id) key
    in `seenEdges` before being recorded, as an explicit, visible guarantee
    rather than an implicit one. Every raw candidate is then run through
    isLegalSPR (the same legality check applySPR itself asserts) anyway,
    as cheap defense-in-depth, even though seeding past pruneDad's own
    edges already means nothing structurally illegal should reach it.
 */
vector<GraftCandidate> findGraftPositions(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, int maxRadius) {
    vector<GraftCandidate> legalCandidates;
    set<pair<int,int> > seenEdges;

    struct QueueItem {
        PhyloNode *node;
        PhyloNode *cameFrom;
        int dist;
    };
    queue<QueueItem> q;
    FOR_NEIGHBOR_IT(pruneDad, pruneNode, it0)
        q.push({(PhyloNode*) (*it0)->node, pruneDad, 0});

    while (!q.empty()) {
        QueueItem cur = q.front();
        q.pop();
        if (cur.dist >= maxRadius)
            continue;

        FOR_NEIGHBOR_IT(cur.node, cur.cameFrom, it) {
            PhyloNode *next = (PhyloNode*) (*it)->node;
            int edgeRadius = cur.dist + 1;

            pair<int,int> key(min(cur.node->id, next->id), max(cur.node->id, next->id));
            if (seenEdges.insert(key).second) {
                SPRMove move;
                move.prune_node = pruneNode;
                move.prune_dad = pruneDad;
                move.regraft_node = next;
                move.regraft_dad = cur.node;
                move.radius = edgeRadius;
                move.screening_score = 0.0;
                move.exact_score = 0.0;
                move.candidate_id = (int) legalCandidates.size();
                move.generation = 0;

                if (tree.isLegalSPR(move))
                    legalCandidates.push_back({next, cur.node, edgeRadius});
            }

            q.push({next, cur.node, edgeRadius});
        }
    }

    return legalCandidates;
}

void collectEdgesExceptRoot(PhyloNode *node, PhyloNode *dad, EdgeRegistry &reg) {
    FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNode *child = (PhyloNode*) (*it)->node;
        reg.addEdge(node, child);
        collectEdgesExceptRoot(child, node, reg);
    }
}

/**
    build an EdgeRegistry over `tree`: one O(n) DFS, starting on the far
    side of tree.root's own single edge so that edge (the "1 edge at the
    top" -- tree.root is always this framework's arbitrary root leaf, never
    a real prune/graft target) is the one edge structurally never added,
    exactly matching what isLegalSPR always rejects.
 */
void buildEdgeRegistry(PhyloTree &tree, EdgeRegistry &reg) {
    reg.slots.clear();
    reg.slotOfKey.clear();
    PhyloNode *root = (PhyloNode*) tree.root;
    PhyloNode *belowRoot = (PhyloNode*) root->neighbors[0]->node;
    collectEdgesExceptRoot(belowRoot, root, reg);
}

/**
    true if tree.root lies within the subtree rooted at `node`, viewed away
    from `awayFrom` -- i.e. `node`'s side of the (node,awayFrom) edge.
    Used only for genuinely rooted trees, to orient a prune so the tree's
    real root always ends up on the "dad" (remaining tree) side, never the
    "node" (pruned branch) side; a plain recursive walk, not O(1), but only
    ever invoked when tree.rooted, which --hillclimb's own trees never are.
 */
bool subtreeContainsRoot(PhyloTree &tree, PhyloNode *node, PhyloNode *awayFrom) {
    if (node == (PhyloNode*) tree.root)
        return true;
    FOR_NEIGHBOR_IT(node, awayFrom, it)
        if (subtreeContainsRoot(tree, (PhyloNode*) (*it)->node, node))
            return true;
    return false;
}

/**
    true if `node` is a bare cherry when viewed away from `awayFrom` --
    its two OTHER neighbors (node is always degree 3 here, so exactly two
    remain once `awayFrom` is excluded) are both leaves. Used by
    choosePrune to keep a cherry from ever ending up as the "remaining
    tree" side of a prune: prune_dad's own two OTHER neighbors are exactly
    the "siblings" chooseGraft's and findGraftPositions' searches seed
    their walk from (see either one's own comment) -- if BOTH of those are
    leaves, neither has any further edge to walk onto at all, so the
    search comes back completely empty no matter how large the pruned
    subtree itself is. This is not about the pruned subtree being
    "illegal" to move; it is that there is nowhere left to move it TO,
    since everything reachable without re-entering the subtree being
    pruned is just those two leaves.
 */
bool isCherryAwayFrom(PhyloNode *node, PhyloNode *awayFrom) {
    FOR_NEIGHBOR_IT(node, awayFrom, it)
        if (!(*it)->node->isLeaf())
            return false;
    return true;
}

/**
    pick a uniformly random SPR-eligible edge from `reg` in O(1) (a single
    array-index pick, no tree traversal) and resolve it into a legal
    (outNode,outDad) prune pair. Requires init_random() to already have
    been called once.

    weightPrune (default false, the "weightprune" flag): instead of a
    uniform O(1) array-index pick, draws the slot with probability
    proportional to that edge's OWN branch length -- built as one
    cumulative-weight array over the whole registry (the same
    cumulative-sum technique chooseGraftDistance uses above, just applied
    to real branch lengths instead of a flat per-distance weight), then
    every retry attempt below draws from that same array via
    upper_bound (O(log E) per draw, O(E) once up front to build it,
    replacing choosePrune's normal O(1) per-pick cost). The intuition:
    long branches are the ones most likely to be hiding a misplaced
    subtree that ordinary short-hop NNI/SPR moves haven't already tried
    moving -- see e.g. Whelan & Money (2010)-style long-branch-attraction
    intuition for why a long branch is a natural place to concentrate
    search effort. Falls back to the ordinary uniform pick if every edge
    in the registry has zero total length (only possible in a degenerate
    or freshly-collapsed tree, where "proportional to length" has no
    well-defined meaning). EXPERIMENTAL.

    Every internal node is degree 3 and every leaf is degree 1 in the
    fully-bifurcating trees this tool works with, and applySPR/rollbackSPR
    never change any node's degree (see applySPRTracked below) -- so for a
    leaf-internal edge exactly one endpoint is a valid prune_dad, while for
    an internal-internal edge either endpoint would do.

    For an unrooted tree (the common case for --hillclimb's own BioNJ/
    Yule-Harding trees), which of the two degree-3 endpoints becomes dad is
    a coin flip, UNLESS exactly one endpoint is a bare cherry when viewed
    away from the other (isCherryAwayFrom) -- in that case the cherry is
    always forced to be the "node" (pruned/relocated) side, never "dad"
    (the remaining tree side chooseGraft/findGraftPositions search from),
    since leaving a cherry as the remaining side guarantees an empty
    search regardless of how large the OTHER side is (see
    isCherryAwayFrom's own comment). The edge itself is still picked
    uniformly at random from the whole registry -- only this ONE edge's
    two possible orientations are no longer equally likely once one side
    is a cherry. If NEITHER side is a cherry, or -- only possible in a
    tiny tree -- BOTH are, this falls back to the fair coin flip, since
    there's no orientation left that avoids the dead end anyway in that
    second case.

    For a genuinely rooted tree (tree.rooted), rootedness takes priority
    over the cherry check above: tree.root is a real topological root, and
    it must always end up on the "dad" (the tree that remains) side, never
    the "node" (pruned branch) side -- otherwise the move would effectively
    prune away the tree's own root. subtreeContainsRoot decides which side
    that is; if degree alone would force the opposite orientation (root's
    side isn't the degree-3 one), this edge has no valid direction at all
    and is skipped, not returned. (--hillclimb's own trees are never
    rooted, so in practice the cherry rule above is what actually applies.)

    @return false if the registry is empty (no edge to prune at all, e.g. a
    2-leaf tree), or if every attempted pick has no valid orientation
    (bounded retries, see below)
 */
bool choosePrune(PhyloTree &tree, EdgeRegistry &reg, PhyloNode* &outNode, PhyloNode* &outDad,
        PruneWeighting weightPrune) {
    if (reg.slots.empty())
        return false;

    // weightPrune: build the cumulative-weight array ONCE, outside the
    // retry loop below -- the registry itself never changes across
    // attempts within one choosePrune() call, so every attempt can safely
    // draw from the same distribution instead of rebuilding it each time.
    vector<double> cumWeight;
    double totalWeight = 0.0;
    if (weightPrune != PRUNE_UNIFORM) {
        vector<double> len(reg.slots.size());
        for (size_t i = 0; i < reg.slots.size(); i++) {
            pair<PhyloNode*, PhyloNode*> &e = reg.slots[i];
            len[i] = e.first->findNeighbor(e.second)->length;
        }

        // PRUNE_SHORT is the mirror image of PRUNE_LONG: weight by the
        // RECIPROCAL length rather than the length itself, so the short
        // edges -- the poorly resolved ones the data barely supports, and
        // so the ones most likely to be sitting in the wrong place -- are
        // the ones proportionally more likely to be pruned.
        //
        // A raw 1/len would be unbounded: a zero-length branch (common
        // after a topology change leaves a placeholder) would take
        // essentially the whole probability mass, collapsing the draw onto
        // one edge. Softening it with a floor tied to the MEDIAN length
        // keeps the bias scale-free -- it means the same thing on a tree
        // whose branches are ~0.001 as on one whose branches are ~1 -- and
        // caps how far any single edge can dominate.
        double floorLen = 0.0;
        if (weightPrune == PRUNE_SHORT) {
            vector<double> sorted = len;
            sort(sorted.begin(), sorted.end());
            const size_t n = sorted.size();
            double median = (n % 2 == 0) ? 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]) : sorted[n / 2];
            floorLen = 0.1 * median;
            if (floorLen <= 0.0)
                floorLen = 1e-9; // degenerate tree: every length zero or below
        }

        cumWeight.resize(reg.slots.size());
        for (size_t i = 0; i < reg.slots.size(); i++) {
            totalWeight += (weightPrune == PRUNE_SHORT) ? 1.0 / (len[i] + floorLen) : len[i];
            cumWeight[i] = totalWeight;
        }
        if (totalWeight <= 0.0)
            weightPrune = PRUNE_UNIFORM; // nothing to weight by -- fall back to uniform
    }

    // bounded retry: for a rooted tree, only edges adjacent to the root
    // can ever fail the orientation check below, a small fraction of the
    // registry, so this succeeds within the first few attempts in
    // practice; the bound just guarantees termination
    for (int attempt = 0; attempt < (int) reg.slots.size(); attempt++) {
        int slotIdx;
        if (weightPrune != PRUNE_UNIFORM) {
            double x = random_double() * totalWeight;
            slotIdx = (int) (upper_bound(cumWeight.begin(), cumWeight.end(), x) - cumWeight.begin());
            if (slotIdx >= (int) reg.slots.size())
                slotIdx = (int) reg.slots.size() - 1; // floating-point fallback; should only trigger on rounding
        } else {
            slotIdx = random_int((int) reg.slots.size());
        }
        pair<PhyloNode*, PhyloNode*> &edge = reg.slots[slotIdx];
        bool aIsDad = edge.first->degree() == 3;
        bool bIsDad = edge.second->degree() == 3;
        if (!aIsDad && !bIsDad)
            continue;

        bool firstHasRoot = tree.rooted && subtreeContainsRoot(tree, edge.first, edge.second);

        if (aIsDad && bIsDad) {
            if (tree.rooted) {
                if (firstHasRoot) { outDad = edge.first; outNode = edge.second; }
                else { outDad = edge.second; outNode = edge.first; }
            } else {
                bool firstIsCherry = isCherryAwayFrom(edge.first, edge.second);
                bool secondIsCherry = isCherryAwayFrom(edge.second, edge.first);
                if (firstIsCherry && !secondIsCherry) {
                    outDad = edge.second; outNode = edge.first;
                } else if (secondIsCherry && !firstIsCherry) {
                    outDad = edge.first; outNode = edge.second;
                } else if (random_int(2) == 0) {
                    outDad = edge.first; outNode = edge.second;
                } else {
                    outDad = edge.second; outNode = edge.first;
                }
            }
        } else if (aIsDad) {
            if (tree.rooted && !firstHasRoot)
                continue; // edge.first as dad would strand root on the pruned side
            outDad = edge.first; outNode = edge.second;
        } else {
            if (tree.rooted && firstHasRoot)
                continue; // edge.second as dad would strand root on the pruned side
            outDad = edge.second; outNode = edge.first;
        }
        return true;
    }
    return false;
}

/**
    resolve one specific SPR-eligible edge (a,b) into a legal (node,dad)
    prune pair, using exactly the same per-edge orientation rule as
    choosePrune() above (a cherry, if only one side is one, is always
    forced to be the "node"/pruned side; otherwise a coin flip) -- except
    the coin flip's outcome is supplied by the caller (`preferA`) instead
    of being drawn internally via random_int(). Needed by
    runBranchLengthCompare(), which must resolve the SAME edge -- named by
    the same slot index -- in two SEPARATE EdgeRegistry/PhyloTree
    instances that started topologically identical (see that function's
    own comment for why): calling choosePrune() independently on each tree
    would consume two different random_int() draws for what is meant to
    be one shared coin flip, letting the two trees disagree on
    orientation even though the exact same edge (by slot index) was
    picked in both.
    @return false if this edge has no legal orientation at all (only
    possible for a genuinely rooted tree, and only for edges next to
    tree.root -- see choosePrune's own comment; this tool's own trees are
    never rooted in practice, so this is effectively unreachable here)
 */
bool resolvePruneOrientationForSlot(PhyloTree &tree, PhyloNode *a, PhyloNode *b, bool preferA,
        PhyloNode* &outNode, PhyloNode* &outDad) {
    bool aIsDad = a->degree() == 3;
    bool bIsDad = b->degree() == 3;
    if (!aIsDad && !bIsDad)
        return false;

    bool firstHasRoot = tree.rooted && subtreeContainsRoot(tree, a, b);

    if (aIsDad && bIsDad) {
        if (tree.rooted) {
            if (firstHasRoot) { outDad = a; outNode = b; }
            else { outDad = b; outNode = a; }
        } else {
            bool firstIsCherry = isCherryAwayFrom(a, b);
            bool secondIsCherry = isCherryAwayFrom(b, a);
            if (firstIsCherry && !secondIsCherry) {
                outDad = b; outNode = a;
            } else if (secondIsCherry && !firstIsCherry) {
                outDad = a; outNode = b;
            } else if (preferA) {
                outDad = a; outNode = b;
            } else {
                outDad = b; outNode = a;
            }
        }
    } else if (aIsDad) {
        if (tree.rooted && !firstHasRoot)
            return false; // edge.first as dad would strand root on the pruned side
        outDad = a; outNode = b;
    } else {
        if (tree.rooted && firstHasRoot)
            return false; // edge.second as dad would strand root on the pruned side
        outDad = b; outNode = a;
    }
    return true;
}

/**
    relative weight for picking graft distance `d` (1..r) in chooseGraft.
    Currently flat (every distance in range equally likely) -- change only
    this function's body to favor closer or farther distances instead
    (e.g. return 1.0 / d to favor close graft points); chooseGraftDistance's
    sampling logic itself never needs to change.
 */
double graftDistanceWeight(int d) {
    return 1.0;
}

/**
    pick a random graft distance in [1,r], weighted by graftDistanceWeight.
 */
int chooseGraftDistance(int r) {
    vector<double> weight(r);
    double total = 0.0;
    for (int d = 1; d <= r; d++) {
        weight[d - 1] = graftDistanceWeight(d);
        total += weight[d - 1];
    }
    double x = random_double() * total;
    double cum = 0.0;
    for (int d = 1; d <= r; d++) {
        cum += weight[d - 1];
        if (x < cum)
            return d;
    }
    return r; // floating-point fallback; should only trigger on rounding
}

/**
    if `raw` (one of `from`'s real neighbors) is pruneDad, returns
    pruneDad's OTHER non-pruneNode neighbor instead of pruneDad itself --
    collapsing the walk's view of pruneDad exactly the way applySPR will
    collapse it for real once the prune actually happens (pruneDad gets
    suppressed, and its two remaining neighbors become directly connected
    to each other). pruneDad's own edges are never a legal regraft target
    (they don't survive the prune -- see isLegalSPR), so rather than let
    the walk stop ON pruneDad and special-case every consumer of
    dad/node/grandDad to tolerate that, this makes pruneDad invisible to
    the walk entirely: "passing through" it is folded into a single
    ordinary hop, as if `from` were already directly connected to whatever
    lies on pruneDad's far side. Returns `raw` unchanged if it isn't
    pruneDad in the first place.

    `from` must actually be adjacent to pruneDad whenever raw==pruneDad
    (true by construction everywhere this is called: it's only invoked on
    a real neighbor of `from` that turned out to be pruneDad, so `from` is
    one of pruneDad's own two non-pruneNode neighbors) -- pruneDad's three
    neighbors are then pruneNode, `from`, and exactly one more, which is
    what gets returned.
 */
PhyloNode* collapsePruneDad(PhyloNode *from, PhyloNode *raw, PhyloNode *pruneDad, PhyloNode *pruneNode) {
    if (raw != pruneDad)
        return raw;
    FOR_NEIGHBOR_IT(pruneDad, from, it)
        if ((*it)->node != pruneNode)
            return (PhyloNode*) (*it)->node;
    ASSERT(0 && "pruneDad always has exactly one non-pruneNode neighbor besides `from`");
    return nullptr;
}

/**
    true if `a` and `b` are adjacent in the collapsed view chooseGraft
    operates on (see collapsePruneDad) -- either really adjacent, or both
    directly incident to pruneDad, which collapsePruneDad papers over as a
    single edge directly between them.
 */
bool virtuallyAdjacent(PhyloNode *a, PhyloNode *b, PhyloNode *pruneDad) {
    return a->isNeighbor(b) || (a->isNeighbor(pruneDad) && b->isNeighbor(pruneDad));
}

/**
    pick a single random SPR regraft target for the (pruneNode,pruneDad)
    edge choosePrune() already selected, via a directed random walk instead
    of enumerating every legal candidate the way findGraftPositions does --
    O(distance walked) rather than O(candidates within radius r).

    The walk operates entirely on a COLLAPSED view of the tree: pruneDad's
    two non-pruneNode edges are treated as a single, ordinary edge directly
    between its two siblings (see collapsePruneDad/virtuallyAdjacent) --
    exactly the topology that will exist once the prune actually happens
    and pruneDad gets suppressed for real. pruneDad itself is therefore
    never visible to the walk at all: it never becomes `dad`, `node`, or
    `grandDad`, so none of the tier logic below needs to know pruneDad
    exists, or reason about how many steps are left before the walk has to
    stop touching it -- every candidate this produces is a real, ordinary
    edge in the tree, full stop. (An earlier version of this function
    tried to allow dad/node to transiently equal pruneDad, gated by how
    many steps remained -- correctly reasoned in isolation, but each fix
    for one illegal case kept exposing another, since pruneDad's
    "everything funnels through one deterministic exit" shape doesn't
    behave like an ordinary node under tier 1's forced-priority walk.
    Collapsing it away entirely sidesteps that whole family of cases
    instead of enumerating them.)

    First picks a target distance d in [1,r] via chooseGraftDistance (see
    graftDistanceWeight to change how distance is weighted). Then walks d
    steps starting from "the edge above the pruned edge" -- the
    (sibling1,sibling2) edge that exists in the collapsed view above, which
    is exactly where a graft-target search should begin (mirrors
    findGraftPositions: pruneDad's own edges are never legal targets).

    Step 1 is a special case: pick uniformly among every real edge incident
    to sibling1 or sibling2 (excluding their shared edge to pruneDad --
    that's the illegal "start edge" itself, never a candidate at all; it
    can't reappear as a candidate later either, since collapsePruneDad
    never returns pruneDad).

    Steps 2..d: only 3 pointers are ever tracked -- grandDad, dad, node --
    where (dad,node) is the current edge and (grandDad,dad) is the edge the
    walk was on immediately before this one (both in the collapsed sense).
    There's no history beyond that single edge: grandDad is overwritten on
    every move (including tier 2's sideways ones), always to whichever node
    was just displaced, so it's always "the edge just come from" for
    whatever the new current edge is -- never a deeper record of how the
    walk got there. Prefer, in order:
      1. node's other edges (walk further outward; up to 2 candidates,
         each passed through collapsePruneDad): grandDad=dad, dad=node,
         node=<pick>
      2. dad's remaining edge, other than the one leading back to grandDad
         (a sideways step onto dad's other, not-yet-explored branch; at
         most 1 candidate after collapsing, and -- since dad is always
         degree 3 here -- effectively always exactly 1): grandDad=node,
         node=<pick> (dad unchanged)
      3. the edge just come from, (grandDad,dad) -- this is not a special
         "undo" move, just the lowest-priority option, tried only when
         neither of the above has anything to offer: dad=grandDad,
         node=dad, grandDad=node (using each variable's value from before
         this reassignment). Only usable when grandDad is actually still
         (virtually) adjacent to dad, which holds right after step 1 and
         right after a tier-1 move, but not in general right after an
         earlier tier-3 move itself (grandDad then refers to whatever node
         just got displaced, two hops from the new dad, not one) -- so two
         tier-3 moves never fire back to back; that's checked directly
         (virtuallyAdjacent(dad, grandDad, pruneDad)) rather than assumed,
         since without it two-in-a-row would fabricate a "current edge"
         between nodes that were never actually (virtually) adjacent
    Ties within whichever tier is chosen are broken uniformly at random.

    Two things are excluded from ever being selected as dad or node at all,
    at every step: the tree's arbitrary root leaf (isLegalSPR always
    rejects grafting onto its edge, regardless of distance -- and since
    that leaf could be anywhere in the tree, not just near the prune point,
    it has to be filtered out of every tier's candidates, not assumed
    unreachable), and pruneNode itself (tier 3 can reach back through the
    collapsed pruneDad -- see above -- and must never continue on into the
    very subtree being pruned). pruneDad itself needs no such check here:
    collapsePruneDad guarantees it can never be produced as a candidate in
    the first place, so there's nothing to filter.

    With all that excluded up front, every candidate this walk can ever
    produce is legal by construction, unlike findGraftPositions which
    instead runs isLegalSPR per candidate after the fact (still ASSERTed
    below anyway, as a cheap defense-in-depth check -- this is exactly how
    an earlier version of this function, which forgot the root-leaf
    exclusion, was caught producing an illegal candidate during testing).

    outDistance, if non-null, receives the walk length d that was chosen
    (not necessarily the true hop-distance of the final edge from the
    prune point -- tier 2/3 moves don't always change hop-distance the way
    a pure tier-1 walk would, and neither does collapsing past pruneDad --
    just the number of random-walk steps taken, for logging purposes).

    @return false if there is no edge to graft onto at all (sibling1 and
    sibling2 are both leaves, e.g. a 3-leaf tree)
 */
bool chooseGraft(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, int r,
        PhyloNode* &outNode, PhyloNode* &outDad, int *outDistance) {
    int d = chooseGraftDistance(r);
    if (outDistance)
        *outDistance = d;

    // grafting onto the edge incident to the tree's arbitrary root leaf is
    // always illegal (see isLegalSPR's last check) regardless of distance
    // from the prune point; since that leaf could be anywhere in the tree,
    // every candidate-building step below excludes it, exactly as if it
    // had no edge to offer at all (it's a real leaf otherwise, so this is
    // the one exclusion beyond the usual tier logic that pruneDad's own
    // collapsing doesn't already take care of -- see pruneNode above)
    PhyloNode *root = (PhyloNode*) tree.root;

    // step 1 (special case): every real edge off sibling1 or sibling2,
    // excluding the edge back to pruneDad. sibling1/sibling2 are also kept
    // around so grandDad can be initialized to "whichever sibling wasn't
    // picked" -- in the collapsed view this function operates in, that's
    // the true far end of "the edge before dad", not pruneDad itself.
    vector<pair<PhyloNode*, PhyloNode*> > firstStep; // (dad=near, node=far)
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    FOR_NEIGHBOR_IT(pruneDad, pruneNode, it) {
        PhyloNode *sibling = (PhyloNode*) (*it)->node;
        if (!sibling1) sibling1 = sibling; else sibling2 = sibling;
        FOR_NEIGHBOR_IT(sibling, pruneDad, it2)
            if ((*it2)->node != root)
                firstStep.push_back(make_pair(sibling, (PhyloNode*) (*it2)->node));
    }
    if (firstStep.empty())
        return false;

    pair<PhyloNode*, PhyloNode*> chosen = firstStep[random_int((int) firstStep.size())];
    PhyloNode *dad = chosen.first;
    PhyloNode *node = chosen.second;
    // grandDad is the other end of the edge the walk was just on -- in the
    // collapsed view, that's dad's sibling on the OTHER side of pruneDad
    // (see this function's own comment above), not pruneDad itself
    PhyloNode *grandDad = (dad == sibling1) ? sibling2 : sibling1;

    for (int step = 2; step <= d; step++) {
        // a candidate produced by actually collapsing pruneDad (raw ==
        // pruneDad, as opposed to collapsePruneDad just handing back an
        // ordinary unrelated neighbor) pairs two of pruneDad's own
        // siblings together -- real nodes, but NOT actually adjacent to
        // each other yet, since pruneDad hasn't been suppressed at the
        // time chooseGraft runs (that only happens later, for real, inside
        // applySPR). Fine as a mid-walk waypoint -- the very next step
        // immediately continues past it to a real edge -- but on the very
        // last step, the walk has to STOP here, and stopping on an edge
        // that doesn't exist yet is exactly what isLegalSPR's own
        // isNeighbor check would (correctly) reject. So collapsed
        // candidates are only excluded on the final step, same as
        // root/pruneNode are excluded always.
        bool isLastStep = (step == d);
        vector<PhyloNode*> tier1, tier2;
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNode *raw = (PhyloNode*) (*it)->node;
            PhyloNode *cand = collapsePruneDad(node, raw, pruneDad, pruneNode);
            if (cand != root && cand != pruneNode && !(isLastStep && raw == pruneDad))
                tier1.push_back(cand);
        }
        FOR_NEIGHBOR_IT(dad, node, it) {
            PhyloNode *raw = (PhyloNode*) (*it)->node;
            PhyloNode *cand = collapsePruneDad(dad, raw, pruneDad, pruneNode);
            if (cand != grandDad && cand != root && cand != pruneNode && !(isLastStep && raw == pruneDad))
                tier2.push_back(cand);
        }

        if (!tier1.empty()) {
            grandDad = dad;
            dad = node;
            node = tier1[random_int((int) tier1.size())];
        } else if (!tier2.empty()) {
            grandDad = node;
            node = tier2[random_int((int) tier2.size())];
        } else if (virtuallyAdjacent(dad, grandDad, pruneDad)
                && (!isLastStep || dad->isNeighbor(grandDad))) {
            // with pruneDad collapsed away, tier 3 is only reached in one
            // irreducible case: `node` is a genuine leaf (tier 1 has
            // nothing to extend to) AND dad's one remaining real neighbor
            // (tier 2's only possible candidate -- dad is always degree 3,
            // and node/grandDad already account for the other two slots)
            // is specifically the root leaf. Nothing else can close off
            // both tiers at once: every internal node still has a real
            // tier-1 option (only root is ever excludable), and tier 2's
            // one slot doesn't care whether grandDad itself is a leaf --
            // that's a property of a DIFFERENT node than the one tier 2
            // actually examines.
            //
            // tier 3 (the edge just come from): only a real candidate if
            // grandDad is still (virtually) adjacent to dad -- which is
            // true right after step 1, and right after a tier-1 move, but
            // NOT in general right after an earlier tier-3 move itself
            // (that reassigns grandDad to whatever node just got
            // displaced, which sits two hops away from the new dad, not
            // one) -- requiring the check here, rather than assuming
            // tier 3 always has a valid destination, is what stops two
            // tier-3 moves in a row from producing a nonexistent "edge"
            // between unrelated nodes (caught by the isLegalSPR ASSERT
            // below during testing). On the last step specifically,
            // virtual adjacency alone isn't enough either, for the same
            // not-a-real-edge-yet reason as the collapsed candidates
            // above -- real adjacency is required there.
            PhyloNode *oldDad = dad, *oldNode = node;
            dad = grandDad;
            node = oldDad;
            grandDad = oldNode;
        }
        // else: no legal move at all this step (only possible in a tiny
        // tree); nothing to do
    }

    outDad = dad;
    outNode = node;

    SPRMove move;
    move.prune_node = pruneNode;
    move.prune_dad = pruneDad;
    move.regraft_node = outNode;
    move.regraft_dad = outDad;
    ASSERT(tree.isLegalSPR(move));

    return true;
}

/**
    <radius>'s meaning for chooseGraftByDistance below: a PERCENTAGE of the
    tree's own total branch length (MTree::treeLength(), summed over every
    edge), not an edge count and not an absolute distance. budget =
    (radiusPercent / 100.0) * tree.treeLength(). This is what makes it
    "normalized": the same <radius> value (say, 10) means "wander up to
    10% of the tree's own total branch length away from the prune point"
    regardless of whether this tree's branch lengths are tiny
    substitutions-per-site values or something on a completely different
    scale -- unlike a fixed absolute-length radius, which would need a
    different value per dataset to mean anything comparable.

    radiusPercent is a double, not an int, even though every value that
    ever reaches it from the command line is a whole number (atoi'd from
    <radius>): "learnradius" (see its own comment on runHillClimb), when
    composed with "distradius", draws a genuinely fractional percentage
    from its own fitted distribution, and needs that fraction to survive
    all the way to the actual budget rather than being rounded away first
    -- the same "continuous" requirement learnRadiusContinuous's own
    comment describes.
 */
double distanceRadiusBudget(PhyloTree &tree, double radiusPercent) {
    return (radiusPercent / 100.0) * tree.treeLength();
}

/**
    alternative to chooseGraft() -- draws a single random SPR regraft
    target the same way (a directed random walk from the collapsed view
    above (pruneNode,pruneDad), never a full enumeration), but governed by
    actual summed BRANCH LENGTH along the walk instead of a fixed hop
    count. Selected via the "distradius" flag (see parseHillClimbFlags);
    everything about legality/candidate-tier priority is identical to
    chooseGraft (same collapsePruneDad/virtuallyAdjacent collapsed view,
    same tier1 forward / tier2 sideways / tier3 backward priority order,
    same root-leaf and pruneNode exclusions, same final isLegalSPR
    ASSERT) -- only the STOPPING rule differs.

    Budget and the up-front deduction for "the branch(es) above the start
    branch": before the walk takes a single step, `remaining` is seeded to
    distanceRadiusBudget(tree, radiusPercent) minus the combined length of
    pruneDad's own two other edges (to sibling1 and sibling2) -- these sit
    directly above the pruned edge and, in the collapsed view this walk
    operates on, are already a single edge of exactly that summed length
    (see collapsePruneDad); step 1 immediately jumps past them to whichever
    real edge it draws, so that distance is spent up front rather than
    charged per-direction.

    Each subsequent hop updates `remaining` by whatever its own tier
    implies, mirroring what real outward distance from the prune point
    that hop actually represents:
      - tier 1 (forward, extends past `node`): pure spend -- charges the
        new edge's own length (or, for a hop that resolves through
        collapsePruneDad, that edge's length PLUS the real edge beyond
        pruneDad it silently continues through -- see the loop body).
      - tier 2 (sideways, pivots at `dad` onto its other branch): the edge
        being abandoned (dad's previous `node`) was never actually
        continued past, so its length is REFUNDED back into `remaining`
        before the new edge's length is charged -- net effect, `remaining`
        only ever reflects the length of edges genuinely still part of the
        current path.
      - tier 3 (backward, undoes the last hop): the edge being left is
        refunded the same way, but the edge walked BACK onto was already
        paid for earlier in the path (at step 1, via the up-front
        deduction above, if backing all the way up to the very first
        collapsed edge; otherwise at whatever earlier tier-1 hop first
        reached `dad`) -- so nothing new is charged, `remaining` just goes
        back to what it was before that earlier charge.
    This is the "sideways/backward steps re-add that edge's distance"
    behavior: those two tiers don't represent real progress away from the
    prune point, so they don't get to spend budget as if they did.

    Stopping: the walk stops the moment `remaining` first drops to zero or
    below AND the edge it is currently on is a real, already-existing edge
    -- i.e. it stops ON the edge that exhausts the budget, not the one
    before it (no "would exceed the budget, so don't take this hop"
    lookahead; whichever edge crosses zero remaining IS the answer). The
    "AND real edge" qualifier matters because a hop that resolves through
    collapsePruneDad lands on an edge that doesn't exist yet (pruneDad
    hasn't actually been suppressed at the time this runs) -- isLegalSPR
    would reject stopping there, so if the budget runs out exactly on such
    a hop, the walk is forced to keep going (ignoring the exhausted budget
    just this once) until it reaches real ground again.

    Because tier2/tier3 refund distance rather than spend it, nothing here
    guarantees `remaining` decreases monotonically the way the hop-count
    version's step counter does -- on a tree with many very-short/
    zero-length branches near the prune point, a long enough run of
    sideways/backward moves could in principle keep refunding as fast as
    it spends, never actually exhausting the budget. maxHops below is a
    hard cap against exactly that: sized off the tree's own edge count
    (generous, practically unreachable on an ordinary walk) rather than
    guessed as a fixed constant, purely so this function is guaranteed to
    terminate regardless of how pathological the branch lengths are. If
    it's ever hit while still sitting on a not-yet-real edge, the walk
    falls back to the most recent REAL edge it was on (always available:
    step 1's own candidate is always real -- see below) rather than
    return something illegal.

    outHops, if non-null, receives the number of hops actually taken (like
    chooseGraft's outDistance, purely for logging -- this is a count of
    walk steps, not a distance).

    outPercentUsed, if non-null, receives the equivalent radiusPercent that
    would have been JUST enough budget to reach the returned edge -- i.e.
    100 * (the actual summed branch length spent, budget minus whatever
    `remaining` stood at when the walk stopped) / tree.treeLength() -- in
    the SAME units as radiusPercent itself, unlike outHops. This is what
    "learnradius" (see its own comment on runHillClimb) records into its
    own history when composed with "distradius": bestDistance/outHops is a
    walk-step count with no fixed relationship to radiusPercent's own
    percentage scale (a handful of short branches or a few long ones can
    both add up to the same budget), so it would be the wrong thing to feed
    back into a distribution that's meant to predict FUTURE radiusPercent
    values.

    @return false if there is no edge to graft onto at all (same
    3-leaf-tree case chooseGraft documents).
 */
bool chooseGraftByDistance(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, double radiusPercent,
        PhyloNode* &outNode, PhyloNode* &outDad, int *outHops, double *outPercentUsed) {
    PhyloNode *root = (PhyloNode*) tree.root;
    double budget = distanceRadiusBudget(tree, radiusPercent);

    vector<pair<PhyloNode*, PhyloNode*> > firstStep; // (dad=near, node=far)
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    double siblingEdgeLength = 0.0;
    FOR_NEIGHBOR_IT(pruneDad, pruneNode, it) {
        PhyloNode *sibling = (PhyloNode*) (*it)->node;
        siblingEdgeLength += (*it)->length;
        if (!sibling1) sibling1 = sibling; else sibling2 = sibling;
        FOR_NEIGHBOR_IT(sibling, pruneDad, it2)
            if ((*it2)->node != root)
                firstStep.push_back(make_pair(sibling, (PhyloNode*) (*it2)->node));
    }
    if (firstStep.empty())
        return false;

    double remaining = budget - siblingEdgeLength;

    pair<PhyloNode*, PhyloNode*> chosen = firstStep[random_int((int) firstStep.size())];
    PhyloNode *dad = chosen.first;
    PhyloNode *node = chosen.second;
    PhyloNode *grandDad = (dad == sibling1) ? sibling2 : sibling1;
    double lastHopLength = dad->findNeighbor(node)->length;
    remaining -= lastHopLength;
    // step 1's own candidate is always a real, already-existing edge (see
    // chooseGraft's own comment: sibling1/sibling2's neighbors, excluding
    // pruneDad, are never collapsed candidates) -- so this is always a
    // legal place to fall back to if the walk below never reaches another
    // one before hitting maxHops
    PhyloNode *lastRealDad = dad, *lastRealNode = node;
    bool currentIsVirtual = false;

    int hops = 1;
    const int maxHops = std::max(1000, 50 * (2 * (int) tree.leafNum - 3));
    bool stopped = (remaining <= 0.0); // already exhausted after step 1 alone

    struct HopCandidate { PhyloNode *node; double length; bool virtualHop; };

    while (!stopped && hops < maxHops) {
        vector<HopCandidate> tier1, tier2;
        FOR_NEIGHBOR_IT(node, dad, it) {
            PhyloNode *raw = (PhyloNode*) (*it)->node;
            PhyloNode *cand = collapsePruneDad(node, raw, pruneDad, pruneNode);
            if (cand == root || cand == pruneNode)
                continue;
            bool virt = (raw == pruneDad);
            double len = (*it)->length + (virt ? pruneDad->findNeighbor(cand)->length : 0.0);
            tier1.push_back(HopCandidate{cand, len, virt});
        }
        FOR_NEIGHBOR_IT(dad, node, it) {
            PhyloNode *raw = (PhyloNode*) (*it)->node;
            PhyloNode *cand = collapsePruneDad(dad, raw, pruneDad, pruneNode);
            if (cand == grandDad || cand == root || cand == pruneNode)
                continue;
            bool virt = (raw == pruneDad);
            double len = (*it)->length + (virt ? pruneDad->findNeighbor(cand)->length : 0.0);
            tier2.push_back(HopCandidate{cand, len, virt});
        }

        bool moved = true;
        if (!tier1.empty()) {
            HopCandidate pick = tier1[random_int((int) tier1.size())];
            grandDad = dad;
            dad = node;
            node = pick.node;
            remaining -= pick.length;
            lastHopLength = pick.length;
            currentIsVirtual = pick.virtualHop;
        } else if (!tier2.empty()) {
            HopCandidate pick = tier2[random_int((int) tier2.size())];
            remaining += lastHopLength; // sideways: refund the abandoned edge
            grandDad = node;
            node = pick.node;
            remaining -= pick.length;
            lastHopLength = pick.length;
            currentIsVirtual = pick.virtualHop;
        } else if (virtuallyAdjacent(dad, grandDad, pruneDad)) {
            bool backEdgeReal = dad->isNeighbor(grandDad);
            double backLen = backEdgeReal
                ? dad->findNeighbor(grandDad)->length
                // still pruneDad's own two (as-yet uncollapsed) siblings --
                // same combined-length convention as siblingEdgeLength above
                : pruneDad->findNeighbor(dad)->length + pruneDad->findNeighbor(grandDad)->length;

            remaining += lastHopLength; // backward: refund the edge being left
            // no charge for backLen -- (grandDad,dad) was already paid for
            // earlier in the path (see this function's own comment)
            PhyloNode *oldDad = dad, *oldNode = node;
            dad = grandDad;
            node = oldDad;
            grandDad = oldNode;
            lastHopLength = backLen;
            currentIsVirtual = !backEdgeReal;
        } else {
            moved = false; // no legal move at all this step (tiny tree)
        }

        if (!moved)
            break;
        hops++;
        if (!currentIsVirtual) {
            lastRealDad = dad;
            lastRealNode = node;
        }
        if (remaining <= 0.0 && !currentIsVirtual)
            stopped = true;
        // else: still have budget, or ran out but landed on an edge that
        // doesn't exist yet -- either way, keep going (maxHops bounds it)
    }

    outDad = currentIsVirtual ? lastRealDad : dad;
    outNode = currentIsVirtual ? lastRealNode : node;
    if (outHops)
        *outHops = hops;
    if (outPercentUsed) {
        double spent = budget - remaining;
        double treeLen = tree.treeLength();
        *outPercentUsed = (treeLen > 0.0) ? (100.0 * spent / treeLen) : 0.0;
    }

    SPRMove move;
    move.prune_node = pruneNode;
    move.prune_dad = pruneDad;
    move.regraft_node = outNode;
    move.regraft_dad = outDad;
    ASSERT(tree.isLegalSPR(move));

    return true;
}

/**
    apply an SPR move and keep `reg` in sync with the result, in O(1).

    Every applySPR call destroys exactly 3 edges and creates exactly 3 new
    ones (see PhyloTree::applySPR in phylotree.cpp): pruning detaches dad1
    from node1 by directly connecting dad1's other two neighbors --
    "sibling1" and "sibling2" -- to each other (destroying edges
    dad1-sibling1 and dad1-sibling2, creating sibling1-sibling2); then dad1
    itself (never deleted -- applySPR keeps the same node object and just
    repoints its now-spare neighbor slots) is spliced into the regraft edge
    (destroying node2-dad2, creating dad1-dad2 and dad1-node2). No node's
    degree ever changes anywhere in this -- dad1 stays degree 3 throughout,
    and sibling1/sibling2/node2/dad2 each just swap which node one
    particular neighbor slot points at -- so the registry never needs a
    slot added or removed, only these same 3 existing slots repointed to
    their new edge.

    The one exception: if dad1 happens to be currently adjacent to the
    tree's root leaf, one of "sibling1"/"sibling2" above IS that root leaf.
    dad1's edge to it is never in the registry to begin with (buildEdgeRegistry
    excludes root's one edge by design -- see its comment), and the edge
    this move creates in its place -- from root to dad1's OTHER sibling --
    becomes the new such excluded edge, so it must not be tracked either.
    In that case only 2 of the usual 3 edges have real slots at all, so
    replaceEdges2 (not replaceEdges3) is used instead.

    sibling1/sibling2 are captured here, before the move, since afterward
    dad1 is no longer adjacent to them.
 */
void applySPRTracked(PhyloTree &tree, EdgeRegistry &reg, const SPRMove &move, TrackedSPR &t) {
    t.move = move;
    PhyloNode *dad1 = move.prune_dad;
    PhyloNode *node1 = move.prune_node;
    PhyloNode *node2 = move.regraft_node;
    PhyloNode *dad2 = move.regraft_dad;
    PhyloNode *root = (PhyloNode*) tree.root;

    t.sibling1 = t.sibling2 = nullptr;
    FOR_NEIGHBOR_IT(dad1, node1, it) {
        if (!t.sibling1)
            t.sibling1 = (PhyloNode*) (*it)->node;
        else
            t.sibling2 = (PhyloNode*) (*it)->node;
    }

    t.sib1IsRoot = (t.sibling1 == root);
    t.sib2IsRoot = (t.sibling2 == root);
    ASSERT(!(t.sib1IsRoot && t.sib2IsRoot)); // root can't be both of dad1's siblings

    t.slotDad1Sib1 = t.sib1IsRoot ? -1 : reg.slotOf(dad1, t.sibling1);
    t.slotDad1Sib2 = t.sib2IsRoot ? -1 : reg.slotOf(dad1, t.sibling2);
    t.slotNode2Dad2 = reg.slotOf(node2, dad2);

    tree.applySPR(move, t.rollback);

    if (!t.sib1IsRoot && !t.sib2IsRoot) {
        reg.replaceEdges3(t.slotDad1Sib1, dad1, dad2,
                          t.slotDad1Sib2, dad1, node2,
                          t.slotNode2Dad2, t.sibling1, t.sibling2);
    } else {
        int nonRootSlot = t.sib1IsRoot ? t.slotDad1Sib2 : t.slotDad1Sib1;
        reg.replaceEdges2(nonRootSlot, dad1, dad2,
                          t.slotNode2Dad2, dad1, node2);
    }
}

/**
    undo an applySPRTracked call: rolls back the tree exactly as
    tree.rollbackSPR always did, and restores `reg`'s tracked slots to
    exactly what they held beforehand (see applySPRTracked for why this is
    sometimes 2 slots, not 3).
 */
void rollbackSPRTracked(PhyloTree &tree, EdgeRegistry &reg, const TrackedSPR &t) {
    tree.rollbackSPR(t.rollback);

    PhyloNode *dad1 = t.move.prune_dad;
    if (!t.sib1IsRoot && !t.sib2IsRoot) {
        reg.replaceEdges3(t.slotDad1Sib1, dad1, t.sibling1,
                          t.slotDad1Sib2, dad1, t.sibling2,
                          t.slotNode2Dad2, t.move.regraft_node, t.move.regraft_dad);
    } else {
        int nonRootSlot = t.sib1IsRoot ? t.slotDad1Sib2 : t.slotDad1Sib1;
        PhyloNode *nonRootSibling = t.sib1IsRoot ? t.sibling2 : t.sibling1;
        reg.replaceEdges2(nonRootSlot, dad1, nonRootSibling,
                          t.slotNode2Dad2, t.move.regraft_node, t.move.regraft_dad);
    }
}

/**
    invalidate cached partial likelihoods after a topology change (an
    applySPR or rollbackSPR call). A plain clearAllPartialLH() is not
    enough here: it only marks existing per-edge likelihood buffers as
    stale, relying on IQ-TREE's incremental buffer-reuse bookkeeping
    (PhyloTree::reorientPartialLh / mem_slots) to still correctly track
    which buffer belongs to which edge. That bookkeeping assumes topology
    changes only ever touch edges adjacent to an already-evaluated branch,
    which is true for NNI but not for SPR (a regraft can jump to a
    distant, never-yet-traversed edge), and silently building on a stale
    slot assignment there crashes with either a "no buffer to reorient"
    assertion or a buffer-aliasing assertion depending on the memory mode.
    A full delete+reinitialize sidesteps the incremental bookkeeping
    entirely by rebuilding every buffer assignment from scratch to match
    the current topology; heavier per candidate, but unambiguously
    correct, matching SPR_IMPLEMENTATION_PLAN's guidance to start from
    full invalidation before ever optimizing cache reuse.
 */
void resetLikelihoodBuffers(PhyloTree &tree) {
    tree.deleteAllPartialLh();
    tree.initializeAllPartialLh();
}

void allocateSPRLocalLhCache(PhyloTree &tree, SPRLocalLhCache &cache) {
    cache.lhCount = tree.getPartialLhSize();
    cache.scaleCount = tree.getScaleNumSize();
}

void growSPRLocalLhCache(SPRLocalLhCache &cache, size_t count) {
    while (cache.lh.size() < count) {
        cache.lh.push_back(aligned_alloc<double>(cache.lhCount));
        cache.scaleNum.push_back(aligned_alloc<UBYTE>(cache.scaleCount));
    }
}

void freeSPRLocalLhCache(SPRLocalLhCache &cache) {
    for (size_t i = 0; i < cache.lh.size(); i++) {
        aligned_free(cache.lh[i]);
        aligned_free(cache.scaleNum[i]);
    }
}

/**
    swap in a dedicated scratch buffer (from `cache`) for every affected
    direction, marking each dirty so the next
    computeLikelihoodBranch call is forced to recompute it rather than
    trust whatever content is already sitting in the scratch buffer from a
    previous candidate's use of the same physical slot. Must be called
    AFTER tree.applySPR() -- dad1/dad2/node2/sibling1/sibling2 name the
    POST-move adjacency (dad1 is now adjacent to dad2, node2, and
    unchanged node1; sibling1 is now adjacent to sibling2 directly).
 */
void beginLocalSPRInvalidation(SPRLocalLhCache &cache, SPRLocalInvalidationState &state,
        const vector<PhyloNeighbor*> &affected) {
    growSPRLocalLhCache(cache, affected.size());
    state.nei = affected;
    state.savedLh.resize(affected.size());
    state.savedScaleNum.resize(affected.size());
    state.savedScaleFactor.resize(affected.size());
    state.savedComputed.resize(affected.size());
    for (size_t i = 0; i < affected.size(); i++) {
        double *lh = cache.lh[i];
        UBYTE *sn = cache.scaleNum[i];
        double sf = 0.0;
        int computed = 0;
        state.nei[i]->swapPartialLhState(lh, sn, sf, computed);
        // swapPartialLhState mutated lh/sn/sf/computed in place to hold
        // whatever was on `state.nei[i]` immediately before this call --
        // exactly what endLocalSPRInvalidationDiscard needs to restore.
        state.savedLh[i] = lh;
        state.savedScaleNum[i] = sn;
        state.savedScaleFactor[i] = sf;
        state.savedComputed[i] = computed;
    }
}

/**
    Compute the candidate likelihood anchored on one of the newly split
    target edges. The traversal follows the dirty removal-to-insertion path;
    unaffected side-subtree partials remain cached.
 */
double computeLocalSPRLikelihood(PhyloTree &tree, PhyloNode *dad1, PhyloNode *dad2) {
    return tree.computeLikelihoodBranch((PhyloNeighbor*) dad1->findNeighbor(dad2), dad1);
}

/**
    undo a beginLocalSPRInvalidation call: every affected direction
    goes back to exactly its pre-call (partial_lh, scale_num,
    lh_scale_factor, partial_lh_computed), including back to null/
    uncomputed if that's what it was. A provably perfect round trip --
    swapPartialLhState is its own inverse given the same saved values.
 */
void endLocalSPRInvalidationDiscard(SPRLocalInvalidationState &state) {
    for (size_t i = 0; i < state.nei.size(); i++) {
        double *lh = state.savedLh[i];
        UBYTE *sn = state.savedScaleNum[i];
        double sf = state.savedScaleFactor[i];
        int computed = state.savedComputed[i];
        state.nei[i]->swapPartialLhState(lh, sn, sf, computed);
    }
}

/**
    dad1's other two neighbors, BEFORE an SPR move is applied to
    (node1, dad1) -- these are exactly the two nodes that end up directly
    connected to each other once dad1 is bypassed (see applySPR's own
    logic in phylotree.cpp: sibling1->updateNeighbor(dad1, sibling2, ...)
    and vice versa). Needed by scoreTrialSPRMove's branch-length
    reoptimization below, to name the third touched edge (the merged
    sibling1-sibling2 edge) that SPRMove itself doesn't otherwise
    identify; the test tool has no other way to name these two nodes once
    the move has already been applied -- SPRRollback records enough
    per-Neighbor state to UNDO a move, but never exposes which NODE owned
    which saved Neighbor.
 */
void findSPRSiblings(PhyloNode *node1, PhyloNode *dad1, PhyloNode *&sibling1, PhyloNode *&sibling2) {
    sibling1 = sibling2 = nullptr;
    FOR_NEIGHBOR_DECLARE(dad1, node1, it) {
        if (!sibling1)
            sibling1 = (PhyloNode*) (*it)->node;
        else
            sibling2 = (PhyloNode*) (*it)->node;
    }
}

/**
    apply `move` as a trial, score it, and always roll it back (the
    caller decides afterward whether to actually keep it for real, via
    applySPRTracked). If reoptimizeBranchLengths is set, four edges around
    the move -- the two split halves of the target edge at the new
    attachment point (prune_dad-regraft_dad, prune_dad-regraft_node), the
    merged edge left behind at the vacated attachment point
    (sibling1-sibling2), and the pruned subtree's own upper edge
    (prune_dad-prune_node, unchanged by applySPR itself -- see its own
    comment in phylotree.cpp -- but now sitting in a different part of the
    tree, with different neighbors on prune_dad's other two sides, so its
    old length is no longer necessarily optimal either) -- are each
    re-optimized via optimizeOneBranch() (Newton-Raphson, same as
    PhyloTree::getBestNNIForBran does for the branches it touches) before
    the final score is taken, instead of just trusting applySPR's own
    naive placeholder lengths for the first three (half the target edge's
    length split evenly, and the sum of the two vacated edges' lengths --
    see applySPR's own comment in phylotree.cpp) and prune_dad-prune_node's
    own pre-move length for the fourth. This can only ever improve (or
    leave unchanged) the likelihood for a given topology, since it's
    searching the exact same space optimizeOneBranch always searches, and
    turns the previously-reported score for a topology from "whatever the
    naive placeholder/pre-move lengths happen to give" into an actual
    (locally) likelihood-optimal set of lengths for those four edges --
    closer to what a real ML search would report for the same topology.

    Deliberately built on the tool's original resetLikelihoodBuffers()
    (full delete+reinit) rather than the selective, incremental
    invalidation tried and abandoned above scoreTrialSPRMove's own
    definition: optimizeOneBranch's default clearLH=true already performs
    its own correct, cascading invalidation internally (the exact
    mechanism NNI itself relies on), so as long as this starts from a
    known-fully-valid state (one plain computeLikelihood() call,
    immediately after applySPR) there is no need to hand-manage buffer
    invalidation here at all, and thus no repeat of the LM_PER_NODE
    buffer-pool crash that ended that attempt.
 */
/**
    optimizeOneBranch() asserts its branch's CURRENT length is >= 0 before
    searching for a better one; BioNJ (a distance-based method) can and
    does sometimes produce a zero or slightly negative branch length
    estimate (a well-known artifact of neighbor-joining-family methods on
    noisy distance matrices), and applySPR's own naive length arithmetic
    (summing/halving whatever length was already there) only preserves
    that sign, never corrects it. Clamp both directions of (a,b) up to
    the model's own minimum branch length first so the precondition
    always holds, regardless of what the naive placeholder length was.
 */
void clampBranchLengthForOptimization(PhyloNode *a, PhyloNode *b, double minLen) {
    Neighbor *ab = a->findNeighbor(b);
    Neighbor *ba = b->findNeighbor(a);
    if (ab->length < minLen)
        ab->length = minLen;
    if (ba->length < minLen)
        ba->length = minLen;
}

/**
    clamp every branch length in the tree up to `minLen`, both directions
    of every edge. optimizeAllBranches() calls optimizeOneBranch on
    literally every edge in the tree, and optimizeOneBranch asserts its
    branch's CURRENT length is >= 0 before searching for a better one --
    BioNJ (a distance-based method) can and does sometimes produce a zero
    or slightly negative estimate on some edge (see
    clampBranchLengthForOptimization's comment above, which this reuses).
    Unlike that SPR-move-specific helper, which only needs to clamp the 3
    edges a given move touches, a whole-tree optimization pass needs every
    edge safe first, since any of them could be the offending one.
 */
void clampAllBranchLengthsForOptimization(PhyloTree &tree, double minLen) {
    NodeVector nodes1, nodes2;
    tree.getBranches(nodes1, nodes2);
    for (size_t i = 0; i < nodes1.size(); i++)
        clampBranchLengthForOptimization((PhyloNode*) nodes1[i], (PhyloNode*) nodes2[i], minLen);
}

/**
    re-optimize (Newton-Raphson) the 4 edges around an SPR move --
    prune_dad-regraft_dad and prune_dad-regraft_node (the two split halves
    of the target edge at the new attachment point), the merged
    sibling1-sibling2 edge left behind at the vacated attachment point,
    and prune_dad-prune_node (the pruned subtree's own upper edge: never
    touched by applySPR itself, since prune_dad stays adjacent to
    prune_node throughout the move, but now sitting between a different
    pair of neighbors on prune_dad's other two sides, so its pre-move
    length is no longer necessarily the locally optimal one either) -- in
    place, on whatever tree/topology is CURRENTLY applied when this is
    called. Used both for scoring a trial candidate (scoreTrialSPRMove,
    where the tree is rolled back afterward regardless) and, separately,
    on the real, permanently-kept tree once a candidate has actually won
    a step (see the main loop in runHillClimb) -- calling this again on
    the real tree is what makes the optimized lengths actually persist,
    rather than being found once for scoring purposes and then discarded
    by scoreTrialSPRMove's own rollback.

    Silences cout for the duration: Newton-Raphson search on a candidate
    that's a poor topological fit (routine for "fast" mode's unfiltered
    random draws, which -- unlike exhaustive mode -- never screens out
    implausible candidates before reaching here) can still hit a non-
    finite likelihood derivative even with the safe/scaled kernel
    (Params::lk_safe_scaling, enabled by runHillClimb whenever this
    feature is on). IQ-TREE's own handling for that is non-fatal (the bad
    derivative is zeroed and the NR search just stops where it is), but
    it also prints a "WARNING: Numerical underflow..." line straight to
    cout every time, bypassing this tool's own "quiet" flag entirely
    since it's the library's own diagnostic output, not this tool's.
 */
void reoptimizeSPREdges(PhyloTree &tree, PhyloNode *dad1, PhyloNode *dad2, PhyloNode *node2, PhyloNode *node1,
        PhyloNode *sibling1, PhyloNode *sibling2) {
    const int maxNRStep = 10; // matches NNI_MAX_NR_STEP's default, utils/pllnni.cpp
    double minLen = Params::getInstance().min_branch_length;

    ostringstream suppressedOptimizeOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedOptimizeOutput.rdbuf());

    clampBranchLengthForOptimization(dad1, dad2, minLen);
    tree.optimizeOneBranch(dad1, dad2, true, maxNRStep);

    clampBranchLengthForOptimization(dad1, node2, minLen);
    tree.optimizeOneBranch(dad1, node2, true, maxNRStep);

    // prune_dad-prune_node itself: never destroyed/recreated by applySPR
    // (dad1 stays adjacent to node1 throughout -- see applySPR's own
    // comment in phylotree.cpp), so unlike the two edges just above there
    // is no naive placeholder length to replace here, only a pre-move one
    // that's no longer necessarily optimal now that dad1's other two
    // neighbors have changed
    clampBranchLengthForOptimization(dad1, node1, minLen);
    tree.optimizeOneBranch(dad1, node1, true, maxNRStep);

    clampBranchLengthForOptimization(sibling1, sibling2, minLen);
    tree.optimizeOneBranch(sibling1, sibling2, true, maxNRStep);

    cout.rdbuf(realCoutBuf);
}

double scoreTrialSPRMoveFullReset(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths) {
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    if (reoptimizeBranchLengths)
        findSPRSiblings(move.prune_node, move.prune_dad, sibling1, sibling2);

    SPRRollback rollback;
    tree.applySPR(move, rollback);
    resetLikelihoodBuffers(tree);
    double score = tree.computeLikelihood();

    if (reoptimizeBranchLengths) {
        reoptimizeSPREdges(tree, move.prune_dad, move.regraft_dad, move.regraft_node, move.prune_node,
                sibling1, sibling2);

        resetLikelihoodBuffers(tree);
        double reoptScore = tree.computeLikelihood();

        // guard against the rare case where repeated underflow left the
        // branch lengths in a degenerate state and the final likelihood
        // itself (not just an intermediate derivative) comes out non-
        // finite -- never let a corrupted score win a "keep the best"
        // comparison; the un-reoptimized `score` computed above is still
        // a valid, finite fallback
        if (std::isfinite(reoptScore))
            score = reoptScore;
    }

    tree.rollbackSPR(rollback);
    resetLikelihoodBuffers(tree);
    return score;
}

/**
    path of nodes from `node` to `goal`, never passing through `blocked`
    (the pruned subtree's own root, which must never be crossed since it's
    leaving the removal-to-insertion path entirely) -- used by
    collectSPRAffectedBeforeApply below to name every persistent edge whose
    cached partial likelihood needs invalidating.
 */
bool findSPRNodePath(PhyloNode *node, PhyloNode *dad, PhyloNode *goal, PhyloNode *blocked,
        vector<PhyloNode*> &path) {
    if (node == blocked)
        return false;
    path.push_back(node);
    if (node == goal)
        return true;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (findSPRNodePath((PhyloNode*) (*it)->node, node, goal, blocked, path))
            return true;
    }
    path.pop_back();
    return false;
}

void addUniqueSPRAffected(vector<PhyloNeighbor*> &affected, PhyloNeighbor *nei) {
    // Tip likelihoods come directly from the alignment and never use an
    // internal partial-likelihood buffer.
    if (nei->node->isLeaf())
        return;
    if (find(affected.begin(), affected.end(), nei) == affected.end())
        affected.push_back(nei);
}

vector<PhyloNeighbor*> collectSPRAffectedBeforeApply(const SPRMove &move) {
    vector<PhyloNode*> pathToDad, pathToNode;
    ASSERT(findSPRNodePath(move.prune_dad, nullptr, move.regraft_dad, move.prune_node, pathToDad));
    ASSERT(findSPRNodePath(move.prune_dad, nullptr, move.regraft_node, move.prune_node, pathToNode));
    const vector<PhyloNode*> &path = pathToDad.size() < pathToNode.size() ? pathToDad : pathToNode;

    vector<PhyloNeighbor*> affected;
    addUniqueSPRAffected(affected,
            (PhyloNeighbor*) move.prune_node->findNeighbor(move.prune_dad));
    // path[0]-path[1] is removed when prune_dad is suppressed. The merged
    // sibling edge is added after apply; every later edge persists, and its
    // direction back toward the removal point is the cached subtree that
    // changes when the pruned clade moves across it.
    for (size_t i = 2; i < path.size(); i++)
        addUniqueSPRAffected(affected,
                (PhyloNeighbor*) path[i]->findNeighbor(path[i - 1]));
    return affected;
}

void addSPRChangedEdgeDirections(vector<PhyloNeighbor*> &affected,
        PhyloNode *dad1, PhyloNode *dad2, PhyloNode *node2,
        PhyloNode *sibling1, PhyloNode *sibling2) {
    PhyloNode *a[6] = {dad1, dad2, dad1, node2, sibling1, sibling2};
    PhyloNode *b[6] = {dad2, dad1, node2, dad1, sibling2, sibling1};
    for (int i = 0; i < 6; i++)
        addUniqueSPRAffected(affected, (PhyloNeighbor*) a[i]->findNeighbor(b[i]));
}

double scoreTrialSPRMove(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths,
        SPRLocalLhCache *localCache) {
    if (!localCache || reoptimizeBranchLengths)
        return scoreTrialSPRMoveFullReset(tree, move, reoptimizeBranchLengths);

    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    findSPRSiblings(move.prune_node, move.prune_dad, sibling1, sibling2);
    vector<PhyloNeighbor*> affected = collectSPRAffectedBeforeApply(move);

    SPRRollback rollback;
    tree.applySPR(move, rollback);
    addSPRChangedEdgeDirections(affected, move.prune_dad, move.regraft_dad, move.regraft_node,
            sibling1, sibling2);

    SPRLocalInvalidationState state;
    beginLocalSPRInvalidation(*localCache, state, affected);
    double score = computeLocalSPRLikelihood(tree, move.prune_dad, move.regraft_dad);
    endLocalSPRInvalidationDiscard(state);

    tree.rollbackSPR(rollback);
    return score;
}

/**
    build a short, descriptive identifier for one runHillClimb call, used
    by "record" (see appendRecordRow) to tell repeated runs in the same
    spreadsheet apart: a timestamp (so identical-flag runs are still
    distinguishable) followed by the flags/parameters that were actually
    used, in the same order printHillClimbFlags prints them -- only the
    ones that are non-default are included, e.g.
    "20260730-121553_r10_s30_fast5_reopt_fullreopt100x5_gtr_investigate3" or
    "20260730-124501_r8_s40_fast_investigate3_alternate".
 */
string buildRunId(int radius, int maxSteps, bool randomStart,
        bool useFastSelection, int numCandidates, bool reoptimizeBranchLengths,
        int fullReoptEveryNSteps, int fullReoptRounds, bool fullReoptInitialFit, bool useGtrModel,
        bool investigateFlag, int investigateRadius, bool alternateFlag, bool shrinkFlag,
        int shrinkStallThreshold, bool learnradiusFlag, int learnradiusN, bool sweepFlag, int sweepCount,
        int findoptEveryNSteps, bool iqtreeStart, int iqtreeStartPoolSize, bool weightpruneFlag) {
    time_t now = time(nullptr);
    char timestamp[32];
    strftime(timestamp, sizeof(timestamp), "%Y%m%d-%H%M%S", localtime(&now));

    ostringstream id;
    id << timestamp << "_r" << radius;
    id << "_s" << maxSteps;
    if (randomStart)
        id << "_random";
    if (iqtreeStart)
        id << "_iqtreestart" << iqtreeStartPoolSize;
    if (useFastSelection)
        id << "_fast" << (numCandidates > 1 ? to_string(numCandidates) : "");
    if (reoptimizeBranchLengths)
        id << "_reopt";
    if (fullReoptEveryNSteps > 0) {
        id << "_fullreopt" << fullReoptRounds << "x" << fullReoptEveryNSteps;
        if (fullReoptInitialFit)
            id << "init";
    }
    if (useGtrModel)
        id << "_gtr";
    if (investigateFlag)
        id << "_investigate" << investigateRadius;
    if (alternateFlag)
        id << "_alternate";
    if (shrinkFlag)
        id << "_shrink" << shrinkStallThreshold;
    if (learnradiusFlag)
        id << "_learnradius" << learnradiusN;
    if (sweepFlag)
        id << "_sweep" << sweepCount;
    if (findoptEveryNSteps > 0)
        id << "_findopt" << findoptEveryNSteps;
    if (weightpruneFlag)
        id << "_weightprune";
    return id.str();
}

/**
    build a short, filename-safe suffix identifying which OPTIONAL
    search-mode flags were used for this run -- fast, reopt,
    investigate -- for recordSpreadsheetPath to route
    differently-tagged runs into separate files instead of this tool's old
    behavior of lumping every run under one model into a single
    record_<model>.csv regardless of which search algorithm actually
    produced it. A plain steepest-descent run's "candidates evaluated" and
    per-step cost have essentially nothing in common with a very different
    search mode's, for instance, so mixing them into one spreadsheet made
    side-by-side comparison misleading rather than useful.

    randomStart is deliberately NOT included: it only affects the
    STARTING topology, not the search mechanics themselves, so a
    random-start and a BioNJ-start run under otherwise-identical flags are
    still a fair, meaningful comparison in the same file. useGtrModel is
    not included either: it already changes modelName itself (to
    "GTR+FO" instead of "JC"), which recordSpreadsheetPath already routes
    to a different file on that basis alone. Numeric sub-parameters
    (numCandidates, fullReoptEveryNSteps, radius, max-steps) are also left
    out of the filename on purpose -- those still vary meaningfully
    *within* a given search-mode file, distinguished via each row's own
    run_id, rather than needing yet another file per value.

    investigateRadius is the one exception to that "numeric parameters
    stay out of the filename" rule: it's included in the tag itself (e.g.
    "_investigate1", "_investigate3"), not just in each row's run_id,
    since it can change an investigation step's cost and behavior
    dramatically -- a radius-1 (NNI-equivalent) refinement and a much
    larger one are different enough searches that mixing their
    trajectories into one file would be misleading.

    alternateFlag is included too (bare "_alternate", no number -- it has
    no numeric parameter of its own): a run that spends half its steps on
    an NNI-equivalent search interleaved with the other mode has a
    meaningfully different cost/behavior profile than a run of that same
    other mode alone.

    shrinkFlag does NOT follow investigateRadius' exception: it's included
    bare ("_shrink", no threshold number) rather than with
    shrinkStallThreshold appended. Different threshold values still shrink
    the SAME way -- a stagnation counter narrowing the same radius
    schedule, just triggered sooner or later -- so they belong side by
    side in one file, distinguished by run_id, the same way different
    fullReoptEveryNSteps values already are. What changes the
    cost/behavior profile enough to need its own file is shrinkFlag itself
    (adaptive narrowing vs. a fixed radius), not the particular threshold
    chosen.

    sweepFlag follows investigateRadius' exception, not alternateFlag's/
    shrinkFlag's convention: sweepCount is included in the tag itself (e.g.
    "_sweep5", "_sweep20"), not just in each row's run_id. Unlike a plain
    on/off flag, sweepFlag now runs as an added post-processing phase whose
    own cost scales directly with sweepCount (each of its N targeted
    positions gets its own exhaustive, whole-tree regraft search) -- a
    "sweep 3" run and a "sweep 50" run add a very different amount of extra
    work on top of the same step loop, different enough to mix misleadingly
    in one file.

    findoptEveryNSteps (like fullReoptEveryNSteps just above) follows
    shrinkFlag's convention: it's included bare ("_findopt", the interval
    itself left to run_id via buildRunId). Unlike the old, single-shot
    "finalreopt" this replaces -- which DID change what the trajectory's
    last row meant, since it overwrote curScore for real -- findopt never
    touches curScore or the main tree at all (see maybeRunFindopt's
    comment), so mixing findopt and non-findopt runs in one file would no
    longer misrepresent anything; it still earns its own file for a
    simpler reason: findopt interleaves EXTRA rows into the trajectory
    (each one a hypothetical "what if we fully refit right now" reading,
    not a real accepted step) that would otherwise look like genuine
    search progress to anything plotting the file's logL column against
    row order.

    fullReoptEveryNSteps follows shrinkFlag's convention, not
    investigateRadius'/sweepFlag's: it's included bare ("_fullreopt", no M
    or N) rather than with either number appended -- both are still left
    to run_id (via buildRunId) the same as any other numeric
    sub-parameter. What earns "fullreopt" its own file at all is that it's
    now a fully independent flag from reoptimizeBranchLengths (unlike its
    old life as "reopt"'s own optional trailing number): a run using it
    alone, with reoptimizeBranchLengths off, has a periodic whole-tree
    branch-length refit that a plain run has no equivalent of at all, so
    it belongs in its own file rather than either the untagged default one
    or "_reopt"'s. Composes with reoptimizeBranchLengths freely, the same
    way sweepFlag composes with everything else -- a run using both gets
    "_reopt_fullreopt".

    fullReoptInitialFit is NOT included here at all, the same way
    randomStart isn't: like randomStart, it only affects the STARTING
    point (whether the initial tree gets one whole-tree branch-length fit
    up front before the step loop begins), not the search's own per-step
    mechanics -- so a "fullreopt M N" run and a "fullreopt M N true" run
    under otherwise-identical flags are still a fair, meaningful
    comparison in the same file, distinguished by run_id (via buildRunId,
    which DOES include it) rather than needing yet another file.

    learnradiusFlag follows shrinkFlag's own convention exactly (it's
    mutually exclusive with shrinkFlag in the first place -- see
    parseHillClimbFlags): included bare ("_learnradius", no window size)
    rather than with learnradiusN appended. Different window sizes still
    learn the SAME way -- a sliding-window Gamma fit determining the same
    radius schedule, just refit over a shorter or longer history -- so they
    belong side by side in one file, distinguished by run_id, the same way
    shrinkStallThreshold's different values do.

    useDistanceRadius ("distradius") is included too (bare "_distradius",
    tucked right after "_fast" since it's specifically "fast" mode's own
    candidate-draw mechanism, with no effect otherwise -- see its own
    comment on parseHillClimbFlags): it swaps chooseGraft's fixed-hop-count
    random walk for chooseGraftByDistance's summed-branch-length one, a
    genuinely different candidate-generation mechanism with its own
    cost/behavior profile (see chooseGraftByDistance's own comment), not
    just a different value of an existing parameter -- exactly the kind of
    difference every other tag in this function already exists to keep
    separated. No number of its own to append (radiusPercent still varies
    *within* the file via run_id's own "_r" field, the same as a plain hop
    radius does), so it follows shrinkFlag's/learnradiusFlag's bare-tag
    convention, not investigateRadius'/sweepCount's.

    noTrueTree ("notree") is included last, independent of every search-mode
    tag above it (it composes with any of them, unlike e.g. distradius which
    only means anything under "fast"): unlike randomStart/useGtrModel/
    fullReoptInitialFit -- deliberately EXCLUDED above because they only
    affect the STARTING point or are already captured via modelName, never
    the recorded trajectory's own shape -- noTrueTree changes what a
    row's own "gap to true tree" column (appendRecordRow's trueTreeLogl -
    logL) actually MEANS: with no ground-truth tree at all, that column is
    written as NaN on every row (see appendRecordRow's own comment) instead
    of a real, meaningful gap value. Mixing noTrueTree and normal runs in
    one file would make that column's own data silently inconsistent --
    real numbers on some rows, NaN on others -- exactly the kind of
    misleading mix this function exists to prevent, so it earns its own
    file the same way findopt's extra, differently-shaped rows do.

    weightpruneFlag ("weightprune") is included too, bare ("_weightprune",
    no number -- it has none of its own), the same way alternateFlag is:
    it replaces choosePrune's uniform edge pick with one weighted by each
    edge's own branch length (see choosePrune's own comment), a genuinely
    different prune-edge selection bias with its own trajectory shape, not
    just a different value of an existing parameter. Independent of every
    other tag above it -- it composes freely with any of them, since it
    only ever changes which edge gets pruned, never how the resulting
    graft search itself proceeds.
 */
string buildRecordTag(bool useFastSelection, bool useDistanceRadius, bool reoptimizeBranchLengths,
        int fullReoptEveryNSteps, bool investigateFlag, int investigateRadius, bool alternateFlag,
        bool shrinkFlag, bool learnradiusFlag, bool sweepFlag, int sweepCount, int findoptEveryNSteps,
        bool noTrueTree, bool weightpruneFlag, bool tunnelFlag, double tunnelTolerance,
        bool slackFlag, double slackDelta, bool slackAnneal) {
    ostringstream tag;
    if (useFastSelection)
        tag << "_fast";
    if (useDistanceRadius)
        tag << "_distradius";
    if (reoptimizeBranchLengths)
        tag << "_reopt";
    if (fullReoptEveryNSteps > 0)
        tag << "_fullreopt";
    if (investigateFlag)
        tag << "_investigate" << investigateRadius;
    if (alternateFlag)
        tag << "_alternate";
    if (shrinkFlag)
        tag << "_shrink";
    if (learnradiusFlag)
        tag << "_learnradius";
    if (sweepFlag)
        tag << "_sweep" << sweepCount;
    if (findoptEveryNSteps > 0)
        tag << "_findopt";
    if (noTrueTree)
        tag << "_notree";
    if (weightpruneFlag)
        tag << "_weightprune";
    // Tolerance is part of the identity, not just the flag: "tunnel 0.5" and
    // "tunnel 20" are different searches and must not share a record file.
    // '.' would break the "_"-delimited tag convention, so it becomes 'p'.
    if (slackFlag) {
        ostringstream d;
        d << slackDelta;
        string t = d.str();
        for (size_t i = 0; i < t.size(); i++)
            if (t[i] == '.') t[i] = 'p';
        tag << "_slack" << t << (slackAnneal ? "anneal" : "");
    }
    if (tunnelFlag) {
        ostringstream tol;
        tol << tunnelTolerance;
        string t = tol.str();
        for (size_t i = 0; i < t.size(); i++)
            if (t[i] == '.') t[i] = 'p';
        tag << "_tunnel" << t;
    }
    return tag.str();
}

/**
    append one row -- this run's id, how many candidates have been
    evaluated (scoreTrialSPRMove calls) so far, CPU-clock seconds elapsed
    since this run started, the current (just-updated) logL, and how far
    that logL still is below the true AliSim tree's own logL
    (trueTreeLogl - logL; positive as long as the search hasn't caught up
    to -- or, since trueTreeLogl isn't necessarily the ceiling, possibly
    surpassed -- the true tree, since curScore should generally trend
    toward trueTreeLogl as the search progresses) -- to a
    model-and-search-mode-specific CSV spreadsheet in the repo root, one
    file per (model, recordTag) pair (record_<model><recordTag>.csv, with
    any non-alphanumeric character in the model name, e.g. GTR+FO's '+',
    replaced by '_' so the filename is always valid; recordTag from
    buildRecordTag, e.g. "_fast_reopt" or "_investigate3"). Writes a header
    row first if the file is new or empty. Deliberately APPENDS rather
    than overwriting (unlike output.txt's own plain-overwrite convention):
    every run tagged with its own runId, so repeated runs using the same
    search mode accumulate in one file, side by side, for comparing their
    convergence trajectories against each other later (e.g. in a
    spreadsheet application, filtered/pivoted by run_id).

    In a "_findopt"-tagged file specifically, true_minus_current does NOT
    mean the same thing on every row: the search's own rows (regular
    per-step/periodic/final writes) compare logL against trueTreeLogl fit
    under the main run's own model (plain, unfit JC unless useGtrModel),
    while findopt's own rows compare against trueTreeLoglForFindopt, a
    separately-fit GTR+FO (BRLEN_FIX) reference -- since findopt's scratch
    refit always runs under GTR+FO regardless of useGtrModel (see
    maybeRunFindopt's comment). This keeps each row's own diff meaningful
    for its own logL's model, but means the column isn't directly
    comparable row-to-row across the two kinds of writes without knowing
    which produced which.
 */
string recordSpreadsheetPath(const string &modelName, const string &recordTag) {
    string sanitizedModel = modelName;
    for (char &c : sanitizedModel)
        if (!isalnum((unsigned char) c))
            c = '_';
    return "record_" + sanitizedModel + recordTag + ".csv";
}

/**
    companion path to recordSpreadsheetPath, for appendRecordRow's own
    per-row topology dump (see its comment) -- same (model, recordTag)
    naming scheme, "topology_" prefix and ".nwk" extension instead of
    "record_"/".csv".
 */
string topologySpreadsheetPath(const string &modelName, const string &recordTag) {
    string sanitizedModel = modelName;
    for (char &c : sanitizedModel)
        if (!isalnum((unsigned char) c))
            c = '_';
    return "topology_" + sanitizedModel + recordTag + ".nwk";
}

void appendRecordRow(const string &modelName, const string &recordTag, const string &runId,
        long candidatesEvaluated, double timeElapsedSec, double logL, double trueTreeLogl,
        bool recordTopology, PhyloTree &treeForTopology) {
    string path = recordSpreadsheetPath(modelName, recordTag);

    ifstream check(path.c_str());
    bool needsHeader = !check.good() || check.peek() == ifstream::traits_type::eof();
    check.close();

    ofstream out(path.c_str(), ios::app);
    if (needsHeader)
        out << "run_id,candidates,time_elapsed,logL,true_minus_current" << endl;
    out << runId << "," << candidatesEvaluated << "," << timeElapsedSec << ","
        << setprecision(12) << logL << "," << (trueTreeLogl - logL) << endl;

    if (!recordTopology)
        return;

    // Companion file, only written when the separate "recordtopology" flag
    // is ALSO on (plain "record" alone writes only the CSV above, exactly
    // as it did before this flag existed -- the Newick dump costs real
    // time on large trees/long runs and most "record" uses have no need
    // for it): this row's tree TOPOLOGY (no branch lengths -- downstream LP
    // branch-length fitting refits those from scratch against a distance
    // matrix, so shipping applySPR's placeholder/reoptimized lengths here
    // would be pointless). One Newick line per CSV row, appended right here
    // so the two files stay positionally aligned (line N of this file <->
    // row N of the CSV, both counted after the CSV's own header) without
    // needing a shared join key -- PROVIDED "recordtopology" was given on
    // every run that contributed to this CSV; mixing runs with and without
    // it breaks that row<->line correspondence, since a without-it run
    // still appends CSV rows but no topology lines.
    stringstream topologyLine;
    treeForTopology.printTree(topologyLine, WT_SORT_TAXA);
    ofstream topologyOut(topologySpreadsheetPath(modelName, recordTag).c_str(), ios::app);
    topologyOut << topologyLine.str() << endl;
}

/**
    path for "trajectory"'s own topology dump -- one file per runHillClimb
    call, named off runId (already filename-safe: buildRunId only ever
    emits alphanumerics, '_', and '-') rather than off (modelName,
    recordTag) the way topologySpreadsheetPath is. recordTag deliberately
    groups MULTIPLE runs of the same search mode into one shared file (see
    buildRecordTag's comment) -- fine for topologySpreadsheetPath, since
    its lines stay positionally aligned with record's own per-run-tagged
    CSV rows, but "trajectory" has no such CSV to align against and no
    per-line run marker of its own, so sharing one file across runs would
    make it impossible to tell, after the fact, where one run's trajectory
    ends and the next one's begins. Keying off runId's own per-run
    timestamp instead keeps every run's trajectory in its own file.
 */
string trajectoryTopologyPath(const string &runId) {
    return "trajectory_" + runId + ".nwk";
}

/**
    "trajectory": append the current tree's topology (no branch lengths,
    same WT_SORT_TAXA convention as appendRecordRow's own topology dump --
    see its comment for why branch lengths are pointless here) as one more
    Newick line to this run's own trajectory_<run-id>.nwk. Called exactly
    twice per event worth recording: once right after the starting tree is
    built (before the step loop, capturing the post-initial-tree topology),
    and once more after every ACCEPTED step (the `if (improved)` branch in
    runHillClimb's main loop) -- never on a reverted step, and never from
    maybeRunPeriodicFullReopt/maybeRunFindopt, since neither of those ever
    changes the topology, only branch lengths and/or model parameters.
    Independent of "record"/"recordtopology": unlike recordTopology, which
    only ever adds a companion column to record's own CSV rows, this flag
    stands on its own and writes even when "record" was never given.
 */
void appendTrajectoryTopology(const string &runId, PhyloTree &tree) {
    stringstream topologyLine;
    tree.printTree(topologyLine, WT_SORT_TAXA);
    ofstream out(trajectoryTopologyPath(runId).c_str(), ios::app);
    out << topologyLine.str() << endl;
}

/**
    if fullReoptEveryNSteps is set and this SUCCESSFUL-step count is due, run one
    full-tree ML refit on the CURRENT tree -- branch lengths only
    (tree.optimizeAllBranches(fullReoptRounds)), or jointly with the
    model's own rate/frequency parameters
    (tree.getModelFactory()->optimizeParameters()) under useGtrModel --
    updating curScore and, if recordProgress is on, appending a row for
    it. fullReoptRounds (the "M" in "M rounds every N steps") only affects
    the non-gtr path: optimizeAllBranches takes an explicit sweep-count
    ceiling, but ModelFactory::optimizeParameters has no equivalent
    parameter of its own, converging by epsilon instead regardless. See
    fullReoptEveryNSteps' and useGtrModel's comments on runHillClimb for
    why. Called only right after an ACCEPTED (improved) step, once
    successfulSteps -- the running count of accepted moves so far,
    including this one -- is itself an exact multiple of
    fullReoptEveryNSteps; a stretch of rejected candidates between two
    accepted moves never counts toward that total, so this fires after
    every Nth SUCCESSFUL move, not every Nth step attempted.
 */
void maybeRunPeriodicFullReopt(PhyloTree &tree, long successfulSteps, int fullReoptEveryNSteps, int fullReoptRounds,
        bool useGtrModel, bool quiet, bool recordProgress, bool recordTopology, const string &modelName,
        const string &recordTag, const string &runId, long candidatesEvaluated, double cpuClockStart,
        double trueTreeLogl, double &curScore) {
    if (!(fullReoptEveryNSteps > 0 && successfulSteps % fullReoptEveryNSteps == 0))
        return;

    clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
    // The clamp just MUTATED branch lengths, so every cached partial
    // likelihood computed from the old ones is stale. optimizeAllBranches
    // tracks its running tree_lh incrementally and cross-checks it against
    // a fresh computeLikelihood at the end, asserting the two agree to
    // within 1.0 -- starting it from a stale cache makes them disagree by
    // however much the clamp moved things. On a degenerate alignment (many
    // near-zero branches, e.g. SARS-CoV-2 where most sites are constant)
    // the clamp touches enough branches at once to blow straight past that
    // tolerance and abort the run.
    tree.clearAllPartialLH();
    ostringstream suppressedFullReoptOutput;
    streambuf *realCoutBufFullReopt = cout.rdbuf(suppressedFullReoptOutput.rdbuf());
    double fullScore = useGtrModel
        ? tree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, Params::getInstance().modelEps)
        : tree.optimizeAllBranches(fullReoptRounds);
    cout.rdbuf(realCoutBufFullReopt);
    if (std::isfinite(fullScore))
        curScore = fullScore;
    if (!quiet)
        cout << "         (periodic full re-optimization: logL -> " << curScore << ")" << endl;
    if (recordProgress)
        appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getCPUTime() - cpuClockStart, curScore,
                trueTreeLogl, recordTopology, tree);
}

/**
    if findoptEveryNSteps is set (see runHillClimb's comment for how its
    default, "the total number of steps", is resolved) and this step index
    is due, run one full-tree ML refit -- exactly the same
    tree.optimizeAllBranches() / ModelFactory::optimizeParameters() call
    maybeRunPeriodicFullReopt makes -- but PURELY as a diagnostic: unlike
    that function (or the old, single-shot "finalreopt" this replaces),
    findopt never keeps its result, and never even touches the real tree in
    the first place. It builds a throwaway scratch copy of the CURRENT tree
    (initClonedTree, re-parsed from its own Newick text -- the same
    approach --branchlength-compare's B/C/D trees use, sharing the same
    Alignment pointer rather than copying it), runs the refit on THAT copy,
    reads off the logL it reaches, and lets the copy fall out of scope --
    PhyloTree's own destructor frees the model/model_factory/site_rate it
    built along the way (never the shared alignment; ownership of aln stays
    with the caller, exactly like trueTree's own sharing of it elsewhere in
    runHillClimb). The real tree, its model, and curScore are never mutated
    at all, so the search that continues after this call is completely
    unaffected by it -- curScore itself is passed by value, not reference,
    read only to print alongside findopt's own result for comparison.

    Because building and refitting the scratch copy is real work with no
    counterpart in the actual search's own cost budget, its own CPU-clock
    cost is excluded from the run's timing entirely ("pausing the timer"):
    t0 is captured before any of it starts, the elapsed search time up to
    (but not including) this call is what gets recorded for findopt's own
    CSV row, and cpuClockStart itself is pushed forward by however long the
    whole thing took -- every timing figure computed afterward (later rows,
    and the final "time elapsed" summary) is therefore exactly as if this
    call had taken zero time.

    The scratch tree's own Params::lk_safe_scaling is forced on for the
    duration of its clone+refit, regardless of the main run's own
    reoptimizeBranchLengths/fullReopt* setting -- a real branch-length and
    model-parameter search from possibly-naive lengths is exactly the
    scenario the plain kernel can hit a FATAL numerical underflow on (an
    outError, not the recoverable "WARNING: Numerical underflow" some other
    code paths print -- see scoreTrialSPRMove's comment for the general
    reasoning), and unlike reoptimizeBranchLengths findopt has no other flag
    guaranteeing it's already been turned on. It's restored to whatever it
    was immediately after, since it's read from the shared Params singleton
    dynamically rather than fixed onto the main tree at some earlier
    setLikelihoodKernel() call, so this window can't leak into the main
    tree's own likelihood evaluations at all.

    The scratch tree always fits under the richest available model for the
    alignment's (auto-detected) sequence type -- GTR+FO for DNA, LG+FO for
    protein, see modelNameFor's comment -- regardless of whether the main
    search itself is running under useGtrModel (plain JC/LG) or not:
    findopt is meant to answer "how good could this topology's branch
    lengths AND model get under the richest model available", not "what
    would this step's own model reach", so it always builds its clone with
    modelNameFor(aln->seq_type, true) and always refits via
    ModelFactory::optimizeParameters(), never plain optimizeAllBranches().

    trueTreeLogl here is NOT necessarily the same value runHillClimb prints
    / records for curScore's own comparisons -- the caller passes
    trueTreeLoglForFindopt, a separate GTR+FO-under-BRLEN_FIX fit of the
    true simulated tree built specifically so it stays model-matched with
    this function's own always-GTR+FO reading (see runHillClimb's own
    comment on trueTreeLoglForFindopt for why the plain trueTreeLogl would
    be an unfair reference whenever the main run isn't already useGtrModel).
 */
void maybeRunFindopt(PhyloTree &tree, int step, int findoptEveryNSteps, bool quiet,
        bool recordProgress, bool recordTopology, const string &modelName, const string &recordTag,
        const string &runId, long candidatesEvaluated, double &cpuClockStart, double trueTreeLogl, double curScore,
        Alignment *aln, Params &params, const string &refitModelName) {
    if (!(findoptEveryNSteps > 0 && (step + 1) % findoptEveryNSteps == 0))
        return;

    double t0 = getCPUTime();
    double searchTimeSoFar = t0 - cpuClockStart;

    // the scratch refit below is a real Newton-Raphson branch-length AND
    // model-parameter search that can start from whatever naive/
    // unreoptimized lengths the main search left the tree at -- exactly the
    // scenario the plain (non-scaled) kernel isn't built to handle without
    // a fatal numerical underflow (see reoptimizeSPREdges' comment on
    // runHillClimb for why reoptimizeBranchLengths needs the same thing).
    // Force it on for this scratch tree's own setup, regardless of whether
    // the main run's own reoptimizeBranchLengths/fullReopt* already enabled
    // it, then restore whatever it was immediately after -- lk_safe_scaling
    // is read off the shared Params singleton dynamically, not fixed for
    // the main tree at its own long-since-past setLikelihoodKernel() call,
    // so this window can't affect the main tree at all.
    bool origSafeScaling = params.lk_safe_scaling;
    params.lk_safe_scaling = true;
    ostringstream suppressedFindoptOutput;
    streambuf *realCoutBufFindopt = cout.rdbuf(suppressedFindoptOutput.rdbuf());

    // newickOf(tree) is NOT usable here -- it deliberately rounds every
    // branch length to 1 decimal place for readable terminal output (see
    // its own comment), which would start the scratch tree from a badly
    // corrupted approximation of the real tree instead of the real tree
    // itself. Print a FULL-precision Newick by hand instead, exactly like
    // runBranchLengthCompare's own startNewick does for its B/C/D clones.
    int savedPrecision = params.numeric_precision;
    params.numeric_precision = 15;
    ostringstream fullPrecisionNewick;
    tree.printTree(fullPrecisionNewick, WT_BR_LEN);
    params.numeric_precision = savedPrecision;

    // always the richest available model for the scratch refit (GTR+FO for
    // DNA, LG+FO for protein -- see modelNameFor's comment), independent of
    // modelName (which reflects the MAIN search's own model, plain JC/LG
    // unless useGtrModel is set)
    PhyloTree scratchTree;
    // `refitModelName` empty keeps this tool's own long-standing choice --
    // the richest model it knows for this sequence type, regardless of what
    // the main search is using. IQ-TREE passes its OWN model name instead,
    // since -m may well name something richer than GTR+FO/LG+FO and a
    // headroom figure measured under a weaker model would read as negative.
    initClonedTree(scratchTree, fullPrecisionNewick.str(), aln, params,
            refitModelName.empty() ? modelNameFor(aln->seq_type, true) : refitModelName);
    clampAllBranchLengthsForOptimization(scratchTree, Params::getInstance().min_branch_length);
    double findoptScore =
        scratchTree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, Params::getInstance().modelEps);
    cout.rdbuf(realCoutBufFindopt);
    params.lk_safe_scaling = origSafeScaling;

    if (!quiet)
        cout << "         (findopt: scratch whole-tree refit -> logL " << findoptScore
             << ", main tree unaffected, still at " << curScore << ")" << endl;
    if (recordProgress && std::isfinite(findoptScore))
        appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, searchTimeSoFar, findoptScore,
                trueTreeLogl, recordTopology, scratchTree);

    cpuClockStart += (getCPUTime() - t0); // pause the timer: this scratch pass never counts
}

/**
    "shrink": a stagnation counter, one option among several discussed for
    deciding WHEN to narrow the search radius from step-history data
    rather than a fixed schedule (see shrinkFlag's comment on runHillClimb
    for the fuller discussion and the other options that were considered
    but not implemented).

    Call once per step, right after that step's own accept/reject outcome
    (`improved`) is known. improved resets the stall counter to 0 (the search is still
    making progress at the CURRENT radius, no reason to narrow yet); a
    non-improving step increments it. Once the count reaches
    shrinkStallThreshold (the "N" the user gave, or its default), shrinkCurrentRadius
    drops by 1 (never below the floor of 1) and the counter resets, so it
    takes a fresh run of stalls before the next shrink.

    The threshold is held CONSTANT for the whole run, not scaled down as
    the step budget runs low -- an earlier version scaled it toward a
    floor of 1 near the end (reasoning: little budget left, so narrow
    eagerly), but that made the last stretch of a long run shrink through
    several radii in quick succession on just a stall or two each,
    skipping past radii that might still have had useful candidates left
    to find with a fair, full-length look. A constant threshold gives
    every radius the same chance to prove itself stalled, at any point in
    the run.

    "Skipped" steps (no legal candidate found at all) are deliberately
    NOT fed into this counter -- neither incrementing nor resetting it --
    since they never produced a genuine accept/reject outcome to begin
    with; counting them as stalls would conflate "no legal move existed"
    with "moves existed but didn't help", which call for different
    responses.
 */
void maybeShrinkRadius(bool improved, int shrinkStallThreshold, bool quiet,
        int &shrinkStallCount, int &shrinkCurrentRadius) {
    if (improved) {
        shrinkStallCount = 0;
        return;
    }
    shrinkStallCount++;
    if (shrinkStallCount >= shrinkStallThreshold && shrinkCurrentRadius > 1) {
        shrinkCurrentRadius--;
        if (!quiet)
            cout << "         (shrink: " << shrinkStallCount << " consecutive stall(s) -- radius reduced to "
                 << shrinkCurrentRadius << ")" << endl;
        shrinkStallCount = 0;
    }
}

/**
    standard-normal draw via Box-Muller, built on this tool's own
    random_double() (uniform on [0,1)) rather than <random>'s
    std::normal_distribution, purely so every source of randomness in this
    file draws from the same seeded stream (init_random, seeded once in
    runHillClimb) instead of mixing in a second, independently-seeded
    generator. Used only as sampleStandardGamma's own building block below.
 */
double sampleStandardNormal() {
    double u1 = random_double();
    while (u1 <= 0.0) // log(0) below would be -inf; redraw the rare zero
        u1 = random_double();
    double u2 = random_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * acos(-1.0) * u2);
}

/**
    draw from a standard Gamma(shape, 1) distribution via Marsaglia & Tsang's
    (2000) squeeze method -- the standard rejection-sampling algorithm for
    shape >= 1 (their "A Simple Method for Generating Gamma Variables"),
    with the well-known shape < 1 case handled by boosting to shape+1 and
    correcting with an independent U^(1/shape) draw (the same reduction
    <random>'s own std::gamma_distribution implementations use). learnRadius
    below scales this standard draw by its own fitted scale parameter to get
    an actual Gamma(shape, scale) sample.
 */
double sampleStandardGamma(double shape) {
    if (shape < 1.0) {
        double boosted = sampleStandardGamma(shape + 1.0);
        double u = random_double();
        while (u <= 0.0)
            u = random_double();
        return boosted * pow(u, 1.0 / shape);
    }
    double d = shape - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do {
            x = sampleStandardNormal();
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        double u = random_double();
        double x2 = x * x;
        if (u < 1.0 - 0.0331 * x2 * x2)
            return d * v;
        if (log(u) < 0.5 * x2 + d * (1.0 - v + log(v)))
            return d * v;
    }
}

/**
    "learnradius": records one more data point into the sliding window of
    the last (up to) windowSize successfully-accepted "normal" moves' own
    achieved radii, evicting the oldest once the window is full. Call ONLY
    for a move that (a) actually improved curScore and (b) was not an
    "investigate" refinement of a previous move -- see learnradiusFlag's
    comment on runHillClimb for why investigation steps are excluded.
    achievedRadius is deliberately a double, not an int, even though every
    current caller passes an integer hop count (bestDistance/walkLength):
    the window and the moment-fitting in learnRadius below are written to
    operate on a continuous scale throughout, so this same machinery keeps
    working unmodified if a future caller ever feeds it a genuinely
    fractional (e.g. substitution-distance) achieved radius instead.
 */
void recordLearnRadiusSample(deque<double> &window, int windowSize, double achievedRadius) {
    window.push_back(achievedRadius);
    while ((int) window.size() > windowSize)
        window.pop_front();
}

/**
    "learnradius"'s STARTING distribution -- what learnRadiusContinuous
    below draws from whenever there isn't (yet, or ever, on a losing
    weighted coin flip) enough history to trust the Gamma fit instead; see
    learnRadiusContinuous's own comment for how the two are mixed, and
    learnradiusFlag's comment on runHillClimb for the full picture and the
    reasoning behind this specific shape.

    A trapezoidal density on [0, maxPath]: FLAT (uniform) from 0 up to
    `center` -- the <radius> the user actually gave, reinterpreted, under
    "learnradius", as the MIDDLE of this starting spread rather than a hard
    ceiling -- then ramping down LINEARLY from `center` to `maxPath`,
    reaching a density of exactly 0 at maxPath itself. maxPath is never
    actually drawn (probability zero, not just unlikely): it's a
    structural limit (every edge the tree has, or 100% of its own branch
    length -- see this function's caller in runHillClimb for which,
    depending on "distradius"), not a value <radius> should realistically
    steer toward.

    The flat and ramp pieces are normalized as ONE continuous density (the
    ramp's own peak, right at `center`, sits at the SAME height the flat
    piece already has there -- no seam/discontinuity), making this a
    proper trapezoid rather than two independently-scaled pieces glued
    together. Derivation of the flat/ramp split probability: let h be that
    shared height at `center` (h is also the flat piece's own constant
    height everywhere on [0,center]). flat area = h*center (a
    center-by-h rectangle); ramp area = 0.5*h*(maxPath-center) (a right
    triangle, base maxPath-center, height h). Requiring both areas to sum
    to 1 (a valid density integrates to 1) gives h = 2/(center+maxPath),
    and therefore
        P(flat) = h*center = 2*center / (center+maxPath)
        P(ramp) = 1 - P(flat) = (maxPath-center) / (center+maxPath)
    Sampling: flip that weighted coin, then draw uniformly within the flat
    piece or via the standard closed-form inverse CDF for a left-mode/
    right-zero triangular distribution within the ramp piece (X = maxPath -
    (maxPath-center)*sqrt(U), U ~ Uniform(0,1) -- the textbook inverse CDF
    for a triangular distribution whose mode sits at its own left edge).

    Falls back to a plain Uniform(0, center) (skipping the ramp piece
    entirely) if maxPath <= center -- a degenerate/tiny-tree edge case
    where there is no room left for a ramp at all (e.g. <radius> already at
    or past the tree's own edge count).
 */
double sampleStartingRadius(double center, double maxPath) {
    if (center <= 0.0)
        return 0.0;
    if (maxPath <= center)
        return random_double() * center;

    double flatProb = (2.0 * center) / (center + maxPath);
    if (random_double() < flatProb)
        return random_double() * center;

    double u = random_double();
    return maxPath - (maxPath - center) * sqrt(u);
}

/**
    "learnradius": draw this step's search radius from a distribution
    fitted to the sliding window recordLearnRadiusSample has been filling in
    -- see learnradiusFlag's comment on runHillClimb for the full picture;
    this is where the mixture it describes is actually implemented.

    The window holds real-valued (not bucketed/rounded) observations
    throughout, and every statistic computed from it here -- sample mean,
    sample variance, the Gamma shape/scale method-of-moments fit, and
    sampleStandardGamma's own continuous draw -- stays in double precision
    all the way out to this function's own return value: nothing here
    rounds or bucket the values. This is what "continuous" means in
    learnradiusFlag's own comment: findGraftPositions/chooseGraft need an
    INT hop count (rounding happens in runHillClimb's own step loop, the
    only place that needs it), but chooseGraftByDistance's radiusPercent is
    itself a double, and the whole point of composing "learnradius" with
    "distradius" is to feed it a genuinely fractional percentage rather than
    one rounded down to a whole number first -- see this function's own
    caller in runHillClimb's step loop for exactly where that split happens.
    Nothing about the fitting itself assumes the window's values are hop
    counts specifically either: the same code works whether they happen to
    be small integers (plain hop-count radius) or fine-grained percentages
    (distradius), without caring which -- maybeFinalizeLearnRadiusExcursion
    and the "fast"+"distradius" branch in runHillClimb's step loop are what
    decide which kind actually goes into the window to begin with.

    Distribution: a two-component mixture,
        gammaWeight * Gamma(shape, scale) + (1 - gammaWeight) * sampleStartingRadius(center, maxPath)
    where gammaWeight = min(kMaxGammaWeight, window.size() / windowSize) -- mostly
    the "proportional part of the distribution" learnradiusFlag's comment
    describes: 0 with an empty window (pure starting distribution, matching
    "initially... uniform between 1 and the max radius" -- see
    sampleStartingRadius's own comment for why that's now a trapezoid
    centered on <radius>, not a flat Uniform(1,radius)), growing linearly
    as the window fills -- EXCEPT it now caps out at kMaxGammaWeight rather
    than reaching a full 1.0: even once windowSize successful moves have
    accumulated, a small kMaxGammaWeight-complement chance of a fresh
    sampleStartingRadius draw always remains, on every single call, for as
    long as the run continues. Each call independently flips a weighted
    coin between the two components rather than blending their outputs,
    which is what makes the RESULT follow that mixture distribution
    (blending two samples' VALUES would not).

    shape/scale come from the standard method-of-moments Gamma fit (mean^2/
    variance, variance/mean) to the window's current contents, but shape is
    additionally capped at kMaxGammaShape (scale recomputed as mean/shape
    afterward, to keep the fitted MEAN unchanged even though the capped
    shape no longer matches the window's own raw sample variance exactly).
    A Gamma's own coefficient of variation is 1/sqrt(shape), so this is a
    floor on how tight the fitted spread can ever get -- shape would
    otherwise grow without bound (and the fit collapse toward a near-point
    mass at the sample mean) if the window ever filled with nearly-identical
    achieved radii, which is exactly the outcome a search that has honed in
    on one effective radius would naturally tend to keep producing on its
    own: every additional near-identical sample would shrink the sample
    variance further, which would shrink the fitted spread further, which
    would make the NEXT drawn radius even more likely to land close enough
    to keep reinforcing the same narrow cluster -- a self-tightening
    feedback loop with nothing in it to ever loosen back up, left
    unchecked. Capping shape breaks that loop directly: no matter how
    tightly clustered the window's own data gets, every draw keeps a real,
    bounded-below chance of landing meaningfully away from the current
    mean, so if a genuinely different radius starts paying off, its own
    achieved values can still enter the window and pull the fitted mean
    toward it, rather than the fit having already locked itself out of
    ever producing a large enough excursion to discover that in the first
    place. When variance is unusable (fewer than 2 samples -- gated by this
    function's own caller below, kept here as an explicit guard against
    division by zero rather than relying on that alone -- or the raw
    sample variance is ~0), shape falls back to kMaxGammaShape outright
    instead of blowing up toward infinity or collapsing to a bare
    `drawn = mean` point value the way this function's own earlier design
    did.
    kMaxGammaWeight and kMaxGammaShape are both EXPERIMENTAL -- picked as
    round, defensible starting points (10% minimum ongoing exploration;
    roughly 22% minimum coefficient of variation), not empirically tuned.

    Unlike this function's own earlier design, the Gamma-fit component is
    NOT clamped down to `center` (<radius>) at all -- real accumulated
    history gets a genuine chance to push the learned radius past it, all
    the way out to `maxPath`, the same structural ceiling
    sampleStartingRadius's own ramp piece reaches zero density at. <radius>
    only shapes the STARTING guess now (sampleStartingRadius's own
    "middle"); it is never a hard cap on what the search can actually learn
    to use, since the accumulated window can legitimately show that a
    larger radius keeps paying off.
 */
double learnRadiusContinuous(const deque<double> &window, int windowSize, double center, double maxPath) {
    const double kMaxGammaWeight = 0.9;
    const double kMaxGammaShape = 20.0;

    double gammaWeight = (windowSize > 0)
        ? min(kMaxGammaWeight, ((double) window.size()) / (double) windowSize) : 0.0;
    double drawn;
    if (window.size() >= 2 && random_double() < gammaWeight) {
        double mean = 0.0;
        for (double v : window)
            mean += v;
        mean /= (double) window.size();
        double variance = 0.0;
        for (double v : window)
            variance += (v - mean) * (v - mean);
        variance /= (double) (window.size() - 1);
        double shape = kMaxGammaShape;
        if (variance > 1e-9 && mean > 0.0)
            shape = min(kMaxGammaShape, (mean * mean) / variance);
        double scale = (mean > 0.0) ? (mean / shape) : 0.0;
        drawn = (mean > 0.0) ? (scale * sampleStandardGamma(shape)) : sampleStartingRadius(center, maxPath);
    } else {
        drawn = sampleStartingRadius(center, maxPath);
    }
    if (drawn < 1.0)
        drawn = 1.0;
    if (drawn > maxPath)
        drawn = maxPath;
    return drawn;
}

/**
    plain, unrestricted hop distance between two REAL, already-existing
    edges of the CURRENT tree, in findGraftPositions' own "radius"
    convention (see its comment: radius 1 is an edge incident to either
    endpoint of the seed edge, an NNI-equivalent distance; radius 2 one hop
    further, and so on). Unlike findGraftPositions, this performs no
    pruneNode-subtree exclusion and no isLegalSPR filtering -- both edges
    given are assumed to already be real tree edges (not a hypothetical
    prune/regraft pair), so a plain, unrestricted BFS is all that's needed,
    seeded from the (seedA,seedB) edge exactly the way findGraftPositions
    seeds from a collapsed sibling1/sibling2 view (see its comment) --
    which is what makes radius 1 mean the same thing in both places.

    Used only by learnradiusFlag's own "investigate"-excursion bookkeeping
    (maybeFinalizeLearnRadiusExcursion below, and learnradiusFlag's comment
    on runHillClimb) to measure how far an entire chain of investigate
    refinements ended up from where it started, once the whole chain is
    done -- in the SAME units findGraftPositions/chooseGraft already use
    elsewhere in this file, so the result is directly comparable to a
    plain (non-chained) move's own radius.

    outLength, if non-null, receives the actual summed branch length along
    that same path (not just its hop count) -- there is exactly one path
    between any two edges of a tree, so the hop-count BFS below is already
    tracing it; this just also accumulates each hop's own edge length along
    the way. Used by maybeFinalizeLearnRadiusExcursion to convert an
    excursion's displacement into a radiusPercent-equivalent when composed
    with "distradius", the same way chooseGraftByDistance's own
    outPercentUsed does for a single, non-chained move.

    @return the hop distance, or -1 if the target edge is never reached
    (should not happen for two edges of the same connected tree; guarded
    rather than risking an infinite loop on any unexpected disconnection)
 */
int edgeHopDistance(PhyloNode *seedA, PhyloNode *seedB, PhyloNode *targetA, PhyloNode *targetB,
        double *outLength) {
    pair<int,int> targetKey(min(targetA->id, targetB->id), max(targetA->id, targetB->id));
    pair<int,int> seedKey(min(seedA->id, seedB->id), max(seedA->id, seedB->id));
    if (targetKey == seedKey) {
        if (outLength)
            *outLength = 0.0;
        return 0;
    }

    set<pair<int,int> > seenEdges;
    seenEdges.insert(seedKey);

    struct QueueItem {
        PhyloNode *node;
        PhyloNode *cameFrom;
        int dist;
        double lengthSoFar;
    };
    queue<QueueItem> q;
    q.push({seedA, seedB, 0, 0.0});
    q.push({seedB, seedA, 0, 0.0});

    while (!q.empty()) {
        QueueItem cur = q.front();
        q.pop();

        FOR_NEIGHBOR_IT(cur.node, cur.cameFrom, it) {
            PhyloNode *next = (PhyloNode*) (*it)->node;
            int edgeRadius = cur.dist + 1;
            double edgeLength = cur.lengthSoFar + (*it)->length;

            pair<int,int> key(min(cur.node->id, next->id), max(cur.node->id, next->id));
            if (key == targetKey) {
                if (outLength)
                    *outLength = edgeLength;
                return edgeRadius;
            }
            if (seenEdges.insert(key).second)
                q.push({next, cur.node, edgeRadius, edgeLength});
        }
    }
    return -1;
}

/**
    "learnradius" + "investigate" together (see learnradiusFlag's comment
    on runHillClimb): close out the CURRENTLY OPEN excursion, if there is
    one, by measuring its net start-to-finish displacement and feeding that
    single number into the learnradius window -- then mark the excursion
    closed. A no-op (via excursionOpen's own guard) whenever there is
    nothing open to close, so every caller below can call this
    unconditionally on every path where an investigation attempt has just
    ended (whether it failed to improve, after its own rollback, or found
    no legal candidate at all) without first re-deriving whether an
    excursion happens to be open.

    currentDad is the CURRENT position of the excursion's own prune_dad --
    i.e. pruneDad/investigatePruneDad as the caller's own local scope has
    it at the moment of the call, already reflecting any rollback that
    needed to happen first. Paired with excursionFinalNode (the last
    successful step's own regraft_node, remembered across the whole
    excursion), this reconstructs exactly one half of that last step's own
    (now again real, since nothing has touched it since) graft edge -- see
    the comment where excursionFinalNode is first set, on runHillClimb's
    own step loop, for why using pruneDad's own THIRD edge (to pruneNode
    itself) here instead would be off by one hop.

    useDistanceRadius (the "distradius" flag): investigation refinements
    are always hop-based regardless of it (investigatingThisStep always
    goes through findGraftPositions/investigateRadius -- see
    investigateFlag's comment on runHillClimb), so edgeHopDistance's own
    hop count is what an excursion is made of either way. But the value
    THIS function feeds back into the window has to match whatever units
    the window's OTHER entries are in, so that a future
    learnRadius(Continuous) draw from it means the same thing every time
    it's used: with "distradius" also on, edgeHopDistance's outLength (the
    actual summed branch length of the excursion's own path) is converted
    to a radiusPercent-equivalent (100 * length / tree.treeLength()) the
    same way chooseGraftByDistance's own outPercentUsed converts a single,
    non-chained move's spent budget; without it, the plain hop count is
    recorded as-is, exactly as before.
 */
void maybeFinalizeLearnRadiusExcursion(PhyloTree &tree, bool useDistanceRadius, bool &excursionOpen,
        PhyloNode *anchorA, PhyloNode *anchorB, PhyloNode *excursionFinalNode, PhyloNode *currentDad,
        deque<double> &window, int windowSize) {
    if (!excursionOpen)
        return;
    double length = 0.0;
    int dist = edgeHopDistance(anchorA, anchorB, excursionFinalNode, currentDad, &length);
    if (dist >= 1) {
        double achieved = dist;
        if (useDistanceRadius) {
            double treeLen = tree.treeLength();
            achieved = (treeLen > 0.0) ? (100.0 * length / treeLen) : 0.0;
        }
        recordLearnRadiusSample(window, windowSize, achieved);
    }
    excursionOpen = false;
}

/**
    "sweep" (see sweepFlag's comment on runHillClimb for the full picture):
    score how well internal edge (p,q) -- i.e. p and q are each other's
    neighbor and both are internal (degree-3) nodes, so each has exactly
    two OTHER neighbors -- supports its OWN current grouping, versus the
    two alternative NNI rearrangements around that same edge. This is
    exactly adjacent_subtree_compatibility.pdf's section 8 statistic
    ("when compatibility means 'should be siblings'"), S_AB, chosen over
    the note's more general per-PAIR statistic C_bar (section 4) because
    C_bar requires isolating a subtree's OWN standalone fit from the rest
    of the tree, which -- properly done, conditioned on everything outside
    the pair per equation 13 -- has no safe, public IQ-TREE API to compute
    without reaching into internal, kernel-specific partial-likelihood
    buffer layouts (computePartialLikelihood's only exposed overload takes
    a TraversalInfo& and a thread id, plainly meant for the library's own
    batched internals, not a single ad-hoc external call). S_AB, by
    contrast, is expressible purely as REAL, whole-tree computeLikelihood()
    evaluations of REAL, already-legal topologies -- exactly the same kind
    of trial this tool already performs everywhere else (scoreTrialSPRMove
    itself) -- so it needs nothing beyond the existing, already-proven
    apply/score/rollback machinery.

    Let p's two other neighbors be A,B and q's two other neighbors be C,D
    (findSPRSiblings gives each pair). The tree's CURRENT topology, at this
    edge, is the quartet resolution "AB|CD" (A,B together on p's side, C,D
    together on q's side) -- exactly the "siblings" whose compatibility is
    in question. The two single-NNI alternatives swap one of {A,B} with
    one of {C,D}: "AC|BD" (realized here as the SPR move prune B from p,
    regraft onto edge q-C) and "AD|BC" (prune B from p, regraft onto edge
    q-D) -- picking A instead of B as the one pruned would just relabel
    which alternative is "first", not change the SET of two alternatives
    considered, so the choice is arbitrary.

    S = curScore - logmeanexp(ell_AC|BD, ell_AD|BC), matching equation (19)
    exactly (curScore standing in for ell_AB|CD, the tree's own current,
    already-known likelihood -- no need to recompute it). A LOW (very
    negative) S means the current AB|CD grouping is poorly supported
    relative to the alternatives -- these two "siblings" are a bad fit --
    which is exactly the ranking sweepFlag uses to pick its N worst
    positions. Each alternative's likelihood is scored via the ordinary
    scoreTrialSPRMove (apply, evaluate, always roll back) with branch-length
    reoptimization deliberately OFF regardless of this run's own
    reoptimizeBranchLengths setting -- this ranking pass is a cheap,
    approximate screen to choose WHICH N positions get sweepFlag's
    expensive, exhaustive, whole-tree treatment (which DOES honor
    reoptimizeBranchLengths); it doesn't need its own answer to be as
    precise as that.

    Also reports (outPruneNode, outPruneDad) = (siblingB, p) -- the exact
    prune point the two scored NNI alternatives above were built from --
    for the caller to reuse AS THE prune position for sweepFlag's own
    follow-up exhaustive whole-tree search. This matters: an earlier
    version of this function only reported the score, and the caller
    re-derived a prune point from the edge (p,q) itself via
    orientPruneEdgeForSweep -- which, since both p and q are always
    degree-3 and the tree is never rooted here, deterministically prunes
    q's WHOLE subtree away from p instead. That is a different, much
    coarser move than the one actually scored (which only prunes siblingB,
    a single one of p's OTHER neighbors) -- so the resulting search never
    even considered the two specific candidates (regraft onto edge q-C or
    q-D) that made this edge's score so low, and reliably found nothing to
    keep. Reporting the SAME (siblingB, p) pair used for scoring guarantees
    the follow-up findGraftPositions call, seeded at p's other neighbors
    (including q), enumerates both of those exact candidates again -- so
    the exhaustive search can never do WORSE than the better of the two
    already-scored alternatives.

    @return false if this edge has no quartet to test at all (p or q is a
    leaf) or either NNI alternative turns out illegal (only possible near
    the tree's arbitrary root leaf; --hillclimb's own trees are never
    rooted, so in practice this only happens if p or q itself IS the root
    leaf, which can't occur here since edgeRegistry excludes the root's one
    edge, but neighboring edges could still involve it as one of A/B/C/D --
    not as p or q). When false, `outScore`/`outPruneNode`/`outPruneDad` are
    left unmodified and this edge is simply excluded from the ranking, the
    same way choosePrune/chooseGraft skip whatever isn't legal rather than
    treating it as an error.
 */
bool computeSiblingCompatibilityScore(PhyloTree &tree, PhyloNode *p, PhyloNode *q, double curScore, double &outScore,
        PhyloNode *&outPruneNode, PhyloNode *&outPruneDad) {
    if (p->degree() != 3 || q->degree() != 3)
        return false;

    PhyloNode *siblingA, *siblingB, *siblingC, *siblingD;
    findSPRSiblings(q, p, siblingA, siblingB); // p's two other neighbors
    findSPRSiblings(p, q, siblingC, siblingD); // q's two other neighbors

    SPRMove moveAC;
    moveAC.prune_node = siblingB;
    moveAC.prune_dad = p;
    moveAC.regraft_node = siblingC;
    moveAC.regraft_dad = q;
    moveAC.radius = 1;
    moveAC.screening_score = 0.0;
    moveAC.exact_score = 0.0;
    moveAC.candidate_id = 0;
    moveAC.generation = -1;

    SPRMove moveAD = moveAC;
    moveAD.regraft_node = siblingD;

    if (!tree.isLegalSPR(moveAC) || !tree.isLegalSPR(moveAD))
        return false;

    double ellAC = scoreTrialSPRMove(tree, moveAC, false);
    double ellAD = scoreTrialSPRMove(tree, moveAD, false);

    // logmeanexp(ellAC, ellAD) = log( (exp(ellAC) + exp(ellAD)) / 2 ),
    // shifted by the larger of the two first (the standard log-sum-exp
    // stabilization) so neither exp() call can overflow
    double hi = max(ellAC, ellAD);
    double logMeanExp = hi + log(0.5 * (exp(ellAC - hi) + exp(ellAD - hi)));

    outScore = curScore - logMeanExp;
    outPruneNode = siblingB;
    outPruneDad = p;
    return true;
}

// A second attempt at speeding up per-candidate trial scoring (beyond the
// standalone-tree subtree-likelihood attempt below) was tried and reverted:
// instead of a full resetLikelihoodBuffers() (delete+reinit of every
// partial-likelihood buffer) around each trial apply/rollback, clear only
// the 7 PhyloNeighbor directions an SPR move can possibly touch (the merged
// edge left behind at the vacated attachment point, plus the two split
// halves of the target edge at the new attachment point) via the PUBLIC
// clearPartialLh(), then evaluate computeLikelihoodBranch() anchored right
// at the graft point instead of at the tree's fixed root -- mathematically
// identical for a reversible, stationary model, and exactly the technique
// IQ-TREE's own NNI search already uses safely (see getBestNNIForBran,
// phylotree.cpp). Unlike the subtree-likelihood attempt, this never
// fabricated any data or bypassed buffer bookkeeping -- it only marked real
// data as "needs recomputing", the same operation NNI performs.
// It still crashed (phylotree.cpp's reorientPartialLh assertion), for a
// related but distinct reason: IQ-TREE's default memory mode (LM_PER_NODE)
// gives each NODE a small shared pool of buffers, not one per possible
// neighbor-direction, on the assumption that at most one direction is
// "away and uncomputed" at a time. An SPR move's dad1 sits at the junction
// of two of the three touched edges (dad1-dad2 and dad1-node2), so clearing
// both simultaneously asks dad1 for two uncomputed buffers at once --
// apparently more than its pool provides, with no sibling buffer available
// to steal via reorientPartialLh's usual takeover logic. NNI never hits
// this because it only ever touches one swapped branch (or processes
// adjacent branches strictly one at a time). Making this safe would need
// either a strictly sequential, NNI-style one-direction-at-a-time
// restructuring, or a deeper look at IQ-TREE's mem_slots sizing -- judged
// not worth pursuing further after a second confirmed crash in this same
// class of internal buffer-management constraint.

// A faster, direct-on-the-big-tree alternative to the standalone-tree
// subtree likelihood below was attempted and reverted. The idea: for a
// reversible, stationary model, freq . P(t) = freq for any branch length t
// (the stationary distribution is a fixed point of the transition matrix),
// so substituting the model's own state frequencies for "everything beyond
// dad" on one side of the (lca, lcaDad) branch, then calling the tree's
// ordinary computeLikelihoodBranch(), should mathematically collapse to
// exactly the local clade's own likelihood -- verified correct on the
// first call (a sane, correctly-scaled logL). It crashed on a second call
// to the same node, however: IQ-TREE's default memory mode (LM_PER_NODE,
// see PhyloTree::reorientPartialLh in phylotree.cpp) treats partial-
// likelihood buffers as a SHARED POOL PER NODE, reassigned ("reoriented")
// between neighbor-directions on demand, not one dedicated buffer per
// direction. Permanently pinning one direction's buffer with substituted
// content (as this approach requires) starves that pool, and a later,
// legitimate computation elsewhere at the same node can find no buffer
// left to reorient into, triggering an assertion failure. Making this
// safe would require also driving IQ-TREE's `mem_slots` lock/unlock/
// takeover bookkeeping to protect the substituted buffer from reorientation
// -- a further layer of undocumented internals with the same risk profile,
// judged not worth pursuing further here.


/*==========================================================================
    Re-entrant driver -- see search_experiments/sprsearch.h for the design.
 *========================================================================*/

AcceptDist::AcceptDist()
        : enabled(false), temperature(0.5), shape(1.0), anneal(false), tempFloor(0.0) {
}

double AcceptDist::effectiveTemperature(double progress) const {
    if (!anneal)
        return temperature;
    if (progress < 0.0) progress = 0.0;
    if (progress > 1.0) progress = 1.0;
    double t = temperature * (1.0 - progress);
    return (t < tempFloor) ? tempFloor : t;
}

double AcceptDist::probability(double delta, double progress) const {
    if (!enabled)
        return (delta > 0.0) ? 1.0 : 0.0;
    if (delta >= 0.0)
        return 1.0;                      // an improvement is never refused
    double t = effectiveTemperature(progress);
    if (t <= 0.0)
        return 0.0;                      // frozen: strict hill-climbing
    double x = -delta / t;               // > 0, since delta < 0
    if (shape != 1.0)
        x = pow(x, shape);
    if (x > 700.0)
        return 0.0;                      // exp() would underflow anyway
    return exp(-x);
}

bool AcceptDist::accept(double delta, double progress) const {
    // Improvements short-circuit WITHOUT drawing, so enabling this rule
    // does not perturb the random stream for moves that would have been
    // kept regardless -- an improve-only run stays reproducible.
    if (delta > 0.0)
        return true;
    if (!enabled)
        return false;
    double p = probability(delta, progress);
    if (p <= 0.0)
        return false;
    return random_double() < p;
}

SPRSearchOptions::SPRSearchOptions()
        : radius(6), useFastSelection(false), numCandidates(1), useDistanceRadius(false),
          weightpruneFlag(PRUNE_UNIFORM), alternateFlag(false), investigateFlag(false), investigateRadius(1),
          shrinkFlag(false), shrinkStallThreshold(10), learnradiusFlag(false), learnradiusN(20),
          reoptimizeBranchLengths(false), fullReoptEveryNSteps(0), fullReoptRounds(100),
          escapeFlag(false), escapeSpan(0), escapeTries(0),
          tunnelFlag(false), tunnelTolerance(0.0), tunnelTries(3),
          slackFlag(false), slackDelta(0.0), slackAnneal(false), acceptProgressBase(0.0),
          useGtrModel(false), sweepFlag(false), sweepCount(10), findoptFlag(false),
          findoptEveryNSteps(0), quiet(false), recordProgress(false), recordTopology(false),
          trajectoryFlag(false), noTrueTree(true), stepsPerPass(0) {
}

SPRSearchState::SPRSearchState()
        : cpuClockStart(0.0), trueTreeLogl(std::numeric_limits<double>::quiet_NaN()),
          trueTreeLoglForFindopt(std::numeric_limits<double>::quiet_NaN()),
          candidatesEvaluated(0), tunnelEntered(0), tunnelProbes(0), tunnelCommitted(0),
          slackAccepted(0), slackRejected(0), slackKicks(0), slackLastDelta(0.0),
          acceptedDownhill(0), offeredDownhill(0), lastTemperature(0.0),
          bestSeenScore(-DBL_MAX), bestSeenRestored(false), bestSeenRestores(0),
          successfulSteps(0), stepsRun(0),
          learnRadiusMaxPath(0.0), learnRadiusExcursionOpen(false), learnRadiusAnchorA(nullptr),
          learnRadiusAnchorB(nullptr), learnRadiusFinalNode(nullptr), investigateNext(false),
          investigatePruneNode(nullptr), investigatePruneDad(nullptr),
          shrinkStallCount(0), shrinkCurrentRadius(0) {
}

void initSPRSearchState(SPRSearchState &st, const SPRSearchOptions &opt, const EdgeRegistry &reg,
        const string &runId, const string &recordTag, const string &modelName,
        double cpuClockStart, double trueTreeLogl, double trueTreeLoglForFindopt) {
    st.runId = runId;
    st.recordTag = recordTag;
    st.modelName = modelName;
    st.cpuClockStart = cpuClockStart;
    st.trueTreeLogl = trueTreeLogl;
    st.trueTreeLoglForFindopt = trueTreeLoglForFindopt;
    st.shrinkCurrentRadius = opt.radius;
    // "learnradius"'s own structural ceiling, computed once rather than
    // per-step: an SPR move always destroys exactly 3 edges and creates 3
    // more, so the registry's edge count never changes across the run.
    // Under "distradius" <radius> is a PERCENTAGE, whose natural ceiling
    // is 100 (a walk with a 100%-of-tree-length budget can already reach
    // anywhere); the plain hop-count case's ceiling is the edge count,
    // since no simple path can use more hops than the tree has edges.
    st.learnRadiusMaxPath = opt.useDistanceRadius ? 100.0 : (double) reg.slots.size();
}

void beginSPRSearchPass(SPRSearchState &st) {
    // Every PhyloNode* below points into whatever tree object the PREVIOUS
    // pass ran on. IQ-TREE rebuilds its tree from a Newick string on every
    // perturbation, so those are dangling by the time the next pass starts
    // -- drop them and let this pass rebuild its own. Everything else in
    // the state (the learnradius window, the shrink counters, the running
    // totals, the record identity) is deliberately left alone: that is the
    // half that is supposed to accumulate across the whole run.
    st.learnRadiusExcursionOpen = false;
    st.learnRadiusAnchorA = nullptr;
    st.learnRadiusAnchorB = nullptr;
    st.learnRadiusFinalNode = nullptr;
    st.investigateNext = false;
    st.investigatePruneNode = nullptr;
    st.investigatePruneDad = nullptr;
}

double runSPRSteps(PhyloTree &tree, EdgeRegistry &reg, const SPRSearchOptions &opt, SPRSearchState &st,
        double curScore, int maxSteps, Alignment *aln, Params &params,
        SPRLocalLhCache *localCache, const string &stepLabelPrefix) {
    // --accept-dist: seed the high-water mark from where this call
    // starts, so a stage that only ever loses ground still reports the tree
    // it began with rather than the endpoint of a downhill walk.
    if (opt.acceptDist.enabled && curScore > st.bestSeenScore) {
        st.bestSeenScore = curScore;
        st.bestSeenTree = tree.getTreeString();
    }

    int step = 0;
    for (; step < maxSteps; step++) {
        // Annealing progress for the acceptance rule. Measured against the
        // whole run's budget (stepsPerPass), not this block's maxSteps: the
        // continuous stage calls this in blocks, and a per-block schedule
        // would restart the cooling curve every block.
        // Prefer this pass's own step budget when there is one; fall back
        // to whatever progress the caller reports otherwise, so an in-loop
        // refinement (which has no step budget) still cools.
        double acceptProgress = opt.acceptProgressBase;
        if (opt.acceptDist.enabled && opt.acceptDist.anneal && opt.stepsPerPass > 0) {
            acceptProgress = (double) (st.stepsRun + step) / (double) opt.stepsPerPass;
        }
        if (acceptProgress > 1.0)
            acceptProgress = 1.0;
        if (acceptProgress < 0.0)
            acceptProgress = 0.0;
        int stepRadius = opt.radius;
        // stepRadiusContinuous mirrors stepRadius, except it's allowed to
        // stay fractional -- only actually diverges from (double)stepRadius
        // when "learnradius" is composed with "distradius" (see
        // learnRadiusContinuous's own comment): chooseGraftByDistance's
        // radiusPercent wants the genuinely continuous draw, everything
        // else (findGraftPositions/chooseGraft's int radius, this step's
        // own printed label) wants stepRadius, rounded
        double stepRadiusContinuous = (double) opt.radius;
        if (opt.shrinkFlag)
            // "shrink" replaces the fixed <radius> with its own
            // stagnation-driven value
            stepRadius = st.shrinkCurrentRadius;
        else if (opt.learnradiusFlag) {
            // "learnradius" replaces the fixed <radius> with a fresh draw
            // from the Gamma/starting-distribution mixture fit to recent
            // successes -- mutually exclusive with "shrink"
            // (parseHillClimbFlags rejects giving both), so this is never
            // reached when shrinkFlag is on. <radius> is passed as the
            // MIDDLE of the starting distribution here, not a ceiling --
            // the draw can legitimately land anywhere up to
            // learnRadiusMaxPath once real history supports it; see
            // learnRadiusContinuous's own comment
            stepRadiusContinuous = learnRadiusContinuous(st.learnRadiusWindow, opt.learnradiusN, (double) opt.radius,
                    st.learnRadiusMaxPath);
            stepRadius = (int) llround(stepRadiusContinuous);
            if (stepRadius < 1)
                stepRadius = 1;
            if (stepRadius > (int) st.learnRadiusMaxPath)
                stepRadius = (int) st.learnRadiusMaxPath;
        }

        bool investigatingThisStep = opt.investigateFlag && st.investigateNext;
        // default to "not investigating next" -- only re-armed below, and
        // only if THIS step's own candidate is actually accepted; every
        // early skip/rejection path below falls through with this staying
        // false, so investigation naturally stops the first time a
        // refinement attempt doesn't pan out
        if (opt.investigateFlag)
            st.investigateNext = false;

        // "alternate": every other step forces a radius-1 (NNI-equivalent)
        // search instead of <radius>, regardless of which candidate-
        // selection path (fast or exhaustive) is actually active this
        // step -- investigate takes priority when it's active, since
        // that's testing a specific
        // already-accepted move at investigateRadius, not a step-parity
        // toggle. See alternateFlag's comment on runHillClimb.
        bool nniStepThisTime = opt.alternateFlag && !investigatingThisStep && (step % 2 != 0);
        int effectiveRadius = nniStepThisTime ? 1 : stepRadius;
        double effectiveRadiusContinuous = nniStepThisTime ? 1.0 : stepRadiusContinuous;

        PhyloNode *pruneNode, *pruneDad;
        if (investigatingThisStep) {
            pruneNode = st.investigatePruneNode;
            pruneDad = st.investigatePruneDad;
        } else if (!choosePrune(tree, reg, pruneNode, pruneDad, opt.weightpruneFlag)) {
            if (!opt.quiet)
                cout << stepLabelPrefix << "step " << (st.stepsRun + step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }

        // stepRadiusText: stepRadius's own printed form -- plain integer,
        // except when "learnradius" and "distradius" are both active, where
        // stepRadius (rounded) would hide the actual fractional percentage
        // chooseGraftByDistance is really about to use (effectiveRadiusContinuous)
        string stepRadiusText;
        if (opt.learnradiusFlag && opt.useDistanceRadius) {
            ostringstream fmt;
            fmt << fixed << setprecision(2) << stepRadiusContinuous;
            stepRadiusText = fmt.str();
        } else
            stepRadiusText = to_string(stepRadius);

        string stepLabel;
        if (investigatingThisStep)
            stepLabel = "(investigate)";
        else if (opt.alternateFlag)
            stepLabel = nniStepThisTime ? "(nni)" : "(spr, radius " + stepRadiusText + ")";
        else
            stepLabel = "(radius " + stepRadiusText + ")";

        PhyloNode *bestNode, *bestDad;
        int bestDistance;
        double bestScore;
        // bestAchievedRadiusForLearning: what learnradiusFlag's own window
        // actually records for this step, once accepted -- bestDistance
        // (a hop count) everywhere EXCEPT "fast"+"distradius", where it's
        // the equivalent radiusPercent instead (chooseGraftByDistance's own
        // outPercentUsed), to stay in the same units as effectiveRadiusContinuous
        // itself. See chooseGraftByDistance's outPercentUsed comment.
        double bestAchievedRadiusForLearning = 0.0;

        if (opt.useFastSelection && !investigatingThisStep) {
            // O(distance) proposal: draw numCandidates independent
            // candidates from chooseGraft's random walk, ALL from this
            // same step's prune position, each applied/scored/rolled-back
            // just like each of the exhaustive branch's many candidates
            // below, and keep the best of the group. numCandidates == 1
            // (plain "fast", the default) is exactly the original single-
            // candidate behavior; "fast N" trades away some of fast
            // mode's cheapness for a chance at a better move each step,
            // without paying to score every candidate within the radius
            // the way the exhaustive (non-fast) branch does.
            bool haveBest = false;
            bestScore = -DBL_MAX;
            for (int c = 0; c < opt.numCandidates; c++) {
                PhyloNode *candNode, *candDad;
                int walkLength;
                double walkPercentUsed = 0.0;
                // useDistanceRadius ("distradius") swaps in
                // chooseGraftByDistance, which reinterprets
                // effectiveRadiusContinuous as a percentage of the tree's
                // total branch length and walks by actual summed branch
                // length instead of hop count -- see its own comment for
                // the full mechanics
                bool found = opt.useDistanceRadius
                    ? chooseGraftByDistance(tree, pruneNode, pruneDad, effectiveRadiusContinuous, candNode, candDad,
                            &walkLength, &walkPercentUsed)
                    : chooseGraft(tree, pruneNode, pruneDad, effectiveRadius, candNode, candDad, &walkLength);
                if (!found)
                    continue; // this draw found no legal target; try the next one

                SPRMove move;
                move.prune_node = pruneNode;
                move.prune_dad = pruneDad;
                move.regraft_node = candNode;
                move.regraft_dad = candDad;
                move.radius = walkLength;
                move.screening_score = 0.0;
                move.exact_score = 0.0;
                move.candidate_id = c;
                move.generation = step;

                double score = scoreTrialSPRMove(tree, move, opt.reoptimizeBranchLengths,
                        localCache);
                st.candidatesEvaluated++;

                if (!haveBest || score > bestScore) {
                    bestScore = score;
                    bestNode = candNode;
                    bestDad = candDad;
                    bestDistance = walkLength;
                    bestAchievedRadiusForLearning = opt.useDistanceRadius ? walkPercentUsed : (double) walkLength;
                    haveBest = true;
                }
            }
            if (!haveBest) {
                if (!opt.quiet)
                    cout << stepLabelPrefix << "step " << (st.stepsRun + step + 1) << " " << stepLabel << ": prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft target found in " << opt.numCandidates << " draw(s); skipping." << endl;
                // finding no legal candidate at all also ends any
                // currently-open "learnradius"+"investigate" excursion --
                // see maybeFinalizeLearnRadiusExcursion's comment
                maybeFinalizeLearnRadiusExcursion(tree, opt.useDistanceRadius, st.learnRadiusExcursionOpen, st.learnRadiusAnchorA, st.learnRadiusAnchorB,
                        st.learnRadiusFinalNode, pruneDad, st.learnRadiusWindow, opt.learnradiusN);
                continue;
            }
        } else {
            // investigatingThisStep forces investigateRadius (default 1,
            // findGraftPositions' own nearest-legal-target tier -- an
            // NNI-equivalent move) regardless of stepRadius/useFastSelection;
            // otherwise effectiveRadius already reflects "alternate"'s own
            // step-parity NNI toggle, or just equals stepRadius unchanged
            int graftRadius = investigatingThisStep ? opt.investigateRadius : effectiveRadius;
            vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, graftRadius);
            if (candidates.empty()) {
                if (!opt.quiet)
                    cout << stepLabelPrefix << "step " << (st.stepsRun + step + 1) << " " << stepLabel << ": prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft candidates; skipping." << endl;
                // finding no legal candidate at all also ends any
                // currently-open "learnradius"+"investigate" excursion --
                // see maybeFinalizeLearnRadiusExcursion's comment
                maybeFinalizeLearnRadiusExcursion(tree, opt.useDistanceRadius, st.learnRadiusExcursionOpen, st.learnRadiusAnchorA, st.learnRadiusAnchorB,
                        st.learnRadiusFinalNode, pruneDad, st.learnRadiusWindow, opt.learnradiusN);
                continue;
            }

            // score every candidate on the SAME tree object via apply ->
            // score -> rollback; never allocate a new tree per candidate.
            bestScore = -DBL_MAX;
            GraftCandidate bestCandidate = candidates[0];
            for (size_t i = 0; i < candidates.size(); i++) {
                const GraftCandidate &c = candidates[i];
                SPRMove move;
                move.prune_node = pruneNode;
                move.prune_dad = pruneDad;
                move.regraft_node = c.node;
                move.regraft_dad = c.dad;
                move.radius = c.radius;
                move.screening_score = 0.0;
                move.exact_score = 0.0;
                move.candidate_id = (int) i;
                move.generation = step;

                double score = scoreTrialSPRMove(tree, move, opt.reoptimizeBranchLengths,
                        localCache);
                st.candidatesEvaluated++;

                if (i == 0 || score > bestScore) {
                    bestScore = score;
                    bestCandidate = c;
                }
            }
            bestNode = bestCandidate.node;
            bestDad = bestCandidate.dad;
            bestDistance = bestCandidate.radius;
            // findGraftPositions is always hop-based, "distradius" or not
            // (it "has no distance-based counterpart" -- see distradius'
            // own comment), so this branch's own achieved radius is
            // always just bestDistance, same as before
            bestAchievedRadiusForLearning = (double) bestDistance;
        }

        // apply the winning candidate for real, once, to decide whether to
        // keep it
        SPRMove bestMove;
        bestMove.prune_node = pruneNode;
        bestMove.prune_dad = pruneDad;
        bestMove.regraft_node = bestNode;
        bestMove.regraft_dad = bestDad;
        bestMove.radius = bestDistance;
        bestMove.screening_score = 0.0;
        bestMove.exact_score = bestScore;
        bestMove.candidate_id = 0;
        bestMove.generation = step;

        TrackedSPR bestTracked;
        applySPRTracked(tree, reg, bestMove, bestTracked);
        bool recomputedAppliedTopology = false;

        if (opt.reoptimizeBranchLengths) {
            resetLikelihoodBuffers(tree);
            // scoreTrialSPRMove's own reoptimization (used to pick this
            // candidate as the step's best) always gets rolled back along
            // with everything else once scoring is done -- redo it here,
            // on the tree as just permanently applied for real, so the
            // optimized lengths actually persist into the tree this step
            // keeps (or are discarded along with everything else below if
            // this candidate turns out not to improve on curScore after all)
            reoptimizeSPREdges(tree, bestMove.prune_dad, bestMove.regraft_dad, bestMove.regraft_node,
                    bestMove.prune_node, bestTracked.sibling1, bestTracked.sibling2);
            resetLikelihoodBuffers(tree);
            double realScore = tree.computeLikelihood();
            if (std::isfinite(realScore))
                bestScore = realScore;
            recomputedAppliedTopology = true;
        } else if (bestScore > curScore || opt.acceptDist.enabled) {
            // The local score is sufficient to reject a non-improving
            // proposal without touching the baseline cache. For a proposed
            // improvement, recompute the winner once from scratch before
            // committing it. This is both a correctness boundary and the
            // fully-populated baseline cache for the next step.
            resetLikelihoodBuffers(tree);
            double realScore = tree.computeLikelihood();
            if (std::isfinite(realScore))
                bestScore = realScore;
            recomputedAppliedTopology = true;
        }

        bool improved = bestScore > curScore;
        // --accept-dist: a losing move may still be kept, with probability
        // set by the rule. `improved` keeps meaning a genuine gain, since
        // shrink's stall counter and learnradius's excursion bookkeeping
        // both key off real progress and would be corrupted by counting a
        // downhill step as one.
        bool accepted = improved;
        if (!improved && opt.acceptDist.enabled && std::isfinite(bestScore)) {
            st.offeredDownhill++;
            st.lastTemperature = opt.acceptDist.effectiveTemperature(acceptProgress);
            if (opt.acceptDist.accept(bestScore - curScore, acceptProgress)) {
                accepted = true;
                st.acceptedDownhill++;
            }
        }
        if (opt.shrinkFlag)
            maybeShrinkRadius(improved, opt.shrinkStallThreshold, opt.quiet,
                    st.shrinkStallCount, st.shrinkCurrentRadius);
        if (!opt.quiet)
            cout << stepLabelPrefix << "step " << (st.stepsRun + step + 1) << " " << stepLabel << ": prune {"
                 << describeEdgeCompact(pruneNode, pruneDad) << "}"
                 << " -> graft {" << describeEdgeCompact(bestNode, bestDad) << "}"
                 << " (" << ((opt.useFastSelection && !investigatingThisStep) ? "d=" : "distance ") << bestDistance << ")"
                 << ", logL " << bestScore << " (cur " << curScore << ")"
                 << (improved ? " [kept]" : (accepted ? " [kept: accept-dist]" : " [reverted]")) << endl;

        if (accepted) {
            curScore = bestScore;
            // --accept-dist walks downhill, so the walk's endpoint is
            // routinely worse than the best tree it passed through. Keep
            // that high-water mark; runSPRSteps restores it on the way out.
            if (opt.acceptDist.enabled && curScore > st.bestSeenScore) {
                st.bestSeenScore = curScore;
                st.bestSeenTree = tree.getTreeString();
            }
            if (improved && opt.learnradiusFlag) {
                if (opt.investigateFlag) {
                    // see learnradiusFlag's comment on runHillClimb: a
                    // "investigate"-chained excursion is measured as ONE
                    // net displacement once it's fully done, not per
                    // individual attempt
                    if (!investigatingThisStep) {
                        // this move starts a fresh excursion -- remember
                        // where it started: the edge applySPRTracked's own
                        // TrackedSPR leaves behind at pruneDad's pre-move
                        // position (its sibling1/sibling2), exactly the
                        // same collapsed view findGraftPositions/
                        // chooseGraft themselves already use for radius 1
                        st.learnRadiusExcursionOpen = true;
                        st.learnRadiusAnchorA = bestTracked.sibling1;
                        st.learnRadiusAnchorB = bestTracked.sibling2;
                    }
                    // every accepted step in the excursion (the initiating
                    // one included) moves its own current resting point
                    st.learnRadiusFinalNode = bestNode;
                } else {
                    // investigateFlag is off -- no excursion to chain,
                    // this accepted move IS the whole (length-1)
                    // "excursion" already; bestAchievedRadiusForLearning is
                    // already in the right units either way (hop count, or
                    // radiusPercent-equivalent under "distradius")
                    recordLearnRadiusSample(st.learnRadiusWindow, opt.learnradiusN, bestAchievedRadiusForLearning);
                }
            }
            if (improved && opt.investigateFlag) {
                // this move -- whether it came from a fresh choosePrune or
                // from investigating a previous one -- just improved the
                // tree, so refine IT one real hop further next step; see
                // investigateFlag's comment on runHillClimb
                st.investigateNext = true;
                st.investigatePruneNode = pruneNode;
                st.investigatePruneDad = pruneDad;
            }
            if (opt.recordProgress)
                appendRecordRow(st.modelName, st.recordTag, st.runId, st.candidatesEvaluated, getCPUTime() - st.cpuClockStart,
                        curScore, st.trueTreeLogl, opt.recordTopology, tree);
            // "trajectory": one more Newick line for this just-accepted
            // move's resulting topology -- see appendTrajectoryTopology's
            // comment on why only accepted moves (never a periodic
            // fullreopt/findopt refit, which never touch the topology)
            // ever add a line here.
            if (opt.trajectoryFlag)
                appendTrajectoryTopology(st.runId, tree);
            // "fullreopt": counts this accepted move toward
            // fullReoptEveryNSteps' own successful-step cadence, then
            // fires the periodic sweep the moment that count crosses the
            // next multiple -- see maybeRunPeriodicFullReopt's comment for
            // why this is gated on successfulSteps rather than the loop's
            // own step index.
            st.successfulSteps++;
            maybeRunPeriodicFullReopt(tree, st.successfulSteps, opt.fullReoptEveryNSteps, opt.fullReoptRounds, opt.useGtrModel,
                    opt.quiet, opt.recordProgress, opt.recordTopology, st.modelName, st.recordTag, st.runId, st.candidatesEvaluated,
                    st.cpuClockStart, st.trueTreeLogl, curScore);
        } else if (opt.tunnelFlag && std::isfinite(bestScore)
                   && (curScore - bestScore) <= opt.tunnelTolerance) {
            // "tunnel": this move is worse, but only marginally. Keep it
            // applied and probe random follow-ups one at a time, looking
            // for a PAIR that beats curScore. Nothing is kept unless it
            // strictly improves on where the tree already was, so a
            // tolerated first move can never leak into the result on its
            // own.
            st.tunnelEntered++;
            bool committed = false;
            for (int probe = 0; probe < opt.tunnelTries && !committed; probe++) {
                PhyloNode *pNode, *pDad;
                if (!choosePrune(tree, reg, pNode, pDad, opt.weightpruneFlag))
                    break;
                PhyloNode *gNode, *gDad;
                bool found = opt.useDistanceRadius
                    ? chooseGraftByDistance(tree, pNode, pDad, (double) opt.radius,
                            gNode, gDad, nullptr, nullptr)
                    : chooseGraft(tree, pNode, pDad, opt.radius, gNode, gDad, nullptr);
                if (!found)
                    continue;

                SPRMove probeMove;
                probeMove.prune_node = pNode;
                probeMove.prune_dad = pDad;
                probeMove.regraft_node = gNode;
                probeMove.regraft_dad = gDad;
                probeMove.radius = 0;
                probeMove.screening_score = 0.0;
                probeMove.exact_score = 0.0;
                probeMove.candidate_id = 0;
                probeMove.generation = step;

                TrackedSPR probeTracked;
                applySPRTracked(tree, reg, probeMove, probeTracked);
                if (opt.reoptimizeBranchLengths) {
                    resetLikelihoodBuffers(tree);
                    reoptimizeSPREdges(tree, probeMove.prune_dad, probeMove.regraft_dad,
                            probeMove.regraft_node, probeMove.prune_node,
                            probeTracked.sibling1, probeTracked.sibling2);
                }
                resetLikelihoodBuffers(tree);
                double pairScore = tree.computeLikelihood();
                st.candidatesEvaluated++;
                st.tunnelProbes++;

                if (std::isfinite(pairScore) && pairScore > curScore) {
                    // the pair as a whole is an improvement -- commit both
                    if (!opt.quiet)
                        cout << stepLabelPrefix << "  tunnel: probe " << (probe + 1)
                             << " cleared the barrier, logL " << pairScore
                             << " (was " << curScore << ", via "
                             << (curScore - bestScore) << " downhill)" << endl;
                    curScore = pairScore;
                    committed = true;
                    st.tunnelCommitted++;
                    if (opt.recordProgress)
                        appendRecordRow(st.modelName, st.recordTag, st.runId, st.candidatesEvaluated,
                                getCPUTime() - st.cpuClockStart, curScore, st.trueTreeLogl,
                                opt.recordTopology, tree);
                    if (opt.trajectoryFlag)
                        appendTrajectoryTopology(st.runId, tree);
                    st.successfulSteps++;
                    maybeRunPeriodicFullReopt(tree, st.successfulSteps, opt.fullReoptEveryNSteps,
                            opt.fullReoptRounds, opt.useGtrModel, opt.quiet, opt.recordProgress,
                            opt.recordTopology, st.modelName, st.recordTag, st.runId,
                            st.candidatesEvaluated, st.cpuClockStart, st.trueTreeLogl, curScore);
                } else {
                    // PARTIAL rollback: undo only this probe, leaving the
                    // tolerated first move in place so the next probe
                    // starts from the same position. Rollbacks are LIFO, so
                    // the probe must come off before bestTracked can.
                    rollbackSPRTracked(tree, reg, probeTracked);
                }
            }
            if (!committed) {
                rollbackSPRTracked(tree, reg, bestTracked);
                maybeFinalizeLearnRadiusExcursion(tree, opt.useDistanceRadius, st.learnRadiusExcursionOpen,
                        st.learnRadiusAnchorA, st.learnRadiusAnchorB, st.learnRadiusFinalNode, pruneDad,
                        st.learnRadiusWindow, opt.learnradiusN);
            }
            // Either path left the tree topologically settled but its
            // likelihood buffers dirty from the probing above; the next
            // step's local scoring reads that cache as its baseline.
            resetLikelihoodBuffers(tree);
            tree.computeLikelihood();
        } else {
            rollbackSPRTracked(tree, reg, bestTracked);
            if (recomputedAppliedTopology) {
                resetLikelihoodBuffers(tree);
                tree.computeLikelihood();
            }
            // this investigation attempt just failed to improve -- if an
            // excursion was open, it's over now (tree just rolled back to
            // exactly where the last successful move in it left things);
            // see maybeFinalizeLearnRadiusExcursion's comment
            maybeFinalizeLearnRadiusExcursion(tree, opt.useDistanceRadius, st.learnRadiusExcursionOpen, st.learnRadiusAnchorA, st.learnRadiusAnchorB,
                    st.learnRadiusFinalNode, pruneDad, st.learnRadiusWindow, opt.learnradiusN);
        }

        // maybeRunPeriodicFullReopt is no longer called unconditionally
        // here -- it now fires from inside the `if (improved)` branch
        // above, gated on successfulSteps rather than this loop's own step
        // index; see its own comment for why. findopt stays on the
        // original per-step-attempt cadence: it's a pure diagnostic that
        // never touches curScore or the tree, so there's no "successful
        // step" concept for it to dilute.
        maybeRunFindopt(tree, st.stepsRun + step, opt.findoptEveryNSteps, opt.quiet, opt.recordProgress, opt.recordTopology, st.modelName, st.recordTag,
                st.runId, st.candidatesEvaluated, st.cpuClockStart, st.trueTreeLoglForFindopt, curScore, aln, params);
    }

    // a "learnradius"+"investigate" excursion still open when the loop
    // above ended (<max-steps> reached, or "no degree-3 node left to
    // prune") rather than via a failed investigation attempt -- finalize
    // it here instead of discarding it; investigatePruneDad (rather than
    // the now out-of-scope loop-local pruneDad) is exactly where it
    // currently sits, since it's updated every time an excursion continues
    maybeFinalizeLearnRadiusExcursion(tree, opt.useDistanceRadius, st.learnRadiusExcursionOpen, st.learnRadiusAnchorA, st.learnRadiusAnchorB,
            st.learnRadiusFinalNode, st.investigatePruneDad, st.learnRadiusWindow, opt.learnradiusN);
    st.stepsRun += step;
    return curScore;
}

double runSPRSweep(PhyloTree &tree, EdgeRegistry &reg, const SPRSearchOptions &opt, SPRSearchState &st,
        double curScore) {
    // "sweep": post-processing phase, run strictly AFTER the step loop
    // above has finished (with whatever mix of fast/investigate/
    // alternate/shrink shaped those steps) -- see sweepFlag's comment on
    // runHillClimb for the full rationale. Unlike every flag above, this
    // never competes for a STEP's own selection logic, so it needs no
    // incompatibility check with any of them.
    if (opt.sweepFlag) {
        // Step 1: score every internal edge (both endpoints degree 3 -- an
        // edge to a leaf has no quartet to test, see
        // computeSiblingCompatibilityScore) by how well its CURRENT
        // grouping is supported relative to the two alternative NNI
        // rearrangements around it. Lower (more negative) means less
        // compatible -- these are the "least-compatible siblings" sweepCount
        // targets.
        vector<pair<double, pair<PhyloNode*, PhyloNode*> > > ranked;
        for (size_t slot = 0; slot < reg.slots.size(); slot++) {
            PhyloNode *p = reg.slots[slot].first;
            PhyloNode *q = reg.slots[slot].second;
            double score;
            PhyloNode *pruneNode, *pruneDad;
            if (computeSiblingCompatibilityScore(tree, p, q, curScore, score, pruneNode, pruneDad))
                ranked.push_back(make_pair(score, make_pair(pruneNode, pruneDad)));
        }
        sort(ranked.begin(), ranked.end(),
                [](const pair<double, pair<PhyloNode*, PhyloNode*> > &x,
                   const pair<double, pair<PhyloNode*, PhyloNode*> > &y) {
                    return x.first < y.first;
                });

        int sweepTargets = min(opt.sweepCount, (int) ranked.size());
        if (!opt.quiet)
            cout << endl << "sweep: " << ranked.size() << " internal edge(s) scored ("
                 << "adjacent_subtree_compatibility.pdf, sec. 8) -- targeting the " << sweepTargets
                 << " least-compatible sibling pair(s)" << endl;

        // Step 2: visit those sweepTargets positions one at a time, each
        // with an exhaustive, whole-tree regraft search (reusing exactly
        // the same findGraftPositions-at-full-radius + scoreTrialSPRMove +
        // accept/reject machinery the ORIGINAL one-pass-over-every-edge
        // "sweep" design used per edge -- see sweepFlag's comment on
        // runHillClimb). This is a FIXED list, computed once up front, not
        // recomputed after each improvement -- exactly the same "not
        // checking from the start after each improvement" limitation the
        // original design had, so this is not guaranteed to reach the
        // SPR-optimal tree either.
        int sweepFullRadius = (int) reg.slots.size();
        for (int k = 0; k < sweepTargets; k++) {
            PhyloNode *pruneNode = ranked[k].second.first;
            PhyloNode *pruneDad = ranked[k].second.second;

            // an earlier position in THIS pass may have moved a subtree
            // that disturbed this specific prune point (e.g. pruneNode sat
            // inside a subtree relocated by an earlier accepted move in
            // this same pass, or WAS itself relocated) -- skip it rather
            // than act on a stale pair of node pointers. pruneDad's own
            // degree never changes (applySPR only repositions it, see
            // applySPRTracked's comment), so only adjacency needs checking.
            if (!pruneNode->isNeighbor(pruneDad))
                continue;

            string sweepLabel = "sweep " + to_string(k + 1) + "/" + to_string(sweepTargets)
                    + " (compat " + to_string(ranked[k].first) + ")";

            vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, sweepFullRadius);
            if (candidates.empty()) {
                if (!opt.quiet)
                    cout << sweepLabel << ": prune {" << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft candidates; skipping." << endl;
                continue;
            }

            double bestScore = -DBL_MAX;
            GraftCandidate bestCandidate = candidates[0];
            for (size_t i = 0; i < candidates.size(); i++) {
                const GraftCandidate &c = candidates[i];
                SPRMove move;
                move.prune_node = pruneNode;
                move.prune_dad = pruneDad;
                move.regraft_node = c.node;
                move.regraft_dad = c.dad;
                move.radius = c.radius;
                move.screening_score = 0.0;
                move.exact_score = 0.0;
                move.candidate_id = (int) i;
                move.generation = -1; // not part of the step loop above

                double score = scoreTrialSPRMove(tree, move, opt.reoptimizeBranchLengths);
                st.candidatesEvaluated++;

                if (i == 0 || score > bestScore) {
                    bestScore = score;
                    bestCandidate = c;
                }
            }

            SPRMove bestMove;
            bestMove.prune_node = pruneNode;
            bestMove.prune_dad = pruneDad;
            bestMove.regraft_node = bestCandidate.node;
            bestMove.regraft_dad = bestCandidate.dad;
            bestMove.radius = bestCandidate.radius;
            bestMove.screening_score = 0.0;
            bestMove.exact_score = bestScore;
            bestMove.candidate_id = 0;
            bestMove.generation = -1;

            TrackedSPR bestTracked;
            applySPRTracked(tree, reg, bestMove, bestTracked);
            resetLikelihoodBuffers(tree);

            if (opt.reoptimizeBranchLengths) {
                reoptimizeSPREdges(tree, bestMove.prune_dad, bestMove.regraft_dad, bestMove.regraft_node,
                        bestMove.prune_node, bestTracked.sibling1, bestTracked.sibling2);
                resetLikelihoodBuffers(tree);
                double realScore = tree.computeLikelihood();
                if (std::isfinite(realScore))
                    bestScore = realScore;
            }

            bool improvedSweep = bestScore > curScore;
            if (!opt.quiet)
                cout << sweepLabel << ": prune {" << describeEdgeCompact(pruneNode, pruneDad) << "}"
                     << " -> graft {" << describeEdgeCompact(bestCandidate.node, bestCandidate.dad) << "}"
                     << " (distance " << bestCandidate.radius << ")"
                     << ", logL " << bestScore << " (cur " << curScore << ")"
                     << (improvedSweep ? " [kept]" : " [reverted]") << endl;

            if (improvedSweep) {
                curScore = bestScore;
                if (opt.recordProgress)
                    appendRecordRow(st.modelName, st.recordTag, st.runId, st.candidatesEvaluated, getCPUTime() - st.cpuClockStart,
                            curScore, st.trueTreeLogl, opt.recordTopology, tree);
            } else {
                rollbackSPRTracked(tree, reg, bestTracked);
                resetLikelihoodBuffers(tree);
            }
        }
    }
    return curScore;
}

void appendFinalRecordRow(const SPRSearchState &st, const SPRSearchOptions &opt, PhyloTree &tree, double curScore) {
    // unconditional final row: the step loop only calls appendRecordRow on
    // an ACCEPTED step (or a periodic full-reopt event), and the sweep phase
    // only on an accepted sweep position -- so if the very last event
    // (whichever phase it came from) was rejected, the CSV's last row would
    // otherwise be stale relative to when the run actually stopped. This
    // guarantees one row always reflects the true final state (curScore as
    // it stands once EVERYTHING is done), regardless of whether that last
    // event was ever kept. findopt is NOT part of this guarantee -- unlike
    // the old, curScore-mutating "finalreopt" it replaced, findopt never
    // changes curScore at all, so it has no result for this final row to
    // reflect; its own periodic checkpoints each write their own row.
    if (!opt.recordProgress)
        return;
    appendRecordRow(st.modelName, st.recordTag, st.runId, st.candidatesEvaluated,
            getCPUTime() - st.cpuClockStart, curScore, st.trueTreeLogl, opt.recordTopology, tree);
}

/*==========================================================================
    --spr-refine / --nni-refine spec parsing
 *========================================================================*/

namespace {

/**
    split a flag spec on whitespace. The spec arrives as ONE command-line
    argument ("radius 10 fast quiet"), so the shell has already done no
    splitting of its own.
 */
vector<string> tokenizeSpec(const string &spec) {
    vector<string> out;
    istringstream in(spec);
    string tok;
    while (in >> tok)
        out.push_back(tok);
    return out;
}

/**
    read the OPTIONAL positive-integer argument that may follow a flag,
    consuming it only if the next token really is one (so "fast quiet"
    leaves "quiet" alone, while "fast 5" takes the 5). Returns whether a
    value was consumed; `out` is left untouched if not.
 */
bool takeOptionalInt(const vector<string> &tok, size_t &i, int &out) {
    if (i + 1 >= tok.size())
        return false;
    const string &next = tok[i + 1];
    char *end = nullptr;
    long v = strtol(next.c_str(), &end, 10);
    if (end == next.c_str() || *end != '\0' || v < 1)
        return false;
    out = (int) v;
    i++;
    return true;
}

/**
    read a REQUIRED positive-integer argument. Unlike takeOptionalInt, a
    missing or non-numeric token here is an error, not a flag without an
    argument.
 */
bool takeRequiredInt(const vector<string> &tok, size_t &i, const string &flag, int &out, string &err) {
    if (i + 1 >= tok.size()) {
        err = "'" + flag + "' needs a number after it";
        return false;
    }
    const string &next = tok[i + 1];
    char *end = nullptr;
    long v = strtol(next.c_str(), &end, 10);
    if (end == next.c_str() || *end != '\0' || v < 1) {
        err = "'" + flag + "' needs a positive whole number, got '" + next + "'";
        return false;
    }
    out = (int) v;
    i++;
    return true;
}

} // namespace

bool parsePerturbSpec(const string &spec, SPRSearchOptions &opt, string &err) {
    // A perturbation is a handful of random moves, not a search, so only
    // the flags that describe a MOVE mean anything: how far it reaches
    // (radius, distradius) and which edge it prunes (weightprune). There
    // record -- the refinement stage that follows owns all of that -- so
    // every other flag is refused rather than silently ignored. The one
    // exception is "slack", which DOES score each draw, precisely so a
    // kick's damage can be bounded; see SPRSearchOptions::slackFlag.
    opt = SPRSearchOptions();
    opt.noTrueTree = true;
    // a perturbation always draws its target by random walk; the
    // exhaustive enumeration path has no meaning when nothing is scored
    opt.useFastSelection = true;
    opt.numCandidates = 1;

    vector<string> tok = tokenizeSpec(spec);
    for (size_t i = 0; i < tok.size(); i++) {
        const string &t = tok[i];
        if (t == "radius") {
            if (!takeRequiredInt(tok, i, t, opt.radius, err))
                return false;
            continue;
        }
        if (t == "distradius") {
            opt.useDistanceRadius = true;
            continue;
        }
        if (t == "slack") {
            // "slack D [anneal]": bound this kick's damage to D
            // log-likelihood per move rather than applying moves blind. D
            // must be positive; a zero threshold would refuse every
            // downhill draw, which is not a perturbation at all.
            if (i + 1 >= tok.size()) {
                err = "\"slack\" needs a delta, e.g. \"slack 1\" or \"slack 5 anneal\"";
                return false;
            }
            opt.slackFlag = true;
            i++;
            opt.slackDelta = atof(tok[i].c_str());
            if (opt.slackDelta <= 0.0) {
                err = "\"slack\" delta must be positive (got " + tok[i] + ")";
                return false;
            }
            if (i + 1 < tok.size() && tok[i + 1] == "anneal") {
                opt.slackAnneal = true;
                i++;
            }
            continue;
        }
        if (t == "weightprune") {
            opt.weightpruneFlag = PRUNE_LONG;
            if (i + 1 < tok.size() && (tok[i + 1] == "long" || tok[i + 1] == "short")) {
                opt.weightpruneFlag = (tok[i + 1] == "short") ? PRUNE_SHORT : PRUNE_LONG;
                i++;
            }
            if (opt.weightpruneFlag == PRUNE_SHORT)
                cout << "WARNING: 'weightprune short' is deprecated -- it measured markedly WORSE"
                        " than both an unweighted prune and 'weightprune long' (see PruneWeighting"
                        " in search_experiments/sprsearch.h). Use 'weightprune long', or drop the flag." << endl;
            continue;
        }
        err = "'" + t + "' does not describe an SPR MOVE, so it means nothing to a"
                " perturbation; --spr-perturb takes only radius/distradius/weightprune"
                " (the kick's strength is --perturb, as it is for the NNI kick)";
        return false;
    }
    return true;
}

int doRandomSPRs(PhyloTree &tree, const SPRSearchOptions &opt, int numMoves,
        SPRSearchState *st, double annealScale) {
    // The SPR counterpart of IQTree::doRandomNNIs: apply `numMoves` random,
    // legal SPR moves and keep none of the scoring machinery -- a kick is
    // supposed to make the tree worse, so nothing here is evaluated,
    // compared, or rolled back. Branch lengths are left exactly as applySPR
    // sets them, matching what the NNI kick does (it re-scores but never
    // re-optimizes; see IQTree::doTreePerturbation).
    EdgeRegistry reg;
    buildEdgeRegistry(tree, reg);

    // "slack": the kick stops being blind. Every proposal is scored and
    // kept only if it costs at most `delta` log-likelihood; anything worse
    // is rolled back and the slot redrawn. This bounds how far a kick can
    // damage the tree, instead of letting numMoves blind moves carry it
    // arbitrarily far from the basin the refinement has to climb back up.
    double delta = 0.0;
    double curScore = 0.0;
    const int maxDrawsPerMove = 8;
    if (opt.slackFlag) {
        delta = opt.slackDelta * (opt.slackAnneal ? annealScale : 1.0);
        if (delta < 0.0)
            delta = 0.0;
        resetLikelihoodBuffers(tree);
        curScore = tree.computeLikelihood();
        if (st) {
            st->slackKicks++;
            st->slackLastDelta = delta;
        }
    }

    int applied = 0;
    int redraws = 0;
    for (int i = 0; i < numMoves; i++) {
        PhyloNode *pruneNode, *pruneDad;
        if (!choosePrune(tree, reg, pruneNode, pruneDad, opt.weightpruneFlag))
            break; // no degree-3 node left to prune from

        PhyloNode *graftNode, *graftDad;
        bool found = opt.useDistanceRadius
            ? chooseGraftByDistance(tree, pruneNode, pruneDad, (double) opt.radius,
                    graftNode, graftDad, nullptr, nullptr)
            : chooseGraft(tree, pruneNode, pruneDad, opt.radius, graftNode, graftDad, nullptr);
        if (!found)
            continue; // this draw found no legal target; try another prune point

        SPRMove move;
        move.prune_node = pruneNode;
        move.prune_dad = pruneDad;
        move.regraft_node = graftNode;
        move.regraft_dad = graftDad;
        move.radius = 0;
        move.screening_score = 0.0;
        move.exact_score = 0.0;
        move.candidate_id = 0;
        move.generation = i;

        // applySPRTracked rather than a bare applySPR: the registry has to
        // stay in sync for the NEXT move's choosePrune to draw from a
        // truthful edge list. The TrackedSPR it fills is discarded -- a
        // perturbation never rolls back.
        TrackedSPR tracked;
        applySPRTracked(tree, reg, move, tracked);

        if (opt.slackFlag) {
            resetLikelihoodBuffers(tree);
            double newScore = tree.computeLikelihood();
            if (std::isfinite(newScore) && newScore > curScore - delta) {
                curScore = newScore;
                if (st) st->slackAccepted++;
                redraws = 0;
            } else {
                // costs more than the threshold allows -- undo it and
                // redraw this slot, up to a bounded number of attempts so a
                // tight delta cannot spin here forever
                rollbackSPRTracked(tree, reg, tracked);
                if (st) st->slackRejected++;
                if (redraws < maxDrawsPerMove) {
                    redraws++;
                    i--;            // slot not filled yet; try again
                }
                else {
                    redraws = 0;    // give up on this slot, move on
                }
                continue;
            }
        }
        applied++;
    }
    if (opt.slackFlag) {
        // leave the buffers consistent with the tree the caller now holds
        resetLikelihoodBuffers(tree);
        tree.computeLikelihood();
    }
    return applied;
}

bool parseRefineSpec(const string &spec, bool sprMode, SPRSearchOptions &opt, string &err) {
    opt = SPRSearchOptions();
    // IQ-TREE never has a ground-truth tree to measure a gap against, so
    // the record CSV's "gap to true tree" column is always NaN here --
    // exactly what --hillclimb's own "notree" mode already does.
    opt.noTrueTree = true;
    // A periodic full re-optimization under "fullreopt" fits model
    // parameters as well as branch lengths, matching what doNNISearch
    // already does whenever an iteration improves on the best score. (In
    // spr_topology_test this same field means "search under the richer
    // model"; the model itself is -m's job here, so only the
    // optimize-parameters-too half of it carries over.)
    opt.useGtrModel = true;

    vector<string> tok = tokenizeSpec(spec);
    // flags that describe a whole standalone run rather than one
    // refinement pass. IQ-TREE already owns every one of these jobs, and
    // silently ignoring a starting-tree or model flag would be worse than
    // refusing it -- the user would reasonably believe it took effect.
    const char *whole_run[][2] = {
        {"random", "IQ-TREE builds its own starting trees; see -t and --ninit"},
        {"iqtreestart", "IQ-TREE builds its own starting trees; see -t and --ninit"},
        {"starttree", "IQ-TREE builds its own starting trees; see -t"},
        {"model", "the model is IQ-TREE's own; see -m"},
        {"gtr", "the model is IQ-TREE's own; see -m"},
        {"notree", "always implied here -- IQ-TREE has no ground-truth tree to compare against"},
    };

    for (size_t i = 0; i < tok.size(); i++) {
        const string &t = tok[i];

        bool rejected = false;
        for (size_t w = 0; w < sizeof(whole_run) / sizeof(whole_run[0]); w++) {
            if (t == whole_run[w][0]) {
                err = "'" + t + "' configures a whole standalone spr_topology_test run, not a"
                        " refinement pass: " + whole_run[w][1];
                rejected = true;
                break;
            }
        }
        if (rejected)
            return false;

        // --- accepted in both modes: these describe what gets reported,
        // --- not how a step is chosen, so they carry over to NNI unchanged
        if (t == "quiet") {
            opt.quiet = true;
            continue;
        }
        if (t == "record") {
            opt.recordProgress = true;
            continue;
        }
        if (t == "trajectory") {
            opt.trajectoryFlag = true;
            continue;
        }
        if (t == "findopt") {
            opt.findoptFlag = true;
            // in IQ-TREE the cadence counts PERTURBATIONS, not SPR steps:
            // one diagnostic refit every N kicks. Default: every kick.
            opt.findoptEveryNSteps = 1;
            takeOptionalInt(tok, i, opt.findoptEveryNSteps);
            continue;
        }

        if (!sprMode) {
            err = "'" + t + "' shapes an SPR step, so it only means something under"
                    " --spr-refine; --nni-refine takes only quiet/record/trajectory/findopt";
            return false;
        }

        // --- SPR-only from here down ---
        if (t == "radius") {
            if (!takeRequiredInt(tok, i, t, opt.radius, err))
                return false;
            continue;
        }
        if (t == "steps") {
            if (!takeRequiredInt(tok, i, t, opt.stepsPerPass, err))
                return false;
            continue;
        }
        if (t == "fast") {
            opt.useFastSelection = true;
            takeOptionalInt(tok, i, opt.numCandidates);
            continue;
        }
        if (t == "reopt") {
            opt.reoptimizeBranchLengths = true;
            continue;
        }
        if (t == "fullreopt") {
            // "fullreopt M N": both required, same as --hillclimb. Its
            // optional third "true" token is an up-front whole-tree fit of
            // the STARTING tree, which is a starting-tree concern -- see
            // the whole_run list above -- so it isn't accepted here.
            if (!takeRequiredInt(tok, i, "fullreopt", opt.fullReoptRounds, err))
                return false;
            if (!takeRequiredInt(tok, i, "fullreopt", opt.fullReoptEveryNSteps, err))
                return false;
            continue;
        }
        if (t == "investigate") {
            opt.investigateFlag = true;
            takeOptionalInt(tok, i, opt.investigateRadius);
            continue;
        }
        if (t == "alternate") {
            opt.alternateFlag = true;
            continue;
        }
        if (t == "shrink") {
            opt.shrinkFlag = true;
            takeOptionalInt(tok, i, opt.shrinkStallThreshold);
            continue;
        }
        if (t == "learnradius") {
            opt.learnradiusFlag = true;
            takeOptionalInt(tok, i, opt.learnradiusN);
            continue;
        }
        if (t == "sweep") {
            opt.sweepFlag = true;
            takeOptionalInt(tok, i, opt.sweepCount);
            continue;
        }
        if (t == "slack") {
            err = "\"slack\" bounds how far a PERTURBATION may go downhill, "
                  "so it belongs to --spr-perturb rather than the refinement spec";
            return false;
        }
        if (t == "acceptdist") {
            // "acceptdist [temp T] [shape S] [anneal] [floor F]" -- the same
            // rule --accept-dist configures, spellable inside a refine spec
            // so a single string can carry the whole search configuration.
            opt.acceptDist.enabled = true;
            while (i + 1 < tok.size()) {
                const string &n = tok[i + 1];
                if (n == "temp" && i + 2 < tok.size()) {
                    opt.acceptDist.temperature = atof(tok[i + 2].c_str());
                    i += 2;
                } else if (n == "shape" && i + 2 < tok.size()) {
                    opt.acceptDist.shape = atof(tok[i + 2].c_str());
                    i += 2;
                } else if (n == "floor" && i + 2 < tok.size()) {
                    opt.acceptDist.tempFloor = atof(tok[i + 2].c_str());
                    i += 2;
                } else if (n == "anneal") {
                    opt.acceptDist.anneal = true;
                    i += 1;
                } else {
                    break;
                }
            }
            if (opt.acceptDist.temperature <= 0.0) {
                err = "\"acceptdist\" temperature must be positive";
                return false;
            }
            if (opt.acceptDist.shape <= 0.0) {
                err = "\"acceptdist\" shape must be positive";
                return false;
            }
            continue;
        }
        if (t == "tunnel") {
            // "tunnel TOL [K]": TOL is required (a tolerance of 0 would
            // make the flag a no-op that still costs probe evaluations),
            // K optional.
            if (i + 1 >= tok.size()) {
                err = "\"tunnel\" needs a tolerance, e.g. \"tunnel 2\" or \"tunnel 2 5\"";
                return false;
            }
            opt.tunnelFlag = true;
            i++;
            opt.tunnelTolerance = atof(tok[i].c_str());
            if (opt.tunnelTolerance <= 0.0) {
                err = "\"tunnel\" tolerance must be positive (got " + tok[i] + ")";
                return false;
            }
            takeOptionalInt(tok, i, opt.tunnelTries);
            if (opt.tunnelTries < 1) {
                err = "\"tunnel\" probe count must be at least 1";
                return false;
            }
            continue;
        }
        if (t == "escape") {
            // both numbers required, same convention as "fullreopt M N"
            opt.escapeFlag = true;
            if (!takeRequiredInt(tok, i, "escape", opt.escapeSpan, err))
                return false;
            if (!takeRequiredInt(tok, i, "escape", opt.escapeTries, err))
                return false;
            continue;
        }
        if (t == "distradius") {
            opt.useDistanceRadius = true;
            continue;
        }
        if (t == "weightprune") {
            // bare "weightprune" keeps its original meaning (bias toward
            // LONG edges); an explicit "long"/"short" token after it picks
            // the direction. See PruneWeighting.
            opt.weightpruneFlag = PRUNE_LONG;
            if (i + 1 < tok.size() && (tok[i + 1] == "long" || tok[i + 1] == "short")) {
                opt.weightpruneFlag = (tok[i + 1] == "short") ? PRUNE_SHORT : PRUNE_LONG;
                i++;
            }
            if (opt.weightpruneFlag == PRUNE_SHORT)
                cout << "WARNING: 'weightprune short' is deprecated -- it measured markedly WORSE"
                        " than both an unweighted prune and 'weightprune long' (see PruneWeighting"
                        " in search_experiments/sprsearch.h). Use 'weightprune long', or drop the flag." << endl;
            continue;
        }

        err = "unknown flag '" + t + "' (see search_experiments/spr_topology_test_usage.txt for the"
                " full vocabulary --spr-refine shares with --hillclimb)";
        return false;
    }

    // same incompatibility --hillclimb's own parser enforces: both replace
    // the step's radius, by different and irreconcilable mechanisms
    if (opt.shrinkFlag && opt.learnradiusFlag) {
        err = "'shrink' and 'learnradius' both replace the step radius; give only one";
        return false;
    }
    return true;
}

string buildRefineRecordTag(const SPRSearchOptions &opt, bool sprMode) {
    if (!sprMode) {
        // an NNI-refined run has no SPR step flags to encode, so its tag is
        // just the refiner plus whichever reporting flags are on
        string tag = "_nni";
        if (opt.findoptEveryNSteps > 0)
            tag += "_findopt";
        return tag;
    }
    // buildRecordTag's own encoding, minus "notree" (always on here, so it
    // would say nothing) and with the refiner named up front
    return "_spr" + string(opt.escapeFlag ? "_escape" : "") + buildRecordTag(opt.useFastSelection, opt.useDistanceRadius, opt.reoptimizeBranchLengths,
            opt.fullReoptEveryNSteps, opt.investigateFlag, opt.investigateRadius, opt.alternateFlag,
            opt.shrinkFlag, opt.learnradiusFlag, opt.sweepFlag, opt.sweepCount, opt.findoptEveryNSteps,
            false, opt.weightpruneFlag != PRUNE_UNIFORM, opt.tunnelFlag, opt.tunnelTolerance,
            opt.slackFlag, opt.slackDelta, opt.slackAnneal)
            + (opt.weightpruneFlag == PRUNE_SHORT ? "short" : "");
}

} // namespace sprsearch
