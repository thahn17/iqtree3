/***************************************************************************
 *   Standalone test/exploration executable for the modern SPR API       *
 *   (isLegalSPR / applySPR / rollbackSPR, declared in phylotree.h).     *
 *                                                                        *
 *   This is a plain executable with its own main() -- it does not run a *
 *   real IQ-TREE analysis pipeline. Most commands are topology-only and *
 *   never touch an alignment or model; the exception is --likelihood,   *
 *   which does load a real alignment to evaluate the tree.              *
 *   Build target: spr_topology_test (see root CMakeLists.txt).          *
 *   Full command reference: tree/spr_topology_test_usage.txt.           *
 *                                                                        *
 *   The self-test's hardcoded start tree mirrors                        *
 *   test_scripts/test_data/spr/six_taxa.start.tree so the fixture file  *
 *   and this test describe the same topology; the test does not read    *
 *   the file itself, to stay runnable from any working directory.       *
 ***************************************************************************/

#include "phylotree.h"
#include "sprsearch.h"
#include "alignment/alignment.h"
#include "model/modelfactory.h"
#include "utils/timeutil.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// This target deliberately does not link the `main` library (it defines
// int main(), which would conflict with the one below). A handful of
// symbols are nonetheless referenced from lower-level libraries this test
// does link (utils/tree/alignment/model) even though they're only ever
// defined in main/*.cpp; this test never exercises the code paths that
// call them, so trivial stubs are enough to satisfy the linker.
void printCopyright(ostream &out) {}
string detectSeqTypeName(string model_name) { return ""; }
void reportRate(ostream &out, PhyloTree &tree) {}
const char *aa_model_names_rax[] = {"LG", "WAG", "JTT", "JTTDCMut", "DCMut", "VT", "PMB", "Blosum62", "Dayhoff",
        "mtREV", "mtART", "mtZOA", "mtMAM",
        "HIVb", "HIVw", "FLU", "rtREV", "cpREV"};

namespace {

// The SPR search machinery this file used to own outright now lives in
// tree/sprsearch.{h,cpp}, so IQTree::doSPRSearch can share it (see that
// header). Pulling the namespace in wholesale keeps every call site
// below spelled exactly as it was before the move.
using namespace sprsearch;

int g_failures = 0;

void expect(bool cond, const string &msg) {
    if (cond) {
        cout << "  [PASS] " << msg << endl;
    } else {
        cout << "  [FAIL] " << msg << endl;
        g_failures++;
    }
}

PhyloNode* requireLeaf(PhyloTree &tree, const string &name) {
    PhyloNode *node = (PhyloNode*) tree.findLeafName(name);
    if (!node) {
        cout << "  [FAIL] could not find leaf '" << name << "' in tree" << endl;
        g_failures++;
    }
    return node;
}

string newickOf(PhyloTree &tree) {
    // limit branch lengths to 1 decimal place so the printed Newick stays
    // readable; this is purely a display setting for this test tool
    Params::getInstance().numeric_precision = 1;
    stringstream ss;
    tree.printTree(ss, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
    return ss.str();
}

typedef unordered_map<PhyloNode*, PhyloNode*> ParentMap;
typedef unordered_map<PhyloNode*, int> DepthMap;

/**
    populate parent/depth maps for every node in the tree, treating
    tree.root (an arbitrary leaf) as the ultimate ancestor. This lets the
    most recent common ancestor (MRCA) of any set of leaves be found, which
    in turn lets an internal edge be addressed from the command line as
    "the edge above the MRCA of these leaves" instead of only ever a
    pendant edge above a single leaf.
 */
void buildAncestry(PhyloNode *node, PhyloNode *dad, int depth, ParentMap &parent, DepthMap &nodeDepth) {
    parent[node] = dad;
    nodeDepth[node] = depth;
    FOR_NEIGHBOR_IT(node, dad, it)
        buildAncestry((PhyloNode*) (*it)->node, node, depth + 1, parent, nodeDepth);
}

PhyloNode* findLCA(PhyloNode *a, PhyloNode *b, ParentMap &parent, DepthMap &nodeDepth) {
    while (nodeDepth[a] > nodeDepth[b])
        a = parent[a];
    while (nodeDepth[b] > nodeDepth[a])
        b = parent[b];
    while (a != b) {
        a = parent[a];
        b = parent[b];
    }
    return a;
}

/**
    split a comma-separated leaf list, trimming whitespace around each name
 */
vector<string> splitLeafList(const string &spec) {
    vector<string> names;
    stringstream ss(spec);
    string item;
    while (getline(ss, item, ',')) {
        size_t start = item.find_first_not_of(" \t");
        size_t stop = item.find_last_not_of(" \t");
        if (start != string::npos)
            names.push_back(item.substr(start, stop - start + 1));
    }
    return names;
}

/**
    resolve an edge spec -- a comma-separated leaf name list -- to the pair
    of currently-adjacent nodes forming that edge:
    - a single leaf name resolves to that leaf's own pendant edge
    - two or more leaf names resolve to the stem edge above their MRCA,
      i.e. a genuinely internal edge
    @return true on success (outNode/outDad set); false on failure (err set)
 */
bool resolveEdgeSpec(PhyloTree &tree, const string &spec, ParentMap &parent, DepthMap &nodeDepth,
        PhyloNode* &outNode, PhyloNode* &outDad, string &err) {
    vector<string> names = splitLeafList(spec);
    if (names.empty()) {
        err = "empty leaf list";
        return false;
    }

    vector<PhyloNode*> leaves;
    for (size_t i = 0; i < names.size(); i++) {
        PhyloNode *leaf = (PhyloNode*) tree.findLeafName(names[i]);
        if (!leaf) {
            err = "leaf '" + names[i] + "' not found in tree";
            return false;
        }
        if (!leaf->isLeaf()) {
            err = "'" + names[i] + "' is not a leaf name in this tree";
            return false;
        }
        leaves.push_back(leaf);
    }

    PhyloNode *node = leaves[0];
    for (size_t i = 1; i < leaves.size(); i++)
        node = findLCA(node, leaves[i], parent, nodeDepth);

    PhyloNode *dad = parent[node];
    if (!dad) {
        err = "'" + spec + "' resolves to the tree's arbitrary root leaf, which has no edge "
              "above it; pick a different leaf set";
        return false;
    }

    outNode = node;
    outDad = dad;
    return true;
}

/**
    collect the names of every leaf in the subtree rooted at `node`, viewed
    away from `dad` -- i.e. the leaf set on the far side of the (node,dad)
    edge from dad. This is the exact inverse of resolveEdgeSpec: feeding the
    comma-joined result back in as an edge spec resolves to this same edge.
 */
} // namespace

/**
    read a tree from `treeArg`, which is either a path to a Newick file or a
    literal Newick string, into `tree`
 */
void readTreeArg(PhyloTree &tree, const string &treeArg) {
    ifstream check(treeArg.c_str());
    if (check.good()) {
        check.close();
        bool is_rooted = false;
        tree.readTree(treeArg.c_str(), is_rooted);
    } else {
        tree.read_TreeString(treeArg, false);
    }
}

/**
    load a tree, prune the edge identified by pruneSpec, and regraft it onto
    the edge identified by regraftSpec, printing the before/after Newick and
    the legality verdict. Does not roll back -- this is a one-shot "what
    would this move do?" check, not a round-trip test.

    Each spec is a comma-separated list of leaf names: a single name
    addresses that leaf's own pendant edge; two or more names address the
    internal edge directly above the MRCA of those leaves (so e.g. "B,D"
    targets the branch leading to the (B,D) clade as a whole, not either
    leaf individually).

    @return 0 if the move was legal and applied, 1 otherwise
 */
int runManualSPR(const string &treeArg, const string &pruneSpec, const string &regraftSpec) {
    PhyloTree tree;
    readTreeArg(tree, treeArg);

    cout << "input tree  : " << newickOf(tree) << endl;

    ParentMap parent;
    DepthMap nodeDepth;
    buildAncestry((PhyloNode*) tree.root, nullptr, 0, parent, nodeDepth);

    PhyloNode *pruneNode, *pruneDad, *regraftNode, *regraftDad;
    string err;

    if (!resolveEdgeSpec(tree, pruneSpec, parent, nodeDepth, pruneNode, pruneDad, err)) {
        cerr << "error resolving prune edge '" << pruneSpec << "': " << err << endl;
        return 2;
    }
    if (!resolveEdgeSpec(tree, regraftSpec, parent, nodeDepth, regraftNode, regraftDad, err)) {
        cerr << "error resolving regraft edge '" << regraftSpec << "': " << err << endl;
        return 2;
    }

    SPRMove move;
    move.prune_node = pruneNode;
    move.prune_dad = pruneDad;
    move.regraft_node = regraftNode;
    move.regraft_dad = regraftDad;
    move.radius = 0;
    move.screening_score = 0.0;
    move.exact_score = 0.0;
    move.candidate_id = 0;
    move.generation = 0;

    cout << "prune edge  : above {" << pruneSpec << "}"
         << (pruneNode->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;
    cout << "regraft edge: above {" << regraftSpec << "}"
         << (regraftNode->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;

    if (!tree.isLegalSPR(move)) {
        cout << "RESULT      : ILLEGAL MOVE" << endl;
        cout << "  (common causes: the two edges are the same or already adjacent; the regraft" << endl;
        cout << "   edge lies inside the subtree being pruned; or the far end of the prune edge" << endl;
        cout << "   is not a normal degree-3 internal node, e.g. it's the tree's arbitrary root)" << endl;
        return 2;
    }

    SPRRollback rollback;
    tree.applySPR(move, rollback);

    cout << "result tree : " << newickOf(tree) << endl;
    return 2;
}

/**
    load a tree, prune the edge identified by pruneSpec, and print every
    distinct, legal graft position within `radius` hops of the prune point
    (see findGraftPositions for exactly how radius/legality/dedup work).
    Does not mutate or apply anything -- this only enumerates candidates.
    @return 0 on success (even if zero candidates were found), 1 on error
 */
int runListGrafts(const string &treeArg, const string &pruneSpec, int radius) {
    PhyloTree tree;
    readTreeArg(tree, treeArg);

    cout << "input tree : " << newickOf(tree) << endl;

    ParentMap parent;
    DepthMap nodeDepth;
    buildAncestry((PhyloNode*) tree.root, nullptr, 0, parent, nodeDepth);

    PhyloNode *pruneNode, *pruneDad;
    string err;
    if (!resolveEdgeSpec(tree, pruneSpec, parent, nodeDepth, pruneNode, pruneDad, err)) {
        cerr << "error resolving prune edge '" << pruneSpec << "': " << err << endl;
        return 2;
    }

    cout << "prune edge : above {" << pruneSpec << "}"
         << (pruneNode->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;
    cout << "radius     : " << radius << endl;
    cout << endl;

    vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, radius);

    if (candidates.empty()) {
        cout << "no legal graft positions found within radius " << radius << endl;
        return 2;
    }

    cout << "found " << candidates.size() << " legal graft position(s):" << endl;
    for (size_t i = 0; i < candidates.size(); i++) {
        const GraftCandidate &c = candidates[i];
        cout << "  [" << (i + 1) << "] radius " << c.radius << ": above {" << describeEdge(c.node, c.dad) << "}"
             << (c.node->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;
    }
    return 2;
}

/**
    load a tree and a DNA or protein alignment (sequence type
    auto-detected from content -- see modelNameFor's comment), and evaluate
    the log-likelihood of that exact tree (topology and branch lengths as
    given) against that alignment under a plain fixed-parameter model (JC
    for DNA, LG for protein) with no rate heterogeneity. Does NOT optimize
    branch lengths or model parameters -- this reports the likelihood of
    the tree exactly as given, not the best achievable likelihood for that
    topology. The alignment's sequence names must match the tree's leaf
    names exactly (case-sensitive).
    @return 0 on success, 1 if the tree/alignment couldn't be read or
    matched, 2 if the alignment's sequence type isn't DNA or protein
 */
int runLikelihood(const string &treeArg, const string &alignmentFile) {
    ifstream check(alignmentFile.c_str());
    if (!check.good()) {
        cerr << "error: cannot open alignment file '" << alignmentFile << "'" << endl;
        return 2;
    }
    check.close();

    // Params is a process-wide singleton read by Alignment/ModelFactory/
    // PhyloTree internals; nothing else in this executable touches it, so
    // it's safe to just reset it to library defaults here.
    Params &params = Params::getInstance();
    params.setDefault();

    InputType intype;
    // nullptr sequence_type (rather than a hardcoded "DNA") makes the
    // Alignment constructor auto-detect DNA vs. protein vs. other types
    // from raw character content (Alignment::detectSequenceType, called
    // unconditionally inside buildPattern before any override is applied)
    Alignment *aln = new Alignment((char*) alignmentFile.c_str(), nullptr, intype, "");
    if (!isSupportedSeqType(aln->seq_type)) {
        cerr << "error: alignment '" << alignmentFile << "' auto-detected as a sequence type this tool "
                "doesn't support (only DNA and protein are handled)" << endl;
        delete aln;
        return 2;
    }

    PhyloTree tree;
    tree.setParams(&params);
    // tree structure must exist before setAlignment, which looks up each
    // alignment sequence name against the tree's already-built leaf nodes
    readTreeArg(tree, treeArg);
    tree.setAlignment(aln);
    // PhyloTree::init() deliberately leaves num_threads at 0, expecting the
    // normal analysis pipeline (which this tool bypasses) to set a real
    // value later; getBufferPartialLhSize() asserts num_threads > 0
    tree.setNumThreads(1);
    // init() also calls setLikelihoodKernel() before any alignment exists,
    // which leaves computeLikelihoodBranchPointer (and friends) null on
    // this build (see the "no alignment specified yet" branch in
    // PhyloTree::setLikelihoodKernel, phylotreesse.cpp); re-run it now
    // that aln is set so it picks the real SSE likelihood kernel instead
    tree.setLikelihoodKernel(LK_SSE2);

    string modelName = modelNameFor(aln->seq_type, false);
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);

    tree.initializeAllPartialLh();
    double logl = tree.computeLikelihood();

    cout << "tree          : " << newickOf(tree) << endl;
    cout << "alignment     : " << alignmentFile << " (" << aln->getNSeq() << " sequences, "
         << aln->getNSite() << " sites, " << (aln->seq_type == SEQ_PROTEIN ? "protein" : "DNA")
         << ", auto-detected)" << endl;
    cout << "model         : " << modelName << " (fixed, no branch length or parameter optimization)" << endl;
    cout << "log-likelihood: " << logl << endl;

    delete aln;
    return 2;
}

/**
    replicate the SHAPE of real IQ-TREE's own starting-tree preprocessing
    (see IQTree::computeInitialTree/initCandidateTreeSet, tree/iqtree.cpp)
    as this standalone tool's own "iqtreestart" --hillclimb starting-tree
    option, instead of the tool's long-standing plain BioNJ estimate:
    generate `poolSize` independent randomized-stepwise-addition parsimony
    trees (PhyloTree::computeParsimonyTree -- IQ-TREE's own native,
    non-PLL parsimony kernel; real IQ-TREE's actual DEFAULT start-tree
    method is PLL's own randomized stepwise-addition parsimony kernel
    instead, which needs a separate PLL instance this standalone tool has
    no other reason to ever set up -- so this is genuinely the SAME KIND
    of preprocessing real IQ-TREE does, and one of its own selectable
    methods ("-t PARS"), but not byte-for-byte its true default), fit each
    one's branch lengths AND (under "gtr") its GTR+FO rate/frequency
    parameters via ModelFactory::optimizeParameters, and return the single
    best-scoring candidate's full-precision Newick string.

    Deliberately stops there rather than also topology-refining the
    winner itself: runHillClimb's own step loop -- already governed by
    <radius>, useFastSelection/numCandidates, reoptimizeBranchLengths,
    and every other --hillclimb flag the caller actually passed -- takes
    over immediately afterward on whatever tree this function returns,
    the same way it already does starting from a plain BioNJ or random
    tree. An earlier version of this function also NNI-hill-climbed the
    pool's best candidates before returning one, via a second,
    separately-parameterized search (fixed radius 1, no fast/shrink/
    alternate/reopt awareness) -- that duplicated runHillClimb's own step
    loop while ignoring the flags the user actually gave, so it was
    removed: whatever refinement iqtreestart's candidate needs now comes
    entirely from the same step loop, controlled by the same parameters,
    as every other starting tree this tool builds.

    Scores every candidate under `modelName` -- so a "gtr" run of
    --hillclimb evaluates this whole pool under GTR+FO, the same model the
    main search itself goes on to use, not JC regardless of that flag.

    Deliberately smaller-scale than real IQ-TREE's own default (roughly
    100 pool trees): this is meant to hand runHillClimb's OWN SPR search a
    better-informed starting point than plain BioNJ, not to reproduce
    IQ-TREE's own tree search wholesale on top of it, so poolSize defaults
    small enough that this preprocessing step doesn't dominate a typical
    --hillclimb run's own cost.

    Every intermediate tree is round-tripped through a FULL-precision
    Newick string (params.numeric_precision temporarily raised to 15),
    never newickOf() -- see maybeRunFindopt's own comment on why: newickOf
    deliberately rounds every branch length to 1 decimal place for
    readable terminal output, which would silently corrupt every
    candidate's branch lengths across this whole pipeline.

    aln->orderPatternByNumChars(PAT_VARIANT) is called once up front if not
    already done -- computeParsimonyTree relies on it (see
    IQTree::computeInitialTree's own identical guard immediately before
    its own first call), and skipping it was a real, previously-shipped
    bug elsewhere in this codebase's history (computeParsimony() silently
    returning 0 for every tree -- see the "Parsimony prescreen" entry in
    this tool's usage doc's "Retired experiments" section).

    The model's own rate/frequency parameters (under "gtr") are fit ONCE,
    on the first pool candidate, then held FIXED for every other
    candidate, which is only ever branch-length-reoptimized
    (PhyloTree::optimizeAllBranches) against that shared, already-fit
    model -- this is the same shape real IQ-TREE's own
    IQTree::initCandidateTreeSet uses (fit via ModelFinder once, then loop
    readTreeString()+optimizeBranches() over the whole candidate pool on
    one reused tree object, tree/iqtree.cpp), replacing an earlier version
    of this function that instead called
    ModelFactory::optimizeParameters() -- a full model+branch-length
    refit -- independently on every candidate. That earlier approach was
    both far slower (poolSize independent GTR+FO fits instead of one) and
    less accurate as a pool ranking, since each candidate's score partly
    reflected how well ITS OWN independently-refit model happened to land
    rather than a shared, consistently-fit yardstick.

    This requires reusing a SINGLE PhyloTree/ModelFactory pair across the
    whole pool (readTreeString() to rebind each new topology onto it)
    rather than one fresh pair per candidate (as initClonedTree/
    runBranchLengthCompare's independent copies do): PhyloTree's own
    destructor unconditionally deletes its model/site_rate/model_factory,
    so handing the same ModelFactory to multiple independently-destructed
    PhyloTree objects would double-free it.
 */
string buildIQTreeStyleStartTree(Alignment *aln, Params &params, const string &modelName, int poolSize) {
    if (aln->ordered_pattern.empty())
        aln->orderPatternByNumChars(PAT_VARIANT);

    // one scratch tree, reused across every computeParsimonyTree() call --
    // each call rebuilds its topology from scratch via randomized stepwise
    // addition (see PhyloTree::computeParsimonyTree's own comment), the
    // same reuse-one-scratch-tree-in-a-loop pattern real IQ-TREE's own
    // initCandidateTreeSet uses (tree/iqtree.cpp) -- so each candidate's
    // Newick must be captured before the next call overwrites it
    PhyloTree parsScratch;
    parsScratch.setParams(&params);
    parsScratch.setNumThreads(1);
    parsScratch.setParsimonyKernel(LK_SSE2);

    int savedPrecision = params.numeric_precision;

    vector<string> poolNewicks(poolSize);
    for (int i = 0; i < poolSize; i++) {
        parsScratch.computeParsimonyTree(nullptr, aln, randstream);
        params.numeric_precision = 15;
        ostringstream fullPrecisionNewick;
        parsScratch.printTree(fullPrecisionNewick, WT_BR_LEN);
        params.numeric_precision = savedPrecision;
        poolNewicks[i] = fullPrecisionNewick.str();
    }

    // same temporary lk_safe_scaling window maybeRunFindopt's own scratch
    // refit uses -- branch-length optimization can legitimately push a
    // branch length toward an extreme value, which the plain
    // (non-scaled) kernel isn't built to handle without a fatal numerical
    // underflow (see reoptimizeSPREdges' comment on runHillClimb for the
    // same concern elsewhere)
    bool origSafeScaling = params.lk_safe_scaling;
    params.lk_safe_scaling = true;

    // candidate #0 both seeds the shared model/rate fit AND is itself the
    // first scored pool member, exactly like real IQ-TREE fitting its
    // model on whatever tree computeInitialTree hands it before the
    // per-candidate loop starts
    PhyloTree tree;
    initClonedTree(tree, poolNewicks[0], aln, params, modelName);
    clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
    // 10x looser than params.modelEps -- matches real IQ-TREE's own
    // one-time initial fit (main/phyloanalysis.cpp: initEpsilon =
    // params.min_iterations == 0 ? params.modelEps : params.modelEps*10,
    // called via IQTree::ensureModelParametersAreSet before its own
    // candidate-set loop starts). Real IQ-TREE only takes the tighter,
    // un-multiplied branch when min_iterations==0, i.e. an explicit "-n 0"
    // (no tree search at all) run; this tool's own maxSteps is a
    // different, unrelated knob (governs runHillClimb's OWN SPR loop, not
    // any NNI iteration count here), so there's no equivalent "-n 0" case
    // to special-case -- the looser epsilon is what a normal run uses
    double bestScore = tree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, params.modelEps * 10);
    params.numeric_precision = 15;
    ostringstream firstPolished;
    tree.printTree(firstPolished, WT_BR_LEN);
    params.numeric_precision = savedPrecision;
    string bestNewick = firstPolished.str();

    for (int i = 1; i < poolSize; i++) {
        // read_TreeString(..., false) -- not the no-argument
        // readTreeString() PhyloTree also has -- to match candidate #0's
        // own initClonedTree call above: readTreeString() auto-detects
        // rootedness from each candidate's own Newick text and updates
        // this->rooted in place, but read_TreeString's is_rooted argument
        // is a local copy that's never written back to this->rooted (see
        // MTree::read_TreeString, tree/mtree.cpp) -- forcing false every
        // time here keeps this->rooted permanently false across the whole
        // loop, exactly as it already implicitly was for candidate #0
        // (never having been set to anything else since this PhyloTree's
        // construction). Without this, a candidate whose own Newick text
        // happens to parse as bifurcating-at-the-root can flip
        // this->rooted to true, which makes PhyloTree::setRootNode take
        // its early "already rooted" return instead of re-deriving root
        // from the alignment's own first sequence name -- leaving root
        // pointed at a node freeNode() already freed, so the very next
        // setAlignment()'s findLeafName() traversal starts from a stale
        // pointer and silently can't reach an entire clade.
        tree.read_TreeString(poolNewicks[i], false);
        tree.setAlignment(aln);
        // full buffer reset, not just clearAllPartialLH(): each pool
        // candidate is an ENTIRELY different topology (not a local SPR/NNI
        // move), so the incremental per-edge buffer-reuse bookkeeping
        // clearAllPartialLH() relies on can't be trusted to still match
        // reality -- same reasoning as every other resetLikelihoodBuffers()
        // call site in this file, see its own comment
        resetLikelihoodBuffers(tree);
        clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
        // params.brlen_num_traversal (default 1), not optimizeAllBranches'
        // own default of 100 -- matches real IQ-TREE's own per-candidate
        // ranking pass (IQTree::optimizeBranches, tree/iqtree.cpp, called
        // as optimizeBranches(params->brlen_num_traversal) from
        // initCandidateTreeSet): candidates only need ONE cheap traversal
        // to be ranked against each other, not near-convergence -- the
        // pool's whole point is picking a starting point for runHillClimb's
        // OWN subsequent search, not delivering a final branch-length fit
        double score = tree.optimizeAllBranches(params.brlen_num_traversal);

        if (score > bestScore) {
            bestScore = score;
            params.numeric_precision = 15;
            ostringstream polishedNewick;
            tree.printTree(polishedNewick, WT_BR_LEN);
            params.numeric_precision = savedPrecision;
            bestNewick = polishedNewick.str();
        }
    }

    params.lk_safe_scaling = origSafeScaling;

    return bestNewick;
}

/**
    greedy, randomized SPR hill-climbing search.

    Takes the tree AliSim used to simulate an alignment (see the AliSim
    command in tree/spr_topology_test_usage.txt for how to produce this
    pair); the alignment file is derived automatically from trueTreeArg by
    AliSim's own default naming convention (<prefix>.treefile paired with
    <prefix>.fa).

    If noTrueTree is set, there's no ground-truth tree at all: trueTreeArg
    is instead a real alignment file path directly (any file/sequence type
    Alignment can auto-detect -- FASTA/NEXUS/PHYLIP, DNA/protein), and the
    final RF-distance-to-true-tree and true-tree-logL reference output are
    both skipped entirely (see every noTrueTree-guarded block further down
    this function). Everything else -- starting tree method, radius,
    selection, model, all other flags -- works exactly the same either way.

    Builds a BioNJ tree from that alignment as the starting "estimate" tree
    (or, if randomStart is true, a random Yule-Harding topology over the
    same taxa instead -- see randomStart below), then repeats up to
    maxSteps times. Each step first picks a prune edge via choosePrune() --
    an O(1) pick from the EdgeRegistry built (and incrementally kept in
    sync) alongside the tree, no traversal needed -- then picks a regraft
    target one of two ways, depending on useFastSelection:

    useFastSelection = false (default, the original --hillclimb behavior):
      1. enumerate every legal regraft candidate within the step's current
         radius (findGraftPositions) -- see shrinkFlag below for how that
         can narrow from `radius` as steps progress
      2. score every candidate by applySPR + computeLikelihood +
         rollbackSPR on the SAME tree object -- no candidate ever gets its
         own copy of the tree
      3. apply the single best-scoring candidate (still via applySPR, on
         that same tree object)
      4. if its likelihood beats the current tree, keep the move; otherwise
         roll it back
    This is a steepest-descent-within-radius search: exhaustive but O(radius)
    candidates scored per step.

    useFastSelection = true: instead of enumerating and scoring every
    candidate in the radius, chooseGraft() picks a single random regraft
    target via a weighted-distance random walk (O(distance walked), no
    enumeration at all -- see chooseGraft's own comment), which is applied,
    scored, and kept or rolled back exactly like the single best candidate
    above. This turns the search into a proposal-based random walk (accept
    if it improves, revert otherwise) rather than steepest descent, trading
    the guarantee of finding the best move in the radius for an O(1)-ish
    per-step cost instead of O(radius).

    Either way, the search also stops early if a chosen prune edge has zero
    legal candidates within the radius (useFastSelection=false) or no edge
    to graft onto at all (useFastSelection=true, only possible on a 3-leaf
    tree), or if the tree runs out of degree-3 nodes to prune from (only
    possible on very small trees).

    randomStart (default false, i.e. the original BioNJ-based --hillclimb
    behavior) replaces the BioNJ estimate tree with a random Yule-Harding
    topology over the same taxa (same engine PhyloTree::generateRandomTree
    uses elsewhere in IQ-TREE for e.g. "-t RANDOM{yh/N}" input trees), so
    the search starts from an arbitrary rather than a distance-based
    starting point. Branch lengths on this random topology come from
    generateRandomTree's own random assignment, not from the alignment.

    iqtreeStart (default false; mutually exclusive with randomStart --
    parseHillClimbFlags rejects giving both) replaces the BioNJ estimate
    tree with buildIQTreeStyleStartTree's own result instead: a pool of
    iqtreeStartPoolSize (default 20, "iqtreestart N" to change it)
    independent randomized-stepwise-addition parsimony trees
    (PhyloTree::computeParsimonyTree), each given a light branch-length
    polish and scored, with the single best-scoring one returned as the
    starting tree. This mirrors the SHAPE of real IQ-TREE's own
    preprocessing (IQTree::computeInitialTree/initCandidateTreeSet,
    tree/iqtree.cpp) -- generate a parsimony pool and pick from it before
    ever handing off to the real search -- at a deliberately smaller scale
    (real IQ-TREE's own default is closer to 100 pool trees; see
    buildIQTreeStyleStartTree's own comment for why this isn't
    byte-for-byte the same method either: real IQ-TREE's actual default
    start-tree method is PLL's own parsimony kernel, not IQ-TREE's native
    one used here, since standing up a separate PLL instance for this
    alone isn't worth it in a tool that otherwise never touches PLL).
    Deliberately does NOT also NNI/SPR-refine the pool's winner itself
    before handing it off -- the step loop just below, governed by
    <radius>, useFastSelection/numCandidates, reoptimizeBranchLengths, and
    every other --hillclimb flag actually given, does that refinement,
    exactly the same way it already does starting from a plain BioNJ or
    random tree; iqtreeStart only ever changes where that loop starts
    from, never how it searches from there.
    Scored under `useGtrModel`'s own modelName throughout (GTR+FO when the
    "gtr" flag is also given, JC otherwise) -- the SAME model the main
    search itself goes on to use, so a "gtr iqtreestart" run's starting
    tree is picked by GTR+FO likelihood, not JC's.

    iqtreeStart also unconditionally triggers the same up-front whole-tree
    fit reoptimizeBranchLengths/fullReoptInitialFit trigger below (see that
    condition's own comment), even if neither of those flags was actually
    given: Newick (what buildIQTreeStyleStartTree hands back) can carry a
    topology and branch lengths, but has no way to carry a substitution
    model's rate/frequency parameter VALUES, so `tree`'s own freshly
    constructed ModelFactory always starts at GTR+FO's arbitrary un-fit
    defaults regardless of how well buildIQTreeStyleStartTree's internal
    fit ranked the pool -- without this, curScore would be capped at
    roughly JC-level likelihood under "gtr iqtreestart" no matter how good
    the chosen topology actually is.

    quiet (default false) suppresses the one printed line per step (prune
    edge, graft target, logL, kept/reverted -- or the "no candidates"/
    "stopping" messages in their place). This is not just cosmetic: each
    line ends with endl, which flushes -- with max-steps in the thousands,
    an interactive console that's slow to render that much output (a
    classic Windows console/PowerShell window is far more prone to this
    than Windows Terminal or a Unix-style pipe) can end up stalling on
    those flushes. All timing in this tool is measured via getCPUTime()
    (process CPU time), not wall-clock time, specifically so a console
    stalled waiting on the terminal to render doesn't inflate the reported
    time -- that dead time is spent blocked, not executing, so it's mostly
    excluded already -- but quiet still removes the per-step writes'
    genuine CPU cost (formatting/streaming all that text) entirely; the
    setup header and the final summary (finished/final tree/RF distance/
    time elapsed) are unaffected either way.

    fullReoptEveryNSteps (default 0, disabled) and fullReoptRounds
    (default 100, only meaningful when fullReoptEveryNSteps is actually
    set) together enable a periodic full tree.optimizeAllBranches(
    fullReoptRounds) sweep every fullReoptEveryNSteps SUCCESSFUL (accepted)
    steps -- steps whose candidate was actually kept, not the total number
    of step attempts, so a stretch of rejected candidates between two
    accepted moves never counts toward the interval and never delays it
    (see maybeRunPeriodicFullReopt's own comment for the exact mechanics).
    This is now a fully independent flag from
    reoptimizeBranchLengths -- it used to only be reachable as "reopt"'s
    own optional trailing number, sharing reoptimizeBranchLengths' 4-edge
    per-move reoptimization automatically; the two are now separate
    concerns that compose freely (either alone, both together, or
    neither), since a periodic whole-tree sweep is useful even without
    paying for a 4-edge NR search on every single candidate along the way,
    and vice versa. The idea behind the periodic sweep itself: when
    reoptimizeBranchLengths is ALSO on, its own reoptimizeSPREdges only
    ever touches the (up to) 4 edges around a given SPR move, so every
    OTHER edge's length is whatever it was left at by the initial sweep
    (or an earlier periodic sweep) -- stale with respect to however much
    the tree's topology has shifted since then; when
    reoptimizeBranchLengths is OFF, branch lengths never change between
    periodic sweeps at all except through them. Either way, a periodic
    full sweep catches that drift at the cost of a whole-tree NR pass
    (roughly 2*(numTaxa-3) edges, several times more expensive than a
    single reoptimizeSPREdges call, since every edge needs re-examining,
    not just the ones a single move touched) instead of 4 edges.
    EXPERIMENTAL -- IQ-TREE's own NNI search loop (IQTree::optimizeNNI,
    iqtree.cpp) calls a full (but single-round-only, my_iterations=1)
    optimizeAllBranches() after every batch of applied NNIs, i.e. N=1 at
    NNI's own granularity; since this tool applies one SPR move per step
    rather than a conflict-free batch per round, N=1 here is a much higher
    relative frequency of full sweeps than NNI's own precedent, and larger
    N trades some of that staleness back for speed. fullReoptRounds caps
    each sweep's own internal convergence loop (optimizeAllBranches' own
    my_iterations, which otherwise defaults to up to 100 -- see its
    comment in phylotree.h) independently of that N/frequency question.
    See the informal comparison in spr_topology_test_usage.txt's
    "fullreopt" section for actual numbers.

    fullReoptInitialFit (default false, only meaningful when
    fullReoptEveryNSteps is set) is fullReopt's own optional trailing
    "true" modifier ("fullreopt M N true" on the command line): when set,
    it makes fullReoptEveryNSteps alone (with reoptimizeBranchLengths OFF)
    also trigger the one-time whole-tree ML branch-length fit on the
    starting tree that reoptimizeBranchLengths always triggers regardless
    of this flag (see the curScore setup block above) -- otherwise, with
    reoptimizeBranchLengths off and fullReoptInitialFit left at its
    default, "fullreopt" alone leaves the starting tree exactly as
    BioNJ/generateRandomTree produced it, only touching branch lengths
    once the first periodic checkpoint is reached. Deliberately excluded
    from buildRecordTag's own tag (see its comment) -- like randomStart,
    it only affects the STARTING point, not per-step search mechanics --
    but included in buildRunId, the same as randomStart.

    useGtrModel (default false) switches the search's substitution model
    from a plain JC (this tool's long-standing default) to GTR+FO
    (ML-estimated rate ratios and frequencies) -- relevant here because
    sim.fa is simulated under GTR{2,4,1,1,4,2}+F{0.3,0.2,0.2,0.3} (see this
    tool's usage doc), so JC is a real model misspecification here, not
    just a simplification: JC has no free rate-matrix or frequency
    parameters at all, so there is nothing for a periodic full reopt to
    "keep in sync" on the model side, only branch lengths. When
    useGtrModel is on and reoptimizeBranchLengths and/or
    fullReoptEveryNSteps is also on, the initial whole-tree fit and every
    periodic sweep both call ModelFactory::optimizeParameters() instead of
    plain PhyloTree::optimizeAllBranches() -- jointly re-fitting the
    model's own rate/frequency parameters alongside branch lengths, not
    just branch lengths -- mirroring real IQ-TREE's sNNI search loop
    (IQTree::optimizeModelParameters, called periodically between NNI
    rounds, not every single accepted move). This is the natural
    combination to test whether "rarely reoptimizing" pays off more when
    there is an actual rate matrix (not just branch lengths) that can
    drift out of sync with the shifting topology. See the informal
    comparison in spr_topology_test_usage.txt's "fullreopt" section
    (searched for "useGtrModel") for what was found.

    On completion, prints the Robinson-Foulds distance between the
    original AliSim tree and the final tree to the terminal, and writes
    both trees plus the RF distance to output.txt (repo root, overwritten
    each run).

    recordProgress (default false) appends this run's convergence
    trajectory to record_<modelName><recordTag>.csv (repo root, '+'
    sanitized to '_' in the model portion; recordTag from buildRecordTag,
    e.g. "_fast_reopt" or "_investigate3" -- see appendRecordRow and
    recordSpreadsheetPath) every time curScore actually changes -- an
    accepted step (the `if (improved)` branch below) or a periodic
    full-reopt change (the fullReoptEveryNSteps block below), NOT every
    candidate considered or every step attempted, so the recorded
    trajectory tracks the search's actual progress rather than its full
    per-candidate volume. Each row is (this run's id from buildRunId,
    candidates evaluated so far, CPU-clock seconds since this run
    started, current logL, and trueTreeLogl - current logL -- how far
    curScore still is below the true AliSim tree's own logL under this
    same model); an initial row (0 candidates, ~setup time, starting
    curScore) is written right after setup so the trajectory has a proper
    starting point. Repeated runs using the same search mode append into
    the same file rather than overwriting it, so their trajectories
    accumulate side by side for comparison; runs using a meaningfully
    different search mode (fast/reopt/investigate) land in a
    separate, appropriately-tagged file instead of being mixed in.

    trajectoryFlag (default false, the "trajectory" command-line flag) is
    independent of recordProgress/recordTopology above: it appends the
    tree's current topology (Newick, no branch lengths) to its own
    trajectory_<run-id>.nwk exactly twice per event worth recording --
    once right after the starting tree is built (the post-initial-tree
    topology, before the step loop's first move) and once more after every
    ACCEPTED step, never a reverted one and never a periodic
    fullreopt/findopt refit (neither ever changes the topology). Unlike
    recordTopology, which only ever adds a companion column to
    recordProgress' own CSV rows and therefore requires it, trajectoryFlag
    needs nothing else turned on. See appendTrajectoryTopology's and
    trajectoryTopologyPath's comments for why it gets its own per-run
    filename (keyed off runId) rather than sharing one file across runs
    the way recordProgress/recordTopology do.

    investigateFlag (default false) changes what happens the step
    immediately AFTER any accepted move (whether that move came from
    "fast" or the exhaustive scan): instead of drawing a fresh prune point
    via choosePrune, that next step re-prunes the EXACT SAME (node,dad)
    pair the accepted move used (dad is never deleted by applySPR, only
    repositioned, so it's always still a valid degree-3 prune point) and
    exhaustively scores every legal regraft candidate within
    investigateRadius real hops of there -- findGraftPositions(tree, node,
    dad, investigateRadius), where a real hop is this tool's own "radius"
    convention as findGraftPositions defines it (radius 1 = the nearest
    possible legal target, equivalent to an NNI; see its own comment).
    If the best of those candidates improves, it's kept (exactly like any
    other accepted step) and the step after THAT one investigates again,
    from the new position -- repeating for as long as each investigation
    keeps improving. As soon as an investigation step fails to improve (or
    finds no legal candidate at all), investigation stops and the
    following step goes back to normal choosePrune-based selection, until
    another move -- from investigation or from normal search -- gets
    accepted and investigation begins again on it.

    investigateRadius (default 1) is how many real hops each investigation
    attempt searches -- 1 (the default) means only the nearest possible
    legal targets (NNI-equivalent moves); larger values search further out
    from the just-accepted move's own position, at proportionally more
    cost per investigation attempt (this is a real exhaustive scan, same
    as the non-"fast" branch, just anchored at a remembered prune point
    and a fixed radius instead of a fresh one and stepRadius). Included in
    both the run's own runId (via buildRunId) and, deliberately, in the
    record spreadsheet's own FILENAME (via buildRecordTag) -- unlike most
    other numeric parameters in this tool, which only ever appear in
    run_id -- since investigateRadius can change an investigation step's
    cost and behavior dramatically (1 vs. a much larger value are
    different enough searches that mixing their trajectories in one file
    would be misleading.

    Each investigation attempt still counts as one of <max-steps>, same as
    any other step -- this does not add extra steps beyond what was
    asked for, it only changes what a given step slot spends its effort
    on. reoptimizeBranchLengths still applies exactly as it does for the
    exhaustive scan, since investigation reuses that same scoring path,
    just against investigateRadius and a remembered prune point instead
    of stepRadius and a fresh one.

    alternateFlag (default false) toggles every other step between a
    plain SPR search at <radius> (even step indices: 0, 2, 4, ...) and an
    NNI-equivalent search forced to radius 1 (odd step indices), where
    "radius 1" means findGraftPositions' own nearest-legal-target tier --
    an actual NNI move, per its comment. This substitutes for <radius> (or
    the current step's chooseGraft/findGraftPositions radius argument)
    wherever it's used for candidate generation -- "fast"'s chooseGraft
    draws and the exhaustive scan's findGraftPositions call alike -- so it
    composes with both without needing its own incompatibility check.
    investigateFlag still takes priority
    on a step it's actively investigating (investigateRadius wins there,
    since that's refining one specific already-accepted move, not a
    step-parity toggle); "alternate" only governs steps investigate isn't
    currently overriding. Step lines print "(spr, radius N)" or "(nni)" in
    place of the usual "(radius N)" while this is active, so the two
    interleaved searches stay visually distinguishable in the log.

    shrinkFlag (default false) replaces the step's own radius -- the fixed
    <radius> otherwise used for every step (overridden/ignored whenever
    this is on) -- with a value that only ever narrows, and only in
    response to actual search stagnation: a persistent shrinkCurrentRadius
    (starts at <radius>, floor of 1) that maybeShrinkRadius decrements by 1 once
    enough CONSECUTIVE non-improving steps (shrinkStallThreshold, held
    CONSTANT for the whole run -- see maybeShrinkRadius' own comment for
    why a fixed threshold beat scaling it down late in the run, and the
    fuller discussion of why a stagnation counter was picked over the
    other options considered: an acceptance-rate moving average, a
    diminishing-logL-gain threshold, or just a nonlinear budget-shaped
    schedule with no history input at all) have piled up, then resets the
    counter. Every actually-improving step (from either candidate-selection
    path -- fast or exhaustive) resets the stall count to 0 regardless of
    radius.
    Composes with "alternate"/"investigate" the same way "alternate" itself
    does: it only changes what stepRadius currently IS, which every other
    flag already reads through that same variable.
    EXPERIMENTAL, including shrinkStallThreshold's own default -- picked
    as a starting point, not empirically tuned.

    learnradiusFlag (default false; mutually exclusive with shrinkFlag --
    parseHillClimbFlags rejects giving both, since they're two different
    mechanisms competing to set the SAME stepRadius) replaces the step's own
    radius with one drawn fresh each step from a distribution fit to recent
    search history, rather than either a fixed <radius> or shrinkFlag's
    one-directional narrowing schedule. learnRadiusContinuous (see its own
    comment, just above maybeShrinkRadius' sibling functions
    sampleStandardNormal/sampleStandardGamma) does the actual drawing; this
    paragraph covers the policy, that function's comment covers the
    mechanics.

    learnRadiusWindow (declared just below investigateNext/
    shrinkCurrentRadius, persisting across step-loop iterations the same
    way) holds the last (up to) learnradiusN successfully-ACCEPTED moves'
    own achieved radii. "Successful" means the move improved curScore (the
    same `improved` flag shrinkFlag's own maybeShrinkRadius call uses). A
    step forced to radius 1 by "alternate"'s own NNI-parity toggle feeds the
    window exactly like any other accepted move: unlike investigation
    (next paragraph), "alternate" is just a different way of picking THIS
    step's radius, not a separate refinement phase layered on top of one
    already-accepted move.

    With investigateFlag OFF, each accepted move is its own, immediate,
    length-1 data point, recorded as bestAchievedRadiusForLearning the
    moment it's accepted -- see "distradius" below for what that actually
    is; without "distradius" it's just bestDistance, the same value each
    step's own log line already reports as "d=" (fast) or "distance "
    (exhaustive).

    With investigateFlag ON, an accepted move's own bestDistance/
    bestAchievedRadiusForLearning is NOT what gets recorded: investigation
    re-prunes and re-grafts the SAME (node,dad) pair repeatedly, each
    attempt measured from wherever the PREVIOUS attempt in the chain left
    it, at investigateRadius rather than the step's own learned radius --
    so any one attempt's own achieved radius describes a local,
    investigateRadius-scaled hop, not a radius the next FRESH excursion
    could sensibly be started at. Instead, the whole chain -- the
    initiating normal move plus every subsequent investigate refinement
    that keeps improving -- is tracked as ONE excursion and recorded as ONE
    number once it's fully done: the net displacement (edgeHopDistance)
    from where the excursion STARTED (learnRadiusAnchorA/B, the
    sibling1/sibling2 edge applySPRTracked's own TrackedSPR leaves behind
    at the initiating move's pre-move position -- see applySPRTracked's
    comment) to where it FINALLY settled (learnRadiusFinalNode paired with
    the excursion's own current pruneDad, reconstructing the last
    successful step's own regraft edge -- see
    maybeFinalizeLearnRadiusExcursion's comment for why pruneDad's OWN edge
    to pruneNode would be off by one hop here instead). Both are declared
    alongside learnRadiusWindow and persist the same way. The excursion
    closes -- and gets measured and recorded, via
    maybeFinalizeLearnRadiusExcursion -- the moment an investigation
    attempt fails to improve OR finds no legal candidate at all (i.e. the
    same point investigateFlag's own comment says investigation itself
    stops); every call site that can end a chain (both "no legal
    candidate(s)" skip branches, and the normal reject/rollback path) calls
    it unconditionally, since it no-ops on its own excursionOpen guard
    whenever there's nothing open to close. A chain still open when
    <max-steps> or "no degree-3 node left to prune" ends the loop entirely
    is likewise finalized once, right after the loop, rather than
    discarded.

    "distradius": findGraftPositions/chooseGraft (everything except
    "fast"+"distradius") are always hop-based, so what goes into the window
    there is always an integer hop count, same as before. "fast"+
    "distradius" is different -- chooseGraftByDistance's own bestDistance
    is a walk-STEP count (outHops), which has no fixed relationship to
    radiusPercent's own percentage scale (see its own comment: a handful of
    long branches and many short ones can spend the same budget in wildly
    different step counts), so recording it directly would feed the fitted
    distribution numbers on a completely different scale than
    sampleStartingRadius's own trapezoid (centered on <radius>, itself a
    single-digit-to-low-double-digit percentage) and the radiusPercent
    draws the whole mixture is meant to produce -- learnRadiusWindow would
    fill with values like 20-50 while <radius> itself stays in single
    digits, corrupting the Gamma fit's own mean/variance estimate (and,
    before this function's own anti-collapse safeguards existed, would
    have driven every draw to clamp at whatever ceiling was in force).
    Avoided by recording chooseGraftByDistance's own outPercentUsed instead
    (the equivalent radiusPercent that would have been JUST enough budget
    to reach the accepted candidate) -- see bestAchievedRadiusForLearning
    and outPercentUsed's own comment. The "investigate" excursion path gets
    the same treatment: edgeHopDistance's own outLength (real summed branch
    length along the excursion's path, not just its hop count) is converted
    to a radiusPercent-equivalent by maybeFinalizeLearnRadiusExcursion the
    same way, so the window stays in one consistent unit -- hop counts
    without "distradius", radiusPercent-equivalents with it -- regardless
    of whether "investigate" is also composed in.

    Each step, learnRadiusContinuous(learnRadiusWindow, learnradiusN,
    radius, learnRadiusMaxPath) draws stepRadiusContinuous from a
    two-component mixture: a Gamma distribution fit (method of moments) to
    the window's current contents, weighted by window.size()/learnradiusN,
    mixed with sampleStartingRadius(radius, learnRadiusMaxPath) weighted by
    the complement -- see sampleStartingRadius's own comment for that
    component's shape. With an empty window this is 100%
    sampleStartingRadius -- "initially, the distribution should be uniform
    between 1 and the max radius" -- and the starting-distribution share
    shrinks in direct proportion as successful moves accumulate. Unlike an
    earlier version of this design, that share does NOT reach 0% outright:
    it bottoms out at a small, permanent minimum (see
    learnRadiusContinuous's own comment for kMaxGammaWeight) once
    learnradiusN moves have landed, so a fresh sampleStartingRadius draw --
    and therefore a genuine chance to notice a radius far from wherever the
    fit has converged -- keeps happening for the rest of the run, not just
    while the window is still filling up. From then on the window slides,
    continuously refitting the SAME Gamma distribution over just the most
    recent learnradiusN successes rather than the run's entire history, so
    the learned radius can keep adapting if the search's own favored radius
    drifts over a long run; the Gamma fit's own shape parameter is
    additionally capped (see learnRadiusContinuous's own comment for
    kMaxGammaShape) so it can never itself narrow all the way down to a
    near-point-mass no matter how tightly clustered the window's recent
    data happens to be -- between that cap and kMaxGammaWeight's own floor,
    learnradius's whole distribution is guaranteed to keep moving in
    response to new data for as long as the run continues, never
    permanently locking onto one value. Recalculation is implicit and
    free: learnRadiusContinuous refits from the window's raw contents on
    every call rather than storing separate running shape/scale parameters
    that would need their own incremental-update logic, since learnradiusN
    is always small enough (a handful to a few dozen samples) that
    repeating the fit from scratch each step is negligible next to a
    single step's own likelihood evaluations.

    learnRadiusMaxPath (declared alongside learnRadiusWindow, computed once
    before the step loop since it never changes across a run) is the
    structural ceiling both sampleStartingRadius's own ramp piece and the
    Gamma-fit component are clamped to: edgeRegistry.slots.size() (no
    simple path in the tree can use more hops than the tree has edges --
    the same bound sweepFlag's own exhaustive full-tree search already
    uses) normally, or 100.0 under "distradius" (100% of the tree's own
    total branch length; a budget that size can already reach anywhere, so
    there is nothing meaningful beyond it).

    stepRadiusContinuous stays fractional all the way into
    chooseGraftByDistance's own radiusPercent (effectiveRadiusContinuous,
    the "alternate"-aware view of it) when "distradius" is active --
    exactly the "continuous" requirement learnRadiusContinuous's own
    comment describes, so a draw of e.g. 4.37 is used as 4.37%, not rounded
    down to 4 first. stepRadius (plain int, rounded via llround and
    clamped) is what findGraftPositions/chooseGraft's own int radius
    parameter gets instead, since a fractional hop count has no meaning
    there. The step's own printed label shows stepRadiusContinuous (2
    decimal places) rather than stepRadius specifically when "learnradius"
    and "distradius" are both active, so the log doesn't silently hide the
    fraction actually being used.
    EXPERIMENTAL, including learnradiusN's own default -- picked as a
    starting point, not empirically tuned.

    sweepFlag (default false) adds a POST-PROCESSING phase that runs
    strictly AFTER the step loop above finishes -- with whatever mix of
    <max-steps> steps, fast/investigate/alternate/shrink shaped
    them -- rather than replacing any of a step's own selection logic the
    way it originally did. Unlike every flag above, it therefore never
    competes with them for what a given STEP does, and composes freely with
    all of them (including "investigate", no longer mutually exclusive with
    it).

    Once the loop is done, it ranks every internal edge of the tree (both
    endpoints degree 3 -- an edge to a leaf has no quartet to test) by
    computeSiblingCompatibilityScore: how well that edge's OWN current
    grouping (its "siblings" -- p's two other neighbors, on one side, q's
    two other neighbors on the other) is supported relative to the two
    single-NNI alternative regroupings around that same edge, following
    adjacent_subtree_compatibility.pdf's section 8 ("when compatibility
    means 'should be siblings'": S_AB = ell(AB|CD) -
    logmeanexp{ell(AC|BD), ell(AD|BC)}, computed here via real trial
    SPR/NNI moves and the tool's own scoreTrialSPRMove, not the note's more
    general per-pair statistic C_bar from section 4 -- see
    computeSiblingCompatibilityScore's own comment for exactly why). A low
    (very negative) score means that edge's current grouping is poorly
    supported -- these are the LEAST-compatible "siblings".

    sweepCount (default 10) is how many of the least-compatible sibling
    pairs, worst-first, get targeted. Each targeted position is the exact
    (node,dad) prune pair computeSiblingCompatibilityScore itself already
    scored (see its own comment for why reusing that SAME pair, rather than
    re-deriving one from the edge some other way, is essential -- otherwise
    the search below isn't guaranteed to even consider the two alternatives
    that made the edge score low to begin with) and searched exhaustively
    across the WHOLE tree -- findGraftPositions(tree, node, dad,
    edgeRegistry.slots.size()), a radius guaranteed to reach every legal
    target since no simple path in the tree can use more hops than the tree
    has edges -- keeping the best regraft if it improves, otherwise
    reverting, exactly like the original one-pass-per-edge design did per
    edge (reoptimizeBranchLengths still applies here, the same way). The
    ranking is a FIXED list, computed once from the tree
    as the step loop left it, not recomputed after each sweep move
    improves -- so, same as the original design, this is NOT guaranteed to
    reach the SPR-optimal tree (an earlier-processed position's move can
    change what the later ones would have scored, or even disturb one of
    the later positions' own edges entirely -- silently skipped if so, see
    the source). It trades that completeness guarantee for a fixed,
    bounded-size, targeted cost: roughly sweepCount * (2n) candidate
    evaluations for the exhaustive searches themselves (compare to the
    original design's roughly (2n)^2, a full pass over every edge), plus a
    cheap O(n) pre-pass (2 trial NNI evaluations per internal edge) to rank
    them. May optionally take a positive integer, e.g. "sweep 20" -- see
    parseHillClimbFlags' comment for the parsing convention. EXPERIMENTAL,
    including sweepCount's own default -- picked as a starting point, not
    empirically tuned.

    findoptFlag (default false) and findoptEveryNSteps (default 0, meaning
    "use the total number of steps", i.e. maxSteps -- resolved right below,
    since maxSteps isn't known yet when parseHillClimbFlags itself runs)
    together enable "findopt": every findoptEveryNSteps loop iterations (by
    step index, the same (step+1) % N idiom fullReoptEveryNSteps uses), run
    ONE whole-tree ML refit on a scratch clone of the current tree, ALWAYS
    under GTR+FO regardless of the main search's own useGtrModel setting
    (ModelFactory::optimizeParameters, jointly reoptimizing branch lengths
    and the model's rate/frequency parameters together) -- PURELY to see
    what logL the richest available model and a full refit would reach
    from here. Unlike every other reopt-family flag in this tool
    (reoptimizeSPREdges, fullReoptEveryNSteps' own periodic sweep, and this
    flag's own previous life as the single-shot, curScore-mutating
    "finalreopt"), findopt makes NO lasting change whatsoever: it never
    touches the real tree, model, or curScore at all, working entirely on a
    throwaway scratch clone that's discarded once its logL is read off --
    see maybeRunFindopt's own comment for the full mechanics. The real
    search that continues after a findopt check point
    is therefore completely unaffected by it; it exists purely to answer
    "how far is curScore from what this exact topology could achieve", read
    off the "findopt: ..." line it prints (or, with recordProgress, its own
    row in the CSV -- see buildRecordTag's comment for why this still earns
    its own file despite no longer mutating curScore the way "finalreopt"
    used to). Because a findopt check point does real, non-search work,
    its CPU-clock cost is deliberately excluded from the run's own timing
    ("pausing the timer" around it) -- see maybeRunFindopt's comment for
    exactly how.

    All timing in this function (cpuClockStart, every recorded CSV row's
    time_elapsed, and the final "time elapsed" summary) is measured via
    getCPUTime() -- process CPU time, not wall-clock -- so repeated runs
    on the same device give consistent readings regardless of unrelated
    system load, background processes, or a slow-to-render console (see
    quiet's own comment for that last one specifically).

    @return 0 on success, 1 if the tree/alignment couldn't be read
 */
int runHillClimb(const string &trueTreeArg, int radius, int maxSteps,
        bool randomStart = false, bool useFastSelection = false, bool quiet = false,
        int numCandidates = 1, bool reoptimizeBranchLengths = false, int fullReoptEveryNSteps = 0,
        int fullReoptRounds = 100, bool fullReoptInitialFit = false, bool useGtrModel = false,
        bool recordProgress = false, bool recordTopology = false,
        bool investigateFlag = false, int investigateRadius = 1, bool alternateFlag = false,
        bool shrinkFlag = false, int shrinkStallThreshold = 10,
        bool learnradiusFlag = false, int learnradiusN = 20,
        bool sweepFlag = false, int sweepCount = 10,
        bool findoptFlag = false, int findoptEveryNSteps = 0,
        bool iqtreeStart = false, int iqtreeStartPoolSize = 20, bool noTrueTree = false,
        bool useDistanceRadius = false, PruneWeighting weightpruneFlag = PRUNE_UNIFORM, bool trajectoryFlag = false,
        bool userStartTree = false, const string &startTreePath = "",
        bool useModelOverride = false, const string &modelOverrideSpec = "") {
    double cpuClockStart = getCPUTime();
    if (findoptFlag && findoptEveryNSteps <= 0)
        findoptEveryNSteps = maxSteps; // "default = total number of steps"

    // AliSim's own default output naming: <prefix>.treefile + <prefix>.fa --
    // UNLESS noTrueTree, in which case there is no ground-truth tree at
    // all and trueTreeArg IS the alignment path directly (any real
    // alignment, not just an AliSim-simulated one; sequence type and file
    // format are both auto-detected below, same as every other command in
    // this file -- see modelNameFor's comment for sequence type, and
    // Alignment's own InputType detection for file format, e.g. NEXUS vs.
    // FASTA vs. PHYLIP).
    string alnFile = trueTreeArg;
    if (!noTrueTree) {
        const string suffix = ".treefile";
        if (alnFile.size() > suffix.size()
                && alnFile.compare(alnFile.size() - suffix.size(), suffix.size(), suffix) == 0)
            alnFile = alnFile.substr(0, alnFile.size() - suffix.size()) + ".fa";
        else
            alnFile += ".fa";
    }

    ifstream alnCheck(alnFile.c_str());
    if (!alnCheck.good()) {
        cerr << "error: could not find alignment '" << alnFile << "'" << endl;
        if (!noTrueTree)
            cerr << "  (derived from the tree argument by replacing '.treefile' with '.fa',"
                    " AliSim's own default output naming; generate a pair with e.g." << endl
                 << "   iqtree3 --alisim <prefix> -m \"GTR{2,4,1,1,4,2}+F{0.3,0.2,0.2,0.3}\""
                    " -t \"RANDOM{yh/100}\" --length 10000)" << endl;
        return 2;
    }
    alnCheck.close();

    // Params is a process-wide singleton read by Alignment/ModelFactory/
    // PhyloTree internals; nothing else in this executable touches it, so
    // it's safe to just reset it to library defaults here, before loading
    // the alignment below (moved up from further down in this function,
    // since the alignment now needs to be loaded before the settings
    // summary prints, to report its auto-detected sequence type).
    Params &params = Params::getInstance();
    params.setDefault();
    // computeBioNJ writes/reads temporary files alongside this prefix
    // (<prefix>.mldist, <prefix>.bionj) as part of how it builds the tree
    params.out_prefix = (char*) "spr_hillclimb_tmp";

    // Load the alignment now (before the settings summary below), with its
    // own narrow cout suppression, so the summary can report the
    // auto-detected sequence type and this run's actual model name.
    // Loading itself is fast even for large alignments -- the slow part
    // this tool avoids blocking the terminal on before printing anything
    // is BioNJ/iqtreeStart preprocessing further below, not this read.
    // nullptr sequence_type (rather than a hardcoded "DNA") makes the
    // Alignment constructor auto-detect DNA vs. protein vs. other types
    // from raw character content (Alignment::detectSequenceType, called
    // unconditionally inside buildPattern before any override is applied).
    InputType intype;
    Alignment *aln;
    {
        ostringstream suppressedAlignmentLoad;
        streambuf *realCoutBufAln = cout.rdbuf(suppressedAlignmentLoad.rdbuf());
        aln = new Alignment((char*) alnFile.c_str(), nullptr, intype, "");
        cout.rdbuf(realCoutBufAln);
    }
    if (!isSupportedSeqType(aln->seq_type)) {
        cerr << "error: alignment '" << alnFile << "' auto-detected as a sequence type this tool "
                "doesn't support (only DNA and protein are handled)" << endl;
        delete aln;
        return 2;
    }

    // GTR+FO (ML-estimated rates and frequencies) for DNA, or LG/LG+FO for
    // protein (see modelNameFor's comment) instead of a plain fixed model:
    // sim.fa was simulated under GTR{2,4,1,1,4,2}+F{0.3,0.2,0.2,0.3} (see
    // the AliSim command in this tool's usage doc), so searching a DNA
    // alignment under plain JC -- which has no free rate/frequency
    // parameters to fit at all -- is a deliberately misspecified model, not
    // just a simplification. See useGtrModel's comment on runHillClimb for
    // why this matters specifically for whether periodic re-optimization
    // is worth its cost.
    string modelName = modelNameFor(aln->seq_type, useGtrModel);
    string recordTag = buildRecordTag(useFastSelection, useDistanceRadius, reoptimizeBranchLengths,
            fullReoptEveryNSteps, investigateFlag, investigateRadius, alternateFlag, shrinkFlag, learnradiusFlag,
            sweepFlag, sweepCount, findoptEveryNSteps, noTrueTree, weightpruneFlag != PRUNE_UNIFORM);

    // print this run's own settings BEFORE building the starting tree,
    // not after: constructing it (BioNJ is near-instant, but iqtreeStart's
    // own parsimony-pool preprocessing can take tens of seconds -- see its
    // own comment) used to leave the terminal completely silent until it
    // finished, which looks indistinguishable from a hang on a slow run.
    // Only the "start tree" line itself (needs the actual tree/curScore)
    // and its own timing still have to wait until after that setup below.
    cout << "sequence type   : " << (aln->seq_type == SEQ_PROTEIN ? "protein" : "DNA")
         << " (auto-detected)" << endl;
    if (noTrueTree)
        cout << "true tree       : none -- no RF distance or true-tree logL reference this run" << endl;
    cout << "radius          : " << radius << endl;
    cout << "max steps       : " << maxSteps << endl;
    if (iqtreeStart)
        cout << "iqtreestart     : best of " << iqtreeStartPoolSize << " parsimony pool tree(s)" << endl;
    if (useFastSelection)
        cout << "selection       : fast (choosePrune/chooseGraft proposal, not exhaustive)"
             << (numCandidates > 1 ? ", " + to_string(numCandidates) + " candidates/step" : "") << endl;
    if (useFastSelection && useDistanceRadius)
        cout << "distradius      : radius " << radius << " read as " << radius
             << "% of the tree's total branch length, walked by summed branch"
                " length instead of hop count" << endl;
    if (weightpruneFlag != PRUNE_UNIFORM)
        cout << "weightprune     : prune edge chosen with probability proportional to "
             << (weightpruneFlag == PRUNE_SHORT ? "the RECIPROCAL of its own branch length"
                                                : "its own branch length")
             << " instead of uniformly, experimental" << endl;
    if (reoptimizeBranchLengths)
        cout << "branch lengths  : re-optimized (Newton-Raphson, like NNI) on the 4 edges around "
                "each candidate (the 3 it changes, plus the pruned subtree's own upper edge) "
                "before scoring, experimental" << endl;
    if (fullReoptEveryNSteps > 0)
        cout << "full reopt      : whole-tree "
             << (useGtrModel ? "optimizeParameters() (model + branch lengths)"
                             : "optimizeAllBranches(" + to_string(fullReoptRounds) + " round(s))")
             << " sweep every " << fullReoptEveryNSteps << " successful step(s)"
             << (fullReoptInitialFit ? ", plus one up front on the starting tree" : "")
             << ", experimental" << endl;
    if (useGtrModel)
        cout << "model           : " << modelName << " (ML-estimated rates/frequencies), experimental" << endl;
    if (recordProgress)
        cout << "record          : appending to " << recordSpreadsheetPath(modelName, recordTag)
             << ", experimental" << endl;
    if (investigateFlag)
        cout << "investigate     : refine each accepted move within " << investigateRadius
             << " real hop(s) next step, experimental" << endl;
    if (alternateFlag)
        cout << "alternate       : toggling every other step between SPR (radius " << radius
             << ") and NNI (radius 1), experimental" << endl;
    if (shrinkFlag)
        cout << "shrink          : radius starts at " << radius << ", narrows by 1 (floor 1) after "
             << shrinkStallThreshold << " consecutive non-improving steps, experimental" << endl;
    if (learnradiusFlag)
        cout << "learnradius     : radius drawn each step from a Gamma fit (method-of-moments, capped spread) "
                "to the last " << learnradiusN << " successful move(s)' own radius, mixed with a starting "
                "distribution centered on " << radius << " (flat below it, tapering to zero at "
             << (useDistanceRadius ? "100%" : "the tree's own edge count")
             << ") in proportion to how many of those " << learnradiusN << " slot(s) are filled so far, "
                "never reaching 100% Gamma, experimental" << endl;
    if (sweepFlag)
        cout << "sweep           : after all steps finish, exhaustive whole-tree regraft search on the "
             << sweepCount << " least-compatible sibling pair(s) (adjacent_subtree_compatibility.pdf), "
                "experimental" << endl;
    if (findoptEveryNSteps > 0)
        cout << "findopt         : every " << findoptEveryNSteps << " step(s), scratch whole-tree "
                "optimizeParameters() (richest available model + branch lengths, always, regardless of "
                "'gtr') refit on a rolled-back copy (main tree unaffected), experimental" << endl;

    // the original AliSim tree, kept as a separate plain tree purely for
    // the final RF-distance comparison -- never touched by any SPR move.
    // Left default-constructed/empty when noTrueTree (trueTreeArg is an
    // alignment path in that case, not a tree to read) -- every use of it
    // and trueTreeNewick further down is itself guarded by !noTrueTree.
    PhyloTree trueTree;
    string trueTreeNewick;
    if (!noTrueTree) {
        readTreeArg(trueTree, trueTreeArg);
        trueTreeNewick = newickOf(trueTree);
    }

    // seeded here, before any random tree generation, so that both a
    // randomStart topology and every step's random prune-edge choice come
    // from the same seeded sequence; time(nullptr) alone has only 1-second
    // resolution, which repeats the exact same "random" sequence across
    // rapid successive runs (e.g. a test script invoking this back to
    // back), so getRealTime()'s sub-second precision is mixed in too
    init_random((int) (time(nullptr) * 1000 + (long) (getRealTime() * 1000) % 1000));

    // the distance/BioNJ/model setup below goes through several library
    // code paths (computeDist, computeBioNJ, ModelFactory) that print
    // their own progress noise (composition test, distance matrix,
    // RapidNJ progress, ...) unconditionally on cout; none of it is
    // useful for this test tool, so silence cout for the duration of the
    // setup and restore it before printing our own summary. (The
    // alignment itself was already loaded, with its own narrow
    // suppression, further up -- before the settings summary.)
    ostringstream suppressedSetupOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedSetupOutput.rdbuf());

    // modelName and aln were already computed/loaded above, before this
    // run's settings summary was printed
    PhyloTree tree(aln);
    tree.setParams(&params);
    if (randomStart) {
        // generateRandomTree requires tree.aln (set by the PhyloTree(aln)
        // constructor above) and tree.params (set just above); it builds a
        // random Yule-Harding topology, renames its leaves to match aln's
        // sequence names, and calls setAlignment(aln) on itself internally
        // (see PhyloTree::generateRandomTree / readTreeStringSeqName) --
        // this is the same engine IQ-TREE's own "-t RANDOM{yh/N}" uses
        tree.generateRandomTree(YULE_HARDING);
    } else if (iqtreeStart) {
        // see buildIQTreeStyleStartTree's own comment for the full
        // preprocessing this replaces plain BioNJ with -- the step loop
        // below (radius/useFastSelection/reoptimizeBranchLengths/etc.)
        // does all of this tree's own topology refinement, same as it
        // already does for a plain BioNJ or random start
        string startNewick = buildIQTreeStyleStartTree(aln, params, modelName, iqtreeStartPoolSize);
        tree.read_TreeString(startNewick, false);
        tree.setAlignment(aln);
    } else if (userStartTree) {
        // caller-supplied starting tree (path to a Newick file, or a literal
        // Newick string) -- topology AND branch lengths are trusted exactly
        // as given, e.g. the final tree of a real iqtree3 NNI search already
        // run on this same alignment. No BioNJ/random/iqtreestart tree-build
        // happens on this path, and (see the reoptimizeBranchLengths gate
        // just below) no up-front full branch-length/model re-fit either --
        // this is deliberately the "skip the initial tree find entirely"
        // starting-tree mode
        readTreeArg(tree, startTreePath);
        tree.setAlignment(aln);
    } else {
        // no tree file is read here: computeDist + computeBioNJ build the
        // starting tree structure directly from the alignment's own distance
        // matrix, instead of the readTree-then-setAlignment order used by
        // runLikelihood/runManualSPR. This first pass uses cheap,
        // model-free pairwise distances: PhyloTree::computeDist falls back
        // to Alignment::computeDist whenever no model/rate is attached yet
        // (tree/phylotree.cpp, computeDist(seq1,seq2,...)'s own
        // "if (!model_factory || !site_rate) return initial_dist" guard),
        // which is unconditionally true here since `tree` has no model at
        // all yet.
        tree.computeDist(params, aln, tree.dist_matrix, tree.var_matrix);
        tree.computeBioNJ(params);
        // re-map leaf ids to match the alignment's sequence order/names, same
        // as runLikelihood does after reading a tree from file
        tree.setAlignment(aln);

        // real IQ-TREE's own two-round bootstrap (main/phyloanalysis.cpp:
        // ensureModelParametersAreSet fits a model on exactly this kind of
        // cheap-distance BIONJ tree first, then computeMLDist reuses that
        // fit to compute true ML pairwise distances, then computeBioNJ
        // runs a SECOND time on the improved matrix -- confirmed by
        // tracing a real "-t BIONJ" run's own log, which prints
        // "Constructing BIONJ tree" and "Estimate model parameters" TWICE
        // each). Without this, this tool's own BIONJ topology stayed stuck
        // at whatever the cheap distance formula alone could produce --
        // RF 76 from the true tree on the checked-in sim.treefile/sim.fa,
        // versus real IQ-TREE's own considerably closer BIONJ topology on
        // the identical alignment -- since the model-fit block below only
        // ever touches branch lengths/model parameters on the topology
        // already fixed above, never the topology itself. This block's own
        // ModelFactory is deliberately a throwaway: it exists only to make
        // the SECOND computeDist call below take computeDist's ML-fitting
        // path instead of its cheap fallback; the common setup further
        // below discards it and fits a fresh one against the resulting,
        // better topology -- see the `delete tree.getModelFactory()` just
        // before that fresh construction
        ModelsBlock *prelimModelsBlock = readModelsDefinition(params);
        string prelimModelName = modelName; // ModelFactory's ctor wants string&, not const string&
        tree.setModelFactory(new ModelFactory(params, prelimModelName, &tree, prelimModelsBlock));
        delete prelimModelsBlock;
        tree.setModel(tree.getModelFactory()->model);
        tree.setRate(tree.getModelFactory()->site_rate);
        tree.setNumThreads(1);
        tree.setLikelihoodKernel(LK_SSE2);
        tree.initializeAllPartialLh();
        clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
        bool prelimSafeScaling = params.lk_safe_scaling;
        params.lk_safe_scaling = true;
        if (useGtrModel)
            tree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, params.modelEps * 10);
        else
            tree.optimizeAllBranches(100);
        params.lk_safe_scaling = prelimSafeScaling;

        // now that a model/rate is attached, this recomputes true ML
        // pairwise distances instead of the cheap fallback above (same
        // function, same internal guard -- see this block's own comment)
        tree.computeDist(params, aln, tree.dist_matrix, tree.var_matrix);
        tree.computeBioNJ(params);
        tree.setAlignment(aln);
    }

    // built once, O(n); kept in sync thereafter by applySPRTracked/
    // rollbackSPRTracked, never rebuilt -- see EdgeRegistry's comment
    EdgeRegistry edgeRegistry;
    buildEdgeRegistry(tree, edgeRegistry);

    tree.setNumThreads(1);
    if (reoptimizeBranchLengths || !randomStart)
        // real Newton-Raphson branch-length search (see scoreTrialSPRMove)
        // can legitimately push a branch length toward an extreme value
        // while searching, which the plain (non-scaled) SSE kernel isn't
        // built to handle without numerical underflow in the likelihood
        // derivative -- exactly the scenario IQ-TREE's own "-safe" option
        // exists for; the naive fixed-length scoring path never searches
        // branch lengths at all, so it doesn't need this. Both iqtreeStart
        // and the plain BioNJ default need it too, since both now trigger
        // the up-front full fit below (see that condition's own comment)
        // -- only randomStart still skips it, since it's the only path
        // that doesn't
        params.lk_safe_scaling = true;
    tree.setLikelihoodKernel(LK_SSE2);

    // modelName was already computed above, before the starting-tree branch.
    // delete first: the plain-BioNJ branch above may have left a throwaway
    // preliminary ModelFactory attached to `tree` (used only to make its
    // own second computeDist call ML-fit pairwise distances -- see that
    // branch's own comment); setModelFactory is a plain pointer
    // reassignment with no ownership transfer of whatever it's
    // overwriting, so skipping this would leak it. A no-op for the other
    // two starting-tree branches, which never attach a model this early
    delete tree.getModelFactory();
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    // useModelOverride ('model <spec>') feeds ModelFactory a literal, fully
    // resolved model string (e.g. "LG+F{0.081,0.055,...}", 20 comma-separated
    // amino acid frequencies in the same A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,
    // Y,V order alignment/alignment.cpp's symbols_protein and IQ-TREE's own
    // .iqtree report use) INSTEAD of modelName's own short "LG"/"LG+FO" label
    // -- modelName itself is untouched and keeps naming record/run-id files
    // (see recordSpreadsheetPath/buildRunId below), only the actual model
    // ModelFactory builds changes. A fixed "+F{...}" spec like that (as
    // opposed to "+FO", which would still re-optimize from scratch) means the
    // tree's model starts, and -- for userStartTree with no further
    // optimization -- STAYS at exactly those externally supplied values, e.g.
    // ones a real iqtree3 run already ML-fit on this same alignment, instead
    // of either this tool's own un-fit "LG+FO" starting parameters or having
    // to refit them itself with no such external input
    string modelNameForFactory = useModelOverride ? modelOverrideSpec : modelName;
    tree.setModelFactory(new ModelFactory(params, modelNameForFactory, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.initializeAllPartialLh();

    double curScore;
    // userStartTree skips this whole up-front fit unconditionally -- its
    // whole point is a starting tree whose branch lengths and model are
    // already trustworthy (e.g. real iqtree3 NNI output), so it's treated
    // like randomStart here even though it isn't randomStart itself
    if (!userStartTree
            && (reoptimizeBranchLengths || (fullReoptEveryNSteps > 0 && fullReoptInitialFit) || !randomStart)) {
        // optimizeAllBranches() needs a model/rate already assigned and
        // valid partial-likelihood buffers (both just set up above by
        // initializeAllPartialLh()), since -- unlike a plain length
        // assignment -- it actually evaluates and re-optimizes the
        // likelihood along the way. BioNJ's distance-based lengths (or
        // generateRandomTree's own random assignment) are never otherwise
        // trusted as meaningful in this tool, so rather than let
        // reoptimizeBranchLengths' per-candidate NR search (or, if
        // fullReoptInitialFit asked for it, fullReoptEveryNSteps' own
        // first periodic sweep) start from that arbitrary point, run one
        // full round of ML branch-length optimization -- the same routine
        // IQ-TREE's own search uses to
        // refine a freshly built tree (see IQTree::doTreeSearch's
        // `curScore = optimizeAllBranches(1)`) -- across every edge of the
        // fixed starting topology first. Every edge must be clamped first:
        // optimizeAllBranches calls optimizeOneBranch on literally every
        // edge, and BioNJ can leave a zero/negative estimate on any of
        // them (see clampBranchLengthForOptimization's comment), not just
        // ones next to a particular SPR move.
        //
        // iqtreeStart ALSO needs this unconditionally (not just when
        // reopt/fullreopt-init happen to be given too): buildIQTreeStyleStartTree
        // already fits its OWN internal candidates' model/branch-length
        // parameters to rank the pool fairly, but Newick has no way to
        // carry a substitution model's rate/frequency parameter VALUES
        // across the round-trip into `tree` -- only topology and branch
        // lengths survive it. Without this, `tree`'s own freshly
        // constructed ModelFactory would sit at GTR+FO's arbitrary un-fit
        // starting parameters regardless of how good the chosen pool
        // candidate's fit was, capping curScore at roughly JC-level
        // likelihood no matter how good the topology is (real IQ-TREE
        // never has this problem: it fits its model once via ModelFinder
        // and keeps reusing that SAME live ModelFactory object throughout
        // preprocessing, never round-tripping through Newick in between).
        //
        // The plain BioNJ default (!randomStart, no other flag given)
        // needs it for the identical reason: real IQ-TREE's own "-t BIONJ"
        // path always fits the model on the BioNJ tree before its search
        // ever starts (IQTree::ensureModelParametersAreSet, called from
        // main/phyloanalysis.cpp before doTreeSearch) -- skipping that
        // step here left this tool's own default starting logL
        // catastrophically far from the true tree's (tens of thousands of
        // log-lik units, not the ~100-unit gap a real, unfit-model-free
        // BioNJ start should show), since GTR+FO's un-fit rate/frequency
        // parameters dominate the score far more than the raw BioNJ
        // branch lengths' own imprecision does. NOT extended to
        // randomStart: that path has no distance-based lengths to trust
        // or distrust in the first place (generateRandomTree assigns
        // arbitrary lengths outright), and fixing it is a separate,
        // unasked-for change
        clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
        if (useGtrModel)
            // GTR+FO's rate ratios and frequencies start at arbitrary
            // (roughly uniform) values -- unlike JC, which has no such
            // parameters to fit -- so a plain optimizeAllBranches() call
            // would leave them untouched and meaningless. optimizeParameters
            // internally alternates model-parameter optimization
            // (model->optimizeParameters/site_rate->optimizeParameters)
            // with its own bounded rounds of optimizeAllBranches() until
            // jointly converged -- the same routine real IQ-TREE analyses
            // use to fit a model to data (IQTree::optimizeModelParameters)
            curScore = tree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, params.modelEps);
        else
            // fullReoptEveryNSteps' own fullReoptRounds ceiling applies
            // here too when it's the one that's on (reoptimizeBranchLengths
            // off): 100 (this call's own default when fullReoptEveryNSteps
            // is 0, i.e. only reoptimizeBranchLengths is on) otherwise.
            curScore = tree.optimizeAllBranches(fullReoptEveryNSteps > 0 ? fullReoptRounds : 100);
    } else {
        curScore = tree.computeLikelihood();
    }

    // logL of the original AliSim tree (topology and branch lengths exactly
    // as simulated) against this same alignment/model, purely as a
    // reference point for how the final hill-climbed tree's logL compares --
    // trueTree gets its own ModelFactory since a PhyloTree's model holds a
    // pointer back to that exact tree, so it can't be shared with `tree`.
    // Skipped entirely when noTrueTree: there's no ground-truth topology to
    // fit a model against, so trueTreeLogl stays NaN (self-documenting --
    // every print/record site that uses it just shows/writes NaN, rather
    // than this needing its own separate on/off plumbing at each of them).
    double trueTreeLogl = std::numeric_limits<double>::quiet_NaN();
    if (!noTrueTree) {
        trueTree.setParams(&params);
        trueTree.setAlignment(aln);
        trueTree.setNumThreads(1);
        trueTree.setLikelihoodKernel(LK_SSE2);
        ModelsBlock *trueTreeModelsBlock = readModelsDefinition(params);
        trueTree.setModelFactory(new ModelFactory(params, modelName, &trueTree, trueTreeModelsBlock));
        delete trueTreeModelsBlock;
        trueTree.setModel(trueTree.getModelFactory()->model);
        trueTree.setRate(trueTree.getModelFactory()->site_rate);
        trueTree.initializeAllPartialLh();
        if (useGtrModel)
            // fit GTR+FO's rate/frequency parameters (only) against the TRUE
            // simulated topology and branch lengths (BRLEN_FIX -- those lengths
            // are exactly as AliSim generated them, not to be touched), so this
            // reference logL reflects a properly-fit model rather than
            // GTR+FO's arbitrary un-optimized starting parameters
            trueTreeLogl = trueTree.getModelFactory()->optimizeParameters(BRLEN_FIX, false, params.modelEps);
        else
            trueTreeLogl = trueTree.computeLikelihood();
    }

    // findopt's own scratch refit (see maybeRunFindopt's comment) always
    // fits under the richest available model (GTR+FO for DNA, LG+FO for
    // protein), regardless of the main run's own model -- so comparing its
    // reading against the plain trueTreeLogl above would be unfair whenever
    // the main run ISN'T already under that richer model (useGtrModel ==
    // false, trueTreeLogl there being an unfit JC/LG likelihood with no
    // free parameters to begin with): any gap would then be partly genuine
    // search progress and partly just a free lunch from a strictly richer
    // model. Build a second, richer-model-under-BRLEN_FIX reference off a
    // scratch clone of the SAME true topology and branch lengths (same
    // clone-from-full-precision-Newick approach maybeRunFindopt uses for
    // the main tree, so trueTree's own state -- still needed below for
    // computeRFDist -- is never touched), purely so findopt's own
    // comparisons stay apples-to-apples. Skipped entirely when findopt
    // isn't even in use (findoptEveryNSteps == 0), when noTrueTree (no
    // trueTree to clone from -- findopt's own diagnostic refit still runs,
    // it just has no "gap to true tree" reference to report, same as
    // trueTreeLogl itself), since it'd otherwise be wasted work;
    // trueTreeLogl itself is reused as-is when useGtrModel already made it
    // a richer-model fit in the first place.
    double trueTreeLoglForFindopt = trueTreeLogl;
    if (!noTrueTree && findoptEveryNSteps > 0 && !useGtrModel) {
        int savedPrecisionTrueTree = params.numeric_precision;
        params.numeric_precision = 15;
        ostringstream trueTreeFullPrecisionNewick;
        trueTree.printTree(trueTreeFullPrecisionNewick, WT_BR_LEN);
        params.numeric_precision = savedPrecisionTrueTree;

        PhyloTree trueTreeGtrScratch;
        initClonedTree(trueTreeGtrScratch, trueTreeFullPrecisionNewick.str(), aln, params,
                modelNameFor(aln->seq_type, true));
        trueTreeLoglForFindopt =
            trueTreeGtrScratch.getModelFactory()->optimizeParameters(BRLEN_FIX, false, params.modelEps);
    }

    cout.rdbuf(realCoutBuf);

    // cpuClockStart was set at the very top of this function, before any
    // of the (possibly slow, always cout-suppressed) setup above -- so
    // this is exactly the CPU-clock cost of building/scoring the starting
    // tree, before the step loop below ever runs its first candidate.
    // Every other timing in this tool is CPU-clock-based for the same
    // reason (see runHillClimb's own comment on `quiet`), so this reuses
    // that same clock/baseline rather than starting a separate one.
    double startTreeCpuTime = getCPUTime() - cpuClockStart;

    cout << (randomStart ? "random start tree: " : (iqtreeStart ? "iqtree start tree: "
            : (userStartTree ? "user start tree  : " : "BioNJ start tree : ")))
         << newickOf(tree) << " (logL = " << curScore << ")" << endl;
    cout << "start tree time : " << fixed << setprecision(2) << startTreeCpuTime << " sec (CPU)" << endl;

    // Everything the step loop needs that outlives an individual step now
    // lives in SPRSearchOptions/SPRSearchState (tree/sprsearch.h) rather
    // than in ~40 locals right here, so IQ-TREE's own stochastic search can
    // re-enter the same loop once per perturbation (see IQTree::doSPRSearch).
    // This command still enters it exactly once, with exactly the values
    // these locals always held.
    SPRSearchOptions opt;
    opt.radius = radius;
    opt.useFastSelection = useFastSelection;
    opt.numCandidates = numCandidates;
    opt.useDistanceRadius = useDistanceRadius;
    opt.weightpruneFlag = weightpruneFlag;
    opt.alternateFlag = alternateFlag;
    opt.investigateFlag = investigateFlag;
    opt.investigateRadius = investigateRadius;
    opt.shrinkFlag = shrinkFlag;
    opt.shrinkStallThreshold = shrinkStallThreshold;
    opt.learnradiusFlag = learnradiusFlag;
    opt.learnradiusN = learnradiusN;
    opt.reoptimizeBranchLengths = reoptimizeBranchLengths;
    opt.fullReoptEveryNSteps = fullReoptEveryNSteps;
    opt.fullReoptRounds = fullReoptRounds;
    opt.useGtrModel = useGtrModel;
    opt.sweepFlag = sweepFlag;
    opt.sweepCount = sweepCount;
    opt.findoptFlag = findoptFlag;
    opt.findoptEveryNSteps = findoptEveryNSteps;
    opt.quiet = quiet;
    opt.recordProgress = recordProgress;
    opt.recordTopology = recordTopology;
    opt.trajectoryFlag = trajectoryFlag;
    opt.noTrueTree = noTrueTree;

    string runId = buildRunId(radius, maxSteps, randomStart, useFastSelection,
            numCandidates, reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds, fullReoptInitialFit,
            useGtrModel, investigateFlag, investigateRadius, alternateFlag, shrinkFlag,
            shrinkStallThreshold, learnradiusFlag, learnradiusN, sweepFlag, sweepCount, findoptEveryNSteps,
            iqtreeStart, iqtreeStartPoolSize, weightpruneFlag != PRUNE_UNIFORM);
    SPRSearchState searchState;
    initSPRSearchState(searchState, opt, edgeRegistry, runId, recordTag, modelName,
            cpuClockStart, trueTreeLogl, trueTreeLoglForFindopt);
    if (recordProgress)
        appendRecordRow(modelName, recordTag, runId, searchState.candidatesEvaluated,
                getCPUTime() - cpuClockStart, curScore, trueTreeLogl, recordTopology, tree);
    // "trajectory": captures the post-initial-tree topology -- the tree
    // exactly as it stands right here, after whatever up-front fit
    // (reopt/fullreopt-init/the plain BioNJ-model fit) already ran above,
    // but before the step loop's first SPR move -- as this run's first
    // trajectory_<run-id>.nwk line. See appendTrajectoryTopology's comment.
    if (trajectoryFlag) {
        cout << "trajectory      : appending to " << trajectoryTopologyPath(runId)
             << " (start + every accepted move), experimental" << endl;
        appendTrajectoryTopology(runId, tree);
    }

    // Reusable path-partial scratch buffers. Not used for branch-length
    // reoptimization trials: those deliberately keep the full-reset path.
    SPRLocalLhCache sprLocalCache;
    bool haveSprLocalCache = !reoptimizeBranchLengths;
    if (haveSprLocalCache)
        allocateSPRLocalLhCache(tree, sprLocalCache);

    curScore = runSPRSteps(tree, edgeRegistry, opt, searchState, curScore, maxSteps, aln, params,
            haveSprLocalCache ? &sprLocalCache : nullptr);

    // "sweep": post-processing phase, run strictly AFTER the step loop
    // above has finished (with whatever mix of fast/investigate/
    // alternate/shrink shaped those steps) -- see sweepFlag's comment on
    // runHillClimb for the full rationale. A no-op unless sweepFlag.
    curScore = runSPRSweep(tree, edgeRegistry, opt, searchState, curScore);

    appendFinalRecordRow(searchState, opt, tree, curScore);

    if (haveSprLocalCache)
        freeSPRLocalLhCache(sprLocalCache);

    if (!quiet)
        cout << endl;
    cout << "=== finished after " << searchState.stepsRun << " step(s) ===" << endl;
    string finalTreeNewick = newickOf(tree);
    cout << "final tree (logL = " << curScore << "): " << finalTreeNewick << endl;

    // RF distance and the true-tree logL reference both need an actual
    // ground-truth tree, which doesn't exist when noTrueTree; rf stays -1
    // (already this function's own "no distance computed" sentinel, same
    // value used further down if rfdist ever comes back empty)
    int rf = -1;
    if (!noTrueTree) {
        cout << "AliSim true tree logL: " << trueTreeLogl << endl;

        stringstream finalTreeStream;
        finalTreeStream << finalTreeNewick;
        finalTreeStream.seekg(0, ios::beg);
        vector<double> rfdist;
        trueTree.computeRFDist(finalTreeStream, rfdist);
        rf = rfdist.empty() ? -1 : (int) rfdist[0];

        cout << "RF distance to the original AliSim tree: " << rf << endl;
    }

    ofstream out("output.txt");
    if (out.good()) {
        if (!noTrueTree) {
            out << "AliSim true tree : " << trueTreeNewick << endl;
            out << "RF distance      : " << rf << endl;
        }
        out << "Final result tree: " << finalTreeNewick << endl;
        out.close();
        cout << "Results written to output.txt" << endl;
    } else {
        cerr << "warning: could not write to output.txt" << endl;
    }

    cout << "time elapsed    : " << fixed << setprecision(2) << (getCPUTime() - cpuClockStart)
         << " sec" << endl;

    delete aln;
    return 2;
}

/**
    apply the SAME sequence of random SPR moves to FOUR separate copies of
    the same starting tree, differing only in how each move's branch
    lengths are set:
      - A keeps applySPR's own naive placeholder lengths (half the target
        edge's length split evenly at the new attachment point, and the
        sum of the two vacated lengths merged back together -- see
        applySPR's own comment in phylotree.cpp) for every move, exactly
        like every OTHER command in this file that doesn't pass "reopt".
      - B re-optimizes (Newton-Raphson, reoptimizeSPREdges) just the 3
        edges each move actually changes, immediately after applying it.
      - C re-optimizes EVERY edge in the whole tree, one full sweep
        (optimizeAllBranches(1)), immediately after applying it.
      - D does the same as C, but up to 10 full sweeps
        (optimizeAllBranches(10)) instead of 1.
    All four are otherwise identical: same starting topology and branch
    lengths (see below), same alignment, same fixed-parameter model (JC for
    DNA, LG for protein -- sequence type auto-detected, see modelNameFor's
    comment), and -- critically -- the same move at every step (same prune
    edge, same regraft target), so any divergence between their logL
    trajectories can only come from how branch lengths were set, never
    from the four trees taking different topological paths.

    This is deliberately NOT a search: every chosen candidate is
    applied unconditionally, every step, on ALL FOUR trees, regardless of
    whether it improves anything -- there is no likelihood comparison to
    accept or reject a move against, and nothing is ever rolled back.
    This is what "the same random moves on [multiple] trees" needs: a
    hill-climb's own accept/reject decision is itself a function of each
    tree's current logL, which would already differ between trees by the
    very quantity this test is trying to isolate, so a search design here
    would let the trees' topologies drift apart after their very first
    accept/reject disagreement.

    STARTING POINT: a single random Yule-Harding tree is built and fully
    branch-length-optimized (clampAllBranchLengthsForOptimization +
    optimizeAllBranches(), the same full ML fit --hillclimb's own "reopt"
    start uses) to become tree A; trees B, C, and D are then each built
    (initClonedTree) by parsing tree A's own full-precision Newick text
    (see the high numeric_precision set just before printTree below --
    newickOf()'s usual 1-decimal rounding is a display-only nicety that
    would bias this comparison). This -- rather than four independent
    generateRandomTree() calls, even with the same seed -- is what
    guarantees all four start out not just topologically identical but
    with bit-for-bit identical branch lengths too.

    KEEPING THE FOUR TREES IN LOCKSTEP: choosePrune() and chooseGraft(),
    this tool's usual random-candidate engine, cannot be called
    independently on each tree -- all draw from the SAME global RNG
    stream (random_int()/random_double()), so calling either one again
    for B/C/D would consume different draws than the call already used
    for A, letting the trees' candidates disagree even though they
    started identical. Instead, each step here:
      1. Draws ONE shared random prune-edge slot index, valid in every
         tree's own edge registry -- all four are built over identical
         topologies (see "starting point" above), and buildEdgeRegistry's
         traversal order depends only on tree structure, never on any
         per-tree state, so slot i names structurally "the same" edge in
         every registry for as long as the four topologies stay
         identical.
      2. Resolves that one edge's prune orientation once via
         resolvePruneOrientationForSlot() (see its own comment), using
         ONE shared coin-flip bit instead of four independent
         choosePrune() calls.
      3. Enumerates every legal regraft target within <radius> hops via
         findGraftPositions() -- separately on each tree, but this is
         safe (unlike chooseGraft) because findGraftPositions is a pure,
         deterministic BFS over each tree's own current structure: no
         randomness at all, so it can't desynchronize the trees by
         itself. If the four trees' candidate counts ever disagree, that
         is a synchronization bug, not a legitimate outcome -- this is
         checked explicitly and the run aborts rather than silently
         continuing on a nonsense comparison.
      4. Draws ONE shared random index into that (equal-length) candidate
         list, and applies the structurally corresponding candidate to
         each tree.

    This means the per-step candidate distribution is NOT exactly
    chooseGraft's own radius-weighted directed walk (see
    graftDistanceWeight) -- it is a uniform pick over every legal target
    findGraftPositions enumerates within <radius> hops instead. That is a
    deliberate, acceptable difference for this test (chooseGraft's own
    weighting was never a requirement here, just a convenient existing
    engine elsewhere), not an attempt to reproduce chooseGraft's exact
    behavior.

    Writes step, then logl_before/logl_after/logl_diff for each of a, b,
    c, d in turn (12 data columns total) to branchlength_compare_data.csv
    (repo root, overwritten each run) -- the absolute before/after values
    are included, not just the diff, since the diff alone doesn't show
    where each tree's likelihood actually sits relative to the others.
    Also prints a shorter step,logl_diff_a,logl_diff_b,logl_diff_c,
    logl_diff_d line per step to stdout.

    C and D are meaningfully more expensive per step than A/B (a full
    optimizeAllBranches() sweep touches every edge in the tree, not just
    the 4 reoptimizeSPREdges touches), so this command is much slower than
    a plain "reopt" hill-climb at the same step count. The final summary
    reports each tree's own average apply+score time per step (excludes
    the prune/graft selection work all four trees pay alike), so the
    relative cost of naive vs. 4-edge reopt vs. 1 vs. 10 full sweeps is
    visible directly, not just inferred from total wall time.

    @return 0 on success, 1 if the trees desynchronize (a bug, should
    never happen in practice -- see point 3 above), 2 if the alignment
    couldn't be found
 */
int runBranchLengthCompare(const string &trueTreeArg, int radius, int maxSteps) {
    double cpuClockStart = getCPUTime();

    // AliSim's own default output naming: <prefix>.treefile + <prefix>.fa
    string alnFile = trueTreeArg;
    const string suffix = ".treefile";
    if (alnFile.size() > suffix.size()
            && alnFile.compare(alnFile.size() - suffix.size(), suffix.size(), suffix) == 0)
        alnFile = alnFile.substr(0, alnFile.size() - suffix.size()) + ".fa";
    else
        alnFile += ".fa";

    ifstream alnCheck(alnFile.c_str());
    if (!alnCheck.good()) {
        cerr << "error: could not find alignment '" << alnFile << "'" << endl;
        cerr << "  (derived from the tree argument by replacing '.treefile' with '.fa',"
                " AliSim's own default output naming; generate a pair with e.g." << endl;
        cerr << "   iqtree3 --alisim <prefix> -m \"GTR{2,4,1,1,4,2}+F{0.3,0.2,0.2,0.3}\""
                " -t \"RANDOM{yh/100}\" --length 10000)" << endl;
        return 2;
    }
    alnCheck.close();

    Params &params = Params::getInstance();
    params.setDefault();

    init_random((int) (time(nullptr) * 1000 + (long) (getRealTime() * 1000) % 1000));

    // silence the library setup noise, same as runHillClimb
    ostringstream suppressedSetupOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedSetupOutput.rdbuf());

    InputType intype;
    // nullptr sequence_type auto-detects DNA vs. protein from content, same
    // as runLikelihood/runHillClimb (see modelNameFor's comment)
    Alignment *aln = new Alignment((char*) alnFile.c_str(), nullptr, intype, "");
    if (!isSupportedSeqType(aln->seq_type)) {
        cout.rdbuf(realCoutBuf);
        cerr << "error: alignment '" << alnFile << "' auto-detected as a sequence type this tool "
                "doesn't support (only DNA and protein are handled)" << endl;
        delete aln;
        return 2;
    }

    string modelName = modelNameFor(aln->seq_type, false);
    // real Newton-Raphson branch-length search runs on every tree here
    // (A's own one-time starting optimization, B's every-step
    // reoptimizeSPREdges, C/D's every-step optimizeAllBranches) -- see
    // runHillClimb's own comment on why this needs the safe/scaled kernel
    params.lk_safe_scaling = true;

    // tree A: random start, then one full ML branch-length fit so ALL
    // FOUR trees start from a real, non-arbitrary baseline (not
    // Yule-Harding's own random lengths) -- see this function's own
    // comment
    PhyloTree treeA(aln);
    treeA.setParams(&params);
    treeA.generateRandomTree(YULE_HARDING);
    treeA.setNumThreads(1);
    treeA.setLikelihoodKernel(LK_SSE2);
    ModelsBlock *modelsBlockA = readModelsDefinition(params);
    treeA.setModelFactory(new ModelFactory(params, modelName, &treeA, modelsBlockA));
    delete modelsBlockA;
    treeA.setModel(treeA.getModelFactory()->model);
    treeA.setRate(treeA.getModelFactory()->site_rate);
    treeA.initializeAllPartialLh();
    clampAllBranchLengthsForOptimization(treeA, params.min_branch_length);
    double curLoglA = treeA.optimizeAllBranches();

    // full-precision starting Newick, WITHOUT WT_SORT_TAXA -- sorting
    // would rewrite tree A's own internal neighbor order into taxon-name
    // order on the way out, so a tree parsed back from this text would
    // start with a DIFFERENT neighbor order than tree A itself, silently
    // breaking the "slot i means the same edge in every registry"
    // invariant this function's lockstep design depends on from step 1
    int savedPrecision = params.numeric_precision;
    params.numeric_precision = 15;
    ostringstream startNewick;
    treeA.printTree(startNewick, WT_BR_LEN);
    params.numeric_precision = savedPrecision;

    // B/C/D are each a fresh, independent parse of tree A's own starting
    // Newick text -- see this function's own "starting point" comment
    PhyloTree treeB, treeC, treeD;
    initClonedTree(treeB, startNewick.str(), aln, params, modelName);
    initClonedTree(treeC, startNewick.str(), aln, params, modelName);
    initClonedTree(treeD, startNewick.str(), aln, params, modelName);
    double curLoglB = treeB.computeLikelihood();
    double curLoglC = treeC.computeLikelihood();
    double curLoglD = treeD.computeLikelihood();

    EdgeRegistry edgeRegistryA, edgeRegistryB, edgeRegistryC, edgeRegistryD;
    buildEdgeRegistry(treeA, edgeRegistryA);
    buildEdgeRegistry(treeB, edgeRegistryB);
    buildEdgeRegistry(treeC, edgeRegistryC);
    buildEdgeRegistry(treeD, edgeRegistryD);

    cout.rdbuf(realCoutBuf);

    cout << "start tree (random Yule-Harding, fully branch-length-optimized): " << newickOf(treeA) << endl;
    cout << "  (logL = " << curLoglA << ", identical for all four copies at the start)" << endl;
    cout << "radius           : " << radius << endl;
    cout << "max steps        : " << maxSteps << endl;
    cout << "tree A           : applySPR's own naive placeholder branch lengths every step (default)" << endl;
    cout << "tree B           : re-optimized (Newton-Raphson) on the 4 edges around each move every step" << endl;
    cout << "tree C           : re-optimized on EVERY edge, 1 full optimizeAllBranches() sweep every step" << endl;
    cout << "tree D           : re-optimized on EVERY edge, up to 10 full optimizeAllBranches() sweeps every step"
         << endl;
    cout << "every candidate is accepted unconditionally on ALL FOUR trees -- this is a random walk, not a search"
         << endl;
    cout << endl;

    ofstream csv("branchlength_compare_data.csv");
    if (csv.good())
        csv << "step,"
               "logl_before_a,logl_after_a,logl_diff_a,"
               "logl_before_b,logl_after_b,logl_diff_b,"
               "logl_before_c,logl_after_c,logl_diff_c,"
               "logl_before_d,logl_after_d,logl_diff_d" << endl;
    else
        cerr << "warning: could not write to branchlength_compare_data.csv" << endl;

    cout << "step,logl_diff_a,logl_diff_b,logl_diff_c,logl_diff_d" << endl;

    // per-tree cumulative time spent applying+scoring a move (excludes
    // the shared prune/graft selection above, which all four trees pay
    // alike) -- reported as a per-step average in the final summary, to
    // show how much C/D's full-tree sweeps actually cost relative to
    // A/B's naive/4-edge handling
    double timeA = 0.0, timeB = 0.0, timeC = 0.0, timeD = 0.0;
    int measuredSteps = 0;

    int step = 0;
    for (; step < maxSteps; step++) {
        if (edgeRegistryA.slots.empty()) {
            cout << "step " << (step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }
        if (edgeRegistryA.slots.size() != edgeRegistryB.slots.size()
                || edgeRegistryA.slots.size() != edgeRegistryC.slots.size()
                || edgeRegistryA.slots.size() != edgeRegistryD.slots.size()) {
            cerr << "error: the four trees' edge registries desynchronized (" << edgeRegistryA.slots.size() << "/"
                 << edgeRegistryB.slots.size() << "/" << edgeRegistryC.slots.size() << "/"
                 << edgeRegistryD.slots.size() << " slots, a/b/c/d) at step " << (step + 1)
                 << " -- this is a bug, aborting." << endl;
            return 1;
        }

        int slot = random_int((int) edgeRegistryA.slots.size());
        bool preferA = random_int(2) == 0; // one shared coin flip, see resolvePruneOrientationForSlot

        pair<PhyloNode*, PhyloNode*> &edgeA = edgeRegistryA.slots[slot];
        pair<PhyloNode*, PhyloNode*> &edgeB = edgeRegistryB.slots[slot];
        pair<PhyloNode*, PhyloNode*> &edgeC = edgeRegistryC.slots[slot];
        pair<PhyloNode*, PhyloNode*> &edgeD = edgeRegistryD.slots[slot];

        PhyloNode *pruneNodeA, *pruneDadA, *pruneNodeB, *pruneDadB, *pruneNodeC, *pruneDadC, *pruneNodeD, *pruneDadD;
        bool okA = resolvePruneOrientationForSlot(treeA, edgeA.first, edgeA.second, preferA, pruneNodeA, pruneDadA);
        bool okB = resolvePruneOrientationForSlot(treeB, edgeB.first, edgeB.second, preferA, pruneNodeB, pruneDadB);
        bool okC = resolvePruneOrientationForSlot(treeC, edgeC.first, edgeC.second, preferA, pruneNodeC, pruneDadC);
        bool okD = resolvePruneOrientationForSlot(treeD, edgeD.first, edgeD.second, preferA, pruneNodeD, pruneDadD);
        if (!okA || !okB || !okC || !okD)
            continue; // unreachable in practice -- none of these trees is ever rooted (see the function comment)

        vector<GraftCandidate> candidatesA = findGraftPositions(treeA, pruneNodeA, pruneDadA, radius);
        vector<GraftCandidate> candidatesB = findGraftPositions(treeB, pruneNodeB, pruneDadB, radius);
        vector<GraftCandidate> candidatesC = findGraftPositions(treeC, pruneNodeC, pruneDadC, radius);
        vector<GraftCandidate> candidatesD = findGraftPositions(treeD, pruneNodeD, pruneDadD, radius);
        if (candidatesA.size() != candidatesB.size() || candidatesA.size() != candidatesC.size()
                || candidatesA.size() != candidatesD.size()) {
            cerr << "error: the four trees found different candidate counts (" << candidatesA.size() << "/"
                 << candidatesB.size() << "/" << candidatesC.size() << "/" << candidatesD.size()
                 << ", a/b/c/d) for structurally the same prune edge at step " << (step + 1)
                 << " -- this is a bug, aborting." << endl;
            return 1;
        }
        if (candidatesA.empty()) {
            cout << "step " << (step + 1) << ": prune {" << describeEdgeCompact(pruneNodeA, pruneDadA)
                 << "} -- no legal graft target; skipping." << endl;
            continue;
        }

        int candIdx = random_int((int) candidatesA.size());

        SPRMove moveA;
        moveA.prune_node = pruneNodeA;
        moveA.prune_dad = pruneDadA;
        moveA.regraft_node = candidatesA[candIdx].node;
        moveA.regraft_dad = candidatesA[candIdx].dad;
        moveA.radius = candidatesA[candIdx].radius;
        moveA.screening_score = 0.0;
        moveA.exact_score = 0.0;
        moveA.candidate_id = 0;
        moveA.generation = step;

        SPRMove moveB = moveA;
        moveB.prune_node = pruneNodeB;
        moveB.prune_dad = pruneDadB;
        moveB.regraft_node = candidatesB[candIdx].node;
        moveB.regraft_dad = candidatesB[candIdx].dad;
        moveB.radius = candidatesB[candIdx].radius;

        SPRMove moveC = moveA;
        moveC.prune_node = pruneNodeC;
        moveC.prune_dad = pruneDadC;
        moveC.regraft_node = candidatesC[candIdx].node;
        moveC.regraft_dad = candidatesC[candIdx].dad;
        moveC.radius = candidatesC[candIdx].radius;

        SPRMove moveD = moveA;
        moveD.prune_node = pruneNodeD;
        moveD.prune_dad = pruneDadD;
        moveD.regraft_node = candidatesD[candIdx].node;
        moveD.regraft_dad = candidatesD[candIdx].dad;
        moveD.radius = candidatesD[candIdx].radius;

        // captured before the move is applied -- see findSPRSiblings' own
        // comment on why this can't be recovered afterward; only B needs
        // this (C/D's optimizeAllBranches sweeps the whole tree via its
        // own internal traversal, not these two specific nodes)
        PhyloNode *sibling1B = nullptr, *sibling2B = nullptr;
        findSPRSiblings(pruneNodeB, pruneDadB, sibling1B, sibling2B);

        TrackedSPR trackedA, trackedB, trackedC, trackedD;
        double t0 = getCPUTime();
        applySPRTracked(treeA, edgeRegistryA, moveA, trackedA);
        resetLikelihoodBuffers(treeA);
        double newLoglA = treeA.computeLikelihood();
        timeA += getCPUTime() - t0;

        // reoptimizeSPREdges'/optimizeAllBranches' optimizeOneBranch calls
        // need valid partial likelihood buffers for the NEW topology
        // already in place before they run -- see scoreTrialSPRMove's own
        // identical reset-then-reoptimize order, which this mirrors for
        // B, C, and D alike
        t0 = getCPUTime();
        applySPRTracked(treeB, edgeRegistryB, moveB, trackedB);
        resetLikelihoodBuffers(treeB);
        reoptimizeSPREdges(treeB, moveB.prune_dad, moveB.regraft_dad, moveB.regraft_node, moveB.prune_node,
                sibling1B, sibling2B);
        resetLikelihoodBuffers(treeB);
        double newLoglB = treeB.computeLikelihood();
        timeB += getCPUTime() - t0;

        t0 = getCPUTime();
        applySPRTracked(treeC, edgeRegistryC, moveC, trackedC);
        resetLikelihoodBuffers(treeC);
        double newLoglC = treeC.optimizeAllBranches(1);
        timeC += getCPUTime() - t0;

        t0 = getCPUTime();
        applySPRTracked(treeD, edgeRegistryD, moveD, trackedD);
        resetLikelihoodBuffers(treeD);
        double newLoglD = treeD.optimizeAllBranches(10);
        timeD += getCPUTime() - t0;
        measuredSteps++;

        double loglDiffA = newLoglA - curLoglA;
        double loglDiffB = newLoglB - curLoglB;
        double loglDiffC = newLoglC - curLoglC;
        double loglDiffD = newLoglD - curLoglD;

        cout << (step + 1) << "," << loglDiffA << "," << loglDiffB << "," << loglDiffC << "," << loglDiffD << endl;
        if (csv.good())
            csv << (step + 1) << ","
                << curLoglA << "," << newLoglA << "," << loglDiffA << ","
                << curLoglB << "," << newLoglB << "," << loglDiffB << ","
                << curLoglC << "," << newLoglC << "," << loglDiffC << ","
                << curLoglD << "," << newLoglD << "," << loglDiffD << endl;

        curLoglA = newLoglA;
        curLoglB = newLoglB;
        curLoglC = newLoglC;
        curLoglD = newLoglD;
    }
    if (csv.good())
        csv.close();

    cout << endl;
    cout << "=== finished after " << step << " step(s) ===" << endl;
    cout << "final tree A (naive lengths, logL = " << curLoglA << "): " << newickOf(treeA) << endl;
    cout << "final tree B (4-edge reopt, logL = " << curLoglB << ")" << endl;
    cout << "final tree C (1 full sweep/step, logL = " << curLoglC << ")" << endl;
    cout << "final tree D (10 full sweeps/step, logL = " << curLoglD << ")" << endl;
    cout << "per-step data written to branchlength_compare_data.csv" << endl;
    if (measuredSteps > 0)
        cout << "avg time/step    : A=" << fixed << setprecision(1) << (timeA / measuredSteps * 1000.0)
             << "ms  B=" << (timeB / measuredSteps * 1000.0) << "ms  C=" << (timeC / measuredSteps * 1000.0)
             << "ms  D=" << (timeD / measuredSteps * 1000.0) << "ms  (" << measuredSteps << " step(s) measured)"
             << endl;
    cout << "time elapsed     : " << fixed << setprecision(2) << (getCPUTime() - cpuClockStart)
         << " sec" << endl;

    delete aln;
    return 0;
}

int runSelfTest() {
    // mirrors test_scripts/test_data/spr/six_taxa.start.tree: groups (A,C) and (B,D)
    string newick = "((A:0.10,C:0.10):0.10,(B:0.10,D:0.10):0.10,(E:0.10,F:0.10):0.10);";

    cout << "=== SPR topology test ===" << endl;
    cout << "start tree: " << newick << endl;

    PhyloTree tree;
    tree.read_TreeString(newick, false);

    string original = newickOf(tree);
    cout << "parsed as : " << original << endl;

    PhyloNode *C = requireLeaf(tree, "C");
    PhyloNode *D = requireLeaf(tree, "D");
    PhyloNode *A = requireLeaf(tree, "A");
    if (!C || !D || !A)
        return 2;

    // C's only neighbor is its dad (the (A,C) cherry node); D's only
    // neighbor is its dad (the (B,D) cherry node)
    PhyloNode *dadOfC = (PhyloNode*) C->neighbors[0]->node;
    PhyloNode *dadOfD = (PhyloNode*) D->neighbors[0]->node;

    // capture dadOfC's other two neighbors (A and the deep trifurcation
    // node) before the move, so we can check they end up directly
    // connected to each other once dadOfC is bypassed
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    FOR_NEIGHBOR_DECLARE(dadOfC, C, it)
    {
        if (!sibling1)
            sibling1 = (PhyloNode*) (*it)->node;
        else
            sibling2 = (PhyloNode*) (*it)->node;
    }

    // --- Move 1: prune C from (A,C), regraft it onto the (B,D) edge ------
    SPRMove move;
    move.prune_node = C;
    move.prune_dad = dadOfC;
    move.regraft_node = D;
    move.regraft_dad = dadOfD;
    move.radius = 1;
    move.screening_score = 0.0;
    move.exact_score = 0.0;
    move.candidate_id = 0;
    move.generation = 0;

    cout << endl << "-- legality checks --" << endl;
    expect(tree.isLegalSPR(move), "prune C from (A,C), regraft onto (B,D) edge is legal");

    // illegal: regraft target incident to the node being suppressed
    SPRMove identityMove = move;
    identityMove.regraft_node = A;
    identityMove.regraft_dad = dadOfC;
    expect(!tree.isLegalSPR(identityMove), "regrafting back onto the node being suppressed is rejected");

    // illegal: regraft target two hops inside the pruned subtree. Prune the
    // deep trifurcation node (X) away from dadOfC; X's subtree contains the
    // (B,D) cherry, so targeting the D-dadOfD edge must be rejected even
    // though neither endpoint is X or dadOfC themselves.
    PhyloNode *deepNode = sibling1->isLeaf() ? sibling2 : sibling1;
    expect(!deepNode->isLeaf(), "found the deep (non-leaf) sibling of C's dad to use as prune_node");
    SPRMove insideSubtree;
    insideSubtree.prune_node = deepNode;
    insideSubtree.prune_dad = dadOfC;
    insideSubtree.regraft_node = D;
    insideSubtree.regraft_dad = dadOfD;
    expect(!tree.isLegalSPR(insideSubtree), "regrafting two hops inside the pruned subtree is rejected");

    cout << endl << "-- apply --" << endl;
    SPRRollback rollback;
    tree.applySPR(move, rollback);

    string afterApply = newickOf(tree);
    cout << "after move: " << afterApply << endl;

    expect(afterApply != original, "topology actually changed after applySPR");
    expect(dadOfC->isNeighbor(C), "prune_dad still connects to the pruned leaf");
    expect(dadOfC->isNeighbor(D), "prune_dad now connects to the regraft leaf");
    expect(dadOfC->isNeighbor(dadOfD), "prune_dad now connects to the regraft dad");
    expect(dadOfC->degree() == 3, "prune_dad is still a valid degree-3 internal node");
    expect(!dadOfD->isNeighbor(D), "regraft dad no longer connects directly to the regraft leaf");
    expect(sibling1->isNeighbor(sibling2), "the two original siblings of the pruned edge are now directly connected");
    expect(!sibling1->isNeighbor(dadOfC) && !sibling2->isNeighbor(dadOfC),
            "the original siblings no longer connect through prune_dad");

    cout << endl << "-- rollback --" << endl;
    tree.rollbackSPR(rollback);
    string afterRollback = newickOf(tree);
    cout << "after rollback: " << afterRollback << endl;

    expect(afterRollback == original, "tree exactly matches the original topology and branch lengths after rollback");
    expect(dadOfC->isNeighbor(sibling1) && dadOfC->isNeighbor(sibling2), "prune_dad's original neighbors are restored");
    expect(!dadOfC->isNeighbor(D), "prune_dad no longer connects to the regraft leaf after rollback");

    cout << endl;
    if (g_failures == 0) {
        cout << "ALL CHECKS PASSED" << endl;
        return 0;
    } else {
        cout << g_failures << " CHECK(S) FAILED" << endl;
        return 1;
    }
}

void printUsage(const char *prog) {
    cerr << "Usage:" << endl;
    cerr << "  " << prog << endl;
    cerr << "      run the built-in self-test (default, no arguments)" << endl;
    cerr << endl;
    cerr << "  " << prog << " <tree.nwk | \"(newick,string);\"> <prune-edge> <regraft-edge>" << endl;
    cerr << "      prune <prune-edge> and regraft it onto <regraft-edge>, printing the tree" << endl;
    cerr << "      before and after the move." << endl;
    cerr << endl;
    cerr << "  " << prog << " --list-grafts <tree.nwk | \"(newick,string);\"> <prune-edge> <radius>" << endl;
    cerr << "      prune <prune-edge> and list every distinct, legal graft position within" << endl;
    cerr << "      <radius> hops of the prune point, without applying any of them. radius 1 is" << endl;
    cerr << "      the nearest possible legal target (equivalent to an NNI)." << endl;
    cerr << endl;
    cerr << "  " << prog << " --likelihood <tree.nwk | \"(newick,string);\"> <alignment.fasta>" << endl;
    cerr << "      evaluate the log-likelihood of the given tree (topology and branch" << endl;
    cerr << "      lengths as given, no optimization) against a DNA or protein alignment" << endl;
    cerr << "      (sequence type auto-detected) under a plain fixed model (JC for DNA, LG" << endl;
    cerr << "      for protein). Sequence names in the alignment must match the tree's" << endl;
    cerr << "      leaf names exactly." << endl;
    cerr << endl;
    cerr << "  " << prog << " --hillclimb <alisim-tree.treefile> <radius> <max-steps> [random] [iqtreestart [N]] [starttree <path>] [model <spec>] [fast [N]] [quiet] [reopt] [fullreopt M N] [gtr] [record] [recordtopology] [trajectory] [investigate [N]] [alternate] [shrink [N]] [learnradius [N]] [sweep [N]] [findopt [N]] [notree] [distradius] [weightprune]" << endl;
    cerr << "      greedy randomized SPR search: build a BioNJ start tree from the" << endl;
    cerr << "      alignment AliSim simulated from <alisim-tree.treefile> (found by" << endl;
    cerr << "      replacing '.treefile' with '.fa'), then repeatedly prune a random edge," << endl;
    cerr << "      evaluate every legal regraft within <radius> hops via applySPR/" << endl;
    cerr << "      rollbackSPR on one tree object, and keep the best if it improves the" << endl;
    cerr << "      likelihood, for up to <max-steps> rounds. Prints the RF distance to the" << endl;
    cerr << "      original AliSim tree and writes both trees + the RF distance to" << endl;
    cerr << "      output.txt (skipped with 'notree', see below). Eighteen optional trailing" << endl;
    cerr << "      flags, in any order:" << endl;
    cerr << "        random     start from a random Yule-Harding topology instead of the" << endl;
    cerr << "                   default BioNJ estimate tree" << endl;
    cerr << "        iqtreestart [N]  start from buildIQTreeStyleStartTree's own result instead:" << endl;
    cerr << "                   a pool of N (default 20) randomized-stepwise-addition parsimony" << endl;
    cerr << "                   trees (PhyloTree::computeParsimonyTree), each given a light" << endl;
    cerr << "                   branch-length polish and scored, keeping the single best -- mirrors" << endl;
    cerr << "                   the SHAPE of real IQ-TREE's own preprocessing" << endl;
    cerr << "                   (IQTree::computeInitialTree/initCandidateTreeSet) at a smaller" << endl;
    cerr << "                   scale, scored under 'gtr' when that flag is also given. Does NOT" << endl;
    cerr << "                   also NNI/SPR-refine the winner itself -- the step loop below" << endl;
    cerr << "                   (radius/fast/reopt/etc., whatever was actually given) does that," << endl;
    cerr << "                   same as it already does from a plain BioNJ or random start. Mutually" << endl;
    cerr << "                   exclusive with 'random'. EXPERIMENTAL" << endl;
    cerr << "        fast [N]   pick each step's prune edge via choosePrune() and draw N (default" << endl;
    cerr << "                   1) independent regraft targets via chooseGraft() -- an" << endl;
    cerr << "                   O(1)/O(distance) random proposal -- from that same prune position," << endl;
    cerr << "                   score each by real likelihood, and keep the best of the group," << endl;
    cerr << "                   applied and kept-or-reverted directly, instead of enumerating and" << endl;
    cerr << "                   scoring every candidate in the radius. 'fast' alone (N omitted) is" << endl;
    cerr << "                   the original single-candidate behavior; the N after 'fast' is only" << endl;
    cerr << "                   consumed if it actually parses as a positive integer, so 'fast" << endl;
    cerr << "                   quiet' still works (quiet is not mistaken for a candidate count)" << endl;
    cerr << "        quiet      suppress the one printed line per step; only the setup header" << endl;
    cerr << "                   and final summary (final tree, RF distance, time elapsed) are" << endl;
    cerr << "                   printed. With max-steps in the thousands, this also avoids the" << endl;
    cerr << "                   per-line flush stalling on a slow interactive console" << endl;
    cerr << "        reopt      re-optimize (Newton-Raphson) the 4 edges around an SPR move -- the" << endl;
    cerr << "                   3 it actually changes, plus the pruned subtree's own upper edge" << endl;
    cerr << "                   (which now sits between a different pair of neighbors, even though" << endl;
    cerr << "                   applySPR itself never touches it) -- before scoring a candidate," << endl;
    cerr << "                   the same way IQ-TREE's own NNI search re-optimizes the branches it" << endl;
    cerr << "                   touches -- instead of trusting applySPR's naive placeholder" << endl;
    cerr << "                   lengths (half the target edge split evenly, the two vacated edges" << endl;
    cerr << "                   summed) for the first 3, or that upper edge's own pre-move length" << endl;
    cerr << "                   for the 4th. When reopt (or" << endl;
    cerr << "                   fullreopt, below) is on, the whole starting tree is also" << endl;
    cerr << "                   ML-optimized once up front (PhyloTree::optimizeAllBranches(), the" << endl;
    cerr << "                   same full-tree sweep IQ-TREE's own search uses on a freshly built" << endl;
    cerr << "                   tree) instead of trusting BioNJ/random-start lengths as the" << endl;
    cerr << "                   search's departure point, and an accepted candidate's re-optimized" << endl;
    cerr << "                   lengths are kept in the tree instead of being discarded after" << endl;
    cerr << "                   scoring. Can only ever improve or leave unchanged a candidate's" << endl;
    cerr << "                   reported likelihood for its topology. A bare flag, no numeric" << endl;
    cerr << "                   argument -- independent of fullreopt below, composes freely either" << endl;
    cerr << "                   way." << endl;
    cerr << "                   EXPERIMENTAL and noticeably slower per candidate" << endl;
    cerr << "                   (real branch-length search plus the safe/scaled likelihood" << endl;
    cerr << "                   kernel this needs -- see scoreTrialSPRMove's comment in the" << endl;
    cerr << "                   source for why); the up-front optimizeAllBranches() sweep alone" << endl;
    cerr << "                   already brings the starting tree's logL close to the true" << endl;
    cerr << "                   simulated tree's own (-7.17e+05 vs. -7.03e+05 on sim.treefile," << endl;
    cerr << "                   before any SPR moves at all)" << endl;
    cerr << "        fullreopt M N  independent of 'reopt' (it used to only be reachable as" << endl;
    cerr << "                   'reopt's own optional trailing number): run one full" << endl;
    cerr << "                   optimizeAllBranches(M) sweep over every edge every N SUCCESSFUL" << endl;
    cerr << "                   (accepted) steps -- not every N steps attempted -- on top of" << endl;
    cerr << "                   whatever 'reopt' itself is or isn't doing per candidate. BOTH M" << endl;
    cerr << "                   (round-count ceiling) and N (successful-step interval) are" << endl;
    cerr << "                   REQUIRED -- unlike every other numeric flag here, there's no" << endl;
    cerr << "                   sensible single-number default for 'M rounds every N successful" << endl;
    cerr << "                   steps'. (mirrors IQ-TREE's own NNI" << endl;
    cerr << "                   loop, which does a full single-round sweep after every batch of" << endl;
    cerr << "                   applied moves). EXPERIMENTAL, and in informal testing on" << endl;
    cerr << "                   sim.treefile showed no measurable logL/RF improvement over plain" << endl;
    cerr << "                   'reopt' at N=1, 5, or 10 (30-100 steps) despite costing up to" << endl;
    cerr << "                   ~3-4x more wall time at N=1 -- see fullReoptEveryNSteps' comment" << endl;
    cerr << "                   in the source for the numbers" << endl;
    cerr << "        fullreopt M N true  optional third token (a literal word, never confused" << endl;
    cerr << "                   with a third required number): also ML-optimize the whole" << endl;
    cerr << "                   starting tree once up front, the same up-front fit 'reopt' always" << endl;
    cerr << "                   does regardless of this modifier. Without it, 'fullreopt' alone" << endl;
    cerr << "                   (no 'reopt') leaves the starting tree's branch lengths untouched" << endl;
    cerr << "                   until the first periodic checkpoint. See fullReoptInitialFit's" << endl;
    cerr << "                   comment in the source" << endl;
    cerr << "        gtr        search under GTR+FO (ML-estimated rates/frequencies) instead of" << endl;
    cerr << "                   JC for a DNA alignment, or LG+FO instead of plain LG for a" << endl;
    cerr << "                   protein one (sequence type auto-detected) -- relevant since" << endl;
    cerr << "                   sim.fa is simulated under a real GTR+F model, so JC is a genuine" << endl;
    cerr << "                   misspecification for the DNA case, not just a simplification." << endl;
    cerr << "                   Combined with 'fullreopt M N', periodic sweeps" << endl;
    cerr << "                   also re-fit the model's own rate/frequency parameters via" << endl;
    cerr << "                   ModelFactory::optimizeParameters(), not just branch lengths." << endl;
    cerr << "                   EXPERIMENTAL, and showed the same result as plain 'fullreopt M N':" << endl;
    cerr << "                   no measurable logL/RF benefit from periodic re-fitting over" << endl;
    cerr << "                   fitting the model once up front -- see useGtrModel's comment" << endl;
    cerr << "                   in the source for the numbers" << endl;
    cerr << "        record     append this run's convergence trajectory to a" << endl;
    cerr << "                   model-and-search-mode-specific CSV spreadsheet," << endl;
    cerr << "                   record_<model><tag>.csv (repo root; '+' in the model name" << endl;
    cerr << "                   sanitized to '_'; <tag> encodes which of fast/reopt/fullreopt/" << endl;
    cerr << "                   investigate were also given, e.g. record_JC_fast_reopt.csv, so" << endl;
    cerr << "                   different search modes never mix in the same file), one row" << endl;
    cerr << "                   every time curScore actually changes (an accepted step or a" << endl;
    cerr << "                   periodic fullreopt change), PLUS one unconditional final row" << endl;
    cerr << "                   once the run ends regardless of whether its last step was ever" << endl;
    cerr << "                   accepted, so the file's last row always reflects the true final" << endl;
    cerr << "                   state: this run's id (timestamp + flags used), candidates" << endl;
    cerr << "                   evaluated so far, CPU-clock seconds elapsed, the current logL," << endl;
    cerr << "                   and how far that logL still is below the true AliSim tree's own" << endl;
    cerr << "                   logL under this model (true_minus_current)." << endl;
    cerr << "                   Repeated runs using the same search mode APPEND (never" << endl;
    cerr << "                   overwrite) so their trajectories accumulate side by side in the" << endl;
    cerr << "                   same file for later comparison -- see appendRecordRow's and" << endl;
    cerr << "                   buildRecordTag's comments in the source" << endl;
    cerr << "        recordtopology  requires 'record'. Additionally appends this row's tree" << endl;
    cerr << "                   TOPOLOGY (no branch lengths) as one more Newick line to a" << endl;
    cerr << "                   companion file, topology_<model><tag>.nwk (same naming as" << endl;
    cerr << "                   'record's own CSV, '.nwk' instead of '.csv'), positionally" << endl;
    cerr << "                   aligned with the CSV's own rows (line N of this file <-> row N" << endl;
    cerr << "                   of the CSV). Off by default -- plain 'record' alone writes only" << endl;
    cerr << "                   the CSV, exactly as before this flag existed -- since the extra" << endl;
    cerr << "                   Newick line costs real time on large trees/long runs and most" << endl;
    cerr << "                   'record' uses have no need for it. See appendRecordRow's comment" << endl;
    cerr << "                   in the source" << endl;
    cerr << "        trajectory  independent of 'record'/'recordtopology' -- needs neither. Writes" << endl;
    cerr << "                   the tree's topology (no branch lengths) to its own" << endl;
    cerr << "                   trajectory_<run-id>.nwk, one Newick line right after the starting" << endl;
    cerr << "                   tree is built (post-initial-tree), then one more after every" << endl;
    cerr << "                   ACCEPTED step -- never a reverted one, and never a periodic" << endl;
    cerr << "                   fullreopt/findopt refit (neither ever changes the topology). Each" << endl;
    cerr << "                   run gets its own file (keyed by run-id) rather than sharing one" << endl;
    cerr << "                   across runs the way record/recordtopology do, since there is no" << endl;
    cerr << "                   per-line run marker to tell separate runs' trajectories apart" << endl;
    cerr << "                   within one file. See appendTrajectoryTopology's comment in the" << endl;
    cerr << "                   source" << endl;
    cerr << "        investigate N  the step right after any accepted move re-prunes that SAME" << endl;
    cerr << "                   (node,dad) pair and exhaustively scores every legal regraft" << endl;
    cerr << "                   candidate within N real hops of there (not stepRadius), keeping" << endl;
    cerr << "                   the best if it improves; repeats on the step after THAT one if it" << endl;
    cerr << "                   does, until a refinement attempt fails to improve, at which point" << endl;
    cerr << "                   normal choosePrune-based selection resumes until another move is" << endl;
    cerr << "                   accepted. radius 1 (real hops, not stepRadius -- see --list-grafts)" << endl;
    cerr << "                   is the nearest possible legal target (an NNI); N omitted defaults" << endl;
    cerr << "                   to 1. N is included in 'record's own spreadsheet filename (e.g." << endl;
    cerr << "                   record_JC_investigate3.csv), not just each row's run_id, since it" << endl;
    cerr << "                   can change an investigation step's cost/behavior substantially." << endl;
    cerr << "                   Does not add extra steps -- each investigation attempt still" << endl;
    cerr << "                   consumes one of <max-steps>. EXPERIMENTAL --" << endl;
    cerr << "                   see investigateFlag's comment in the source" << endl;
    cerr << "        alternate  toggle every other step between a plain SPR search at <radius>" << endl;
    cerr << "                   (even step indices) and an NNI-equivalent search forced to" << endl;
    cerr << "                   radius 1 (odd step indices) -- 'radius 1' meaning" << endl;
    cerr << "                   findGraftPositions' own nearest-legal-target tier (see" << endl;
    cerr << "                   --list-grafts), a real NNI move. Substitutes for <radius>" << endl;
    cerr << "                   wherever it feeds candidate generation -- 'fast' and the exhaustive" << endl;
    cerr << "                   scan alike -- so it composes with both, no incompatibility check" << endl;
    cerr << "                   needed. 'investigate' still takes priority on a step it's" << endl;
    cerr << "                   actively refining (investigateRadius wins there); 'alternate'" << endl;
    cerr << "                   only governs steps investigate isn't currently overriding. Step" << endl;
    cerr << "                   lines print '(spr, radius N)' / '(nni)' in place of the usual" << endl;
    cerr << "                   '(radius N)' while active. EXPERIMENTAL -- see alternateFlag's" << endl;
    cerr << "                   comment in the source" << endl;
    cerr << "        shrink N   replace the step's own radius (fixed <radius>, overridden" << endl;
    cerr << "                   whenever this is on) with a value that only ever narrows, and" << endl;
    cerr << "                   only in response to actual search stagnation: starts at <radius>" << endl;
    cerr << "                   (floor 1), and drops by 1 once N consecutive non-improving steps" << endl;
    cerr << "                   have piled up, then the stall count resets. N is held CONSTANT for" << endl;
    cerr << "                   the whole run, not scaled down late in it -- every radius gets the" << endl;
    cerr << "                   same fair, full-length chance to prove itself stalled before" << endl;
    cerr << "                   narrowing past it. Every actually-improving step (fast or exhaustive)" << endl;
    cerr << "                   resets the stall count to 0 regardless of radius. N omitted" << endl;
    cerr << "                   defaults to 10. Composes with 'alternate'/'investigate' the same" << endl;
    cerr << "                   way 'alternate' itself does -- it only changes what the step's" << endl;
    cerr << "                   radius currently IS, which every other flag already reads through." << endl;
    cerr << "                   EXPERIMENTAL, including the default threshold -- see" << endl;
    cerr << "                   maybeShrinkRadius' comment in the source" << endl;
    cerr << "        learnradius N   replace the step's own radius (mutually exclusive with" << endl;
    cerr << "                   'shrink' -- both replace the same value, so giving both is a parse" << endl;
    cerr << "                   error) with one drawn fresh each step from a Gamma distribution" << endl;
    cerr << "                   (method-of-moments fit, spread capped so it can never fully collapse)" << endl;
    cerr << "                   over the last N successfully-accepted moves' own achieved radius," << endl;
    cerr << "                   mixed with a STARTING distribution centered on <radius> -- flat below" << endl;
    cerr << "                   it, tapering linearly to zero at a structural ceiling (every edge the" << endl;
    cerr << "                   tree has, or 100% of its own branch length under 'distradius') -- in" << endl;
    cerr << "                   proportion to how many of those N slots are filled so far. <radius> is" << endl;
    cerr << "                   the MIDDLE of that starting spread here, not a cap: the Gamma fit is" << endl;
    cerr << "                   free to drift past it, all the way to that same structural ceiling," << endl;
    cerr << "                   once real history supports it. The starting-distribution share never" << endl;
    cerr << "                   reaches 0% either, even once N successes have accumulated -- a small," << endl;
    cerr << "                   permanent minimum stays in force for the rest of the run, together" << endl;
    cerr << "                   with the capped Gamma spread guaranteeing the whole distribution keeps" << endl;
    cerr << "                   responding to new data rather than ever locking onto one value." << endl;
    cerr << "                   Composed with 'investigate', each entry is not any one investigation" << endl;
    cerr << "                   attempt's own local radius, but the NET distance, start to finish, of" << endl;
    cerr << "                   the whole excursion (the initiating move plus every subsequent" << endl;
    cerr << "                   investigate refinement that kept improving). Composed with" << endl;
    cerr << "                   'distradius', draws stay fractional percentages (see distradius' own" << endl;
    cerr << "                   entry above) instead of being rounded to hop-count integers first," << endl;
    cerr << "                   and what feeds the history is the equivalent radiusPercent actually" << endl;
    cerr << "                   used, not a walk-step/hop count -- so the two stay in the same unit" << endl;
    cerr << "                   regardless of which mode is active. N omitted defaults to 20." << endl;
    cerr << "                   EXPERIMENTAL, including the default window size -- see" << endl;
    cerr << "                   learnRadiusContinuous'/sampleStartingRadius'/learnradiusFlag's" << endl;
    cerr << "                   comments in the source" << endl;
    cerr << "        sweep N    AFTER every step above finishes, rank every internal edge of the" << endl;
    cerr << "                   tree by how well its own current 'siblings' are supported relative" << endl;
    cerr << "                   to the two single-NNI alternative regroupings around that same edge" << endl;
    cerr << "                   (adjacent_subtree_compatibility.pdf, sec. 8), then run an exhaustive," << endl;
    cerr << "                   whole-tree regraft search on the N LEAST-compatible sibling pairs" << endl;
    cerr << "                   found, one at a time, keeping each move that improves (same" << endl;
    cerr << "                   exhaustive-search machinery, and the same reopt handling, the" << endl;
    cerr << "                   tool's original one-pass-per-edge 'sweep' used). This" << endl;
    cerr << "                   ranking is a FIXED list, computed once from the tree as the step" << endl;
    cerr << "                   loop left it, not recomputed after each sweep move improves -- so," << endl;
    cerr << "                   like the original design, NOT guaranteed to reach the SPR-optimal" << endl;
    cerr << "                   tree. No longer touches a step's own selection logic at all (that" << endl;
    cerr << "                   was the ORIGINAL design), so it now composes freely with every other" << endl;
    cerr << "                   flag, including 'investigate'. May optionally be" << endl;
    cerr << "                   immediately followed by a positive integer -- same parsing special" << endl;
    cerr << "                   case as 'fast N'/'shrink N'; N omitted defaults to 10. EXPERIMENTAL" << endl;
    cerr << "                   -- see sweepFlag's and computeSiblingCompatibilityScore's comments" << endl;
    cerr << "                   in the source" << endl;
    cerr << "        findopt [N]  every N steps (default: the total number of steps, i.e. once,"  << endl;
    cerr << "                   effectively at the end), run ONE whole-tree ML refit on a scratch" << endl;
    cerr << "                   clone of the current tree, ALWAYS under the richest available model" << endl;
    cerr << "                   (GTR+FO for DNA, LG+FO for protein) regardless of whether" << endl;
    cerr << "                   the main search itself is using 'gtr' or not (ModelFactory::" << endl;
    cerr << "                   optimizeParameters, jointly refitting branch lengths and the model's" << endl;
    cerr << "                   rate/frequency parameters together) -- PURELY as a diagnostic: the" << endl;
    cerr << "                   scratch clone is discarded once its logL is read off, so the main" << endl;
    cerr << "                   tree, its model, and the real search that continues past this" << endl;
    cerr << "                   point are completely unaffected -- curScore itself is never written" << endl;
    cerr << "                   to. Reports what logL a full refit would reach from here, printed" << endl;
    cerr << "                   alongside curScore for comparison (or, with 'record', its own row in" << endl;
    cerr << "                   the CSV -- 'findopt' is folded into the record filename's tag so" << endl;
    cerr << "                   these diagnostic-only rows don't mix into a file expecting only" << endl;
    cerr << "                   genuine search progress). May optionally be immediately followed by" << endl;
    cerr << "                   a positive integer -- same parsing special case as 'fast N'/" << endl;
    cerr << "                   'shrink N'. Because the refit+restore round trip is real work with" << endl;
    cerr << "                   no counterpart in the actual search's own cost, its CPU-clock cost" << endl;
    cerr << "                   is excluded from the run's own timing entirely ('pausing the timer'" << endl;
    cerr << "                   around it, both for its own CSV row and every later one)." << endl;
    cerr << "                   EXPERIMENTAL -- see maybeRunFindopt's comment in the source" << endl;
    cerr << "        starttree <path>  start from this Newick tree (a file path, or a literal" << endl;
    cerr << "                   Newick string) exactly as given -- topology AND branch lengths" << endl;
    cerr << "                   are trusted as-is, e.g. the final tree of a real iqtree3 NNI run" << endl;
    cerr << "                   already done on this same alignment. Skips BioNJ/random/" << endl;
    cerr << "                   iqtreestart entirely, and (unlike the plain BioNJ default) also" << endl;
    cerr << "                   skips the up-front full branch-length/model re-fit -- this tree is" << endl;
    cerr << "                   used exactly as supplied, with no 'initial tree find' of any kind." << endl;
    cerr << "                   Mutually exclusive with 'random'/'iqtreestart'" << endl;
    cerr << "        model <spec>  build the tree's ModelFactory from this literal model string" << endl;
    cerr << "                   (e.g. 'LG+F{0.081,0.055,...}', 20 comma-separated frequencies in" << endl;
    cerr << "                   the same A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V order IQ-TREE's" << endl;
    cerr << "                   own .iqtree report prints its 'pi(X) = ...' lines in) INSTEAD of" << endl;
    cerr << "                   this run's own short 'LG'/'LG+FO' label. Doesn't affect that label" << endl;
    cerr << "                   itself -- record/run-id filenames still go by 'gtr' as usual -- only" << endl;
    cerr << "                   which actual model gets built. A '+F{...}' spec (fixed values, not" << endl;
    cerr << "                   '+FO', which would still re-optimize from scratch) is how real" << endl;
    cerr << "                   ML-fit parameters from an external run -- e.g. iqtree3's own -- can" << endl;
    cerr << "                   be carried over: combine with 'starttree' and 'notree' so neither" << endl;
    cerr << "                   the topology/branch-lengths NOR the model get re-derived, e.g.:" << endl;
    cerr << "                     iqtree3 -s real.fa -m LG+FO --prefix nni_run" << endl;
    cerr << "                     spr_topology_test --hillclimb real.fa 10 10000 fast quiet notree \\" << endl;
    cerr << "                         starttree nni_run.treefile model \"LG+F{<20 pi(X) values" << endl;
    cerr << "                         from nni_run.iqtree, same order>}\" record" << endl;
    cerr << "        notree     no ground-truth tree: <alisim-tree.treefile> is instead a real" << endl;
    cerr << "                   alignment file path directly (any format/sequence type Alignment" << endl;
    cerr << "                   can auto-detect -- FASTA/NEXUS/PHYLIP, DNA/protein; no '.fa'" << endl;
    cerr << "                   derivation happens). Skips the final RF-distance and true-tree-" << endl;
    cerr << "                   logL output entirely (output.txt still gets the final tree, just" << endl;
    cerr << "                   not the true tree or RF line); 'record' still works, writing NaN" << endl;
    cerr << "                   for the 'gap to true tree' column since there's no true tree to" << endl;
    cerr << "                   compare against. Every other flag works the same as usual" << endl;
    cerr << "        distradius only affects 'fast' mode's own candidate draws: replaces" << endl;
    cerr << "                   chooseGraft's directed random walk (fixed number of edge hops) with" << endl;
    cerr << "                   chooseGraftByDistance, which reinterprets <radius> as a PERCENTAGE" << endl;
    cerr << "                   of the tree's own total branch length (normalized -- the same value" << endl;
    cerr << "                   means the same thing regardless of this dataset's branch-length" << endl;
    cerr << "                   scale) and walks by actual summed branch length instead of hop" << endl;
    cerr << "                   count, stopping ON whichever edge exhausts that budget (not the one" << endl;
    cerr << "                   before it). Sideways/backward hops refund the distance of the edge" << endl;
    cerr << "                   they abandon rather than spend it, since neither is real outward" << endl;
    cerr << "                   progress from the prune point; a generous hop-count cap guarantees" << endl;
    cerr << "                   termination regardless. No effect outside 'fast' mode. Composes with" << endl;
    cerr << "                   'learnradius': draws stay fractional percentages (not rounded to" << endl;
    cerr << "                   whole numbers first), and what feeds learnradius' own history is the" << endl;
    cerr << "                   equivalent radiusPercent actually spent, not a walk-step count, so" << endl;
    cerr << "                   the two stay in the same unit. EXPERIMENTAL -- see" << endl;
    cerr << "                   chooseGraftByDistance's comment in the source" << endl;
    cerr << "        weightprune  replace choosePrune's uniform random edge pick with one" << endl;
    cerr << "                   weighted by each edge's own branch length -- long branches are" << endl;
    cerr << "                   proportionally more likely to be chosen as the prune edge than" << endl;
    cerr << "                   short ones, instead of every edge being equally likely regardless" << endl;
    cerr << "                   of length. A bare flag, no numeric argument of its own; independent" << endl;
    cerr << "                   of every other flag here (only changes WHICH edge gets pruned, never" << endl;
    cerr << "                   how the graft search that follows proceeds), so it composes freely" << endl;
    cerr << "                   with 'fast'/exhaustive, 'distradius', 'learnradius', 'investigate'," << endl;
    cerr << "                   etc. EXPERIMENTAL -- see choosePrune's comment in the source" << endl;
    cerr << endl;
    cerr << "  " << prog << " --branchlength-compare <alisim-tree.treefile> <radius> <max-steps>" << endl;
    cerr << "      NOT a search: applies the SAME sequence of random SPR moves to FOUR" << endl;
    cerr << "      separate copies of one starting tree (see runBranchLengthCompare's own" << endl;
    cerr << "      comment for how the copies are kept in lockstep without reusing" << endl;
    cerr << "      choosePrune()/chooseGraft(), which would desync them). Copy A keeps" << endl;
    cerr << "      applySPR's own naive placeholder branch lengths every move (the default" << endl;
    cerr << "      everywhere else in this tool); copy B re-optimizes (Newton-Raphson) just" << endl;
    cerr << "      the 4 edges around every move; copy C re-optimizes EVERY edge, one" << endl;
    cerr << "      full optimizeAllBranches() sweep after every move; copy D does the same" << endl;
    cerr << "      as C but up to 10 full sweeps. Every candidate is accepted unconditionally" << endl;
    cerr << "      on all four copies, so any divergence in their logL trajectories comes" << endl;
    cerr << "      only from branch-length handling, never from different topological paths." << endl;
    cerr << "      Prints/writes branchlength_compare_data.csv with each copy's logL" << endl;
    cerr << "      before/after/diff per step. C/D are much more expensive per step than" << endl;
    cerr << "      A/B (a full sweep touches every edge, not just 3). EXPERIMENTAL." << endl;
    cerr << endl;
    cerr << "  In the move/list-grafts forms, an edge is a comma-separated leaf name list:" << endl;
    cerr << "    a single leaf, e.g. C         -> that leaf's own pendant edge" << endl;
    cerr << "    two or more leaves, e.g. B,D  -> the internal edge above their MRCA" << endl;
    cerr << endl;
    cerr << "  Examples:" << endl;
    cerr << "    " << prog << " tree.nwk C D                      (leaf onto leaf)" << endl;
    cerr << "    " << prog << " tree.nwk C \"B,D\"                  (leaf onto an internal edge)" << endl;
    cerr << "    " << prog << " tree.nwk \"B,D\" \"E,F\"              (internal edge onto internal edge)" << endl;
    cerr << "    " << prog << " --list-grafts tree.nwk C 3         (list candidates within 3 hops of C)" << endl;
    cerr << "    " << prog << " --likelihood tree.nwk aln.fasta    (evaluate likelihood)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 3 20      (hill-climb search)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 3 20 random   (random Yule-Harding start tree)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 3 20 fast     (O(1)/O(distance) proposal search)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 3 20 fast 5   (5 proposals/step, keep the best)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 3 20 random fast   (both, in either order)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 6 18000 fast quiet (many steps, no per-step spam)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 reopt        (re-optimize branch lengths, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fullreopt 100 5" << endl;
    cerr << "                                                        (full 100-round whole-tree sweep every 5" << endl;
    cerr << "                                                         SUCCESSFUL steps, independent of 'reopt'," << endl;
    cerr << "                                                         experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fullreopt 100 5 true" << endl;
    cerr << "                                                        (same, plus one up-front whole-tree fit" << endl;
    cerr << "                                                         on the starting tree, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fast reopt gtr record" << endl;
    cerr << "                                                        (append convergence trajectory to" << endl;
    cerr << "                                                         record_GTR_FO_fast_reopt.csv, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fast trajectory" << endl;
    cerr << "                                                        (write topology to trajectory_<run-id>.nwk" << endl;
    cerr << "                                                         at the start and after every accepted move," << endl;
    cerr << "                                                         independent of 'record', experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 20 fast investigate" << endl;
    cerr << "                                                        (refine each accepted move one hop" << endl;
    cerr << "                                                         further next step, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 20 fast alternate" << endl;
    cerr << "                                                        (toggle SPR/NNI every other step," << endl;
    cerr << "                                                         experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 5000 fast shrink" << endl;
    cerr << "                                                        (radius narrows on stagnation instead" << endl;
    cerr << "                                                         of a fixed schedule, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 5000 fast shrink 15  (same, stall threshold 15 instead of 10)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 5000 fast learnradius" << endl;
    cerr << "                                                        (radius drawn each step from a Gamma fit to" << endl;
    cerr << "                                                         the last 20 successful moves' own radius," << endl;
    cerr << "                                                         mutually exclusive with shrink, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 5000 fast learnradius 30  (same, window of 30 instead of 20)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 5000 fast distradius learnradius" << endl;
    cerr << "                                                        (fractional percentage drawn each step," << endl;
    cerr << "                                                         e.g. \"radius 4.37\", experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 500 fast sweep" << endl;
    cerr << "                                                        (after 500 fast steps, exhaustive whole-tree" << endl;
    cerr << "                                                         search on the 10 least-compatible sibling" << endl;
    cerr << "                                                         pairs, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 500 fast sweep 20  (same, 20 sibling pairs instead of 10)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 6 500 fast quiet findopt record" << endl;
    cerr << "                                                        (cheap naive-length search, with a" << endl;
    cerr << "                                                         non-destructive scratch refit reading" << endl;
    cerr << "                                                         recorded at the end, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 6 500 fast quiet findopt 100 record" << endl;
    cerr << "                                                        (same, but checked every 100 steps" << endl;
    cerr << "                                                         instead of just once at the end)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 8 20 fast weightprune" << endl;
    cerr << "                                                        (prune edge chosen proportional to its" << endl;
    cerr << "                                                         own branch length, not uniformly)" << endl;
    cerr << "    " << prog << " --branchlength-compare sim.treefile 6 50" << endl;
    cerr << "                                                        (naive vs 4-edge-reopt vs full-sweep-x1 vs" << endl;
    cerr << "                                                         full-sweep-x10 branch lengths, same moves," << endl;
    cerr << "                                                         experimental, slow -- see command's own doc)" << endl;
    cerr << endl;
    cerr << "  Full reference: tree/spr_topology_test_usage.txt" << endl;
}

/**
    parse --hillclimb's trailing optional flags: the literal words
    "random", "fast", and "quiet", in any order, each at most
    once. argv[fromIndex..argc-1]
    must consist of exactly these (in any combination); anything else
    (typos, duplicates, unrelated tokens) is treated as a parse failure so
    main() falls through to printUsage() rather than silently ignoring a
    misspelled flag.

    "fast" may optionally be immediately followed by a positive integer,
    e.g. "fast 5" -- how many independent chooseGraft() candidates to draw
    per step instead of the default 1 (plain "fast", unchanged from
    before). This is a parsing special case, not a separate flag: the
    token right after "fast" is only consumed as its candidate count if it
    actually parses as a positive integer, so "fast quiet" still parses
    "quiet" as its own flag rather than erroring. See numCandidates'
    comment on runHillClimb.

    "reopt" re-optimizes (Newton-Raphson) the 4 edges around each candidate
    (the 3 it actually changes, plus the pruned subtree's own upper edge)
    before it's scored, the same way IQ-TREE's own NNI search does for the
    branches it touches, instead of trusting applySPR's naive placeholder
    lengths (or that upper edge's own pre-move length) -- a bare flag, no
    numeric argument. See scoreTrialSPRMove's comment on runHillClimb.

    "fullreopt" is a fully independent flag from "reopt" (it used to only
    be reachable as "reopt"'s own optional trailing number): it takes TWO
    required positive integers immediately after it, e.g. "fullreopt 100
    5" -- M (fullReoptRounds, the round-count ceiling passed to
    tree.optimizeAllBranches) and N (fullReoptEveryNSteps, how often, in
    SUCCESSFUL/accepted steps, not total step attempts) -- and periodically
    runs a full tree.optimizeAllBranches(M) sweep every N accepted steps,
    on top of whatever "reopt" itself is or isn't doing per move. Unlike
    every other numeric flag in this parser, both numbers are REQUIRED,
    not optional: there's no sensible default for either half of "M
    rounds every N successful steps" the way e.g. "fast" defaults its
    candidate count to 1. See fullReoptEveryNSteps' comment on runHillClimb
    and maybeRunPeriodicFullReopt's comment for exactly how the successful-
    step count is tracked.

    "fullreopt M N" may optionally be followed by a THIRD token, the
    literal word "true" (not a number, so it's never confused with a
    would-be third numeric argument), e.g. "fullreopt 100 5 true" -- turns
    on fullReoptInitialFit, which makes "fullreopt" alone (without
    "reopt") also ML-optimize the whole starting tree once up front,
    exactly like "reopt" always does regardless of this modifier. Omitted,
    "fullreopt" alone leaves the starting tree's branch lengths untouched
    until the first periodic checkpoint. See fullReoptInitialFit's comment
    on runHillClimb.

    "gtr" switches the search's model from JC to GTR+FO; only meaningful
    combined with "reopt" and/or "fullreopt" (JC has no model parameters
    for a GTR fit to replace). See useGtrModel's comment on runHillClimb.

    "record" appends this run's convergence trajectory (candidates
    evaluated, wall time, logL, every time curScore actually changes) to a
    model-specific CSV spreadsheet, record_<model>.csv. See
    recordProgress's comment on runHillClimb and appendRecordRow.

    "investigate" replaces the step immediately after an accepted move
    with an attempt to refine that SAME move one or more real hops
    further, trying again on the step after that if the refinement is
    ALSO accepted, and so on until a refinement attempt fails to improve.
    May optionally be immediately followed by a positive integer, e.g.
    "investigate 3" -- same parsing special case as "fast N"/"shrink N":
    how many real hops (in findGraftPositions' own
    convention, where 1 is the nearest possible legal target -- an NNI)
    the refinement search covers each time. "investigate" alone (N
    omitted) defaults to 1. See investigateFlag's and
    investigateRadius' comments on runHillClimb.

    "alternate" toggles every other step between a plain SPR search at
    <radius> (even step indices) and an NNI-equivalent one forced to
    radius 1 (odd step indices), regardless of which candidate-selection
    path (fast or exhaustive) is active that step; a bare
    flag, no numeric argument. See alternateFlag's comment on runHillClimb.

    "shrink" replaces the step's own radius with a value that only ever
    narrows (starting at <radius>, floor 1), narrowing by 1 once N
    consecutive non-improving steps have piled up (N held constant for
    the whole run). May optionally be immediately followed by a positive
    integer, e.g. "shrink 15" -- same parsing special case as "fast N"/
    "investigate N": the stall threshold. "shrink" alone (N omitted)
    defaults to 10. See shrinkFlag's and maybeShrinkRadius' comments on
    runHillClimb/in the source.

    "learnradius" replaces the step's own radius with one drawn fresh each
    step from a Gamma-vs-starting-distribution mixture fit to the last N
    successfully-accepted moves' own achieved radii -- an adaptive
    alternative to both the fixed <radius> and "shrink"'s one-directional
    narrowing. <radius> is reinterpreted as the MIDDLE of the starting
    spread, not a ceiling: the learned (Gamma) component is free to drift
    anywhere up to a structural maximum (every edge the tree has, or 100%
    of its own branch length under "distradius") once real history
    supports it -- see sampleStartingRadius's and learnRadiusContinuous's
    own comments in the source for the full distribution shape and the
    safeguards (a capped mixture weight, a capped Gamma shape) that keep it
    from ever permanently narrowing onto one value even after it's fully
    converged. Composed with "investigate", each entry fed into that
    history is not any one investigation attempt's own local radius, but
    the NET distance, start to finish, of the whole excursion (the
    initiating move plus every subsequent investigate refinement that kept
    improving) -- see learnradiusFlag's comment on runHillClimb for why.
    Composed with "distradius", draws stay genuinely fractional
    percentages (not rounded to whole numbers before use), and what's fed
    back into the history is the equivalent radiusPercent actually used
    rather than a walk-step/hop count, so the fitted distribution stays in
    one consistent unit -- see learnradiusFlag's own comment for why
    recording the raw hop/step count there would have been wrong. Mutually
    exclusive with "shrink" (checked below,
    the same way randomStart/iqtreeStart are): both exist to replace the
    SAME step's-own-radius value via two different, incompatible
    mechanisms, so giving both is rejected as a parse failure rather than
    silently letting one win. May optionally be immediately followed by a
    positive integer, e.g. "learnradius 30" -- same parsing special case as
    "fast N"/"shrink N": the sliding window size N. "learnradius" alone (N
    omitted) defaults to 20. See learnradiusFlag's comment on runHillClimb
    and learnRadiusContinuous's own comment in the source for the full
    mechanics.

    "sweep" no longer touches the step loop's own selection logic at all --
    it runs as an added phase AFTER every step above has finished, so it
    composes freely with every other flag (including "investigate", no
    longer mutually exclusive with it). It ranks every
    internal edge by how well its current grouping is supported relative to
    its two NNI alternatives (adjacent_subtree_compatibility.pdf, section
    8), then runs an exhaustive, whole-tree regraft search on the N
    LEAST-compatible sibling pairs found, one at a time, keeping each move
    that improves. May optionally be immediately followed by a positive
    integer, e.g. "sweep 20" -- same parsing special case as "fast N"/
    "shrink N": how many of the least-compatible sibling pairs to target.
    "sweep" alone (N omitted) defaults to 10. See sweepFlag's and
    computeSiblingCompatibilityScore's comments on/near runHillClimb for the
    full mechanics and why this is not guaranteed to reach the SPR-optimal
    tree.

    "findopt" (replacing the old, single-shot, curScore-mutating
    "finalreopt") runs a non-destructive, whole-tree ML refit every N steps
    on a throwaway scratch clone of the current tree, ALWAYS under the
    richest available model for the alignment's (auto-detected) sequence
    type -- GTR+FO for DNA, LG+FO for protein -- regardless of whether the
    main search itself is using "gtr" or not (ModelFactory::
    optimizeParameters, jointly refitting branch lengths and the model's
    rate/frequency parameters together) -- the scratch clone is discarded
    once its logL is read off, so it never actually touches the main tree,
    its model, or curScore, only reports what that full refit would find.
    May optionally be immediately followed by a positive
    integer, e.g.
    "findopt 50" -- same parsing special case as "fast N"/"shrink N": N
    omitted defaults to the TOTAL number of steps (maxSteps) -- resolved in
    runHillClimb itself, since maxSteps isn't available yet here -- which
    makes bare "findopt" check exactly once, effectively at the end,
    mirroring "finalreopt"'s old one-shot behavior. See findoptFlag's and
    maybeRunFindopt's comments on/near runHillClimb.

    "notree" means there's no ground-truth tree at all -- the first
    positional argument is a real alignment file directly (any format
    Alignment can auto-detect: FASTA, NEXUS, PHYLIP, ...), not an
    AliSim .treefile with a same-prefix .fa alignment derived from it.
    Disables the final RF-distance-to-true-tree and true-tree-logL
    reference output entirely (see noTrueTree's uses on runHillClimb) --
    everything else about the search (starting tree method, radius,
    selection, model, etc.) is unaffected.

    "distradius" only affects "fast" mode's own candidate draws --
    chooseGraft's directed random walk is replaced with
    chooseGraftByDistance, which reinterprets <radius> as a PERCENTAGE of
    the tree's own total branch length (normalized, so the same value
    means the same thing regardless of this dataset's particular branch-
    length scale) rather than a fixed number of edge hops, and walks by
    actual summed branch length instead of hop count. See
    chooseGraftByDistance's own comment for the full mechanics (why
    sideways/backward hops refund distance rather than spend it, and why
    the walk stops ON the edge that exhausts the budget rather than the
    one before it). Has no effect outside "fast" mode (the exhaustive
    scan's findGraftPositions has no distance-based counterpart). Composes
    with "learnradius" (see its own comment on runHillClimb): draws stay
    genuinely fractional percentages rather than being rounded to whole
    numbers first, and what "learnradius" records into its own history is
    chooseGraftByDistance's own outPercentUsed (an equivalent radiusPercent)
    rather than its outHops walk-step count, so the two stay in the same
    unit.

    "weightprune" replaces choosePrune's uniform random edge pick with one
    weighted by each edge's own branch length -- a bare flag, no numeric
    argument of its own. Independent of every other flag: it only ever
    changes WHICH edge gets pruned, never how the resulting graft search
    (fast or exhaustive, distradius or not, learnradius or not) proceeds
    from there, so it composes freely with all of them. See choosePrune's
    own comment for the mechanics.

    "trajectory" is a bare flag, no numeric argument, independent of every
    other flag here (including "record"/"recordtopology" -- it needs
    neither): it appends the tree's topology to its own
    trajectory_<run-id>.nwk right after the starting tree is built, then
    again after every accepted step, never a reverted one. See
    trajectoryFlag's comment on runHillClimb and appendTrajectoryTopology's
    comment for the mechanics.
    @return false if any trailing argument isn't recognized
 */
bool parseHillClimbFlags(int argc, char **argv, int fromIndex, bool &randomStart, bool &useFastSelection,
        bool &quiet, int &numCandidates, bool &reoptimizeBranchLengths,
        int &fullReoptEveryNSteps, int &fullReoptRounds, bool &fullReoptInitialFit, bool &useGtrModel,
        bool &recordProgress, bool &recordTopology, bool &investigateFlag, int &investigateRadius,
        bool &alternateFlag, bool &shrinkFlag, int &shrinkStallThreshold,
        bool &learnradiusFlag, int &learnradiusN,
        bool &sweepFlag, int &sweepCount,
        bool &findoptFlag, int &findoptEveryNSteps, bool &iqtreeStart, int &iqtreeStartPoolSize,
        bool &noTrueTree, bool &useDistanceRadius, PruneWeighting &weightpruneFlag, bool &trajectoryFlag,
        bool &userStartTree, string &startTreePath,
        bool &useModelOverride, string &modelOverrideSpec) {
    randomStart = false;
    iqtreeStart = false;
    iqtreeStartPoolSize = 20;
    useFastSelection = false;
    quiet = false;
    numCandidates = 1;
    reoptimizeBranchLengths = false;
    fullReoptEveryNSteps = 0;
    fullReoptRounds = 100;
    fullReoptInitialFit = false;
    useGtrModel = false;
    recordProgress = false;
    recordTopology = false;
    investigateFlag = false;
    investigateRadius = 1;
    alternateFlag = false;
    shrinkFlag = false;
    shrinkStallThreshold = 10;
    learnradiusFlag = false;
    learnradiusN = 20;
    sweepFlag = false;
    sweepCount = 10;
    findoptFlag = false;
    findoptEveryNSteps = 0;
    noTrueTree = false;
    useDistanceRadius = false;
    weightpruneFlag = PRUNE_UNIFORM;
    trajectoryFlag = false;
    userStartTree = false;
    startTreePath = "";
    useModelOverride = false;
    modelOverrideSpec = "";
    for (int i = fromIndex; i < argc; i++) {
        string arg = argv[i];
        if (arg == "random" && !randomStart)
            randomStart = true;
        else if (arg == "iqtreestart" && !iqtreeStart) {
            iqtreeStart = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    iqtreeStartPoolSize = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "fast" && !useFastSelection) {
            useFastSelection = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    numCandidates = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "quiet" && !quiet)
            quiet = true;
        else if (arg == "gtr" && !useGtrModel)
            useGtrModel = true;
        else if (arg == "record" && !recordProgress)
            recordProgress = true;
        else if (arg == "recordtopology" && !recordTopology)
            recordTopology = true;
        else if (arg == "alternate" && !alternateFlag)
            alternateFlag = true;
        else if (arg == "sweep" && !sweepFlag) {
            sweepFlag = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    sweepCount = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "shrink" && !shrinkFlag) {
            shrinkFlag = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    shrinkStallThreshold = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "learnradius" && !learnradiusFlag) {
            learnradiusFlag = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    learnradiusN = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "investigate" && !investigateFlag) {
            investigateFlag = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    investigateRadius = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "reopt" && !reoptimizeBranchLengths) {
            reoptimizeBranchLengths = true;
        } else if (arg == "fullreopt" && fullReoptEveryNSteps == 0) {
            // unlike every other numeric flag here, BOTH trailing integers
            // are required -- "M rounds every N steps" has no sensible
            // single-number default the way "fast"/"shrink"/etc. do
            if (i + 2 >= argc)
                return false;
            char *endM = nullptr;
            long m = strtol(argv[i + 1], &endM, 10);
            char *endN = nullptr;
            long n = strtol(argv[i + 2], &endN, 10);
            if (endM == argv[i + 1] || *endM != '\0' || m < 1)
                return false;
            if (endN == argv[i + 2] || *endN != '\0' || n < 1)
                return false;
            fullReoptRounds = (int) m;
            fullReoptEveryNSteps = (int) n;
            i += 2; // consume both numeric arguments
            // optional trailing "true" -- a bare modifier, not a third
            // required number -- turns on the one-time up-front whole-tree
            // fit that reoptimizeBranchLengths always does regardless
            if (i + 1 < argc && string(argv[i + 1]) == "true") {
                fullReoptInitialFit = true;
                i++; // consume the modifier too
            }
        } else if (arg == "findopt" && !findoptFlag) {
            findoptFlag = true;
            if (i + 1 < argc) {
                char *end = nullptr;
                long n = strtol(argv[i + 1], &end, 10);
                if (end != argv[i + 1] && *end == '\0' && n >= 1) {
                    findoptEveryNSteps = (int) n;
                    i++; // consume the numeric argument too
                }
            }
        } else if (arg == "starttree" && !userStartTree) {
            userStartTree = true;
            if (i + 1 >= argc)
                return false; // requires a tree-file-or-Newick-string argument
            startTreePath = argv[i + 1];
            i++; // consume the path argument too
        } else if (arg == "model" && !useModelOverride) {
            useModelOverride = true;
            if (i + 1 >= argc)
                return false; // requires a model spec string
            modelOverrideSpec = argv[i + 1];
            i++; // consume the spec argument too
        } else if (arg == "notree" && !noTrueTree)
            noTrueTree = true;
        else if (arg == "distradius" && !useDistanceRadius)
            useDistanceRadius = true;
        else if (arg == "weightprune" && weightpruneFlag == PRUNE_UNIFORM) {
            // bare "weightprune" keeps its original meaning (bias toward
            // LONG edges); an explicit "long"/"short" after it picks the
            // direction -- see PruneWeighting in tree/sprsearch.h
            weightpruneFlag = PRUNE_LONG;
            if (i + 1 < argc) {
                string next = argv[i + 1];
                if (next == "long" || next == "short") {
                    weightpruneFlag = (next == "short") ? PRUNE_SHORT : PRUNE_LONG;
                    i++;
                }
            }
        }
        else if (arg == "trajectory" && !trajectoryFlag)
            trajectoryFlag = true;
        else
            return false;
    }
    // random and iqtreestart are alternative STARTING-tree methods -- both
    // replace the tool's plain BioNJ estimate with something else, so
    // giving both at once has no well-defined meaning
    if (randomStart && iqtreeStart)
        return false;
    // starttree is a third, mutually exclusive alternative starting-tree
    // method alongside random/iqtreestart (all three replace the default
    // BioNJ estimate)
    if ((randomStart || iqtreeStart) && userStartTree)
        return false;
    // shrink and learnradius are both alternative STEP RADIUS schedules --
    // both replace the fixed <radius> with their own per-step value, via
    // two different, incompatible mechanisms, so giving both at once has no
    // well-defined meaning either
    if (shrinkFlag && learnradiusFlag)
        return false;
    // recordtopology only means anything alongside 'record' -- it adds a
    // companion Newick dump to 'record's own CSV rows, so it has nothing to
    // attach to on its own
    if (recordTopology && !recordProgress)
        return false;
    return true;
}

int main(int argc, char **argv) {
    if (argc == 1)
        return runSelfTest();
    // flag-based modes must be checked before the plain positional "apply
    // one move" form below, since e.g. "--likelihood <tree> <aln>" is also
    // exactly 4 arguments and would otherwise be misread as a move command
    if (argc == 5 && string(argv[1]) == "--list-grafts")
        return runListGrafts(argv[2], argv[3], atoi(argv[4]));
    if (argc == 5 && string(argv[1]) == "--branchlength-compare")
        return runBranchLengthCompare(argv[2], atoi(argv[3]), atoi(argv[4]));
    if (argc >= 5 && string(argv[1]) == "--hillclimb") {
        bool randomStart, useFastSelection, quiet, reoptimizeBranchLengths, fullReoptInitialFit, useGtrModel;
        bool recordProgress, recordTopology, investigateFlag, alternateFlag, shrinkFlag, learnradiusFlag;
        bool sweepFlag, findoptFlag;
        bool iqtreeStart, noTrueTree, useDistanceRadius, trajectoryFlag;
        PruneWeighting weightpruneFlag;
        bool userStartTree;
        string startTreePath;
        bool useModelOverride;
        string modelOverrideSpec;
        int numCandidates, fullReoptEveryNSteps, fullReoptRounds, investigateRadius;
        int shrinkStallThreshold, learnradiusN, sweepCount, findoptEveryNSteps, iqtreeStartPoolSize;
        if (parseHillClimbFlags(argc, argv, 5, randomStart, useFastSelection, quiet, numCandidates,
                reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds, fullReoptInitialFit, useGtrModel,
                recordProgress, recordTopology, investigateFlag, investigateRadius, alternateFlag, shrinkFlag,
                shrinkStallThreshold, learnradiusFlag, learnradiusN, sweepFlag, sweepCount, findoptFlag,
                findoptEveryNSteps, iqtreeStart, iqtreeStartPoolSize, noTrueTree, useDistanceRadius,
                weightpruneFlag, trajectoryFlag, userStartTree, startTreePath, useModelOverride,
                modelOverrideSpec)) {
            return runHillClimb(argv[2], atoi(argv[3]), atoi(argv[4]), randomStart, useFastSelection, quiet,
                    numCandidates, reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds,
                    fullReoptInitialFit, useGtrModel, recordProgress, recordTopology, investigateFlag,
                    investigateRadius, alternateFlag,
                    shrinkFlag, shrinkStallThreshold, learnradiusFlag, learnradiusN,
                    sweepFlag, sweepCount, findoptFlag, findoptEveryNSteps,
                    iqtreeStart, iqtreeStartPoolSize, noTrueTree, useDistanceRadius, weightpruneFlag,
                    trajectoryFlag, userStartTree, startTreePath, useModelOverride, modelOverrideSpec);
        }
    }
    if (argc == 4 && string(argv[1]) == "--likelihood")
        return runLikelihood(argv[2], argv[3]);
    if (argc == 4)
        return runManualSPR(argv[1], argv[2], argv[3]);
    printUsage(argv[0]);
    return 2;
}
