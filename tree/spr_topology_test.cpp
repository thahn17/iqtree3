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
#include "alignment/alignment.h"
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

// defined later (near runBranchLengthCompare, which builds its own B/C/D
// scratch trees the same way); forward-declared here, in the SAME (global)
// scope as that later definition, so maybeRunFindopt below can reuse it too
// without duplicating the model/alignment setup boilerplate -- declaring it
// INSIDE the anonymous namespace below instead would create a second,
// internal-linkage entity distinct from the real one, ambiguous with it at
// every call site from here to the end of the file
void initClonedTree(PhyloTree &t, const string &newickStr, Alignment *aln, Params &params, string modelName);

namespace {

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
string describeEdgeCompact(PhyloNode *node, PhyloNode *dad, size_t maxNames = 4) {
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
    a single legal SPR regraft target: the edge (node,dad), and how many
    hops away from the prune point it lies
 */
struct GraftCandidate {
    PhyloNode *node;
    PhyloNode *dad;
    int radius;
};

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
bool choosePrune(PhyloTree &tree, EdgeRegistry &reg, PhyloNode* &outNode, PhyloNode* &outDad) {
    if (reg.slots.empty())
        return false;
    // bounded retry: for a rooted tree, only edges adjacent to the root
    // can ever fail the orientation check below, a small fraction of the
    // registry, so this succeeds within the first few attempts in
    // practice; the bound just guarantees termination
    for (int attempt = 0; attempt < (int) reg.slots.size(); attempt++) {
        pair<PhyloNode*, PhyloNode*> &edge = reg.slots[random_int((int) reg.slots.size())];
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
    pick a single random SPR regraft target for the (pruneNode,pruneDad)
    edge choosePrune() already selected, via a directed random walk instead
    of enumerating every legal candidate the way findGraftPositions does --
    O(distance walked) rather than O(candidates within radius r).

    First picks a target distance d in [1,r] via chooseGraftDistance (see
    graftDistanceWeight to change how distance is weighted). Then walks d
    steps starting from "the edge above the pruned edge" -- the
    (sibling1,sibling2) edge that would exist once pruneDad is suppressed,
    even though pruneDad hasn't actually been suppressed yet (this mirrors
    findGraftPositions: pruneDad's own edges are never legal targets, so
    this is exactly where a graft-target search actually begins).

    Step 1 is a special case: pick uniformly among every real edge incident
    to sibling1 or sibling2 (excluding their shared edge to pruneDad --
    that's the illegal "start edge" itself, never a candidate at all).

    Steps 2..d: only 3 pointers are ever tracked -- grandDad, dad, node --
    where (dad,node) is the current edge and (grandDad,dad) is the edge the
    walk was on immediately before this one. There's no history beyond
    that single edge: grandDad is overwritten on every move (including
    tier 2's sideways ones), always to whichever node was just displaced,
    so it's always "the edge just come from" for whatever the new current
    edge is -- never a deeper record of how the walk got there. Prefer, in
    order:
      1. node's other edges (walk further outward; up to 2 candidates):
         grandDad=dad, dad=node, node=<pick>
      2. dad's remaining edge, other than the one leading back to grandDad
         (a sideways step onto dad's other, not-yet-explored branch; at
         most 1 candidate, and -- since dad is always degree 3 here --
         effectively always exactly 1): grandDad=node, node=<pick> (dad
         unchanged)
      3. the edge just come from, (grandDad,dad) -- this is not a special
         "undo" move, just the lowest-priority option, tried only when
         neither of the above has anything to offer: dad=grandDad,
         node=dad, grandDad=node (using each variable's value from before
         this reassignment). Only usable when grandDad is actually still
         adjacent to dad, which holds right after step 1 and right after a
         tier-1 move, but not in general right after an earlier tier-3 move
         itself (grandDad then refers to whatever node just got displaced,
         two hops from the new dad, not one) -- so two tier-3 moves never
         fire back to back; that's checked directly (dad->isNeighbor(
         grandDad)) rather than assumed, since without it two-in-a-row
         would fabricate a "current edge" between nodes that were never
         actually adjacent
    Ties within whichever tier is chosen are broken uniformly at random.

    Two things are excluded from ever being selected as dad or node at all,
    at every step: the tree's arbitrary root leaf (isLegalSPR always
    rejects grafting onto its edge, regardless of distance -- and since
    that leaf could be anywhere in the tree, not just near the prune point,
    it has to be filtered out of every tier's candidates, not assumed
    unreachable), and pruneNode itself (tier 3 can reach back through
    pruneDad -- see below -- and must never continue on into the very
    subtree being pruned).

    pruneDad itself can end up as dad mid-walk (reachable via tier 3 right
    after step 1, since grandDad starts out equal to pruneDad there), but
    never as the walk's FINAL position, since that's the one edge no legal
    regraft target can ever be: tier 3 is skipped as the last move if
    grandDad is pruneDad (that would make dad=pruneDad the result), and
    separately, if dad is already pruneDad going into the last step, tier 2
    is skipped too, since it never changes dad -- only tier 1 is guaranteed
    to move dad away from pruneDad, so it's the only option left in that
    situation.

    With all that excluded up front, every candidate this walk can ever
    produce is legal by construction, unlike findGraftPositions which
    instead runs isLegalSPR per candidate after the fact (still ASSERTed
    below anyway, as a cheap defense-in-depth check -- this is exactly how
    an earlier version of this function, which forgot the root-leaf
    exclusion, was caught producing an illegal candidate during testing).

    outDistance, if non-null, receives the walk length d that was chosen
    (not necessarily the true hop-distance of the final edge from the
    prune point -- tier 2/3 moves don't always change hop-distance the way
    a pure tier-1 walk would -- just the number of random-walk steps
    taken, for logging purposes).

    @return false if there is no edge to graft onto at all (sibling1 and
    sibling2 are both leaves, e.g. a 3-leaf tree)
 */
bool chooseGraft(PhyloTree &tree, PhyloNode *pruneNode, PhyloNode *pruneDad, int r,
        PhyloNode* &outNode, PhyloNode* &outDad, int *outDistance = nullptr) {
    int d = chooseGraftDistance(r);
    if (outDistance)
        *outDistance = d;

    // grafting onto the edge incident to the tree's arbitrary root leaf is
    // always illegal (see isLegalSPR's last check) regardless of distance
    // from the prune point; since that leaf could be anywhere in the tree,
    // every candidate-building step below excludes it, exactly as if it
    // had no edge to offer at all (it's a real leaf otherwise, so this is
    // one of two exclusions beyond the usual tier logic -- see pruneNode
    // below for the other)
    PhyloNode *root = (PhyloNode*) tree.root;

    // step 1 (special case): every real edge off sibling1 or sibling2,
    // excluding the edge back to pruneDad
    vector<pair<PhyloNode*, PhyloNode*> > firstStep; // (dad=near, node=far)
    FOR_NEIGHBOR_IT(pruneDad, pruneNode, it) {
        PhyloNode *sibling = (PhyloNode*) (*it)->node;
        FOR_NEIGHBOR_IT(sibling, pruneDad, it2)
            if ((*it2)->node != root)
                firstStep.push_back(make_pair(sibling, (PhyloNode*) (*it2)->node));
    }
    if (firstStep.empty())
        return false;

    pair<PhyloNode*, PhyloNode*> chosen = firstStep[random_int((int) firstStep.size())];
    // grandDad is the other end of the edge the walk was just on -- always
    // pruneDad itself right after step 1, since that's the edge (sibling,
    // pruneDad) the walk left to get here
    PhyloNode *grandDad = pruneDad;
    PhyloNode *dad = chosen.first;
    PhyloNode *node = chosen.second;

    for (int step = 2; step <= d; step++) {
        // pruneDad is excluded here (in addition to root/pruneNode) since
        // it can only ever legitimately become dad, via the tier-3
        // crossing below, which sets dad directly and never consults these
        // lists; letting it slip in here as a tier-1/tier-2 candidate would
        // instead make it `node`, e.g. whenever node currently sits at one
        // of pruneDad's own siblings and dad is one hop further out --
        // `node`'s own neighbors then legitimately include pruneDad -- and
        // node==pruneDad is exactly as illegal as dad==pruneDad
        vector<PhyloNode*> tier1, tier2;
        FOR_NEIGHBOR_IT(node, dad, it)
            if ((*it)->node != root && (*it)->node != pruneNode && (*it)->node != pruneDad)
                tier1.push_back((PhyloNode*) (*it)->node);

        // on the last step, dad must not still be pruneDad afterward (that
        // would make the final answer pruneDad's own edge); tier 1 always
        // moves dad away from pruneDad (dad becomes node, which by
        // construction is never pruneDad), but tier 2 never changes dad at
        // all -- so if dad is already pruneDad, tier 2 must be skipped on
        // the last step, forcing tier 1 (or, failing that, nothing) to fire
        bool mustLeavePruneDad = (step == d) && (dad == pruneDad);
        if (!mustLeavePruneDad) {
            FOR_NEIGHBOR_IT(dad, node, it)
                if ((*it)->node != grandDad && (*it)->node != root && (*it)->node != pruneNode
                        && (*it)->node != pruneDad)
                    tier2.push_back((PhyloNode*) (*it)->node);
        }

        if (!tier1.empty()) {
            grandDad = dad;
            dad = node;
            node = tier1[random_int((int) tier1.size())];
        } else if (!tier2.empty()) {
            grandDad = node;
            node = tier2[random_int((int) tier2.size())];
        } else if (dad->isNeighbor(grandDad) && !(step == d && grandDad == pruneDad)) {
            // tier 3 (the edge just come from): only a real candidate if
            // grandDad is still actually adjacent to dad -- which is true
            // right after step 1, and right after a tier-1 move, but NOT
            // in general right after an earlier tier-3 move itself (that
            // reassigns grandDad to whatever node just got displaced,
            // which sits two hops away from the new dad, not one) --
            // requiring the check here, rather than assuming tier 3 always
            // has a valid destination, is what stops two tier-3 moves in a
            // row from producing a nonexistent "edge" between unrelated
            // nodes (caught by the isLegalSPR ASSERT below during testing).
            // Also skipped as the last move if it would make the result
            // pruneDad's own edge (grandDad==pruneDad would become
            // dad==pruneDad)
            PhyloNode *oldDad = dad, *oldNode = node;
            dad = grandDad;
            node = oldDad;
            grandDad = oldNode;
        }
        // else: no legal move at all this step (only possible in a tiny
        // tree, or when mustLeavePruneDad forbids the only tier that would
        // otherwise fire); nothing to do
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
    applySPRTracked). If reoptimizeBranchLengths is set, the three edges
    an SPR move actually changes -- the two split halves of the target
    edge at the new attachment point (prune_dad-regraft_dad,
    prune_dad-regraft_node) and the merged edge left behind at the
    vacated attachment point (sibling1-sibling2) -- are each re-optimized
    via optimizeOneBranch() (Newton-Raphson, same as
    PhyloTree::getBestNNIForBran does for the branches it touches) before
    the final score is taken, instead of just trusting applySPR's own
    naive placeholder lengths (half the target edge's length split
    evenly, and the sum of the two vacated edges' lengths -- see
    applySPR's own comment in phylotree.cpp). This can only ever improve
    (or leave unchanged) the likelihood for a given topology, since it's
    searching the exact same space optimizeOneBranch always searches, and
    turns the previously-reported score for a topology from "whatever the
    naive placeholder lengths happen to give" into an actual (locally)
    likelihood-optimal set of lengths for those three edges -- closer to
    what a real ML search would report for the same topology.

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
    re-optimize (Newton-Raphson) the 3 edges an SPR move actually changes
    -- prune_dad-regraft_dad, prune_dad-regraft_node, and the merged
    sibling1-sibling2 edge left behind at the vacated attachment point --
    in place, on whatever tree/topology is CURRENTLY applied when this is
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
void reoptimizeSPREdges(PhyloTree &tree, PhyloNode *dad1, PhyloNode *dad2, PhyloNode *node2,
        PhyloNode *sibling1, PhyloNode *sibling2) {
    const int maxNRStep = 10; // matches NNI_MAX_NR_STEP's default, utils/pllnni.cpp
    double minLen = Params::getInstance().min_branch_length;

    ostringstream suppressedOptimizeOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedOptimizeOutput.rdbuf());

    clampBranchLengthForOptimization(dad1, dad2, minLen);
    tree.optimizeOneBranch(dad1, dad2, true, maxNRStep);

    clampBranchLengthForOptimization(dad1, node2, minLen);
    tree.optimizeOneBranch(dad1, node2, true, maxNRStep);

    clampBranchLengthForOptimization(sibling1, sibling2, minLen);
    tree.optimizeOneBranch(sibling1, sibling2, true, maxNRStep);

    cout.rdbuf(realCoutBuf);
}

double scoreTrialSPRMove(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths) {
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    if (reoptimizeBranchLengths)
        findSPRSiblings(move.prune_node, move.prune_dad, sibling1, sibling2);

    SPRRollback rollback;
    tree.applySPR(move, rollback);
    resetLikelihoodBuffers(tree);
    double score = tree.computeLikelihood();

    if (reoptimizeBranchLengths) {
        reoptimizeSPREdges(tree, move.prune_dad, move.regraft_dad, move.regraft_node, sibling1, sibling2);

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
        int shrinkStallThreshold, bool sweepFlag, int sweepCount, int findoptEveryNSteps) {
    time_t now = time(nullptr);
    char timestamp[32];
    strftime(timestamp, sizeof(timestamp), "%Y%m%d-%H%M%S", localtime(&now));

    ostringstream id;
    id << timestamp << "_r" << radius;
    id << "_s" << maxSteps;
    if (randomStart)
        id << "_random";
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
    if (sweepFlag)
        id << "_sweep" << sweepCount;
    if (findoptEveryNSteps > 0)
        id << "_findopt" << findoptEveryNSteps;
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
 */
string buildRecordTag(bool useFastSelection, bool reoptimizeBranchLengths, int fullReoptEveryNSteps,
        bool investigateFlag, int investigateRadius, bool alternateFlag,
        bool shrinkFlag, bool sweepFlag, int sweepCount, int findoptEveryNSteps) {
    ostringstream tag;
    if (useFastSelection)
        tag << "_fast";
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
    if (sweepFlag)
        tag << "_sweep" << sweepCount;
    if (findoptEveryNSteps > 0)
        tag << "_findopt";
    return tag.str();
}

/**
    append one row -- this run's id, how many candidates have been
    evaluated (scoreTrialSPRMove calls) so far, wall-clock seconds elapsed
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

void appendRecordRow(const string &modelName, const string &recordTag, const string &runId,
        long candidatesEvaluated, double timeElapsedSec, double logL, double trueTreeLogl) {
    string path = recordSpreadsheetPath(modelName, recordTag);

    ifstream check(path.c_str());
    bool needsHeader = !check.good() || check.peek() == ifstream::traits_type::eof();
    check.close();

    ofstream out(path.c_str(), ios::app);
    if (needsHeader)
        out << "run_id,candidates,time_elapsed,logL,true_minus_current" << endl;
    out << runId << "," << candidatesEvaluated << "," << timeElapsedSec << ","
        << setprecision(12) << logL << "," << (trueTreeLogl - logL) << endl;
}

/**
    if fullReoptEveryNSteps is set and this step index is due, run one
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
    why. Called once per step in runHillClimb's main loop, independent of
    how that step's own candidate was chosen.
 */
void maybeRunPeriodicFullReopt(PhyloTree &tree, int step, int fullReoptEveryNSteps, int fullReoptRounds,
        bool useGtrModel, bool quiet, bool recordProgress, const string &modelName, const string &recordTag,
        const string &runId, long candidatesEvaluated, double wallClockStart, double trueTreeLogl,
        double &curScore) {
    if (!(fullReoptEveryNSteps > 0 && (step + 1) % fullReoptEveryNSteps == 0))
        return;

    clampAllBranchLengthsForOptimization(tree, Params::getInstance().min_branch_length);
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
        appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getRealTime() - wallClockStart, curScore,
                trueTreeLogl);
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
    counterpart in the actual search's own cost budget, its own wall-clock
    cost is excluded from the run's timing entirely ("pausing the timer"):
    t0 is captured before any of it starts, the elapsed search time up to
    (but not including) this call is what gets recorded for findopt's own
    CSV row, and wallClockStart itself is pushed forward by however long the
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

    The scratch tree always fits under GTR+FO, regardless of whether the
    main search itself is running under useGtrModel (JC) or not -- findopt
    is meant to answer "how good could this topology's branch lengths AND
    model get under the richest model available", not "what would this
    step's own model reach", so it always builds its clone with modelName
    "GTR+FO" and always refits via ModelFactory::optimizeParameters(), never
    plain optimizeAllBranches().

    trueTreeLogl here is NOT necessarily the same value runHillClimb prints
    / records for curScore's own comparisons -- the caller passes
    trueTreeLoglForFindopt, a separate GTR+FO-under-BRLEN_FIX fit of the
    true simulated tree built specifically so it stays model-matched with
    this function's own always-GTR+FO reading (see runHillClimb's own
    comment on trueTreeLoglForFindopt for why the plain trueTreeLogl would
    be an unfair reference whenever the main run isn't already useGtrModel).
 */
void maybeRunFindopt(PhyloTree &tree, int step, int findoptEveryNSteps, bool quiet,
        bool recordProgress, const string &modelName, const string &recordTag, const string &runId,
        long candidatesEvaluated, double &wallClockStart, double trueTreeLogl, double curScore, Alignment *aln,
        Params &params) {
    if (!(findoptEveryNSteps > 0 && (step + 1) % findoptEveryNSteps == 0))
        return;

    double t0 = getRealTime();
    double searchTimeSoFar = t0 - wallClockStart;

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

    // always GTR+FO for the scratch refit, independent of modelName (which
    // reflects the MAIN search's own model, JC unless useGtrModel is set)
    PhyloTree scratchTree;
    initClonedTree(scratchTree, fullPrecisionNewick.str(), aln, params, "GTR+FO");
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
                trueTreeLogl);

    wallClockStart += (getRealTime() - t0); // pause the timer: this scratch pass never counts
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
    load a tree and a DNA alignment, and evaluate the log-likelihood of
    that exact tree (topology and branch lengths as given) against that
    alignment under a plain JC model with no rate heterogeneity. Does NOT
    optimize branch lengths or model parameters -- this reports the
    likelihood of the tree exactly as given, not the best achievable
    likelihood for that topology. The alignment's sequence names must
    match the tree's leaf names exactly (case-sensitive).
    @return 0 on success, 1 if the tree/alignment couldn't be read or
    matched
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
    Alignment *aln = new Alignment((char*) alignmentFile.c_str(), (char*) "DNA", intype, "");

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

    string modelName = "JC";
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);

    tree.initializeAllPartialLh();
    double logl = tree.computeLikelihood();

    cout << "tree          : " << newickOf(tree) << endl;
    cout << "alignment     : " << alignmentFile << " (" << aln->getNSeq() << " sequences, "
         << aln->getNSite() << " sites)" << endl;
    cout << "model         : " << modelName << " (fixed, no branch length or parameter optimization)" << endl;
    cout << "log-likelihood: " << logl << endl;

    delete aln;
    return 2;
}

/**
    greedy, randomized SPR hill-climbing search.

    Takes the tree AliSim used to simulate an alignment (see the AliSim
    command in tree/spr_topology_test_usage.txt for how to produce this
    pair); the alignment file is derived automatically from trueTreeArg by
    AliSim's own default naming convention (<prefix>.treefile paired with
    <prefix>.fa).

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

    quiet (default false) suppresses the one printed line per step (prune
    edge, graft target, logL, kept/reverted -- or the "no candidates"/
    "stopping" messages in their place). This is not just cosmetic: each
    line ends with endl, which flushes -- with max-steps in the thousands,
    an interactive console that's slow to render that much output (a
    classic Windows console/PowerShell window is far more prone to this
    than Windows Terminal or a Unix-style pipe) can end up stalling on
    those flushes, inflating wall-clock time well beyond actual compute
    time. quiet removes the per-step writes entirely; the setup header and
    the final summary (finished/final tree/RF distance/time elapsed) are
    unaffected either way.

    fullReoptEveryNSteps (default 0, disabled) and fullReoptRounds
    (default 100, only meaningful when fullReoptEveryNSteps is actually
    set) together enable a periodic full tree.optimizeAllBranches(
    fullReoptRounds) sweep every fullReoptEveryNSteps loop iterations (by
    step index, regardless of whether that particular
    step's candidate was kept). This is now a fully independent flag from
    reoptimizeBranchLengths -- it used to only be reachable as "reopt"'s
    own optional trailing number, sharing reoptimizeBranchLengths' 3-edge
    per-move reoptimization automatically; the two are now separate
    concerns that compose freely (either alone, both together, or
    neither), since a periodic whole-tree sweep is useful even without
    paying for a 3-edge NR search on every single candidate along the way,
    and vice versa. The idea behind the periodic sweep itself: when
    reoptimizeBranchLengths is ALSO on, its own reoptimizeSPREdges only
    ever touches the 2-3 edges a given SPR move directly changes, so every
    OTHER edge's length is whatever it was left at by the initial sweep
    (or an earlier periodic sweep) -- stale with respect to however much
    the tree's topology has shifted since then; when
    reoptimizeBranchLengths is OFF, branch lengths never change between
    periodic sweeps at all except through them. Either way, a periodic
    full sweep catches that drift at the cost of a whole-tree NR pass
    (roughly 2*(numTaxa-3) edges, several times more expensive than a
    single reoptimizeSPREdges call, since every edge needs re-examining,
    not just the ones a single move touched) instead of 3 edges.
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
    candidates evaluated so far, wall-clock seconds since this run
    started, current logL, and trueTreeLogl - current logL -- how far
    curScore still is below the true AliSim tree's own logL under this
    same model); an initial row (0 candidates, ~setup time, starting
    curScore) is written right after setup so the trajectory has a proper
    starting point. Repeated runs using the same search mode append into
    the same file rather than overwriting it, so their trajectories
    accumulate side by side for comparison; runs using a meaningfully
    different search mode (fast/reopt/investigate) land in a
    separate, appropriately-tagged file instead of being mixed in.

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
    its wall-clock cost is deliberately excluded from the run's own timing
    ("pausing the timer" around it) -- see maybeRunFindopt's comment for
    exactly how.

    @return 0 on success, 1 if the tree/alignment couldn't be read
 */
int runHillClimb(const string &trueTreeArg, int radius, int maxSteps,
        bool randomStart = false, bool useFastSelection = false, bool quiet = false,
        int numCandidates = 1, bool reoptimizeBranchLengths = false, int fullReoptEveryNSteps = 0,
        int fullReoptRounds = 100, bool fullReoptInitialFit = false, bool useGtrModel = false,
        bool recordProgress = false,
        bool investigateFlag = false, int investigateRadius = 1, bool alternateFlag = false,
        bool shrinkFlag = false, int shrinkStallThreshold = 10, bool sweepFlag = false, int sweepCount = 10,
        bool findoptFlag = false, int findoptEveryNSteps = 0) {
    double wallClockStart = getRealTime();
    if (findoptFlag && findoptEveryNSteps <= 0)
        findoptEveryNSteps = maxSteps; // "default = total number of steps"

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

    // the original AliSim tree, kept as a separate plain tree purely for
    // the final RF-distance comparison -- never touched by any SPR move
    PhyloTree trueTree;
    readTreeArg(trueTree, trueTreeArg);
    string trueTreeNewick = newickOf(trueTree);

    Params &params = Params::getInstance();
    params.setDefault();
    // computeBioNJ writes/reads temporary files alongside this prefix
    // (<prefix>.mldist, <prefix>.bionj) as part of how it builds the tree
    params.out_prefix = (char*) "spr_hillclimb_tmp";

    // seeded here, before any random tree generation, so that both a
    // randomStart topology and every step's random prune-edge choice come
    // from the same seeded sequence; time(nullptr) alone has only 1-second
    // resolution, which repeats the exact same "random" sequence across
    // rapid successive runs (e.g. a test script invoking this back to
    // back), so getRealTime()'s sub-second precision is mixed in too
    init_random((int) (time(nullptr) * 1000 + (long) (getRealTime() * 1000) % 1000));

    // the alignment/distance/BioNJ/model setup below goes through several
    // library code paths (Alignment, computeDist, computeBioNJ, ModelFactory)
    // that print their own progress noise (format detection, composition
    // test, distance matrix, RapidNJ progress, ...) unconditionally on cout;
    // none of it is useful for this test tool, so silence cout for the
    // duration of the setup and restore it before printing our own summary
    ostringstream suppressedSetupOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedSetupOutput.rdbuf());

    InputType intype;
    Alignment *aln = new Alignment((char*) alnFile.c_str(), (char*) "DNA", intype, "");

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
    } else {
        // no tree file is read here: computeDist + computeBioNJ build the
        // starting tree structure directly from the alignment's own distance
        // matrix, instead of the readTree-then-setAlignment order used by
        // runLikelihood/runManualSPR
        tree.computeDist(params, aln, tree.dist_matrix, tree.var_matrix);
        tree.computeBioNJ(params);
        // re-map leaf ids to match the alignment's sequence order/names, same
        // as runLikelihood does after reading a tree from file
        tree.setAlignment(aln);
    }

    // built once, O(n); kept in sync thereafter by applySPRTracked/
    // rollbackSPRTracked, never rebuilt -- see EdgeRegistry's comment
    EdgeRegistry edgeRegistry;
    buildEdgeRegistry(tree, edgeRegistry);

    tree.setNumThreads(1);
    if (reoptimizeBranchLengths)
        // real Newton-Raphson branch-length search (see scoreTrialSPRMove)
        // can legitimately push a branch length toward an extreme value
        // while searching, which the plain (non-scaled) SSE kernel isn't
        // built to handle without numerical underflow in the likelihood
        // derivative -- exactly the scenario IQ-TREE's own "-safe" option
        // exists for; the naive fixed-length scoring path never searches
        // branch lengths at all, so it doesn't need this
        params.lk_safe_scaling = true;
    tree.setLikelihoodKernel(LK_SSE2);

    // GTR+FO (ML-estimated rates and frequencies) instead of a plain JC
    // model: sim.fa was simulated under GTR{2,4,1,1,4,2}+F{0.3,0.2,0.2,0.3}
    // (see the AliSim command in this tool's usage doc), so searching under
    // JC -- which has no free rate/frequency parameters to fit at all --
    // is a deliberately misspecified model, not just a simplification. See
    // useGtrModel's comment on runHillClimb for why this matters
    // specifically for whether periodic re-optimization is worth its cost.
    string modelName = useGtrModel ? "GTR+FO" : "JC";
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.initializeAllPartialLh();

    double curScore;
    if (reoptimizeBranchLengths || (fullReoptEveryNSteps > 0 && fullReoptInitialFit)) {
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
    // pointer back to that exact tree, so it can't be shared with `tree`
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
    double trueTreeLogl;
    if (useGtrModel)
        // fit GTR+FO's rate/frequency parameters (only) against the TRUE
        // simulated topology and branch lengths (BRLEN_FIX -- those lengths
        // are exactly as AliSim generated them, not to be touched), so this
        // reference logL reflects a properly-fit model rather than
        // GTR+FO's arbitrary un-optimized starting parameters
        trueTreeLogl = trueTree.getModelFactory()->optimizeParameters(BRLEN_FIX, false, params.modelEps);
    else
        trueTreeLogl = trueTree.computeLikelihood();

    // findopt's own scratch refit (see maybeRunFindopt's comment) always
    // fits under GTR+FO, regardless of the main run's own model -- so
    // comparing its reading against the plain trueTreeLogl above would be
    // unfair whenever the main run ISN'T already under GTR+FO
    // (useGtrModel == false, trueTreeLogl there being an unfit JC
    // likelihood with no free parameters to begin with): any gap would
    // then be partly genuine search progress and partly just a free lunch
    // from a strictly richer model. Build a second, GTR+FO-under-BRLEN_FIX
    // reference off a scratch clone of the SAME true topology and branch
    // lengths (same clone-from-full-precision-Newick approach
    // maybeRunFindopt uses for the main tree, so trueTree's own state --
    // still needed below for computeRFDist -- is never touched), purely so
    // findopt's own comparisons stay apples-to-apples. Skipped entirely
    // when findopt isn't even in use (findoptEveryNSteps == 0), since it'd
    // otherwise be wasted work; trueTreeLogl itself is reused as-is when
    // useGtrModel already made it a GTR+FO fit in the first place.
    double trueTreeLoglForFindopt = trueTreeLogl;
    if (findoptEveryNSteps > 0 && !useGtrModel) {
        int savedPrecisionTrueTree = params.numeric_precision;
        params.numeric_precision = 15;
        ostringstream trueTreeFullPrecisionNewick;
        trueTree.printTree(trueTreeFullPrecisionNewick, WT_BR_LEN);
        params.numeric_precision = savedPrecisionTrueTree;

        PhyloTree trueTreeGtrScratch;
        initClonedTree(trueTreeGtrScratch, trueTreeFullPrecisionNewick.str(), aln, params, "GTR+FO");
        trueTreeLoglForFindopt =
            trueTreeGtrScratch.getModelFactory()->optimizeParameters(BRLEN_FIX, false, params.modelEps);
    }

    cout.rdbuf(realCoutBuf);

    cout << (randomStart ? "random start tree: " : "BioNJ start tree : ")
         << newickOf(tree) << " (logL = " << curScore << ")" << endl;
    cout << "radius          : " << radius << endl;
    cout << "max steps       : " << maxSteps << endl;
    if (useFastSelection)
        cout << "selection       : fast (choosePrune/chooseGraft proposal, not exhaustive)"
             << (numCandidates > 1 ? ", " + to_string(numCandidates) + " candidates/step" : "") << endl;
    if (reoptimizeBranchLengths)
        cout << "branch lengths  : re-optimized (Newton-Raphson, like NNI) on each candidate's 3 "
                "changed edges before scoring, experimental" << endl;
    if (fullReoptEveryNSteps > 0)
        cout << "full reopt      : whole-tree "
             << (useGtrModel ? "optimizeParameters() (model + branch lengths)"
                             : "optimizeAllBranches(" + to_string(fullReoptRounds) + " round(s))")
             << " sweep every " << fullReoptEveryNSteps << " step(s)"
             << (fullReoptInitialFit ? ", plus one up front on the starting tree" : "")
             << ", experimental" << endl;
    if (useGtrModel)
        cout << "model           : GTR+FO (ML-estimated rates/frequencies), experimental" << endl;
    string recordTag = buildRecordTag(useFastSelection, reoptimizeBranchLengths, fullReoptEveryNSteps,
            investigateFlag, investigateRadius, alternateFlag, shrinkFlag, sweepFlag, sweepCount,
            findoptEveryNSteps);
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
    if (sweepFlag)
        cout << "sweep           : after all steps finish, exhaustive whole-tree regraft search on the "
             << sweepCount << " least-compatible sibling pair(s) (adjacent_subtree_compatibility.pdf), "
                "experimental" << endl;
    if (findoptEveryNSteps > 0)
        cout << "findopt         : every " << findoptEveryNSteps << " step(s), scratch whole-tree "
                "optimizeParameters() (GTR+FO model + branch lengths, always, regardless of 'gtr') "
                "refit on a rolled-back copy (main tree unaffected), experimental" << endl;

    string runId = buildRunId(radius, maxSteps, randomStart, useFastSelection,
            numCandidates, reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds, fullReoptInitialFit,
            useGtrModel, investigateFlag, investigateRadius, alternateFlag, shrinkFlag,
            shrinkStallThreshold, sweepFlag, sweepCount, findoptEveryNSteps);
    long candidatesEvaluated = 0;
    if (recordProgress)
        appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getRealTime() - wallClockStart, curScore,
                trueTreeLogl);

    // "investigate": persists across loop iterations (unlike everything
    // else declared inside the loop body) -- investigateNext is whether
    // THIS iteration should refine the previous iteration's accepted move
    // rather than draw a fresh one; investigatePruneNode/Dad is that
    // move's own (node,dad) pair, still valid to re-prune from since
    // applySPR never deletes dad or changes its degree, only repositions
    // it. See investigateFlag's comment on runHillClimb.
    bool investigateNext = false;
    PhyloNode *investigatePruneNode = nullptr, *investigatePruneDad = nullptr;

    // "shrink": persists across iterations like investigateNext above.
    // shrinkCurrentRadius starts at <radius> and only ever decreases (down
    // to a floor of 1), each decrease triggered by maybeShrinkRadius once
    // enough CONSECUTIVE non-improving steps have piled up; see
    // shrinkFlag's comment on runHillClimb.
    int shrinkStallCount = 0;
    int shrinkCurrentRadius = radius;

    int step = 0;
    for (; step < maxSteps; step++) {
        int stepRadius = radius;
        if (shrinkFlag)
            // "shrink" replaces the fixed <radius> with its own
            // stagnation-driven value
            stepRadius = shrinkCurrentRadius;

        bool investigatingThisStep = investigateFlag && investigateNext;
        // default to "not investigating next" -- only re-armed below, and
        // only if THIS step's own candidate is actually accepted; every
        // early skip/rejection path below falls through with this staying
        // false, so investigation naturally stops the first time a
        // refinement attempt doesn't pan out
        if (investigateFlag)
            investigateNext = false;

        // "alternate": every other step forces a radius-1 (NNI-equivalent)
        // search instead of <radius>, regardless of which candidate-
        // selection path (fast or exhaustive) is actually active this
        // step -- investigate takes priority when it's active, since
        // that's testing a specific
        // already-accepted move at investigateRadius, not a step-parity
        // toggle. See alternateFlag's comment on runHillClimb.
        bool nniStepThisTime = alternateFlag && !investigatingThisStep && (step % 2 != 0);
        int effectiveRadius = nniStepThisTime ? 1 : stepRadius;

        PhyloNode *pruneNode, *pruneDad;
        if (investigatingThisStep) {
            pruneNode = investigatePruneNode;
            pruneDad = investigatePruneDad;
        } else if (!choosePrune(tree, edgeRegistry, pruneNode, pruneDad)) {
            if (!quiet)
                cout << "step " << (step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }

        string stepLabel;
        if (investigatingThisStep)
            stepLabel = "(investigate)";
        else if (alternateFlag)
            stepLabel = nniStepThisTime ? "(nni)" : "(spr, radius " + to_string(stepRadius) + ")";
        else
            stepLabel = "(radius " + to_string(stepRadius) + ")";

        PhyloNode *bestNode, *bestDad;
        int bestDistance;
        double bestScore;

        if (useFastSelection && !investigatingThisStep) {
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
            for (int c = 0; c < numCandidates; c++) {
                PhyloNode *candNode, *candDad;
                int walkLength;
                if (!chooseGraft(tree, pruneNode, pruneDad, effectiveRadius, candNode, candDad, &walkLength))
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

                double score = scoreTrialSPRMove(tree, move, reoptimizeBranchLengths);
                candidatesEvaluated++;

                if (!haveBest || score > bestScore) {
                    bestScore = score;
                    bestNode = candNode;
                    bestDad = candDad;
                    bestDistance = walkLength;
                    haveBest = true;
                }
            }
            if (!haveBest) {
                if (!quiet)
                    cout << "step " << (step + 1) << " " << stepLabel << ": prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft target found in " << numCandidates << " draw(s); skipping." << endl;
                continue;
            }
        } else {
            // investigatingThisStep forces investigateRadius (default 1,
            // findGraftPositions' own nearest-legal-target tier -- an
            // NNI-equivalent move) regardless of stepRadius/useFastSelection;
            // otherwise effectiveRadius already reflects "alternate"'s own
            // step-parity NNI toggle, or just equals stepRadius unchanged
            int graftRadius = investigatingThisStep ? investigateRadius : effectiveRadius;
            vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, graftRadius);
            if (candidates.empty()) {
                if (!quiet)
                    cout << "step " << (step + 1) << " " << stepLabel << ": prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft candidates; skipping." << endl;
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

                double score = scoreTrialSPRMove(tree, move, reoptimizeBranchLengths);
                candidatesEvaluated++;

                if (i == 0 || score > bestScore) {
                    bestScore = score;
                    bestCandidate = c;
                }
            }
            bestNode = bestCandidate.node;
            bestDad = bestCandidate.dad;
            bestDistance = bestCandidate.radius;
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
        applySPRTracked(tree, edgeRegistry, bestMove, bestTracked);
        resetLikelihoodBuffers(tree);

        if (reoptimizeBranchLengths) {
            // scoreTrialSPRMove's own reoptimization (used to pick this
            // candidate as the step's best) always gets rolled back along
            // with everything else once scoring is done -- redo it here,
            // on the tree as just permanently applied for real, so the
            // optimized lengths actually persist into the tree this step
            // keeps (or are discarded along with everything else below if
            // this candidate turns out not to improve on curScore after all)
            reoptimizeSPREdges(tree, bestMove.prune_dad, bestMove.regraft_dad, bestMove.regraft_node,
                    bestTracked.sibling1, bestTracked.sibling2);
            resetLikelihoodBuffers(tree);
            double realScore = tree.computeLikelihood();
            if (std::isfinite(realScore))
                bestScore = realScore;
        }

        bool improved = bestScore > curScore;
        if (shrinkFlag)
            maybeShrinkRadius(improved, shrinkStallThreshold, quiet,
                    shrinkStallCount, shrinkCurrentRadius);
        if (!quiet)
            cout << "step " << (step + 1) << " " << stepLabel << ": prune {"
                 << describeEdgeCompact(pruneNode, pruneDad) << "}"
                 << " -> graft {" << describeEdgeCompact(bestNode, bestDad) << "}"
                 << " (" << ((useFastSelection && !investigatingThisStep) ? "d=" : "distance ") << bestDistance << ")"
                 << ", logL " << bestScore << " (cur " << curScore << ")"
                 << (improved ? " [kept]" : " [reverted]") << endl;

        if (improved) {
            curScore = bestScore;
            if (investigateFlag) {
                // this move -- whether it came from a fresh choosePrune or
                // from investigating a previous one -- just improved the
                // tree, so refine IT one real hop further next step; see
                // investigateFlag's comment on runHillClimb
                investigateNext = true;
                investigatePruneNode = pruneNode;
                investigatePruneDad = pruneDad;
            }
            if (recordProgress)
                appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getRealTime() - wallClockStart,
                        curScore, trueTreeLogl);
        } else {
            rollbackSPRTracked(tree, edgeRegistry, bestTracked);
            resetLikelihoodBuffers(tree);
        }

        maybeRunPeriodicFullReopt(tree, step, fullReoptEveryNSteps, fullReoptRounds, useGtrModel, quiet,
                recordProgress, modelName, recordTag, runId, candidatesEvaluated, wallClockStart, trueTreeLogl,
                curScore);
        maybeRunFindopt(tree, step, findoptEveryNSteps, quiet, recordProgress, modelName, recordTag,
                runId, candidatesEvaluated, wallClockStart, trueTreeLoglForFindopt, curScore, aln, params);
    }

    // "sweep": post-processing phase, run strictly AFTER the step loop
    // above has finished (with whatever mix of fast/investigate/
    // alternate/shrink shaped those steps) -- see sweepFlag's comment on
    // runHillClimb for the full rationale. Unlike every flag above, this
    // never competes for a STEP's own selection logic, so it needs no
    // incompatibility check with any of them.
    if (sweepFlag) {
        // Step 1: score every internal edge (both endpoints degree 3 -- an
        // edge to a leaf has no quartet to test, see
        // computeSiblingCompatibilityScore) by how well its CURRENT
        // grouping is supported relative to the two alternative NNI
        // rearrangements around it. Lower (more negative) means less
        // compatible -- these are the "least-compatible siblings" sweepCount
        // targets.
        vector<pair<double, pair<PhyloNode*, PhyloNode*> > > ranked;
        for (size_t slot = 0; slot < edgeRegistry.slots.size(); slot++) {
            PhyloNode *p = edgeRegistry.slots[slot].first;
            PhyloNode *q = edgeRegistry.slots[slot].second;
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

        int sweepTargets = min(sweepCount, (int) ranked.size());
        if (!quiet)
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
        int sweepFullRadius = (int) edgeRegistry.slots.size();
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
                if (!quiet)
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

                double score = scoreTrialSPRMove(tree, move, reoptimizeBranchLengths);
                candidatesEvaluated++;

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
            applySPRTracked(tree, edgeRegistry, bestMove, bestTracked);
            resetLikelihoodBuffers(tree);

            if (reoptimizeBranchLengths) {
                reoptimizeSPREdges(tree, bestMove.prune_dad, bestMove.regraft_dad, bestMove.regraft_node,
                        bestTracked.sibling1, bestTracked.sibling2);
                resetLikelihoodBuffers(tree);
                double realScore = tree.computeLikelihood();
                if (std::isfinite(realScore))
                    bestScore = realScore;
            }

            bool improvedSweep = bestScore > curScore;
            if (!quiet)
                cout << sweepLabel << ": prune {" << describeEdgeCompact(pruneNode, pruneDad) << "}"
                     << " -> graft {" << describeEdgeCompact(bestCandidate.node, bestCandidate.dad) << "}"
                     << " (distance " << bestCandidate.radius << ")"
                     << ", logL " << bestScore << " (cur " << curScore << ")"
                     << (improvedSweep ? " [kept]" : " [reverted]") << endl;

            if (improvedSweep) {
                curScore = bestScore;
                if (recordProgress)
                    appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getRealTime() - wallClockStart,
                            curScore, trueTreeLogl);
            } else {
                rollbackSPRTracked(tree, edgeRegistry, bestTracked);
                resetLikelihoodBuffers(tree);
            }
        }
    }

    // unconditional final row: the loop above only calls appendRecordRow on
    // an ACCEPTED step (or a periodic full-reopt event), and the sweep phase
    // above only on an accepted sweep position -- so if the very last event
    // (whichever phase it came from) was rejected, the CSV's last row would
    // otherwise be stale relative to when the run actually stopped -- this
    // guarantees one row always reflects the true final state (curScore as
    // it stands once EVERYTHING above, including sweep, is done), regardless
    // of whether that last event was ever kept. findopt is NOT part of this
    // guarantee -- unlike the old, curScore-mutating "finalreopt" this
    // replaced, findopt never changes curScore at all, so it has no result
    // for this final row to reflect; its own periodic checkpoints (inside
    // the loop above) each write their own separate row instead.
    if (recordProgress)
        appendRecordRow(modelName, recordTag, runId, candidatesEvaluated, getRealTime() - wallClockStart, curScore,
                trueTreeLogl);

    if (!quiet)
        cout << endl;
    cout << "=== finished after " << step << " step(s) ===" << endl;
    string finalTreeNewick = newickOf(tree);
    cout << "final tree (logL = " << curScore << "): " << finalTreeNewick << endl;
    cout << "AliSim true tree logL: " << trueTreeLogl << endl;

    stringstream finalTreeStream;
    finalTreeStream << finalTreeNewick;
    finalTreeStream.seekg(0, ios::beg);
    vector<double> rfdist;
    trueTree.computeRFDist(finalTreeStream, rfdist);
    int rf = rfdist.empty() ? -1 : (int) rfdist[0];

    cout << "RF distance to the original AliSim tree: " << rf << endl;

    ofstream out("output.txt");
    if (out.good()) {
        out << "AliSim true tree : " << trueTreeNewick << endl;
        out << "Final result tree: " << finalTreeNewick << endl;
        out << "RF distance      : " << rf << endl;
        out.close();
        cout << "Results written to output.txt" << endl;
    } else {
        cerr << "warning: could not write to output.txt" << endl;
    }

    cout << "time elapsed    : " << fixed << setprecision(2) << (getRealTime() - wallClockStart)
         << " sec" << endl;

    delete aln;
    return 2;
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
    lengths (see below), same alignment, same JC model, and -- critically
    -- the same move at every step (same prune edge, same regraft
    target), so any divergence between their logL trajectories can only
    come from how branch lengths were set, never from the four trees
    taking different topological paths.

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
    the 3 an SPR move changes), so this command is much slower than a
    plain "reopt" hill-climb at the same step count. The final summary
    reports each tree's own average apply+score time per step (excludes
    the prune/graft selection work all four trees pay alike), so the
    relative cost of naive vs. 3-edge reopt vs. 1 vs. 10 full sweeps is
    visible directly, not just inferred from total wall time.

    @return 0 on success, 1 if the trees desynchronize (a bug, should
    never happen in practice -- see point 3 above), 2 if the alignment
    couldn't be found
 */
int runBranchLengthCompare(const string &trueTreeArg, int radius, int maxSteps) {
    double wallClockStart = getRealTime();

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
    Alignment *aln = new Alignment((char*) alnFile.c_str(), (char*) "DNA", intype, "");

    string modelName = "JC";
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
    cout << "tree B           : re-optimized (Newton-Raphson) on the 3 changed edges every step" << endl;
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
    // A/B's naive/3-edge handling
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
        double t0 = getRealTime();
        applySPRTracked(treeA, edgeRegistryA, moveA, trackedA);
        resetLikelihoodBuffers(treeA);
        double newLoglA = treeA.computeLikelihood();
        timeA += getRealTime() - t0;

        // reoptimizeSPREdges'/optimizeAllBranches' optimizeOneBranch calls
        // need valid partial likelihood buffers for the NEW topology
        // already in place before they run -- see scoreTrialSPRMove's own
        // identical reset-then-reoptimize order, which this mirrors for
        // B, C, and D alike
        t0 = getRealTime();
        applySPRTracked(treeB, edgeRegistryB, moveB, trackedB);
        resetLikelihoodBuffers(treeB);
        reoptimizeSPREdges(treeB, moveB.prune_dad, moveB.regraft_dad, moveB.regraft_node, sibling1B, sibling2B);
        resetLikelihoodBuffers(treeB);
        double newLoglB = treeB.computeLikelihood();
        timeB += getRealTime() - t0;

        t0 = getRealTime();
        applySPRTracked(treeC, edgeRegistryC, moveC, trackedC);
        resetLikelihoodBuffers(treeC);
        double newLoglC = treeC.optimizeAllBranches(1);
        timeC += getRealTime() - t0;

        t0 = getRealTime();
        applySPRTracked(treeD, edgeRegistryD, moveD, trackedD);
        resetLikelihoodBuffers(treeD);
        double newLoglD = treeD.optimizeAllBranches(10);
        timeD += getRealTime() - t0;
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
    cout << "final tree B (3-edge reopt, logL = " << curLoglB << ")" << endl;
    cout << "final tree C (1 full sweep/step, logL = " << curLoglC << ")" << endl;
    cout << "final tree D (10 full sweeps/step, logL = " << curLoglD << ")" << endl;
    cout << "per-step data written to branchlength_compare_data.csv" << endl;
    if (measuredSteps > 0)
        cout << "avg time/step    : A=" << fixed << setprecision(1) << (timeA / measuredSteps * 1000.0)
             << "ms  B=" << (timeB / measuredSteps * 1000.0) << "ms  C=" << (timeC / measuredSteps * 1000.0)
             << "ms  D=" << (timeD / measuredSteps * 1000.0) << "ms  (" << measuredSteps << " step(s) measured)"
             << endl;
    cout << "time elapsed     : " << fixed << setprecision(2) << (getRealTime() - wallClockStart)
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
    cerr << "      lengths as given, no optimization) against a DNA alignment under a" << endl;
    cerr << "      plain JC model. Sequence names in the alignment must match the tree's" << endl;
    cerr << "      leaf names exactly." << endl;
    cerr << endl;
    cerr << "  " << prog << " --hillclimb <alisim-tree.treefile> <radius> <max-steps> [random] [fast [N]] [quiet] [reopt] [fullreopt M N] [gtr] [record] [investigate [N]] [alternate] [shrink [N]] [sweep [N]] [findopt [N]]" << endl;
    cerr << "      greedy randomized SPR search: build a BioNJ start tree from the" << endl;
    cerr << "      alignment AliSim simulated from <alisim-tree.treefile> (found by" << endl;
    cerr << "      replacing '.treefile' with '.fa'), then repeatedly prune a random edge," << endl;
    cerr << "      evaluate every legal regraft within <radius> hops via applySPR/" << endl;
    cerr << "      rollbackSPR on one tree object, and keep the best if it improves the" << endl;
    cerr << "      likelihood, for up to <max-steps> rounds. Prints the RF distance to the" << endl;
    cerr << "      original AliSim tree and writes both trees + the RF distance to" << endl;
    cerr << "      output.txt. Twelve optional trailing flags, in any order:" << endl;
    cerr << "        random     start from a random Yule-Harding topology instead of the" << endl;
    cerr << "                   default BioNJ estimate tree" << endl;
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
    cerr << "        reopt      re-optimize (Newton-Raphson) the 3 edges an SPR move actually" << endl;
    cerr << "                   changes before scoring a candidate, the same way IQ-TREE's own" << endl;
    cerr << "                   NNI search re-optimizes the branches it touches -- instead of" << endl;
    cerr << "                   trusting applySPR's naive placeholder lengths (half the target" << endl;
    cerr << "                   edge split evenly, the two vacated edges summed). When reopt (or" << endl;
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
    cerr << "                   optimizeAllBranches(M) sweep over every edge every N steps, on top" << endl;
    cerr << "                   of whatever 'reopt' itself is or isn't doing per candidate. BOTH M" << endl;
    cerr << "                   (round-count ceiling) and N (step interval) are REQUIRED -- unlike" << endl;
    cerr << "                   every other numeric flag here, there's no sensible single-number" << endl;
    cerr << "                   default for 'M rounds every N steps'. (mirrors IQ-TREE's own NNI" << endl;
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
    cerr << "                   JC -- relevant since sim.fa is simulated under a real GTR+F" << endl;
    cerr << "                   model, so JC is a genuine misspecification, not just a" << endl;
    cerr << "                   simplification. Combined with 'fullreopt M N', periodic sweeps" << endl;
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
    cerr << "                   evaluated so far, wall-clock seconds elapsed, the current logL," << endl;
    cerr << "                   and how far that logL still is below the true AliSim tree's own" << endl;
    cerr << "                   logL under this model (true_minus_current)." << endl;
    cerr << "                   Repeated runs using the same search mode APPEND (never" << endl;
    cerr << "                   overwrite) so their trajectories accumulate side by side in the" << endl;
    cerr << "                   same file for later comparison -- see appendRecordRow's and" << endl;
    cerr << "                   buildRecordTag's comments in the source" << endl;
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
    cerr << "                   clone of the current tree, ALWAYS under GTR+FO regardless of whether" << endl;
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
    cerr << "                   no counterpart in the actual search's own cost, its wall-clock cost" << endl;
    cerr << "                   is excluded from the run's own timing entirely ('pausing the timer'" << endl;
    cerr << "                   around it, both for its own CSV row and every later one)." << endl;
    cerr << "                   EXPERIMENTAL -- see maybeRunFindopt's comment in the source" << endl;
    cerr << endl;
    cerr << "  " << prog << " --branchlength-compare <alisim-tree.treefile> <radius> <max-steps>" << endl;
    cerr << "      NOT a search: applies the SAME sequence of random SPR moves to FOUR" << endl;
    cerr << "      separate copies of one starting tree (see runBranchLengthCompare's own" << endl;
    cerr << "      comment for how the copies are kept in lockstep without reusing" << endl;
    cerr << "      choosePrune()/chooseGraft(), which would desync them). Copy A keeps" << endl;
    cerr << "      applySPR's own naive placeholder branch lengths every move (the default" << endl;
    cerr << "      everywhere else in this tool); copy B re-optimizes (Newton-Raphson) just" << endl;
    cerr << "      the 3 changed edges after every move; copy C re-optimizes EVERY edge, one" << endl;
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
    cerr << "                                                         steps, independent of 'reopt', experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fullreopt 100 5 true" << endl;
    cerr << "                                                        (same, plus one up-front whole-tree fit" << endl;
    cerr << "                                                         on the starting tree, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 fast reopt gtr record" << endl;
    cerr << "                                                        (append convergence trajectory to" << endl;
    cerr << "                                                         record_GTR_FO_fast_reopt.csv, experimental)" << endl;
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
    cerr << "    " << prog << " --branchlength-compare sim.treefile 6 50" << endl;
    cerr << "                                                        (naive vs 3-edge-reopt vs full-sweep-x1 vs" << endl;
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

    "reopt" re-optimizes (Newton-Raphson) each candidate's 3 changed edges
    before it's scored, the same way IQ-TREE's own NNI search does for the
    branches it touches, instead of trusting applySPR's naive placeholder
    lengths -- a bare flag, no numeric argument. See scoreTrialSPRMove's
    comment on runHillClimb.

    "fullreopt" is a fully independent flag from "reopt" (it used to only
    be reachable as "reopt"'s own optional trailing number): it takes TWO
    required positive integers immediately after it, e.g. "fullreopt 100
    5" -- M (fullReoptRounds, the round-count ceiling passed to
    tree.optimizeAllBranches) and N (fullReoptEveryNSteps, how often, in
    steps) -- and periodically runs a full tree.optimizeAllBranches(M)
    sweep every N steps, on top of whatever "reopt" itself is or isn't
    doing per move. Unlike every other numeric flag in this parser, both
    numbers are REQUIRED, not optional: there's no sensible default for
    either half of "M rounds every N steps" the way e.g. "fast" defaults
    its candidate count to 1. See fullReoptEveryNSteps' comment on
    runHillClimb.

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
    on a throwaway scratch clone of the current tree, ALWAYS under GTR+FO
    regardless of whether the main search itself is using "gtr" or not
    (ModelFactory::optimizeParameters, jointly refitting branch lengths and
    the model's rate/frequency parameters together) -- the scratch clone is
    discarded once its logL is read off, so it never actually touches the
    main tree, its model, or curScore, only reports what a full GTR+FO
    refit would find. May optionally be immediately followed by a positive
    integer, e.g.
    "findopt 50" -- same parsing special case as "fast N"/"shrink N": N
    omitted defaults to the TOTAL number of steps (maxSteps) -- resolved in
    runHillClimb itself, since maxSteps isn't available yet here -- which
    makes bare "findopt" check exactly once, effectively at the end,
    mirroring "finalreopt"'s old one-shot behavior. See findoptFlag's and
    maybeRunFindopt's comments on/near runHillClimb.
    @return false if any trailing argument isn't recognized
 */
bool parseHillClimbFlags(int argc, char **argv, int fromIndex, bool &randomStart, bool &useFastSelection,
        bool &quiet, int &numCandidates, bool &reoptimizeBranchLengths,
        int &fullReoptEveryNSteps, int &fullReoptRounds, bool &fullReoptInitialFit, bool &useGtrModel,
        bool &recordProgress, bool &investigateFlag, int &investigateRadius,
        bool &alternateFlag, bool &shrinkFlag, int &shrinkStallThreshold, bool &sweepFlag, int &sweepCount,
        bool &findoptFlag, int &findoptEveryNSteps) {
    randomStart = false;
    useFastSelection = false;
    quiet = false;
    numCandidates = 1;
    reoptimizeBranchLengths = false;
    fullReoptEveryNSteps = 0;
    fullReoptRounds = 100;
    fullReoptInitialFit = false;
    useGtrModel = false;
    recordProgress = false;
    investigateFlag = false;
    investigateRadius = 1;
    alternateFlag = false;
    shrinkFlag = false;
    shrinkStallThreshold = 10;
    sweepFlag = false;
    sweepCount = 10;
    findoptFlag = false;
    findoptEveryNSteps = 0;
    for (int i = fromIndex; i < argc; i++) {
        string arg = argv[i];
        if (arg == "random" && !randomStart)
            randomStart = true;
        else if (arg == "fast" && !useFastSelection) {
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
        } else
            return false;
    }
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
        bool recordProgress, investigateFlag, alternateFlag, shrinkFlag, sweepFlag, findoptFlag;
        int numCandidates, fullReoptEveryNSteps, fullReoptRounds, investigateRadius;
        int shrinkStallThreshold, sweepCount, findoptEveryNSteps;
        if (parseHillClimbFlags(argc, argv, 5, randomStart, useFastSelection, quiet, numCandidates,
                reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds, fullReoptInitialFit, useGtrModel,
                recordProgress, investigateFlag, investigateRadius, alternateFlag, shrinkFlag,
                shrinkStallThreshold, sweepFlag, sweepCount, findoptFlag, findoptEveryNSteps)) {
            return runHillClimb(argv[2], atoi(argv[3]), atoi(argv[4]), randomStart, useFastSelection, quiet,
                    numCandidates, reoptimizeBranchLengths, fullReoptEveryNSteps, fullReoptRounds,
                    fullReoptInitialFit, useGtrModel, recordProgress, investigateFlag,
                    investigateRadius, alternateFlag,
                    shrinkFlag, shrinkStallThreshold, sweepFlag, sweepCount, findoptFlag, findoptEveryNSteps);
        }
    }
    if (argc == 4 && string(argv[1]) == "--likelihood")
        return runLikelihood(argv[2], argv[3]);
    if (argc == 4)
        return runManualSPR(argv[1], argv[2], argv[3]);
    printUsage(argv[0]);
    return 2;
}
