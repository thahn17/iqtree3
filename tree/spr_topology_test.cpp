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
    full-precision Newick text for the subtree hanging off `node`, away from
    `dad` (dad may be nullptr, meaning "no direction is excluded" -- correct
    when node is tree.root itself, which has only one neighbor to begin
    with). Unlike newickOf() elsewhere in this file, branch lengths are NOT
    rounded to 1 decimal place -- that rounding is a display-only nicety and
    would bias the actual likelihood computation this text feeds.
 */
void writeFullPrecisionNewick(ostream &out, PhyloNode *node, PhyloNode *dad) {
    if (node->isLeaf()) {
        out << node->name;
        return;
    }
    out << "(";
    bool first = true;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!first)
            out << ",";
        first = false;
        writeFullPrecisionNewick(out, (PhyloNode*) (*it)->node, node);
        out << ":" << (*it)->length;
    }
    out << ")";
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
    edges of (pruneNode,pruneDad), where radius is measured in tree hops
    from pruneDad (pruneDad's own two remaining edges are radius 1, their
    neighbors' edges are radius 2, and so on).

    Implementation: a breadth-first search rooted at pruneDad, never
    crossing into pruneNode's subtree (so no candidate can ever lie inside
    the subtree being pruned), recording one candidate per edge the first
    time it is reached. A tree traversal can only ever visit each edge
    once, so duplicates are structurally impossible here -- but every edge
    is additionally checked against a canonical (min id, max id) key in
    `seenEdges` before being recorded, as an explicit, visible guarantee
    rather than an implicit one. Every raw candidate is then run through
    isLegalSPR (the same legality check applySPR itself asserts), so only
    genuinely legal candidates -- e.g. never one of pruneDad's own original
    edges, which isLegalSPR always rejects as incident to the node being
    suppressed -- are returned.
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
    q.push({pruneDad, pruneNode, 0});

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
    just a coin flip -- tree.root is only this framework's arbitrary
    bookkeeping leaf there (its own edge is already excluded from `reg` by
    buildEdgeRegistry), so it doesn't matter which side it ends up on.

    For a genuinely rooted tree (tree.rooted), it does matter: tree.root is
    a real topological root, and it must always end up on the "dad" (the
    tree that remains) side, never the "node" (pruned branch) side --
    otherwise the move would effectively prune away the tree's own root.
    subtreeContainsRoot decides which side that is; if degree alone would
    force the opposite orientation (root's side isn't the degree-3 one),
    this edge has no valid direction at all and is skipped, not returned.

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
            } else if (random_int(2) == 0) {
                outDad = edge.first; outNode = edge.second;
            } else {
                outDad = edge.second; outNode = edge.first;
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
    invalidate cached partial parsimony vectors after a topology change, for
    exactly the same reason resetLikelihoodBuffers exists: parsimony
    partials use the identical per-Neighbor "already computed" bit-flag
    caching as likelihood partials (see computePartialParsimonyFast's
    `partial_lh_computed & 2` check in phylotreepars.cpp), so the same
    stale-buffer risk after a long-distance SPR jump applies here too.

    tree.initializeAllPartialPars() alone is sufficient (unlike the
    likelihood case, no separate delete step is needed): it reallocates
    central_partial_pars only if resetLikelihoodBuffers has already freed
    it (deleteAllPartialLh frees central_partial_pars too, since the two
    buffer types share that lifecycle), re-derives every Neighbor's
    partial_pars pointer either way, and its own trailing call to
    clearAllPartialLH() resets the "already computed" flag (both the
    likelihood and parsimony bits) on every Neighbor -- exactly the
    invalidation this needs.
 */
void resetParsimonyBuffers(PhyloTree &tree) {
    tree.initializeAllPartialPars();
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

double scoreTrialSPRMove(PhyloTree &tree, const SPRMove &move, bool reoptimizeBranchLengths) {
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    if (reoptimizeBranchLengths)
        findSPRSiblings(move.prune_node, move.prune_dad, sibling1, sibling2);

    SPRRollback rollback;
    tree.applySPR(move, rollback);
    resetLikelihoodBuffers(tree);
    double score = tree.computeLikelihood();

    if (reoptimizeBranchLengths) {
        const int maxNRStep = 10; // matches NNI_MAX_NR_STEP's default, utils/pllnni.cpp
        double minLen = Params::getInstance().min_branch_length;

        clampBranchLengthForOptimization(move.prune_dad, move.regraft_dad, minLen);
        tree.optimizeOneBranch(move.prune_dad, move.regraft_dad, true, maxNRStep);

        clampBranchLengthForOptimization(move.prune_dad, move.regraft_node, minLen);
        tree.optimizeOneBranch(move.prune_dad, move.regraft_node, true, maxNRStep);

        clampBranchLengthForOptimization(sibling1, sibling2, minLen);
        tree.optimizeOneBranch(sibling1, sibling2, true, maxNRStep);

        resetLikelihoodBuffers(tree);
        score = tree.computeLikelihood();
    }

    tree.rollbackSPR(rollback);
    resetLikelihoodBuffers(tree);
    return score;
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
// judged not worth pursuing further here. See runRandomWalkSubtree's own
// comment for the (correct, if slower) standalone-tree method actually
// used, and its measured efficiency numbers.

/**
    how many extra parsimony substitutions a candidate is allowed to have,
    beyond the CURRENT tree's own parsimony score, before it gets rejected
    outright by the prescreen (see usePrescreen on runHillClimb) without
    ever running the real, expensive likelihood check on it.

    EXPERIMENTAL -- two rounds of testing against sim.treefile (~100
    taxa/3000bp, JC model, exhaustive mode), and they don't fully agree,
    which is itself the most honest thing to report here.

    Round 1: tolerance = 0, 1, 2, 3, 5, 8, 15 (radius 6, 50 steps) and
    0, 2, 5, 10, 20 (radius 15, 20 steps, repeated a few times). RF
    distance showed no consistent trend with tolerance (noisy, ~50-76
    throughout). Wall-clock time was consistently WORSE with prescreen on
    than off at every tolerance tried (e.g. radius 15: ~23s with no
    prescreen vs ~26-32s with it).

    Round 2 (radius 10, 30 steps, 3 repeats per value, averaged):
        no prescreen : avg 22.5s  (16.5-27.4s), avg RF 60.7
        tolerance 0  : avg 25.5s  (24.4-26.1s), avg RF 67.3
        tolerance 2  : avg 16.9s  (14.9-18.5s), avg RF 62.0
        tolerance 5  : avg 24.6s  (22.4-26.8s), avg RF 65.3
        tolerance 10 : avg 24.7s  (23.1-25.5s), avg RF 64.0
        tolerance 20 : avg 22.4s  (20.0-26.1s), avg RF 62.7
        tolerance 40 : avg 26.6s  (22.1-30.8s), avg RF 65.3
    Here tolerance 2 was the clear outlier: consistently ~25% faster than
    no-prescreen across all 3 repeats (its slowest rep was still faster
    than every no-prescreen rep but one), for an RF cost of about +1.3
    over baseline -- comfortably inside this search's own run-to-run RF
    noise (~10-15 points). No other tolerance tested, including 0, showed
    a similar speed advantage; most were flat or slower than baseline.

    Given round 1 found no benefit at any tolerance and round 2 found a
    real-looking benefit specifically at 2 (with only 3 repeats per value
    and no fixed RNG seed, so this could still be a lucky draw of which
    prune edges got picked), 2 is kept as the shipped default. It is the
    only value with any empirical support across two rounds of testing,
    but "3 repeats showed a 25% speedup once" is not strong evidence --
    treat this as a reasonable starting point, not a settled result.
    Retest with more repeats (and ideally a fixed seed) before relying on
    it for anything real.

    IMPORTANT CORRECTION, discovered while building runRandomWalk():
    computeParsimony() silently returned 0 for every tree/candidate for
    this entire file's history up to that point, because nsites in
    computeParsimonyBranchFast() is derived from aln->num_parsimony_sites,
    which stays 0 until Alignment::orderPatternByNumChars() is called --
    a call IQ-TREE's real pipeline makes once per alignment (see
    IQTree::doTreeSearch, iqtree.cpp) but this tool never did, for either
    runHillClimb's prescreen or the code that produced the numbers above.
    That means every one of the round 1/round 2 runs above was comparing
    "0 > 0 + tolerance", which is false for every non-negative tolerance
    -- prescreen never rejected a single candidate at any tolerance ever
    tested. All those wall-clock/RF differences were pure RNG noise from
    the search itself, not any effect of tolerance or of prescreen being
    on vs off. This is now fixed (both here and in runHillClimb, which
    now calls aln->orderPatternByNumChars(PAT_VARIANT) once before first
    computing parsimony) -- computeParsimony() returns real, large,
    sensibly-varying scores now. This makes prescreen an actual filter
    for the first time, and means the round 1/round 2 conclusions above
    should be discounted entirely, not treated as evidence for or against
    any tolerance value. See runRandomWalk() for a fresh, unbiased look at
    how parsimony difference actually relates to likelihood/RF difference
    now that the bug is fixed; parsimonyPrescreenTolerance() itself has
    NOT been retested against real (non-zero) parsimony scores yet.
 */
int parsimonyPrescreenTolerance() {
    return 2;
}

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
         radius (findGraftPositions) -- see radiusDecayPerStep below for how
         that can shrink from `radius` as steps progress
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

    radiusDecayPerStep (default 0.0, i.e. constant radius -- the original
    --hillclimb behavior) lets the radius shrink as the search progresses:
    step i (0-indexed) uses radius = max(1, ceil(radius - radiusDecayPerStep
    * i)) instead of the fixed `radius` for every step. The idea is to
    search wide early on, when the start tree is furthest from optimal,
    then narrow to cheaper, more local moves later. Under useFastSelection,
    this instead shrinks chooseGraft's own maximum walk distance.

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

    usePrescreen (default false) adds a cheap parsimony-based auto-reject
    check before each candidate's real (expensive) likelihood evaluation,
    in both useFastSelection and the exhaustive branch alike. Once per
    step, before any candidate is generated, the CURRENT tree's parsimony
    score is computed (computeParsimony() -- Fitch parsimony, plain bitwise
    set operations over site patterns, no substitution-model math at all,
    on the order of 1-2 magnitudes cheaper per call than
    computeLikelihood()). Then, for each candidate considered (the single
    one chooseGraft produces under useFastSelection, or each one
    findGraftPositions returns otherwise): apply it, compute its own
    post-move parsimony score the same cheap way, and roll back. If that
    score is worse than the current tree's by more than
    parsimonyPrescreenTolerance() substitutions, the candidate is rejected
    right there -- it never touches computeLikelihood() at all. Only
    candidates that pass this check get the real likelihood evaluation
    (and, in the exhaustive branch, only those survivors are compared
    against each other to find the best one). If every candidate in a
    step gets rejected this way, that step is skipped exactly like a step
    with zero legal candidates. This trades a chance of discarding a
    candidate that would have actually improved the likelihood (parsimony
    and likelihood don't always agree) for skipping the expensive check on
    candidates that are unlikely to be worth it. See
    parsimonyPrescreenTolerance's own comment for how that cutoff was
    chosen -- it is explicitly experimental.

    On completion, prints the Robinson-Foulds distance between the
    original AliSim tree and the final tree to the terminal, and writes
    both trees plus the RF distance to output.txt (repo root, overwritten
    each run).

    @return 0 on success, 1 if the tree/alignment couldn't be read
 */
int runHillClimb(const string &trueTreeArg, int radius, int maxSteps, double radiusDecayPerStep = 0.0,
        bool randomStart = false, bool useFastSelection = false, bool quiet = false, bool usePrescreen = false,
        int numCandidates = 1, bool reoptimizeBranchLengths = false) {
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

    string modelName = "JC";
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.initializeAllPartialLh();
    if (usePrescreen) {
        // computeParsimony() silently scores every site as 0 substitutions
        // (nsites == 0 in computeParsimonyBranchFast) until aln's pattern
        // order/count is set up this way -- IQ-TREE's own pipeline does
        // this once per alignment (see IQTree::doTreeSearch, iqtree.cpp)
        // before ever calling into parsimony; this tool must too
        if (aln->ordered_pattern.empty())
            aln->orderPatternByNumChars(PAT_VARIANT);
        tree.initializeAllPartialPars();
    }

    double curScore = tree.computeLikelihood();

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
    double trueTreeLogl = trueTree.computeLikelihood();

    cout.rdbuf(realCoutBuf);

    cout << (randomStart ? "random start tree: " : "BioNJ start tree : ")
         << newickOf(tree) << " (logL = " << curScore << ")" << endl;
    cout << "radius          : " << radius << endl;
    if (radiusDecayPerStep > 0.0)
        cout << "radius decay    : " << radiusDecayPerStep << " per step (min 1)" << endl;
    cout << "max steps       : " << maxSteps << endl;
    if (useFastSelection)
        cout << "selection       : fast (choosePrune/chooseGraft proposal, not exhaustive)"
             << (numCandidates > 1 ? ", " + to_string(numCandidates) + " candidates/step" : "") << endl;
    if (usePrescreen)
        cout << "prescreen       : parsimony, tolerance " << parsimonyPrescreenTolerance()
             << " (experimental)" << endl;
    if (reoptimizeBranchLengths)
        cout << "branch lengths  : re-optimized (Newton-Raphson, like NNI) on each candidate's 3 "
                "changed edges before scoring, experimental" << endl;

    int step = 0;
    for (; step < maxSteps; step++) {
        // radiusDecayPerStep == 0.0 keeps this equal to `radius` for every
        // step, i.e. identical to the original constant-radius --hillclimb
        int stepRadius = radiusDecayPerStep > 0.0
            ? max(1, (int) ceil(radius - radiusDecayPerStep * step))
            : radius;

        PhyloNode *pruneNode, *pruneDad;
        if (!choosePrune(tree, edgeRegistry, pruneNode, pruneDad)) {
            if (!quiet)
                cout << "step " << (step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }

        // usePrescreen: the parsimony score of the tree as it stands right
        // now, before this step's candidate(s) are even generated -- the
        // baseline each candidate's own post-move parsimony score gets
        // compared against below. Computed at most once per step,
        // regardless of how many candidates there are.
        int curParsimony = 0;
        if (usePrescreen) {
            resetParsimonyBuffers(tree);
            curParsimony = tree.computeParsimony();
        }

        PhyloNode *bestNode, *bestDad;
        int bestDistance;
        double bestScore;

        if (useFastSelection) {
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
            int drawsWithTarget = 0;
            for (int c = 0; c < numCandidates; c++) {
                PhyloNode *candNode, *candDad;
                int walkLength;
                if (!chooseGraft(tree, pruneNode, pruneDad, stepRadius, candNode, candDad, &walkLength))
                    continue; // this draw found no legal target; try the next one
                drawsWithTarget++;

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

                if (usePrescreen) {
                    // auto-reject this candidate on parsimony alone if
                    // it's clearly worse than the tree already is -- skip
                    // the expensive exact likelihood check entirely when so
                    SPRRollback trialRollback;
                    tree.applySPR(move, trialRollback);
                    resetParsimonyBuffers(tree);
                    int candidateParsimony = tree.computeParsimony();
                    tree.rollbackSPR(trialRollback);
                    resetParsimonyBuffers(tree);

                    if (candidateParsimony > curParsimony + parsimonyPrescreenTolerance())
                        continue; // rejected; never touches computeLikelihood()
                }

                double score = scoreTrialSPRMove(tree, move, reoptimizeBranchLengths);

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
                    cout << "step " << (step + 1) << " (radius " << stepRadius << "): prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << (drawsWithTarget == 0
                             ? " -- no legal graft target found in " + to_string(numCandidates) + " draw(s)"
                             : " -- every drawn candidate rejected by parsimony prescreen ("
                                 + to_string(drawsWithTarget) + "/" + to_string(numCandidates) + " had a target)")
                         << "; skipping." << endl;
                continue;
            }
        } else {
            vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, stepRadius);
            if (candidates.empty()) {
                if (!quiet)
                    cout << "step " << (step + 1) << " (radius " << stepRadius << "): prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- no legal graft candidates; skipping." << endl;
                continue;
            }

            // score every candidate on the SAME tree object via apply ->
            // score -> rollback; never allocate a new tree per candidate.
            // usePrescreen: before spending a full computeLikelihood() on
            // a candidate, first auto-reject it on parsimony alone if it's
            // clearly worse than the tree already is (same rule as the
            // fast branch above, just applied to each of many candidates
            // here instead of the sole one there) -- candidates that
            // survive the prescreen are then compared against each other
            // by real likelihood exactly as before.
            bestScore = -DBL_MAX;
            bool haveBest = false;
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

                if (usePrescreen) {
                    SPRRollback trialRollback;
                    tree.applySPR(move, trialRollback);
                    resetParsimonyBuffers(tree);
                    int candidateParsimony = tree.computeParsimony();
                    tree.rollbackSPR(trialRollback);
                    resetParsimonyBuffers(tree);

                    if (candidateParsimony > curParsimony + parsimonyPrescreenTolerance())
                        continue; // rejected; never touches computeLikelihood()
                }

                double score = scoreTrialSPRMove(tree, move, reoptimizeBranchLengths);

                if (!haveBest || score > bestScore) {
                    bestScore = score;
                    bestCandidate = c;
                    haveBest = true;
                }
            }
            if (!haveBest) {
                if (!quiet)
                    cout << "step " << (step + 1) << " (radius " << stepRadius << "): prune {"
                         << describeEdgeCompact(pruneNode, pruneDad) << "}"
                         << " -- every candidate rejected by parsimony prescreen; skipping." << endl;
                continue;
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

        bool improved = bestScore > curScore;
        if (!quiet)
            cout << "step " << (step + 1) << " (radius " << stepRadius << "): prune {"
                 << describeEdgeCompact(pruneNode, pruneDad) << "}"
                 << " -> graft {" << describeEdgeCompact(bestNode, bestDad) << "}"
                 << " (" << (useFastSelection ? "d=" : "distance ") << bestDistance << ")"
                 << ", logL " << bestScore << " (cur " << curScore << ")"
                 << (improved ? " [kept]" : " [reverted]") << endl;

        if (improved) {
            curScore = bestScore;
        } else {
            rollbackSPRTracked(tree, edgeRegistry, bestTracked);
            resetLikelihoodBuffers(tree);
        }
    }

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
    an unconditional-acceptance SPR random walk. This exists purely to
    gather data on how a candidate's parsimony difference relates to its
    likelihood difference and its RF-distance-to-the-true-tree difference,
    without the confound hill-climbing introduces: runHillClimb only ever
    keeps a move that already improved the likelihood, so any candidate
    ever recorded there (with or without usePrescreen) is a biased sample
    for asking "does a candidate's parsimony difference predict its
    likelihood/RF difference in general?" This command instead accepts
    every chosen candidate unconditionally, every step, so every
    parsimony/logL/RF difference recorded here is a plain random draw,
    not one that already survived a real-likelihood filter.

    Unlike runHillClimb, the starting tree is ALWAYS a random Yule-Harding
    topology (never BioNJ) -- this isn't a search-quality test, so starting
    from a "reasonable" tree isn't the point; a maximally uninformed start
    keeps every subsequent step's candidate similarly uninformed, which is
    what an unbiased sample needs. Candidates are picked exactly like
    fast-mode hill-climbing (choosePrune()/chooseGraft(), same O(1)/
    O(distance) engine) but are applied for keeps immediately via
    applySPRTracked() -- there is no trial apply/rollback pair here, since
    the move is going to be kept regardless of what it scores.

    Per step: recompute parsimony (Fitch, computeParsimony()), likelihood
    (computeLikelihood(), plain JC model, no branch-length/parameter
    optimization), and RF distance to trueTreeArg (the AliSim ground-truth
    tree), and record how each changed from the tree as it stood right
    before this step. A step whose chosen prune edge has no legal graft
    target is skipped (still consumes one of maxSteps), same as the same
    situation in runHillClimb; the walk never stops early otherwise, and
    never rolls a move back.

    Writes every step's raw before/after values and diffs to
    randomwalk_data.csv (repo root, overwritten each run) for offline
    analysis; also prints a shorter step,parsimony_diff,logl_diff,rf_diff
    line per step to stdout.

    RESULTS (sim.treefile, ~100 taxa/3000bp, JC model, radius 10, 2000
    steps, one run -- EXPERIMENTAL, single dataset/single seed):
      - parsimony_diff vs logl_diff: Pearson r = -0.56 (moderate, correct-
        sign correlation -- fewer extra substitutions genuinely tends to
        mean better likelihood). Scanning parsimonyPrescreenTolerance()-
        style thresholds T (reject if parsimony_diff > T) against "did
        likelihood actually improve" gives a WIDE, fairly flat plateau of
        ~69% classification accuracy for any T from about -15 to +2
        (peak 69.4% at T=-14), versus a 53% majority-class baseline --
        a real but moderate lift, not a sharp optimum. T=2 (the value
        runHillClimb ships) sits inside that plateau (68.7%), not
        meaningfully worse than the observed peak.
      - parsimony_diff vs rf_diff: RF distance barely moved at all over
        2000 steps (only 14/1869 recorded steps changed it, all by
        exactly +/-2) -- the random start tree is already at RF=194,
        which is 2*100-6, the theoretical MAXIMUM for 100 taxa, and a
        single uninformed SPR move essentially never lands on or off a
        bipartition shared with the true tree by chance. This makes
        "which threshold best aligns with RF" a degenerate question here:
        >99% of steps have rf_diff==0 regardless of parsimony_diff, so
        no threshold does meaningfully better than the trivial "predict
        no change" baseline. The rare 14 exceptions do show a clean
        pattern for what little it's worth: every rf_diff==-2 (improving)
        step had parsimony_diff < -170 and logl_diff > 0, and every
        rf_diff==+2 (worsening) step had parsimony_diff > 420 and
        logl_diff < 0 -- consistent with, but not adding evidence beyond,
        the parsimony/likelihood correlation above, and far too few
        events (n=14) to support picking a specific numeric threshold.
      Bottom line: parsimony difference is a real, moderate proxy for
      likelihood difference in general, with no sharp best cutoff in the
      range tested; it says essentially nothing usable about RF change in
      a random walk, because RF change is itself too rare an event here to
      learn anything from. Retest with more steps/seeds before treating
      the -15..+2 plateau as settled.

    @return 0 on success, 1 if the tree/alignment couldn't be read
 */
int runRandomWalk(const string &trueTreeArg, int radius, int maxSteps) {
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

    // topology-only -- computeRFDist below never touches branch lengths,
    // alignment, or model, so trueTree needs no further setup at all
    PhyloTree trueTree;
    readTreeArg(trueTree, trueTreeArg);

    Params &params = Params::getInstance();
    params.setDefault();

    init_random((int) (time(nullptr) * 1000 + (long) (getRealTime() * 1000) % 1000));

    // silence the library setup noise (Alignment/generateRandomTree/
    // ModelFactory progress messages), same as runHillClimb
    ostringstream suppressedSetupOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedSetupOutput.rdbuf());

    InputType intype;
    Alignment *aln = new Alignment((char*) alnFile.c_str(), (char*) "DNA", intype, "");

    PhyloTree tree(aln);
    tree.setParams(&params);
    // always a random start -- see this function's own comment for why
    tree.generateRandomTree(YULE_HARDING);

    // built once, O(n); kept in sync thereafter by applySPRTracked, since
    // every move here is committed for real, never rolled back
    EdgeRegistry edgeRegistry;
    buildEdgeRegistry(tree, edgeRegistry);

    tree.setNumThreads(1);
    tree.setLikelihoodKernel(LK_SSE2);

    string modelName = "JC";
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.initializeAllPartialLh();
    // see the matching comment in runHillClimb: computeParsimony() silently
    // scores every site as 0 substitutions until this is called once
    if (aln->ordered_pattern.empty())
        aln->orderPatternByNumChars(PAT_VARIANT);
    tree.initializeAllPartialPars();

    double curLogl = tree.computeLikelihood();
    int curParsimony = tree.computeParsimony();

    auto rfToTrueTree = [&](PhyloTree &t) -> int {
        stringstream ss;
        ss << newickOf(t);
        ss.seekg(0, ios::beg);
        vector<double> rfdist;
        trueTree.computeRFDist(ss, rfdist);
        return rfdist.empty() ? -1 : (int) rfdist[0];
    };
    int curRf = rfToTrueTree(tree);

    cout.rdbuf(realCoutBuf);

    cout << "random start tree: " << newickOf(tree) << endl;
    cout << "  (logL = " << curLogl << ", parsimony = " << curParsimony << ", RF = " << curRf << ")" << endl;
    cout << "radius           : " << radius << endl;
    cout << "max steps        : " << maxSteps << endl;
    cout << "every candidate is accepted unconditionally -- this is a random walk, not a search" << endl;
    cout << endl;

    ofstream csv("randomwalk_data.csv");
    if (csv.good())
        csv << "step,parsimony_before,parsimony_after,parsimony_diff,"
               "logl_before,logl_after,logl_diff,rf_before,rf_after,rf_diff" << endl;
    else
        cerr << "warning: could not write to randomwalk_data.csv" << endl;

    cout << "step,parsimony_diff,logl_diff,rf_diff" << endl;

    int step = 0;
    for (; step < maxSteps; step++) {
        PhyloNode *pruneNode, *pruneDad;
        if (!choosePrune(tree, edgeRegistry, pruneNode, pruneDad)) {
            cout << "step " << (step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }

        PhyloNode *graftNode, *graftDad;
        int walkLength;
        if (!chooseGraft(tree, pruneNode, pruneDad, radius, graftNode, graftDad, &walkLength)) {
            cout << "step " << (step + 1) << ": prune {" << describeEdgeCompact(pruneNode, pruneDad)
                 << "} -- no legal graft target; skipping." << endl;
            continue;
        }

        SPRMove move;
        move.prune_node = pruneNode;
        move.prune_dad = pruneDad;
        move.regraft_node = graftNode;
        move.regraft_dad = graftDad;
        move.radius = walkLength;
        move.screening_score = 0.0;
        move.exact_score = 0.0;
        move.candidate_id = 0;
        move.generation = step;

        // committed immediately -- no trial apply/rollback pair, because
        // this candidate is kept no matter what it scores
        TrackedSPR tracked;
        applySPRTracked(tree, edgeRegistry, move, tracked);

        resetParsimonyBuffers(tree);
        int newParsimony = tree.computeParsimony();
        resetLikelihoodBuffers(tree);
        double newLogl = tree.computeLikelihood();
        int newRf = rfToTrueTree(tree);

        int parsimonyDiff = newParsimony - curParsimony;
        double loglDiff = newLogl - curLogl;
        int rfDiff = newRf - curRf;

        cout << (step + 1) << "," << parsimonyDiff << "," << loglDiff << "," << rfDiff << endl;
        if (csv.good())
            csv << (step + 1) << "," << curParsimony << "," << newParsimony << "," << parsimonyDiff << ","
                << curLogl << "," << newLogl << "," << loglDiff << ","
                << curRf << "," << newRf << "," << rfDiff << endl;

        curParsimony = newParsimony;
        curLogl = newLogl;
        curRf = newRf;
    }
    if (csv.good())
        csv.close();

    cout << endl;
    cout << "=== finished after " << step << " step(s) ===" << endl;
    cout << "final tree (logL = " << curLogl << ", parsimony = " << curParsimony << "): " << newickOf(tree) << endl;
    cout << "final RF distance to the original AliSim tree: " << curRf << endl;
    cout << "per-step data written to randomwalk_data.csv" << endl;
    cout << "time elapsed     : " << fixed << setprecision(2) << (getRealTime() - wallClockStart)
         << " sec" << endl;

    delete aln;
    return 0;
}

/**
    a second unconditional-acceptance SPR random walk (see runRandomWalk's
    own comment for why "accept everything" is the right design for this
    kind of measurement), testing a different cheap proxy: instead of
    parsimony, this scores just the LOCAL subtree spanning the prune and
    graft positions, evaluated as its own small standalone tree under the
    same JC model, and compares how ITS likelihood changes to how the
    WHOLE tree's likelihood changes for the same move.

    Per step:
      1. choosePrune()/chooseGraft() pick a candidate exactly like
         runRandomWalk.
      2. Before applying: find L = LCA(pruneDad, graftDad) in the CURRENT
         topology (ancestry rooted at tree.root, rebuilt fresh every step
         since the topology changes every step). Both attachment points
         are pruneDad/graftDad's own descendants-or-self, so everything
         the move can possibly change is confined to the clade hanging
         off L, away from tree.root -- that clade is "the subtree".
      3. Extract that clade's leaf names and build a standalone small
         PhyloTree from a full-precision Newick fragment (see
         writeFullPrecisionNewick) plus a projected sub-alignment
         (Alignment::extractSubAlignment) covering just those leaves,
         with its own fresh JC ModelFactory. computeLikelihood() on this
         small tree gives "the subtree's own likelihood, evaluated as if
         it were an independent tree" -- genuinely cheaper than the whole
         tree's likelihood when the clade is small, since it never
         touches partial likelihoods for anything outside it.
      4. Apply the candidate for real via applySPRTracked (never rolled
         back, exactly like runRandomWalk).
      5. Re-locate the same clade in the NEW topology: since nothing
         outside L's old clade was touched, the same leaf set is still a
         valid clade after the move; its MRCA is refound via a fresh
         ancestry rebuild and findLCA() on two reference leaves (the
         first and last names collected in step 3), which by
         collectLeafNames' traversal order sit on opposite sides of the
         original split. The SAME sub-alignment is reused (same leaves,
         same sequences); only a new small tree is built for the new
         topology.
      6. Recompute the small tree's likelihood (after) and the whole
         tree's likelihood (after, via computeLikelihood() + resetLikelihoodBuffers
         same as runRandomWalk), and record how much each changed from
         the step's own "before" values.
    Wall-clock time spent specifically inside the small-tree construction
    + likelihood calls, and separately inside the whole-tree likelihood
    calls, is accumulated across the whole run and reported at the end,
    to compare the two approaches' actual cost head to head -- both are
    computed every step here (needed to measure their correlation at
    all), so neither number is what a real prescreen/replacement would
    cost in production use, but their RATIO is exactly what answers
    whether the subtree approach is worth using that way.

    A step whose subtree ends up with fewer than 2 leaves (degenerate,
    only possible right at a 3-taxon corner of the tree) skips the
    subtree-specific measurement for that step but still applies the move
    and updates the whole-tree running likelihood, same as any other
    skip case in this file.

    Writes step,subtree_size,subtree_logl_diff,whole_logl_diff to stdout
    and a fuller version (before/after/diff for both, plus per-step
    timings) to randomwalk_subtree_data.csv (repo root, overwritten each
    run).

    RESULTS (sim.treefile, ~100 taxa/3000bp, JC model, radius 10, 1000
    steps/489 measured, one run -- EXPERIMENTAL, single dataset/seed):
      - CORRELATION IS STRONG: Pearson r = 0.87 between subtree_logl_diff
        and whole_logl_diff, with 82% sign agreement (both improve or
        both worsen together) -- this local proxy tracks the whole
        tree's likelihood change far better than the parsimony proxy did
        in runRandomWalk (r=-0.56 there). Unsurprising in hindsight: the
        subtree literally contains every node whose local topology
        changed, so its own likelihood captures most of what moved,
        unlike parsimony which is a different scoring function entirely.
      - EFFICIENCY IS NOT: despite that strong correlation, this
        implementation of "evaluate just the subtree" is SLOWER than
        evaluating the whole tree, at every subtree size tested --
        exactly backwards from what a useful prescreen/replacement needs:
            subtree leaves   avg subtree ms   avg whole-tree ms   ratio
                 0-10             14.9              11.1          1.3x
                10-20             17.3              11.2          1.5x
                20-35             21.2              10.8          2.0x
                35-55             24.5               9.9          2.5x
                55-80             29.8              10.8          2.8x
                80-101            36.8               9.8          3.7x
        The whole-tree number is roughly flat (~10-11ms) regardless of
        subtree size, as expected (it's always the same ~100-taxon
        computation via resetLikelihoodBuffers()+computeLikelihood()).
        The subtree number, despite operating on far fewer leaves, is
        NEVER cheaper -- because this function rebuilds a brand new
        PhyloTree + projected sub-Alignment + fresh JC ModelFactory from
        scratch every single time (twice per step: before and after), and
        that construction overhead (Newick parsing, pattern
        re-compression, buffer allocation) dominates over the actual
        likelihood math it's supposedly saving on, even for a 5-leaf
        subtree.
      CONCLUSION: NOT worth using as either a prescreen or a full
      replacement for the whole-tree likelihood check, as implemented --
      it would strictly add cost (you'd still pay the ~10-11ms whole-tree
      check on anything the prescreen didn't reject, on top of the
      14-37ms subtree check). The correlation result suggests this idea
      COULD pay off with a cheaper implementation -- e.g. computing just
      the inner partial-likelihood vector at the LCA branch directly on
      the EXISTING tree object (IQ-TREE's own computePartialLikelihood,
      combined with the model's stationary frequencies as a stand-in for
      "everything outside the clade"), instead of constructing independent
      PhyloTree/Alignment/ModelFactory objects -- but that requires
      correctly reimplementing low-level likelihood-summation logic this
      tool currently gets for free from computeLikelihood(), which was
      judged too risky to get right without extensive testing; this
      version deliberately traded some possible efficiency for confidence
      the numbers are actually correct.

    @return 0 on success, 1 if the tree/alignment couldn't be read
 */
int runRandomWalkSubtree(const string &trueTreeArg, int radius, int maxSteps) {
    double wallClockStart = getRealTime();

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

    ostringstream suppressedSetupOutput;
    streambuf *realCoutBuf = cout.rdbuf(suppressedSetupOutput.rdbuf());

    InputType intype;
    Alignment *aln = new Alignment((char*) alnFile.c_str(), (char*) "DNA", intype, "");

    PhyloTree tree(aln);
    tree.setParams(&params);
    tree.generateRandomTree(YULE_HARDING);

    EdgeRegistry edgeRegistry;
    buildEdgeRegistry(tree, edgeRegistry);

    tree.setNumThreads(1);
    tree.setLikelihoodKernel(LK_SSE2);

    string modelName = "JC";
    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, modelName, &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.initializeAllPartialLh();

    double curLogl = tree.computeLikelihood();

    cout.rdbuf(realCoutBuf);

    cout << "random start tree: " << newickOf(tree) << " (logL = " << curLogl << ")" << endl;
    cout << "radius           : " << radius << endl;
    cout << "max steps        : " << maxSteps << endl;
    cout << "every candidate is accepted unconditionally -- this is a random walk, not a search" << endl;
    cout << endl;

    ofstream csv("randomwalk_subtree_data.csv");
    if (csv.good())
        csv << "step,subtree_leaves,subtree_logl_before,subtree_logl_after,subtree_logl_diff,"
               "whole_logl_before,whole_logl_after,whole_logl_diff,subtree_eval_sec,whole_eval_sec" << endl;
    else
        cerr << "warning: could not write to randomwalk_subtree_data.csv" << endl;

    cout << "step,subtree_leaves,subtree_logl_diff,whole_logl_diff" << endl;

    double totalSubtreeEvalTime = 0.0, totalWholeTreeEvalTime = 0.0;
    int subtreeStepsMeasured = 0;

    int step = 0;
    for (; step < maxSteps; step++) {
        PhyloNode *pruneNode, *pruneDad;
        if (!choosePrune(tree, edgeRegistry, pruneNode, pruneDad)) {
            cout << "step " << (step + 1) << ": no degree-3 node left to prune from; stopping." << endl;
            step++;
            break;
        }

        PhyloNode *graftNode, *graftDad;
        int walkLength;
        if (!chooseGraft(tree, pruneNode, pruneDad, radius, graftNode, graftDad, &walkLength)) {
            cout << "step " << (step + 1) << ": prune {" << describeEdgeCompact(pruneNode, pruneDad)
                 << "} -- no legal graft target; skipping." << endl;
            continue;
        }

        // find the clade spanning both attachment points, in the CURRENT
        // (pre-move) topology
        ParentMap parentBefore;
        DepthMap depthBefore;
        buildAncestry((PhyloNode*) tree.root, nullptr, 0, parentBefore, depthBefore);
        PhyloNode *lca = findLCA(pruneDad, graftDad, parentBefore, depthBefore);
        PhyloNode *lcaDad = parentBefore.count(lca) ? parentBefore[lca] : nullptr;

        vector<string> subtreeLeaves;
        collectLeafNames(lca, lcaDad, subtreeLeaves);

        bool measureSubtree = subtreeLeaves.size() >= 2;
        Alignment *subAln = nullptr;
        double subtreeLoglBefore = 0.0, subtreeLoglAfter = 0.0;
        double stepSubtreeTime = 0.0;

        if (measureSubtree) {
            double t0 = getRealTime();
            ostringstream suppressedStepOutput;
            streambuf *realCoutBuf2 = cout.rdbuf(suppressedStepOutput.rdbuf());

            IntVector seqIds;
            for (size_t i = 0; i < subtreeLeaves.size(); i++)
                seqIds.push_back(aln->getSeqID(subtreeLeaves[i]));
            subAln = aln->extractSubAlignment(seqIds, 0);

            ostringstream nwkBefore;
            writeFullPrecisionNewick(nwkBefore, lca, lcaDad);
            nwkBefore << ";";

            PhyloTree smallTreeBefore;
            smallTreeBefore.setParams(&params);
            smallTreeBefore.read_TreeString(nwkBefore.str(), false);
            smallTreeBefore.setAlignment(subAln);
            smallTreeBefore.setNumThreads(1);
            smallTreeBefore.setLikelihoodKernel(LK_SSE2);
            ModelsBlock *smallModelsBlockB = readModelsDefinition(params);
            smallTreeBefore.setModelFactory(new ModelFactory(params, modelName, &smallTreeBefore, smallModelsBlockB));
            delete smallModelsBlockB;
            smallTreeBefore.setModel(smallTreeBefore.getModelFactory()->model);
            smallTreeBefore.setRate(smallTreeBefore.getModelFactory()->site_rate);
            smallTreeBefore.initializeAllPartialLh();
            subtreeLoglBefore = smallTreeBefore.computeLikelihood();

            cout.rdbuf(realCoutBuf2);
            stepSubtreeTime += getRealTime() - t0;
        }

        SPRMove move;
        move.prune_node = pruneNode;
        move.prune_dad = pruneDad;
        move.regraft_node = graftNode;
        move.regraft_dad = graftDad;
        move.radius = walkLength;
        move.screening_score = 0.0;
        move.exact_score = 0.0;
        move.candidate_id = 0;
        move.generation = step;

        TrackedSPR tracked;
        applySPRTracked(tree, edgeRegistry, move, tracked);

        bool subtreeMeasured = false;
        if (measureSubtree) {
            double t0 = getRealTime();
            ostringstream suppressedStepOutput;
            streambuf *realCoutBuf2 = cout.rdbuf(suppressedStepOutput.rdbuf());

            // re-locate the same clade in the new topology. In the common
            // case, nothing outside the original clade was touched, so its
            // MRCA is found fresh via two reference leaves from opposite
            // sides of the original split -- but this is NOT guaranteed in
            // every case (e.g. if the pruned branch was itself one whole
            // side of the LCA split, the LCA node itself can get bypassed
            // by the same "degree drops to 2" mechanic that removes
            // pruneDad), so the resulting leaf set is explicitly checked
            // against the original before trusting it, rather than assumed
            ParentMap parentAfter;
            DepthMap depthAfter;
            buildAncestry((PhyloNode*) tree.root, nullptr, 0, parentAfter, depthAfter);
            PhyloNode *leafX = (PhyloNode*) tree.findLeafName(subtreeLeaves.front());
            PhyloNode *leafY = (PhyloNode*) tree.findLeafName(subtreeLeaves.back());
            PhyloNode *lcaAfter = findLCA(leafX, leafY, parentAfter, depthAfter);
            PhyloNode *lcaAfterDad = parentAfter.count(lcaAfter) ? parentAfter[lcaAfter] : nullptr;

            vector<string> afterLeaves;
            collectLeafNames(lcaAfter, lcaAfterDad, afterLeaves);
            vector<string> beforeSorted = subtreeLeaves, afterSorted = afterLeaves;
            sort(beforeSorted.begin(), beforeSorted.end());
            sort(afterSorted.begin(), afterSorted.end());

            if (afterSorted == beforeSorted) {
                ostringstream nwkAfter;
                writeFullPrecisionNewick(nwkAfter, lcaAfter, lcaAfterDad);
                nwkAfter << ";";

                PhyloTree smallTreeAfter;
                smallTreeAfter.setParams(&params);
                smallTreeAfter.read_TreeString(nwkAfter.str(), false);
                smallTreeAfter.setAlignment(subAln);
                smallTreeAfter.setNumThreads(1);
                smallTreeAfter.setLikelihoodKernel(LK_SSE2);
                ModelsBlock *smallModelsBlockA = readModelsDefinition(params);
                smallTreeAfter.setModelFactory(new ModelFactory(params, modelName, &smallTreeAfter, smallModelsBlockA));
                delete smallModelsBlockA;
                smallTreeAfter.setModel(smallTreeAfter.getModelFactory()->model);
                smallTreeAfter.setRate(smallTreeAfter.getModelFactory()->site_rate);
                smallTreeAfter.initializeAllPartialLh();
                subtreeLoglAfter = smallTreeAfter.computeLikelihood();
                subtreeMeasured = true;
            }

            cout.rdbuf(realCoutBuf2);
            stepSubtreeTime += getRealTime() - t0;

            delete subAln;
            if (subtreeMeasured) {
                totalSubtreeEvalTime += stepSubtreeTime;
                subtreeStepsMeasured++;
            }
        }

        double t1 = getRealTime();
        resetLikelihoodBuffers(tree);
        double newLogl = tree.computeLikelihood();
        double stepWholeTime = getRealTime() - t1;
        totalWholeTreeEvalTime += stepWholeTime;

        double wholeDiff = newLogl - curLogl;

        if (subtreeMeasured) {
            double subtreeDiff = subtreeLoglAfter - subtreeLoglBefore;
            cout << (step + 1) << "," << subtreeLeaves.size() << "," << subtreeDiff << "," << wholeDiff << endl;
            if (csv.good())
                csv << (step + 1) << "," << subtreeLeaves.size() << "," << subtreeLoglBefore << ","
                    << subtreeLoglAfter << "," << subtreeDiff << ","
                    << curLogl << "," << newLogl << "," << wholeDiff << ","
                    << stepSubtreeTime << "," << stepWholeTime << endl;
        }

        curLogl = newLogl;
    }
    if (csv.good())
        csv.close();

    cout << endl;
    cout << "=== finished after " << step << " step(s) ===" << endl;
    cout << "final tree (logL = " << curLogl << "): " << newickOf(tree) << endl;
    cout << "per-step data written to randomwalk_subtree_data.csv" << endl;
    cout << "subtree steps measured : " << subtreeStepsMeasured << " / " << step << endl;
    if (subtreeStepsMeasured > 0) {
        cout << "total time in subtree eval  : " << fixed << setprecision(4) << totalSubtreeEvalTime
             << " sec (avg " << (totalSubtreeEvalTime / subtreeStepsMeasured * 1000.0) << " ms/step)" << endl;
        cout << "total time in whole-tree eval: " << fixed << setprecision(4) << totalWholeTreeEvalTime
             << " sec (avg " << (totalWholeTreeEvalTime / subtreeStepsMeasured * 1000.0) << " ms/step)" << endl;
    }
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
    cerr << "      <radius> hops of the prune point, without applying any of them." << endl;
    cerr << endl;
    cerr << "  " << prog << " --likelihood <tree.nwk | \"(newick,string);\"> <alignment.fasta>" << endl;
    cerr << "      evaluate the log-likelihood of the given tree (topology and branch" << endl;
    cerr << "      lengths as given, no optimization) against a DNA alignment under a" << endl;
    cerr << "      plain JC model. Sequence names in the alignment must match the tree's" << endl;
    cerr << "      leaf names exactly." << endl;
    cerr << endl;
    cerr << "  " << prog << " --hillclimb <alisim-tree.treefile> <radius> <max-steps> [random] [fast [N]] [quiet] [prescreen]" << endl;
    cerr << "      greedy randomized SPR search: build a BioNJ start tree from the" << endl;
    cerr << "      alignment AliSim simulated from <alisim-tree.treefile> (found by" << endl;
    cerr << "      replacing '.treefile' with '.fa'), then repeatedly prune a random edge," << endl;
    cerr << "      evaluate every legal regraft within <radius> hops via applySPR/" << endl;
    cerr << "      rollbackSPR on one tree object, and keep the best if it improves the" << endl;
    cerr << "      likelihood, for up to <max-steps> rounds. Prints the RF distance to the" << endl;
    cerr << "      original AliSim tree and writes both trees + the RF distance to" << endl;
    cerr << "      output.txt. Five optional trailing flags, in any order:" << endl;
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
    cerr << "        prescreen  before each candidate's real likelihood check, auto-reject it" << endl;
    cerr << "                   on a cheap parsimony score alone if that's clearly worse than" << endl;
    cerr << "                   the current tree already is -- see parsimonyPrescreenTolerance()" << endl;
    cerr << "                   in the source. EXPERIMENTAL: tuned informally against one" << endl;
    cerr << "                   dataset, and on that dataset did not show a net speed benefit" << endl;
    cerr << "                   under a plain JC model (see its comment for the actual numbers)" << endl;
    cerr << "                   -- likelier to help with a heavier substitution model" << endl;
    cerr << "        reopt      re-optimize (Newton-Raphson) the 3 edges an SPR move actually" << endl;
    cerr << "                   changes before scoring a candidate, the same way IQ-TREE's own" << endl;
    cerr << "                   NNI search re-optimizes the branches it touches -- instead of" << endl;
    cerr << "                   trusting applySPR's naive placeholder lengths (half the target" << endl;
    cerr << "                   edge split evenly, the two vacated edges summed). Can only ever" << endl;
    cerr << "                   improve or leave unchanged a candidate's reported likelihood for" << endl;
    cerr << "                   its topology. EXPERIMENTAL and noticeably slower per candidate" << endl;
    cerr << "                   (real branch-length search plus the safe/scaled likelihood" << endl;
    cerr << "                   kernel this needs -- see scoreTrialSPRMove's comment in the" << endl;
    cerr << "                   source for why); ~1.3x slower in one informal comparison, for a" << endl;
    cerr << "                   modest logL improvement (-7.60e+05 -> -7.59e+05 at radius 10," << endl;
    cerr << "                   15 steps on sim.treefile)" << endl;
    cerr << endl;
    cerr << "  " << prog << " --hillclimb-decay <alisim-tree.treefile> <radius> <max-steps> <decay> [random] [fast [N]] [quiet] [prescreen] [reopt]" << endl;
    cerr << "      same search as --hillclimb, but the radius shrinks by <decay> each step" << endl;
    cerr << "      (step i uses radius = max(1, ceil(<radius> - <decay> * i)) instead of a" << endl;
    cerr << "      fixed <radius> for every step). Accepts the same trailing 'random', 'fast [N]'," << endl;
    cerr << "      'quiet', 'prescreen', and 'reopt' flags as --hillclimb, in any order." << endl;
    cerr << endl;
    cerr << "  " << prog << " --randomwalk <alisim-tree.treefile> <radius> <max-steps>" << endl;
    cerr << "      NOT a search: an unconditional-acceptance SPR random walk, always" << endl;
    cerr << "      starting from a random Yule-Harding tree, that accepts every" << endl;
    cerr << "      choosePrune()/chooseGraft()-picked candidate regardless of whether it" << endl;
    cerr << "      improves. Exists to gather an unbiased sample of how a candidate's" << endl;
    cerr << "      parsimony difference relates to its likelihood difference and its RF-" << endl;
    cerr << "      distance-to-<alisim-tree.treefile> difference (hill-climbing only ever" << endl;
    cerr << "      keeps improving moves, which biases that question). Prints and writes" << endl;
    cerr << "      to randomwalk_data.csv one row per step: parsimony/logL/RF before, after," << endl;
    cerr << "      and their difference." << endl;
    cerr << endl;
    cerr << "  " << prog << " --randomwalk-subtree <alisim-tree.treefile> <radius> <max-steps>" << endl;
    cerr << "      Same idea as --randomwalk, but tests a different cheap proxy: finds the" << endl;
    cerr << "      LCA of the prune and graft positions, evaluates JUST that clade as its" << endl;
    cerr << "      own small standalone tree (own Newick + projected sub-alignment, same" << endl;
    cerr << "      JC model) before and after the move, and compares its likelihood change" << endl;
    cerr << "      to the whole tree's likelihood change. Prints/writes" << endl;
    cerr << "      randomwalk_subtree_data.csv with both series plus per-step timings for" << endl;
    cerr << "      each. EXPERIMENTAL result on one dataset (see runRandomWalkSubtree's" << endl;
    cerr << "      comment for the numbers): the subtree score correlates STRONGLY with the" << endl;
    cerr << "      whole tree's likelihood change (r=0.87, much stronger than the parsimony" << endl;
    cerr << "      proxy's r=-0.56), but this implementation is SLOWER than just checking the" << endl;
    cerr << "      whole tree at every subtree size tested, because rebuilding a standalone" << endl;
    cerr << "      PhyloTree/Alignment/ModelFactory every step costs more than the likelihood" << endl;
    cerr << "      math it saves -- not worth using as a prescreen or replacement as built." << endl;
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
    cerr << "    " << prog << " --hillclimb sim.treefile 15 20 prescreen    (parsimony auto-reject, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb sim.treefile 10 20 reopt        (re-optimize branch lengths, experimental)" << endl;
    cerr << "    " << prog << " --hillclimb-decay sim.treefile 6 20 0.5  (hill-climb, shrinking radius)" << endl;
    cerr << "    " << prog << " --hillclimb-decay sim.treefile 6 20 0.5 random  (same, random start tree)" << endl;
    cerr << "    " << prog << " --hillclimb-decay sim.treefile 6 20 0.5 fast    (same, proposal search)" << endl;
    cerr << "    " << prog << " --randomwalk sim.treefile 10 300   (unbiased parsimony/logL/RF data)" << endl;
    cerr << "    " << prog << " --randomwalk-subtree sim.treefile 10 300  (LCA-subtree vs whole-tree logL)" << endl;
    cerr << endl;
    cerr << "  Full reference: tree/spr_topology_test_usage.txt" << endl;
}

/**
    parse the trailing optional flags shared by --hillclimb and
    --hillclimb-decay: the literal words "random", "fast", "quiet", and
    "prescreen", in any order, each at most once. argv[fromIndex..argc-1]
    must consist of exactly these (in any combination); anything else
    (typos, duplicates, unrelated tokens) is treated as a parse failure so
    main() falls through to printUsage() rather than silently ignoring a
    misspelled flag. "prescreen" combines with "fast" too -- in that case
    the auto-reject check applies to each of fast mode's candidates
    individually, exactly as it does for each candidate in the exhaustive
    (non-fast) branch; see usePrescreen's comment on runHillClimb.

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
    lengths. See scoreTrialSPRMove's comment on runHillClimb.
    @return false if any trailing argument isn't recognized
 */
bool parseHillClimbFlags(int argc, char **argv, int fromIndex, bool &randomStart, bool &useFastSelection,
        bool &quiet, bool &usePrescreen, int &numCandidates, bool &reoptimizeBranchLengths) {
    randomStart = false;
    useFastSelection = false;
    quiet = false;
    usePrescreen = false;
    numCandidates = 1;
    reoptimizeBranchLengths = false;
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
        else if (arg == "prescreen" && !usePrescreen)
            usePrescreen = true;
        else if (arg == "reopt" && !reoptimizeBranchLengths)
            reoptimizeBranchLengths = true;
        else
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
    if (argc == 5 && string(argv[1]) == "--randomwalk")
        return runRandomWalk(argv[2], atoi(argv[3]), atoi(argv[4]));
    if (argc == 5 && string(argv[1]) == "--randomwalk-subtree")
        return runRandomWalkSubtree(argv[2], atoi(argv[3]), atoi(argv[4]));
    if (argc >= 5 && string(argv[1]) == "--hillclimb") {
        bool randomStart, useFastSelection, quiet, usePrescreen, reoptimizeBranchLengths;
        int numCandidates;
        if (parseHillClimbFlags(argc, argv, 5, randomStart, useFastSelection, quiet, usePrescreen, numCandidates,
                reoptimizeBranchLengths))
            return runHillClimb(argv[2], atoi(argv[3]), atoi(argv[4]), 0.0, randomStart, useFastSelection, quiet,
                    usePrescreen, numCandidates, reoptimizeBranchLengths);
    }
    if (argc >= 6 && string(argv[1]) == "--hillclimb-decay") {
        bool randomStart, useFastSelection, quiet, usePrescreen, reoptimizeBranchLengths;
        int numCandidates;
        if (parseHillClimbFlags(argc, argv, 6, randomStart, useFastSelection, quiet, usePrescreen, numCandidates,
                reoptimizeBranchLengths))
            return runHillClimb(argv[2], atoi(argv[3]), atoi(argv[4]), atof(argv[5]), randomStart, useFastSelection,
                    quiet, usePrescreen, numCandidates, reoptimizeBranchLengths);
    }
    if (argc == 4 && string(argv[1]) == "--likelihood")
        return runLikelihood(argv[2], argv[3]);
    if (argc == 4)
        return runManualSPR(argv[1], argv[2], argv[3]);
    printUsage(argv[0]);
    return 2;
}
