/***************************************************************************
 *   Standalone topology-only test for the modern SPR API                *
 *   (isLegalSPR / applySPR / rollbackSPR, declared in phylotree.h).     *
 *                                                                        *
 *   This is a plain executable with its own main(); it never touches    *
 *   an alignment, a model, or the likelihood machinery, so it exercises *
 *   applySPR()/rollbackSPR() without running a real IQ-TREE analysis.   *
 *   Build target: spr_topology_test (see root CMakeLists.txt).          *
 *                                                                        *
 *   The hardcoded start tree mirrors                                    *
 *   test_scripts/test_data/spr/six_taxa.start.tree so the fixture file  *
 *   and this test describe the same topology; the test does not read    *
 *   the file itself, to stay runnable from any working directory.       *
 ***************************************************************************/

#include "phylotree.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
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
        return 1;
    }
    if (!resolveEdgeSpec(tree, regraftSpec, parent, nodeDepth, regraftNode, regraftDad, err)) {
        cerr << "error resolving regraft edge '" << regraftSpec << "': " << err << endl;
        return 1;
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
        return 1;
    }

    SPRRollback rollback;
    tree.applySPR(move, rollback);

    cout << "result tree : " << newickOf(tree) << endl;
    return 0;
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
        return 1;
    }

    cout << "prune edge : above {" << pruneSpec << "}"
         << (pruneNode->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;
    cout << "radius     : " << radius << endl;
    cout << endl;

    vector<GraftCandidate> candidates = findGraftPositions(tree, pruneNode, pruneDad, radius);

    if (candidates.empty()) {
        cout << "no legal graft positions found within radius " << radius << endl;
        return 0;
    }

    cout << "found " << candidates.size() << " legal graft position(s):" << endl;
    for (size_t i = 0; i < candidates.size(); i++) {
        const GraftCandidate &c = candidates[i];
        cout << "  [" << (i + 1) << "] radius " << c.radius << ": above {" << describeEdge(c.node, c.dad) << "}"
             << (c.node->isLeaf() ? " (pendant edge)" : " (internal edge)") << endl;
    }
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
        return 1;

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
    cerr << "  In both forms, an edge is a comma-separated leaf name list:" << endl;
    cerr << "    a single leaf, e.g. C         -> that leaf's own pendant edge" << endl;
    cerr << "    two or more leaves, e.g. B,D  -> the internal edge above their MRCA" << endl;
    cerr << endl;
    cerr << "  Examples:" << endl;
    cerr << "    " << prog << " tree.nwk C D                      (leaf onto leaf)" << endl;
    cerr << "    " << prog << " tree.nwk C \"B,D\"                  (leaf onto an internal edge)" << endl;
    cerr << "    " << prog << " tree.nwk \"B,D\" \"E,F\"              (internal edge onto internal edge)" << endl;
    cerr << "    " << prog << " --list-grafts tree.nwk C 3         (list candidates within 3 hops of C)" << endl;
}

int main(int argc, char **argv) {
    if (argc == 1)
        return runSelfTest();
    if (argc == 4)
        return runManualSPR(argv[1], argv[2], argv[3]);
    if (argc == 5 && string(argv[1]) == "--list-grafts")
        return runListGrafts(argv[2], argv[3], atoi(argv[4]));
    printUsage(argv[0]);
    return 1;
}
