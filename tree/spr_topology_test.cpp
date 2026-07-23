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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
    stringstream ss;
    tree.printTree(ss, WT_BR_LEN | WT_SORT_TAXA);
    return ss.str();
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
    load a tree, prune pruneLeafName away from its current position, and
    regraft it directly above regraftLeafName (i.e. so the two become
    siblings), printing the before/after Newick and the legality verdict.
    Does not roll back -- this is a one-shot "what would this move do?"
    check, not a round-trip test.
    @return 0 if the move was legal and applied, 1 otherwise
 */
int runManualSPR(const string &treeArg, const string &pruneLeafName, const string &regraftLeafName) {
    PhyloTree tree;
    readTreeArg(tree, treeArg);

    cout << "input tree : " << newickOf(tree) << endl;

    PhyloNode *pruneLeaf = (PhyloNode*) tree.findLeafName(pruneLeafName);
    if (!pruneLeaf) {
        cerr << "error: leaf '" << pruneLeafName << "' not found in tree" << endl;
        return 1;
    }
    if (!pruneLeaf->isLeaf()) {
        cerr << "error: '" << pruneLeafName << "' is not a leaf name in this tree" << endl;
        return 1;
    }

    PhyloNode *regraftLeaf = (PhyloNode*) tree.findLeafName(regraftLeafName);
    if (!regraftLeaf) {
        cerr << "error: leaf '" << regraftLeafName << "' not found in tree" << endl;
        return 1;
    }
    if (!regraftLeaf->isLeaf()) {
        cerr << "error: '" << regraftLeafName << "' is not a leaf name in this tree" << endl;
        return 1;
    }

    SPRMove move;
    move.prune_node = pruneLeaf;
    move.prune_dad = (PhyloNode*) pruneLeaf->neighbors[0]->node;
    move.regraft_node = regraftLeaf;
    move.regraft_dad = (PhyloNode*) regraftLeaf->neighbors[0]->node;
    move.radius = 0;
    move.screening_score = 0.0;
    move.exact_score = 0.0;
    move.candidate_id = 0;
    move.generation = 0;

    cout << "move       : prune '" << pruneLeafName << "', regraft onto the pendant edge above '"
         << regraftLeafName << "' (so they become siblings)" << endl;

    if (!tree.isLegalSPR(move)) {
        cout << "RESULT     : ILLEGAL MOVE" << endl;
        cout << "  (common causes: the two leaves are already siblings/adjacent; the regraft leaf" << endl;
        cout << "   lies inside the subtree being pruned; or the pruned leaf's parent is not a" << endl;
        cout << "   normal degree-3 internal node, e.g. it's directly on the arbitrary root)" << endl;
        return 1;
    }

    SPRRollback rollback;
    tree.applySPR(move, rollback);

    cout << "result tree: " << newickOf(tree) << endl;
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
    cerr << "  " << prog << " <tree.nwk | \"(newick,string);\"> <prune-leaf> <regraft-leaf>" << endl;
    cerr << "      prune <prune-leaf> and regraft it directly above <regraft-leaf>," << endl;
    cerr << "      printing the tree before and after the move" << endl;
}

int main(int argc, char **argv) {
    if (argc == 1)
        return runSelfTest();
    if (argc == 4)
        return runManualSPR(argv[1], argv[2], argv[3]);
    printUsage(argv[0]);
    return 1;
}
