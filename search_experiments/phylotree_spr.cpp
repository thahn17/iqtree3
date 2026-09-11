/***************************************************************************
 *   search_experiments/phylotree_spr.cpp                                 *
 *                                                                         *
 *   PhyloTree's transactional SPR move API -- isLegalSPR, applySPR,       *
 *   rollbackSPR -- declared in tree/phylotree.h. Kept here, with the rest *
 *   of the search experiments, rather than in tree/phylotree.cpp.         *
 ***************************************************************************/

#include "tree/phylotree.h"

/****************************************************************************
 Subtree Pruning and Regrafting by maximum likelihood (modern API)
 Topology-only "basic" implementation of the transactional move API declared
 in phylotree.h (see SPR_IMPLEMENTATION_PLAN section 4.3). isLegalSPR,
 applySPR and rollbackSPR rewire/restore tree topology and branch lengths
 only; they never touch cached partial likelihoods. Until
 invalidateSPRPartials()/scoreSPR()/commitSPR() are implemented, callers
 must clear partial likelihoods themselves (e.g. clearAllPartialLH())
 before computing likelihood on a tree mutated by applySPR().
 ****************************************************************************/

/**
    true if `node` lies in the subtree rooted at `subtreeRoot`, viewing the
    tree with the branch to `awayFrom` cut off
 */
static bool isInPrunedSubtree(Node *node, Node *subtreeRoot, Node *awayFrom) {
    if (node == subtreeRoot)
        return true;
    for (NeighborVec::iterator it = subtreeRoot->neighbors.begin(); it != subtreeRoot->neighbors.end(); it++)
        if ((*it)->node != awayFrom && isInPrunedSubtree(node, (*it)->node, subtreeRoot))
            return true;
    return false;
}

/**
    record the current (node, length) of `nei` into `rollback`, before
    applySPR() overwrites it
 */
static void saveSPRNeighbor(SPRRollback &rollback, Neighbor *nei) {
    SPRRollback::SavedNeighbor saved;
    saved.nei = nei;
    saved.orig_node = nei->node;
    saved.orig_length = nei->length;
    rollback.saved_neighbors.push_back(saved);
}

bool PhyloTree::isLegalSPR(const SPRMove &move) {
    PhyloNode *node1 = move.prune_node;
    PhyloNode *dad1 = move.prune_dad;
    PhyloNode *node2 = move.regraft_node;
    PhyloNode *dad2 = move.regraft_dad;

    if (!node1 || !dad1 || !node2 || !dad2)
        return false;
    if (!dad1->isNeighbor(node1) || !dad2->isNeighbor(node2))
        return false;
    // dad1 must become degree-2 (and thus suppressible) once node1 is pruned
    if (dad1->degree() != 3)
        return false;
    // the target edge cannot be incident to the node being suppressed, or be
    // the pruned edge itself
    if (node2 == dad1 || dad2 == dad1 || node2 == node1 || dad2 == node1)
        return false;
    // the target edge must lie entirely outside the pruned subtree
    if (isInPrunedSubtree(node2, node1, dad1) || isInPrunedSubtree(dad2, node1, dad1))
        return false;
    // the target edge cannot be incident to the tree's root. `root` is an
    // arbitrary leaf used only as a traversal anchor (see MTree::readTree),
    // not a real attachment point; regrafting onto its pendant edge would
    // rewire that bookkeeping leaf's own neighbor for no biological reason.
    if (node2 == root || dad2 == root)
        return false;
    return true;
}

void PhyloTree::applySPR(const SPRMove &move, SPRRollback &rollback) {
    ASSERT(isLegalSPR(move));

    PhyloNode *node1 = move.prune_node;
    PhyloNode *dad1 = move.prune_dad;
    PhyloNode *node2 = move.regraft_node;
    PhyloNode *dad2 = move.regraft_dad;

    rollback.saved_neighbors.clear();
    rollback.clear_all_partial_lh = true;

    // --- detach the pruned subtree by suppressing dad1 ---------------------
    // dad1 has exactly 3 neighbors: node1 and two others ("siblings"), which
    // become directly connected to each other once dad1 is bypassed. dad1
    // itself is kept around (not deleted/reallocated) so that a pointer to
    // it -- including the tree's root, if dad1 happens to be it -- stays
    // valid; its two now-spare neighbor slots are repurposed below.
    PhyloNode *sibling1 = nullptr, *sibling2 = nullptr;
    FOR_NEIGHBOR_DECLARE(dad1, node1, it)
    {
        if (!sibling1)
            sibling1 = (PhyloNode*) (*it)->node;
        else
            sibling2 = (PhyloNode*) (*it)->node;
    }
    ASSERT(sibling1 && sibling2);

    Neighbor *sibling1_recip = sibling1->findNeighbor(dad1);
    Neighbor *sibling2_recip = sibling2->findNeighbor(dad1);
    double sum_len = sibling1_recip->length + sibling2_recip->length;

    saveSPRNeighbor(rollback, sibling1_recip);
    saveSPRNeighbor(rollback, sibling2_recip);

    sibling1->updateNeighbor(dad1, sibling2, sum_len);
    sibling2->updateNeighbor(dad1, sibling1, sum_len);

    // --- split the target edge and graft dad1 (now spare) into it ----------
    Neighbor *node2_recip = node2->findNeighbor(dad2);
    Neighbor *dad2_recip = dad2->findNeighbor(node2);
    double half_len = node2_recip->length / 2.0;

    saveSPRNeighbor(rollback, node2_recip);
    saveSPRNeighbor(rollback, dad2_recip);

    // repurpose dad1's two now-stale neighbor slots (that used to point at
    // sibling1/sibling2) to point at dad2 and node2 respectively
    bool first = true;
    FOR_NEIGHBOR_IT(dad1, node1, it2)
    {
        saveSPRNeighbor(rollback, *it2);
        if (first) {
            (*it2)->node = dad2;
            (*it2)->length = half_len;
            first = false;
        } else {
            (*it2)->node = node2;
            (*it2)->length = half_len;
        }
    }
    dad2_recip->node = dad1;
    dad2_recip->length = half_len;
    node2_recip->node = dad1;
    node2_recip->length = half_len;
}

void PhyloTree::rollbackSPR(const SPRRollback &rollback) {
    for (vector<SPRRollback::SavedNeighbor>::const_reverse_iterator it = rollback.saved_neighbors.rbegin();
            it != rollback.saved_neighbors.rend(); it++) {
        it->nei->node = it->orig_node;
        it->nei->length = it->orig_length;
    }
}
