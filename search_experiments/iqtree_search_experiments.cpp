/***************************************************************************
 *   search_experiments/iqtree_search_experiments.cpp                     *
 *                                                                         *
 *   Every IQTree member function added by the search experiments          *
 *   (--spr-refine, --spr-perturb, --spr-continuous, --perturb-slack,      *
 *   --accept-dist, "record", "sweep", ...). Declared in tree/iqtree.h     *
 *   between the SEARCH EXPERIMENTS markers; tree/iqtree.cpp only calls    *
 *   the hooks at the bottom of this file, each call site marked           *
 *   "[search-experiments hook]". See search_experiments/README.md.        *
 ***************************************************************************/

#include "tree/iqtree.h"
#include "tree/phylosupertree.h"
#include "utils/timeutil.h"
#include "utils/MPIHelper.h"
#include "utils/pllnni.h"
#include "utils/tools.h"
#include <cfloat>
#include <cmath>

/****************************************************************************
 SPR refinement -- the --spr-refine alternative to doNNISearch
 ****************************************************************************/

void IQTree::initRefinement(Params &params) {
    refine_kick_count = 0;

    // --accept-dist: "temp T", "shape S", "anneal", "floor F", in any
    // order. A bare --accept-dist (empty spec) means every default, which
    // is the point of the defaults -- see sprsearch.h's AcceptDist for why
    // T = 0.5 is the convergence-oriented choice.
    if (params.accept_dist) {
        accept_dist.enabled = true;
        istringstream tokens(params.accept_dist_spec);
        string t;
        while (tokens >> t) {
            if (t == "temp") {
                if (!(tokens >> accept_dist.temperature))
                    outError("--accept-dist: \"temp\" needs a value");
            } else if (t == "shape") {
                if (!(tokens >> accept_dist.shape))
                    outError("--accept-dist: \"shape\" needs a value");
            } else if (t == "floor") {
                if (!(tokens >> accept_dist.tempFloor))
                    outError("--accept-dist: \"floor\" needs a value");
            } else if (t == "anneal") {
                accept_dist.anneal = true;
            } else {
                outError("--accept-dist: unknown token \"" + t + "\" "
                         "(expected temp/shape/anneal/floor)");
            }
        }
        if (accept_dist.temperature <= 0.0)
            outError("--accept-dist: temperature must be positive");
        if (accept_dist.shape <= 0.0)
            outError("--accept-dist: shape must be positive");
        if (accept_dist.tempFloor < 0.0)
            outError("--accept-dist: floor must not be negative");
        cout << "Stochastic acceptance: P(keep) = exp(-(|dL|/T)^" << accept_dist.shape
             << "), T = " << accept_dist.temperature
             << (accept_dist.anneal ? " annealing to " : " constant (floor ")
             << accept_dist.tempFloor << (accept_dist.anneal ? "" : ")")
             << "; perturbation stage suppressed" << endl;
        // the SPR refiner reads the rule off its own options struct
        spr_opt.acceptDist = accept_dist;
    }

    // --spr-perturb is independent of the refinement mode: it can be given
    // on its own (SPR kick, NNI refinement), so it is parsed and validated
    // before the early return below.
    if (params.spr_perturb) {
        string perr;
        if (!sprsearch::parsePerturbSpec(params.spr_perturb_spec, spr_perturb_opt, perr))
            outError("--spr-perturb: " + perr);
        if (rooted || params.is_rooted || params.root != nullptr)
            outError("--spr-perturb does not support rooted trees; SPR moves are written"
                    " against an unrooted topology");
        if (isSuperTree())
            outError("--spr-perturb does not support partition models yet");
        cout << "Perturbation     : random SPR moves"
             << (params.spr_perturb_spec.empty() ? "" : " [" + params.spr_perturb_spec + "]") << endl;
    }

    // --perturb-slack is the kick-agnostic spelling of "slack", so with an
    // SPR kick it has to reach the SPR kick's own settings -- doRandomSPRs
    // reads only those. (doRandomNNIs, which reads the Params copy, never
    // runs under --spr-perturb.) Without this, the combination was silently
    // a no-op.
    if (params.perturb_slack && params.spr_perturb) {
        if (spr_perturb_opt.slackFlag)
            outError("give the kick's slack once: either --perturb-slack or"
                     " --spr-perturb \"... slack D\", not both");
        spr_perturb_opt.slackFlag = true;
        spr_perturb_opt.slackDelta = params.perturb_slack_delta;
        spr_perturb_opt.slackAnneal = params.perturb_slack_anneal;
    }

    // no flag given: stock IQ-TREE, nothing to set up and nothing to report
    if (params.refine_mode == REFINE_NNI && params.refine_spec.empty())
        return;

    bool sprMode = (params.refine_mode == REFINE_SPR);
    string err;
    if (!sprsearch::parseRefineSpec(params.refine_spec, sprMode, spr_opt, err))
        outError(string(sprMode ? "--spr-refine: " : "--nni-refine: ") + err);

    if (sprMode) {
        // The shared SPR machinery works on one plain unrooted PhyloTree
        // whose root is an arbitrary leaf (isLegalSPR/applySPR and
        // buildEdgeRegistry are both written against exactly that), and
        // scores candidates through PhyloTree::computeLikelihood. None of
        // the shapes below is that, and each would fail in its own quiet
        // way rather than loudly, so refuse them here -- before the run
        // spends anything on model fitting or starting trees. The tree
        // itself does not exist yet at this point, so rootedness is read
        // off the parameters that force it rather than off this->rooted;
        // doSPRSearch re-checks the built tree for the same thing.
        if (params.is_rooted || params.root != nullptr)
            outError("--spr-refine does not support rooted trees; the SPR search is written"
                    " against an unrooted topology (drop -o/--root, or refine with NNI)");
        if (isSuperTree())
            outError("--spr-refine does not support partition models yet; its moves are applied to a"
                    " single tree, with no per-partition tree mapping (use NNI for -p/-q analyses)");
        if (params.pll)
            outError("--spr-refine has no PLL implementation; run with -nopll (IQ-TREE's own kernel)");
        if (!sprsearch::isSupportedSeqType(aln->seq_type))
            outError("--spr-refine supports DNA and protein alignments only");
    }

    // A real Newton-Raphson branch-length search can legitimately push a
    // length toward an extreme while searching, which the plain (non-scaled)
    // kernel is not built to handle without numerical underflow in the
    // likelihood derivative -- exactly what IQ-TREE's own "-safe" exists
    // for. spr_topology_test turns this on for the same two flags (see
    // runHillClimb); the naive fixed-length scoring path never searches
    // branch lengths at all, so it does not need it. This has to happen
    // before initSettings' setLikelihoodKernel call, which is where
    // safe_numeric is latched from this flag -- see the call site.
    if (spr_opt.reoptimizeBranchLengths || spr_opt.fullReoptEveryNSteps > 0) {
        params.lk_safe_scaling = true;
        cout << "NOTE: enabling safe numerical scaling (-safe) for --spr-refine's"
                " branch-length re-optimization" << endl;
    }

    // The CPU clock the record CSV's time column is measured against. Set
    // here, at the earliest point the run is known to be happening, so
    // that column covers the whole analysis (candidate-set construction
    // included) rather than just the search loop.
    spr_state.cpuClockStart = getCPUTime();
    // There is no ground-truth tree in a real analysis, so the record
    // CSV's "gap to true tree" column is NaN throughout -- the same thing
    // --hillclimb's own "notree" mode writes.
    spr_state.trueTreeLogl = std::numeric_limits<double>::quiet_NaN();
    spr_state.trueTreeLoglForFindopt = std::numeric_limits<double>::quiet_NaN();
    spr_state.shrinkCurrentRadius = spr_opt.radius;

    cout << "Refinement stage : " << (sprMode ? "SPR hill-climbing" : "NNI (default)");
    if (!params.refine_spec.empty())
        cout << " [" << params.refine_spec << "]";
    cout << endl;
}

void IQTree::finalizeRefineIdentity() {
    // Deferred out of initRefinement because the record CSV is named after
    // the model, and getModelName() reads model/model_factory -- neither of
    // which exists yet when initSettings runs (the model is built and fitted
    // afterwards, and under ModelFinder is not even chosen until then). By
    // the first refinement iteration it is settled, and it cannot change
    // after that, so this runs exactly once.
    if (!spr_state.runId.empty())
        return;

    bool sprMode = (params->refine_mode == REFINE_SPR);
    // buildRunId is spr_topology_test's own flag-encoding scheme, built
    // around ITS positional arguments (start-tree method, total step
    // budget); none of that describes an IQ-TREE run, so the id here is
    // the refiner plus what actually has to make it unique per run.
    ostringstream runId;
    runId << (sprMode ? "iqtree_spr" : "iqtree_nni")
          << "_" << (long) time(nullptr) << "_" << params->ran_seed;
    spr_state.runId = runId.str();
    spr_state.recordTag = sprsearch::buildRefineRecordTag(spr_opt, sprMode);
    spr_state.modelName = getModelName();

    if (spr_opt.recordProgress)
        cout << "Refinement record: appending to "
             << sprsearch::recordSpreadsheetPath(spr_state.modelName, spr_state.recordTag) << endl;
    if (spr_opt.trajectoryFlag)
        cout << "Refinement trajectory: appending to "
             << sprsearch::trajectoryTopologyPath(spr_state.runId) << endl;
}

void IQTree::doRandomSPRs() {
    // Mirror doRandomNNIs' own strength rule exactly, so --perturb means
    // the same thing whichever kick is in use: floor((ntaxa - 3) * initPS)
    // random moves, with a floor of 1 once the tree is big enough to have
    // any.
    int numMoves = (int) floor((leafNum - 3) * Params::getInstance().initPS);
    if (leafNum >= 4 && numMoves == 0)
        numMoves = 1;

    // "slack anneal": decay the kick's permitted damage toward 0 over the
    // run, so early kicks explore and the last one is strictly
    // non-worsening. Read off the shared clock (see updateSearchProgress),
    // the same one "--accept-dist ... anneal" cools against.
    double annealScale = spr_perturb_opt.slackAnneal ? (1.0 - experiment_progress) : 1.0;
    int applied = sprsearch::doRandomSPRs(*this, spr_perturb_opt, numMoves, &spr_state, annealScale);
    if (verbose_mode >= VB_MAX)
        cout << "Tree perturbation: number of random SPR performed = " << applied << endl;

    // An SPR kick needs more cleanup than the NNI one does, in two ways
    // an NNI never triggers.
    //
    // First, branch lengths. applySPR splits the target edge in half to
    // make room for the regrafted subtree (see PhyloTree::applySPR), so a
    // chain of random SPRs can halve the same region repeatedly and drive
    // lengths toward zero -- whereas an NNI swap leaves every length
    // untouched. On a weak-signal alignment, where many branches already
    // sit near zero, that produces lengths the subsequent Newton-Raphson
    // search cannot work with: it was tripping optimizeNNI's own
    // "curScore > appliedNNIs.at(0).newloglh - 0.1" assertion, which
    // asserts an applied NNI's predicted score actually held. Clamp to the
    // same floor the SPR scoring path clamps to.
    sprsearch::clampAllBranchLengthsForOptimization(*this, params->min_branch_length);
    // the clamp mutated branch lengths -- invalidate before anything reads
    // a cached partial (see the same note in maybeRunPeriodicFullReopt)
    clearAllPartialLH();

    // Second, the partial-likelihood buffers. An NNI rearranges neighbors
    // in place, so clearing the computed flags is enough; an SPR moves a
    // whole subtree, leaving the per-neighbor buffer assignment itself
    // stale. Rebuild them outright, the same reset the shared SPR scoring
    // path uses between moves.
    setAlignment(aln);
    setRootNode(params->root);

    // THE important one. Branch ids (PhyloNeighbor::id) are assigned by a
    // traversal in MTree::initializeTree, and IQ-TREE indexes by them --
    // saveBranchLengths/restoreBranchLengths write into a vector at
    // neighbor->id, so the ids must describe the CURRENT topology.
    //
    // An NNI swap preserves the branch set (it only changes which subtrees
    // hang off two existing edges), so the NNI kick leaves ids valid and
    // never had to think about this. An SPR destroys three edges and
    // creates three new ones, leaving ids that no longer correspond to the
    // tree. optimizeNNI's revert path then calls restoreBranchLengths and
    // writes every saved length onto the WRONG branch, producing a tree
    // scoring far below the NNI it just applied -- which is what tripped
    // its "curScore > appliedNNIs.at(0).newloglh - 0.1" assertion (measured
    // 76.8 logL below the prediction). Renumber before anything reads an id.
    initializeTree();

    deleteAllPartialLh();
    initializeAllPartialLh();
    clearAllPartialLH();
    if (isSuperTree())
        ((PhyloSuperTree*) this)->mapTrees();
    if (params->pll)
        pllReadNewick(getTreeString());
    resetCurScore();

    // Third difference from the NNI kick, and the one that actually
    // matters: an NNI swap leaves every branch length exactly as the last
    // refinement fitted it, so the perturbed tree still has MEANINGFUL
    // lengths. applySPR instead invents lengths for the edges it rewires
    // (it halves the target edge to make room), so after this many moves a
    // large share of the tree carries placeholder values.
    //
    // Handing that to the NNI refiner breaks one of its own invariants:
    // getBestNNIForBran predicts a post-NNI score from locally optimized
    // lengths, then optimizeNNI asserts that applying the move and running
    // one optimizeAllBranches round actually achieves it
    // ("curScore > appliedNNIs.at(0).newloglh - 0.1"). Starting from
    // placeholder lengths, one round is not enough to get there and the
    // assertion aborts the run.
    //
    // So refit the lengths here. This is not making the SPR kick easier --
    // it is giving it the same property the NNI kick already has for free.
    // A light refit, since applySPR invents lengths for the edges it
    // rewires (it halves the target edge) whereas an NNI swap leaves every
    // length as the last refinement fitted it. One sweep is enough now
    // that branch ids are correct -- the heavy 100-round fit this replaced
    // was compensating for the id bug above, not for the lengths.
    clearAllPartialLH();
    optimizeAllBranches(1);
    resetCurScore();
}

double IQTree::doSPRSearch(int blockCap, int patience) {
    finalizeRefineIdentity();
    // initRefinement could only check the parameters that force a rooted
    // tree; this is the built tree itself
    if (rooted)
        outError("--spr-refine does not support rooted trees; the SPR search is written"
                " against an unrooted topology (drop -o/--root, or refine with NNI)");

    // doTreePerturbation left the tree scrambled and re-scored but not
    // re-optimized; that score is this pass's starting point, exactly as
    // curScore is for doNNISearch.
    curScore = computeLikelihood();
    double curBestScore = getBestScore();

    sprsearch::EdgeRegistry reg;
    sprsearch::buildEdgeRegistry(*this, reg);
    spr_state.learnRadiusMaxPath = spr_opt.useDistanceRadius ? 100.0 : (double) reg.slots.size();
    // Every PhyloNode* the previous pass left behind points into a tree
    // object readTreeString has since replaced -- see beginSPRSearchPass.
    sprsearch::beginSPRSearchPass(spr_state);
    // --accept-dist's best-seen tree is per REFINEMENT: what this pass
    // walked through, not the best tree of any earlier iteration (that one
    // is already in the candidate set). Left run-wide, the restore below
    // handed most iterations an older tree back instead of their own result.
    spr_state.bestSeenScore = -DBL_MAX;
    spr_state.bestSeenTree.clear();

    sprsearch::SPRLocalLhCache localCache;
    bool haveLocalCache = !spr_opt.reoptimizeBranchLengths;
    if (haveLocalCache)
        sprsearch::allocateSPRLocalLhCache(*this, localCache);

    // "findopt"'s cadence counts KICKS here, not SPR steps, so the step
    // loop must not fire it -- recordRefineIteration does, once per
    // completed iteration. See parseRefineSpec.
    sprsearch::SPRSearchOptions stepOpt = spr_opt;
    stepOpt.findoptEveryNSteps = 0;

    string stepPrefix = "  [kick " + convertIntToString(refine_kick_count + 1) + "] ";

    if (spr_opt.stepsPerPass > 0) {
        // "steps N": an explicit budget, honoured exactly as asked
        curScore = sprsearch::runSPRSteps(*this, reg, stepOpt, spr_state, curScore,
                spr_opt.stepsPerPass, aln, *params, haveLocalCache ? &localCache : nullptr, stepPrefix);
    } else {
        // No explicit budget: climb to a local optimum, the way the NNI
        // refinement this replaces does (optimizeNNI loops until no
        // positive NNI is left). A FIXED step count would quietly
        // handicap SPR here -- the perturbation applies
        // floor((ntaxa-3) * initPS) random NNIs, 48 of them on a
        // 100-taxon tree at the default strength, and a pass that runs
        // out of steps before undoing them hands the candidate set a
        // tree far worse than the one it was kicked from. So run the
        // shared step loop in blocks and stop once a whole block fails
        // to improve, which is the same "nothing left that helps" rule,
        // with a cap so a pathological case still terminates.
        int blockSteps = max((int) leafNum - 3, 1);
        int stepCap = blockCap * blockSteps;
        // `patience` consecutive non-improving blocks before declaring
        // convergence. One is right for the NNI-equivalent case, where a
        // block exhaustively examines every candidate and so genuinely
        // proves there is nothing left. Under "fast N" a block only samples
        // N random regrafts per step, which on a large tree is a tiny
        // fraction of the neighbourhood -- a single quiet block there means
        // "this sample found nothing", not "nothing exists". Requiring
        // several in a row keeps a sparse sampler from mistaking bad luck
        // for a local optimum.
        int quietBlocks = 0;
        for (int done = 0; done < stepCap; done += blockSteps) {
            double before = curScore;
            curScore = sprsearch::runSPRSteps(*this, reg, stepOpt, spr_state, curScore, blockSteps,
                    aln, *params, haveLocalCache ? &localCache : nullptr, stepPrefix);
            if (curScore <= before + params->loglh_epsilon) {
                if (++quietBlocks >= patience)
                    break;
            } else {
                quietBlocks = 0;
            }
        }
    }

    if (haveLocalCache)
        sprsearch::freeSPRLocalLhCache(localCache);

    // Hand IQ-TREE's own optimizer a coherent likelihood cache. The SPR
    // scoring path deliberately swaps scratch buffers in and out around
    // each trial move and restores them on rollback -- an internal
    // optimization that is invisible to, and not accounted for by,
    // optimizeAllBranches, which tracks tree_lh incrementally and asserts
    // it still matches a fresh computeLikelihood to within 1.0. Crossing
    // back out of the SPR machinery is exactly the boundary at which to
    // re-establish that, the same way the step loop recomputes in full
    // before committing an improving move.
    // applySPR halves the target edge on every move, so a long run drives
    // branches toward zero -- on a degenerate alignment (SARS-CoV-2, where
    // most sites are constant) enough of them get there that IQ-TREE's own
    // optimizer cannot track its incremental tree_lh accurately and trips
    // its consistency assertion. Clamp to the same floor the periodic
    // re-optimization uses before handing the tree over.
    // --accept-dist walks downhill by design, so the stage's endpoint is
    // routinely worse than the best tree it passed through. Restore that
    // high-water mark before the model-reopt tail, so what gets reported is
    // the best tree actually found rather than wherever the walk stopped.
    if (spr_opt.acceptDist.enabled && !spr_state.bestSeenTree.empty()
            && spr_state.bestSeenScore > curScore) {
        double drifted = curScore;
        readTreeString(spr_state.bestSeenTree);
        initializeTree();
        deleteAllPartialLh();
        initializeAllPartialLh();
        clearAllPartialLH();
        curScore = computeLikelihood();
        spr_state.bestSeenRestored = true;
        spr_state.bestSeenRestores++;
        if (!spr_opt.quiet)
            cout << "  [accept-dist] restored best-seen tree: " << curScore
                 << " (walk ended at " << drifted << ")" << endl;
    }

    sprsearch::clampAllBranchLengthsForOptimization(*this, params->min_branch_length);
    clearAllPartialLH();
    curScore = computeLikelihood();

    // Same tail as doNNISearch: a genuinely better tree earns a fresh
    // model-parameter fit (the sNNI algorithm's own rule), and the same
    // checkpoint save.
    if (!on_refine_btree && curScore > curBestScore + params->modelEps) {
        optimizeModelParameters(false, params->modelEps * 10);
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            getModelFactory()->saveCheckpoint();
        }
    }
    return curScore;
}

double IQTree::doContinuousSPRStage() {
    cout << "--------------------------------------------------------------------" << endl;
    cout << "|        CONTINUOUS SPR HILL-CLIMB (--spr-continuous)              |" << endl;
    cout << "--------------------------------------------------------------------" << endl;

    // Start from the best tree the candidate-set phase produced -- the
    // same handoff --hillclimb's own "starttree" flag makes, just without
    // the round trip through a Newick file and a model string.
    readTreeString(candidateTrees.getBestTreeStrings()[0]);
    double startScore = candidateTrees.getBestScore();

    // A far larger block cap than a per-kick refinement gets: this single
    // stage is the whole search, not one of forty iterations, so it should
    // keep climbing until it genuinely cannot improve rather than hand
    // back after a fixed slice of work.
    if (spr_opt.escapeFlag) {
        curScore = runEscapingSPRStage();
        recordRefineIteration();
        double escScore = curScore;
        if (escScore > startScore) {
            addTreeToCandidateSet(getTreeString(), escScore, false, MPIHelper::getInstance().getProcessID());
            printBestCandidateTree();
        }
        curScore = escScore;
        cout << "Continuous SPR stage: logL " << startScore << " -> " << escScore
             << " (" << spr_state.stepsRun << " steps, " << spr_state.candidatesEvaluated
             << " candidates evaluated)" << endl;
    // "tunnel" accounting: entered counts marginally-rejected steps that
    // were probed from, committed counts those that actually cleared the
    // barrier. committed near zero means the flag only bought extra
    // evaluations on this data.
    if (spr_opt.tunnelFlag) {
        cout << "tunnel summary  : " << spr_state.tunnelEntered << " entered, "
             << spr_state.tunnelCommitted << " cleared, "
             << spr_state.tunnelProbes << " probe evaluation(s), tolerance "
             << spr_opt.tunnelTolerance << ", tries " << spr_opt.tunnelTries << endl;
    }
        return escScore;
    }

    // The continuous stage is the whole search rather than one of forty
    // iterations, so it gets both a far larger block cap and real patience
    // before concluding it has converged.
    //
    // Patience, not candidates-per-step, is the right knob for buying
    // thoroughness here. "fast N" scores N candidates per step but still
    // applies at most ONE move, so it pays N times over for a single
    // attempt; N cheap steps instead give N chances to actually move. This
    // was measured: raising fast from 5 to 50 halved the remaining gap but
    // cost 8.8x the time. So keep fast at 1 and spend the budget on more
    // attempts, which is what a large patience does.
    curScore = doSPRSearch(200, 20);
    recordRefineIteration();

    // Hold the result in a local across the two calls below:
    // printBestCandidateTree re-reads the best tree through readTreeString,
    // which resets curScore to PhyloTree's -DBL_MAX "not computed"
    // sentinel. Reading curScore afterwards would report that sentinel as
    // this stage's result.
    double finalScore = curScore;
    if (finalScore > startScore) {
        addTreeToCandidateSet(getTreeString(), finalScore, false, MPIHelper::getInstance().getProcessID());
        printBestCandidateTree();
    }
    curScore = finalScore;
    cout << "Continuous SPR stage: logL " << startScore << " -> " << finalScore
         << " (" << spr_state.stepsRun << " steps, " << spr_state.candidatesEvaluated
         << " candidates evaluated)" << endl;
    // "tunnel" accounting: entered counts marginally-rejected steps that
    // were probed from, committed counts those that actually cleared the
    // barrier. committed near zero means the flag only bought extra
    // evaluations on this data.
    if (spr_opt.tunnelFlag) {
        cout << "tunnel summary  : " << spr_state.tunnelEntered << " entered, "
             << spr_state.tunnelCommitted << " cleared, "
             << spr_state.tunnelProbes << " probe evaluation(s), tolerance "
             << spr_opt.tunnelTolerance << ", tries " << spr_opt.tunnelTries << endl;
    }
    return finalScore;
}

double IQTree::runEscapingSPRStage() {
    finalizeRefineIdentity();
    initializeAllPartialLh();
    clearAllPartialLH();
    curScore = computeLikelihood();

    // findopt's cadence belongs to recordRefineIteration, not the step loop
    sprsearch::SPRSearchOptions stepOpt = spr_opt;
    stepOpt.findoptEveryNSteps = 0;

    const int span = spr_opt.escapeSpan;
    const int tries = spr_opt.escapeTries;
    // total step budget: "steps N" if given, else generous
    const long stepCap = spr_opt.stepsPerPass > 0 ? (long) spr_opt.stepsPerPass
            : 200L * max((int) leafNum - 3, 1);
    // consecutive FAILED excursions before giving up -- the same role
    // patience plays in the plain continuous stage
    const int maxFailedEscapes = 20;

    long stepsUsed = 0;
    int failedEscapes = 0, escapesTried = 0, escapesKept = 0;

    while (stepsUsed < stepCap && failedEscapes < maxFailedEscapes) {
        // ---- ordinary hill-climbing for one span ----
        // The registry and the scratch cache are rebuilt each phase: every
        // rollback and every kick replaces node pointers wholesale, so
        // neither can be carried across.
        {
            sprsearch::EdgeRegistry reg;
            sprsearch::buildEdgeRegistry(*this, reg);
            spr_state.learnRadiusMaxPath = spr_opt.useDistanceRadius ? 100.0 : (double) reg.slots.size();
            sprsearch::beginSPRSearchPass(spr_state);
            sprsearch::SPRLocalLhCache cache;
            bool haveCache = !spr_opt.reoptimizeBranchLengths;
            if (haveCache)
                sprsearch::allocateSPRLocalLhCache(*this, cache);
            double before = curScore;
            curScore = sprsearch::runSPRSteps(*this, reg, stepOpt, spr_state, curScore, span,
                    aln, *params, haveCache ? &cache : nullptr, "  [climb] ");
            if (haveCache)
                sprsearch::freeSPRLocalLhCache(cache);
            stepsUsed += span;
            if (curScore > before + params->loglh_epsilon)
                continue; // still making progress -- no reason to kick
        }

        // ---- stalled: remember exactly where we are, then kick ----
        // The incumbent is stored as Newick rather than as node pointers,
        // because that is the only representation that survives the kick
        // and the rebuild that a rollback needs.
        string incumbentTree = getTreeString();
        double incumbentScore = curScore;
        escapesTried++;

        // The negative move. Whichever kick the run is configured for:
        // --spr-perturb's random SPRs, otherwise IQ-TREE's own random-NNI
        // kick, at the same floor((ntaxa-3) * --perturb) strength either
        // way -- so "how hard to kick" stays one setting across the tool.
        if (params->spr_perturb)
            doRandomSPRs();
        else
            doRandomNNIs(false);
        curScore = computeLogL();

        {
            sprsearch::EdgeRegistry reg;
            sprsearch::buildEdgeRegistry(*this, reg);
            spr_state.learnRadiusMaxPath = spr_opt.useDistanceRadius ? 100.0 : (double) reg.slots.size();
            sprsearch::beginSPRSearchPass(spr_state);
            sprsearch::SPRLocalLhCache cache;
            bool haveCache = !spr_opt.reoptimizeBranchLengths;
            if (haveCache)
                sprsearch::allocateSPRLocalLhCache(*this, cache);
            curScore = sprsearch::runSPRSteps(*this, reg, stepOpt, spr_state, curScore, tries,
                    aln, *params, haveCache ? &cache : nullptr, "  [escape] ");
            if (haveCache)
                sprsearch::freeSPRLocalLhCache(cache);
            stepsUsed += tries;
        }

        if (curScore > incumbentScore + params->loglh_epsilon) {
            // the excursion actually paid off -- keep it and carry on
            escapesKept++;
            failedEscapes = 0;
            if (!spr_opt.quiet)
                cout << "  [escape] kept: logL " << incumbentScore << " -> " << curScore << endl;
        } else {
            // it did not: revert to the incumbent exactly, as if the kick
            // had never happened, and count it against the patience budget
            readTreeString(incumbentTree);
            initializeTree();
            deleteAllPartialLh();
            initializeAllPartialLh();
            clearAllPartialLH();
            curScore = computeLikelihood();
            failedEscapes++;
            if (!spr_opt.quiet)
                cout << "  [escape] reverted (reached " << curScore << ", needed > "
                     << incumbentScore << ")" << endl;
        }
    }

    cout << "escape summary  : " << escapesTried << " excursion(s), " << escapesKept
         << " kept, " << stepsUsed << " steps used" << endl;
    return curScore;
}

/**
    One line of key=value pairs per run, to <prefix>.searchstats.txt.

    The same numbers already appear in the .log, but only inside prose that
    a sweep script has to regex its way through -- and this session added
    several search variants whose results otherwise live only there. Writing
    them in a stable, parseable form means a comparison table can be built
    from the files themselves rather than from log scraping that breaks the
    moment a message is reworded.
 */
void IQTree::writeSearchStats() {
    string path = string(params->out_prefix) + ".searchstats.txt";
    ofstream out(path.c_str());
    if (!out.is_open())
        return;
    out.precision(10);
    out << "prefix=" << params->out_prefix
        << " seed=" << params->ran_seed
        << " iterations=" << stop_rule.getCurIt()
        << " best_logl=" << candidateTrees.getBestScore()
        << " cpu_seconds=" << getCPUTime()
        << " wall_seconds=" << (getRealTime() - params->start_real_time);

    out << " refiner=" << (params->refine_mode == REFINE_SPR ? "spr" : "nni");
    out << " kick=" << (params->accept_dist ? "none"
                        : (params->spr_perturb ? "spr" : "nni"));

    out << " accept_dist=" << (accept_dist.enabled ? "on" : "off");
    if (accept_dist.enabled) {
        out << " accept_temp=" << accept_dist.temperature
            << " accept_shape=" << accept_dist.shape
            << " accept_anneal=" << (accept_dist.anneal ? "yes" : "no")
            << " accept_floor=" << accept_dist.tempFloor
            << " accept_last_temp=" << spr_state.lastTemperature
            << " downhill_kept=" << spr_state.acceptedDownhill
            << " downhill_offered=" << spr_state.offeredDownhill
            << " bestseen_restored=" << (spr_state.bestSeenRestored ? "yes" : "no")
            << " bestseen_restores=" << spr_state.bestSeenRestores;
    }

    bool slackOn = Params::getInstance().perturb_slack || spr_perturb_opt.slackFlag;
    out << " slack=" << (slackOn ? "on" : "off");
    if (slackOn) {
        double d = Params::getInstance().perturb_slack
                 ? Params::getInstance().perturb_slack_delta : spr_perturb_opt.slackDelta;
        bool an = Params::getInstance().perturb_slack
                 ? Params::getInstance().perturb_slack_anneal : spr_perturb_opt.slackAnneal;
        out << " slack_delta=" << d
            << " slack_anneal=" << (an ? "yes" : "no")
            << " slack_last_delta=" << spr_state.slackLastDelta
            << " slack_applied=" << spr_state.slackAccepted
            << " slack_refused=" << spr_state.slackRejected
            << " slack_kicks=" << spr_state.slackKicks;
    }

    if (spr_opt.tunnelFlag) {
        out << " tunnel=on tunnel_tol=" << spr_opt.tunnelTolerance
            << " tunnel_tries=" << spr_opt.tunnelTries
            << " tunnel_entered=" << spr_state.tunnelEntered
            << " tunnel_cleared=" << spr_state.tunnelCommitted
            << " tunnel_probes=" << spr_state.tunnelProbes;
    } else {
        out << " tunnel=off";
    }

    out << " candidates_evaluated=" << spr_state.candidatesEvaluated;
    out << endl;
    out.close();
    cout << "Search statistics written to " << path << endl;
}

void IQTree::recordRefineIteration() {
    // Count the completed iteration unconditionally, before the early
    // return below: the "[kick N]" step labels and findopt's cadence read
    // it. (Annealing does not -- it runs on updateSearchProgress's clock.)
    // Everything past this point (the record CSV, trajectory file,
    // findopt) only matters under --spr-refine/--nni-refine, so it stays
    // behind the guard.
    refine_kick_count++;

    if (params->refine_mode == REFINE_NNI && params->refine_spec.empty())
        return;

    finalizeRefineIdentity();

    if (spr_opt.recordProgress)
        sprsearch::appendRecordRow(spr_state.modelName, spr_state.recordTag, spr_state.runId,
                spr_state.candidatesEvaluated, getCPUTime() - spr_state.cpuClockStart, curScore,
                spr_state.trueTreeLogl, spr_opt.recordTopology, *this);
    if (spr_opt.trajectoryFlag)
        sprsearch::appendTrajectoryTopology(spr_state.runId, *this);
    // "findopt": one throwaway whole-tree ML refit every N kicks, read for
    // its logL and discarded -- the real tree, its model and curScore are
    // never touched. maybeRunFindopt's own cadence test is
    // (index + 1) % N, so the count is passed 0-based. Refitting under
    // this run's OWN model (rather than the tool's default GTR+FO/LG+FO)
    // is what keeps the headroom figure meaningful when -m names something
    // richer.
    if (spr_opt.findoptEveryNSteps > 0)
        sprsearch::maybeRunFindopt(*this, refine_kick_count - 1, spr_opt.findoptEveryNSteps,
                spr_opt.quiet, spr_opt.recordProgress, spr_opt.recordTopology, spr_state.modelName,
                spr_state.recordTag, spr_state.runId, spr_state.candidatesEvaluated,
                spr_state.cpuClockStart, spr_state.trueTreeLoglForFindopt, curScore, aln, *params,
                getModelName());
}

/****************************************************************************
 Hooks called from IQ-TREE's own search loop (tree/iqtree.cpp)

 Each call site there is marked "[search-experiments hook]". Every hook is a
 no-op -- no output, no random draw, no change to the tree -- unless its
 experiment flag was given, so stock IQ-TREE runs exactly as before.
 ****************************************************************************/

void IQTree::initSearchExperimentState() {
    // spr_opt/spr_state default-construct themselves; initRefinement fills
    // them in only when --spr-refine/--nni-refine was actually given
    refine_kick_count = 0;
    experiment_progress = 0.0;
    search_loop_start = -1;
}

void IQTree::beginSearchExperiments() {
    search_loop_start = -1;
    experiment_progress = 0.0;
    if (params->spr_continuous)
        doContinuousSPRStage();
}

void IQTree::updateSearchProgress() {
    // The one annealing clock for every experiment that anneals: 0 on the
    // first main-loop iteration, 1 on the last. It used to be the raw
    // iteration count over min_iterations (accept-dist) or the kick count
    // over min_iterations (slack) -- two different clocks, neither of which
    // reached 1: starting-tree iterations counted toward one and not the
    // other, and under the default stopping rule min_iterations is
    // 100 x ntaxa, so neither annealed at all.
    int it = stop_rule.getCurIt();
    if (search_loop_start < 0)
        search_loop_start = it;
    // the stop_rule iteration index the loop will END on
    int lastIt;
    if (params->stop_condition == SC_FIXED_ITERATION) {
        // -n N: the loop runs while curIt < N
        lastIt = params->min_iterations - 1;
    } else {
        // The default rule (and the ones built on the same "unsuccessful
        // iterations" test) stops --nstop iterations after the last
        // improvement. Aim at that projected end; it moves out whenever the
        // search improves, so the schedule warms back up slightly then.
        lastIt = stop_rule.getLastImprovedIteration() + params->unsuccess_iteration;
    }
    int span = lastIt - search_loop_start;
    double p = (span > 0) ? (double) (it - search_loop_start) / (double) span : 1.0;
    if (p < 0.0)
        p = 0.0;
    if (p > 1.0)
        p = 1.0;
    experiment_progress = p;
    if (accept_dist.enabled) {
        spr_opt.acceptDist = accept_dist;
        spr_opt.acceptProgressBase = experiment_progress;
    }
}

bool IQTree::doExperimentPerturbation() {
    if (Params::getInstance().accept_dist) {
        // --accept-dist replaces the kick outright: escaping a local
        // optimum is the acceptance rule's job now, and running a kick as
        // well would make it impossible to say which mechanism did the
        // escaping. The candidate tree already read is left exactly as drawn.
        return true;
    }
    if (Params::getInstance().spr_perturb) {
        // --spr-perturb: random SPR moves instead of random NNIs, the same
        // count the NNI kick would have used so --perturb still means the
        // same thing. See doRandomSPRs.
        doRandomSPRs();
        return true;
    }
    return false;
}

void IQTree::beginNNIKickSlack(sprsearch::KickSlack &slack) {
    // --perturb-slack: bound how far one kick move may drive the tree
    // downhill, the same rule --spr-perturb's "slack" applies to the SPR
    // kick. Without it the NNI kick is blind, exactly as stock IQ-TREE's.
    slack = sprsearch::KickSlack();
    slack.active = Params::getInstance().perturb_slack;
    if (!slack.active)
        return;
    slack.delta = Params::getInstance().perturb_slack_delta;
    if (Params::getInstance().perturb_slack_anneal)
        slack.delta *= (1.0 - experiment_progress);
    clearAllPartialLH();
    slack.curScore = computeLikelihood();
}

sprsearch::KickSlack::Verdict IQTree::judgeNNIKickMove(sprsearch::KickSlack &slack, NNIMove &move) {
    // score what the kick just did; undo it and redraw if it costs more
    // than the bound allows
    clearAllPartialLH();
    double after = computeLikelihood();
    if (std::isfinite(after) && after > slack.curScore - slack.delta) {
        slack.curScore = after;
        slack.kept++;
        slack.draws = 0;
        return sprsearch::KickSlack::KEEP;
    }
    doNNI(move);            // an NNI is its own inverse
    clearAllPartialLH();
    slack.refused++;
    if (slack.draws < sprsearch::KickSlack::MAX_DRAWS) {
        slack.draws++;
        return sprsearch::KickSlack::REDRAW;    // same slot, fresh draw
    }
    slack.draws = 0;
    return sprsearch::KickSlack::GIVE_UP;       // give up on this slot
}

void IQTree::endNNIKickSlack(const sprsearch::KickSlack &slack) {
    if (!slack.active)
        return;
    spr_state.slackAccepted += slack.kept;
    spr_state.slackRejected += slack.refused;
    spr_state.slackKicks++;
    spr_state.slackLastDelta = slack.delta;
}

bool IQTree::acceptNNICandidate(double newLogl, double curLogl) {
    // "record"'s x-axis, kept in the same unit for both refiners:
    // getBestNNIForBran scores BOTH alternative topologies around this
    // branch, so one NNI candidate is two evaluations -- the same thing
    // scoreTrialSPRMove counts one of per trial regraft.
    spr_state.candidatesEvaluated += 2;
    if (newLogl > curLogl)
        return true;
    if (!accept_dist.enabled)
        return false;
    // --accept-dist: an NNI that loses ground may still be admitted, with
    // probability set by the rule. This is the NNI refiner's accept/reject
    // point, the exact counterpart of the SPR step loop's own.
    spr_state.offeredDownhill++;
    spr_state.lastTemperature = accept_dist.effectiveTemperature(experiment_progress);
    if (accept_dist.accept(newLogl - curLogl, experiment_progress)) {
        spr_state.acceptedDownhill++;
        return true;
    }
    return false;
}

void IQTree::beginNNIAcceptWalk(sprsearch::NNIAcceptState &walk) {
    walk = sprsearch::NNIAcceptState();
    walk.bestScore = curScore;
    walk.bestTree = getTreeString();
}

bool IQTree::nniAcceptWalkShouldStop(sprsearch::NNIAcceptState &walk, double newScore, double oldScore) {
    if (newScore > walk.bestScore) {
        walk.bestScore = newScore;
        walk.bestTree = getTreeString();
    }
    if (accept_dist.effectiveTemperature(experiment_progress) > 0.0) {
        // A step is ALLOWED to lose ground, so the usual "stop as soon as a
        // step stops improving" test would end the walk on its first
        // tolerated move. Allow a bounded run of non-improving steps
        // instead. The bound matters: without it the loop runs to MAXSTEPS
        // (= leafNum) every iteration, which measured ~8x the CPU of a plain
        // NNI iteration on a 500-taxon alignment for no benefit.
        if (newScore - oldScore < params->loglh_epsilon)
            return ++walk.stall >= sprsearch::NNIAcceptState::STALL_LIMIT;
        walk.stall = 0;
        return false;
    }
    // annealed all the way to T = 0: a strict hill-climb again
    return newScore - oldScore < params->loglh_epsilon;
}

void IQTree::endNNIAcceptWalk(sprsearch::NNIAcceptState &walk) {
    // The walk is allowed downhill, so where it stopped is routinely worse
    // than the best tree it passed through. Hand back that high-water mark
    // -- the same restore doSPRSearch does for the SPR refiner. Before this
    // existed the NNI refiner handed the candidate set wherever the walk
    // ended.
    if (walk.bestTree.empty() || walk.bestScore <= curScore + params->loglh_epsilon)
        return;
    double drifted = curScore;
    readTreeString(walk.bestTree);
    initializeTree();
    deleteAllPartialLh();
    initializeAllPartialLh();
    clearAllPartialLH();
    if (isSuperTree())
        ((PhyloSuperTree*) this)->mapTrees();
    curScore = computeLikelihood();
    spr_state.bestSeenRestored = true;
    spr_state.bestSeenRestores++;
    if (verbose_mode >= VB_MED)
        cout << "[accept-dist] restored best-seen NNI tree: " << curScore
             << " (walk ended at " << drifted << ")" << endl;
}

void IQTree::finishSearchExperiments() {
    // "sweep": an end-of-run phase in spr_topology_test too, so it stays
    // one here -- a single exhaustive whole-tree regraft pass over the
    // least-compatible sibling pairs of the FINAL tree, not something run
    // per perturbation. Anything it improves goes back into the candidate
    // set so the rest of the pipeline sees it.
    if (params->refine_mode == REFINE_SPR && spr_opt.sweepFlag) {
        curScore = computeLikelihood();
        sprsearch::EdgeRegistry sweepReg;
        sprsearch::buildEdgeRegistry(*this, sweepReg);
        sprsearch::beginSPRSearchPass(spr_state);
        double sweptScore = sprsearch::runSPRSweep(*this, sweepReg, spr_opt, spr_state, curScore);
        if (sweptScore > curScore) {
            curScore = sweptScore;
            addTreeToCandidateSet(getTreeString(), curScore, false, MPIHelper::getInstance().getProcessID());
            printBestCandidateTree();
        }
        readTreeString(candidateTrees.getBestTreeStrings()[0]);
    }
    // One last "record" row reflecting the final state, whether or not the
    // last iteration happened to improve on anything -- see
    // appendFinalRecordRow for why this is unconditional. The score comes
    // from the candidate set, not curScore: readTreeString above swaps in
    // the best tree without recomputing its likelihood, leaving curScore
    // at PhyloTree's own -DBL_MAX "not computed" sentinel.
    if (params->refine_mode == REFINE_SPR || !params->refine_spec.empty()) {
        finalizeRefineIdentity();
        sprsearch::appendFinalRecordRow(spr_state, spr_opt, *this, candidateTrees.getBestScore());
    }
}

void IQTree::reportSearchExperiments() {
    writeSearchStats();

    if (accept_dist.enabled) {
        cout << "accept-dist     : " << spr_state.acceptedDownhill << " downhill move(s) kept of "
             << spr_state.offeredDownhill << " offered; T0 " << accept_dist.temperature
             << " shape " << accept_dist.shape
             << (accept_dist.anneal ? " annealed" : " constant")
             << ", last effective T " << spr_state.lastTemperature
             << ", best-seen restore ";
        if (spr_state.bestSeenRestores > 0)
            cout << "fired in " << spr_state.bestSeenRestores << " refinement(s)";
        else
            cout << "not needed";
        cout << endl << endl;
    }
    // Report the settings the kick actually ran with: --perturb-slack's own
    // values when it was given (NNI kick, or propagated to the SPR kick by
    // initRefinement), otherwise the SPR kick's "slack D" spec. Reading
    // spr_perturb_opt unconditionally printed "delta 0.000 constant" for
    // every --perturb-slack run.
    bool viaFlag = Params::getInstance().perturb_slack;
    if (viaFlag || spr_perturb_opt.slackFlag) {
        double delta = viaFlag ? Params::getInstance().perturb_slack_delta : spr_perturb_opt.slackDelta;
        bool anneal = viaFlag ? Params::getInstance().perturb_slack_anneal : spr_perturb_opt.slackAnneal;
        cout << "slack summary   : " << spr_state.slackAccepted << " kick move(s) applied, "
             << spr_state.slackRejected << " refused, over " << spr_state.slackKicks
             << " kick(s); delta " << delta
             << (anneal ? " annealed" : " constant")
             << ", last effective " << spr_state.slackLastDelta << endl << endl;
    }
}
