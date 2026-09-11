/***************************************************************************
 *   search_experiments/search_params.cpp                                 *
 *                                                                         *
 *   Defaults, parsing and usage text for the search-experiment flags --   *
 *   called from utils/tools.cpp (Params::setDefault, parseArg,            *
 *   usage_iqtree) at the lines marked "[search-experiments hook]".        *
 ***************************************************************************/

#include "search_params.h"
#include "utils/tools.h"
#include <cstring>

void SearchExperimentParams::setSearchExperimentDefaults() {
    refine_mode = REFINE_NNI;
    refine_spec = "";
    spr_perturb = false;
    spr_perturb_spec = "";
    spr_continuous = false;
    accept_dist = false;
    accept_dist_spec = "";
    perturb_slack = false;
    perturb_slack_delta = 0.0;
    perturb_slack_anneal = false;
}

bool parseSearchExperimentArg(int argc, char *argv[], int &cnt, SearchExperimentParams &params) {
    // --spr-refine / --nni-refine: which hill-climber runs after
    // each perturbation, plus that climber's own settings. The
    // argument is passed through verbatim in spr_topology_test
    // --hillclimb's own trailing-flag vocabulary rather than being
    // exploded into ~25 separate IQ-TREE flags, so the two tools
    // can never drift apart; sprsearch::parseRefineSpec does the
    // actual parsing (and the rejecting), called from
    // IQTree::initSettings. An empty string is legal and means
    // "this mode, all defaults".
    if (strcmp(argv[cnt], "--spr-refine") == 0) {
        cnt++;
        if (cnt >= argc)
            throw "Use --spr-refine \"<flags>\" (e.g. --spr-refine \"radius 10 fast quiet\")";
        params.refine_mode = REFINE_SPR;
        params.refine_spec = argv[cnt];
        return true;
    }

    // --spr-perturb: swap the perturbation half of the search
    // from random NNIs to random SPRs. Orthogonal to
    // --spr-refine/--nni-refine, which choose the REFINEMENT half.
    if (strcmp(argv[cnt], "--spr-continuous") == 0) {
        cnt++;
        if (cnt >= argc)
            throw "Use --spr-continuous \"<flags>\" (e.g. --spr-continuous \"radius 10 fast 5 distradius investigate\")";
        params.refine_mode = REFINE_SPR;
        params.spr_continuous = true;
        params.refine_spec = argv[cnt];
        return true;
    }

    if (strcmp(argv[cnt], "--spr-perturb") == 0) {
        params.spr_perturb = true;
        // the spec is optional here, unlike --spr-refine's: a bare
        // --spr-perturb is the common case (all defaults)
        if (cnt + 1 < argc && argv[cnt + 1][0] != '-') {
            cnt++;
            params.spr_perturb_spec = argv[cnt];
        }
        return true;
    }

    // --accept-dist: stochastic acceptance in the REFINEMENT,
    // and no perturbation at all. See Params::accept_dist.
    if (strcmp(argv[cnt], "--accept-dist") == 0) {
        params.accept_dist = true;
        // the spec is optional: a bare --accept-dist means all
        // defaults (temp 0.5, shape 1, no annealing)
        if (cnt + 1 < argc && argv[cnt + 1][0] != '-') {
            cnt++;
            params.accept_dist_spec = argv[cnt];
        }
        return true;
    }

    // --perturb-slack: the kick-agnostic spelling of --spr-perturb's
    // "slack", so the same bound can be applied to the NNI kick.
    if (strcmp(argv[cnt], "--perturb-slack") == 0) {
        cnt++;
        if (cnt >= argc)
            throw "Use --perturb-slack <delta> [anneal]";
        params.perturb_slack = true;
        params.perturb_slack_delta = convert_double(argv[cnt]);
        if (params.perturb_slack_delta <= 0.0)
            throw "--perturb-slack delta must be positive";
        if (cnt + 1 < argc && strcmp(argv[cnt + 1], "anneal") == 0) {
            cnt++;
            params.perturb_slack_anneal = true;
        }
        return true;
    }

    if (strcmp(argv[cnt], "--nni-refine") == 0) {
        cnt++;
        if (cnt >= argc)
            throw "Use --nni-refine \"<flags>\" (e.g. --nni-refine \"record\")";
        params.refine_mode = REFINE_NNI;
        params.refine_spec = argv[cnt];
        return true;
    }
    return false;
}

void usageSearchExperiments(std::ostream &out) {
    out
    << "  --spr-refine \"FLAGS\" Refine each perturbed tree by SPR hill-climbing instead" << endl
    << "                       of NNI. FLAGS uses spr_topology_test --hillclimb's own" << endl
    << "                       vocabulary, e.g. \"radius 10 fast quiet learnradius 30\"" << endl
    << "                       (see search_experiments/spr_topology_test_usage.txt)" << endl
    << "  --nni-refine \"FLAGS\" Keep the default NNI refinement, but enable the" << endl
    << "                       refiner-independent FLAGS, e.g. \"record findopt 25\"" << endl
    << "  --spr-continuous \"FLAGS\" Skip the perturb/refine loop: build the candidate set," << endl
    << "                       then run ONE long SPR hill-climb from its best tree" << endl
    << "  --spr-perturb [\"FLAGS\"] Perturb with random SPR moves instead of random NNIs." << endl
    << "                       FLAGS may set radius/distradius/weightprune [long|short]" << endl
    << "  --perturb-slack D [anneal]  Bound one kick move's damage to D logL (any kick)." << endl
    << "  --accept-dist [\"SPEC\"]  Stochastic acceptance in refinement; no perturbation." << endl
    << "                       SPEC: temp T | shape S | anneal | floor F" << endl;
}
