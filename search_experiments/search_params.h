/***************************************************************************
 *   search_experiments/search_params.h                                   *
 *                                                                         *
 *   Command-line settings for the search-experiment flags. IQ-TREE's      *
 *   Params (utils/tools.h) inherits SearchExperimentParams, so every      *
 *   field is still read as params.<field> / Params::getInstance().<field> *
 *   -- but none of them is declared in IQ-TREE's own Params.              *
 *   See search_experiments/README.md.                                     *
 ***************************************************************************/

#ifndef SEARCH_EXPERIMENTS_SEARCH_PARAMS_H
#define SEARCH_EXPERIMENTS_SEARCH_PARAMS_H

#include <string>
#include <ostream>

/**
        Which hill-climber refines the tree after each perturbation of the
        stochastic search. See Params::refine_mode.
 */
enum RefineMode {
    REFINE_NNI, REFINE_SPR
};

struct SearchExperimentParams {

    /**
     *  Which algorithm hill-climbs back up after each perturbation ("kick")
     *  of the stochastic search -- the refinement stage that follows
     *  doTreePerturbation() in IQTree::doTreeSearch. REFINE_NNI is
     *  IQ-TREE's own long-standing NNI search (doNNISearch), unchanged and
     *  still the default; REFINE_SPR runs the SPR hill-climber shared with
     *  spr_topology_test --hillclimb instead (IQTree::doSPRSearch, see
     *  search_experiments/sprsearch.h).
     */
    RefineMode refine_mode;

    /**
     *  Whether the PERTURBATION half of the stochastic search applies
     *  random SPR moves instead of random NNIs (--spr-perturb). Entirely
     *  orthogonal to refine_mode: an SPR kick can be paired with either
     *  refiner, and the kick's STRENGTH is still --perturb/initPS either
     *  way (the same floor((ntaxa-3) * initPS) count of random moves).
     *
     *  It exists because the default NNI kick and the NNI refiner are
     *  matched in a way that flatters NNI: undoing N random NNIs is
     *  precisely what an exhaustive NNI scan is built for. An SPR kick
     *  removes that asymmetry.
     */
    bool spr_perturb;

    /**
     *  --spr-perturb's own flag string -- the subset of the shared SPR
     *  vocabulary that describes a MOVE rather than a search: radius,
     *  distradius, weightprune [long|short]. Parsed by
     *  sprsearch::parseRefineSpec in "perturbation" mode.
     */
    std::string spr_perturb_spec;

    /**
     *  --spr-continuous: skip the perturb/refine loop entirely. Build the
     *  candidate tree set as usual (IQ-TREE's own NNI machinery), then run
     *  ONE long SPR hill-climb from the best candidate tree and stop --
     *  the in-process equivalent of running iqtree3 for its starting
     *  candidate set and handing that tree to spr_topology_test
     *  --hillclimb. The bet is that SPR's own long-range moves escape
     *  local optima directly, making the stochastic restart loop
     *  unnecessary rather than something SPR has to keep recovering from.
     *
     *  Uses refine_spec for its settings, so it implies REFINE_SPR.
     */
    bool spr_continuous;

    /**
     *  --spr-refine/--nni-refine's own flag string, in exactly the
     *  vocabulary spr_topology_test --hillclimb takes its trailing flags
     *  in (e.g. "radius 10 fast quiet learnradius 30"). Stored raw here
     *  and parsed by sprsearch::parseRefineSpec, called from
     *  IQTree::initSettings -- the parser lives with the search it
     *  configures so the two can never drift apart. Empty means neither
     *  flag was given.
     */
    std::string refine_spec;

    /**
     *  --accept-dist "<spec>": replace the refinement's strict improve-only
     *  rule with a stochastic one -- see sprsearch.h's AcceptDist. A move
     *  losing |d| log-likelihood is kept with probability
     *  exp(-(|d|/T)^shape); an improvement is always kept.
     *
     *  Giving this flag also SUPPRESSES the perturbation stage entirely.
     *  The two are alternative ways to leave a local optimum -- a kick
     *  jumps out and climbs back, whereas this walks out one tolerated
     *  step at a time -- and running both at once would make it impossible
     *  to attribute an escape to either. See Params::accept_dist_spec for
     *  the spec vocabulary.
     *
     *  Works with the NNI refiner and the SPR refiner alike: the rule is
     *  applied at each one's own accept/reject point.
     */
    bool accept_dist;
    std::string accept_dist_spec;

    /**
     *  --perturb-slack D [anneal]: bound how far a single PERTURBATION move
     *  may drive the tree downhill, whichever kick is in use. This is the
     *  kick-agnostic spelling of --spr-perturb's own "slack" flag, and the
     *  only way to get the behaviour with IQ-TREE's default NNI kick.
     */
    bool perturb_slack;
    double perturb_slack_delta;
    bool perturb_slack_anneal;

    /** every field above at its "flag not given" value (from Params::setDefault) */
    void setSearchExperimentDefaults();
};

/**
 *  Parse one search-experiment flag at argv[cnt], advancing cnt past any
 *  argument it takes. Returns false, touching nothing, if argv[cnt] is not
 *  one of these flags. Throws const char* on a malformed flag, exactly like
 *  parseArg itself (which calls this).
 */
bool parseSearchExperimentArg(int argc, char *argv[], int &cnt, SearchExperimentParams &params);

/** the TREE SEARCH ALGORITHM usage lines for these flags (from usage_iqtree) */
void usageSearchExperiments(std::ostream &out);

#endif
