/***************************************************************************
 *   Standalone executable that samples random 4-taxon quartets from a   *
 *   real alignment and, for each one, finds which of the 3 possible     *
 *   unrooted quartet topologies has the highest likelihood.             *
 *                                                                        *
 *   Does almost none of this itself: PhyloTree::computeQuartetLikelihoods*
 *   (tree/quartet.cpp) -- IQ-TREE's own likelihood-mapping engine, the   *
 *   same code "--lmap" uses -- already draws random quartets with all 4 *
 *   taxa guaranteed distinct, builds each of the 3 possible topologies  *
 *   ((0,1)|(2,3), (0,2)|(1,3), (0,3)|(1,2)) via readTreeStringSeqName,  *
 *   and branch-length-optimizes each one's likelihood via               *
 *   optimizeAllBranches. This tool is just the setup (load an alignment,*
 *   fit a shared starting model/tree for computeQuartetLikelihoods to   *
 *   reuse) and the output (pick each quartet's best-scoring topology,   *
 *   write the list to a file) around that one call.                     *
 *
 *   This is a plain executable with its own main() -- it does not run a *
 *   real IQ-TREE analysis pipeline, and shares no code with (only the   *
 *   same house style as) tree/spr_topology_test.cpp, its nearest sibling*
 *   in this directory.                                                  *
 *   Build target: quartet_topology_test (see root CMakeLists.txt).      *
 ***************************************************************************/

#include "phylotree.h"
#include "alignment/alignment.h"
#include "model/modelfactory.h"
#include "utils/timeutil.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// This target deliberately does not link the `main` library (it defines
// int main(), which would conflict with the one below). A handful of
// symbols are nonetheless referenced from lower-level libraries this tool
// does link (utils/tree/alignment/model) even though they're only ever
// defined in main/*.cpp; this tool never exercises the code paths that
// call them, so trivial stubs are enough to satisfy the linker -- same
// stubs, for the same reason, as tree/spr_topology_test.cpp's own.
void printCopyright(ostream &out) {}
string detectSeqTypeName(string model_name) { return ""; }
void reportRate(ostream &out, PhyloTree &tree) {}
const char *aa_model_names_rax[] = {"LG", "WAG", "JTT", "JTTDCMut", "DCMut", "VT", "PMB", "Blosum62", "Dayhoff",
        "mtREV", "mtART", "mtZOA", "mtMAM",
        "HIVb", "HIVw", "FLU", "rtREV", "cpREV"};

/**
    true iff `seqType` is one this tool's model handling actually
    supports. Alignment content can auto-detect to other SeqType values
    (SEQ_BINARY, SEQ_MORPH, SEQ_CODON, ...) this tool has no model logic
    for; call this right after loading an alignment and bail out with a
    clear error instead of silently mishandling one of those. Same
    definition as spr_topology_test.cpp's own isSupportedSeqType, kept
    separate since these are two independent translation units.
 */
bool isSupportedSeqType(SeqType seqType) {
    return seqType == SEQ_DNA || seqType == SEQ_PROTEIN;
}

/**
    the shared starting model every sampled quartet's own 3 candidate
    topologies get scored under (computeQuartetLikelihoods reuses
    `this`'s model/rate object directly -- see its own comment in
    tree/quartet.cpp -- rather than refitting frequencies/rate parameters
    per quartet, only re-optimizing each candidate topology's own branch
    lengths). Deliberately richer than spr_topology_test.cpp's own
    modelNameFor default (which only ever offers plain LG/JC or LG+FO/
    GTR+FO, no rate heterogeneity at all): this tool has no other model
    to fall back on if this one is a poor fit, and real biological data
    is routinely rate-heterogeneous enough that a from-scratch model
    lacking +G materially changes which topology looks best (see
    PARTIAL_LIKELIHOOD_ATTEMPT.md-adjacent findings elsewhere in this
    project) -- so a Gamma-rate term is included by default here, not
    left as something the caller has to remember to ask for. Still fully
    overridable via --model for a literal external spec (e.g. carrying
    over a real iqtree3 run's own fitted values the same way
    spr_topology_test.cpp's "model <spec>" flag and
    test_scripts/extract_iqtree_model.py do).
 */
string defaultModelFor(SeqType seqType) {
    if (seqType == SEQ_PROTEIN)
        return "LG+FO+G4";
    return "GTR+FO+G4";
}

void printUsage(const char *progName) {
    cerr << "Usage: " << progName << " <alignment> <num-quartets> [options]" << endl;
    cerr << endl;
    cerr << "Samples <num-quartets> random 4-taxon quartets from <alignment> (all 4 taxa" << endl;
    cerr << "in a quartet always distinct), finds which of the 3 possible unrooted" << endl;
    cerr << "quartet topologies has the highest likelihood for each one, and writes the" << endl;
    cerr << "list to a file." << endl;
    cerr << endl;
    cerr << "  <alignment>       PHYLIP/FASTA/NEXUS/CLUSTAL/MSF alignment file (DNA or protein," << endl;
    cerr << "                    sequence type auto-detected)" << endl;
    cerr << "  <num-quartets>    how many quartets to sample. May sample the same quartet" << endl;
    cerr << "                    more than once (with replacement across quartets --" << endl;
    cerr << "                    only the 4 taxa WITHIN one quartet are guaranteed" << endl;
    cerr << "                    distinct from each other), same as --lmap" << endl;
    cerr << endl;
    cerr << "  out <path>        output file (default: <alignment>.quartets.tsv)" << endl;
    cerr << "  seed <N>          random seed, for reproducible sampling (default: time-based)" << endl;
    cerr << "  model <spec>      literal IQ-TREE model spec (e.g. \"LG+F{0.081,...}+G4{1.2}\")" << endl;
    cerr << "                    used as the shared model every quartet's 3 candidate" << endl;
    cerr << "                    topologies are scored under, INSTEAD of the default" << endl;
    cerr << "                    LG+FO+G4/GTR+FO+G4 fit fresh from this alignment -- the same" << endl;
    cerr << "                    fixed-parameter-spec convention spr_topology_test.cpp's own" << endl;
    cerr << "                    \"model <spec>\" flag and test_scripts/extract_iqtree_model.py use" << endl;
    cerr << "  threads <N>       number of threads for the shared model's own initial fit and" << endl;
    cerr << "                    (OpenMP builds only) quartet likelihood computation (default: 1)." << endl;
    cerr << "                    N>1 has NOT been verified safe -- concurrent quartets share the" << endl;
    cerr << "                    same fitted model object, which segfaulted under OpenMP's own" << endl;
    cerr << "                    default (one thread per CPU core) in testing; only 1 is known-good" << endl;
    cerr << "  quiet             suppress the underlying computeQuartetLikelihoods progress dots" << endl;
}

/**
    build a reasonable shared starting tree/model for computeQuartetLikelihoods
    to reuse across every sampled quartet: a BioNJ topology (computeDist +
    computeBioNJ, IQ-TREE's own standard cheap-distance starting tree) with
    its branch lengths and model parameters (frequencies, rate
    heterogeneity) jointly ML-fit via ModelFactory::optimizeParameters.
    The BioNJ topology itself is thrown away in every sense that matters
    to this tool -- computeQuartetLikelihoods never touches `tree`'s own
    topology, it only reuses its fitted model/rate object and builds each
    quartet's own 3 candidate trees from scratch -- so unlike
    spr_topology_test.cpp's own starting-tree code, there's no need for
    BioNJ's second-pass ML-distance refinement (computeDist -> computeBioNJ
    -> fit -> computeDist again -> computeBioNJ again): that exists purely
    to improve topology quality, and this tool never scores this
    particular topology at all.
 */
void fitSharedModel(PhyloTree &tree, Alignment *aln, Params &params, const string &modelSpec, int threads) {
    tree.setParams(&params);
    // params.out_prefix is a raw char*, left nullptr by Params::setDefault()
    // (a normal analysis pipeline always sets it from -s/--prefix before
    // reaching here); PhyloTree::computeDist does an unconditional
    // "dist_file = params.out_prefix" (tree/phylotree.cpp), which throws
    // "basic_string: construction from null is not valid" without this --
    // confirmed by hitting exactly that crash first. Same placeholder-
    // literal-cast-to-char* workaround spr_topology_test.cpp's own
    // "spr_hillclimb_tmp" uses; the resulting stray distance-file name is
    // never read back by anything here
    params.out_prefix = (char*) "quartet_topology_test_tmp";
    // cheap, model-free pairwise distances (PhyloTree::computeDist falls
    // back to Alignment::computeDist whenever no model/rate is attached
    // yet -- unconditionally true here, `tree` has no model at all).
    // Needs tree.aln already set (via the PhyloTree(aln) constructor the
    // caller used) even though setAlignment() below hasn't run yet --
    // computeDist/computeBioNJ crash (std::string built from a null
    // char*) without it; confirmed by testing the default constructor
    // first and hitting exactly that
    tree.computeDist(params, aln, tree.dist_matrix, tree.var_matrix);
    tree.computeBioNJ(params);
    // re-map leaf ids to match the alignment's own sequence order/names
    tree.setAlignment(aln);
    // BioNJ (a distance-based method) can and does produce a zero or
    // slightly negative branch length on some edge; optimizeOneBranch
    // (called on every edge by optimizeParameters below) asserts its
    // branch's CURRENT length is >= 0 before searching for a better one,
    // so this has to run first -- confirmed by hitting exactly that
    // assertion without it. Same built-in fix
    // computeQuartetLikelihoods itself calls on each of ITS OWN 3
    // candidate quartet trees, right after building them (tree/quartet.cpp)
    tree.wrapperFixNegativeBranch(true);

    ModelsBlock *modelsBlock = readModelsDefinition(params);
    tree.setModelFactory(new ModelFactory(params, const_cast<string&>(modelSpec), &tree, modelsBlock));
    delete modelsBlock;
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.setNumThreads(threads);
    tree.setLikelihoodKernel(LK_SSE2);
    tree.initializeAllPartialLh();
    tree.getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, params.modelEps);
}

int runQuartetSearch(const string &alignmentFile, int64_t numQuartets, const string &outPath, int seed,
        bool useModelOverride, const string &modelOverrideSpec, int threads, bool quiet) {
    ifstream check(alignmentFile.c_str());
    if (!check.good()) {
        cerr << "error: cannot open alignment file '" << alignmentFile << "'" << endl;
        return 2;
    }
    check.close();

    if (numQuartets <= 0) {
        cerr << "error: <num-quartets> must be a positive integer, got " << numQuartets << endl;
        return 2;
    }

    Params &params = Params::getInstance();
    params.setDefault();
    // setDefault() leaves SSE at LK_AVX512, assuming a full-featured
    // build; this target compiles with -D__NOAVX__ (see its CMakeLists.txt
    // block), which doesn't wire up a real function pointer for that
    // kernel -- computeQuartetLikelihoods calls
    // quartet_tree->setLikelihoodKernel(params->SSE) internally on each of
    // its own per-quartet trees (tree/quartet.cpp), so left at the
    // AVX512 default this null-function-pointer-calls and segfaults
    // (confirmed via gdb backtrace: crashes inside
    // PhyloTree::optimizeAllBranches's computeLikelihoodBranch call).
    // LK_SSE2 is the same kernel every other standalone tool in this
    // directory (e.g. spr_topology_test.cpp's own tree.setLikelihoodKernel
    // calls) hardcodes for exactly this reason
    params.SSE = LK_SSE2;

    InputType intype;
    // nullptr sequence_type makes the Alignment constructor auto-detect
    // DNA vs. protein vs. other types from raw character content
    Alignment *aln = new Alignment((char*) alignmentFile.c_str(), nullptr, intype, "");
    if (!isSupportedSeqType(aln->seq_type)) {
        cerr << "error: alignment '" << alignmentFile << "' auto-detected as a sequence type this tool "
                "doesn't support (only DNA and protein are handled)" << endl;
        delete aln;
        return 2;
    }
    if (aln->getNSeq() < 4) {
        cerr << "error: alignment '" << alignmentFile << "' has only " << aln->getNSeq()
             << " sequence(s); quartets need at least 4" << endl;
        delete aln;
        return 2;
    }
    // computeQuartetLikelihoods itself guards exactly this condition for
    // each quartet's own sub-alignment (tree/quartet.cpp: "if
    // (quartet_aln->ordered_pattern.empty()) ...orderPatternByNumChars");
    // the parsimony-tree/branch-length code this tool's own shared-model
    // setup goes through below needs the SAME done on the main alignment
    // first -- confirmed by hitting PhyloTree::computePartialParsimonyFast's
    // own "!aln->ordered_pattern.empty()" assertion without it
    if (aln->ordered_pattern.empty())
        aln->orderPatternByNumChars(PAT_VARIANT);

    // seeded here, before computeQuartetLikelihoods' own sampling: the
    // OpenMP build path inside it seeds each thread from
    // Params::ran_seed directly (params->ran_seed + omp_get_thread_num()),
    // while the non-OpenMP path uses the global randstream init_random()
    // sets up -- neither is set by the other, so both need setting for
    // reproducible --seed behavior regardless of how this was built
    int actualSeed = (seed >= 0) ? seed : (int) (time(nullptr) * 1000 + (long) (getRealTime() * 1000) % 1000);
    params.ran_seed = actualSeed;
    init_random(actualSeed, !quiet);

    string modelSpec = useModelOverride ? modelOverrideSpec : defaultModelFor(aln->seq_type);

    PhyloTree tree(aln);
    fitSharedModel(tree, aln, params, modelSpec, threads);

    if (!quiet) {
        cout << "sequence type   : " << (aln->seq_type == SEQ_PROTEIN ? "protein" : "DNA")
             << " (auto-detected)" << endl;
        cout << "alignment       : " << alignmentFile << " (" << aln->getNSeq() << " sequences, "
             << aln->getNSite() << " sites)" << endl;
        cout << "model           : " << modelSpec
             << (useModelOverride ? " (externally supplied, fixed)" : " (fit fresh on this alignment)") << endl;
        cout << "quartets        : " << numQuartets << " (random, all 4 taxa in each guaranteed distinct)" << endl;
        cout << "seed            : " << actualSeed << endl;
        cout << "out             : " << outPath << endl;
        cout << endl;
    }

    // groups.numGroups == 0 is computeQuartetLikelihoods' own "not
    // initialized" sentinel (see QuartetGroups' comment in phylotree.h):
    // it fills in a single group containing every taxon in the alignment,
    // i.e. plain uniform random sampling across all of them -- exactly
    // what this tool wants, no A/B/C/D clustering
    QuartetGroups groups;
    groups.numGroups = 0;

    // params.lmap_num_quartets, not a local variable, is what
    // computeQuartetLikelihoods itself reads to decide how many quartets
    // to draw -- and, via the `params.lmap_num_quartets == uniqueQuarts`
    // check inside it, whether to enumerate every unique quartet
    // (only possible if this equals the alignment's own C(n,4) exactly)
    // or draw this many at random with replacement across quartets
    // (the normal case for any real-sized alignment, and always the
    // requested behavior here)
    params.lmap_num_quartets = numQuartets;

    // computeQuartetLikelihoods parallelizes its own quartet loop via a
    // bare "#pragma omp parallel" (tree/quartet.cpp), governed by
    // OpenMP's own thread count -- entirely separate from
    // tree.setNumThreads() above, which only affects the shared tree's
    // OWN likelihood-kernel threading, not this. Left uncontrolled (an
    // OpenMP build's own default is usually one thread per CPU core),
    // this segfaulted every time it was tested: multiple threads each
    // build their own quartet_tree but all read the SAME shared
    // model_factory/getModel()/getRate() from `tree` (set up once, above,
    // specifically so every quartet is scored under one consistent
    // model) concurrently, and something in that sharing path isn't
    // actually safe for concurrent first-use despite the per-quartet
    // PhyloTree/buffers themselves being independent. Pinning this to
    // exactly `threads` (default 1, i.e. no OpenMP parallelism at all)
    // is the only configuration this was ever actually tested safe under
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif

    vector<QuartetInfo> quartetInfo;
    tree.computeQuartetLikelihoods(quartetInfo, groups);

    ofstream out(outPath.c_str());
    if (!out.is_open()) {
        cerr << "error: cannot open output file '" << outPath << "' for writing" << endl;
        delete aln;
        return 2;
    }
    out << "# alignment=" << alignmentFile << " num_quartets=" << numQuartets
        << " model=" << modelSpec << " seed=" << actualSeed << endl;
    out << "taxon1\ttaxon2\ttaxon3\ttaxon4\tbest_topology\tbest_logl\tlogl_12_34\tlogl_13_24\tlogl_14_23" << endl;

    // pairing for logl[k] follows computeQuartetLikelihoods' own qc[] table
    // (tree/quartet.cpp): k=0 -> (seqID[0],seqID[1])|(seqID[2],seqID[3]),
    // k=1 -> (seqID[0],seqID[2])|(seqID[1],seqID[3]),
    // k=2 -> (seqID[0],seqID[3])|(seqID[1],seqID[2])
    static const int pairA[3][2] = {{0, 1}, {0, 2}, {0, 3}};
    static const int pairB[3][2] = {{2, 3}, {1, 3}, {1, 2}};

    int64_t skipped = 0;
    for (size_t qid = 0; qid < quartetInfo.size(); qid++) {
        const QuartetInfo &q = quartetInfo[qid];
        string names[4];
        for (int i = 0; i < 4; i++)
            names[i] = aln->getSeqName(q.seqID[i]);

        // logl[k] == -1.0 for all 3 k is computeQuartetLikelihoods' own
        // "this quartet's sub-alignment had no usable partition" sentinel
        // (see its own "kept_partitions.size() == 0" branch) -- skip it
        // rather than reporting a fabricated winner
        if (q.logl[0] == -1.0 && q.logl[1] == -1.0 && q.logl[2] == -1.0) {
            skipped++;
            continue;
        }

        int best = 0;
        if (q.logl[1] > q.logl[best]) best = 1;
        if (q.logl[2] > q.logl[best]) best = 2;

        out << names[0] << "\t" << names[1] << "\t" << names[2] << "\t" << names[3] << "\t"
            << "(" << names[pairA[best][0]] << "," << names[pairA[best][1]] << ")|("
            << names[pairB[best][0]] << "," << names[pairB[best][1]] << ")\t"
            << fixed << setprecision(6) << q.logl[best] << "\t"
            << q.logl[0] << "\t" << q.logl[1] << "\t" << q.logl[2] << endl;
    }
    out.close();

    if (!quiet) {
        cout << "wrote " << (quartetInfo.size() - skipped) << " quartet(s) to " << outPath;
        if (skipped > 0)
            cout << " (" << skipped << " skipped: no usable alignment partition for that quartet)";
        cout << endl;
    }

    delete aln;
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 2;
    }

    string alignmentFile = argv[1];
    char *end = nullptr;
    long long numQuartets = strtoll(argv[2], &end, 10);
    if (end == argv[2] || *end != '\0') {
        cerr << "error: <num-quartets> must be an integer, got '" << argv[2] << "'" << endl;
        printUsage(argv[0]);
        return 2;
    }

    string outPath = alignmentFile + ".quartets.tsv";
    int seed = -1;
    bool useModelOverride = false;
    string modelOverrideSpec;
    int threads = 1;
    bool quiet = false;

    for (int i = 3; i < argc; i++) {
        string arg = argv[i];
        if (arg == "out" && i + 1 < argc) {
            outPath = argv[++i];
        } else if (arg == "seed" && i + 1 < argc) {
            seed = atoi(argv[++i]);
        } else if (arg == "model" && i + 1 < argc) {
            useModelOverride = true;
            modelOverrideSpec = argv[++i];
        } else if (arg == "threads" && i + 1 < argc) {
            threads = atoi(argv[++i]);
        } else if (arg == "quiet") {
            quiet = true;
        } else {
            cerr << "error: unrecognized option '" << arg << "'" << endl;
            printUsage(argv[0]);
            return 2;
        }
    }

    return runQuartetSearch(alignmentFile, (int64_t) numQuartets, outPath, seed, useModelOverride,
            modelOverrideSpec, threads, quiet);
}
