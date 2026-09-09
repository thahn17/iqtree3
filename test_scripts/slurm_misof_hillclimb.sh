#!/bin/bash
# Run from the directory this job was submitted from, so record/output
# files land somewhere predictable regardless of SLURM's own
# working-directory defaults
cd "$SLURM_SUBMIT_DIR" || exit 1

# Sanity check: fail loudly and immediately if a binary isn't there, instead
# of the job burning its walltime allocation on a "command not found" error
# buried in the log.
if [ ! -x build/spr_topology_test ]; then
    echo "ERROR: build/spr_topology_test not found or not executable." >&2
    echo "Build it first (interactively, not as part of this job) with:" >&2
    echo "  cmake -S . -B build -DBOOST_ROOT=... -DEIGEN3_INCLUDE_DIR=..." >&2
    echo "  cmake --build build --target spr_topology_test" >&2
    exit 1
fi
if [ ! -x build/iqtree3 ]; then
    echo "ERROR: build/iqtree3 not found or not executable." >&2
    echo "Build it first with: cmake --build build --target iqtree3" >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 not found on PATH (needed for misof_subsample.py and" >&2
    echo "extract_iqtree_model.py)." >&2
    exit 1
fi

ALIGNMENT="misofproteinalignment.nex"
if [ ! -f "$ALIGNMENT" ]; then
    echo "ERROR: $ALIGNMENT not found in $SLURM_SUBMIT_DIR." >&2
    exit 1
fi
SUBSAMPLE_SCRIPT="test_scripts/misof_subsample.py"
if [ ! -f "$SUBSAMPLE_SCRIPT" ]; then
    echo "ERROR: $SUBSAMPLE_SCRIPT not found." >&2
    exit 1
fi
EXTRACT_MODEL_SCRIPT="test_scripts/extract_iqtree_model.py"
if [ ! -f "$EXTRACT_MODEL_SCRIPT" ]; then
    echo "ERROR: $EXTRACT_MODEL_SCRIPT not found." >&2
    exit 1
fi

# --- adjust these to your actual experiment ---
TARGET_SITES=10000
MAX_STEPS=10000
NNI_MODEL="MFP"        # ModelFinder Plus: let iqtree3 pick the best-fit protein model itself
                       # (any AA matrix x frequency-handling x rate-heterogeneity combination --
                       # LG/WAG/JTT/Q.insect/..., +F/+FO/+FQ, +I/+G/+R -- scored by BIC, not a
                       # single model fixed up front) for every partial's own NNI run below, then
                       # continue straight into that same run's tree search under whichever model
                       # won. extract_iqtree_model.py (see its own comment) reads back WHATEVER
                       # ModelFinder actually selected -- base matrix, fitted frequencies, and
                       # fitted +I/+G/+R rate-heterogeneity parameters, in any combination -- and
                       # that fitted model (not just spr_topology_test's own un-fit LG default,
                       # which has neither rate heterogeneity nor fitted frequencies) is what gets
                       # carried into the SPR search afterward, so its own topology comparisons
                       # are scored under the same model ModelFinder judged best for this data
                       # rather than a model picked without looking at it. Confirmed in practice
                       # to badly overrun even a 3h PER_RUN_CAP under -T 1 on a 144-taxon/10000-site
                       # partial (ModelFinder's own per-candidate-model branch-length fits are the
                       # cost, not the later tree search) -- NNI_MSET and NNI_THREADS below exist
                       # specifically to bring that back down to something that actually finishes
NNI_MSET="LG,WAG,JTT,Q.insect,Q.pfam" # restrict ModelFinder to a handful of protein matrices
                       # known to fit real (esp. insect phylogenomic) data well, instead of its
                       # full default candidate set. On its own this was NOT enough: a real test
                       # run still printed "ModelFinder will test up to 220 protein models" with
                       # --mset alone (5 matrices x ~44 rate/frequency combinations each -- --mset
                       # only restricts the matrix axis, not the rate-heterogeneity one) -- see
                       # NNI_MRATE below for the other, equally necessary axis
NNI_MRATE="G,I+G"      # restrict ModelFinder's rate-heterogeneity search to plain Gamma and
                       # Invariant+Gamma, instead of its full default sweep (E, I, G, I+G, and
                       # FreeRate R2 through R10 -- nine separate FreeRate fits alone, the single
                       # most expensive part of the 220-model default). Combined with NNI_MSET,
                       # this cuts ModelFinder down to 5 matrices x 2 rate models = 10 total fits
                       # instead of 220 -- still real, data-driven selection between two of the
                       # most commonly-favored rate treatments, just not an exhaustive one. Widen
                       # this (e.g. add ",R" for FreeRate) once a run is confirmed to finish
                       # comfortably inside PER_RUN_CAP
NNI_THREADS="${SLURM_CPUS_PER_TASK:-4}" # ModelFinder parallelizes across candidate models far
                       # more effectively than a single NNI tree search does (the previous fixed
                       # "-T 1" was sized for that single-model NNI case, not for MFP testing many
                       # models), so this is the other main lever on wall-clock time. Picks up
                       # however many CPUs this SLURM job actually requested
                       # (--cpus-per-task); falls back to 4 outside SLURM (e.g. local testing).
                       # Deliberately a concrete count rather than "-T AUTO": AUTO makes iqtree3
                       # run its own thread-count benchmarking pass first, which on a
                       # near-zero-pattern-compression alignment has been seen (elsewhere in this
                       # project, on the full un-subsampled alignment) to attempt a many-GB
                       # allocation baseline it doesn't need for these much smaller partials
COMMON_SPR_FLAGS="fast quiet notree record investigate 10 fullreopt 1 1"   # starttree/model are appended separately below, since they take their own arguments
RADII=(10)             # one radius=10 SPR pass per newly-sampled partial alignment. A radius=1
                       # pass used to also run here (against the same starting tree, not chained
                       # after radius=10); dropped after record_LG_fast_notree.csv showed it never
                       # once found an improving move across all 7 partials it was tried on (every
                       # radius=1 pass ended at exactly its own starting logL), while radius=10
                       # found one on 1 of those same 7 -- expected, since each starting tree is
                       # already a genuine iqtree3 NNI-search optimum, so NNI-equivalent radius=1
                       # moves from it have little left to find. Halves this loop's own wall-clock
                       # cost per partial (each radius gets its own PER_RUN_CAP budget) for a pass
                       # that had a 0-for-7 track record
PER_RUN_CAP="2h"      # wall-time cap (via `timeout`) on any single iqtree3-NNI or spr_topology_test
                       # call, so one stuck pass can't eat the whole job's walltime allocation;
                       # each partial alignment is only TARGET_SITES bp (unlike the full
                       # ~595k-bp misofproteinalignment.nex), so this is far smaller than a cap
                       # sized for the full-alignment case would need to be
PARTIAL_DIR="misof_partials"          # every random partial alignment is kept here, never deleted
IQTREE_RUN_DIR="misof_iqtree_runs"    # every iqtree3 NNI run's own output files (.treefile,
                                       # .iqtree report, ...), same story -- same naming
                                       # convention as PARTIAL_DIR: ${IQTREE_RUN_DIR}/<run_id>.*
NNI_LOG="misof_iqtree_nni_log.csv"    # one row per iqtree3 NNI run: run_id,partial_alignment,tree_file,logL,timestamp
MODEL_LOG="misof_model_specs_used.csv" # run_id,model_spec -- one row per run_id, written once its
                                       # NNI_MODEL spec is extracted from that run's own .iqtree
                                       # report (the "model tuning" file NNI_MODEL's fit produces),
                                       # independently of spr_topology_test's own record CSV/run_id
                                       # tag (which only marks "gtr" was passed, not the actual
                                       # fitted frequency/alpha values used) -- see
                                       # extract_iqtree_model.py's own comment for why those values
                                       # need carrying over explicitly rather than re-estimating fresh
# -----------------------------------------------
mkdir -p "$PARTIAL_DIR" "$IQTREE_RUN_DIR"
[ -s "$NNI_LOG" ] || echo "run_id,partial_alignment,tree_file,logL,timestamp" >> "$NNI_LOG"
[ -s "$MODEL_LOG" ] || echo "run_id,model_spec" >> "$MODEL_LOG"

# Sanity check #2: confirm THIS binary's own --hillclimb actually recognizes
# every flag word this script is about to pass it, before ever entering the
# loop below -- see the header comment history in git blame for the incident
# this guards against (a stale binary's parseHillClimbFlags rejecting the
# whole flag list falls through to printUsage()+return 2 near-instantly, so
# an unguarded loop would spin doing nothing but writing usage text to the
# log for the entire walltime allocation).
usage_output="$(./build/spr_topology_test --hillclimb 2>&1)"
for flag in $COMMON_SPR_FLAGS starttree model gtr; do
    if ! grep -q -- "$flag" <<< "$usage_output"; then
        echo "ERROR: build/spr_topology_test's own --hillclimb usage text doesn't mention" >&2
        echo "'$flag' -- this binary is almost certainly built from source that predates" >&2
        echo "that flag. Rebuild before resubmitting: cmake --build build --target spr_topology_test" >&2
        exit 1
    fi
done

# Runs iqtree3's own NNI search on a partial alignment under NNI_MODEL (a
# real ModelFinder selection + ML fit: whichever AA matrix and fitted
# frequency/rate-heterogeneity parameters BIC actually favors for this
# data, not spr_topology_test's own fixed plain-LG default -- this
# selection+fit IS the "model tuning" step, its .iqtree report the file
# that tuning produces), records
# its final logL to $NNI_LOG, and (on success) echoes the path to the
# resulting .treefile on stdout for the caller to feed into
# spr_topology_test. Returns non-zero (nothing echoed) if the run failed or
# its logL/treefile couldn't be found -- the caller is expected to skip
# this iteration's spr passes rather than abort the whole job over one bad
# run.
run_iqtree_nni() {
    local partial_nex="$1"
    local run_id="$2"
    local prefix="${IQTREE_RUN_DIR}/${run_id}"

    if ! timeout "$PER_RUN_CAP" ./build/iqtree3 -s "$partial_nex" -m "$NNI_MODEL" --mset "$NNI_MSET" \
            --mrate "$NNI_MRATE" -nstop 100 -T "$NNI_THREADS" --prefix "$prefix" --redo >&2; then
        echo "  -- iqtree3 NNI run failed or hit the ${PER_RUN_CAP} cap (run_id=$run_id); see ${prefix}.log" >&2
        return 1
    fi

    local treefile="${prefix}.treefile"
    if [ ! -f "$treefile" ]; then
        echo "  -- iqtree3 exited 0 but ${treefile} is missing (run_id=$run_id)" >&2
        return 1
    fi

    local logl
    logl="$(grep -m1 "Log-likelihood of the tree:" "${prefix}.iqtree" \
        | sed -E 's/.*Log-likelihood of the tree:[[:space:]]*(-?[0-9.]+).*/\1/')"
    if [ -z "$logl" ]; then
        echo "  -- could not find logL in ${prefix}.iqtree (run_id=$run_id)" >&2
        return 1
    fi

    echo "${run_id},${partial_nex},${treefile},${logl},$(date -Iseconds)" >> "$NNI_LOG"
    echo "  -- iqtree3 NNI logL=${logl} -- ${treefile}" >&2
    echo "$treefile"
}

# Extract whichever model ModelFinder selected under NNI_MODEL -- base
# matrix, fitted frequencies, fitted +I/+G/+R rate heterogeneity, in
# whatever combination BIC picked -- from a run_iqtree_nni run's .iqtree
# report, log the resulting spec to
# $MODEL_LOG, and echo it on stdout for the caller to pass into run_spr.
# Returns non-zero (nothing echoed) on any failure -- same "skip, don't
# abort the job" contract as run_iqtree_nni.
extract_model_spec() {
    local run_id="$1"
    local iqtree_report="${IQTREE_RUN_DIR}/${run_id}.iqtree"

    if [ ! -s "$iqtree_report" ]; then
        echo "  -- $run_id has no .iqtree report ($iqtree_report missing or empty); can't" >&2
        echo "     build a fitted $NNI_MODEL spec" >&2
        return 1
    fi

    local model_spec
    model_spec="$(python3 "$EXTRACT_MODEL_SCRIPT" "$iqtree_report" 2>&1)"
    if [ $? -ne 0 ] || [ -z "$model_spec" ]; then
        echo "  -- $run_id: couldn't extract a model spec from $iqtree_report:" >&2
        echo "     $model_spec" >&2
        return 1
    fi

    echo "  -- model spec: $model_spec" >&2
    # double-quoted: model_spec's own "+F{f1,f2,...}" frequency list is
    # full of internal commas that would otherwise fragment this into many
    # CSV columns instead of staying inside one
    echo "${run_id},\"${model_spec}\"" >> "$MODEL_LOG"
    echo "$model_spec"
}

# One spr_topology_test --hillclimb pass, starting from the already-fit
# iqtree3 NNI tree/branch-lengths (via 'starttree', so no BioNJ/random/
# iqtreestart start-tree construction happens) directly against the same
# partial alignment (via 'notree', so <partial_nex> IS the alignment argument,
# not an AliSim .treefile), scored under that same NNI run's own fitted
# NNI_MODEL (via 'model "<spec>"', built by extract_model_spec above from
# the matching .iqtree report) instead of spr_topology_test's plain-LG
# default. 'gtr' is passed alongside 'model' purely so record/run_id
# filenames self-document that a richer-than-default model was used (see
# 'model <spec>'s own comment in tree/spr_topology_test_usage.txt -- 'gtr'
# alone would NOT get the real fitted values, and 'model' alone would build
# the right model but leave record's own filenames looking like a
# plain-LG run). 'record' appends this run's trajectory to
# record_LG+FO_fast_notree_gtr.csv (or similar, see recordSpreadsheetPath in
# tree/spr_topology_test.cpp) as it goes, so a pass that hits the
# ${PER_RUN_CAP} cap still leaves its partial progress recorded.
run_spr() {
    local radius="$1"
    local partial_nex="$2"
    local start_treefile="$3"
    local model_spec="$4"
    echo "  === spr radius=$radius -- $(date) ==="
    timeout "$PER_RUN_CAP" ./build/spr_topology_test --hillclimb "$partial_nex" "$radius" "$MAX_STEPS" \
            $COMMON_SPR_FLAGS starttree "$start_treefile" model "$model_spec" gtr
    local status=$?
    if [ "$status" -eq 124 ]; then
        echo "  -- radius=$radius hit the ${PER_RUN_CAP} cap (progress up to that point is already recorded)"
    else
        # --hillclimb's own success path returns 2, not 0 (tree/spr_topology_test.cpp,
        # runHillClimb's final "return 2;") -- this is the NORMAL outcome, not a failure
        echo "  -- radius=$radius exited with status $status"
    fi
    return 0
}

# Repeat, until SLURM's own --time kills the job: draw a fresh random
# TARGET_SITES-bp partial alignment from misofproteinalignment.nex, run
# iqtree3's own NNI search on it under NNI_MODEL (the "model tuning" step
# -- ModelFinder selects whichever model actually fits this data best,
# then that model's own frequency/rate-heterogeneity parameters get
# ML-fit, not spr_topology_test's own fixed, un-fit plain-LG default),
# then feed that NNI run's final tree (topology + branch lengths already
# fit) AND its ModelFinder-selected, fitted model into spr_topology_test
# for one radius=10 SPR pass (see RADII's own comment for why just the
# one pass/radius) starting from that SAME tree/model rather than
# rebuilding its own BioNJ estimate or falling back to plain LG. No
# explicit stop condition here by design, same as the rest of
# this codebase's SLURM scripts (see e.g. test_scripts/record_iqtree_nni.sh's
# own comments) -- every iteration generates its own freshly-timestamped
# run_id, so a job resubmitted after hitting SLURM's own --time limit just
# keeps generating new partials rather than needing to resume a specific
# one; nothing here is ever re-processed, so (unlike the old fixed-backlog
# consumer version of this script) there's no done-log to check against.
ITER=0
while true; do
    ITER=$((ITER + 1))
    RUN_ID="${SLURM_JOB_ID:-local}_$(date +%Y%m%d-%H%M%S)_it$(printf '%04d' "$ITER")"
    PARTIAL_NEX="${PARTIAL_DIR}/${RUN_ID}.nex"

    echo "=========================================="
    echo "=== iteration $ITER (run_id=$RUN_ID) -- $(date) ==="
    echo "=========================================="

    echo "  -- sampling ${TARGET_SITES} random sites from ${ALIGNMENT} -> ${PARTIAL_NEX}"
    if ! python3 "$SUBSAMPLE_SCRIPT" "$ALIGNMENT" "$TARGET_SITES" "$PARTIAL_NEX"; then
        echo "  -- subsampling failed (run_id=$RUN_ID); skipping this iteration" >&2
        continue
    fi

    TREEFILE="$(run_iqtree_nni "$PARTIAL_NEX" "$RUN_ID")"
    if [ -z "$TREEFILE" ]; then
        echo "  -- no usable iqtree3 NNI tree this iteration; skipping the spr passes" >&2
        continue
    fi

    MODEL_SPEC="$(extract_model_spec "$RUN_ID")"
    if [ -z "$MODEL_SPEC" ]; then
        echo "  -- no usable fitted model this iteration; skipping the spr passes" >&2
        continue
    fi

    for radius in "${RADII[@]}"; do
        run_spr "$radius" "$PARTIAL_NEX" "$TREEFILE" "$MODEL_SPEC"
    done
done
