#!/bin/bash
# Run from the directory this job was submitted from, so record/output
# files land somewhere predictable regardless of SLURM's own
# working-directory defaults
cd "$SLURM_SUBMIT_DIR" || exit 1

# Sanity check: fail loudly and immediately if the binary isn't there,
# instead of the job burning its walltime allocation on a "command not
# found" error buried in the log.
if [ ! -x build/spr_topology_test ]; then
    echo "ERROR: build/spr_topology_test not found or not executable." >&2
    echo "Build it first (interactively, not as part of this job) with:" >&2
    echo "  cmake -S . -B build -DBOOST_ROOT=... -DEIGEN3_INCLUDE_DIR=..." >&2
    echo "  cmake --build build --target spr_topology_test" >&2
    exit 1
fi

# --- adjust these to your actual experiment ---
MAX_STEPS=10000
COMMON_SPR_FLAGS="fast quiet notree record"   # starttree is checked separately below since it takes its own argument
RADII=(10 1)          # radius=10 pass, then radius=1 pass, on each already NNI-fit partial alignment
PER_RUN_CAP="2h"      # wall-time cap (via `timeout`) on any single spr_topology_test call, so one
                       # stuck pass can't eat the whole job's walltime allocation; each partial
                       # alignment is far smaller than the full ~595k-bp misofproteinalignment.nex,
                       # so this is comfortably above what a single pass over one of them should need
PARTIAL_DIR="misof_partials"          # partial alignments -- NOT written here; this script only
                                       # consumes whatever the data-generating version of this script
                                       # (see git history) already left behind
IQTREE_RUN_DIR="misof_iqtree_runs"    # matching iqtree3 NNI output (incl. .treefile) for each partial
                                       # above, same naming convention: ${IQTREE_RUN_DIR}/<run_id>.treefile
SPR_DONE_LOG="misof_spr_done_log.csv" # one run_id per line, appended after that run_id's spr passes
                                       # (every radius in RADII) finish successfully -- lets a job
                                       # resubmitted after hitting SLURM's own --time limit pick up
                                       # where the last one left off instead of re-recording (and
                                       # duplicating, in record's own CSV) work already done
# -----------------------------------------------
touch "$SPR_DONE_LOG"

if [ ! -d "$PARTIAL_DIR" ]; then
    echo "ERROR: $PARTIAL_DIR not found -- run the data-generating version of this script" >&2
    echo "(see git history for slurm_misof_hillclimb.sh) first to create partial alignments" >&2
    echo "and their iqtree3 NNI trees for this script to consume." >&2
    exit 1
fi

# Sanity check #2: confirm THIS binary's own --hillclimb actually recognizes
# every flag word this script is about to pass it, before ever entering the
# loop below -- see the header comment history in git blame for the incident
# this guards against (a stale binary's parseHillClimbFlags rejecting the
# whole flag list falls through to printUsage()+return 2 near-instantly, so
# an unguarded loop would spin doing nothing but writing usage text to the
# log for the entire walltime allocation).
usage_output="$(./build/spr_topology_test --hillclimb 2>&1)"
for flag in $COMMON_SPR_FLAGS starttree; do
    if ! grep -q -- "$flag" <<< "$usage_output"; then
        echo "ERROR: build/spr_topology_test's own --hillclimb usage text doesn't mention" >&2
        echo "'$flag' -- this binary is almost certainly built from source that predates" >&2
        echo "that flag. Rebuild before resubmitting: cmake --build build --target spr_topology_test" >&2
        exit 1
    fi
done

# One spr_topology_test --hillclimb pass, starting from the already-fit
# iqtree3 NNI tree/branch-lengths/model (via 'starttree', so no BioNJ/random/
# iqtreestart start-tree construction happens) directly against the same
# partial alignment (via 'notree', so <partial_nex> IS the alignment argument,
# not an AliSim .treefile). 'record' appends this run's trajectory to
# record_LG_fast_notree.csv (or similar, see recordSpreadsheetPath in
# tree/spr_topology_test.cpp) as it goes, so a pass that hits the
# ${PER_RUN_CAP} cap still leaves its partial progress recorded.
run_spr() {
    local radius="$1"
    local partial_nex="$2"
    local start_treefile="$3"
    echo "  === spr radius=$radius -- $(date) ==="
    timeout "$PER_RUN_CAP" ./build/spr_topology_test --hillclimb "$partial_nex" "$radius" "$MAX_STEPS" \
            $COMMON_SPR_FLAGS starttree "$start_treefile"
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

# Walk every partial alignment already sitting in $PARTIAL_DIR (sorted, so
# repeated submissions of this same job process them in the same order),
# one at a time, running the full RADII sweep of spr passes -- via
# run_spr's own 'starttree' -- against its matching already-NNI-fit tree.
# "Matching" here means exactly the naming convention the data-generating
# version of this script uses: a partial alignment
# ${PARTIAL_DIR}/<run_id>.nex pairs with the NNI tree
# ${IQTREE_RUN_DIR}/<run_id>.treefile (run_iqtree_nni's own $prefix there).
# A partial with no such treefile -- or an empty one -- has no valid
# matching end topology to start from, so it's skipped (logged, not
# fatal) rather than aborting the whole job over one bad/incomplete pair.
# Ends -- no "while true" here -- the moment every partial alignment
# currently in $PARTIAL_DIR has either been processed or skipped as
# invalid; re-running new data through this script means re-submitting
# the job after the data-generating version has produced more of it.
mapfile -t PARTIALS < <(find "$PARTIAL_DIR" -maxdepth 1 -name '*.nex' | sort)
if [ "${#PARTIALS[@]}" -eq 0 ]; then
    echo "ERROR: no partial alignments (*.nex) found in $PARTIAL_DIR" >&2
    exit 1
fi

ITER=0
PROCESSED=0
for PARTIAL_NEX in "${PARTIALS[@]}"; do
    ITER=$((ITER + 1))
    RUN_ID="$(basename "$PARTIAL_NEX" .nex)"
    TREEFILE="${IQTREE_RUN_DIR}/${RUN_ID}.treefile"

    echo "=========================================="
    echo "=== entry $ITER/${#PARTIALS[@]} (run_id=$RUN_ID) -- $(date) ==="
    echo "=========================================="

    if grep -qxF "$RUN_ID" "$SPR_DONE_LOG"; then
        echo "  -- $RUN_ID already recorded on an earlier submission of this job; skipping"
        continue
    fi
    if [ ! -s "$TREEFILE" ]; then
        echo "  -- $RUN_ID has no valid matching end topology ($TREEFILE missing or empty); skipping" >&2
        continue
    fi

    for radius in "${RADII[@]}"; do
        run_spr "$radius" "$PARTIAL_NEX" "$TREEFILE"
    done
    echo "$RUN_ID" >> "$SPR_DONE_LOG"
    PROCESSED=$((PROCESSED + 1))
done

echo "=========================================="
echo "=== done: $PROCESSED/${#PARTIALS[@]} partial alignment(s) processed, none left -- $(date) ==="
echo "=========================================="
