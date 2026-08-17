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
    echo "ERROR: python3 not found on PATH (needed for misof_subsample.py)." >&2
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

# --- adjust these to your actual experiment ---
TARGET_SITES=10000
MAX_STEPS=10000
COMMON_SPR_FLAGS="fast quiet notree record"   # starttree is checked separately below since it takes its own argument
RADII=(10 1)          # radius=10 pass, then radius=1 pass, on each iqtree-NNI'd partial alignment
PER_RUN_CAP="2h"      # wall-time cap (via `timeout`) on any single iqtree3-NNI or spr_topology_test
                       # call, so one stuck pass can't eat the whole job's walltime allocation;
                       # each partial alignment is only TARGET_SITES bp (unlike the full
                       # ~595k-bp misofproteinalignment.nex), so this is far smaller than the
                       # 9h cap the old fixed-alignment version of this script used
PARTIAL_DIR="misof_partials"          # every random partial alignment is kept here, never deleted
IQTREE_RUN_DIR="misof_iqtree_runs"    # every iqtree3 NNI run's own output files, same story
NNI_LOG="misof_iqtree_nni_log.csv"    # one row per iqtree3 NNI run: run_id,partial_alignment,tree_file,logL,timestamp
# -----------------------------------------------
mkdir -p "$PARTIAL_DIR" "$IQTREE_RUN_DIR"
[ -s "$NNI_LOG" ] || echo "run_id,partial_alignment,tree_file,logL,timestamp" >> "$NNI_LOG"

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

# Runs iqtree3's own NNI search on a partial alignment, records its final
# logL to $NNI_LOG, and (on success) echoes the path to the resulting
# .treefile on stdout for the caller to feed into spr_topology_test.
# Returns non-zero (nothing echoed) if the run failed or its logL/treefile
# couldn't be found -- the caller is expected to skip this iteration's spr
# passes rather than abort the whole job over one bad run.
run_iqtree_nni() {
    local partial_nex="$1"
    local run_id="$2"
    local prefix="${IQTREE_RUN_DIR}/${run_id}"

    if ! timeout "$PER_RUN_CAP" ./build/iqtree3 -s "$partial_nex" -m LG -nstop 100 -T 1 \
            --prefix "$prefix" --redo >&2; then
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
}

# Repeat, until SLURM's own --time kills the job: draw a fresh random
# TARGET_SITES-bp partial alignment from misofproteinalignment.nex, run
# iqtree3's own NNI search on it, then feed that NNI run's final tree
# (topology + branch lengths + model already fit) into spr_topology_test
# twice -- once at radius=10, once at radius=1 -- each starting from that
# SAME tree rather than rebuilding its own BioNJ estimate. No explicit stop
# condition here by design, same as the rest of this codebase's SLURM
# scripts (see e.g. test_scripts/record_iqtree_nni.sh's own comments).
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

    for radius in "${RADII[@]}"; do
        run_spr "$radius" "$PARTIAL_NEX" "$TREEFILE"
    done
done
