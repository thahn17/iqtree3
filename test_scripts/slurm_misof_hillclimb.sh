#!/bin/bash
# Run from the directory this job was submitted from, so record_*.csv and
# other outputs land somewhere predictable regardless of SLURM's own
# working-directory defaults
cd "$SLURM_SUBMIT_DIR" || exit 1

# Sanity check: fail loudly and immediately if the binary isn't there,
# instead of the job burning its walltime allocation on a "command not
# found" error buried in the log. Unlike a replicate sweep over AliSim-
# simulated data, this job runs a single fixed real alignment
# (misofproteinalignment.nex) with no ground-truth tree via --hillclimb's
# "notree" mode, so there's no per-pass iqtree3/AliSim regeneration step --
# only spr_topology_test is needed.
if [ ! -x build/spr_topology_test ]; then
    echo "ERROR: build/spr_topology_test not found or not executable." >&2
    echo "Build it first (interactively, not as part of this job) with:" >&2
    echo "  cmake -S . -B build -DBOOST_ROOT=... -DEIGEN3_INCLUDE_DIR=..." >&2
    echo "  cmake --build build --target spr_topology_test" >&2
    exit 1
fi

ALIGNMENT="misofproteinalignment.nex"
if [ ! -f "$ALIGNMENT" ]; then
    echo "ERROR: $ALIGNMENT not found in $SLURM_SUBMIT_DIR." >&2
    exit 1
fi

# --- adjust these to your actual experiment ---
MAX_STEPS=10000
COMMON_FLAGS="fast quiet notree findopt 20 record"
# Sanity check #2: confirm THIS binary's own --hillclimb actually recognizes
# every flag word in COMMON_FLAGS, before ever entering the "while true"
# loop below. A binary built from source that predates one of these flags
# (e.g. "notree"/"findopt") doesn't error out per-call -- parseHillClimbFlags
# just rejects the whole trailing flag list, main() falls through to
# printUsage()+return 2, and since that happens near-instantly (no alignment
# ever gets touched), the while-loop below would otherwise spin on it as
# fast as the OS allows for the entire walltime allocation, writing nothing
# but repeated usage text to the log and never producing a record_*.csv --
# exactly what happened the first time this ran against a stale binary.
# `--hillclimb` with no further args always falls through to that same
# printUsage() (argc is too small to match the real --hillclimb dispatch),
# so its output is a reliable, near-instant proxy for "does this binary's
# own help text even mention the flags I'm about to pass it" -- no alignment
# load, no BioNJ build, nothing slow.
usage_output="$(./build/spr_topology_test --hillclimb 2>&1)"
for flag in $COMMON_FLAGS; do
    case "$flag" in
        ''|*[!0-9]*) ;;   # not purely numeric (a real flag word) -- check it
        *) continue ;;    # purely numeric (e.g. findopt's "20") -- skip
    esac
    if ! grep -q -- "$flag" <<< "$usage_output"; then
        echo "ERROR: build/spr_topology_test's own --hillclimb usage text doesn't" >&2
        echo "mention '$flag' (COMMON_FLAGS=\"$COMMON_FLAGS\") -- this binary is almost" >&2
        echo "certainly built from source that predates that flag. Copy the CURRENT" >&2
        echo "tree/spr_topology_test.cpp here and rebuild before resubmitting:" >&2
        echo "  cmake --build build --target spr_topology_test" >&2
        exit 1
    fi
done
PER_RUN_CAP="9h"   # each individual --hillclimb call is capped at this wall
                   # time via `timeout`; on a 595k-character/144-taxon real
                   # alignment there's no guarantee 10000 steps finishes in
                   # any bounded time, so this keeps any one pass from
                   # eating the whole job's walltime allocation. "record"
                   # appends its trajectory as it goes (see appendRecordRow
                   # in tree/spr_topology_test.cpp), so a pass that gets cut
                   # off here still leaves its partial progress recorded.
RADII=(10 1)       # radius=10 pass, then radius=1 pass, every loop iteration
# -----------------------------------------------

run_one() {
    local radius="$1"
    echo "=========================================="
    echo "=== radius=$radius -- $(date) ==="
    echo "=========================================="
    timeout "$PER_RUN_CAP" ./build/spr_topology_test --hillclimb "$ALIGNMENT" "$radius" "$MAX_STEPS" $COMMON_FLAGS
    local status=$?
    if [ "$status" -eq 124 ]; then
        echo "  -- radius=$radius hit the ${PER_RUN_CAP} cap and was cut off (progress up to that point is already in record_*.csv)"
    else
        # NOTE: --hillclimb's own success path returns 2, not 0 (see
        # runHillClimb's final "return 2;" in tree/spr_topology_test.cpp) --
        # this predates and is unrelated to the "notree"/protein work, it's
        # just this tool's own long-standing convention. So status==2 here
        # is the NORMAL outcome, not a failure -- only used for the log.
        echo "  -- radius=$radius exited with status $status"
    fi
}

# Keep alternating radius=10 and radius=1 passes, each individually capped
# at $PER_RUN_CAP, every pass appending to the same record_LG_fast.csv
# (recordProgress/"record" always appends, never truncates). There is no
# explicit stop condition here by design -- if both passes finish (or get
# cut off) within their caps, the loop just starts another radius=10/
# radius=1 pair, continuing to add both types to the recording until
# SLURM's own --time above is what finally kills the job.
while true; do
    for radius in "${RADII[@]}"; do
        run_one "$radius"
    done
done
