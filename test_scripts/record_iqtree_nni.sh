#!/bin/bash
# Records IQ-TREE's own built-in stochastic NNI tree search (the real
# pipeline -- main/phyloanalysis.cpp + tree/iqtree.cpp -- run via the full
# iqtree3.exe, NOT spr_topology_test's standalone SPR API) into a CSV with
# the same columns spr_topology_test's --hillclimb "record" flag writes:
# run_id,candidates,time_elapsed,logL,true_minus_current
# (see recordSpreadsheetPath/appendRecordRow in tree/spr_topology_test.cpp).
#
# Bash port of record_iqtree_nni.ps1 (same logic) -- use this one instead of
# the .ps1 version if PowerShell script execution is disabled on your system
# (e.g. "running scripts is disabled on this system", an ExecutionPolicy
# restriction); this runs fine under Git Bash / MSYS2 bash with no such
# restriction, and needs no interpreter flag or execution policy at all.
#
# Two things distinguish this from just eyeballing an `iqtree3 -s ... -m ...`
# log by hand:
#
# 1) "Own start": the search NEVER receives the AliSim ground-truth tree as
#    a starting tree or constraint (no -t/-g) -- it builds its own
#    BIONJ/parsimony start, exactly like a real analysis would. The true
#    tree is only ever used in a SEPARATE `-te` (fixed-tree) evaluation run,
#    purely to compute the true_minus_current reference column, the same
#    "known answer, only used to score the attempt" role sim.treefile plays
#    in spr_topology_test's own --hillclimb.
#
# 2) Full wall-clock timing. IQ-TREE prints "Iteration N / LogL: ... / Time:
#    HhMmSs" per iteration and "Total wall-clock time used: ..." at the end,
#    but BOTH are measured from params.start_real_time, which
#    phyloanalysis.cpp only sets AFTER alignment reading, PLL setup, and the
#    initial BIONJ/parsimony distance computation have already run (see the
#    runPhyloAnalysis call order around "params.start_real_time =
#    getRealTime();") -- i.e. IQ-TREE's own clock starts only once a start
#    tree already exists. This script instead times the whole external
#    process (spawn to exit) and adds the (externalTotal - internalTotal)
#    difference as a constant offset to every per-iteration timestamp, so
#    time_elapsed in the CSV really does span from process launch to
#    completion.
#
# Usage (mirrors spr_topology_test --hillclimb's calling convention):
#     test_scripts/record_iqtree_nni.sh <alisim-tree.treefile> <nstop> [gtr] [quiet]
#
# <alisim-tree.treefile>  path to the AliSim ground-truth tree; the
#                          alignment is derived by replacing ".treefile"
#                          with ".fa", same convention as --hillclimb.
# <nstop>                 passed straight through as iqtree3's own -nstop
#                          (number of unsuccessful iterations before the
#                          stochastic NNI search gives up) -- the closest
#                          built-in-NNI analogue of --hillclimb's <max-steps>.
# gtr                     search under GTR+FO (ML-estimated rates/freqs)
#                          instead of plain JC.
# quiet                   only print the run summary, not one line per
#                          recorded iteration (does not affect what's
#                          written to the CSV).
#
# Env overrides (all optional): IQTREE_BIN (default build/iqtree3.exe),
# WORK_DIR (default ${TMPDIR:-/tmp}/iqtree_nni_runs), CSV_PATH_OVERRIDE
# (default record_<model>_iqtree_nni_ownstart.csv).
#
# Examples:
#     test_scripts/record_iqtree_nni.sh sim.treefile 100
#     test_scripts/record_iqtree_nni.sh sim.treefile 100 gtr quiet

set -euo pipefail

TREE_FILE="${1:-}"
MAX_STEPS="${2:-}"
if [ -z "$TREE_FILE" ] || [ -z "$MAX_STEPS" ]; then
    echo "Usage: $0 <alisim-tree.treefile> <nstop> [gtr] [quiet]" >&2
    exit 1
fi
shift 2

USE_GTR=0
QUIET=0
for arg in "$@"; do
    case "$arg" in
        gtr) USE_GTR=1 ;;
        quiet) QUIET=1 ;;
        *) echo "Unknown flag: $arg" >&2; exit 1 ;;
    esac
done

IQTREE_BIN="${IQTREE_BIN:-build/iqtree3.exe}"
WORK_DIR="${WORK_DIR:-${TMPDIR:-/tmp}/iqtree_nni_runs}"
CSV_PATH_OVERRIDE="${CSV_PATH_OVERRIDE:-}"

[ -f "$TREE_FILE" ] || { echo "Tree file not found: $TREE_FILE" >&2; exit 1; }
[ -f "$IQTREE_BIN" ] || { echo "iqtree3 binary not found at $IQTREE_BIN -- build it first: cmake --build build --target iqtree3" >&2; exit 1; }

ALIGNMENT="${TREE_FILE%.treefile}.fa"
if [ "$ALIGNMENT" = "$TREE_FILE" ]; then
    echo "Expected a .treefile path (alignment is derived by replacing .treefile with .fa): $TREE_FILE" >&2
    exit 1
fi
[ -f "$ALIGNMENT" ] || { echo "Alignment not found next to tree file (expected $ALIGNMENT)" >&2; exit 1; }

if [ "$USE_GTR" = 1 ]; then MODEL="GTR+FO"; else MODEL="JC"; fi

mkdir -p "$WORK_DIR"

TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
RUN_ID="${TIMESTAMP}_iqtree_nni_ownstart_nstop${MAX_STEPS}"
PREFIX="${WORK_DIR}/${RUN_ID}"

# iqtree3 tees everything it prints into its own <prefix>.log file
# regardless of what its stdout is connected to, so unless quiet was
# requested we just let stdout inherit the terminal directly -- that's what
# gives live progress instead of the whole run going silent until it exits.
# Post-run parsing below always reads the authoritative <prefix>.log file,
# not this stdout stream, so this only affects what you see live.
run_iqtree() {
    if [ "$QUIET" = 1 ]; then
        "$IQTREE_BIN" "$@" > /dev/null 2>&1
    else
        "$IQTREE_BIN" "$@"
    fi
}

# --- Step 1: evaluate the true AliSim tree's own logL under the same model.
# This run is NEVER fed into the search below; it exists only to compute
# the true_minus_current reference column.
[ "$QUIET" = 1 ] || echo "Evaluating true tree's likelihood (reference only -- not used by the search)..."

TRUEEVAL_PREFIX="${PREFIX}_trueeval"
if ! run_iqtree -s "$ALIGNMENT" -te "$TREE_FILE" -m "$MODEL" --show-lh -T 1 \
        --prefix "$TRUEEVAL_PREFIX" --redo; then
    echo "trueeval run failed; see ${TRUEEVAL_PREFIX}.log" >&2
    exit 1
fi

TRUE_LOGL="$(grep -m1 "Log-likelihood of the tree:" "${TRUEEVAL_PREFIX}.iqtree" \
    | sed -E 's/.*Log-likelihood of the tree:[[:space:]]*(-?[0-9.]+).*/\1/')"
if [ -z "$TRUE_LOGL" ]; then
    echo "Could not find true tree log-likelihood in ${TRUEEVAL_PREFIX}.iqtree" >&2
    exit 1
fi
[ "$QUIET" = 1 ] || echo "True tree logL: $TRUE_LOGL"

# --- Step 2: the actual search, own BIONJ/parsimony start (no -t/-g), timed
# externally from just before the process is spawned to just after it exits.
[ "$QUIET" = 1 ] || echo "Running iqtree3's own NNI search (own start, model=$MODEL, nstop=$MAX_STEPS)..."

MAIN_PREFIX="${PREFIX}_run"
EXTERNAL_START="$(date +%s)"
if ! run_iqtree -s "$ALIGNMENT" -m "$MODEL" -v -T 1 --prefix "$MAIN_PREFIX" --redo -nstop "$MAX_STEPS"; then
    echo "NNI search run failed; see ${MAIN_PREFIX}.log" >&2
    exit 1
fi
EXTERNAL_END="$(date +%s)"
EXTERNAL_TOTAL=$((EXTERNAL_END - EXTERNAL_START))

# --- Step 3: parse the log for per-iteration progress and IQ-TREE's own
# (partial) total, then correct every timestamp by the difference.
LOG_FILE="${MAIN_PREFIX}.log"

INTERNAL_TOTAL="$(grep -o "Total wall-clock time used: [0-9.]*" "$LOG_FILE" | tail -1 | grep -o "[0-9.]*$")"
if [ -z "$INTERNAL_TOTAL" ]; then
    echo "Could not find 'Total wall-clock time used' in $LOG_FILE" >&2
    exit 1
fi

OFFSET="$(awk -v e="$EXTERNAL_TOTAL" -v i="$INTERNAL_TOTAL" 'BEGIN{printf "%.3f", e-i}')"
[ "$QUIET" = 1 ] || echo "External process time: ${EXTERNAL_TOTAL}s -- IQ-TREE-reported time: ${INTERNAL_TOTAL}s -- pre-search-timer overhead added back: ${OFFSET}s"

# --- Step 4: append to the record spreadsheet, same format as
# recordSpreadsheetPath/appendRecordRow in tree/spr_topology_test.cpp.
SANITIZED_MODEL="$(echo "$MODEL" | sed -E 's/[^A-Za-z0-9]/_/g')"
if [ -n "$CSV_PATH_OVERRIDE" ]; then
    CSV_PATH="$CSV_PATH_OVERRIDE"
else
    CSV_PATH="record_${SANITIZED_MODEL}_iqtree_nni_ownstart.csv"
fi

if [ ! -s "$CSV_PATH" ]; then
    echo "run_id,candidates,time_elapsed,logL,true_minus_current" >> "$CSV_PATH"
fi

ROW_COUNT=0
while IFS= read -r line; do
    if [[ "$line" =~ ^(Iteration|Bootstrap)\ ([0-9]+)\ /\ LogL:\ (-?[0-9.]+)\ /\ Time:\ ([0-9]+)h:([0-9]+)m:([0-9]+)s ]]; then
        CANDIDATES="${BASH_REMATCH[2]}"
        LOGL="${BASH_REMATCH[3]}"
        H="${BASH_REMATCH[4]}"; M="${BASH_REMATCH[5]}"; S="${BASH_REMATCH[6]}"
        INTERNAL_SEC=$((10#$H * 3600 + 10#$M * 60 + 10#$S))
        TIME_ELAPSED="$(awk -v s="$INTERNAL_SEC" -v o="$OFFSET" 'BEGIN{printf "%.3f", s+o}')"
        TRUE_MINUS_CURRENT="$(awk -v t="$TRUE_LOGL" -v l="$LOGL" 'BEGIN{printf "%.6f", t-l}')"
        echo "${RUN_ID},${CANDIDATES},${TIME_ELAPSED},${LOGL},${TRUE_MINUS_CURRENT}" >> "$CSV_PATH"
        ROW_COUNT=$((ROW_COUNT + 1))
        [ "$QUIET" = 1 ] || echo "  candidates=${CANDIDATES} time_elapsed=${TIME_ELAPSED}s logL=${LOGL} true_minus_current=${TRUE_MINUS_CURRENT}"
    fi
done < "$LOG_FILE"

if [ "$ROW_COUNT" -eq 0 ]; then
    echo "No 'Iteration N / LogL: ... / Time: ...' lines found in $LOG_FILE -- was -v honored?" >&2
    exit 1
fi

echo "Appended ${ROW_COUNT} rows to ${CSV_PATH} (run_id=${RUN_ID})"
echo "Raw iqtree3 outputs kept under ${WORK_DIR} (prefix ${RUN_ID}*) for inspection."
