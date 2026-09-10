#!/bin/bash
#
# Run test_scripts/run_acceptdist_bootstrap.sh repeatedly until the job's
# walltime is nearly gone, accumulating independent bootstrap trials.
#
# No #SBATCH header on purpose -- supply the resource request on the command
# line (sbatch -t 24:00:00 -c 8 --mem=8G test_scripts/slurm_acceptdist_bootstrap.sh)
# or prepend your own directives. One header line IS worth adding if your site
# allows it, because it enables the graceful-save path below:
#
#     #SBATCH --signal=B:USR1@600
#
# which makes SLURM send SIGUSR1 ten minutes before the kill, so the job can
# finish its bookkeeping instead of dying mid-write.
#
# WHY A LOOP. Each invocation of the inner script draws a NEW random BOOTSEED,
# so it generates a fresh bootstrap alignment and a fresh set of paired blocks.
# Running it repeatedly accumulates independent samples rather than repeating
# the same ones, and the inner script's pooled summary spans every trial it
# finds. A long job therefore produces a progressively stronger comparison
# rather than one big fixed-size result.
#
# EVERYTHING IS WRITTEN STRAIGHT TO A PERSISTENT DIRECTORY under the submit
# directory -- not to node-local scratch with a copy-back at the end. These
# runs are CPU-bound, not I/O-bound, so there is nothing to gain from scratch,
# and a job killed between the last run and the copy-back would lose the whole
# allocation's work. Results are durable the moment each run finishes.
#
set -u

cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1

# ------------------------------------------------------------ sanity checks --
# Fail loudly and immediately rather than burning the allocation on an error
# buried thousands of lines into a log.
INNER="test_scripts/run_acceptdist_bootstrap.sh"
IQ="${IQ:-$PWD/build/iqtree3}"
[ -x "$IQ" ] || IQ="$PWD/build/iqtree3.exe"

if [ ! -f "$INNER" ]; then
    echo "ERROR: $INNER not found (run this from the repo root, or set SLURM_SUBMIT_DIR)." >&2
    exit 1
fi
if [ ! -x "$IQ" ]; then
    echo "ERROR: iqtree3 not found or not executable at $IQ" >&2
    echo "Build it first (interactively, not as part of this job) with:" >&2
    echo "  cmake --build build --target iqtree3" >&2
    exit 1
fi
if ! command -v python3 >/dev/null 2>&1 && ! command -v python >/dev/null 2>&1; then
    echo "ERROR: no python on PATH (the inner script needs it to generate" >&2
    echo "bootstrap alignments and to build the summary tables)." >&2
    exit 1
fi

SRCALN="${SRCALN:-$PWD/covid_data/covid_500.fasta}"
if [ ! -f "$SRCALN" ]; then
    echo "ERROR: source alignment $SRCALN not found." >&2
    echo "Set SRCALN=/path/to/alignment.fasta to use a different one." >&2
    exit 1
fi

# ---------------------------------------------------------------- settings --
RESULTS="${RESULTS:-$PWD/acceptdist_results}"
REPS="${REPS:-1}"
NIT="${NIT:-15}"
SEEDS="${SEEDS:-777 1234 2025}"
# One thread per core the job actually got, capped: these runs are sequential
# and IQ-TREE's scaling past a handful of threads is poor on this workload.
THREADS="${THREADS:-$(( ${SLURM_CPUS_PER_TASK:-2} > 4 ? 4 : ${SLURM_CPUS_PER_TASK:-2} ))}"

# Seconds held back at the end for the final summary and log compression.
RESERVE="${RESERVE:-900}"
# Safety factor applied to the longest trial seen so far, before deciding
# there is room for another one.
SAFETY="${SAFETY:-1.15}"
# Used only for the FIRST trial, before any duration has been measured.
FIRST_TRIAL_ESTIMATE="${FIRST_TRIAL_ESTIMATE:-21600}"   # 6h

mkdir -p "$RESULTS" || exit 1
MANIFEST="$RESULTS/trials.tsv"
[ -f "$MANIFEST" ] || printf 'trial\tbootseed\tstarted\tfinished\tstatus\n' > "$MANIFEST"

# ------------------------------------------------------------ time budget ---
# Remaining seconds of walltime, or empty if it cannot be determined (in which
# case the loop falls back to MAXTRIALS and never guesses).
remaining_seconds () {
    local now end
    now=$(date +%s)
    # digits-only test first: [ -gt ] on a non-numeric value is a hard error,
    # not something 2>/dev/null suppresses
    case "${SLURM_JOB_END_TIME:-}" in
        ''|*[!0-9]*) : ;;
        *) echo $(( SLURM_JOB_END_TIME - now )); return ;;
    esac
    if [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        end=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
              | tr ' ' '\n' | sed -n 's/^EndTime=//p' | head -1)
        if [ -n "$end" ]; then
            end=$(date -d "$end" +%s 2>/dev/null)
            [ -n "$end" ] && { echo $(( end - now )); return; }
        fi
    fi
    echo ""
}

MAXTRIALS="${MAXTRIALS:-0}"        # 0 = unlimited (bounded by walltime instead)

# ------------------------------------------------------- graceful shutdown --
# SLURM sends SIGUSR1 early only if --signal was requested; SIGTERM always
# arrives at the kill. Either way, stop launching new trials and finish up.
STOPPING=0
finish_up () {
    echo
    echo "### wrapping up ($(date '+%F %T')) ###"
    # compress logs, which dominate the output size and are rarely read again
    find "$RESULTS" -name '*.log' -size +1M -exec gzip -f {} \; 2>/dev/null
    {
        echo "pooled summary written $(date '+%F %T')"
        echo "job=${SLURM_JOB_ID:-none} host=$(hostname)"
        echo
    } > "$RESULTS/POOLED_SUMMARY.txt"
    # the inner script's own pooled table, over every trial present
    BASEDIR="$RESULTS" SUMMARY_ONLY=1 bash "$INNER" --summary-only \
        >> "$RESULTS/POOLED_SUMMARY.txt" 2>/dev/null \
        || pooled_fallback >> "$RESULTS/POOLED_SUMMARY.txt"
    echo "results in: $RESULTS"
    echo "  trials.tsv            one row per trial attempted"
    echo "  POOLED_SUMMARY.txt    comparison pooled across every trial"
    echo "  t<seed>/              per-trial runs (.searchstats.txt, .treefile, .log.gz)"
    sed -n '1,80p' "$RESULTS/POOLED_SUMMARY.txt"
}

# Standalone pooled table, used if the inner script has no --summary-only mode.
# Reads only the .searchstats.txt files, so it cannot be broken by a reworded
# log message.
pooled_fallback () {
    "${PY:-python3}" - "$RESULTS" <<'PYEOF'
import glob, os, re, sys, collections
base = sys.argv[1]
rows = {}
for f in glob.glob(os.path.join(base, "t*", "*.searchstats.txt")):
    kv = dict(p.split("=", 1) for p in open(f).read().split() if "=" in p)
    m = re.match(r"^(t\w+)_b(\d+)_s(\d+)_(nni|spr)_(.+)$", kv.get("prefix", ""))
    if m:
        rows[(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5))] = kv
if not rows:
    print("no completed runs found")
    raise SystemExit
trials = sorted({k[0] for k in rows})
blocks = sorted({(k[0], k[1], k[2]) for k in rows})
arms = ["ctrl", "ad05c", "ad05a", "ad2a", "ad2s2a"]
print("trials=%d  paired blocks=%d" % (len(trials), len(blocks)))
for fam in ("nni", "spr"):
    fb = [b for b in blocks if any((b[0], b[1], b[2], fam, a) in rows for a in arms)]
    if not fb:
        continue
    print()
    print("--- %s family, pooled over %d block(s) ---" % (fam.upper(), len(fb)))
    print("%-10s %11s %10s %12s" % ("arm", "beats ctrl", "blocks", "mean cpu_s"))
    for a in arms:
        cpus = [float(rows[(b[0], b[1], b[2], fam, a)]["cpu_seconds"])
                for b in fb if (b[0], b[1], b[2], fam, a) in rows]
        if a == "ctrl":
            if cpus:
                print("%-10s %11s %10d %12.0f" % ("ctrl", "-", len(cpus), sum(cpus)/len(cpus)))
            continue
        better = tot = 0
        for b in fb:
            ka, kc = (b[0], b[1], b[2], fam, a), (b[0], b[1], b[2], fam, "ctrl")
            if ka in rows and kc in rows:
                tot += 1
                if float(rows[ka]["best_logl"]) > float(rows[kc]["best_logl"]):
                    better += 1
        if tot:
            print("%-10s %11d %10d %12.0f"
                  % (a, better, tot, sum(cpus)/len(cpus) if cpus else float("nan")))
PYEOF
}

on_signal () {
    echo
    echo "!! signal received at $(date '+%F %T') -- no new trials will start."
    echo "!! the trial in flight keeps running; it is checkpointed and can be"
    echo "!! resumed in a later job with its BOOTSEED (see trials.tsv)."
    STOPPING=1
}
trap on_signal USR1 TERM INT

# ------------------------------------------------------------------- main ---
echo "=================================================================="
echo " accept-dist bootstrap sweep"
echo " job=${SLURM_JOB_ID:-none}  host=$(hostname)  started=$(date '+%F %T')"
echo " results=$RESULTS"
echo " REPS=$REPS NIT=$NIT THREADS=$THREADS SEEDS='$SEEDS'"
R=$(remaining_seconds)
if [ -n "$R" ]; then
    echo " walltime remaining: ${R}s ($(( R / 3600 ))h $(( (R % 3600) / 60 ))m)"
else
    echo " walltime remaining: unknown (no SLURM time info; using MAXTRIALS=$MAXTRIALS)"
fi
echo "=================================================================="
echo

LONGEST=0
NTRIALS=0

while : ; do
    [ "$STOPPING" = "1" ] && { echo "stopping: signal received"; break; }
    if [ "$MAXTRIALS" != "0" ] && [ "$NTRIALS" -ge "$MAXTRIALS" ]; then
        echo "stopping: MAXTRIALS=$MAXTRIALS reached"; break
    fi

    R=$(remaining_seconds)
    if [ -n "$R" ]; then
        # Room for another trial? Use the longest one seen so far, times a
        # safety factor, plus the reserve held back for wrap-up. Before any
        # trial has finished, fall back to FIRST_TRIAL_ESTIMATE.
        if [ "$LONGEST" -gt 0 ]; then
            NEED=$(awk -v l="$LONGEST" -v s="$SAFETY" 'BEGIN{printf "%d", l*s}')
        else
            NEED="$FIRST_TRIAL_ESTIMATE"
        fi
        if [ "$R" -lt $(( NEED + RESERVE )) ]; then
            echo "stopping: ${R}s left, need ~$(( NEED + RESERVE ))s for another trial"
            break
        fi
    elif [ "$MAXTRIALS" = "0" ]; then
        echo "ERROR: cannot determine remaining walltime and MAXTRIALS=0 (unlimited)." >&2
        echo "Set MAXTRIALS=<n> so the loop has a bound." >&2
        break
    fi

    BOOTSEED=$(( (RANDOM * 32768 + RANDOM) % 1000000000 ))
    TRIAL="t${BOOTSEED}"
    NTRIALS=$(( NTRIALS + 1 ))
    T0=$SECONDS
    echo "################ trial $NTRIALS: $TRIAL ################"
    printf '%s\t%s\t%s\t\t%s\n' "$TRIAL" "$BOOTSEED" "$(date '+%F %T')" "started" >> "$MANIFEST"

    BASEDIR="$RESULTS" BOOTSEED="$BOOTSEED" REPS="$REPS" NIT="$NIT" \
        SEEDS="$SEEDS" THREADS="$THREADS" SRCALN="$SRCALN" IQ="$IQ" \
        bash "$INNER"
    RC=$?

    DUR=$(( SECONDS - T0 ))
    [ "$DUR" -gt "$LONGEST" ] && LONGEST="$DUR"
    printf '%s\t%s\t\t%s\t%s\n' "$TRIAL" "$BOOTSEED" "$(date '+%F %T')" \
        "finished rc=$RC dur=${DUR}s" >> "$MANIFEST"
    echo "################ trial $TRIAL done: rc=$RC dur=${DUR}s ################"
    echo

    # The generated alignment is ~15 MB and is fully reproducible from the
    # BOOTSEED recorded above, so compress it rather than keeping it raw.
    gzip -f "$RESULTS/$TRIAL/aln/"*.fasta 2>/dev/null
done

finish_up
echo
echo "=== JOB COMPLETE: $NTRIALS trial(s) attempted ==="
