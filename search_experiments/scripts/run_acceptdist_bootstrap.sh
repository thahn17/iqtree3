#!/usr/bin/env bash
#
# Compare --accept-dist against the perturb/refine controls, on bootstrap
# replicates of an alignment (COVID by default).
#
# TWO FAMILIES, so each accept-dist arm is judged against a control that
# differs ONLY in how moves are accepted -- not in which refiner is running:
#
#   NNI family   control = stock IQ-TREE          (NNI kick  + NNI refinement)
#                arms    = --accept-dist          (no kick   + NNI refinement)
#
#   SPR family   control = --spr-perturb + --spr-refine "...weightprune"
#                                                  (SPR kick  + SPR refinement)
#                arms    = --spr-refine "...weightprune" --accept-dist
#                                                  (no kick   + SPR refinement)
#
# Both controls are the perturb+refine COMBINATION, because that pairing is
# what was actually being run before this feature existed -- comparing against
# a refinement-only run would flatter accept-dist by removing a stage the real
# baseline has.
#
# TRIALS. The script generates its own bootstrap alignments, seeded by
# BOOTSEED. Change BOOTSEED and you get a completely independent trial: new
# alignments, a new output directory, and new run prefixes. Nothing from a
# previous trial is reused or overwritten, so trials can be run repeatedly and
# compared. Within one trial the alignments are written once and then left
# alone, so re-running the same trial resumes rather than regenerating.
#
# Resumable: a finished run is detected from its own log and skipped, and an
# unfinished one restarts WITHOUT -redo so IQ-TREE resumes from its .ckp.gz.
# Killing this script and re-running it loses at most the run in flight.
#
# Usage:
#   bash test_scripts/run_acceptdist_bootstrap.sh
#       new randomly-seeded trial: 1 fresh bootstrap alignment x 3 search
#       seeds = 3 paired blocks. Run it again to ADD another trial; the
#       pooled summary at the end covers every trial found so far.
#
#   BOOTSEED=12345 bash test_scripts/run_acceptdist_bootstrap.sh
#       resume/repeat a SPECIFIC trial (the script prints its BOOTSEED at
#       startup for exactly this purpose)
#
#   REPS=3 NIT=30 bash test_scripts/run_acceptdist_bootstrap.sh
#   SEEDS="777 1234" bash test_scripts/run_acceptdist_bootstrap.sh
#   SRCALN=/path/to/other.fasta bash ...                # different source data
#
set -u

# --summary-only / SUMMARY_ONLY=1: print the pooled table over whatever trials
# already exist and exit WITHOUT running or generating anything. The wrapper
# job calls this during wrap-up, where accidentally starting a fresh trial
# would be actively harmful.
SUMMARY_ONLY="${SUMMARY_ONLY:-0}"
for _a in "$@"; do
  [ "$_a" = "--summary-only" ] && SUMMARY_ONLY=1
done

# ----------------------------------------------------------------- config --
IQ="${IQ:-C:/Users/tinst/Desktop/iqtree3/build/iqtree3.exe}"
SRCALN="${SRCALN:-C:/Users/tinst/Desktop/iqtree3/covid_data/covid_500.fasta}"

# Trial identity. BOOTSEED seeds the site resampling; every alignment, output
# path and run prefix is keyed to it, so two BOOTSEEDs never collide.
#
# RANDOM BY DEFAULT, deliberately: each invocation draws a fresh trial, so
# running this script repeatedly ACCUMULATES independent bootstrap samples
# instead of re-running the same ones. The value is printed below and can be
# passed back explicitly (BOOTSEED=<value>) to resume an interrupted trial --
# without that, a re-run would start a new trial rather than resuming.
BOOTSEED="${BOOTSEED:-$(( RANDOM * 32768 + RANDOM ))}"
TRIAL="${TRIAL:-t${BOOTSEED}}"

BASEDIR="${BASEDIR:-$(pwd)/acceptdist_bootstrap}"
OUTDIR="${OUTDIR:-$BASEDIR/$TRIAL}"
ALNDIR="${ALNDIR:-$OUTDIR/aln}"

# One fresh alignment per invocation, crossed with three search seeds, gives
# 3 paired blocks per run -- enough to see a direction without a single
# invocation running for a day. Re-run the script to add more.
REPS="${REPS:-1}"                    # bootstrap replicates per trial
NIT="${NIT:-15}"                     # -n, iterations per run
SEEDS="${SEEDS:-777 1234 2025}"      # IQ-TREE seeds; varies the SEARCH
THREADS="${THREADS:-2}"      # keep low: runs are sequential, RAM is the limit

# Shared by both families so only the accept rule differs within a family.
INITF="--ninit 10 --ntop 5"
MODEL="-m GTR+F"
SPRSPEC="radius 10 fast quiet fullreopt 100 20 weightprune"
PERTURBSPEC="radius 6"

# The accept-dist variants under test. name|spec
#  - default temperature, no annealing: the convergence-oriented default
#  - default temperature, annealed
#  - hotter and annealed: more exploratory
#  - shape 2: near-flat inside T then a sharp cliff, instead of Boltzmann
ADVARIANTS=(
  "ad05c|temp 0.5"
  "ad05a|temp 0.5 anneal"
  "ad2a|temp 2 anneal"
  "ad2s2a|temp 2 shape 2 anneal"
)

if [ "$SUMMARY_ONLY" = "1" ]; then
  cd "$BASEDIR" 2>/dev/null || { echo "no results directory at $BASEDIR"; exit 0; }
else
mkdir -p "$OUTDIR" "$ALNDIR"

# ------------------------------------------------- bootstrap alignments ----
# Written once per trial and then left alone: regenerating them mid-trial
# would silently invalidate every run already completed against them.
echo "=================================================================="
echo " TRIAL $TRIAL   (BOOTSEED=$BOOTSEED)"
echo " to resume THIS trial after an interruption, re-run with:"
echo "     BOOTSEED=$BOOTSEED bash test_scripts/run_acceptdist_bootstrap.sh"
echo " running it without BOOTSEED starts a NEW trial instead"
echo "=================================================================="
echo "Bootstrap alignments (source=$(basename "$SRCALN"))"
NEEDGEN=0
for r in $(seq 1 "$REPS"); do
  [ -f "$ALNDIR/boot$r.fasta" ] || NEEDGEN=1
done

if [ "$NEEDGEN" = "1" ]; then
  python - "$SRCALN" "$ALNDIR" "$REPS" "$BOOTSEED" <<'PYEOF'
import os, random, sys

src, outdir, nrep, bootseed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

names, seqs, cur = [], [], []
with open(src) as fh:
    for line in fh:
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
            names.append(line[1:])
        elif line:
            cur.append(line)
if cur:
    seqs.append("".join(cur))

n = len(seqs[0])
assert all(len(s) == n for s in seqs), "ragged alignment"
print("  source: %d taxa x %d sites" % (len(seqs), n))

for r in range(1, nrep + 1):
    path = os.path.join(outdir, "boot%d.fasta" % r)
    if os.path.exists(path):
        print("  keep   boot%d.fasta (already present)" % r)
        continue
    # Seed is BOOTSEED + r, so trial 1000 and trial 2000 draw disjoint
    # resamplings rather than shifted copies of the same ones.
    rng = random.Random(bootseed + r)
    cols = [rng.randrange(n) for _ in range(n)]
    with open(path, "w") as out:
        for nm, s in zip(names, seqs):
            out.write(">" + nm + "\n")
            out.write("".join([s[c] for c in cols]) + "\n")
    print("  wrote  boot%d.fasta (seed %d, %d distinct source columns)"
          % (r, bootseed + r, len(set(cols))))
PYEOF
  if [ $? -ne 0 ]; then
    echo "!! alignment generation failed -- aborting"
    exit 1
  fi
else
  echo "  all $REPS alignments already present, reusing"
fi
echo

cd "$OUTDIR" || exit 1

# ---------------------------------------------------------------- helpers --
have_result () { [ -f "$1.log" ] && grep -q "BEST SCORE FOUND" "$1.log" 2>/dev/null; }
score ()   { grep -h "BEST SCORE FOUND" "$1.log" 2>/dev/null | tail -1 | sed 's/.*: //'; }
cputime () { grep -hE "^Total CPU time used" "$1.log" 2>/dev/null | tail -1 | awk '{print $5}'; }

# run <prefix> <alignment> <seed> <extra args...>
run () {
  local P="$1"; shift
  local A="$1"; shift
  local S="$1"; shift
  if have_result "$P"; then
    echo "  SKIP   $P  logL=$(score "$P")  cpu=$(cputime "$P")s"
    return
  fi
  [ -f "$P.ckp.gz" ] && echo "  RESUME $P (from checkpoint)"
  local t0=$SECONDS
  "$IQ" -s "$A" -st DNA $MODEL -nt "$THREADS" $INITF -n "$NIT" -seed "$S" \
        "$@" --prefix "$P" > /dev/null 2>&1
  local rc=$?
  echo "  DONE   $P  exit=$rc  logL=$(score "$P")  cpu=$(cputime "$P")s  wall=$((SECONDS-t0))s"
}

# ------------------------------------------------------------------- main --
echo "accept-dist vs perturb/refine controls"
echo "  trial=$TRIAL  replicates=$REPS  iterations=$NIT  seeds='$SEEDS'  threads=$THREADS"
echo "  output=$OUTDIR"
echo

for r in $(seq 1 "$REPS"); do
  A="$ALNDIR/boot$r.fasta"
  [ -f "$A" ] || { echo "!! missing $A -- skipping replicate $r"; continue; }

  for S in $SEEDS; do
    echo "######## replicate $r / seed $S ########"
    TAG="${TRIAL}_b${r}_s${S}"

    # ---- NNI family ----------------------------------------------------
    echo " [nni] control: stock IQ-TREE (NNI kick + NNI refinement)"
    run "${TAG}_nni_ctrl" "$A" "$S"
    for v in "${ADVARIANTS[@]}"; do
      name="${v%%|*}"; spec="${v#*|}"
      echo " [nni] accept-dist: $spec"
      run "${TAG}_nni_${name}" "$A" "$S" --accept-dist "$spec"
    done

    # ---- SPR family ----------------------------------------------------
    echo " [spr] control: --spr-perturb + --spr-refine (SPR kick + SPR refinement)"
    run "${TAG}_spr_ctrl" "$A" "$S" --spr-perturb "$PERTURBSPEC" --spr-refine "$SPRSPEC"
    for v in "${ADVARIANTS[@]}"; do
      name="${v%%|*}"; spec="${v#*|}"
      echo " [spr] accept-dist: $spec"
      run "${TAG}_spr_${name}" "$A" "$S" --spr-refine "$SPRSPEC" --accept-dist "$spec"
    done

    echo "--- replicate $r / seed $S complete ---"
    echo
  done
done

fi   # end of the acting section; the summaries below also run under
     # --summary-only

# ---------------------------------------------------------------- summary --
# Built from the .searchstats.txt files rather than by scraping the logs, so a
# reworded log message cannot silently break the table.
if [ "$SUMMARY_ONLY" != "1" ]; then
echo "=== SUMMARY (trial $TRIAL) ==="
python - <<'PYEOF'
import glob, re, collections

rows = {}
for f in sorted(glob.glob("*.searchstats.txt")):
    kv = dict(p.split("=", 1) for p in open(f).read().split() if "=" in p)
    m = re.match(r"^(t\w+)_b(\d+)_s(\d+)_(nni|spr)_(.+)$", kv.get("prefix", ""))
    if m:
        rows[(int(m.group(2)), int(m.group(3)), m.group(4), m.group(5))] = kv

if not rows:
    print("no .searchstats.txt files yet")
    raise SystemExit

# a "block" is one paired comparison: same alignment, same seed
blocks = sorted({(k[0], k[1]) for k in rows})
arms = ["ctrl", "ad05c", "ad05a", "ad2a", "ad2s2a"]

for fam in ("nni", "spr"):
    present = [a for a in arms if any((r, s, fam, a) in rows for r, s in blocks)]
    if not present:
        continue
    print()
    print("--- %s family (logL; higher is better; * = best in row) ---" % fam.upper())
    print("%-12s" % "rep/seed" + "".join("%16s" % a for a in present))
    wins = collections.Counter()
    for (r, s) in blocks:
        cells, best, bestarm = [], None, None
        for a in present:
            kv = rows.get((r, s, fam, a))
            v = float(kv["best_logl"]) if kv else None
            cells.append(v)
            if v is not None and (best is None or v > best):
                best, bestarm = v, a
        if bestarm:
            wins[bestarm] += 1
        print("%-12s" % ("b%d/s%d" % (r, s)) + "".join(
            ("%16s" % "-") if v is None else
            (("%15.3f*" % v) if a == bestarm else ("%16.3f" % v))
            for a, v in zip(present, cells)))

    def meancpu(a):
        vals = [float(rows[(r, s, fam, a)]["cpu_seconds"])
                for (r, s) in blocks if (r, s, fam, a) in rows]
        return sum(vals) / len(vals) if vals else float("nan")
    print("%-12s" % "mean cpu_s" + "".join("%16.0f" % meancpu(a) for a in present))
    print("%-12s" % "wins" + "".join("%16d" % wins[a] for a in present))

    # head-to-head against the control, which is the comparison that matters
    if "ctrl" in present:
        print("  vs control:")
        for a in present:
            if a == "ctrl":
                continue
            better = sum(1 for (r, s) in blocks
                         if (r, s, fam, a) in rows and (r, s, fam, "ctrl") in rows
                         and float(rows[(r, s, fam, a)]["best_logl"])
                           > float(rows[(r, s, fam, "ctrl")]["best_logl"]))
            total = sum(1 for (r, s) in blocks
                        if (r, s, fam, a) in rows and (r, s, fam, "ctrl") in rows)
            if total:
                print("    %-8s beats control in %d/%d block(s)" % (a, better, total))
PYEOF
fi
echo
echo "=== POOLED ACROSS ALL TRIALS IN $BASEDIR ==="
python - "$BASEDIR" <<'PYEOF2'
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
print("trials=%d  paired blocks=%d  (%s)" % (len(trials), len(blocks), ", ".join(trials)))

for fam in ("nni", "spr"):
    fb = [b for b in blocks if any((b[0], b[1], b[2], fam, a) in rows for a in arms)]
    if not fb:
        continue
    print()
    print("--- %s family, pooled ---" % fam.upper())
    print("%-10s %10s %10s %12s" % ("arm", "beats ctrl", "of blocks", "mean cpu_s"))
    for a in arms:
        if a == "ctrl":
            continue
        better = tot = 0
        cpus = []
        for b in fb:
            ka, kc = (b[0], b[1], b[2], fam, a), (b[0], b[1], b[2], fam, "ctrl")
            if ka in rows:
                cpus.append(float(rows[ka]["cpu_seconds"]))
            if ka in rows and kc in rows:
                tot += 1
                if float(rows[ka]["best_logl"]) > float(rows[kc]["best_logl"]):
                    better += 1
        if tot:
            print("%-10s %10d %10d %12.0f"
                  % (a, better, tot, sum(cpus) / len(cpus) if cpus else float("nan")))
    cc = [float(rows[(b[0], b[1], b[2], fam, "ctrl")]["cpu_seconds"])
          for b in fb if (b[0], b[1], b[2], fam, "ctrl") in rows]
    if cc:
        print("%-10s %10s %10d %12.0f" % ("ctrl", "-", len(cc), sum(cc) / len(cc)))
PYEOF2
echo
if [ "$SUMMARY_ONLY" = "1" ]; then
  echo "=== SUMMARY ONLY: nothing was run and no trial was created ==="
else
  echo "=== ALL DONE (trial $TRIAL) ==="
fi
