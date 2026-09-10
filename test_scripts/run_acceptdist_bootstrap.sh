#!/usr/bin/env bash
#
# Compare --accept-dist against the perturb/refine controls, on bootstrap
# replicates of the COVID alignment.
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
# Resumable: a finished run is detected from its own log and skipped, and an
# unfinished one restarts WITHOUT -redo so IQ-TREE resumes from its .ckp.gz.
# Killing this script and re-running it loses at most the run in flight.
#
# Usage:
#   bash test_scripts/run_acceptdist_bootstrap.sh
#   REPS=5 NIT=30 bash test_scripts/run_acceptdist_bootstrap.sh
#   OUTDIR=/path/to/results bash test_scripts/run_acceptdist_bootstrap.sh
#
set -u

# ----------------------------------------------------------------- config --
IQ="${IQ:-C:/Users/tinst/Desktop/iqtree3/build/iqtree3.exe}"
BOOTDIR="${BOOTDIR:-C:/Users/tinst/AppData/Local/Temp/claude/c--Users-tinst-Desktop-iqtree3/ee861b91-f8cd-46b1-8e13-e75f82d98841/scratchpad/boot}"
OUTDIR="${OUTDIR:-$(pwd)/acceptdist_bootstrap}"
REPS="${REPS:-3}"            # bootstrap replicates to use (boot1..bootN)
NIT="${NIT:-15}"             # -n, iterations per run
SEED="${SEED:-777}"
THREADS="${THREADS:-2}"      # keep low: runs are sequential, RAM is the limit

# Shared by both families so only the accept rule differs within a family.
INITF="--ninit 10 --ntop 5"
MODEL="-m GTR+F"
SPRSPEC="radius 10 fast quiet fullreopt 100 20 weightprune"
PERTURBSPEC="radius 6"

# The accept-dist variants under test. Name|spec
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

mkdir -p "$OUTDIR"
cd "$OUTDIR" || exit 1

# ---------------------------------------------------------------- helpers --
have_result () { [ -f "$1.log" ] && grep -q "BEST SCORE FOUND" "$1.log" 2>/dev/null; }
score ()   { grep -h "BEST SCORE FOUND" "$1.log" 2>/dev/null | tail -1 | sed 's/.*: //'; }
cputime () { grep -hE "^Total CPU time used" "$1.log" 2>/dev/null | tail -1 | awk '{print $5}'; }

# run <prefix> <alignment> <extra args...>
run () {
  local P="$1"; shift
  local A="$1"; shift
  if have_result "$P"; then
    echo "  SKIP   $P  logL=$(score "$P")  cpu=$(cputime "$P")s"
    return
  fi
  [ -f "$P.ckp.gz" ] && echo "  RESUME $P (from checkpoint)"
  local t0=$SECONDS
  "$IQ" -s "$A" -st DNA $MODEL -nt "$THREADS" $INITF -n "$NIT" -seed "$SEED" \
        "$@" --prefix "$P" > /dev/null 2>&1
  local rc=$?
  echo "  DONE   $P  exit=$rc  logL=$(score "$P")  cpu=$(cputime "$P")s  wall=$((SECONDS-t0))s"
}

# ------------------------------------------------------------------- main --
echo "accept-dist vs perturb/refine controls"
echo "  replicates=$REPS  iterations=$NIT  seed=$SEED  threads=$THREADS"
echo "  output=$OUTDIR"
echo

for r in $(seq 1 "$REPS"); do
  A="$BOOTDIR/boot$r.fasta"
  if [ ! -f "$A" ]; then
    echo "!! missing $A -- skipping replicate $r"
    continue
  fi
  echo "######## bootstrap replicate $r ########"

  # ---- NNI family -------------------------------------------------------
  echo " [nni] control: stock IQ-TREE (NNI kick + NNI refinement)"
  run "b${r}_nni_ctrl" "$A"

  for v in "${ADVARIANTS[@]}"; do
    name="${v%%|*}"; spec="${v#*|}"
    echo " [nni] accept-dist: $spec"
    run "b${r}_nni_${name}" "$A" --accept-dist "$spec"
  done

  # ---- SPR family -------------------------------------------------------
  echo " [spr] control: --spr-perturb + --spr-refine (SPR kick + SPR refinement)"
  run "b${r}_spr_ctrl" "$A" --spr-perturb "$PERTURBSPEC" --spr-refine "$SPRSPEC"

  for v in "${ADVARIANTS[@]}"; do
    name="${v%%|*}"; spec="${v#*|}"
    echo " [spr] accept-dist: $spec"
    run "b${r}_spr_${name}" "$A" --spr-refine "$SPRSPEC" --accept-dist "$spec"
  done

  echo "--- replicate $r complete ---"
  echo
done

# ---------------------------------------------------------------- summary --
# Built from the .searchstats.txt files rather than by scraping the logs, so
# a reworded log message cannot silently break the table.
echo "=== SUMMARY ==="
python - <<'PYEOF'
import glob, os, re, collections

rows = {}
for f in sorted(glob.glob("*.searchstats.txt")):
    kv = dict(p.split("=", 1) for p in open(f).read().split() if "=" in p)
    m = re.match(r"b(\d+)_(nni|spr)_(.+)$", kv.get("prefix", ""))
    if m:
        rows[(int(m.group(1)), m.group(2), m.group(3))] = kv

if not rows:
    print("no .searchstats.txt files yet")
    raise SystemExit

reps = sorted({k[0] for k in rows})
arms = ["ctrl", "ad05c", "ad05a", "ad2a", "ad2s2a"]

for fam in ("nni", "spr"):
    present = [a for a in arms if any((r, fam, a) in rows for r in reps)]
    if not present:
        continue
    print()
    print("--- %s family (logL; higher is better) ---" % fam.upper())
    print("%-6s" % "rep" + "".join("%16s" % a for a in present))
    wins = collections.Counter()
    for r in reps:
        cells, best, bestarm = [], None, None
        for a in present:
            kv = rows.get((r, fam, a))
            v = float(kv["best_logl"]) if kv else None
            cells.append(v)
            if v is not None and (best is None or v > best):
                best, bestarm = v, a
        if bestarm:
            wins[bestarm] += 1
        print("%-6d" % r + "".join(
            ("%16s" % "-") if v is None else
            (("%15.3f*" % v) if a == bestarm else ("%16.3f" % v))
            for a, v in zip(present, cells)))
    print("%-6s" % "cpu_s" + "".join(
        "%16s" % (("%.0f" % (sum(float(rows[(r, fam, a)]["cpu_seconds"]) for r in reps
                                 if (r, fam, a) in rows) /
                             max(1, sum(1 for r in reps if (r, fam, a) in rows)))))
        for a in present))
    print("%-6s" % "wins" + "".join("%16d" % wins[a] for a in present))
    ctrl_wins = wins.get("ctrl", 0)
    total = sum(wins.values())
    print("  control won %d of %d replicate(s); best arm: %s"
          % (ctrl_wins, total, max(wins, key=wins.get) if wins else "-"))
PYEOF
echo
echo "=== ALL DONE ==="
