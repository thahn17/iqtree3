# Search-strategy experiments — CLI reference

Every search variant added in this line of work, how to run it, and what it
writes. All of them are optional: with none of these flags, `iqtree3` behaves
exactly as it always has.

Each run writes `<prefix>.searchstats.txt` — one line of `key=value` pairs
holding the configuration and the result, so a sweep can be tabulated without
regex-scraping the human-readable `.log`.

---

## 1. SPR as the refinement stage

Replaces IQ-TREE's NNI hill-climb with an SPR search; the perturbation,
candidate set and stopping rule are untouched.

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-refine "radius 10 fast quiet weightprune"
```

`--nni-refine "record"` keeps the default NNI hill-climb but turns on the same
reporting, so NNI and SPR runs land in comparable spreadsheets.

**Recommended settings**, from the bootstrap study in `covid_500`:
`radius 10 fast fullreopt 100 20 weightprune`. Notably **not** `distradius`
(it halved improvement) and **not** `investigate` (a net negative at radius 1
and worse at radius 3).

---

## 2. One continuous SPR stage

Builds the candidate set as usual, then runs a single long SPR hill-climb from
the best candidate instead of the perturb/refine loop.

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F \
        --spr-continuous "radius 10 fast quiet fullreopt 100 20 weightprune steps 30000"
```

---

## 3. SPR perturbation (instead of the random-NNI kick)

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6"
```

Accepts kick-shaping flags only: `radius N`, `distradius`,
`weightprune [long]`, `slack D [anneal]`.

---

## 4. `slack` — bound the kick's damage

Without it a kick is blind: no score is consulted, so its total damage is
unbounded. With it, each proposed kick move is scored and kept only if it
costs at most `D` log-likelihood; a costlier draw is rolled back and redrawn
(up to 8 attempts per slot).

```bash
# SPR kick, constant bound
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6 slack 5"

# SPR kick, bound annealed to 0 over the run
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6 slack 5 anneal"

# NNI kick (IQ-TREE's default kick) -- use the kick-agnostic spelling
iqtree3 -s aln.fasta -st DNA -m GTR+F --perturb-slack 5
iqtree3 -s aln.fasta -st DNA -m GTR+F --perturb-slack 5 anneal
```

`--perturb-slack` and `--spr-perturb "slack ..."` are the same rule; the
former is the only way to apply it to the **NNI** kick.

Reports: `slack summary : N applied, M refused, over K kicks; delta D
<constant|annealed>, last effective E`.

---

## 5. `tunnel` — two-move lookahead in the refinement

On a *marginally* losing SPR step (within `TOL`), keep it applied and probe up
to `K` random follow-ups one at a time, looking for a **pair** that beats the
pre-move score. Failed probes roll back individually, leaving the tolerated
first move in place; if all fail, that move is rolled back too. Nothing is
ever kept unless it strictly improves on where the tree already was.

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F \
        --spr-continuous "radius 10 fast quiet fullreopt 100 20 weightprune tunnel 2 3 steps 30000"
```

`--spr-refine`-family only. Reports `tunnel summary`.

---

## 6. `--accept-dist` — stochastic acceptance, no perturbation

Replaces the refinement's strict improve-only rule with a probabilistic one
**and suppresses the perturbation stage entirely**. A move losing `|dL|` is
kept with probability

```
P = exp(-(|dL| / T)^shape)
```

and any improvement is kept outright (`P = 1`). So `(-inf, 0)` is a tunable
curve and `[0, inf)` is flat at 1.

```bash
# defaults: temp 0.5, shape 1, no annealing
iqtree3 -s aln.fasta -st DNA -m GTR+F --accept-dist

# explicit, annealed, with the NNI refiner
iqtree3 -s aln.fasta -st DNA -m GTR+F --accept-dist "temp 2 shape 1 anneal"

# same rule, SPR refiner
iqtree3 -s aln.fasta -st DNA -m GTR+F \
        --spr-refine "radius 10 fast quiet weightprune" --accept-dist "temp 2 anneal"

# sharper cutoff, annealing to a non-zero floor
iqtree3 -s aln.fasta -st DNA -m GTR+F --accept-dist "temp 3 shape 2 anneal floor 0.1"
```

| token | meaning | default |
|---|---|---|
| `temp T` | scale of a tolerable loss | 0.5 |
| `shape S` | falloff exponent | 1 |
| `anneal` | decay `T` toward `floor` over the run | off |
| `floor F` | lowest `T` annealing may reach | 0 |

`shape = 1` is the classic Boltzmann/Metropolis rule; `shape = 2` stays near 1
well inside `T` then falls off a cliff; `shape → ∞` becomes a hard threshold
at `|dL| = T` — which is exactly what `slack` does, making slack the limiting
case of this same family rather than an unrelated idea.

**Why the default `T` is 0.5.** It targets final convergence, not
exploration: a move 0.5 worse is kept with P = 0.37, one 2 worse with
P = 0.018, one 5 worse with P = 4.5e-5. Small backward steps that let a search
round a barrier get through; genuinely destructive moves effectively never do.
Raise `T` to explore. With `anneal` and the default floor of 0, the run ends
as a strict hill-climb however exploratory it began.

**Safety.** Accepted downhill moves let the current tree drift below the best
one seen, so the best topology is tracked and restored before the result is
reported. The stats line records whether that restore fired
(`bestseen_restored=yes|no`).

**Both refiners.** The rule is applied at each refiner's own accept/reject
point — the SPR step loop's keep/revert test, and the NNI refiner's
`evaluateNNIs` filter — so it composes with `--spr-refine` or the default NNI
refinement with no further configuration.

---

## Reading the results

```bash
cat run.searchstats.txt
```

```
prefix=AD_NNI seed=777 iterations=10 best_logl=-84382.29963 cpu_seconds=534.28
wall_seconds=356.53 refiner=nni kick=none accept_dist=on accept_temp=2
accept_shape=1 accept_anneal=yes accept_floor=0 accept_last_temp=0.2
downhill_kept=3561 downhill_offered=12934 bestseen_restored=no slack=off
tunnel=off candidates_evaluated=31450
```

To tabulate a sweep:

```bash
cat *.searchstats.txt | tr ' ' '\n' | grep -E '^(prefix|best_logl|cpu_seconds)='
```

**Use `cpu_seconds`, not `wall_seconds`, for any timing comparison.** Wall
clock on a laptop is not trustworthy for unattended runs — one run in this
work recorded 24,008 s wall against 2,276 s CPU (9.5% utilisation) because the
machine slept through it, which would have made a like-for-like comparison
meaningless.
