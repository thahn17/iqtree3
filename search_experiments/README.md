# Search experiments

Every search variant added in this line of work: how to run it, what it
writes, and where its code lives. All of it is optional. With none of these
flags, `iqtree3` runs IQ-TREE's own search unchanged. On a regression run it
produces byte-identical output to the build from before this code was moved
here, and the same tree and log-likelihood as the official 3.1.3 release
binary.

Each run writes `<prefix>.searchstats.txt`: one line of `key=value` pairs
holding the configuration and the result, so a sweep can be tabulated without
regex-scraping the human-readable `.log`.

---

## Layout

Everything lives in this directory:

| path | contents |
|---|---|
| `sprsearch.h`, `sprsearch.cpp` | the SPR search machinery shared by IQ-TREE and `spr_topology_test` |
| `iqtree_search_experiments.cpp` | every `IQTree` member the experiments add, including the hooks below |
| `phylotree_spr.cpp` | `PhyloTree::isLegalSPR` / `applySPR` / `rollbackSPR` |
| `search_params.h`, `search_params.cpp` | the flags: their fields, defaults, parsing and `--help` text |
| `spr_topology_test.cpp`, `spr_topology_test_usage.txt` | standalone SPR test program, built as `build/spr_topology_test` |
| `quartet_topology_test.cpp` | standalone quartet test program, built as `build/quartet_topology_test` |
| `CMakeLists.txt` | adds the sources above to IQ-TREE's `tree`/`utils` libraries; builds the two test programs |
| `scripts/` | experiment harnesses: `run_acceptdist_bootstrap.sh`, `slurm_acceptdist_bootstrap.sh`, `slurm_misof_hillclimb.sh`, `record_iqtree_nni.{sh,ps1}` |
| `scripts/analysis/` | post-processing for the record CSVs and their charts |
| `test_data/spr/` | fixtures for `spr_topology_test` |

Run the scripts from the repository root, e.g.
`bash search_experiments/scripts/run_acceptdist_bootstrap.sh`.

### How it attaches to IQ-TREE

IQ-TREE's own files keep only hook calls. Each call site is marked
`[search-experiments hook]`, and every hook does nothing unless its flag was
given: no output, no random draw, no change to the tree.

```bash
grep -n "search-experiments hook" tree/iqtree.cpp utils/tools.cpp
```

- **`tree/iqtree.h`** declares the experiment members in one block, between
  the `SEARCH EXPERIMENTS` and `END SEARCH EXPERIMENTS` markers.
- **`utils/tools.h`** has `class Params : public SearchExperimentParams`.
  The flag fields are declared in `search_params.h`, yet are still read as
  `params.spr_perturb` etc.
- **`tree/phylotree.h`** declares the SPR move API (`SPRMove`,
  `SPRRollback`, `isLegalSPR`, ...). Upstream's dormant `SPRMove` prototype
  was renamed `SPRMoveLegacy` to free the name.

Other edits to upstream files that are not hooks:

- **Root `CMakeLists.txt`:** `add_subdirectory(search_experiments)`, plus
  static MinGW linking for `iqtree3.exe` on Windows.
- **`pll/hardware.c`:** a MinGW `cpuid` build fix.
- **`tree/phylonode.h`:** `swapPartialLhState`, used by the SPR local
  likelihood cache.
- **`tree/quartet.cpp`:** a duplicate-quartet cache in
  `computeQuartetLikelihoods`, the `--lmap` engine `quartet_topology_test`
  uses.

---

## 0. Combining flags — read this first

The sections below are numbered for reference, **not** because they are
mutually exclusive. `--spr-refine` and `--spr-perturb` in particular are not
alternatives you switch between: they are two independent settings on the
*same* loop, and a normal run passes both at once.

Every iteration of IQ-TREE's iterated local search does:

```
doTreePerturbation();                          // <-- --spr-perturb controls this half
if (refine_mode == REFINE_SPR) doSPRSearch();
else                           doNNISearch();  // <-- --spr-refine  controls this half
addTreeToCandidateSet(...);
```

So there are four combinations:

| flags | kick | refinement |
|---|---|---|
| *(none)* | random NNI | NNI |
| `--spr-perturb "radius 6"` | random SPR | NNI |
| `--spr-refine "..."` | random NNI | SPR |
| **both** | random SPR | SPR |

A fully-SPR run is one command line carrying both:

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F \
        --spr-perturb "radius 6" \
        --spr-refine "radius 10 fast quiet fullreopt 100 20 weightprune"
```

**Two flags break this pattern:**

- **`--spr-continuous` replaces the whole loop.** It runs one uninterrupted
  SPR hill-climb from the best candidate tree and never perturbs. Combining it
  with `--spr-perturb` is meaningless, because there is no perturbation stage
  left.
- **`--accept-dist` suppresses the perturbation half** and leaves the
  refinement half intact. It works with `--spr-refine` or with default NNI
  refinement. Pairing it with `--spr-perturb` is pointless because the kick
  never runs.

Which pairs are meaningful:

| combination | valid | note |
|---|---|---|
| `--spr-perturb` + `--spr-refine` | ✅ | the full SPR configuration |
| `--spr-perturb` + `--perturb-slack` | ✅ | bounds the SPR kick |
| `--perturb-slack` alone | ✅ | bounds the **NNI** kick |
| `--accept-dist` + `--spr-refine` | ✅ | stochastic accept, SPR refinement |
| `--accept-dist` alone | ✅ | stochastic accept, NNI refinement |
| `--accept-dist` + `--spr-perturb` | ⚠️ | kick is suppressed; the flag does nothing |
| `--spr-continuous` + `--spr-perturb` | ⚠️ | no perturbation stage exists |
| `--perturb-slack` + `--spr-perturb "... slack D"` | ❌ | rejected: give the kick's bound once |
| `tunnel` in a `--spr-refine`/`--spr-continuous` spec | ✅ | refinement-side only |
| `slack` in a `--spr-refine` spec | ❌ | rejected: it is a kick flag |

---

## 1. SPR as the refinement stage

Replaces IQ-TREE's NNI hill-climb with an SPR search. The perturbation,
candidate set and stopping rule are untouched.

```bash
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-refine "radius 10 fast quiet weightprune"
```

`--nni-refine "record"` keeps the default NNI hill-climb but turns on the same
reporting, so NNI and SPR runs land in comparable spreadsheets.

**Recommended settings**, from the bootstrap study on `covid_500`:
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

The best single result on `covid_500` so far: +55.1 logL over stock NNI at
3.3× its CPU (one seed, 9102026). It is the `cont` arm of the bootstrap
harness below, so it can be replicated.

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
costs at most `D` log-likelihood. A costlier draw is rolled back and redrawn,
up to 8 attempts per slot.

```bash
# SPR kick, constant bound
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6 slack 5"

# SPR kick, bound annealed to 0 over the run
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6 slack 5 anneal"

# NNI kick (IQ-TREE's default kick) -- the kick-agnostic spelling
iqtree3 -s aln.fasta -st DNA -m GTR+F --perturb-slack 5
iqtree3 -s aln.fasta -st DNA -m GTR+F --perturb-slack 5 anneal

# ...which also works with the SPR kick
iqtree3 -s aln.fasta -st DNA -m GTR+F --spr-perturb "radius 6" --perturb-slack 5
```

`--perturb-slack` and `--spr-perturb "slack ..."` are the same rule.
`--perturb-slack` is the only way to apply it to the **NNI** kick. With
`anneal` the bound is `D × (1 − p)` on the shared clock (see *Annealing*
below), so the last kick is strictly non-worsening.

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

`shape = 1` is the classic Boltzmann/Metropolis rule. `shape = 2` stays near 1
well inside `T`, then falls off a cliff. As `shape → ∞` it becomes a hard
threshold at `|dL| = T`, which is exactly what `slack` does. So slack is the
limiting case of this same family, not an unrelated idea.

**Why the default `T` is 0.5.** It targets final convergence, not
exploration. A move 0.5 worse is kept with P = 0.37, one 2 worse with
P = 0.018, and one 5 worse with P = 4.5e-5. Small backward steps that let a
search get around a barrier get through; genuinely destructive moves
effectively never do. Raise `T` to explore. With `anneal` and the default
floor of 0, the run ends as a strict hill-climb however exploratory it began.

**Best-seen restore.** Accepted downhill moves let a refinement end below the
best tree it passed through. Both refiners therefore remember the best tree
*each refinement* passes through, and hand that back instead of wherever the
walk stopped. `bestseen_restores=N` in the stats line counts how many
refinements needed it.

**Both refiners.** The rule is applied at each refiner's own accept/reject
point: the SPR step loop's keep/revert test, and the NNI refiner's
`evaluateNNIs` filter. It works with `--spr-refine` or the default NNI
refinement with no further configuration.

---

## Annealing: one clock for everything that anneals

`slack ... anneal`, `--perturb-slack D anneal` and `--accept-dist ... anneal`
all read the same clock, `p ∈ [0, 1]`: 0 on the first main-loop iteration, 1
on the last. The bound is `D × (1 − p)` and the temperature is
`max(floor, T × (1 − p))`.

| stopping rule | when is "the last iteration"? |
|---|---|
| `-n N` | exactly iteration `N − 1`, so the schedule ends on time |
| default / `--nstop K` | `K` iterations after the last improvement, i.e. the projected stopping point. It moves out whenever the search improves, so the schedule re-warms a little then. |

Starting-tree iterations (`--ntop`) are not part of the schedule.

> **Results from before this change used different schedules.** In runs up to
> commit `a16009cb`:
>
> - Slack annealing counted kicks against `min_iterations`, so with
>   `-n 15`/`-n 20 --ntop 5` it only got down to 30–40 % of `D`.
> - Accept-dist counted all iterations, running from about 67 % to 7 % of `T`.
> - Under the default stopping rule, neither annealed at all.
> - The NNI refiner had no best-seen restore.
> - The SPR refiner's best-seen tree was run-wide, so a refinement could be
>   handed back a tree from an earlier iteration.
>
> Every annealed-slack and accept-dist result recorded before then carries
> that behaviour, and should be re-run before being compared with new ones.

---

## Reading the results

```bash
cat run.searchstats.txt
```

```
prefix=AD_NNI seed=777 iterations=10 best_logl=-84382.29963 cpu_seconds=534.28
wall_seconds=356.53 refiner=nni kick=none accept_dist=on accept_temp=2
accept_shape=1 accept_anneal=yes accept_floor=0 accept_last_temp=0
downhill_kept=3561 downhill_offered=12934 bestseen_restored=yes
bestseen_restores=4 slack=off tunnel=off candidates_evaluated=31450
```

To tabulate a sweep:

```bash
cat *.searchstats.txt | tr ' ' '\n' | grep -E '^(prefix|best_logl|cpu_seconds)='
```

**Use `cpu_seconds`, not `wall_seconds`, for any timing comparison.** Wall
clock on a laptop is not trustworthy for unattended runs. One run in this work
recorded 24,008 s wall against 2,276 s CPU (9.5 % utilisation) because the
machine slept through it, which would have made a like-for-like comparison
meaningless. CPU time is also only comparable on the same machine: a cluster
node's AVX+FMA kernel and a laptop's SSE2 kernel differ several-fold.

---

## Bootstrap harness

`scripts/run_acceptdist_bootstrap.sh` generates its own bootstrap replicates of
an alignment (covid_500 by default) and runs two families, each against the
control it should be judged by:

| family | control | arms |
|---|---|---|
| NNI | stock IQ-TREE | `--accept-dist` × 4 variants; **`cont`** = `--spr-continuous "... weightprune steps 30000"` |
| SPR | `--spr-perturb` + `--spr-refine "... weightprune"` | `--spr-refine "... weightprune" --accept-dist` × 4 variants |

```bash
bash search_experiments/scripts/run_acceptdist_bootstrap.sh              # one new trial
CONT=0 bash search_experiments/scripts/run_acceptdist_bootstrap.sh       # without the continuous arm
BASEDIR=acceptdist_results bash search_experiments/scripts/run_acceptdist_bootstrap.sh --summary-only
```

Each run is independent: `BOOTSEED` is random unless you set it, so each
invocation adds a fresh bootstrap trial. Finished runs are skipped and
unfinished ones resume from their checkpoint.
`scripts/slurm_acceptdist_bootstrap.sh` loops it until the job's walltime is
nearly gone. Submit it from the repository root; results go to
`acceptdist_results/`, or to `$RESULTS` if set.
