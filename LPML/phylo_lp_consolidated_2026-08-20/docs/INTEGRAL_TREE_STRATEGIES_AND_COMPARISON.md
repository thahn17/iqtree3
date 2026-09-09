# Integral-tree recovery strategies: evidence, rationale, and remaining GTR comparison

## Purpose

This document consolidates the experiments used to turn the continuous fractional phylogenetic LP into a small set of integral-tree candidates. The objective is not to interpret one LP extreme point as a literal posterior distribution over trees. The objective is to generate a compact candidate set that contains high-likelihood integral trees, then rank those candidates with the actual alignment likelihood after continuous branch-length optimization.

The main production target is approximately **100 taxa** with alignments from roughly **1,000 bp upward**, often much longer.

## Current model context

The current branch model has two non-JC alternatives that must be distinguished:

1. **Two-point GTR relaxation**: each biological lower arm interpolates between a short and a long transition effect. It is restrictive but often produces concentrated topology marginals that are easy for bottleneck decoding.
2. **Spectral GTR relaxation**: a reversible non-JC transition model is decomposed into its nonstationary spectral decay modes. This covers branch effects from near zero through the stationary/infinite-length endpoint. It is the more physically complete GTR branch-effect relaxation, but its additional freedom can leave topology variables less concentrated.

JC remains simpler because its nonstationary decay rates coincide, so the branch effect effectively has one decay mode.

## Why integral-tree recovery is difficult

Several experiments showed that raw unit activity is not a calibrated probability that an integral tree contains that unit. The LP returns one optimal extreme point of a large relaxation. Interchangeable physical copies, weakly coupled internal topology variables, and branch-effect freedom can all change which extreme point is returned without a comparable change in objective.

The most useful recovery methods therefore share three properties:

- they do **not** greedily remove fractional mass as the primary tree-identification step;
- they preserve the original bottleneck candidates as a safety baseline;
- they rank final candidates by the **actual optimized integral-tree likelihood**, not by the fractional score.

---

# 1. Strategies tested

## 1.1 Greedy residual-mass peeling

### Idea
Construct one complete integral tree, subtract its supported mass from the fractional tree, and repeat.

### Reason it was considered
The fixed-flow construction gives each active subtree a frequency-like activity, so literal extraction initially appeared to provide interpretable mixture proportions.

### What happened
It could remove support shared by several near-identical trees. In exact synthetic mixtures it could explain 100% of root mass while recovering a different mixture of trees. Soft peeling improved some mixture-recovery metrics, but the result depended on extraction order.

### Decision
**Do not use as the main integral-tree identification strategy.** It remains only a descriptive decomposition experiment.

---

## 1.2 Full tree-space NNI/SPR search on the fractional score

### Idea
Treat the fractional-tree match score as an objective and search completed integral trees by NNI/SPR moves.

### Reason it was considered
It avoids committing to a fixed pre-generated candidate set.

### Problem
There is no demonstrated monotonicity of the fractional matching score under NNI/SPR. A better tree may require an intermediate decrease. Search cost also grows poorly with taxon count.

### Decision
Do not use as the main fractional decoder. NNI remains useful **after** a candidate is integral, when moves are accepted only if the real optimized alignment likelihood improves.

---

## 1.3 Root-down / leaves-up hierarchical construction

### Idea
Construct the topology directly in the fixed-flow hierarchy rather than traversing completed-tree space.

- Root-down: choose supported split units recursively.
- Leaves-up: begin with known taxa and merge supported child subtrees.
- Bidirectional/beam variants retain several alternatives when local decisions are ambiguous.

### Reason it helps
A partial-tree support bound is monotone nonincreasing as more required subtrees/leaves are added. This gives meaningful beam pruning that NNI/SPR on the fractional score lacks.

### Small-case performance
On an 8-taxon fractional solution, a root-down beam recovered 100% of the structurally extractable root mass as four trees, with the largest component 37.5%; leaves-up independently recovered the same total. On one 12-taxon case both directions again recovered essentially all root mass as four trees.

### Limitation
Hard support from a single LP extreme point can exclude a plausible topology entirely even when a nearby LP solution supports it. The hierarchy is therefore best used as a **candidate generator**, not as a probability model.

---

## 1.4 G2/local-fragment hierarchy + bottleneck ranking

### Idea
Retain fixed-depth local co-occurrence: a split is scored together with the immediate child decisions below it. Generate trees non-destructively, then rank them by the **weakest supported local split** rather than multiplying many supports.

### Why bottleneck was chosen
Product/log-product scores can punish a good tree for one uncertain split and can reward a hybrid built from many moderately supported pieces. The minimum local support is conservative and performed better on difficult synthetic mixtures.

### Alignment-backed best-known-tree ranks
From `bestknown_hierarchy_summary.csv`:

| Taxa | Alignment | Cases | Worst observed rank of best-known tree |
|---:|---:|---:|---:|
| 20 | 20 bp | 1 | 5 |
| 20 | 50 bp | 5 | 5 |
| 20 | 200 bp | 5 | 1 |
| 20 | 1000 bp | 5 | 1 |
| 50 | 20 bp | 1 | 3 |
| 50 | 50 bp | 3 | 10 |
| 50 | 200 bp | 3 | 1 |
| 50 | 1000 bp | 3 | 1 |
| 100 | 20 bp | 1 | 27 |
| 100 | 50 bp | 3 | 30 |
| 100 | 200 bp | 1 | 1 |
| 100 | 1000 bp | 1 | 1 |

These experiments used a synthetic/local-fragment representation richer than the most aggressively pruned production LP, so they are evidence for the **ranking strategy**, not proof that the production LP always supplies equally good G2 information.

### Decision
**Keep raw bottleneck as the mandatory baseline candidate list.**

---

## 1.5 Root objective profiling and near-optimal-face sampling

### Idea
Do not treat zero activity in one LP optimum as zero plausibility. Profile alternative root types by forcing them and measuring objective loss, or sample several extreme points of the near-optimal face.

### Findings
Root profiling is meaningful and cheap because only O(L) root split types need to be tested. Deeper profiling was not useful in the aggressively pruned LP: in one 8-taxon test, 15 of 18 non-root split units could be forced to activity 1 with essentially zero objective loss.

Averaging or centering sampled near-optimal vertices was not robust; more probes were not monotonically better.

### Decision
Root-profile loss may be retained as a **supplemental root-diversification signal**. Do not recursively profile a large portion of the unit library and do not use a sampled near-optimal polytope barycenter as the construction object.

---

## 1.6 Local-fragment viability as a secondary objective

### Idea
After the primary likelihood LP reaches value z*, solve a second LP on the **same topology variables**:

    maximize summed local-fragment viability
    subject to primary LP constraints
               primary score >= z* - Delta.

The local viability of a split jointly reflects both child compositions and branch effects. This gives all units a likelihood-informed usefulness score while requiring them to coexist on one shared topology.

### Exact local-score small-case results
From `local_fragment_coupling_clean.csv`, with Delta = 2:

| Case | Spectral baseline | Spectral + exact local coupling | Two-point + exact local coupling | Reference |
|---|---:|---:|---:|---:|
| 8 taxa moderate #0 | -310.174 | **-299.979** | -303.759 | -291.779 |
| 8 taxa moderate #1 | -302.641 | **-284.819** | -294.918 | -282.182 |
| 8 taxa wide | -295.612 | **-294.922** | -295.411 | -294.999 |
| 12 taxa moderate | -320.307 | **-305.108** | -326.264 | -270.747 |

In all four completed paired exact-local tests, spectral + local coupling beat two-point + the same coupling. The 12-taxon gap to the reference remained large, so this is not yet a complete recovery solution.

### Enumeration issue
The strongest exact local score enumerated child composition possibilities and is not appropriate at 100 taxa. Saturated fixed-order approximations were tested. At 8 taxa, T=1 (presence-only) was often as robust as or more robust than T=2, including the wide-branch case.

### Decision
Local-fragment coupling is one of the most promising ways to make the spectral relaxation produce useful topology structure, but the scalable score construction still needs the full 100-taxon end-to-end test.

---

## 1.7 Enumeration-free local overlap score

### Idea
For child base-frequency vectors fL and fR, use a reference transition matrix Pbar and score

    sum_i pi_i min((Pbar fL)_i, (Pbar fR)_i).

It uses four small auxiliary score variables per unit/site and no composition-table enumeration.

### Finding
It does not numerically reproduce the exact enumerated local upper score; in some regimes the two are only moderately correlated. However it can still act as a useful topology-selection heuristic. It should be viewed as a separate fixed-cost usefulness score rather than an approximation guaranteed to preserve the old score numerically.

### Interior-vs-extreme secondary solutions
Blending the primary and secondary LP solutions sometimes improved bottleneck recovery:

    x_alpha = (1-alpha) x_primary + alpha x_secondary,

with alpha in {0, 0.5, 1}. In some cases alpha=0.5 was best; in others the primary point alpha=0 was best. Therefore the safe policy is to decode **all three views and union the candidates** rather than choosing a universal alpha.

---

# 2. Unit weighting and normalization

The local child score already uses base **fractions**, not raw counts, so a strong second flow normalization is usually harmful.

Define a unit split m=a+b and

    B = 4ab / m^2.

A useful coupled size/balance family is

    w = m^{-q} B^{-q} = (4ab/m)^{-q}.

This couples flow and balance; weighting balance alone can explode for large unbalanced splits.

## 2.1 20-30 taxon long-alignment results

From `unit_weighting_large_tree_summary.csv`, mild negative exponents often improved the known-tree rank, but no single exponent dominated on every replicate. This motivated a **multi-view** policy rather than replacing raw bottleneck.

## 2.2 100 taxa, 10,000 bp, moderate branches

From `unit_weighting_100taxa_summary.csv`, three replicates:

| beta=gamma | Mean true-tree rank | Worst rank | Mean RF of rank-1 tree |
|---:|---:|---:|---:|
| **0.00** | **3.67** | **5** | 5.33 |
| -0.05 | 4.67 | 5 | 5.33 |
| -0.10 | 4.33 | 7 | **4.00** |
| -0.15 | 4.00 | **5** | 6.00 |
| -0.20 | 16.00 | **40** | 9.33 |
| -0.30 | 7.00 | 10 | 6.00 |

Balance-only weighting was catastrophic: beta=0, gamma=-0.2 had mean true-tree rank about 108 and mean RF 196. Strong size-only alternatives were also poor.

### Decision at 100 taxa
- **Primary:** raw unweighted bottleneck.
- **Secondary structural view:** mild q about 0.10-0.15.
- Do not use q around 0.20 or stronger as the sole score.
- Never use a strong balance-only penalty.

## 2.3 Whole-alignment pairwise silhouette

For split A|B, a useful size-independent score is

    mean D(A,B) - max(mean D(A,A), mean D(B,B)).

It uses an O(L^2) alignment summary and is independent of the number of alignment columns after preprocessing.

At 20-30 taxa it improved average known-tree rank and RF distance, but not on every replicate. At 100 taxa with very wide/saturated branches it remained only a complementary view.

### Decision
Use pairwise/silhouette ranking only to **add candidates**. Do not replace raw bottleneck or put the pairwise coefficient directly into the secondary LP objective by default.

---

# 3. Plain spectral GTR versus old two-point GTR

This is the most important unresolved comparison.

## 3.1 Without local coupling: two-point is currently the more stable tree finder

From `spectral_bottleneck_paired_summary.csv`, 11 direct paired cases using the same bottleneck hierarchy:

- spectral wins: **2**
- spectral losses: **9**
- mean spectral-minus-two-point final optimized log-likelihood: **-7.70**
- median: **-1.56**

Moderate branch cases:
- 1 win / 6 cases
- mean difference: **-10.32**
- median: **-5.18**

Wide branch cases:
- 1 win / 5 cases
- mean difference: **-4.54**
- median: **-1.24**

Therefore **plain spectral + hard bottleneck is not as reliable as the old two-point GTR + bottleneck pipeline**.

## 3.2 Why spectral is still attractive

Spectral is the more complete non-JC branch-effect relaxation. It covers the full reversible-model decay curve from zero length toward the stationary/infinite-length limit instead of restricting each edge to one short/long interval.

The issue is not the spectral decomposition itself; the extra branch freedom can leave the pruned LP's internal topology variables less concentrated, which hurts a decoder that reads one optimal extreme point as hard support.

## 3.3 With local coupling: spectral becomes much more competitive

The small exact-local results above reverse the plain comparison: spectral + local coupling beat two-point + the same coupling in all four completed paired tests.

This is the strongest evidence that the right long-term approach may be:

    spectral GTR likelihood relaxation
    + cheap fixed-order local-unit coherence
    + multi-view bottleneck candidate generation
    + exact continuous branch optimization.

But the evidence is currently small-scale.

---

# 4. Current provisional policy

Until the 100-taxon end-to-end comparison is complete:

1. Keep **raw bottleneck** as a mandatory candidate list.
2. Add a mild coupled size/balance view with q approximately 0.10-0.15.
3. Add a whole-alignment pairwise/silhouette view when inexpensive.
4. If a secondary local-coherence LP is used, decode alpha in {0, 0.5, 1} and union the candidates.
5. Final ranking is always by the real GTR/J C likelihood after continuous branch-length optimization.
6. Do **not** remove the old two-point GTR result as a regression baseline yet.

A union policy has an important safety property: if the original candidate list is retained, adding alternative rankings cannot reduce the best final likelihood after exact rescoring.

---

# 5. Future definitive comparison: two-point GTR vs current spectral methods

The remaining comparison should be performed directly at the intended use case rather than inferred from small trees.

## Required dataset grid

At minimum:

- taxa: **100**
- alignment lengths: **10,000; 20,000; 50,000 bp**
- branch regimes:
  - moderate heterogeneous branches;
  - very wide heterogeneous branches;
  - one or more near-clock/uniform-branch cases;
- at least 3 independent simulated replicates per regime when runtime permits.

The same sampled alignment, taxon order, GTR+F parameters, and final branch optimizer must be used for both methods.

## Pipelines to compare

### Baseline A: old two-point GTR

    two-point GTR LP
    -> raw bottleneck/G2 hierarchy
    -> top-K candidates
    -> continuous branch optimization
    -> final log-likelihood ranking.

### Current B: spectral primary only

    spectral GTR LP
    -> raw bottleneck
    -> top-K
    -> continuous branch optimization.

This is expected to remain less stable than A in some regimes and is included to quantify the value of the coupling layer.

### Current C: spectral + multi-view recovery

    spectral primary LP
    -> secondary local-coherence LP with small Delta
    -> decode alpha = 0, 0.5, 1
    -> raw bottleneck + q=0.10/0.15 structural view
    -> optional pairwise silhouette view
    -> union candidates
    -> continuous branch optimization.

A scalable non-enumerating local score must be used for the 100-taxon test.

## Candidate budgets

Report at least:

- K = 20
- K = 40
- K = 80

For the multi-view method, report the **actual unique union size**, because the lists overlap substantially.

## Primary comparison metrics

The decisive metric is:

    best optimized integral-tree log likelihood found within the fixed total time/candidate budget.

Also report:

- whether the known generating/best-known topology was present;
- rank of the generating/best-known topology under each structural view;
- RF/clade distance of the top structurally ranked candidate;
- LP build/preprocessing time;
- LP solve time;
- secondary solve time;
- candidate-generation time;
- branch-optimization/rescoring time;
- total wall time;
- peak memory;
- number of variables/rows/nonzeros.

## No-regression criterion

The new spectral pipeline should replace the two-point baseline only if, at 100 taxa:

1. its best final optimized likelihood is **never worse** than the two-point baseline within the tested budget, or any losses are small enough to be recovered by a cheap retained two-point fallback;
2. its median/mean final likelihood is at least as good;
3. the improvement persists in both moderate and wide branch regimes;
4. the added preprocessing and solve time remain proportional to the total computation rather than dominating it.

## Provisional verdict today

**As a branch-effect relaxation:** spectral is better/more complete.

**As an input to the current hard bottleneck decoder without coupling:** two-point is better/more stable based on the 11 paired tests.

**With exact local-fragment coupling in small tests:** spectral is better in the completed paired examples, suggesting that the previous spectral failures were largely topology-underconstraint rather than a bad branch model.

**At 100 taxa:** the unit-weighting and pairwise-ranking components have been tested, but a definitive full LP + branch-rescoring comparison has not yet completed. Therefore the two-point method should remain the explicit regression baseline/fallback until that experiment is done.

---

# 6. Preprocessing modes retained in the project

The project also retains the preprocessing tradeoff experiments. In the saturated-count benchmark for a balanced 100-taxon 25/25/25/25 pattern:

- T=1 preprocessing: about **1.02 s**, root upper 6.137e-8;
- T=2 preprocessing: about **8.40 s**, root upper 5.064e-8.

T=2 remains the preferred accuracy/runtime compromise in that experiment. Exact composition/hull routines are retained as small-tree oracles, not as the intended 100-taxon preprocessing path.

---

# 7. Files supporting this document

Key result files in `research/results/`:

- `bestknown_hierarchy_summary.csv`
- `spectral_bottleneck_paired_summary.csv`
- `spectral_bottleneck_benchmark.csv`
- `local_fragment_coupling_clean.csv`
- `spectral_secondary_weighting_screen.csv`
- `unit_weighting_large_tree_summary.csv`
- `unit_weighting_100taxa_summary.csv`
- `unit_candidate_union.csv`
- `preprocessing_solve_tradeoff_summary.csv`
- `near_optimal_face_test_summary.csv`

Key research notes in `docs/`:

- `TREE_IDENTIFICATION_HIERARCHY.md`
- `TREE_RECOVERY_EVALUATION.md`
- `SPECTRAL_BRANCH_TESTS.md`
- `LOCAL_FRAGMENT_COUPLING_TESTS.md`
- `LOCAL_SCORE_SCALING_EXPERIMENTS.md`
- `UNIT_WEIGHTING_LARGE_TREE_TESTS.md`
- `UNIT_WEIGHTING_100TAXA_EXPERIMENTS.md`
- `NEAR_OPTIMAL_FACE_TESTS.md`
- `MODERATE_GTR_DIAGNOSTIC.md`

