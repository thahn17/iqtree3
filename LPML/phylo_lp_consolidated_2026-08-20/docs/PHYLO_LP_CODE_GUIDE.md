# Variable-topology phylogenetic LP: code modes and post-solve tree decomposition

## Preprocessing modes

`DEFAULT_PRODUCTION_PREPROCESS_MODE = "balanced"` in `multisite_pruned_general.py`.
The **currently wired runnable core** varies the root-envelope preprocessing budget as follows:

| Mode | Root envelope | Quadratic root facets | Root-plane preprocessing | Intended use |
|---|---|---:|---|---|
| `fast` | 4 shifted valid majorants (+ global plane) | 1 | `O(KN)` | minimum build latency |
| **`balanced` (default)** | **8 shifted valid majorants (+ global plane)** | **3** | **`O(KN)`** | normal production runs |
| `accuracy` | 16 shifted valid majorants (+ global plane) | 6 | `O(KN)` | spend more rows for a tighter root envelope |
| `oracle` | exact upper convex hull | 12 | potentially much larger / Qhull-dependent | small-L validation only |

The previously benchmarked **saturated-count** policies (`T=1,2,3`) remain a separate
preprocessing experiment and are not silently claimed by this constructor.  Those benchmarks
showed that fixed `T=2` is the preferred future replacement for the current exact composition-box
DP at large `L`, but wiring it through every retained root/moment feature would be a larger core
change than the FASTA/GTR/final-rescoring additions here.  Thus the runner's `--preprocess-mode`
currently controls the root-envelope/facet budget; it does **not** change the composition DP order.

The production criterion remains that preprocessing/assembly should eventually be
`O(nnz(final LP))` or smaller.  Exact composition + exact Qhull is therefore retained only as an
oracle path, not the intended 100-taxon endpoint.

## Branch-length convention

`fixed_flow_unit_generator_variable_lengths.py --variable-length-scope branches` now means **only the `left_arm` and `right_arm` of split units**.  These are the two parent-to-child biological branches.  `upper_arm`, taxon-unit internal edges, and Beneš helper edges do not get independent evolutionary lengths by default.

Branch-effect variables are shared across all alignment columns.  JC-like models use one long-path activity `ell[u,port]` per lower arm.  Non-JC reversible models use up to three homogenized spectral mode activities `q_r[u,port]`; a general four-state GTR model has three distinct nonstationary modes, while degeneracies automatically collapse to one or two unique modes.

## Post-solve fractional-tree representation and integral-tree extraction

The default extractor is now **hierarchical**.  It does not use NNI/SPR and does not
search a pre-generated candidate set.

```python
frac = model.build_fractional_tree(scipy_result)
mix = model.extract_integral_tree_mixture(
    scipy_result,
    max_components=20,
    root_beam=24,
    leaf_beam=20,
    child_branching=6,
    leaf_branching=4,
    merge_branching=48,
)
```

Implementation: `fractional_tree_hierarchical_units.py`.

### Fractional tree used by the decoder

The solved split units are first converted to a hierarchical fractional-tree support
object.  A virtual super-root has mass 1.  Each physical split unit `u` has frequency
`act[u]`; its two lower requests have the generator's known descendant flows and named-
taxon marginal vectors.

For an internal request of flow `k>1`, **any unused count-k split unit is structurally
compatible**.  The child unit's residual activity is its frequency support.  Taxon-
membership overlap is used only as a secondary/tie-break score.  Exact named-taxon
identity is enforced at flow 1.  Therefore a fractional mix of unresolved internal
subtrees is not penalized simply because a single leaf set has not yet been chosen.

A second diagnostic representation is available as

```python
dag = model.build_physical_fractional_dag(scipy_result)
```

which collapses the solved physical Beneš switch activities into request->child
fractional edges.  This is useful for inspecting the actual router solution, but it can
be more restrictive than the projected named-taxon relaxation, so it is not the default
decoding object.

### Root-down and leaves-up search

Two complementary bounded searches operate on the same fractional support:

* **root-down**: choose a root candidate, then recursively choose the most strongly
  supported count-matched child subtrees; at flow 1 choose named taxa;
* **leaves-up**: start from the integral taxa and repeatedly merge partial subtrees into
  the most strongly supported compatible split unit.

For a partial tree define

```text
theta = min(residual activity of every selected split unit,
            residual flow-1 taxon membership of every selected leaf attachment).
```

Extending a partial tree can only keep or lower `theta`, so this is a monotone pruning
bound.  Beam width gives bounded backtracking without traversing the whole rooted-tree
space.  Pure greedy is obtained with beam/branching width 1.

After a complete tree is found, its `theta` is subtracted from every selected split-unit
frequency and selected flow-1 taxon attachment.  A full integral tree removes balanced
request/child frequency at every fixed flow, so the residual frequencies remain
well-defined.  Repeating gives multiple integral components.  Their reported
`proportion` values are **extractable residual frequencies**, not coefficients fitted by
another LP/QP.  Any residual root mass is reported rather than forced into a tree.

### Frequency versus branch length

Variable branch effect remains mandatory in the likelihood LP.  Topology extraction uses

```text
arm_frequency[e] = act[u]
```

plus taxon-state structure, not the magnitude of a branch-length coordinate.  Under the JC/two-point
representation one may still inspect `ell[e]/act[u]` as a conditional long-path intensity.  Under the
non-JC spectral representation the shared coordinates are mode survival activities rather than a
literal scalar length, and the heavily pruned internal modes need not be uniquely identified.  They
are therefore not used to decide topology or to initialize internal fixed-tree branches.  Every
recovered integral topology has its actual scalar branch lengths re-estimated by the final ML scorer.

### Current directional checks

On the main 8-taxon fractional solution:

* root-down greedy extracted 12.5% before getting stuck;
* root-down beam extracted **100%** as four trees, largest component **37.5%**;
* leaves-up beam independently extracted **100%** as four trees;
* the previous NNI/Frank-Wolfe comparison took about 1.3 s and explained about 98.6% of
  its weighted feature norm, whereas hierarchical decoding took milliseconds.

On a 12-taxon test, root-down and leaves-up both extracted essentially **100%** as four
trees, with the largest component about **60%**.  Other random fractional LP points can
leave nonzero residual mass; this is expected when the relaxed marginals lie outside the
convex hull representable by positive-mass integral trees under the chosen hierarchy.

The older NNI/Frank-Wolfe implementation remains available only as
`extract_integral_tree_mixture_tree_space(...)`, and the even older finite-candidate
routine as `extract_integral_tree_mixture_legacy_candidates(...)`.

## Tree-recovery accuracy: ranking versus mixture peeling

The post-solve decoder now distinguishes two tasks that should not be conflated:

1. **Find a high-quality integral tree.**  Do **not** remove root mass first.  Use
   `rank_integral_tree_candidates(...)` to generate a non-destructive structural
   candidate list, then rescore those integral topologies with the actual alignment
   likelihood and fast branch-length optimization.
2. **Describe the fractional LP point as an approximate tree mixture.**  Residual
   peeling remains available, but its proportions are descriptive rather than a
   uniquely identifiable latent mixture.

### Why immediate peeling is not a reliable tree-identification criterion

Synthetic benchmarks were constructed by averaging known integral trees through the
same limited physical split-unit copies as the universal library.  Therefore unit
activities and taxon-membership vectors are exact linear averages of known trees before
controlled perturbation is added.

For exact mixtures of nearby trees, destructive peeling could remove 100% of root mass
while recovering only about 94% of the true topology mass on average.  A 50% soft peel
improved exact true-component recovery to about 98% in the same benchmark.  More
importantly, a strict decoder can reproduce *all* unit/taxon marginals exactly while
still returning a different topology mixture: first-order marginals do not retain
cross-level subtree co-occurrence correlations.

A concrete 6-taxon zero-noise example had true mixture weights 0.45/0.40/0.15 and an
alternative 0.60/0.25/0.15 mixture with the same clade-support and internal taxon-
membership marginals to numerical precision.  Their second-order nested-clade
co-occurrence differs (maximum difference 0.15).  Therefore the current fractional
representation is a valid first-order average of trees, but the latent component mixture
need not be unique.

### Non-destructive candidate recall

Across 6- and 8-taxon synthetic tests (24 cases per mixture/noise cell):

- a single underlying tree was rank 1 in every test through 8% perturbation;
- for near-tree 45/40/15 mixtures, the dominant true tree was in the top 3 in
  95.8% (zero noise), 91.7% (3% noise), and 100% (8% noise) of cases;
- it was in the top 10 in 100% of every tested cell.

This supports likelihood-rescoring a small *non-destructive* candidate list instead of
committing the first structural match.

### Large-taxon direction choice

Activity-first root-down ranking becomes ambiguous as the number of physical units and
possible hybrids grows.  A taxon-aware leaves-up decoder now uses internal named-taxon
membership amounts as soon as a partial subtree's taxon set is known.  It considers only
the best few left/right subtree matches for each physical split unit.

`best_tree_leaves_up_taxon_aware(...)` has a bounded local beam; it does not traverse
NNI/SPR tree space.  In targeted stress cases:

- 100 taxa, far 50/30/20 mixture, 3% perturbation: greedy leaves-up failed, beam 4
  recovered the dominant true tree exactly in about 3.5 s;
- 50 taxa, far mixture, 8% perturbation: beam 4 recovered it exactly in about 0.4 s;
- 20 taxa, near mixture, 3% perturbation: beam 8 recovered it exactly in about 0.05 s.

The recommended candidate generator is therefore bidirectional but non-destructive:

```python
candidates = rank_integral_tree_candidates(
    fractional_tree,
    topk=10,
    root_beam=256,
    leaves_beam=4,
)
```

For production, rescore every returned topology using the true alignment likelihood.
Branch effect remains mandatory inside the LP, but matching uses branch frequency `act[u]`.  JC/two-point conditional intensity and non-JC spectral mode coordinates are ignored for topology identity and branch lengths are reoptimized after the topology is selected.

### Fractional-average diagnostic

A stricter small-tree diagnostic can require every taxon below an extracted internal
port to have enough residual membership to support the component mass.  On exact
synthetic averages it extracted 100% root mass.  With controlled 3% and 8% perturbation,
mean strictly extractable mass fell to roughly 95.7% and 88.9%, respectively.  This is
useful as an **integral-convexity diagnostic**, not as the default large-tree search
criterion: a relaxed LP point can legitimately lie outside the convex hull encoded by
those first-order marginals.

## Non-destructive integral-tree identification (current recommendation)

Residual/root-mass peeling is **not** the primary tree-identification diagnostic.  It is a greedy packing operation and can consume support shared by several near-identical trees.  The production interpretation is instead non-destructive:

1. Keep LP branch effect in the likelihood model, but use lower-arm frequency `a_u` rather than length/mode intensity for topology matching.
2. Build a fractional hierarchy from split-unit activities, taxon-state marginals, and (where available) local parent/child routing support.
3. Compose the hierarchy with a fixed local-fragment model.  The preferred tested form keeps the two immediate child decisions jointly for each parent production (fixed context order; it does not grow with tree height).
4. Generate a top-K list top-down (leaves-up as a complementary beam when local taxon information is more decisive).
5. Rank conservatively by the weakest supported split (bottleneck score), using product/log support only as a tie-breaker.
6. Re-score those K integral topologies by the actual alignment likelihood with the fast integral-tree branch-length optimizer.  Do not assign/removal proportions before this re-scoring.

### Empirical candidate-depth rule

In an alignment-backed JC benchmark with an intentionally harsh 15% structural-contamination component, the target was the best topology found by a multi-start exhaustive-NNI neighborhood search, not merely the winner of a pre-generated candidate pool.  The largest observed hierarchy ranks were:

- 20 taxa / 50 bp: 5
- 50 taxa / 50 bp: 10
- 100 taxa / 50 bp: 30
- 20--100 taxa / 200 or 1000 bp: 1 in every completed case
- representative 20 bp cases: 5 (20 taxa), 3 (50 taxa), 27 (100 taxa)

A conservative one-pass production rule supported by these tests is therefore:

`K_score = 3` when hierarchy uncertainty is low, otherwise `K_score = max(ceil(n_taxa/3), ceil(N_eff))`,

where `N_eff` is an effective-tree-count/entropy diagnostic from the non-destructive hierarchy.  At the target 100 taxa this normally means scoring roughly 30--45 integral candidates in the uncertain regime, which is cheap compared with the LP.

This is empirical, not a mathematical guarantee.  If the hierarchy is exceptionally diffuse, increase K rather than peel mass.

### Fractional-tree composition fidelity

First-order split/clade marginals can recombine into hybrid trees that were not major components of the underlying average.  A fixed local-fragment composition reduces this.  In the controlled averages, retaining the joint pair of immediate child decisions increased probability mass assigned to the actual contributing trees and reduced hybrid leakage by roughly 5--10 percentage points in the difficult large-tree cases.  It also improved 100-taxon/50-bp ML rank (e.g. zero-order rank 4 to local-fragment rank 2 in one 5% contamination group).

The next implementation target for an actual LP solution is to reconstruct these local fragment weights from the shared topology flow plus taxon-compatibility marginals (no new per-site variables), then use the same non-destructive top-K decoder.

## FASTA production runner, GTR+F, final top-K likelihoods

The end-to-end entry point is now:

```bash
python phylo_lp_runner.py alignment.fa \
  --out-prefix run1 \
  --preprocess-mode balanced \
  --top-k 40
```

The five timed phases are:

1. FASTA read, optional random column sampling, and substitution-model setup;
2. fixed-flow unit-library creation (unless `--library` is supplied) plus LP build/preprocessing;
3. continuous LP solve;
4. non-destructive hierarchical integral-tree identification;
5. continuous branch-length optimization and true alignment-likelihood rescoring of every recovered top-K tree.

The runner writes `<prefix>.timings.json`, `<prefix>.topk.tsv`, `<prefix>.best.tree`, and
`<prefix>.summary.json`.  If columns are sampled, it also writes `<prefix>.sampled.fa` and
records the original 1-based column indices in the summary.

### Whole alignment versus random column sample

By default all FASTA columns are used.  To draw a reproducible sample without replacement:

```bash
python phylo_lp_runner.py alignment.fa --max-columns 500 --seed 17
```

If the alignment contains fewer than 500 columns, the whole alignment is used.  Duplicate
sampled site patterns are compressed internally by the LP and by the final likelihood scorer.
The compression is exact at the site-pattern level: a full across-taxon column tuple is stored
once and receives an integer multiplicity weight equal to the number of sampled columns with
that tuple.  The LP objective, log-chord error accounting, and final Felsenstein likelihood all
multiply that pattern's contribution by this weight.  The run summary reports
`unique_exact_site_patterns`, `unique_lp_site_patterns`, and `duplicate_columns_merged_lp`.
Because the LP maps non-ACGT ambiguity codes to uninformative `N` before compression, two
columns that differ only in ambiguity coding can merge in the LP even though the final scorer
continues to distinguish their IUPAC state sets.

Standard DNA ambiguity characters are accepted by the final likelihood scorer.  In the LP,
non-ACGT characters are conservatively treated as uninformative/null leaves with partial vector
`(1,1,1,1)`; the null count is implicit and does not add a fifth composition coordinate.

### JC69, GTR+F, and LP branch-effect defaults

JC69 remains available with:

```bash
--model jc
```

GTR+F uses empirical frequencies from the **selected FASTA columns** and six user-supplied
exchangeabilities in the order `AC,AG,AT,CG,CT,GT`:

```bash
--model gtr+f --gtr-rates 1.0,2.1,0.7,0.8,1.4,0.9
```

The resulting generator is normalized to one expected substitution per unit branch length.
The same `Q` and stationary-frequency vector are used by LP preprocessing and by final
integral-tree likelihood evaluation.

`--branch-relaxation auto` is the default.  It dispatches as follows:

- JC-like transition model: the existing one-scalar/two-point branch representation;
- any non-JC reversible transition model: the **three-mode spectral branch relaxation**.

The old non-JC two-point representation is retained only as a diagnostic/compatibility option:

```bash
--branch-relaxation two-point
```

The spectral relaxation can be requested explicitly with:

```bash
--branch-relaxation spectral --spectral-support-points 7
```

For reversible non-JC `Q`, write

```text
P(t) = 1*pi^T + z1(t) B1 + z2(t) B2 + z3(t) B3
zr(t) = exp(-mu_r t),  0 <= t < infinity.
```

A physical branch therefore traces a one-dimensional curve in `(z1,z2,z3)`.  The LP does
**not** require a scalar branch length.  For an arm with activity `a`, it stores the homogenized
mode activities

```text
q_r = a * z_r.
```

A small 3-D polyhedron outer-approximates the convex hull of all physical branch effects.  Its
facet orientations are obtained from a sampled spectral curve, but every right-hand side is
recomputed as the maximum over the **continuous** curve, so every real branch length from zero
through infinity remains feasible.  Linear transition-entry nonnegativity constraints ensure
that every relaxed spectral vertex is still a stochastic transition matrix.

With the default 7 support points, the tested GTR model had 26 inequalities and 24 vertices in
its shared 3-D branch-effect polytope.  These constraints are generated once from the model and
then homogenized on every active lower arm.  At 100 taxa / 4012 split units / 8024 lower arms:

```text
old two-point branch variables:       8,024
spectral branch variables:            24,072
additional shared variables:          16,048
spectral branch-polytope rows:        about 26 * 8,024 = 208,624
```

Thus the variable increase is small relative to the million-variable shared topology/taxon
layer.  The current prototype's larger cost is root-envelope preprocessing, because a
7-support spectral polytope has 24 vertices per arm instead of two short/long endpoints.
The root-cloud code is vectorized and slope fitting is subsampled while every resulting plane
is still shifted against all branch/composition points; on the 8-taxon/60-column stress test
this reduced spectral LP build time from about 24 s to about 9 s.

#### Spectral branch validation and observed effects

For the tested GTR+F generator, reconstructing `P(t)` from the three spectral modes matched
`expm(Q*t)` to about `1.4e-14` over branch lengths spanning zero through very large values, and
all sampled physical effects satisfied the LP branch polytope.

Two 8-taxon synthetic GTR alignments were used to compare the old two-point proxy with the
new full-range spectral default.  The first deliberately contained true branches from `1e-4`
through `10` substitutions/site; the second used only moderate branches around `0.05--0.24`.
The spectral relaxation produced substantially larger (less restrictive) LP scores, which is
expected because the old proxy excluded most of the physical branch-length range:

```text
wide-branch 60 bp:
  two-point LP score   -396.51
  spectral LP score    -163.91

moderate 80 bp:
  two-point LP score   -446.39
  spectral LP score    -203.79
```

The score increase is **not by itself a regression**: the spectral model contains physical
branch effects that the old interval simply excluded.  A real issue found during this test was
post-LP branch initialization.  Most internal spectral-mode variables in the heavily pruned
root-only model are weakly likelihood-coupled, so using their arbitrary LP values to initialize
fixed-tree branch optimization can start L-BFGS-B in the infinite-branch saturation plateau.
The final scorer now ignores those arbitrary internal spectral values and uses deterministic
finite multistarts (`--branch-multistart 4` by default).

After that correction, the best top-20 tree from the wide-branch spectral run optimized to
about `-482.94`, versus about `-518.43` for the top-10 two-point run and `-544.29` for the known
generating topology after its own branch reoptimization.  On the moderate-length test, however,
the two-point fractional hierarchy still produced a better top-20 candidate (`-639.39`) than
the spectral hierarchy (`-670.91`).  This shows the expected tradeoff: full branch coverage
can make the fractional topology less concentrated even while making the likelihood relaxation
more physically complete.  For that reason `two-point` remains available as a diagnostic, but
`auto` still selects spectral for non-JC as the production default.

### Final continuous branch-length optimization

Variable branch effect remains mandatory inside the LP.  For recovered integral candidates,
branch lengths are optimized as actual scalar `t` values; the spectral LP coordinates are not
reported as literal lengths.

`phylo_branch_opt.py` implements Felsenstein pruning plus an analytic gradient and bounded
L-BFGS-B on log branch lengths.  A deterministic finite multistart set is used by default to
avoid saturation-start artifacts:

```bash
--branch-multistart 4
```

The default lower bound is `1e-12`, effectively zero.  The default upper bound is `auto`: from
the nonzero eigenvalues of the normalized substitution generator it chooses

```text
t_max = max(10, -log(1e-12) / spectral_gap(Q))
```

so that even the slowest nonstationary CTMC mode is at most `1e-12` at the upper boundary.
For JC69 this is about 20.72 substitutions/site; extreme GTR+F models automatically receive a
larger range.  `--branch-max X` can still impose an explicit user bound.  Thus the final
integral-tree optimization spans effectively zero through a model-dependent numerical infinity.

The final `.topk.tsv` is sorted by these optimized alignment log likelihoods, not by the
fractional structural score.  `<prefix>.best.tree` is the best rescored topology in Newick
format with optimized branch lengths.

The structural search is non-destructive.  It uses bottleneck support as the primary hierarchy
score and local/taxon support as the secondary score.  If fewer than the requested `K` complete
positive-support trees exist in the searched fractional hierarchy, all recovered trees are
rescored and reported; zero-support trees are not fabricated merely to fill K.

### Timing fields

`<prefix>.timings.json` contains:

```text
01_input_model_s
02_library_lp_build_s
03_lp_solve_s
04_tree_identification_s
05_branch_opt_rescore_s
total_s
```

The same timings are printed while the run executes.  The intent is to keep preprocessing
visible and comparable to solve time rather than hiding expensive root/composition work in a
single opaque build step.
