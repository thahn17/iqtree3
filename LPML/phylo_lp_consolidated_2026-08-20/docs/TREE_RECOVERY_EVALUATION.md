# Fractional-tree recovery evaluation

## Main conclusions

- Root mass explained is not the same as correct topology recovery.
- Immediate destructive peeling can create hybrid components by consuming marginal
  support shared by several near-identical trees.
- The fractional representation is an exact *first-order* average of unit activities,
  taxon memberships, and conditional branch-length intensity, but it does not retain
  enough cross-level co-occurrence information to make the underlying tree mixture
  unique.
- For selecting an integral tree, non-destructive candidate generation followed by
  exact alignment-likelihood rescoring is safer than interpreting peel proportions.
- At 50--100 taxa, taxon-aware leaves-up construction is materially stronger than
  activity-first root-down ranking; a small local beam cures several greedy failures
  without NNI/SPR traversal.

## Small synthetic benchmark

Known mixtures were passed through the same limited-copy split-unit structure used by
this project, then optionally perturbed by 3% or 8% while preserving fixed-flow totals.

For single-tree inputs, the true tree was rank 1 in every 6/8-taxon case through 8%
perturbation.  For 45/40/15 near-tree mixtures the dominant true tree was always in the
top 10 and was usually in the top 3.

Destructive peeling often extracted 100% root mass but did not exactly recover the
original mixture.  On zero-noise near-tree mixtures its mean exact true-topology mass
was 94.4%; soft 50% peeling raised that to 98.1%.  Strict internal taxon capacities can
reconstruct all first-order marginals exactly, yet a different component mixture may
remain indistinguishable.

## Non-identifiability example

A 6-taxon zero-noise average was generated from:

- 0.45 `((0,((2,5),3)),(1,4))`
- 0.40 `((0,(2,(3,5))),(1,4))`
- 0.15 `(0,((1,4),((2,5),3)))`

A strict decoder recovered an alternative exact-marginal decomposition:

- 0.60 `((0,((2,5),3)),(1,4))`
- 0.25 `((0,(2,(3,5))),(1,4))`
- 0.15 `(0,((1,4),(2,(3,5))))`

Clade-support and internal taxon-membership RMSE were ~1e-17, yet topology-distribution
TV distance was 0.30.  Nested-clade co-occurrence differed by as much as 0.15.  This is
a structural information loss, not a search failure.

## Large-tree targeted tests

A taxon-aware leaves-up decoder with per-port top-match lists scales much better than
the original all-pairs implementation.  In targeted perturbed mixtures:

| Case | beam 1 | beam 4 | beam 8 |
|---|---|---|---|
| 100 taxa, far mixture, 3% noise | RF 0.357, wrong | exact, ~3.5 s | exact, ~7 s |
| 50 taxa, far mixture, 8% noise | RF 0.188, wrong | exact, ~0.4 s | exact, ~0.8 s |
| 20 taxa, near mixture, 3% noise | RF 0.056, wrong | same local optimum | exact, ~0.05 s |

These are targeted stress cases rather than a large Monte-Carlo claim, but they show
that bounded leaves-up backtracking can repair local mistakes without searching the
whole rooted-tree graph.

## Recommended workflow

1. Solve the LP with variable lower-arm branch lengths.
2. Convert the solution to the hierarchical fractional-tree representation.
3. Generate 5--20 integral candidates non-destructively from both root-down and
   taxon-aware leaves-up search.
4. Optimize branch lengths and evaluate the true alignment likelihood on every
   candidate; choose/rank by that likelihood.
5. Only if a mixture description is desired, run residual peeling separately and
   report its unexplained root mass plus a warning that component proportions may be
   non-unique.
