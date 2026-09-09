# Fractional-tree composition and integral-tree identification

## Main conclusion

Do not use greedy residual-mass removal to identify the best integral topology. It is a packing decomposition, not a likelihood ranking, and shared support between near-identical trees can be consumed by the first selected tree.

The useful workflow is non-destructive:

- preserve variable lower-arm branch lengths in the LP;
- normalize them away from topology frequency when matching (`frequency=a_u`, `conditional length=ell/a_u`);
- compose the fractional tree as a fixed-order local hierarchy;
- produce a top-K list without subtracting anything;
- rank conservatively by bottleneck split support;
- re-score all K trees with the actual alignment likelihood and fast integral branch-length optimization.

## Composition strategies tested

### G0: first-order split grammar
Uses `P(split | current clade)`.

### G1: parent-conditioned grammar
Uses `P(child split | parent split, current clade)`.

### G2: local-fragment grammar (preferred)
Keeps the two immediate child decisions jointly for each current split.  This retains one fixed layer of split co-occurrence without a context depth that grows with tree height.

G2 consistently reduced hybrid recombination relative to G0.  Example mean original-component probability mass:

- 100 taxa, 50 bp, 5% contamination: 0.571 (G0) -> 0.648 (G2)
- 100 taxa, 50 bp, 15% contamination: 0.326 -> 0.384
- 50 taxa, 50 bp, 15% contamination: 0.417 -> 0.502

This does not eliminate hybrid mass; deeper cross-level co-occurrence is genuinely absent from first-order LP marginals.  The important point is that a fixed local fragment improves composition without requiring depth to grow with tree height.

## Ranking metric

Product/log support is useful for ordinary cases but can over-penalize one weak local decision and allow many moderately supported hybrid decisions to outrank the ML tree.  For candidate selection the safer ranking found here is lexicographic:

1. maximize the minimum supported split (bottleneck);
2. break ties by total log split support;
3. finally evaluate true alignment likelihood.

On the hardest earlier 50-taxon/20-bp cases, G2 product ranks of 139 and 363 became bottleneck ranks 8 and 15.

## Alignment-backed best-known-ML benchmark

Alignments were simulated under JC69 with a common known branch length.  For each replicate an offline multi-start search exhaustively scored every NNI neighbor at each hill-climbing step; the best visited topology was used as the target.  The fractional average used the high-likelihood visited trees plus a deliberate 15% structural-contamination component, making this harsher than a pure likelihood average.

Observed maximum target ranks under G2 + bottleneck:

| taxa | bp | cases | maximum rank |
|---:|---:|---:|---:|
| 20 | 20 | 1 | 5 |
| 50 | 20 | 1 | 3 |
| 100 | 20 | 1 | 27 |
| 20 | 50 | 5 | 5 |
| 50 | 50 | 3 | 10 |
| 100 | 50 | 3 | 30 |
| 20 | 200 | 5 | 1 |
| 50 | 200 | 3 | 1 |
| 100 | 200 | 1 | 1 |
| 20 | 1000 | 5 | 1 |
| 50 | 1000 | 3 | 1 |
| 100 | 1000 | 1 | 1 |

At 50 bp the best-known topology often differed from the generating topology; at 200--1000 bp the generating topology usually became the best-known tree and the hierarchy collapsed to rank 1.

The data therefore support a taxa-dependent safety list in weak-signal cases, but **not** a list size that grows with alignment length.  Once enough sites are present, more sites make the hierarchy easier, not harder.

## Practical top-K rule

The empirical 50-bp maxima are covered by `ceil(n_taxa/3)`.  A safer uncertainty-aware rule is

`K = max(3, ceil(n_taxa/3), ceil(N_eff))`

when the hierarchy is diffuse, and `K=3` when it is concentrated.  `N_eff` should be estimated from the normalized non-destructive hierarchy/candidate weights, not from peeled proportions.

For 100 taxa this suggests scoring roughly 30--45 integral trees in an uncertain alignment, then choosing by the actual optimized alignment likelihood.

This is an empirical safety rule, not a formal probability guarantee.
