# Local score approximation / bottleneck scaling experiments

## Main findings

1. The exact enumeration-based local-fragment upper score is not well reproduced by the simplest fixed-state saturated DP. A T=1 fixed-state DP was cheap, but did not reliably reproduce the tree-selection effect of the enumerated score.

2. A much cheaper alignment-conditioned overlap score requires no composition enumeration. For a unit with left/right normalized base-count vectors f_L,f_R and a reference transition matrix Pbar,

   v(u,s) = sum_i pi_i min((Pbar f_L)_i, (Pbar f_R)_i).

   It uses four auxiliary variables per unit / selected alignment pattern and two inequalities per variable. Score construction was around 0.02 s at 8 taxa and 0.06 s at 12 taxa in the tested models.

3. This overlap score is a useful replacement *heuristic*, not a numerically faithful approximation of the old optimistic local upper. On split-normalized exact small-tree clouds, overlap alone gave R^2 about 0.64 / 0.40 and Spearman about 0.74 / 0.57 at 8 / 12 taxa. Adding presence information improved this to about 0.76 / 0.50 and 0.89 / 0.64. Richer simple features can reach about R^2 0.95 at 8 taxa, but only about 0.69 at 12 taxa and require nonlinear concentration / max / L1 features.

4. Fully maximizing a secondary unit score is not always best for bottleneck decoding. For the same primary and secondary LP vertices, decoding the midpoint x_0.5 = 0.5 x_primary + 0.5 x_secondary sometimes produced better final integral trees than decoding either endpoint. In a wide-branch 8-taxon case, fully maximizing the overlap score produced a very poor candidate while the primary/midpoint retained the near-optimal candidate. In a 12-taxon short-alignment case the fully maximized point was better. Therefore no single alpha is robust.

   Recommended low-cost strategy: form candidates from alpha in {0, 0.5, 1}, union the candidate topologies, and perform the usual exact fixed-tree branch optimization / likelihood rescoring. This needs only the original primary solve plus one secondary solve; the extra alpha values are convex combinations and require no additional LP optimization.

5. On larger, known-tree simulations (20 and 30 taxa) with 1k--30k bp, child-flow-normalized local scores were compared on hundreds of NNI-nearby / random candidate trees. Raw bottleneck of the per-unit local score was the best overall ranking statistic. Additional per-flow z normalization was consistently worse.

   Same-tree / same-candidate prefix experiment:

   - 20 taxa: true-tree bottleneck rank 9 at 1,000 bp, 2 at 5,000 bp, 2 at 10,000 bp, 2 at 30,000 bp.
   - 30 taxa: true-tree bottleneck rank 1 at 1,000, 5,000, 10,000, and 30,000 bp.

   Among summed unit scores, mild weightings were preferable to aggressive ones. 1/sqrt(flow) was strongest in the 20-taxon prefix case; sqrt(split-balance) was strongest in the 30-taxon prefix case. Combining both penalties was worse. Equal weighting remained competitive. 1/flow was inconsistent and should not be the default.

6. On five independent 5,000-bp 20/30-taxon candidate pools, raw bottleneck had the best average true-tree rank (mean 21, median 13) among tested metrics. Its top-ranked tree was on average only about 2.4 rooted clade differences from the known tree, even when the exact known topology was not ranked first.

## Recommended experimental design

- Keep the primary likelihood LP unchanged.
- Add an optional secondary local-score pass using the enumeration-free overlap score.
- Normalize child counts by child flow before scoring; this removes the dominant tree-size scaling.
- Use equal or at most 1/sqrt(flow) unit weights. Avoid per-flow z normalization and avoid strong 1/flow weighting as defaults.
- Do not trust only the fully maximized secondary extreme point. Generate bottleneck candidates from alpha = 0, 0.5, 1 between primary and secondary solutions and union them before exact final likelihood rescoring.
- Longer alignments improve/stabilize the bottleneck score without requiring the score state dimension to grow with alignment length.
- The overlap score should be described as a scalable unit-usefulness surrogate, not as an accurate numeric approximation to the old enumerated optimistic local upper.
