# 100-taxon unit-weighting and spectral-secondary experiments

## Main conclusions

1. Raw bottleneck remains the safest single ranking at 100 taxa.
2. A mild coupled size/balance view is useful as a complementary ranking:
   `w=(m B)^(-q)=(4ab/m)^(-q)`, with `q` around 0.1--0.15 for the 100-taxon moderate tests. Stronger `q=0.2` had a rank-40 failure in one replicate.
3. Balance-only weighting is unsafe at 100 taxa. `B^-0.2` put the generating topology near rank 100--116 and sometimes selected an RF-maximal tree.
4. Size-only weighting is also much less stable than raw bottleneck.
5. For 100 taxa / 10k bp moderate GTR+F simulations, raw bottleneck put the generating topology at ranks 3, 5, 3 across three replicates. At 20k bp it was rank 2; at 30k bp it was rank 3 in the additional run.
6. Under a deliberately extreme wide-branch 100-taxon / 10k-bp simulation, raw rank was 41, while the coupled `q=0.2`/`0.3` views moved it to rank 25. The top pairwise-score tree was nevertheless poor (RF near maximum), so these views must only add candidates and must be followed by real likelihood rescoring.
7. Whole-alignment pairwise preprocessing is cheap: about 0.38--0.47 s for 100x10k, 0.83 s for 100x20k, and 1.15 s for 100x30k in the current Python benchmark.
8. A GTR-centered pair covariance did not materially outperform simple Hamming separation in the extreme saturation test. It retained similar ranks, so the difficulty is lack of topology signal after saturation rather than the particular distance normalization.

## Secondary spectral-LP coupling screening

The secondary LP score uses child base fractions and a local overlap score. The best simple structural weighting tested was mild and did not require the pairwise alignment modifier inside the LP objective.

For one 8-taxon moderate case:

- unweighted secondary, best alpha blend: about -299.11;
- `m^-0.1 B^-0.1`, no pairwise modifier: about -293.11;
- same structural weight plus pairwise modifier: about -294.09.

Thus the whole-alignment pairwise signal appears more useful as a separate candidate-ranking view than as a multiplier inside the secondary LP objective.

The secondary optimum should not replace the primary spectral solution. Candidate extraction should use at least alpha in `{0, 0.5, 1}` and union the resulting trees. On completed 8-taxon screening pairs, the best candidate from that union met or exceeded the two-point/no-spectral candidate under the same (fast, one-start) final branch-optimization settings. This is screening evidence only; production comparisons should retain multi-start branch optimization.

## 100-taxon candidate policy suggested by the tests

Use several non-destructive views and union them before final branch optimization:

1. mandatory raw bottleneck top K;
2. mild effective-child-size bottleneck using `q=0.1` (optionally `q=0.15`);
3. whole-alignment pairwise silhouette top K_add;
4. for heavily diffuse/saturated solutions, optionally add `q=0.2` as a diversification view, never as the sole ranking.

The final winner is always chosen by optimized GTR+F likelihood.

## Scaling of the secondary LP score at 100 taxa

With 4,012 split units and A fixed alignment-score patterns, the overlap score layer adds `4 A U` variables and `8 A U` rows. At A=8 this is 128,384 variables and 256,768 rows. Structural weighting changes only objective coefficients and adds no variables or rows.

The whole-alignment pairwise summary is O(L^2 S) preprocessing and O(L^2) stored data, independent of the number of unique site patterns.

## Full 100-taxon LP status

The current Python research constructor is not yet fast enough for a full 100-taxon paired spectral/two-point solve in this interactive benchmark. A 20-taxon structural paired test also crossed a 120-second cap. Therefore the 100-taxon results here test the unit-weight/ranking layer directly on long alignments; they do not yet constitute a full 100-taxon LP end-to-end validation.
