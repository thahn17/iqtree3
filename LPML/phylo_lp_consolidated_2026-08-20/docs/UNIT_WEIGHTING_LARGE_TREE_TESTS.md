# Unit weighting tests on larger trees

Tests used simulated GTR+F alignments at 20 and 30 taxa, primarily 5,000--20,000 bp, with 180--300 nearby/random candidate topologies per replicate. Unit composition scores use child base fractions, so the base local score is already normalized for subtree size.

## Main findings

1. **Raw strict bottleneck remains the only individually robust default.** No tested size/balance normalization dominated it on every replicate.
2. A mild alternative bottleneck based on unit deficiency with
   `w = m^-0.2 * B^-0.2` or `m^-0.3 * B^-0.3`, where `B=4ab/m^2`, often improves the known-tree rank and RF distance, especially on difficult 20-taxon cases, but can worsen individual 30-taxon cases.
3. Strong flow standardization, per-flow z-scores, `1/m`, and low-tail/CVaR replacements were less stable than strict bottleneck.
4. Binary-entropy balance weights were not competitive with the simple power factor `B^gamma`.
5. A whole-alignment pairwise silhouette-like unit score,
   `cross(A,B) - max(within(A), within(B))`, is a useful complementary view and is automatically normalized for child sizes. At 10k bp across six 20/30-taxon tests, its mean known-tree rank was 10.0 versus 18.5 for the overlap bottleneck, but it was not better on every replicate.
6. Longer alignments make all of these statistics much less noisy, but do not select one universally optimal normalization. On individual 20k-bp tests, the preferred view differed between 20 and 30 taxa.

## Baseline-preserving policy

Do **not** replace the original bottleneck list. Preserve all top-K candidates from the raw/no-normalization bottleneck and only add candidates from complementary rankings. This guarantees that exact branch-length likelihood rescoring cannot be worse than the baseline candidate result.

The current best complementary views are:

- mild size/balance bottleneck: `w = m^-0.2 B^-0.2` and/or `m^-0.3 B^-0.3`;
- pairwise silhouette bottleneck: `cross - max(within-left, within-right)`.

For a raw top-K of 40, adding the top 20 candidates from each complementary view increased the candidate union from 40 to about 55.5 on average in six 20/30-taxon 10k-bp tests. Known-tree recall increased from 4/6 to 5/6 and mean minimum RF in the union improved from 0.67 to 0.33. From the observed structural ranks, allowing up to 40 additions from each complementary view would have included the known tree in all six cases while still preserving every raw-bottleneck candidate.

## Relation to the old two-point/no-spectral baseline

The large-tree weighting benchmark isolates the scoring behavior and does not by itself prove that spectral+weighted recovery has reached the old two-point/no-spectral likelihood baseline. Until that is demonstrated broadly, the strict consistency guarantee requires retaining the old two-point bottleneck candidate list as an auxiliary fallback. On the saved 8--12 taxon paired timing tests, producing that fallback costs about 0.49x the spectral LP build+solve time on average (range about 0.44--0.60x).

The target is to remove that fallback once the local coupling/weighting scheme matches it consistently on 20+ taxa.
