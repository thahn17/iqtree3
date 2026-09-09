# Moderate GTR spectral counterexample diagnostic

Alignment: `spectral_gtr_8_moderate.fa`, 8 taxa, 80 sites, GTR+F rates AC/AG/AT/CG/CT/GT = 1,2.3,0.6,0.85,1.5,0.95.

## Main results

- The production run used K=20. Increasing requested K to 50,100,200,500,1000 did not recover the known better topology. The hard-support decoder generated only 24 distinct positive-support trees; the better topology was absent.
- A likelihood NNI climb from the generating topology found a better-known topology with optimized logL -597.6503942: `((((0,4),3),7),((1,6),(2,5)))` in zero-based taxon indices. The generating topology optimized to -597.8761060.
- The best spectral top-20 candidate had logL -670.9083095.
- The better-known topology is not simply rank 21+: it has a 4+4 root while the single unrestricted spectral LP optimum assigns root activities (1+7,2+6,3+5,4+4) = (0,1,0,0).
- Root-profile solves show spectral is actually much kinder to the correct root than two-point: forcing 4+4 loses only 4.4529 LP objective units under spectral, versus 34.9737 under two-point. Spectral root-profile drops are approximately: 1+7 1.1948, 2+6 0, 3+5 4.0452, 4+4 4.4529.
- Therefore raw root activity from one LP extreme point is not a calibrated probability. A soft profile score based on objective loss is more appropriate for tree identification.

## Branch-range tests

Hard-bounding the spectral physical branch curve does not solve the miss. Full spectral makes the 4+4 root *more* competitive than bounded variants. At [0,0.5], forcing 4+4 costs ~19.50; at [0,0.35] it costs ~24.94; at [0.04,0.35] ~24.97. Narrow ranges [0.05,0.25], [0.08,0.25], and [0.10,0.20] still put all root mass on 1+7.

The best-known topology's optimized lengths are mostly moderate: after collapsing the artificial degree-2 root into one unrooted edge, branch CV is ~0.522, positive-edge CV ~0.417, positive log-SD ~0.530, one near-zero edge, and max/median ~1.68. The best spectral hierarchy candidate has CV ~0.989, positive-edge CV ~0.607, positive log-SD ~0.782, four near-zero unrooted edges, and max/median ~3.73.

## Post-processing normalization

1. Root/unit activities are marginals of one optimal LP vertex, not posterior probabilities. Zero activity can mean “not used by this extreme optimum,” not “structurally implausible.”
2. Physical copy activity is also a poor abstract-clade score when interchangeable unit copies split support. A copy-aggregated taxon-support experiment helps in other cases, but here even the root-4+4-constrained LP gives the best-known topology zero support on several internal splits; under that experimental aggregate metric it ranks about 4987/7875 among 4+4-root trees. The root-only pruning therefore leaves internal structure underidentified in this example.
3. Spectral q/a should not be interpreted as branch length on internal arms in the pruned model. Many internal mode variables are not likelihood-coupled and sit at arbitrary vertices. The final branch optimizer correctly ignores them as initialization. Root modes are likelihood-coupled and here correspond to near-zero effects, not huge lengths.
4. For reversible models the two branch lengths incident to the artificial degree-2 root are non-identifiable separately; only their sum is an unrooted branch length. Uniformity diagnostics should collapse them before computing dispersion.

## Search/refinement

Likelihood-monotone NNI refinement helps but does not fully cure the hierarchy miss. Starting from the best spectral candidate improves -670.91 -> -632.49; another top candidate reaches -620.33. Starting from the best candidate from a root-4+4-profiled LP reaches -621.27. These are still below -597.65, so the candidate hierarchy needs a softer support notion, not merely a larger K or a few NNI steps.

## Recommended next adjustment

Keep full spectral as non-JC default. For identification, replace hard root activity with a profile-loss normalization, e.g. root weight proportional to exp(-(z* - z_r*)/tau). In this case tau=2 gives approximate root weights (1+7,2+6,3+5,4+4) = (0.307,0.559,0.074,0.060), rather than (0,1,0,0). The same principle should be used cautiously lower in the hierarchy, preferably with cheap local/taxon support rather than thousands of auxiliary LP solves.

A branch-uniformity parameter is useful as an explicit optional search prior/diagnostic, but hard length bounds are not the fix. A scale-free choice is the expected log-length SD (or equivalently a dispersion strength) computed after collapsing the root edge. It should guide candidate expansion/refinement; final reported ML likelihood should remain unpenalized unless the user explicitly requests a branch-length prior.
