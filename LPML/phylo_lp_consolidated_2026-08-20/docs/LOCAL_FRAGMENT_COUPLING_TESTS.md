# Local-fragment coupling / unit-score experiments

## Construction tested
For selected anchor site patterns `s` and every split unit `u` with child flows `a,b`, introduce a normalized local-fragment viability score `v[s,u]`. Its piecewise-linear upper score depends jointly on:

- left child A/C/G/T composition from the shared named-taxon flow,
- right child A/C/G/T composition,
- the shared left branch-effect coordinates,
- the shared right branch-effect coordinates.

Thus the score is a fixed-depth local fragment (G2-like) quantity rather than a marginal activity score.

First solve the ordinary likelihood LP, obtaining raw LP log-objective `z*`. Then solve a second LP

    maximize sum_s w_s sum_u v[s,u]
    subject to all original constraints
               primary_log_objective >= z* - Delta.

All unit scores are therefore optimized on one shared fractional topology. Delta=0 is a lexicographic tie-break on the optimal face; positive Delta permits a controlled trade of fractional likelihood for better local topology coherence.

## Main exact local-score results
All final numbers below are continuous-branch optimized integral-tree likelihoods from the non-destructive bottleneck candidate list.

| case | spectral baseline | spectral local score, Delta=2 | two-point local score, Delta=2 | generating topology |
|---|---:|---:|---:|---:|
| 8 taxa, moderate, rep 0 | -310.17 (another equivalent vertex gave -304.13) | **-299.98** | -303.76 | -291.78 |
| 8 taxa, moderate, rep 1 | -302.64 | **-284.82** | -294.92 | -282.18 |
| 8 taxa, wide, rep 0 | -295.61 | **-294.92** | -295.41 | -295.00 |
| 12 taxa, moderate, rep 0 | -320.31 | **-305.11** | -326.26 | -270.75 |

The 12-taxon case remains far from the reference topology, so the coupling is a substantial improvement rather than a complete solution.

## Delta behavior
Delta=0 was inconsistent: it improved some moderate cases substantially but worsened the wide 8-taxon case and slightly worsened the 12-taxon moderate case. Delta=2 was much more robust in the completed tests. With 30-40 columns this corresponds to roughly 0.05-0.067 LP log units per column; a production parameter should be normalized by alignment size / known relaxation error rather than kept as an absolute 2.

## Raw unit-score objectives
Maximizing `sum_u act[u]` is essentially useless: an integral rooted binary L-taxon tree uses L-1 split units, so the sum is effectively fixed. Weighting each unit by its activity in the first LP vertex also selected worse tree basins in the tested moderate case. Unit scores need alignment-conditioned local information.

## Scaling tests / rejected compromises
- Score only low-flow units (m<=3 or m<=4): failed; the 12-taxon moderate case worsened badly and one m<=4 run produced no complete hard-support candidate.
- Random fixed-size local-fragment clouds: cheap but unstable across replicates.
- Score only physical GTR branch points instead of the spectral branch polytope: much faster but materially weaker on moderate cases.
- Exact score clouds: strongest, but spectral preprocessing currently enumerates too many composition/branch combinations.

## Fixed-order saturated child-composition test
A local score using child base counts truncated to 0/1 (T=1) or 0/1/>=2 (T=2) was tested. The LP representation uses perspective convex-hull variables for the saturated counts, so its state dimension does not grow with tree height.

8-taxon results with Delta=2:

| case | T=1 | T=2 | baseline |
|---|---:|---:|---:|
| moderate rep 0 | **-299.01** | -299.98 | about -304 to -310 |
| moderate rep 1 | **-298.25** | -302.44 | -302.64 |
| wide rep 0 | **-295.15** | -302.53 | -295.61 |

T=1 was more robust in this very small sample; T=2 matched the strongest exact-count result on moderate rep 0 but was less stable elsewhere. The prototype still used exact composition enumeration to construct the saturated score table, so its measured preprocessing time is *not* the intended scaling. A real saturated DP would make the score-table state count fixed in L.

## Current recommendation
Do not add raw unit-activity maximization. Keep the primary spectral likelihood LP unchanged, then use an optional second-pass local-fragment objective with a small likelihood slack. The strongest evidence supports an all-flow local score; the most promising scalable approximation is a fixed-order saturated child-composition score, probably starting with T=1 and testing T=2 as an accuracy option. This should be implemented with the saturated DP directly before it is promoted to the production default.
