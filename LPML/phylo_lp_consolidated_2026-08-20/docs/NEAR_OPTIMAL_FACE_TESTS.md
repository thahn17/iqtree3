# Near-optimal-face / objective-profile construction tests

## Conclusion
The production root-only LP contains useful objective-profile information at the root, but most internal split-unit directions are deliberately projected out of the likelihood objective.  Therefore deep objective profiling is not a useful tree-construction metric unless internal likelihood coupling is restored.  Do not scan/profile a large fraction of units.

## Exact internal profile flatness (8 taxa, spectral GTR, 60 sites)
For every non-root split unit we forced `act[u] = 1` and re-solved the likelihood LP.

- flow 2: 4/4 units had objective drop <= 1e-7
- flow 3: 2/2
- flow 4: 4/4
- flow 5: 2/2
- flow 6: 3/3
- flow 7: 0/3 (drop about 3.4064)
- total: 15/18 internal units could be made fully active at essentially zero objective loss.

This is not a weakness of the profiler: those coordinates are genuinely flat because the aggressively pruned model profiles most site likelihood information to the root.

## Reduced-cost test
At one 8-taxon optimum, all four root activity variables had zero HiGHS bound marginals/reduced costs, while exact forced-root objective drops were approximately 4.3743, 0, 0.5372, and 1.8657.  Degeneracy therefore makes single-basis reduced costs unusable as substitutes for profile solves.

## Root-only profiles
Root profiling requires only floor(L/2) auxiliary solves.  In small tests:
- 8 taxa: 4 root profiles took about 0.53 s versus a 0.82 s base solve; candidate likelihood did not improve in that replicate.
- 12 taxa: 6 root profiles took about 1.43 s versus a 3.04 s base solve; best rescored candidate improved about -207.29 -> -204.99, still below the generating topology (-185.40).

Root profiles are therefore useful as a soft root prior, but not sufficient for internal topology recovery.

## Near-optimal-face extreme-point sketch
Add `objective >= z* - Delta`, optimize 4--32 diverse topology directions, and summarize the returned extreme points.

Barycentric averaging was bad: on an 8-taxon case the original best rescored candidate was about -371.42 while face averages gave about -402 to -406 despite making the previously-zero true root type nonzero.

Upper/quantile summaries were somewhat better but inconsistent:
- 8-taxon case A: single vertex best -508.86; 90th-percentile / max support sketch about -488.46; generating tree -465.57.
- 8-taxon case B: single -504.67; 4 probes -480.14; 24 probes -473.82; 32 probes worsened to -488.76; generating tree -460.11.
- 12-taxon case: single -219.50; 90th-percentile sketch -223.11; max sketch -225.27; generating tree -185.12.

The method is therefore not stable enough for production construction, and probe counts comparable to many units are explicitly not justified.

## One-solve central face point
Inside `objective >= z* - Delta`, add auxiliary `s_u <= act_u`, `s_u <= 1-act_u` and maximize normalized sum of `s_u`.

It sometimes improved candidate quality (8 taxa -454.49 -> -442.84; 12 taxa -257.62 -> -254.75) but a Delta grid was inconsistent; adding root taxon-membership centering often worsened results and some settings yielded no complete hard-support candidates.

## Polytope interpretation
The near-optimal polytope already has an exact H-representation: the original LP constraints plus one objective-floor inequality.  Explicit vertex/facet enumeration is generally exponential and is not returned by SciPy's HiGHS wrapper.  A support-function sketch from auxiliary objectives is a compact approximation, but the experiments above show that the current internal face is too flat for that sketch to identify integral trees reliably.

## Recommended use
1. Keep the ordinary fractional solution / G2 hierarchy for internal construction.
2. Optionally profile the O(L) root split types and use profile loss as a soft root prior instead of treating zero activity in one extreme point as impossible.
3. Do not objective-profile internal physical units in the current root-only model.
4. If deep objective-profile construction is desired later, first retain a small fixed amount of internal likelihood information in the main LP; otherwise the missing information is not present for a second pass to recover.
