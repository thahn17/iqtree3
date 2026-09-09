# Spectral non-JC branch relaxation test

`--branch-relaxation auto` now selects the 3-mode spectral relaxation for non-JC reversible models and the legacy one-mode/two-point representation for JC-like models.

The spectral branch curve is

`P(t) = 1*pi^T + sum_r exp(-mu_r t) B_r`.

The LP stores three homogenized mode activities per lower arm and uses an outer polyhedron whose support planes are valid for the full continuous `t in [0,infinity)` curve.  The default 7-support test GTR model produced 26 inequalities and 24 vertices in this 3-D polytope.  Direct reconstruction of `P(t)` matched `scipy.linalg.expm(Q*t)` to `1.39e-14` maximum absolute error over the validation grid.

## 8-taxon stress alignment (60 bp; simulated true branches 1e-4 through 10)

| branch LP | build s | solve s | LP score | best rescored tree |
|---|---:|---:|---:|---:|
| two-point | 1.71 | 1.94 | -396.5106 | -518.4262 (top 10, 4-start branch optimization) |
| spectral-7 | 8.61--9.17 | 2.89--3.07 | -163.9140 | -482.9365 (top 20, 4-start branch optimization) |

The known generating topology optimized to -544.2910, so the spectral hierarchy contained a substantially better likelihood topology in this replicate.

## 8-taxon moderate alignment (80 bp; simulated branches about 0.05--0.24)

| branch LP | build s | solve s | LP score | best rescored tree |
|---|---:|---:|---:|---:|
| two-point | 1.57 | 2.17 | -446.3896 | -639.3872 |
| spectral-7 | 8.31 | 3.66 | -203.7871 | -670.9083 |

The known generating topology optimized to -598.8922.  Thus full branch coverage can loosen the fractional topology enough to worsen structural candidate recovery on some moderate-length cases, even though it is the physically safer relaxation.

## Important bug found by the test

Using the LP spectral modes as fixed-tree optimizer initialization was wrong in the pruned/root-only model: most internal mode variables are weakly likelihood-coupled and can sit at saturation endpoints.  Starting L-BFGS-B there produced a false plateau around -643/-873.  The final tree scorer now ignores arbitrary internal spectral modes and uses deterministic finite multistarts by default.
