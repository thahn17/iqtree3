# Consolidated phylogenetic LP research package

This folder consolidates the current runnable pipeline snapshot and the research code/results used to develop the fractional phylogenetic LP and integral-tree recovery strategy.

## Start here

- `RUN_CURRENT_SOLVER.md` — tested environment, library generation, and current fast/standard launch commands.
- `requirements.txt` — tested Python dependencies.
- `docs/INTEGRAL_TREE_STRATEGIES_AND_COMPARISON.md` — consolidated strategy/evidence document and the planned definitive two-point-vs-spectral comparison.
- `pipeline_v3/` — last explicitly self-contained end-to-end pipeline bundle, including FASTA input, JC/GTR+F support, spectral non-JC branch relaxation, top-K recovery, continuous branch optimization, and timing/output code.
- `research/code/` — later experimental scripts and working code snapshots. These are included for reproducibility/history and are not all production-clean.
- `research/results/` — CSV outputs from the experiments, including 20–30 and 100-taxon weighting studies.
- `docs/` — detailed research notes and benchmark reports.
- `source_archives/` — original self-contained v3 archive.

## Important status note

The research after v3 explored several alternative topology-recovery objectives and weightings. Some newest top-level scratch files referenced helper modules from the self-contained v3 snapshot or were written as benchmark harnesses rather than as a polished runnable application. For that reason, `pipeline_v3/` is kept intact as the known self-contained baseline, while later work is kept under `research/`.

## Provisional algorithmic status

- Raw bottleneck remains the safest single integral-tree ranking.
- Spectral GTR is the more complete branch-effect relaxation but, without topology coupling, often gives weaker hard-support tree recovery than the old two-point GTR relaxation.
- Small exact-local-fragment experiments show that coupling can reverse this and make spectral outperform two-point.
- At 100 taxa, raw bottleneck and mild coupled size/balance views have been tested extensively as ranking heuristics; the definitive full spectral-vs-two-point LP comparison remains future work.

See the strategy document for the data and proposed no-regression comparison protocol.
