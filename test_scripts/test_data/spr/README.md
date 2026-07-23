# SPR test fixtures

This directory is reserved for deterministic, small inputs used while implementing
subtree-pruning-and-regrafting (SPR) topology search.  It is deliberately separate
from `test_scripts/test_data`, whose files are consumed by the cross-platform
runtime/regression drivers (`test_iqtree.sh` and `test_iqtree.ps1`).

## Current fixture

`six_taxa` is a six-taxon DNA alignment with three clearly separated close pairs.
The supplied start tree deliberately groups `A` with `C` and `B` with `D`; the
reference topology instead groups `A,B`, `C,D`, and `E,F`.

| File | Purpose |
| --- | --- |
| `six_taxa.phy` | Small sequential PHYLIP alignment. |
| `six_taxa.start.tree` | Intentionally suboptimal user start tree. |
| `six_taxa.reference.tree` | Reference topology for structural assertions; it is not an exact likelihood oracle. |
| `six_taxa.ab.constraint.tree` | Constraint requiring the `A,B` split, for a later constraint-compatibility test. |
| `manifest.tsv` | Machine-readable fixture metadata and intended future checks. |

## Conventions

- Inputs are immutable.  Tests must write all run products under a temporary or
  caller-provided output directory, never into this directory.
- A fixture must not be added to `test_iqtree.sh` or `test_iqtree.ps1` until the
  SPR CLI and its deterministic result contract are implemented.
- Tests that evaluate a temporary SPR candidate must assert topology, branch-length,
  and likelihood restoration after rollback.  A final accepted move must be checked
  by a fresh likelihood computation, not a screening score.
- The reference topology documents the intended local rearrangement only.  Exact
  likelihood values are intentionally omitted until the production SPR evaluator is
  implemented and validated across supported platforms.

## Manual baseline command

Once an SPR option is implemented, use a caller-owned output prefix, for example:

```text
iqtree3 -s test_scripts/test_data/spr/six_taxa.phy \
  -t test_scripts/test_data/spr/six_taxa.start.tree \
  --prefix <temporary-output>/six_taxa
```

Add the eventual SPR-specific options only after their names and semantics are
stabilized.  Keep the seed and thread count explicit in automated tests.
