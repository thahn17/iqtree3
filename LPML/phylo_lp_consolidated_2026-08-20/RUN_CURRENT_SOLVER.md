# Running the current phylogenetic LP code

## Environment

The tested environment is Python 3.12 with the versions in `requirements.txt`.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

SciPy's `linprog(method="highs")` supplies the HiGHS dual-simplex LP solver; no separate solver installation is required.

## Generate a fixed-flow unit library

Generate a library whose taxon count matches the intended run. For example, for 20 taxa:

```bash
python pipeline_v3/fixed_flow_unit_generator_variable_lengths.py \
  --leaf-count 20 \
  --out-dir run/lib20 \
  --variable-length-paths \
  --variable-length-scope branches
```

This creates `run/lib20/library.json` and the supporting CSV files.

## Current slightly-faster profile

The slightly-faster profile is the two-point relaxation with bottleneck decoding and exact final tree rescoring. The research benchmark below simulates a matched alignment and runs the actual unit LP:

```bash
python research/code/benchmark_actual_100taxa_lp.py \
  --library run/lib20/library.json \
  --n 20 \
  --sites 1000 \
  --regime moderate \
  --mode two-point \
  --fingerprint-columns 4 \
  --saturated-threshold 2 \
  --root-planes 4 \
  --log-floor 1e-9 \
  --solver-method highs \
  --top-k 100 \
  --out-prefix run/two_point
```

## Current standard profile components

The standard profile unions projected and physical two-point/spectral candidate views, then performs exact full-alignment likelihood rescoring. Run the projected views first:

```bash
python research/code/benchmark_actual_100taxa_lp.py \
  --library run/lib20/library.json --n 20 --sites 1000 --regime moderate \
  --mode two-point --fingerprint-columns 4 --saturated-threshold 2 \
  --root-planes 4 --log-floor 1e-9 --projected-routing \
  --solver-method highs --top-k 100 --out-prefix run/projected_two_point

python research/code/benchmark_actual_100taxa_lp.py \
  --library run/lib20/library.json --n 20 --sites 1000 --regime moderate \
  --mode spectral --spectral-support-points 5 --fingerprint-columns 4 \
  --saturated-threshold 2 --root-planes 4 --log-floor 1e-9 \
  --projected-routing --solver-method highs --top-k 100 \
  --out-prefix run/projected_spectral
```

If either projected view returns fewer than 10--20 complete candidates, rerun that command without `--projected-routing` and retain candidates from both runs. Candidate CSVs are written beside each output prefix.

`benchmark_actual_100taxa_lp.py` is a controlled simulation harness. The self-contained FASTA application remains `pipeline_v3/phylo_lp_runner.py`; see `docs/FASTA_GTR_USAGE.txt`. The saturated projected/physical standard profile is still research code rather than a single production FASTA wrapper.

## Numerical settings

- Keep `--solver-method highs` for production runs.
- Use `--log-floor 1e-9`, matching HiGHS's small-matrix cutoff.
- Do not use IPM as the sole candidate-generating solve.
- Exact final likelihood rescoring, rather than the LP objective, selects the returned tree.
