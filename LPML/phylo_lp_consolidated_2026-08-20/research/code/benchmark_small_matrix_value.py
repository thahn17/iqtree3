#!/usr/bin/env python3
"""Test HiGHS coefficient cutoff/tolerance effects on phylogeny LP decoding."""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path[:0] = [str(HERE / "v3_experiments"), str(ROOT / "pipeline_v3"), str(HERE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import lp_exact_speedups as exact  # noqa: E402
import multisite_saturated_100taxa as saturated  # noqa: E402


def decode(model, result, top_k):
    candidates, _ = saturated.rank_fingerprint_candidates(
        model,
        result,
        topk=top_k,
        root_beam=192,
        leaf_beam=48,
        child_branching=16,
        leaf_branching=10,
        merge_branching=128,
        max_expansions=800000,
    )
    return tuple(bench.key(candidate.topology) for candidate in candidates)


def run(args):
    rng = random.Random(972000 + args.n + args.rep)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    tree = bench.random_tree(args.n, rng)
    _, columns, _ = bench.simulate(tree, args.sites, simulation_model, rng, args.regime)
    fitted = bench.empirical_model(columns)
    bench.configure(fitted)
    model = saturated.FingerprintSaturatedModel(
        args.library,
        columns,
        fingerprint_columns=args.fingerprint_columns,
        anchor_weighting="uniform",
        substitution_model=fitted,
        branch_relaxation=args.mode,
        spectral_support_points=5,
        saturation_threshold=2,
        root_planes=4,
    )
    compiled = exact.CompiledSparseLP(model.lp)
    variants = [
        ("default", 1e-8, None),
        ("keep_1e-12", 1e-8, {"small_matrix_value": 1e-12}),
        ("keep_1e-12_tol_1e-9", 1e-9, {"small_matrix_value": 1e-12}),
        ("keep_1e-12_tol_1e-10", 1e-10, {"small_matrix_value": 1e-12}),
    ]
    rows = []
    base_value = None
    base_x = None
    base_keys = None
    for name, tolerance, options in variants:
        value, result, timing = compiled.solve(
            model.obj,
            feasibility_tolerance=tolerance,
            extra_options=options,
        )
        keys = decode(model, result, args.top_k)
        if base_value is None:
            base_value, base_x, base_keys = value, result.x.copy(), keys
        base_set, key_set = set(base_keys), set(keys)
        union = base_set | key_set
        rows.append(
            {
                "variant": name,
                "solve_s": timing.solve_s,
                "objective": value,
                "objective_abs_diff": abs(value - base_value),
                "x_inf_diff": float(np.max(np.abs(result.x - base_x))),
                "positive_below_1e-9": int(np.count_nonzero((result.x > 0) & (result.x < 1e-9))),
                "positive_below_1e-12": int(np.count_nonzero((result.x > 0) & (result.x < 1e-12))),
                "candidates": len(keys),
                "candidate_jaccard": 1.0 if not union else len(base_set & key_set) / len(union),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--sites", type=int, default=500)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--rep", type=int, default=0)
    parser.add_argument("--mode", choices=["two-point", "spectral"], required=True)
    parser.add_argument("--fingerprint-columns", type=int, default=4)
    parser.add_argument("--top-k", type=int, default=50)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows = run(args)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
