#!/usr/bin/env python3
"""Benchmark exact LP solve speedups on real 15--25 taxon unit libraries."""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PIPELINE = ROOT / "pipeline_v3"
V3_EXPERIMENTS = HERE / "v3_experiments"
sys.path[:0] = [str(V3_EXPERIMENTS), str(PIPELINE), str(HERE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import lp_exact_speedups as exact  # noqa: E402
import multisite_saturated_100taxa as saturated  # noqa: E402


def candidate_keys(model, result, top_k):
    candidates, _ = saturated.rank_fingerprint_candidates(
        model,
        result,
        topk=top_k,
        root_beam=128,
        leaf_beam=32,
        child_branching=12,
        leaf_branching=8,
        merge_branching=96,
        max_expansions=500000,
    )
    return tuple(bench.key(candidate.topology) for candidate in candidates)


def run(args):
    rng = random.Random(912000 + 1000 * args.n + args.rep)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    true_tree = bench.random_tree(args.n, rng)
    _, columns, _ = bench.simulate(true_tree, args.sites, simulation_model, rng, args.regime)
    fitted_model = bench.empirical_model(columns)
    bench.configure(fitted_model)

    started = time.perf_counter()
    model = saturated.FingerprintSaturatedModel(
        args.library,
        columns,
        fingerprint_columns=args.fingerprint_columns,
        fingerprint_offset=args.fingerprint_offset,
        anchor_weighting=args.anchor_weighting,
        substitution_model=fitted_model,
        branch_relaxation=args.mode,
        spectral_support_points=args.spectral_support_points,
        saturation_threshold=args.saturation_threshold,
        root_planes=args.root_planes,
    )
    build_s = time.perf_counter() - started
    objective_scale = max(float(sum(model.weights)), 1.0)

    rows = []
    baseline_value, baseline_result, baseline_internal_s = model.lp.solve(model.obj, True)
    baseline_keys = candidate_keys(model, baseline_result, args.top_k)
    rows.append(
        {
            "variant": "baseline_recompile",
            "method": "highs",
            "row_scale": 0,
            "objective_scale": 1.0,
            "edge_strategy": "default",
            "compile_s": np.nan,
            "solve_s": float(baseline_internal_s),
            "objective": float(baseline_value),
            "objective_abs_diff": 0.0,
            "x_inf_diff": 0.0,
            "candidate_jaccard": 1.0,
            "first_candidate_same": 1,
            "candidate_count": len(baseline_keys),
        }
    )

    direct = exact.DirectHighsLP(model.lp)
    direct_value, direct_result, direct_timing = direct.solve(model.obj)
    direct_keys = candidate_keys(model, direct_result, args.top_k)
    base_set, direct_set = set(baseline_keys), set(direct_keys)
    direct_union = len(base_set | direct_set)
    rows.append(
        {
            "variant": "direct_highs_cached_csc",
            "method": "highs",
            "row_scale": 0,
            "objective_scale": 1.0,
            "edge_strategy": "default",
            "compile_s": direct_timing.compile_s,
            "solve_s": direct_timing.solve_s,
            "objective": float(direct_value),
            "objective_abs_diff": abs(float(direct_value - baseline_value)),
            "x_inf_diff": float(np.max(np.abs(direct_result.x - baseline_result.x))),
            "candidate_jaccard": 1.0 if direct_union == 0 else len(base_set & direct_set) / direct_union,
            "first_candidate_same": int(
                bool(direct_keys and baseline_keys and direct_keys[0] == baseline_keys[0])
            ),
            "candidate_count": len(direct_keys),
        }
    )

    if args.profile == "full":
        variants = [
            ("compiled", False, 1.0, "highs", None),
            ("objective_normalized", False, objective_scale, "highs", None),
            ("row_and_objective_scaled", True, objective_scale, "highs", None),
            ("scaled_dual_simplex", True, objective_scale, "highs-ds", None),
            ("scaled_ipm", True, objective_scale, "highs-ipm", None),
            ("scaled_dantzig", True, objective_scale, "highs-ds", "dantzig"),
            ("scaled_devex", True, objective_scale, "highs-ds", "devex"),
        ]
    else:
        variants = [
            ("compiled", False, 1.0, "highs", None),
            ("unscaled_ipm", False, 1.0, "highs-ipm", None),
            ("unscaled_devex", False, 1.0, "highs-ds", "devex"),
            ("row_scaled_only", True, 1.0, "highs", None),
            ("row_and_objective_scaled", True, objective_scale, "highs", None),
        ]

    compiled_by_scale = {}
    for name, row_scale, scale, method, edge_strategy in variants:
        if row_scale not in compiled_by_scale:
            compiled_by_scale[row_scale] = exact.CompiledSparseLP(
                model.lp, row_scale=row_scale, deduplicate=args.deduplicate
            )
        compiled = compiled_by_scale[row_scale]
        value, result, timing = compiled.solve(
            model.obj,
            maximize=True,
            objective_scale=scale,
            method=method,
            feasibility_tolerance=args.tolerance,
            simplex_dual_edge_weight_strategy=edge_strategy,
        )
        keys = candidate_keys(model, result, args.top_k)
        base_set, key_set = set(baseline_keys), set(keys)
        union_size = len(base_set | key_set)
        rows.append(
            {
                "variant": name,
                "method": method,
                "row_scale": int(row_scale),
                "objective_scale": scale,
                "edge_strategy": edge_strategy or "default",
                "compile_s": timing.compile_s,
                "solve_s": timing.solve_s,
                "objective": float(value),
                "objective_abs_diff": abs(float(value - baseline_value)),
                "x_inf_diff": float(np.max(np.abs(result.x - baseline_result.x))),
                "candidate_jaccard": 1.0 if union_size == 0 else len(base_set & key_set) / union_size,
                "first_candidate_same": int(bool(keys and baseline_keys and keys[0] == baseline_keys[0])),
                "candidate_count": len(keys),
            }
        )

    metadata = {
        "n": args.n,
        "sites": args.sites,
        "regime": args.regime,
        "rep": args.rep,
        "mode": args.mode,
        "fingerprint_columns": args.fingerprint_columns,
        "anchor_weighting": args.anchor_weighting,
        "split_units": len(model.splits),
        "variables": len(model.lp.names),
        "inequalities": len(model.lp.le),
        "equalities": len(model.lp.eq),
        "unique_patterns_seen": len(set(columns)),
        "patterns_in_lp": len(model.patterns),
        "build_s": build_s,
        "objective_scale": objective_scale,
        "deduplicate": args.deduplicate,
    }
    return metadata, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--sites", type=int, default=500)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--rep", type=int, default=0)
    parser.add_argument("--mode", choices=["two-point", "spectral"], required=True)
    parser.add_argument("--fingerprint-columns", type=int, default=2)
    parser.add_argument("--fingerprint-offset", type=int, default=0)
    parser.add_argument("--anchor-weighting", choices=["uniform", "nearest"], default="uniform")
    parser.add_argument("--spectral-support-points", type=int, default=5)
    parser.add_argument("--saturation-threshold", type=int, default=2)
    parser.add_argument("--root-planes", type=int, default=4)
    parser.add_argument("--top-k", type=int, default=30)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--profile", choices=["full", "confirm"], default="full")
    parser.add_argument("--deduplicate", action="store_true")
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()
    metadata, rows = run(args)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    with prefix.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(metadata), flush=True)
    for row in rows:
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
