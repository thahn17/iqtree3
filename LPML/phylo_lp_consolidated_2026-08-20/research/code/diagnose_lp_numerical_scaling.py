#!/usr/bin/env python3
"""Report coefficient, bound, and solution scales for a real phylogeny LP."""
from __future__ import annotations

import argparse
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path[:0] = [str(HERE / "v3_experiments"), str(ROOT / "pipeline_v3"), str(HERE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import lp_exact_speedups as exact  # noqa: E402
import multisite_saturated_100taxa as saturated  # noqa: E402


def quantiles(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {}
    return {
        str(q): float(np.quantile(values, q))
        for q in (0, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 1)
    }


def run(args):
    rng = random.Random(972000 + args.n + args.rep)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    tree = bench.random_tree(args.n, rng)
    _, columns, _ = bench.simulate(tree, args.sites, simulation_model, rng, args.regime)
    fitted = bench.empirical_model(columns)
    bench.configure(fitted)
    model_class = (
        saturated.ProjectedFingerprintSaturatedModel
        if args.projected_routing
        else saturated.FingerprintSaturatedModel
    )
    model = model_class(
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
    coefficients = [
        abs(float(value))
        for rows in (model.lp.le, model.lp.eq)
        for row in rows
        for value in row.values()
        if value != 0
    ]
    finite_upper = [float(upper) for _, upper in model.lp.bounds if upper is not None and upper > 0]
    objective = [abs(float(value)) for value in model.obj.values() if value != 0]
    value, result, timing = exact.CompiledSparseLP(model.lp).solve(model.obj)
    positive = result.x[result.x > 0]
    groups = defaultdict(list)
    for name, number in zip(model.lp.names, result.x):
        groups[name.split(":", 1)[0]].append(float(number))
    tiny_by_group = {
        group: {
            "count": len(values),
            "positive": sum(value > 0 for value in values),
            "between_1e-12_1e-8": sum(1e-12 < value < 1e-8 for value in values),
            "at_most_1e-12": sum(0 < value <= 1e-12 for value in values),
        }
        for group, values in groups.items()
    }
    report = {
        "n": args.n,
        "sites": args.sites,
        "mode": args.mode,
        "projected_routing": args.projected_routing,
        "variables": len(model.lp.names),
        "inequalities": len(model.lp.le),
        "equalities": len(model.lp.eq),
        "coefficient_abs_quantiles": quantiles(coefficients),
        "coefficients_below_1e-9": sum(value < 1e-9 for value in coefficients),
        "coefficients_below_1e-12": sum(value < 1e-12 for value in coefficients),
        "finite_positive_upper_bound_quantiles": quantiles(finite_upper),
        "upper_bounds_below_1e-8": sum(value < 1e-8 for value in finite_upper),
        "objective_abs_quantiles": quantiles(objective),
        "solution_positive_quantiles": quantiles(positive),
        "solution_positive_below_1e-9": int(np.count_nonzero((result.x > 0) & (result.x < 1e-9))),
        "solution_positive_below_1e-12": int(np.count_nonzero((result.x > 0) & (result.x < 1e-12))),
        "objective": float(value),
        "solve_s": timing.solve_s,
        "variable_groups": tiny_by_group,
    }
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--sites", type=int, default=1000)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--rep", type=int, default=0)
    parser.add_argument("--mode", choices=["two-point", "spectral"], required=True)
    parser.add_argument("--fingerprint-columns", type=int, default=4)
    parser.add_argument("--projected-routing", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    report = run(args)
    Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
