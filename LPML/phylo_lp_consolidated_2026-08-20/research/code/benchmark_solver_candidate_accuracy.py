#!/usr/bin/env python3
"""Check whether a faster exact solver configuration preserves best rescored trees."""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PIPELINE = ROOT / "pipeline_v3"
V3_EXPERIMENTS = HERE / "v3_experiments"
sys.path[:0] = [str(V3_EXPERIMENTS), str(PIPELINE), str(HERE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import lp_exact_speedups as exact  # noqa: E402
import multisite_saturated_100taxa as saturated  # noqa: E402
from phylo_branch_opt import TreeLikelihood  # noqa: E402


def decode(model, result, top_k):
    candidates, fractional = saturated.rank_fingerprint_candidates(
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
    return candidates, fractional


def run_case(args, rep, mode):
    seed = 932000 + 1000 * args.n + rep
    rng = random.Random(seed)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    true_tree = bench.random_tree(args.n, rng)
    sequences, columns, _ = bench.simulate(
        true_tree, args.sites, simulation_model, rng, args.regime
    )
    fitted_model = bench.empirical_model(columns)
    bench.configure(fitted_model)
    model = saturated.FingerprintSaturatedModel(
        args.library,
        columns,
        fingerprint_columns=args.fingerprint_columns,
        anchor_weighting=args.anchor_weighting,
        substitution_model=fitted_model,
        branch_relaxation=mode,
        spectral_support_points=args.spectral_support_points,
        saturation_threshold=2,
        root_planes=4,
    )
    compiled = exact.CompiledSparseLP(model.lp)

    primary_value, primary_result, primary_timing = compiled.solve(model.obj, method="highs")
    fast_value, fast_result, fast_timing = compiled.solve(
        model.obj,
        method="highs-ds",
        simplex_dual_edge_weight_strategy="devex",
    )
    primary, primary_fractional = decode(model, primary_result, args.top_k)
    fast, fast_fractional = decode(model, fast_result, args.top_k)

    labels = tuple(f"leaf_{index + 1}" for index in range(args.n))
    by_key = {}
    source = {}
    for tag, candidates in (("primary", primary), ("devex", fast)):
        for candidate in candidates:
            marker = bench.key(candidate.topology)
            by_key.setdefault(marker, candidate.topology)
            source.setdefault(marker, set()).add(tag)
    scores = {}
    started = time.perf_counter()
    for marker, topology in by_key.items():
        scores[marker] = TreeLikelihood(
            topology,
            labels,
            sequences,
            fitted_model,
            short_time=bench.pg.base.branch_times,
            long_ratio=2.0,
        ).optimize(multistart=args.multistart, maxiter=150, ftol=1e-8).log_likelihood
    rescore_s = time.perf_counter() - started

    def best(candidates):
        if not candidates:
            return -float("inf"), ""
        pairs = [(scores[bench.key(candidate.topology)], bench.key(candidate.topology)) for candidate in candidates]
        return max(pairs)

    primary_best, primary_key = best(primary)
    fast_best, fast_key = best(fast)
    union_best = max(scores.values()) if scores else -float("inf")
    return {
        "n": args.n,
        "sites": args.sites,
        "regime": args.regime,
        "rep": rep,
        "mode": mode,
        "anchor_weighting": args.anchor_weighting,
        "variables": len(model.lp.names),
        "primary_solve_s": primary_timing.solve_s,
        "devex_solve_s": fast_timing.solve_s,
        "speedup": primary_timing.solve_s / fast_timing.solve_s,
        "objective_abs_diff": abs(primary_value - fast_value),
        "primary_candidates": len(primary),
        "devex_candidates": len(fast),
        "candidate_intersection": len(
            {bench.key(item.topology) for item in primary}
            & {bench.key(item.topology) for item in fast}
        ),
        "primary_best_logL": primary_best,
        "devex_best_logL": fast_best,
        "devex_minus_primary": fast_best - primary_best,
        "best_tree_same": int(primary_key == fast_key),
        "union_best_logL": union_best,
        "rescore_s": rescore_s,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--sites", type=int, default=500)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--modes", nargs="+", choices=["two-point", "spectral"], default=["two-point", "spectral"])
    parser.add_argument("--fingerprint-columns", type=int, default=2)
    parser.add_argument("--anchor-weighting", choices=["uniform", "nearest"], default="nearest")
    parser.add_argument("--spectral-support-points", type=int, default=5)
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--multistart", type=int, default=2)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows = []
    for rep in range(args.reps):
        for mode in args.modes:
            row = run_case(args, rep, mode)
            rows.append(row)
            print(json.dumps(row), flush=True)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
