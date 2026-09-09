#!/usr/bin/env python3
"""Compare physical Beneš routing with its count-wise projected master."""
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
    return saturated.rank_fingerprint_candidates(
        model,
        result,
        topk=top_k,
        root_beam=192,
        leaf_beam=48,
        child_branching=16,
        leaf_branching=10,
        merge_branching=128,
        max_expansions=800000,
    )[0]


def one(args, rep, mode):
    rng = random.Random(952000 + 1000 * args.n + rep)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    true_tree = bench.random_tree(args.n, rng)
    sequences, columns, _ = bench.simulate(true_tree, args.sites, simulation_model, rng, args.regime)
    fitted_model = bench.empirical_model(columns)
    bench.configure(fitted_model)
    common = dict(
        library=args.library,
        columns=columns,
        fingerprint_columns=args.fingerprint_columns,
        anchor_weighting="uniform",
        substitution_model=fitted_model,
        branch_relaxation=mode,
        spectral_support_points=5,
        saturation_threshold=2,
        root_planes=4,
    )
    build_started = time.perf_counter()
    physical = saturated.FingerprintSaturatedModel(**common)
    physical_build_s = time.perf_counter() - build_started
    build_started = time.perf_counter()
    projected = saturated.ProjectedFingerprintSaturatedModel(**common)
    projected_build_s = time.perf_counter() - build_started
    p_value, p_result, p_time = exact.CompiledSparseLP(physical.lp).solve(physical.obj)
    q_value, q_result, q_time = exact.CompiledSparseLP(projected.lp).solve(projected.obj)
    p_candidates = decode(physical, p_result, args.top_k)
    q_candidates = decode(projected, q_result, args.top_k)

    p_keys = {bench.key(candidate.topology) for candidate in p_candidates}
    q_keys = {bench.key(candidate.topology) for candidate in q_candidates}
    topologies = {bench.key(candidate.topology): candidate.topology for candidate in p_candidates + q_candidates}
    labels = tuple(f"leaf_{index + 1}" for index in range(args.n))
    scores = {}
    for marker, topology in topologies.items():
        scores[marker] = TreeLikelihood(
            topology,
            labels,
            sequences,
            fitted_model,
            short_time=bench.pg.base.branch_times,
            long_ratio=2.0,
        ).optimize(multistart=args.multistart, maxiter=150, ftol=1e-8).log_likelihood
    p_best = max((scores[key] for key in p_keys), default=-float("inf"))
    q_best = max((scores[key] for key in q_keys), default=-float("inf"))
    return {
        "n": args.n,
        "sites": args.sites,
        "rep": rep,
        "mode": mode,
        "fingerprint_columns": args.fingerprint_columns,
        "physical_variables": len(physical.lp.names),
        "projected_variables": len(projected.lp.names),
        "physical_build_s": physical_build_s,
        "projected_build_s": projected_build_s,
        "physical_solve_s": p_time.solve_s,
        "projected_solve_s": q_time.solve_s,
        "speedup": p_time.solve_s / q_time.solve_s,
        "objective_abs_diff": abs(p_value - q_value),
        "physical_candidates": len(p_keys),
        "projected_candidates": len(q_keys),
        "physical_contained": int(p_keys <= q_keys),
        "candidate_intersection": len(p_keys & q_keys),
        "physical_best_logL": p_best,
        "projected_best_logL": q_best,
        "projected_minus_physical": q_best - p_best,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--sites", type=int, default=500)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--modes", nargs="+", choices=["two-point", "spectral"], default=["two-point", "spectral"])
    parser.add_argument("--fingerprint-columns", type=int, default=4)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--multistart", type=int, default=2)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows = []
    for rep in range(args.reps):
        for mode in args.modes:
            row = one(args, rep, mode)
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
