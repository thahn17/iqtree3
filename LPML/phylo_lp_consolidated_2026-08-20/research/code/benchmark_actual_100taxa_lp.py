#!/usr/bin/env python3
"""Run the actual fixed-flow LP on a 100-taxon library.

Unlike ``benchmark_100taxa_hybrid.py``, this script instantiates all 4,012 split
units and the complete shared router/taxon-flow layer.  It is deliberately a
single-case runner so large cases can be scheduled and monitored independently.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
import time
import types
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PIPELINE = ROOT / "pipeline_v3"
V3_EXPERIMENTS = HERE / "v3_experiments"
sys.path[:0] = [str(V3_EXPERIMENTS), str(PIPELINE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import fractional_tree_hierarchical_units as decoder  # noqa: E402
import multisite_pruned_general as pg  # noqa: E402
import multisite_saturated_100taxa as saturated  # noqa: E402


def configure_solver(sparse_lp, method):
    if method == "highs":
        return

    def solve(self, objective, maximize=True):
        vector = np.zeros(len(self.names))
        for index, value in objective.items():
            vector[index] = -value if maximize else value
        started = time.perf_counter()
        result = linprog(
            vector,
            A_ub=self._mat(self.le),
            b_ub=np.asarray(self.lerhs) if self.le else None,
            A_eq=self._mat(self.eq),
            b_eq=np.asarray(self.eqrhs) if self.eq else None,
            bounds=self.bounds,
            method=method,
            options={"dual_feasibility_tolerance": 1e-8, "primal_feasibility_tolerance": 1e-8},
        )
        elapsed = time.perf_counter() - started
        if not result.success:
            raise RuntimeError(result.message)
        value = sum(coefficient * result.x[index] for index, coefficient in objective.items())
        return value, result, elapsed

    sparse_lp.solve = types.MethodType(solve, sparse_lp)


def run(args):
    seed = 127000000 + args.rep + (1000000 if args.regime == "wide" else 0)
    rng = random.Random(seed)
    simulation_model = bench.gtr_f(bench.PI, bench.RATES)
    true_tree = bench.random_tree(args.n, rng)
    sequences, columns, _ = bench.simulate(true_tree, args.sites, simulation_model, rng, args.regime)
    model = bench.empirical_model(columns)
    bench.configure(model)

    print(json.dumps({"event": "build_start", "n": args.n, "sites": args.sites, "mode": args.mode}), flush=True)
    started = time.perf_counter()
    if args.saturated_threshold:
        if args.fingerprint_columns and args.projected_routing:
            model_class = saturated.ProjectedFingerprintSaturatedModel
        elif args.fingerprint_columns:
            model_class = saturated.FingerprintSaturatedModel
        else:
            model_class = saturated.SaturatedPrunedModel
        extra = (
            {
                "fingerprint_columns": args.fingerprint_columns,
                "fingerprint_offset": args.fingerprint_offset,
                "anchor_weighting": args.anchor_weighting,
            }
            if args.fingerprint_columns
            else {}
        )
        lp = model_class(
            args.library,
            columns,
            substitution_model=model,
            branch_relaxation=args.mode,
            spectral_support_points=args.spectral_support_points,
            saturation_threshold=args.saturated_threshold,
            root_planes=args.root_planes,
            log_floor=args.log_floor,
            **extra,
        )
    else:
        lp = pg.PrunedGeneralModel(
            args.library,
            columns,
            preprocess_mode=args.preprocess,
            substitution_model=model,
            branch_relaxation=args.mode,
            spectral_support_points=args.spectral_support_points,
        )
    build_s = time.perf_counter() - started

    configure_solver(lp.lp, args.solver_method)
    if args.build_only:
        summary = {
            "n": args.n,
            "sites": args.sites,
            "regime": args.regime,
            "rep": args.rep,
            "mode": args.mode,
            "saturated_threshold": args.saturated_threshold,
            "root_planes": args.root_planes,
            "log_floor": args.log_floor,
            "fingerprint_columns": args.fingerprint_columns,
            "fingerprint_offset": args.fingerprint_offset,
            "anchor_weighting": args.anchor_weighting,
            "projected_routing": args.projected_routing,
            "solver_method": args.solver_method,
            "split_units": len(lp.splits),
            "variables": len(lp.lp.names),
            "inequalities": len(lp.lp.le),
            "equalities": len(lp.lp.eq),
            "build_s": build_s,
            "build_only": True,
        }
        return summary, []
    print(json.dumps({"event": "solve_start", "build_s": build_s, "variables": len(lp.lp.names)}), flush=True)
    started = time.perf_counter()
    lp_score, result, solve_meta = lp.solve()
    solve_s = time.perf_counter() - started

    print(json.dumps({"event": "decode_start", "solve_s": solve_s}), flush=True)
    started = time.perf_counter()
    ranker = saturated.rank_fingerprint_candidates if args.fingerprint_columns else decoder.rank_integral_tree_candidates
    candidates, _ = ranker(
        lp,
        result,
        topk=args.top_k,
        root_beam=args.root_beam,
        leaf_beam=args.leaf_beam,
        child_branching=args.child_branching,
        leaf_branching=args.leaf_branching,
        merge_branching=args.merge_branching,
        max_expansions=args.max_expansions,
    )
    decode_s = time.perf_counter() - started

    true_key = bench.key(true_tree)
    true_rank = 999999
    records = []
    for rank, item in enumerate(candidates, 1):
        if bench.key(item.topology) == true_key:
            true_rank = rank
        records.append(
            {
                "rank": rank,
                "newick": item.newick,
                "bottleneck": float(item.proportion),
                "secondary": float(item.secondary_score),
                "direction": item.direction,
                "is_generating": int(bench.key(item.topology) == true_key),
            }
        )

    summary = {
        "n": args.n,
        "sites": args.sites,
        "regime": args.regime,
        "rep": args.rep,
        "mode": args.mode,
        "preprocess": args.preprocess,
        "saturated_threshold": args.saturated_threshold,
        "root_planes": args.root_planes,
        "log_floor": args.log_floor,
        "fingerprint_columns": args.fingerprint_columns,
        "fingerprint_offset": args.fingerprint_offset,
        "anchor_weighting": args.anchor_weighting,
        "projected_routing": args.projected_routing,
        "solver_method": args.solver_method,
        "split_units": len(lp.splits),
        "variables": len(lp.lp.names),
        "inequalities": len(lp.lp.le),
        "equalities": len(lp.lp.eq),
        "build_s": build_s,
        "solve_s": solve_s,
        "decode_s": decode_s,
        "lp_score": float(lp_score),
        "candidate_count": len(candidates),
        "generating_rank": true_rank,
        "solve_status": str(solve_meta),
        "generating_newick": decoder._nw(true_tree, tuple(f"leaf_{i+1}" for i in range(args.n))) + ";",
    }
    return summary, records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", required=True)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--sites", type=int, required=True)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--rep", type=int, default=0)
    parser.add_argument("--mode", choices=["two-point", "spectral"], required=True)
    parser.add_argument("--preprocess", choices=["fast", "balanced", "accuracy"], default="fast")
    parser.add_argument("--spectral-support-points", type=int, default=5)
    parser.add_argument("--saturated-threshold", type=int, choices=[1, 2, 3], default=2)
    parser.add_argument("--root-planes", type=int, default=4)
    parser.add_argument(
        "--log-floor",
        type=float,
        default=1e-9,
        help="Lowest log-chord likelihood support; 1e-9 matches HiGHS small_matrix_value",
    )
    parser.add_argument("--fingerprint-columns", type=int, default=0)
    parser.add_argument("--fingerprint-offset", type=int, default=0)
    parser.add_argument("--anchor-weighting", choices=["uniform", "nearest"], default="uniform")
    parser.add_argument("--projected-routing", action="store_true")
    parser.add_argument("--solver-method", choices=["highs", "highs-ipm", "highs-ds"], default="highs")
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--root-beam", type=int, default=128)
    parser.add_argument("--leaf-beam", type=int, default=32)
    parser.add_argument("--child-branching", type=int, default=12)
    parser.add_argument("--leaf-branching", type=int, default=8)
    parser.add_argument("--merge-branching", type=int, default=96)
    parser.add_argument("--max-expansions", type=int, default=500000)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()
    summary, records = run(args)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    with prefix.with_suffix(".candidates.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]) if records else ["rank", "newick", "bottleneck", "secondary", "direction", "is_generating"])
        writer.writeheader()
        writer.writerows(records)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
