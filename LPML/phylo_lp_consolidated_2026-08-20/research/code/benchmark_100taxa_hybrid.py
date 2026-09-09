#!/usr/bin/env python3
"""100-taxon out-of-sample calibration for bottleneck retrieval views.

This benchmark intentionally isolates the scalable candidate-retrieval layer.  It does
not claim to be a full 100-taxon spectral/two-point LP solve.  For each simulated GTR+F
alignment it generates one shared candidate pool, computes the whole-alignment local
split score used by the saved 100-taxon experiments, and compares:

* raw bottleneck (q=0);
* a coupled effective-child-size view w=(4ab/m)^(-q);
* a baseline-preserving union of the two rankings.

The generating topology is included in the pool so its structural rank and top-K
recall can be measured.  q must be selected across replicates, never from the held-out
replicate being evaluated.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PIPELINE = ROOT / "pipeline_v3"
V3_EXPERIMENTS = HERE / "v3_experiments"
# The research snapshot contains older modules with hard-coded /mnt/data imports.
# Put the self-contained v3 modules ahead of those snapshots.
sys.path[:0] = [str(V3_EXPERIMENTS), str(PIPELINE), str(HERE)]

import benchmark_spectral_bottleneck as bench  # noqa: E402
import benchmark_100taxa_weighting as bw  # noqa: E402


def _rankings(candidates, data, ref, true_key, q_grid):
    out = {}
    for q in q_grid:
        scored = []
        for tree in candidates:
            value = bw.score(tree, data, ref, -q, -q)
            scored.append((value, bench.key(tree), bw.rf(tree, true_key[1])))
        scored.sort(key=lambda row: row[0], reverse=True)
        keys = [row[1] for row in scored]
        rank = keys.index(true_key[0]) + 1
        out[q] = {
            "keys": keys,
            "true_rank": rank,
            "top_rf": scored[0][2],
        }
    return out


def run_case(n, sites, rep, ncand, regime, q_grid, budgets):
    # The seed is deliberately independent of ``sites``.  Runs at 1k, 5k, and
    # 10k sites therefore use prefixes of the same alignment on the same tree.
    seed = 93000000 + n * 100000 + rep
    if regime == "wide":
        seed += 100000000
    rng = random.Random(seed)
    candidate_rng = random.Random(seed + 700000003)
    sim_model = bench.gtr_f(bench.PI, bench.RATES)
    true = bench.random_tree(n, rng)

    start = time.perf_counter()
    sequences, _, _ = bench.simulate(true, sites, sim_model, rng, regime)
    simulate_s = time.perf_counter() - start

    start = time.perf_counter()
    # A separate seed keeps the candidate pool identical across alignment prefixes.
    candidates = bw.candidates(true, ncand, candidate_rng)
    candidate_s = time.perf_counter() - start

    start = time.perf_counter()
    distances = bw.pair_dist(sequences)
    pair_s = time.perf_counter() - start

    start = time.perf_counter()
    data, ref = bw.build_unit_scores(candidates, distances)
    unit_s = time.perf_counter() - start

    true_key = (bench.key(true), true)
    start = time.perf_counter()
    rankings = _rankings(candidates, data, ref, true_key, q_grid)
    rank_s = time.perf_counter() - start

    rows = []
    raw = rankings[0.0]
    for q in q_grid:
        view = rankings[q]
        for budget in budgets:
            raw_set = set(raw["keys"][:budget])
            q_set = set(view["keys"][:budget])
            union = raw_set | q_set
            # Compact union retains all raw candidates and adds only budget/2 new q-view
            # candidates.  It is useful for runtime estimates but lacks the strict
            # no-regression property of the full union against the q view.
            additions = [key for key in view["keys"] if key not in raw_set]
            compact = raw_set | set(additions[: max(1, budget // 2)])
            rows.append(
                {
                    "n": n,
                    "sites": sites,
                    "regime": regime,
                    "rep": rep,
                    "ncand": len(candidates),
                    "q": q,
                    "budget": budget,
                    "raw_true_rank": raw["true_rank"],
                    "q_true_rank": view["true_rank"],
                    "raw_top_rf": raw["top_rf"],
                    "q_top_rf": view["top_rf"],
                    "raw_recall": int(true_key[0] in raw_set),
                    "q_recall": int(true_key[0] in q_set),
                    "union_recall": int(true_key[0] in union),
                    "union_size": len(union),
                    "compact_recall": int(true_key[0] in compact),
                    "compact_size": len(compact),
                    "simulate_s": simulate_s,
                    "candidate_s": candidate_s,
                    "pair_s": pair_s,
                    "unit_s": unit_s,
                    "rank_s": rank_s,
                    "unique_units": len(data),
                    "score_ref": ref,
                }
            )
    return rows


def write_rows(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--sites", type=int, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--ncand", type=int, default=240)
    parser.add_argument("--regime", choices=["moderate", "wide"], default="moderate")
    parser.add_argument("--q-grid", default="0,0.05,0.10,0.15,0.20,0.25,0.30")
    parser.add_argument("--budgets", default="20,40,80")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    q_grid = tuple(float(x) for x in args.q_grid.split(","))
    if 0.0 not in q_grid:
        raise ValueError("q-grid must include 0")
    budgets = tuple(int(x) for x in args.budgets.split(","))
    rows = run_case(args.n, args.sites, args.rep, args.ncand, args.regime, q_grid, budgets)
    write_rows(args.out, rows)
    timing = {key: rows[0][key] for key in ("simulate_s", "candidate_s", "pair_s", "unit_s", "rank_s")}
    print(json.dumps({"case": [args.n, args.sites, args.regime, args.rep], "timing": timing, "rows": len(rows)}))


if __name__ == "__main__":
    main()
