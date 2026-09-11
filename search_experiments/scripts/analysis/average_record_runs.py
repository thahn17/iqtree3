#!/usr/bin/env python3
"""
Averages repeated runs within a record_<model><tag>.csv-style file (columns
run_id,candidates,time_elapsed,logL,true_minus_current -- see
appendRecordRow in tree/spr_topology_test.cpp) into one approximate
"mean run" per distinct parameter combination.

Only operates on files that have already been split by
split_findopt_csv.py, i.e. files named <base>_regular.csv or
<base>_findopt.csv where <base>.csv also exists (that's how a genuine
split output is told apart from an ordinary record file that just happens
to end in "_findopt.csv" because that's its own recordTag, e.g.
record_GTR_FO_findopt.csv is NOT a split output -- there's no
record_GTR_FO.csv to have been split from).

Grouping "same exact parameters": run_id is
"<timestamp>_<param-tail>" (see buildRunId), e.g.
"20260804-164420_r10_s10000_fast_findopt500" -- runs are grouped by the
param-tail alone, so repeated invocations with identical flags (differing
only by when they were run) land in the same group regardless of their
timestamp.

Averaging "into a single run per time, approximate where time does not
align": within a parameter group, each run_id's own rows already form a
real trajectory over its own time_elapsed values, but different runs
essentially never land on the exact same time_elapsed. So the "mean run"
is built on the union of every run's own time_elapsed values in the group
(the finest grid available); at each such time t, every OTHER run in the
group that has real data bracketing t (or exactly at it) contributes its
own piecewise-linearly-interpolated value there, and the row written is
the average of every run's value at t. A run contributes at t only within
its own observed time range [its first row's time, its last row's time]
-- deliberately no extrapolation before a run starts or after it ends, so
a run that finished early or started late simply doesn't weigh in on
times outside what it actually measured (the printed per-t "n" isn't
constant across the file for this reason).

Usage:
    python3 test_scripts/average_record_runs.py [directory] [prefix]

    [directory]  where to look for <prefix>*_regular.csv / *_findopt.csv
                 (default: current directory)
    [prefix]     filename prefix to scan for (default: zrecord)

For every genuine split output <name>.csv found, writes <name>_mean.csv
(same 5-column header), and prints, per parameter group, how many run_ids
were merged and how many rows resulted.
"""

import csv
import glob
import os
import re
import sys
from bisect import bisect_left
from collections import defaultdict

RUN_ID_RE = re.compile(r"^(\d{8}-\d{6})_(.+)$")


def parse_param_tail(run_id):
    m = RUN_ID_RE.match(run_id)
    return m.group(2) if m else run_id


def interpolate(times, values, t):
    """Linear interpolation of `values` (parallel to sorted `times`) at t.
    Returns None if t is outside [times[0], times[-1]] (no extrapolation)."""
    if t < times[0] or t > times[-1]:
        return None
    i = bisect_left(times, t)
    if i < len(times) and times[i] == t:
        return values[i]
    # times[i-1] < t < times[i]
    t0, t1 = times[i - 1], times[i]
    v0, v1 = values[i - 1], values[i]
    frac = (t - t0) / (t1 - t0)
    return v0 + frac * (v1 - v0)


def average_selected(series, run_ids, mean_run_id):
    """Average an explicit set of run_ids together into one approximate
    trajectory (see module docstring for the time-alignment method).

    series: dict run_id -> (times, candidates_vals, logL_vals,
            true_minus_current_vals), each parallel list sorted by time.
    run_ids: which keys of `series` to average (order doesn't matter).
    mean_run_id: run_id string to stamp on every output row.

    Returns a list of (run_id, candidates, time_elapsed, logL,
    true_minus_current) tuples, one per distinct time in the union of the
    selected runs' own time_elapsed values, each an average across
    whichever of those runs have real data covering that time (see
    interpolate -- no extrapolation past a run's own first/last row).
    """
    grid = sorted({t for run_id in run_ids for t in series[run_id][0]})

    out_rows = []
    for t in grid:
        cands, logls, diffs = [], [], []
        for run_id in run_ids:
            times, c_vals, l_vals, d_vals = series[run_id]
            c = interpolate(times, c_vals, t)
            if c is None:
                continue
            cands.append(c)
            logls.append(interpolate(times, l_vals, t))
            diffs.append(interpolate(times, d_vals, t))
        if not cands:
            continue
        out_rows.append((
            mean_run_id,
            round(sum(cands) / len(cands)),
            t,
            sum(logls) / len(logls),
            sum(diffs) / len(diffs),
        ))
    return out_rows


def average_file(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    try:
        run_id_col = header.index("run_id")
        candidates_col = header.index("candidates")
        time_col = header.index("time_elapsed")
        logl_col = header.index("logL")
        diff_col = header.index("true_minus_current")
    except ValueError:
        print(f"  skipping {path}: unexpected header {header}", file=sys.stderr)
        return

    # group rows by run_id, then run_ids by parameter tail
    runs = defaultdict(list)
    for row in rows:
        runs[row[run_id_col]].append((
            float(row[time_col]),
            float(row[candidates_col]),
            float(row[logl_col]),
            float(row[diff_col]),
        ))

    param_groups = defaultdict(list)
    for run_id in runs:
        param_groups[parse_param_tail(run_id)].append(run_id)

    # sort each run's own rows by time
    for run_id in runs:
        runs[run_id].sort(key=lambda r: r[0])

    out_rows = []
    for param_tail, run_ids in param_groups.items():
        series = {}
        for run_id in run_ids:
            pts = runs[run_id]
            series[run_id] = (
                [p[0] for p in pts],
                [p[1] for p in pts],  # candidates
                [p[2] for p in pts],  # logL
                [p[3] for p in pts],  # true_minus_current
            )

        mean_run_id = f"{param_tail}_mean_of_{len(run_ids)}"
        group_rows = average_selected(series, run_ids, mean_run_id)
        out_rows.extend(group_rows)
        print(f"  {param_tail}: merged {len(run_ids)} run(s) -> {len(group_rows)} rows")

    out_path = path[:-4] + "_mean.csv" if path.endswith(".csv") else path + "_mean.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(out_rows)
    print(f"  wrote {len(out_rows)} rows to {out_path}")


def find_split_outputs(directory, prefix):
    candidates = sorted(
        glob.glob(os.path.join(directory, f"{prefix}*_regular.csv")) +
        glob.glob(os.path.join(directory, f"{prefix}*_findopt.csv"))
    )
    targets = []
    for path in candidates:
        for suffix in ("_regular.csv", "_findopt.csv"):
            if path.endswith(suffix):
                base = path[: -len(suffix)] + ".csv"
                if os.path.isfile(base):
                    targets.append(path)
                break
    return targets


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    prefix = sys.argv[2] if len(sys.argv) > 2 else "zrecord"

    targets = find_split_outputs(directory, prefix)
    if not targets:
        print(f"No split outputs found for prefix '{prefix}' in {directory}", file=sys.stderr)
        sys.exit(1)

    for path in targets:
        print(f"{path}:")
        average_file(path)


if __name__ == "__main__":
    main()
