#!/usr/bin/env python3
"""
Averages an explicit, user-chosen list of run_ids from a record_*.csv-style
file (columns run_id,candidates,time_elapsed,logL,true_minus_current -- see
appendRecordRow in tree/spr_topology_test.cpp) into one approximate "mean
run" -- score (logL) and score-difference (true_minus_current) both
averaged, same time-alignment method average_record_runs.py uses for its
automatic per-parameter-group averaging, but for whatever specific
run_ids you name here, regardless of whether their own parameters
(radius/flags/etc, encoded in the run_id's own text) match each other at
all. Useful for merging runs across different files, or across
intentionally different flag combinations that average_record_runs.py's
automatic exact-parameter-match grouping wouldn't ever put together.

Averaging method (see average_selected in average_record_runs.py): the
output grid is the union of the selected runs' own time_elapsed values;
at each grid time, every selected run with real data spanning that time
contributes a linearly-interpolated value (no extrapolation before a
run's first row or after its last), and the row written is the average
across whichever runs cover that point.

Usage:
    python3 test_scripts/average_selected_runs.py <input.csv> <run_id> [<run_id> ...] [-o output.csv]

    <input.csv>  a record_*.csv (or anything with the same 5-column
                 format, e.g. output from split_findopt_csv.py)
    <run_id>     one or more run_id values to average together, exactly as
                 they appear in the file's run_id column
    -o/--output  output path (default: <input>_selected_mean.csv)

If any given run_id isn't found, the script lists every run_id actually
present in the file and exits without writing anything.

Example:
    python3 test_scripts/average_selected_runs.py record_JC_fast_findopt.csv \\
        20260804-164420_r10_s10000_fast_findopt500 \\
        20260804-165315_r10_s10000_fast_findopt500
"""

import argparse
import csv
import sys
from collections import defaultdict

from average_record_runs import average_selected


def load_series(path):
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
        print(f"{path}: expected the usual 5-column record header, got {header}", file=sys.stderr)
        sys.exit(1)

    runs = defaultdict(list)
    for row in rows:
        runs[row[run_id_col]].append((
            float(row[time_col]),
            float(row[candidates_col]),
            float(row[logl_col]),
            float(row[diff_col]),
        ))

    series = {}
    for run_id, pts in runs.items():
        pts.sort(key=lambda p: p[0])
        series[run_id] = (
            [p[0] for p in pts],
            [p[1] for p in pts],
            [p[2] for p in pts],
            [p[3] for p in pts],
        )
    return header, series


def main():
    parser = argparse.ArgumentParser(description="Average an explicit list of run_ids from a record CSV.")
    parser.add_argument("input_csv")
    parser.add_argument("run_ids", nargs="+", help="run_id values to average together")
    parser.add_argument("-o", "--output", help="output path (default: <input>_selected_mean.csv)")
    args = parser.parse_args()

    header, series = load_series(args.input_csv)

    missing = [r for r in args.run_ids if r not in series]
    if missing:
        print(f"run_id(s) not found in {args.input_csv}:", file=sys.stderr)
        for r in missing:
            print(f"  {r}", file=sys.stderr)
        print(f"Available run_ids ({len(series)}):", file=sys.stderr)
        for r in sorted(series):
            print(f"  {r}", file=sys.stderr)
        sys.exit(1)

    mean_run_id = f"selected_mean_of_{len(args.run_ids)}"
    out_rows = average_selected(series, args.run_ids, mean_run_id)

    if args.output:
        output = args.output
    elif args.input_csv.endswith(".csv"):
        output = args.input_csv[:-4] + "_selected_mean.csv"
    else:
        output = args.input_csv + "_selected_mean.csv"

    with open(output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(out_rows)

    print(f"Averaged {len(args.run_ids)} run(s) -> {len(out_rows)} rows -> {output}")


if __name__ == "__main__":
    main()
