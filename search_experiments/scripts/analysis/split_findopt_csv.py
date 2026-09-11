#!/usr/bin/env python3
"""
Splits a record_<model><tag>.csv file (columns: run_id,candidates,
time_elapsed,logL,true_minus_current -- see appendRecordRow in
tree/spr_topology_test.cpp) written by a --hillclimb run using both
"record" and "findopt" into two files: one with the real search-progress
rows, one with findopt's own diagnostic checkpoint rows.

appendRecordRow writes the exact same 5 columns for both kinds of row --
there's no explicit "source" marker -- so rows are told apart by a
property of how they're produced, not by reading a column directly:

  Primary method (requires "fast N" and no "investigate", both readable
  off run_id, e.g. "..._fast5_findopt50..."): in "fast N" mode every step
  scores EXACTLY N candidates, so cumulative `candidates` after step k is
  exactly (k+1)*N -- deterministic, no exceptions. findopt fires every F
  steps (the number right after "findopt" in run_id), i.e. right after
  step k where (k+1) % F == 0, so a findopt row's candidates count is an
  exact multiple of F*N, while an accepted-move row's candidates count is
  still a multiple of N but lands at an essentially random step number, so
  is only a multiple of F*N by chance. "investigate" breaks this: it adds
  extra, irregular scoreTrialSPRMove calls right after accepted moves, so
  candidates is no longer a clean multiple of N at all.

  Fallback method (used when "fast" isn't in run_id, or "investigate" is,
  so the multiple-of-F*N test doesn't apply): regular rows are only ever
  appended when curScore improves, so within one run_id their logL values
  are strictly non-decreasing; findopt rows are a one-off scratch refit
  that never touches curScore, so they don't participate in that trend.
  This finds the longest non-decreasing subsequence (LNDS) of logL per
  run_id and calls that "regular" -- an approximation, not a certainty.

Either way, the guaranteed final row of each run_id (runHillClimb always
appends one unconditionally at the end, whether or not the last real event
was a findopt checkpoint) is always treated as regular -- see
runHillClimb's own comment on that final appendRecordRow call.

Rows are grouped by run_id first (never compared across different runs),
since neither method's assumptions hold across a run boundary.

Usage:
    python3 test_scripts/split_findopt_csv.py <input.csv> [output_prefix]

    <input.csv>      a record_*.csv containing both regular and findopt
                     rows (i.e. generated with "findopt ... record")
    [output_prefix]  base name for the two output files (default: input
                     path with its extension stripped)

Writes <prefix>_regular.csv and <prefix>_findopt.csv, each with the
original header, and prints per-run_id which method was used and a row
count, so you can sanity-check the split (e.g. against how many findopt
checkpoints you expected from the findoptEveryNSteps/max-steps you
actually ran with).

Example:
    python3 test_scripts/split_findopt_csv.py record_JC_fast_findopt.csv
        -> record_JC_fast_findopt_regular.csv
           record_JC_fast_findopt_findopt.csv
"""

import csv
import re
import sys
from bisect import bisect_right
from collections import defaultdict

FINDOPT_RE = re.compile(r"_findopt(\d+)")
FAST_RE = re.compile(r"_fast(\d*)(?:_|$)")
INVESTIGATE_RE = re.compile(r"_investigate")


def longest_nondecreasing_indices(values):
    """Indices (into `values`, ascending) of one longest non-decreasing
    subsequence -- standard O(n log n) patience-style algorithm
    (bisect_right so equal values extend a run)."""
    tails_val = []
    tails_idx = []
    parent = [-1] * len(values)

    for i, v in enumerate(values):
        pos = bisect_right(tails_val, v)
        if pos == len(tails_val):
            tails_val.append(v)
            tails_idx.append(i)
        else:
            tails_val[pos] = v
            tails_idx[pos] = i
        parent[i] = tails_idx[pos - 1] if pos > 0 else -1

    if not tails_idx:
        return []

    result = []
    cur = tails_idx[-1]
    while cur != -1:
        result.append(cur)
        cur = parent[cur]
    result.reverse()
    return result


def classify_group(run_id, candidates, logls):
    """Returns a list of bool, True = regular, for one run_id's rows (in
    original file order). Picks the multiple-of-F*N method when run_id
    makes it applicable, else falls back to the LNDS method."""
    findopt_m = FINDOPT_RE.search(run_id)
    fast_m = FAST_RE.search(run_id)
    has_investigate = INVESTIGATE_RE.search(run_id) is not None

    if findopt_m and fast_m and not has_investigate:
        f_interval = int(findopt_m.group(1))
        per_step = int(fast_m.group(1)) if fast_m.group(1) else 1
        stride = f_interval * per_step
        is_regular = [not (c > 0 and c % stride == 0) for c in candidates]
        is_regular[-1] = True  # guaranteed final row is always regular
        return is_regular, f"multiple-of-{stride} (fast{per_step} x findopt{f_interval})"

    method = "LNDS fallback"
    if has_investigate:
        method += " (investigate present -- multiple-of-F*N test doesn't apply)"
    elif not fast_m:
        method += " (no 'fast' in run_id -- no fixed per-step candidate count)"

    regular_idx = set(longest_nondecreasing_indices(logls))
    is_regular = [i in regular_idx for i in range(len(logls))]
    is_regular[-1] = True
    return is_regular, method


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.csv> [output_prefix]", file=sys.stderr)
        sys.exit(1)

    in_path = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else in_path.rsplit(".", 1)[0]

    with open(in_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    if not rows:
        print(f"{in_path}: no data rows found", file=sys.stderr)
        sys.exit(1)

    try:
        run_id_col = header.index("run_id")
        candidates_col = header.index("candidates")
        logl_col = header.index("logL")
    except ValueError:
        print(f"{in_path}: expected 'run_id', 'candidates', 'logL' columns, got {header}", file=sys.stderr)
        sys.exit(1)

    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[row[run_id_col]].append(i)

    is_regular = [False] * len(rows)
    for run_id, idxs in groups.items():
        candidates = [int(rows[i][candidates_col]) for i in idxs]
        logls = [float(rows[i][logl_col]) for i in idxs]
        group_regular, method = classify_group(run_id, candidates, logls)
        for local_i, keep in zip(idxs, group_regular):
            is_regular[local_i] = keep
        n_regular = sum(group_regular)
        print(f"run_id={run_id}: {n_regular} regular, {len(idxs) - n_regular} findopt "
              f"(of {len(idxs)} rows) -- method: {method}")

    regular_path = f"{prefix}_regular.csv"
    findopt_path = f"{prefix}_findopt.csv"

    with open(regular_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(row for row, keep in zip(rows, is_regular) if keep)

    with open(findopt_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(row for row, keep in zip(rows, is_regular) if not keep)

    print(f"Wrote {sum(is_regular)} regular rows to {regular_path}")
    print(f"Wrote {len(rows) - sum(is_regular)} findopt rows to {findopt_path}")


if __name__ == "__main__":
    main()
