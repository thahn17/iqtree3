#!/usr/bin/env python3
"""
One-off script: processes the yrecord_*.csv files (batches dated the 7th
or 8th only), splits each complete batch of 10 same-parameter runs into
findopt/regular points, averages each, and writes the results into
PhylogenyTreeLengths.xlsx as new sheets plus new chart series on the
existing "r=10 vs r=1" chart, the existing "Standard/BLO/fullreopt" chart,
and a new "fast r10 regular vs shrink vs investigate" chart.

Incomplete trailing batches (<10 runs) are still split+averaged (to prove
the pipeline handles them) but never written into the workbook.

Not meant to be reusable/general -- run once, inspect the result, discard.
"""

import csv
import re
import sys
from collections import defaultdict

sys.path.insert(0, "test_scripts")
from average_record_runs import average_selected
from split_findopt_csv import classify_group

import openpyxl
from openpyxl.chart import ScatterChart, Series, Reference
from openpyxl.chart.marker import Marker
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.axis import NumericAxis, Scaling
from openpyxl.chart.data_source import NumData, NumVal
from openpyxl.drawing.line import LineProperties

RUN_ID_RE = re.compile(r"^(\d{8})-(\d{6})_(.+)$")
HEADER = ["run_id", "candidates", "time_elapsed", "logL", "true_minus_current"]
QUALIFYING_DATES = ("20260807", "20260808")


def load_series(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)
    run_id_col = header.index("run_id")
    cand_col = header.index("candidates")
    time_col = header.index("time_elapsed")
    logl_col = header.index("logL")
    diff_col = header.index("true_minus_current")
    runs = defaultdict(list)
    for row in rows:
        runs[row[run_id_col]].append((
            float(row[time_col]), float(row[cand_col]), float(row[logl_col]), float(row[diff_col])
        ))
    series = {}
    for rid, pts in runs.items():
        pts.sort(key=lambda p: p[0])
        series[rid] = ([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts], [p[3] for p in pts])
    return series


def find_batches(path):
    """Yields (tail, label, run_ids, complete, series) for every
    contiguous same-parameter chunk of up to 10 runs, filtered to
    QUALIFYING_DATES. `label` is a/b/c/... per distinct tail in this file."""
    series = load_series(path)
    run_ids = sorted(series.keys())
    parsed = [RUN_ID_RE.match(r) for r in run_ids]
    tails = [m.group(3) for m in parsed]
    dates = [m.group(1) for m in parsed]

    counters = defaultdict(int)
    i = 0
    while i < len(run_ids):
        tail = tails[i]
        j = i
        while j < len(run_ids) and tails[j] == tail and (j - i) < 10:
            j += 1
        chunk = run_ids[i:j]
        date = dates[i]
        complete = len(chunk) == 10
        if date in QUALIFYING_DATES:
            label = "abcdefgh"[counters[tail]]
            counters[tail] += 1
            yield tail, label, chunk, complete, series
        i = j


def split_and_average(chunk, series, mean_run_id_base):
    """Returns (regular_rows, findopt_rows). findopt_rows is [] if no run
    in the chunk contributed any findopt-classified rows at all."""
    regular_series = {}
    findopt_series = {}
    findopt_run_ids = []
    for rid in chunk:
        times, cand, logl, diff = series[rid]
        is_regular, _method = classify_group(rid, cand, logl)
        reg_idx = [k for k, r in enumerate(is_regular) if r]
        fo_idx = [k for k, r in enumerate(is_regular) if not r]
        regular_series[rid] = (
            [times[k] for k in reg_idx], [cand[k] for k in reg_idx],
            [logl[k] for k in reg_idx], [diff[k] for k in reg_idx],
        )
        if fo_idx:
            findopt_series[rid] = (
                [times[k] for k in fo_idx], [cand[k] for k in fo_idx],
                [logl[k] for k in fo_idx], [diff[k] for k in fo_idx],
            )
            findopt_run_ids.append(rid)

    regular_rows = average_selected(regular_series, chunk, mean_run_id_base + "_regular")
    findopt_rows = average_selected(findopt_series, findopt_run_ids, mean_run_id_base + "_findopt") if findopt_run_ids else []
    return regular_rows, findopt_rows


def write_sheet(wb, name, rows):
    name = name[:31]
    if name in wb.sheetnames:
        del wb[name]
    ws = wb.create_sheet(name)
    ws.append(HEADER)
    for row in rows:
        ws.append(list(row))
    return ws


def build_numcache(values):
    return NumData(formatCode="General", ptCount=len(values),
                    pt=[NumVal(idx=i, v=v) for i, v in enumerate(values)])


def add_series(chart, ws, rows, title):
    """rows: the same [run_id, candidates, time_elapsed, logL,
    true_minus_current] rows written to `ws`. numCache is built explicitly
    from these values rather than left for openpyxl to populate on save --
    that turned out to be inconsistent across runs (sometimes producing a
    numRef with no cache at all), unlike every pre-existing series in this
    workbook, which all carry a full cache."""
    n_rows = len(rows)
    x_values = [r[2] for r in rows]
    y_values = [r[3] for r in rows]

    xref = Reference(ws, min_col=3, min_row=2, max_row=n_rows + 1)
    yref = Reference(ws, min_col=4, min_row=2, max_row=n_rows + 1)
    s = Series(yref, xref, title=title)
    s.xVal.numRef.numCache = build_numcache(x_values)
    s.yVal.numRef.numCache = build_numcache(y_values)
    # line only, no marker at each point -- matches the minimal
    # <c:marker><c:symbol val="none"/></c:marker> the existing series use
    s.marker = Marker(symbol="none")
    s.graphicalProperties = GraphicalProperties()
    s.graphicalProperties.line = LineProperties(w=19050)
    s.smooth = False
    chart.series.append(s)
    return s


def main():
    processed = []  # (tag, label, complete, regular_rows, findopt_rows)

    plan = [
        ("yrecord_GTR_FO_fast_findopt.csv", "r10_findopt", lambda tail: tail.startswith("r10_")),
        ("yrecord_GTR_FO_fast_findopt.csv", "r1_findopt", lambda tail: tail.startswith("r1_")),
        ("yrecord_GTR_FO_fast_fullreopt_findopt.csv", "r10_fullreopt", lambda tail: True),
        ("yrecord_GTR_FO_fast_investigate3.csv", "r10_investigate3", lambda tail: True),
        ("yrecord_GTR_FO_fast_shrink.csv", "r10_shrink", lambda tail: True),
    ]

    results = defaultdict(list)  # key -> list of dict(label, complete, regular_rows, findopt_rows)
    for path, key, tail_filter in plan:
        for tail, label, chunk, complete, series in find_batches(path):
            if not tail_filter(tail):
                continue
            mean_base = f"{key}_07{label}"
            regular_rows, findopt_rows = split_and_average(chunk, series, mean_base)
            results[key].append(dict(label=label, complete=complete, chunk=chunk,
                                      regular_rows=regular_rows, findopt_rows=findopt_rows))
            print(f"{key} batch 07{label}: {len(chunk)} runs complete={complete} "
                  f"-> {len(regular_rows)} regular rows, {len(findopt_rows)} findopt rows")

    wb = openpyxl.load_workbook("PhylogenyTreeLengths.xlsx")
    ws1 = wb["Sheet1"]
    chart_radius = ws1._charts[0]   # r=10 / r=1
    chart_strategy = ws1._charts[1]  # Standard / BLO / fullreopt

    written_sheets = {}  # key+label+kind -> (ws, rows)

    def persist(key, entry, kind):
        rows = entry["regular_rows"] if kind == "regular" else entry["findopt_rows"]
        if not entry["complete"] or not rows:
            return None
        sheet_name = f"{key}_07{entry['label']}_{kind}"
        ws = write_sheet(wb, sheet_name, rows)
        written_sheets[(key, entry["label"], kind)] = (ws, rows)
        return ws, rows

    for key in results:
        for entry in results[key]:
            persist(key, entry, "regular")
            persist(key, entry, "findopt")

    # --- Chart A (radius comparison): all r10_findopt + r1_findopt series ---
    for key in ("r10_findopt", "r1_findopt"):
        radius_label = "r10" if key == "r10_findopt" else "r1"
        for entry in results[key]:
            if not entry["complete"]:
                continue
            for kind in ("regular", "findopt"):
                got = written_sheets.get((key, entry["label"], kind))
                if got:
                    ws, rows = got
                    title = f"{radius_label} findopt 07{entry['label']} {kind}"
                    add_series(chart_radius, ws, rows, title)

    # --- Chart B (strategy comparison): r10_findopt (both kinds) + r10_fullreopt (both kinds) ---
    for key in ("r10_findopt", "r10_fullreopt"):
        for entry in results[key]:
            if not entry["complete"]:
                continue
            for kind in ("regular", "findopt"):
                got = written_sheets.get((key, entry["label"], kind))
                if got:
                    ws, rows = got
                    label_kind = "regular" if kind == "regular" else "findopt (BLO)"
                    prefix = "r10 findopt" if key == "r10_findopt" else "r10 fullreopt"
                    title = f"{prefix} 07{entry['label']} {label_kind}"
                    add_series(chart_strategy, ws, rows, title)

    # --- Chart C (new): r10_findopt regular + investigate3 + shrink ---
    chart_c = ScatterChart()
    chart_c.scatterStyle = "lineMarker"
    chart_c.title = "fast r10 regular vs shrink vs investigate"
    # explicit, fresh axis objects -- a freshly constructed ScatterChart()
    # can otherwise end up holding onto stray scaling state (max/min) left
    # over from whatever was most recently loaded/parsed in this process,
    # so force a clean auto-scaled axis rather than trust the defaults
    chart_c.x_axis = NumericAxis(axId=10, crossAx=20)
    chart_c.y_axis = NumericAxis(axId=20, crossAx=10)
    chart_c.x_axis.scaling = Scaling(orientation="minMax")
    chart_c.y_axis.scaling = Scaling(orientation="minMax")
    chart_c.x_axis.title = "Time (s)"
    chart_c.y_axis.title = "Log Likelihood"
    chart_c.x_axis.delete = False
    chart_c.y_axis.delete = False
    chart_c.x_axis.axPos = "b"
    chart_c.y_axis.axPos = "l"

    for entry in results["r10_findopt"]:
        if not entry["complete"]:
            continue
        got = written_sheets.get(("r10_findopt", entry["label"], "regular"))
        if got:
            ws, rows = got
            add_series(chart_c, ws, rows, f"r10 fast 07{entry['label']} regular")

    for key, human in (("r10_investigate3", "investigate3"), ("r10_shrink", "shrink500")):
        for entry in results[key]:
            if not entry["complete"]:
                continue
            got = written_sheets.get((key, entry["label"], "regular"))
            if got:
                ws, rows = got
                add_series(chart_c, ws, rows, f"r10 {human} 07{entry['label']}")

    ws1.add_chart(chart_c, "S2")

    wb.save("PhylogenyTreeLengths.xlsx")
    print("\nSaved. Chart A series:", len(chart_radius.series))
    print("Chart B series:", len(chart_strategy.series))
    print("Chart C series:", len(chart_c.series))
    print("New sheets written:", len(written_sheets))


if __name__ == "__main__":
    main()
