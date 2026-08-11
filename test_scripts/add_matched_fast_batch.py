#!/usr/bin/env python3
"""
Adds the r10_findopt batch dated 2026-08-06 (10 runs, tail
'r10_s10000_fast_gtr_findopt500', starting logL -766483.5) as new sheets
'r10_findopt_06a_regular'/'_findopt', and wires the regular series into
chart4 ("fast r10 regular vs shrink vs investigate") as "r10 fast 06a
regular" -- this is the confirmed dataset match for investigate3/shrink's
"07a" batches (same starting logL, verified against every run in the file).
It was excluded from the original build because it's dated the 6th, not
the 7th/8th the original date filter used.

Writes new worksheets and wires them into the workbook (Content_Types,
workbook.xml.rels, workbook.xml <sheets>, and the raw worksheet XML) via
direct zip surgery rather than openpyxl, since any openpyxl resave at this
point would strip queryTables/connections/dxfs and destroy all of the
chart styling fixed by hand in this and earlier scripts.

Only adds two new worksheet parts and touches xl/charts/chart4.xml,
xl/workbook.xml, xl/_rels/workbook.xml.rels, [Content_Types].xml.
"""

import sys
import zipfile

sys.path.insert(0, "test_scripts")
from build_phylogeny_chart import load_series, split_and_average
import restyle_charts as rc

TARGET = "PhylogenyTreeLengths.xlsx"
HEADER = ["run_id", "candidates", "time_elapsed", "logL", "true_minus_current"]
COLS = ["A", "B", "C", "D", "E"]


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def fmt_num(v):
    if isinstance(v, int):
        return str(v)
    return repr(float(v))


def build_row_xml(row_num, values):
    cells = []
    for col, val in zip(COLS, values):
        ref = f"{col}{row_num}"
        if isinstance(val, str):
            cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{esc(val)}</t></is></c>')
        else:
            cells.append(f'<c r="{ref}"><v>{fmt_num(val)}</v></c>')
    return f'<row r="{row_num}" spans="1:5" x14ac:dyDescent="0.3">{"".join(cells)}</row>'


def build_sheet_xml(rows):
    all_rows = [HEADER] + rows
    n = len(all_rows)
    rows_xml = "".join(build_row_xml(i + 1, r) for i, r in enumerate(all_rows))
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" mc:Ignorable="x14ac" '
        'xmlns:x14ac="http://schemas.microsoft.com/office/spreadsheetml/2009/9/ac">'
        f'<dimension ref="A1:E{n}"/>'
        '<sheetViews><sheetView workbookViewId="0"/></sheetViews>'
        '<sheetFormatPr defaultRowHeight="14.5" x14ac:dyDescent="0.35"/>'
        f'<sheetData>{rows_xml}</sheetData>'
        '<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" header="0.3" footer="0.3"/>'
        '</worksheet>'
    )


def build_numcache(values):
    pts = "".join(f'<c:pt idx="{i}"><c:v>{fmt_num(v)}</c:v></c:pt>' for i, v in enumerate(values))
    return f'<c:formatCode>General</c:formatCode><c:ptCount val="{len(values)}"/>{pts}'


def build_series_xml(idx, title, sheet_name, rows, color_key):
    x_values = [r[2] for r in rows]
    y_values = [r[3] for r in rows]
    n = len(rows)
    fill = rc.FILLS[color_key]
    return (
        f'<c:ser><c:idx val="{idx}"/><c:order val="{idx}"/>'
        f'<c:tx><c:v>{esc(title)}</c:v></c:tx>'
        f'<c:spPr><a:ln w="19050"><a:solidFill>{fill}</a:solidFill><a:prstDash val="solid"/></a:ln></c:spPr>'
        '<c:marker><c:symbol val="none"/></c:marker>'
        f"<c:xVal><c:numRef><c:f>'{sheet_name}'!$C$2:$C${n + 1}</c:f><c:numCache>{build_numcache(x_values)}</c:numCache></c:numRef></c:xVal>"
        f"<c:yVal><c:numRef><c:f>'{sheet_name}'!$D$2:$D${n + 1}</c:f><c:numCache>{build_numcache(y_values)}</c:numCache></c:numRef></c:yVal>"
        '<c:smooth val="0"/></c:ser>'
    )


def main():
    series = load_series("yrecord_GTR_FO_fast_findopt.csv")
    r10_runs = sorted(rid for rid in series if rid.startswith("20260806") and "_r10_" in rid)
    assert len(r10_runs) == 10, f"expected 10 matching runs, found {len(r10_runs)}"

    regular_rows, findopt_rows = split_and_average(r10_runs, series, "r10_findopt_06a")
    print(f"regular_rows={len(regular_rows)} findopt_rows={len(findopt_rows)}")

    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    # --- new worksheets ---
    reg_name, fo_name = "r10_findopt_06a_regular", "r10_findopt_06a_findopt"
    entries["xl/worksheets/sheet37.xml"] = build_sheet_xml(regular_rows).encode("utf-8")
    entries["xl/worksheets/sheet38.xml"] = build_sheet_xml(findopt_rows).encode("utf-8")
    order += ["xl/worksheets/sheet37.xml", "xl/worksheets/sheet38.xml"]

    # --- wire into workbook.xml ---
    wb_xml = entries["xl/workbook.xml"].decode("utf-8")
    new_sheets = (
        f'<sheet name="{reg_name}" sheetId="37" r:id="rId42"/>'
        f'<sheet name="{fo_name}" sheetId="38" r:id="rId43"/>'
    )
    assert "</sheets>" in wb_xml
    wb_xml = wb_xml.replace("</sheets>", new_sheets + "</sheets>")
    entries["xl/workbook.xml"] = wb_xml.encode("utf-8")

    # --- wire into workbook.xml.rels ---
    wb_rels = entries["xl/_rels/workbook.xml.rels"].decode("utf-8")
    new_rels = (
        '<Relationship Id="rId42" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet37.xml"/>'
        '<Relationship Id="rId43" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet38.xml"/>'
    )
    assert "</Relationships>" in wb_rels
    wb_rels = wb_rels.replace("</Relationships>", new_rels + "</Relationships>")
    entries["xl/_rels/workbook.xml.rels"] = wb_rels.encode("utf-8")

    # --- wire into [Content_Types].xml ---
    ct = entries["[Content_Types].xml"].decode("utf-8")
    new_ct = (
        '<Override PartName="/xl/worksheets/sheet37.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/worksheets/sheet38.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
    )
    assert "</Types>" in ct
    ct = ct.replace("</Types>", new_ct + "</Types>")
    entries["[Content_Types].xml"] = ct.encode("utf-8")

    # --- add series to chart4.xml ---
    chart4 = entries["xl/charts/chart4.xml"].decode("utf-8")
    n_existing = chart4.count("<c:ser>")
    new_series = build_series_xml(n_existing, "r10 fast 06a regular", reg_name, regular_rows, "accent1")
    marker = "<c:dLbls>"
    idx = chart4.find(marker)
    assert idx != -1
    chart4 = chart4[:idx] + new_series + chart4[idx:]
    entries["xl/charts/chart4.xml"] = chart4.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print(f"Added sheets '{reg_name}' ({len(regular_rows)} rows) and '{fo_name}' ({len(findopt_rows)} rows).")
    print(f"Added chart4 series 'r10 fast 06a regular' (idx={n_existing}, color=accent1).")


if __name__ == "__main__":
    main()
