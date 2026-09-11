#!/usr/bin/env python3
"""
Conservative retry of rebuild_fullreopt.py, after two consecutive attempts
both caused Excel to drop drawing2.xml (all three Sheet1 charts) on open.

Both failed attempts combined several operation *types* in one pass, some
of which had never been used successfully before in this session (deleting
a sheet, renaming a sheet's <sheet name>, removing a whole visible <c:ser>,
editing content inside <c15:filteredScatterSeries>). Since either attempt
could have failed for any of several reasons and there's no local OOXML
schema validator to pin it down, this version uses *only* operation types
already proven safe earlier in this session:
  - adding a brand new sheet (used by add_matched_fast_batch.py, fine)
  - overwriting an existing sheet's full row content in place (used by
    restore_sheet1_stack.py on Sheet1, fine)
  - adding a brand new visible <c:ser> (used repeatedly, fine)
  - refreshing an existing visible series' numCache in place (used by
    restore_sheet1_stack.py on chart2/chart3, fine)
  - renumbering visible <c:idx>/<c:order> (used by unfilter_and_recolor.py
    successfully)

It deliberately avoids, this round:
  - deleting any sheet (the old r10_fullreopt_07a_* sheets are left in
    place, just no longer referenced by any visible series)
  - renaming any sheet's <sheet name> (investigate3/shrink's underlying
    sheets keep their old "07a" names; only the chart legend text changes)
  - removing an existing visible series (the old "r10 fullreopt 07a
    regular" series IS removed, since that's unavoidable for correctness,
    but it's the only removal operation performed)
  - touching anything inside <c15:filteredScatterSeries> at all -- the
    hidden findopt(BLO) sub-series for 07a/07b/07c keep their OLD
    (currently-wrong) cached data untouched; if unhidden they'll show
    stale numbers until a future, more carefully isolated fix corrects
    them without risking the whole drawing again

Only touches xl/charts/chart3.xml, chart4.xml, and adds two new worksheet
parts (no sheet deletions or renames).
"""

import re
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


def build_visible_series_xml(idx, title, sheet_name, rows, color_key):
    x_values, y_values = [r[2] for r in rows], [r[3] for r in rows]
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


def get_series_block(xml, title):
    idx = xml.find(f"<c:tx><c:v>{title}</c:v></c:tx>")
    assert idx != -1, f"title not found: {title}"
    start = xml.rfind("<c:ser>", 0, idx)
    end = xml.find("</c:ser>", idx) + len("</c:ser>")
    return start, end, xml[start:end]


def remove_series(xml, title):
    start, end, _ = get_series_block(xml, title)
    return xml[:start] + xml[end:]


def renumber_visible(xml):
    parts, pos, counter = [], 0, 0
    while True:
        start = xml.find("<c:ser>", pos)
        if start == -1:
            parts.append(xml[pos:])
            break
        end = xml.find("</c:ser>", start) + len("</c:ser>")
        block = xml[start:end]
        block = re.sub(r'<c:idx val="\d+"/>', f'<c:idx val="{counter}"/>', block, count=1)
        block = re.sub(r'<c:order val="\d+"/>', f'<c:order val="{counter}"/>', block, count=1)
        parts.append(xml[pos:start])
        parts.append(block)
        pos = end
        counter += 1
    return "".join(parts), counter


def refresh_visible_numcache(xml, title, sheet_name, rows):
    start, end, block = get_series_block(xml, title)
    x_values, y_values = [r[2] for r in rows], [r[3] for r in rows]
    n = len(rows)
    new_block = block
    new_block = re.sub(
        rf"<c:f>'{sheet_name}'!\$C\$2:\$C\$\d+</c:f><c:numCache>.*?</c:numCache>",
        f"<c:f>'{sheet_name}'!$C$2:$C${n + 1}</c:f><c:numCache>{build_numcache(x_values)}</c:numCache>",
        new_block, count=1, flags=re.DOTALL,
    )
    new_block = re.sub(
        rf"<c:f>'{sheet_name}'!\$D\$2:\$D\$\d+</c:f><c:numCache>.*?</c:numCache>",
        f"<c:f>'{sheet_name}'!$D$2:$D${n + 1}</c:f><c:numCache>{build_numcache(y_values)}</c:numCache>",
        new_block, count=1, flags=re.DOTALL,
    )
    assert new_block != block, f"no substitution happened for {title!r}"
    return xml[:start] + new_block + xml[end:]


def rename_series_title_only(xml, old_title, new_title):
    """Changes only the display title text of an existing visible series --
    leaves its <c:f> sheet references (and therefore the underlying sheet
    name) completely untouched, to avoid the sheet-rename operation that
    hasn't been proven safe."""
    start, end, block = get_series_block(xml, old_title)
    new_block = block.replace(f"<c:v>{old_title}</c:v>", f"<c:v>{new_title}</c:v>", 1)
    assert new_block != block
    return xml[:start] + new_block + xml[end:]


def main():
    fr_series = load_series("yrecord_GTR_FO_fast_fullreopt_findopt.csv")
    groups = {"07b": -648090.947, "07c": -649071.173, "07d": -605880.034}
    run_lists = {label: [] for label in groups}
    for rid in sorted(fr_series):
        times, cand, logl, diff = fr_series[rid]
        tl = logl[0] + diff[0]
        for label, target in groups.items():
            if abs(tl - target) < 0.01:
                run_lists[label].append(rid)
    for label, runs in run_lists.items():
        assert len(runs) == 10, f"{label}: expected 10 runs, got {len(runs)}"

    computed = {}
    for label, runs in run_lists.items():
        regular_rows, findopt_rows = split_and_average(runs, fr_series, f"r10_fullreopt_{label}_true")
        computed[label] = (regular_rows, findopt_rows)
        print(f"{label}: regular={len(regular_rows)} findopt={len(findopt_rows)}")

    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    # ---------- overwrite 07b/07c sheet data in place (no rename, no delete) ----------
    wb_xml = entries["xl/workbook.xml"].decode("utf-8")
    wb_rels = entries["xl/_rels/workbook.xml.rels"].decode("utf-8")
    ct = entries["[Content_Types].xml"].decode("utf-8")

    def sheet_filename(name):
        m = re.search(rf'<sheet name="{re.escape(name)}" sheetId="(\d+)" r:id="(rId\d+)"/>', wb_xml)
        assert m, f"sheet not found: {name}"
        rid = m.group(2)
        m2 = re.search(rf'<Relationship Id="{rid}" [^>]*Target="worksheets/(sheet\d+\.xml)"/>', wb_rels)
        assert m2, f"relationship not found for {rid}"
        return m2.group(1)

    for label in ["07b", "07c"]:
        reg_rows, fo_rows = computed[label]
        reg_fname = sheet_filename(f"r10_fullreopt_{label}_regular")
        fo_fname = sheet_filename(f"r10_fullreopt_{label}_findopt")
        entries[f"xl/worksheets/{reg_fname}"] = build_sheet_xml(reg_rows).encode("utf-8")
        entries[f"xl/worksheets/{fo_fname}"] = build_sheet_xml(fo_rows).encode("utf-8")
        print(f"overwrote data for r10_fullreopt_{label}_regular/_findopt "
              f"({len(reg_rows)}/{len(fo_rows)} rows)")

    # ---------- add new 07d sheets ----------
    max_sheet_num = max(int(re.search(r"sheet(\d+)\.xml", n).group(1))
                         for n in order if re.search(r"xl/worksheets/sheet\d+\.xml$", n))
    max_sheet_id = max(int(m) for m in re.findall(r'sheetId="(\d+)"', wb_xml))
    max_rid = max(int(m) for m in re.findall(r'Id="rId(\d+)"', wb_rels))

    reg_rows, fo_rows = computed["07d"]
    new_sheets = [("r10_fullreopt_07d_regular", reg_rows), ("r10_fullreopt_07d_findopt", fo_rows)]
    for i, (name, rows) in enumerate(new_sheets):
        num, sid, rid = max_sheet_num + 1 + i, max_sheet_id + 1 + i, max_rid + 1 + i
        fname = f"sheet{num}.xml"
        entries[f"xl/worksheets/{fname}"] = build_sheet_xml(rows).encode("utf-8")
        order.append(f"xl/worksheets/{fname}")
        wb_xml = wb_xml.replace("</sheets>", f'<sheet name="{name}" sheetId="{sid}" r:id="rId{rid}"/></sheets>')
        wb_rels = wb_rels.replace(
            "</Relationships>",
            f'<Relationship Id="rId{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/{fname}"/></Relationships>',
        )
        ct = ct.replace(
            "</Types>",
            f'<Override PartName="/xl/worksheets/{fname}" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>',
        )
        print(f"added sheet {name} ({fname}, sheetId={sid}, rId=rId{rid}) with {len(rows)} rows")

    entries["xl/workbook.xml"] = wb_xml.encode("utf-8")
    entries["xl/_rels/workbook.xml.rels"] = wb_rels.encode("utf-8")
    entries["[Content_Types].xml"] = ct.encode("utf-8")

    # ---------- chart3.xml: visible series only, no c15 edits ----------
    chart3 = entries["xl/charts/chart3.xml"].decode("utf-8")
    chart3 = remove_series(chart3, "r10 fullreopt 07a regular")
    for label in ["07b", "07c"]:
        reg_rows, _ = computed[label]
        chart3 = refresh_visible_numcache(chart3, f"r10 fullreopt {label} regular",
                                           f"r10_fullreopt_{label}_regular", reg_rows)

    n_before = chart3.count("<c:ser>")
    reg_rows, fo_rows = computed["07d"]
    new_reg = build_visible_series_xml(n_before, "r10 fullreopt 07d regular",
                                        "r10_fullreopt_07d_regular", reg_rows, "accent5")
    marker = "<c:dLbls>"
    mi = chart3.find(marker)
    assert mi != -1
    chart3 = chart3[:mi] + new_reg + chart3[mi:]

    n_before2 = chart3.count("<c:ser>")
    new_fo = build_visible_series_xml(n_before2, "r10 fullreopt 07d findopt (BLO)",
                                       "r10_fullreopt_07d_findopt", fo_rows, "red")
    mi2 = chart3.find(marker)
    chart3 = chart3[:mi2] + new_fo + chart3[mi2:]

    chart3, n_total = renumber_visible(chart3)
    entries["xl/charts/chart3.xml"] = chart3.encode("utf-8")
    print(f"chart3: removed 07a regular, refreshed 07b/07c, added 07d regular + findopt (BLO) "
          f"(visible total={n_total}); hidden c15 content untouched")

    # ---------- chart4.xml: title-only rename (sheet refs untouched) ----------
    chart4 = entries["xl/charts/chart4.xml"].decode("utf-8")
    chart4 = rename_series_title_only(chart4, "r10 investigate3 07a", "r10 investigate3 06a")
    chart4 = rename_series_title_only(chart4, "r10 shrink500 07a", "r10 shrink500 06a")
    entries["xl/charts/chart4.xml"] = chart4.encode("utf-8")
    print("chart4: renamed investigate3/shrink 07a titles to 06a (sheet references unchanged)")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print("Done.")


if __name__ == "__main__":
    main()
