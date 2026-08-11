#!/usr/bin/env python3
"""
Rebuilds the fullreopt data in chart3 ("Standard/BLO/fullreopt" comparison)
using correct dataset groupings, per the discovery that fullreopt's runs
don't line up with naive "chunks of 10 in file order" batching -- the
actual dataset transitions happen mid-chunk (verified via true_logL =
logL + true_minus_current, which is invariant per dataset and consistent
even across different algorithms/time-points, unlike raw starting logL).

True dataset groups found (10 clean runs each):
- true_logL=-648090.947 -> matches findopt/investigate3/shrink's "07b"
- true_logL=-649071.173 -> matches findopt/investigate3/shrink's "07c"
- true_logL=-605880.034 -> matches findopt/investigate3/shrink's "07d"
- true_logL=-623642.183 (findopt's "07a") has only 3 fullreopt runs
  (incomplete) -- no complete fullreopt series exists for this dataset,
  so the old (mixed-data) "r10 fullreopt 07a" series is removed outright.

Rebuilds:
- deletes sheets r10_fullreopt_07a_regular/_findopt (no longer valid)
- overwrites r10_fullreopt_07b_regular/_findopt and 07c_regular/_findopt
  with the correct 10-run data for those datasets (same sheet names, same
  chart series names -- only the underlying data was wrong)
- adds new sheets + chart series for the previously-missing 07d dataset
- renames r10_investigate3_07a_regular -> ..._06a_regular and
  r10_shrink_07a_regular -> ..._06a_regular, since those already-hidden
  series correspond to the "06"-dated fast dataset, not fullreopt's or
  findopt's own "07a" -- and updates chart4's hidden series titles/formula
  refs to match

Preserves existing look: same colors (accent5 for fullreopt regular, red
for findopt/BLO sub-series), same hidden/visible split (fullreopt findopt
sub-series stay hidden, matching the existing pattern), sharp edges.
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


def get_c15_ser_block(xml, title):
    pos = 0
    while True:
        start = xml.find("<c15:ser>", pos)
        if start == -1:
            return None, None, None
        end = xml.find("</c15:ser>", start) + len("</c15:ser>")
        block = xml[start:end]
        if f"<c:v>{title}</c:v>" in block:
            return start, end, block
        pos = end


def remove_c15_series(xml, title):
    start, end, block = get_c15_ser_block(xml, title)
    assert block is not None, f"hidden series not found: {title}"
    return xml[:start] + xml[end:]


def refresh_hidden_numcache(xml, title, sheet_name, rows):
    start, end, block = get_c15_ser_block(xml, title)
    assert block is not None, f"hidden series not found: {title}"
    x_values, y_values = [r[2] for r in rows], [r[3] for r in rows]
    n = len(rows)
    new_block = block
    new_block = re.sub(
        rf"(<c15:sqref>'{sheet_name}'!\$C\$2:\$C\$)\d+(</c15:sqref>.*?<c:numCache>).*?(</c:numCache>)",
        lambda m: m.group(1) + str(n + 1) + m.group(2) + build_numcache(x_values) + m.group(3),
        new_block, count=1, flags=re.DOTALL,
    )
    new_block = re.sub(
        rf"(<c15:sqref>'{sheet_name}'!\$D\$2:\$D\$)\d+(</c15:sqref>.*?<c:numCache>).*?(</c:numCache>)",
        lambda m: m.group(1) + str(n + 1) + m.group(2) + build_numcache(y_values) + m.group(3),
        new_block, count=1, flags=re.DOTALL,
    )
    assert new_block != block, f"no substitution happened for hidden {title!r}"
    return xml[:start] + new_block + xml[end:]


def build_hidden_series_from_template(template_block, old_title, new_title, old_sheet, new_sheet, rows):
    """Builds a new <c15:ser> block by copying an existing, proven-valid one
    and substituting title/sheet/data -- avoids hand-authoring the tricky
    c15 extension structure from scratch."""
    x_values, y_values = [r[2] for r in rows], [r[3] for r in rows]
    n = len(rows)
    block = template_block.replace(f"<c:v>{old_title}</c:v>", f"<c:v>{new_title}</c:v>")
    block = block.replace(f"'{old_sheet}'", f"'{new_sheet}'")
    block = re.sub(
        r"(\$C\$2:\$C\$)\d+(</c15:sqref>.*?<c:numCache>).*?(</c:numCache>)",
        lambda m: m.group(1) + str(n + 1) + m.group(2) + build_numcache(x_values) + m.group(3),
        block, count=1, flags=re.DOTALL,
    )
    block = re.sub(
        r"(\$D\$2:\$D\$)\d+(</c15:sqref>.*?<c:numCache>).*?(</c:numCache>)",
        lambda m: m.group(1) + str(n + 1) + m.group(2) + build_numcache(y_values) + m.group(3),
        block, count=1, flags=re.DOTALL,
    )
    # fresh idx (doesn't need to be globally unique, just distinct-looking) and uniqueId
    block = re.sub(r'<c:idx val="\d+"/>', '<c:idx val="99"/>', block, count=1)
    block = re.sub(r'<c:order val="\d+"/>', '<c:order val="98"/>', block, count=1)
    block = re.sub(r'<c16:uniqueId val="\{[^}]*\}"/>',
                    '<c16:uniqueId val="{99999999-0000-0000-0000-000000000000}"/>', block)
    return block


def insert_c15_series(xml, new_block):
    marker = "</c15:filteredScatterSeries>"
    idx = xml.find(marker)
    assert idx != -1
    return xml[:idx] + new_block + xml[idx:]


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

    # ---------- worksheets ----------
    # delete 07a sheets
    for name in ["r10_fullreopt_07a_regular", "r10_fullreopt_07a_findopt"]:
        pass  # handled via workbook.xml removal below; file removal happens after we know rIds

    wb_xml = entries["xl/workbook.xml"].decode("utf-8")
    wb_rels = entries["xl/_rels/workbook.xml.rels"].decode("utf-8")
    ct = entries["[Content_Types].xml"].decode("utf-8")

    def sheet_file_for(name):
        m = re.search(rf'<sheet name="{re.escape(name)}" sheetId="(\d+)" r:id="(rId\d+)"/>', wb_xml)
        assert m, f"sheet not found in workbook.xml: {name}"
        rid = m.group(2)
        m2 = re.search(rf'<Relationship Id="{rid}" [^>]*Target="worksheets/(sheet\d+\.xml)"/>', wb_rels)
        assert m2, f"relationship not found for {rid}"
        return m.group(1), rid, m2.group(1)  # sheetId, rId, filename

    # --- delete 07a sheets entirely ---
    for name in ["r10_fullreopt_07a_regular", "r10_fullreopt_07a_findopt"]:
        sheet_id, rid, fname = sheet_file_for(name)
        wb_xml = re.sub(rf'<sheet name="{re.escape(name)}" sheetId="{sheet_id}" r:id="{rid}"/>', "", wb_xml)
        wb_rels = re.sub(rf'<Relationship Id="{rid}" [^>]*Target="worksheets/{fname}"/>', "", wb_rels)
        ct = re.sub(rf'<Override PartName="/xl/worksheets/{fname}"[^>]*/>', "", ct)
        del entries[f"xl/worksheets/{fname}"]
        order.remove(f"xl/worksheets/{fname}")
        print(f"deleted sheet {name} ({fname})")

    # --- overwrite 07b/07c data in place ---
    for label, sheet_kind in [("07b", "regular"), ("07b", "findopt"), ("07c", "regular"), ("07c", "findopt")]:
        pass
    for label in ["07b", "07c"]:
        reg_rows, fo_rows = computed[label]
        _, _, reg_fname = sheet_file_for(f"r10_fullreopt_{label}_regular")
        _, _, fo_fname = sheet_file_for(f"r10_fullreopt_{label}_findopt")
        entries[f"xl/worksheets/{reg_fname}"] = build_sheet_xml(reg_rows).encode("utf-8")
        entries[f"xl/worksheets/{fo_fname}"] = build_sheet_xml(fo_rows).encode("utf-8")
        print(f"overwrote data for r10_fullreopt_{label}_regular/_findopt "
              f"({len(reg_rows)}/{len(fo_rows)} rows)")

    # --- add new 07d sheets ---
    max_sheet_num = max(int(re.search(r"sheet(\d+)\.xml", n).group(1))
                         for n in order if re.search(r"xl/worksheets/sheet\d+\.xml$", n))
    max_sheet_id = max(int(m) for m in re.findall(r'sheetId="(\d+)"', wb_xml))
    max_rid = max(int(m) for m in re.findall(r'Id="rId(\d+)"', wb_rels))

    reg_rows, fo_rows = computed["07d"]
    new_sheets = [
        ("r10_fullreopt_07d_regular", reg_rows),
        ("r10_fullreopt_07d_findopt", fo_rows),
    ]
    new_series_visible_xml = None
    for i, (name, rows) in enumerate(new_sheets):
        num = max_sheet_num + 1 + i
        sid = max_sheet_id + 1 + i
        rid = max_rid + 1 + i
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

    # --- rename investigate3/shrink 07a sheets to 06a ---
    rename_map = {
        "r10_investigate3_07a_regular": "r10_investigate3_06a_regular",
        "r10_shrink_07a_regular": "r10_shrink_06a_regular",
    }
    for old, new in rename_map.items():
        wb_xml = wb_xml.replace(f'<sheet name="{old}" ', f'<sheet name="{new}" ')

    entries["xl/workbook.xml"] = wb_xml.encode("utf-8")
    entries["xl/_rels/workbook.xml.rels"] = wb_rels.encode("utf-8")
    entries["[Content_Types].xml"] = ct.encode("utf-8")

    # ---------- chart3.xml ----------
    chart3 = entries["xl/charts/chart3.xml"].decode("utf-8")

    # remove 07a visible + hidden
    chart3 = remove_series(chart3, "r10 fullreopt 07a regular")
    chart3 = remove_c15_series(chart3, "r10 fullreopt 07a findopt (BLO)")

    # refresh 07b/07c visible + hidden with correct data
    for label in ["07b", "07c"]:
        reg_rows, fo_rows = computed[label]
        chart3 = refresh_visible_numcache(chart3, f"r10 fullreopt {label} regular",
                                           f"r10_fullreopt_{label}_regular", reg_rows)
        chart3 = refresh_hidden_numcache(chart3, f"r10 fullreopt {label} findopt (BLO)",
                                          f"r10_fullreopt_{label}_findopt", fo_rows)

    # add new visible 07d regular series
    marker = "<c:dLbls>"
    mi = chart3.find(marker)
    assert mi != -1
    n_visible_before = chart3.count("<c:ser>")
    new_series = build_visible_series_xml(n_visible_before, "r10 fullreopt 07d regular",
                                           "r10_fullreopt_07d_regular", computed["07d"][0], "accent5")
    chart3 = chart3[:mi] + new_series + chart3[mi:]

    # NOTE: the new 07d findopt (BLO) series is added VISIBLE, not hidden --
    # hand-constructing a new <c15:filteredScatterSeries> entry from a
    # template has twice now caused Excel to drop the entire drawing (all
    # three Sheet1 charts) on open, for reasons that couldn't be pinned down
    # without a real schema validator. Adding it as a normal visible series
    # is the proven-safe path; it can be hidden manually via Excel's own
    # chart filter afterward if desired.
    n_visible_before2 = chart3.count("<c:ser>")
    mi2 = chart3.find(marker)
    new_series2 = build_visible_series_xml(n_visible_before2, "r10 fullreopt 07d findopt (BLO)",
                                            "r10_fullreopt_07d_findopt", computed["07d"][1], "red")
    chart3 = chart3[:mi2] + new_series2 + chart3[mi2:]

    chart3, n_total = renumber_visible(chart3)

    entries["xl/charts/chart3.xml"] = chart3.encode("utf-8")
    print(f"chart3: removed 07a (visible+hidden), refreshed 07b/07c, added 07d "
          f"(visible total={n_total}, hidden total={chart3.count('<c15:ser>')})")

    # ---------- chart4.xml: rename investigate3/shrink 07a -> 06a ----------
    # NOTE: as of this run, these two are VISIBLE (<c:ser>), not hidden --
    # their hidden/visible status apparently changed since they were last
    # checked (likely from the file being opened/interacted with in Excel
    # between turns); "r10 fast 07a regular" is still hidden, untouched here.
    chart4 = entries["xl/charts/chart4.xml"].decode("utf-8")
    for old_title, new_title, old_sheet, new_sheet in [
        ("r10 investigate3 07a", "r10 investigate3 06a",
         "r10_investigate3_07a_regular", "r10_investigate3_06a_regular"),
        ("r10 shrink500 07a", "r10 shrink500 06a",
         "r10_shrink_07a_regular", "r10_shrink_06a_regular"),
    ]:
        start, end, block = get_series_block(chart4, old_title)
        new_block = block.replace(f"<c:v>{old_title}</c:v>", f"<c:v>{new_title}</c:v>")
        new_block = new_block.replace(f"'{old_sheet}'", f"'{new_sheet}'")
        chart4 = chart4[:start] + new_block + chart4[end:]
        print(f"chart4: renamed visible series {old_title!r} -> {new_title!r}")
    entries["xl/charts/chart4.xml"] = chart4.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print("Done.")


if __name__ == "__main__":
    main()
