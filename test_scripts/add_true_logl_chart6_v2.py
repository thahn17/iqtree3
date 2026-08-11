#!/usr/bin/env python3
"""
Retry of add_true_logl_chart6.py. The first attempt used <c:numLit> for
the 3 reference lines' inline 2-point data, which turned out to be an
untested pattern nowhere else in this workbook (every other successful
series addition this session used <c:numRef> pointing at a real sheet
range). That combined with a later over-broad dash edit both preceded
Excel dropping drawing2.xml -- rather than keep guessing which one broke
it, this redo eliminates the untested numLit path entirely and sticks to
the proven-safe pattern: a tiny new worksheet per reference line + numRef.

Adds 3 small worksheets (2 data rows each: one point at tmin, one at tmax,
both y = that dataset's true logL) and wires them into chart6.xml only, as
3 new visible series with a dashed gray line. Also (re)applies dash only
to chart6's other 6 data-series lines, scoped precisely to each <c:ser>
block so gridlines/axis/border lines are never touched.

Only touches: adds 3 new worksheet parts, and edits xl/charts/chart6.xml,
xl/workbook.xml, xl/_rels/workbook.xml.rels, [Content_Types].xml. No other
chart file is read or written.
"""

import re
import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"
HEADER = ["label", "unused_b", "time_elapsed", "true_logL", "unused_e"]
COLS = ["A", "B", "C", "D", "E"]

DATASETS = {
    "07a": (-623642.1826232, 3.33, 147.37),
    "07c": (-649071.173, 3.25, 163.14),
    "07d": (-605880.034, 3.30, 153.60),
}

MAIN_SERIES_TITLES = [
    "r=10", "r10 findopt 07c regular", "r10 findopt 07d regular",
    "r=1", "r1 findopt 07c regular", "r1 findopt 07d regular",
]


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


def build_ref_line_series(idx, label, sheet_name, true_logl, tmin, tmax):
    return (
        f'<c:ser><c:idx val="{idx}"/><c:order val="{idx}"/>'
        f'<c:tx><c:v>true logL {label}</c:v></c:tx>'
        '<c:spPr><a:ln w="12700"><a:solidFill><a:schemeClr val="tx1"><a:lumMod val="50000"/>'
        '<a:lumOff val="50000"/></a:schemeClr></a:solidFill><a:prstDash val="dash"/></a:ln></c:spPr>'
        '<c:marker><c:symbol val="none"/></c:marker>'
        f"<c:xVal><c:numRef><c:f>'{sheet_name}'!$C$2:$C$3</c:f><c:numCache>{build_numcache([tmin, tmax])}</c:numCache></c:numRef></c:xVal>"
        f"<c:yVal><c:numRef><c:f>'{sheet_name}'!$D$2:$D$3</c:f><c:numCache>{build_numcache([true_logl, true_logl])}</c:numCache></c:numRef></c:yVal>"
        '<c:smooth val="0"/></c:ser>'
    )


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    chart6 = entries["xl/charts/chart6.xml"].decode("utf-8")
    n_before = chart6.count("<c:ser>")
    assert n_before == 6, f"expected 6 existing series, found {n_before}"

    wb_xml = entries["xl/workbook.xml"].decode("utf-8")
    wb_rels = entries["xl/_rels/workbook.xml.rels"].decode("utf-8")
    ct = entries["[Content_Types].xml"].decode("utf-8")

    max_sheet_num = max(int(re.search(r"sheet(\d+)\.xml", n).group(1))
                         for n in order if re.search(r"xl/worksheets/sheet\d+\.xml$", n))
    max_sheet_id = max(int(m) for m in re.findall(r'sheetId="(\d+)"', wb_xml))
    max_rid = max(int(m) for m in re.findall(r'Id="rId(\d+)"', wb_rels))

    new_series_xml = []
    for i, (label, (true_logl, tmin, tmax)) in enumerate(DATASETS.items()):
        sheet_name = f"true_logL_{label}"
        num, sid, rid = max_sheet_num + 1 + i, max_sheet_id + 1 + i, max_rid + 1 + i
        fname = f"sheet{num}.xml"
        rows = [["min", 0, tmin, true_logl, 0], ["max", 0, tmax, true_logl, 0]]
        entries[f"xl/worksheets/{fname}"] = build_sheet_xml(rows).encode("utf-8")
        order.append(f"xl/worksheets/{fname}")
        wb_xml = wb_xml.replace("</sheets>", f'<sheet name="{sheet_name}" sheetId="{sid}" r:id="rId{rid}"/></sheets>')
        wb_rels = wb_rels.replace(
            "</Relationships>",
            f'<Relationship Id="rId{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/{fname}"/></Relationships>',
        )
        ct = ct.replace(
            "</Types>",
            f'<Override PartName="/xl/worksheets/{fname}" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>',
        )
        new_series_xml.append(build_ref_line_series(n_before + i, label, sheet_name, true_logl, tmin, tmax))
        print(f"added sheet {sheet_name} ({fname}, sheetId={sid}, rId=rId{rid})")

    entries["xl/workbook.xml"] = wb_xml.encode("utf-8")
    entries["xl/_rels/workbook.xml.rels"] = wb_rels.encode("utf-8")
    entries["[Content_Types].xml"] = ct.encode("utf-8")

    marker = "<c:dLbls>"
    mi = chart6.find(marker)
    assert mi != -1
    chart6 = chart6[:mi] + "".join(new_series_xml) + chart6[mi:]

    # precisely dash the 6 main data-series lines, scoped per <c:ser> block
    def series_blocks(xml):
        parts, pos = [], 0
        while True:
            start = xml.find("<c:ser>", pos)
            if start == -1:
                break
            end = xml.find("</c:ser>", start) + len("</c:ser>")
            parts.append((start, end, xml[start:end]))
            pos = end
        return parts

    parts, pos, n_dashed = [], 0, 0
    for start, end, block in series_blocks(chart6):
        t = re.search(r"<c:tx><c:v>(.*?)</c:v></c:tx>", block).group(1)
        parts.append(chart6[pos:start])
        if t in MAIN_SERIES_TITLES:
            new_block = block.replace('<a:prstDash val="solid"/>', '<a:prstDash val="dash"/>', 1)
            assert new_block != block, f"no dash substitution for {t!r}"
            n_dashed += 1
            parts.append(new_block)
        else:
            parts.append(block)
        pos = end
    parts.append(chart6[pos:])
    chart6 = "".join(parts)
    assert n_dashed == len(MAIN_SERIES_TITLES), f"expected {len(MAIN_SERIES_TITLES)} dashed, got {n_dashed}"

    entries["xl/charts/chart6.xml"] = chart6.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print(f"Added {len(DATASETS)} true-logL reference lines (numRef-based) and dashed "
          f"{n_dashed} main series. Total chart6 series now {n_before + len(DATASETS)}.")


if __name__ == "__main__":
    main()
