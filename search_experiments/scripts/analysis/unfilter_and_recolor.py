#!/usr/bin/env python3
"""
One-off repair #6: found that 14 series across the three charts (8 in
chart2, 3 in chart3, 3 in chart4) had been silently demoted into Excel's
"filtered series" representation (`<c15:filteredScatterSeries><c15:ser>...`
inside a `<c:extLst>`, with the series' `<c:f>` range formula replaced by a
`c15:formulaRef`/`c15:sqref` extension) -- Excel's mechanism for a
series that's still defined but not currently plotted/shown in the legend.
No data was lost (the numCache is intact), but the series don't render.

This was found while checking chart4 before adding a new series: 12 series
were expected but only 9 real <c:ser> elements existed; the missing 3 were
exactly the "07a" batch across fast/investigate3/shrink. In chart2 it was
every "findopt" sub-series; in chart3 every "fullreopt findopt (BLO)"
sub-series. Root cause of why Excel filtered these isn't fully certain, but
the fix is mechanical and lossless: unwrap each c15:ser back into a normal
c:ser (reconstructing `<c:f>` from the c15:sqref range text, keeping the
existing numCache as-is), then force-reapply this workbook's per-type color
scheme (from restyle_charts.py's classify functions, extended to also cover
the 5 native series) and the sharp-edges styling to every series uniformly,
since at least one promoted series (chart4's "r10 fast 07a regular") had
also lost its assigned color along the way.

Only touches xl/charts/chart2.xml, chart3.xml, chart4.xml.
"""

import re
import sys
import zipfile

sys.path.insert(0, "test_scripts")
import restyle_charts as rc

TARGET = "PhylogenyTreeLengths.xlsx"

FILLS = rc.FILLS


def sharp_sppr(color_key):
    return f'<c:spPr><a:ln w="19050"><a:solidFill>{FILLS[color_key]}</a:solidFill><a:prstDash val="solid"/></a:ln></c:spPr>'


def classify_chart2_all(title):
    if title == "r=10":
        return "accent1"
    if title == "r=1":
        return "accent2"
    return rc.classify_chart2(title)


def classify_chart3_all(title):
    if title == "Standard":
        return "accent1"
    if title == "BLO":
        return "red"
    if title == "fullreopt":
        return "accent5"
    return rc.classify_chart3(title)


def classify_chart4_all(title):
    return rc.classify_chart4(title)


def unwrap_filtered(xml):
    """Remove the <c:extLst> wrapping <c15:filteredScatterSeries>, returning
    (xml_without_it, list_of_converted_c_ser_blocks)."""
    m = re.search(
        r"<c:extLst><c:ext uri=\"\{02D57815-91ED-43cb-92C2-25804820EDAC\}\"[^>]*>"
        r"<c15:filteredScatterSeries>(.*?)</c15:filteredScatterSeries>"
        r"</c:ext></c:extLst>",
        xml, re.DOTALL,
    )
    if not m:
        return xml, []

    inner = m.group(1)
    xml_without = xml[:m.start()] + xml[m.end():]

    blocks = []
    pos = 0
    while True:
        start = inner.find("<c15:ser>", pos)
        if start == -1:
            break
        end = inner.find("</c15:ser>", start) + len("</c15:ser>")
        blocks.append(inner[start:end])
        pos = end

    converted = []
    for block in blocks:
        block = block.replace("<c15:ser>", "<c:ser>").replace("</c15:ser>", "</c:ser>")
        # replace each numRef's c15:formulaRef/c15:sqref extLst with a plain <c:f>
        # (attributes on the intervening extLst/ext tags vary -- some redeclare
        # xmlns:c15 inline, some don't -- so match loosely on tag names only)
        def repl(mm):
            return f"<c:f>{mm.group(1)}</c:f>"
        block, n_sub = re.subn(
            r"<c:extLst[^>]*><c:ext[^>]*>"
            r"<c15:formulaRef><c15:sqref>(.*?)</c15:sqref></c15:formulaRef></c:ext></c:extLst>",
            repl, block,
        )
        assert n_sub == 2, f"expected 2 formulaRef substitutions (xVal+yVal), got {n_sub}"
        assert "c15:" not in block, "leftover c15: content after conversion"
        converted.append(block)
    return xml_without, converted


def insert_series(xml, new_blocks):
    if not new_blocks:
        return xml
    marker = "<c:dLbls>"
    idx = xml.find(marker)
    assert idx != -1
    return xml[:idx] + "".join(new_blocks) + xml[idx:]


def renumber_and_recolor(xml, classify):
    parts = []
    pos = 0
    counter = 0
    while True:
        start = xml.find("<c:ser>", pos)
        if start == -1:
            parts.append(xml[pos:])
            break
        end = xml.find("</c:ser>", start) + len("</c:ser>")
        block = xml[start:end]

        block = re.sub(r'<c:idx val="\d+"/>', f'<c:idx val="{counter}"/>', block, count=1)
        block = re.sub(r'<c:order val="\d+"/>', f'<c:order val="{counter}"/>', block, count=1)

        t = re.search(r"<c:tx><c:v>(.*?)</c:v></c:tx>", block)
        title = t.group(1) if t else ""
        color = classify(title)
        if color:
            block = re.sub(r"<c:spPr>.*?</c:spPr>", sharp_sppr(color), block, count=1, flags=re.DOTALL)

        parts.append(xml[pos:start])
        parts.append(block)
        pos = end
        counter += 1
    return "".join(parts), counter


def process(xml, classify):
    xml, restored_blocks = unwrap_filtered(xml)
    xml = insert_series(xml, restored_blocks)
    xml, n_total = renumber_and_recolor(xml, classify)
    return xml, len(restored_blocks), n_total


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    for name, classify in [
        ("xl/charts/chart2.xml", classify_chart2_all),
        ("xl/charts/chart3.xml", classify_chart3_all),
        ("xl/charts/chart4.xml", classify_chart4_all),
    ]:
        xml = entries[name].decode("utf-8")
        xml, n_restored, n_total = process(xml, classify)
        entries[name] = xml.encode("utf-8")
        print(f"{name}: restored {n_restored} filtered series, {n_total} total series now visible, recolored.")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])


if __name__ == "__main__":
    main()
