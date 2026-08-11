#!/usr/bin/env python3
"""
Reverts unfilter_and_recolor.py's un-hiding: the user hid these 14 series
deliberately (via Excel's chart filter) and wants that preserved. Puts them
back into the c15:filteredScatterSeries representation, using their
CURRENT (correctly colored, sharp-edged, correct-data) content rather than
reverting to the pre-fix state, so if they're ever unhidden again they'll
display correctly.

Leaves everything else alone: the newly added "r10 fast 06a regular"
series, colors, sharp-edge styling, all data.

Only touches xl/charts/chart2.xml, chart3.xml, chart4.xml.
"""

import re
import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"

TO_HIDE = {
    "xl/charts/chart2.xml": [
        "r10 findopt 07a findopt", "r10 findopt 07b findopt",
        "r10 findopt 07c findopt", "r10 findopt 07d findopt",
        "r1 findopt 07a findopt", "r1 findopt 07b findopt",
        "r1 findopt 07c findopt", "r1 findopt 07d findopt",
    ],
    "xl/charts/chart3.xml": [
        "r10 fullreopt 07a findopt (BLO)",
        "r10 fullreopt 07b findopt (BLO)",
        "r10 fullreopt 07c findopt (BLO)",
    ],
    "xl/charts/chart4.xml": [
        "r10 fast 07a regular", "r10 investigate3 07a", "r10 shrink500 07a",
    ],
}

EXT_URI = "{02D57815-91ED-43cb-92C2-25804820EDAC}"


def extract_series(xml, titles):
    """Removes each <c:ser> whose title is in `titles` from xml, returns
    (xml_without_them, {title: block})."""
    wanted = set(titles)
    found = {}
    parts = []
    pos = 0
    while True:
        start = xml.find("<c:ser>", pos)
        if start == -1:
            parts.append(xml[pos:])
            break
        end = xml.find("</c:ser>", start) + len("</c:ser>")
        block = xml[start:end]
        t = re.search(r"<c:tx><c:v>(.*?)</c:v></c:tx>", block).group(1)
        if t in wanted:
            found[t] = block
            parts.append(xml[pos:start])
        else:
            parts.append(xml[pos:end])
        pos = end
    assert set(found) == wanted, f"missing titles: {wanted - set(found)}"
    return "".join(parts), found


def to_filtered(block):
    block = block.replace("<c:ser>", "<c15:ser>").replace("</c:ser>", "</c15:ser>")

    def repl(m):
        return (f'<c:extLst><c:ext uri="{EXT_URI}">'
                f"<c15:formulaRef><c15:sqref>{m.group(1)}</c15:sqref></c15:formulaRef>"
                f"</c:ext></c:extLst>")
    block, n = re.subn(r"<c:f>(.*?)</c:f>", repl, block)
    assert n == 2, f"expected 2 <c:f> replacements, got {n}"
    return block


def renumber(xml):
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
    return "".join(parts)


def insert_filtered_block(xml, filtered_ser_blocks):
    wrapper = (
        f'<c:extLst><c:ext uri="{EXT_URI}" xmlns:c15="http://schemas.microsoft.com/office/drawing/2012/chart">'
        f'<c15:filteredScatterSeries>{"".join(filtered_ser_blocks)}</c15:filteredScatterSeries>'
        f"</c:ext></c:extLst>"
    )
    marker = "</c:scatterChart>"
    idx = xml.find(marker)
    assert idx != -1
    return xml[:idx] + wrapper + xml[idx:]


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    for path, titles in TO_HIDE.items():
        xml = entries[path].decode("utf-8")
        xml, found = extract_series(xml, titles)
        xml = renumber(xml)
        filtered_blocks = [to_filtered(found[t]) for t in titles]
        xml = insert_filtered_block(xml, filtered_blocks)
        entries[path] = xml.encode("utf-8")
        print(f"{path}: re-hid {len(titles)} series.")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])


if __name__ == "__main__":
    main()
