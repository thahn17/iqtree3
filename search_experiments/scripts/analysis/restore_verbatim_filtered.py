#!/usr/bin/env python3
"""
Follow-up fix: the hand-crafted <c15:filteredScatterSeries> wrapper written
by rehide_series.py apparently violates the schema somewhere subtle enough
that Excel dropped the whole drawing2.xml part (all three Sheet1 charts)
rather than flag one specific series. Rather than keep guessing at the
exact schema requirements for this Microsoft chart-filter extension,
this replaces that hand-built wrapper with the byte-for-byte original
extracted from PhylogenyTreeLengths.presurgery6.xlsx (the last version
confirmed to open in Excel without error) for all three charts.

This restores the "hidden for a reason" 14 series using Excel's own
previously-validated bytes. One known side effect: chart4's copy of this
wrapper has "r10 fast 07a regular" stored with Excel's auto-assigned
accent2 color rather than the accent1 the rest of that series family uses
-- cosmetic only, invisible unless someone un-hides it via the chart
filter, and not touched here to avoid reintroducing a hand-edit into this
proven-good block.

Only touches xl/charts/chart2.xml, chart3.xml, chart4.xml.
"""

import re
import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"
GOOD_SOURCE = "PhylogenyTreeLengths.presurgery6.xlsx"

EXT_URI_ESC = r"\{02D57815-91ED-43cb-92C2-25804820EDAC\}"
WRAPPER_RE = re.compile(
    rf'<c:extLst><c:ext uri="{EXT_URI_ESC}"[^>]*><c15:filteredScatterSeries>.*?'
    r"</c15:filteredScatterSeries></c:ext></c:extLst>",
    re.DOTALL,
)


def main():
    with zipfile.ZipFile(GOOD_SOURCE) as z:
        good = {name: WRAPPER_RE.search(z.read(f"xl/charts/{name}.xml").decode("utf-8")).group(0)
                for name in ["chart2", "chart3", "chart4"]}

    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    for name in ["chart2", "chart3", "chart4"]:
        path = f"xl/charts/{name}.xml"
        xml = entries[path].decode("utf-8")
        m = WRAPPER_RE.search(xml)
        assert m, f"no current wrapper found in {path}"
        xml = xml[:m.start()] + good[name] + xml[m.end():]
        entries[path] = xml.encode("utf-8")
        print(f"{path}: replaced hand-built wrapper with verbatim original ({len(good[name])} bytes).")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for n in order:
            z.writestr(n, entries[n])


if __name__ == "__main__":
    main()
