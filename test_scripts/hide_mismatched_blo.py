#!/usr/bin/env python3
"""
Hides "r10 findopt 07a findopt (BLO)" from chart3's plot -- its dataset
(true_logL=-623642.183) has no complete matching fullreopt series (fullreopt
only has complete 10-run data for the 07b/07c/07d datasets), so it doesn't
pair with anything in this "Standard/BLO/fullreopt" comparison chart.

Its legend entry is already deleted (idx=2 is in the existing
<c:legendEntry><c:delete val="1"/></c:legendEntry> list alongside its
sibling findopt(BLO) series), so this only needs to make the line itself
invisible (noFill) -- a single, low-risk spPr edit, deliberately avoiding
the <c15:filteredScatterSeries> "proper hide" mechanism that has twice
caused Excel to drop the whole drawing this session.

Only touches xl/charts/chart3.xml.
"""

import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"
TITLE = "r10 findopt 07a findopt (BLO)"
OLD_LINE = '<a:ln w="19050"><a:solidFill><a:srgbClr val="FF0000"/></a:solidFill><a:prstDash val="solid"/></a:ln>'
NEW_LINE = '<a:ln w="19050"><a:noFill/><a:prstDash val="solid"/></a:ln>'


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    chart3 = entries["xl/charts/chart3.xml"].decode("utf-8")
    idx = chart3.find(f"<c:tx><c:v>{TITLE}</c:v></c:tx>")
    assert idx != -1, "series not found"
    start = chart3.rfind("<c:ser>", 0, idx)
    end = chart3.find("</c:ser>", idx) + len("</c:ser>")
    block = chart3[start:end]
    assert OLD_LINE in block, "expected line style not found"
    new_block = block.replace(OLD_LINE, NEW_LINE, 1)
    chart3 = chart3[:start] + new_block + chart3[end:]
    entries["xl/charts/chart3.xml"] = chart3.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print(f"Made {TITLE!r}'s line invisible (noFill).")


if __name__ == "__main__":
    main()
