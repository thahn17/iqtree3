#!/usr/bin/env python3
"""
Follow-up to restyle_charts.py, per user correction:
- "chamfering" means rounded corners (their definition) -- they want SHARP
  edges everywhere, which is the opposite of what restyle_charts.py did
  (it added cap="rnd" + <a:round/> joins to every series, matching what the
  workbook's own pre-existing series happened to use).
- The axes of chart4 should go back to exactly what they were before
  restyle_charts.py touched them (gridlines/title position/label font/the
  axis max values) -- reverting fix_chart4's two <c:valAx> edits verbatim.

This script:
1. Reverts chart4.xml's two <c:valAx> blocks (axId 10 and 20) to their
   pre-restyle_charts.py content, by re-importing the exact old/new strings
   from restyle_charts.py and applying the substitution in reverse. This
   also restores the c:max="200"/c:max="-600000" values that were removed
   as (what looked like) a leftover-contamination bug -- flagged to the
   user separately, not silently dropped.
2. Across chart2.xml, chart3.xml, chart4.xml: strips every `cap="rnd"`
   attribute and every `<a:round/>` line-join element (series lines,
   gridlines, axis lines, chart-area border), and sets
   `<c:roundedCorners val="1"/>` to `val="0"` -- i.e. removes every form of
   corner/edge rounding introduced or already present, leaving default
   sharp (flat cap, miter join, square chart-area corners) throughout.
   Series colors from restyle_charts.py are untouched.

Only touches xl/charts/chart2.xml, chart3.xml, chart4.xml.
"""

import re
import sys
import zipfile

sys.path.insert(0, "test_scripts")
import restyle_charts as rc

TARGET = "PhylogenyTreeLengths.xlsx"


def revert_chart4_axes(xml):
    # x-axis (axId=10): reverse of fix_chart4 step 2
    new_block = (
        '<c:axId val="10"/><c:scaling><c:orientation val="minMax"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="b"/>' + rc.GRIDLINES + rc.X_AXIS_TITLE +
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>' + rc.AXIS_LINE_AND_LABEL_TXPR +
        '<c:crossAx val="20"/>'
    )
    old_block = (
        '<c:axId val="10"/><c:scaling><c:orientation val="minMax"/><c:max val="200"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="b"/><c:majorGridlines/>'
        '<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr/></a:pPr>'
        '<a:r><a:rPr lang="en-US"/><a:t>Time (s)</a:t></a:r></a:p></c:rich></c:tx>'
        '<c:overlay val="1"/></c:title>'
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/><c:crossAx val="20"/>'
    )
    assert new_block in xml, "x-axis styled block not found -- already reverted?"
    xml = xml.replace(new_block, old_block, 1)

    # y-axis (axId=20): reverse of fix_chart4 step 3
    new_block = (
        '<c:axId val="20"/><c:scaling><c:orientation val="minMax"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="l"/>' + rc.GRIDLINES + rc.Y_AXIS_TITLE +
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>' + rc.AXIS_LINE_AND_LABEL_TXPR +
        '<c:crossAx val="10"/>'
    )
    old_block = (
        '<c:axId val="20"/><c:scaling><c:orientation val="minMax"/><c:max val="-600000"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="l"/><c:majorGridlines/>'
        '<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr/></a:pPr>'
        '<a:r><a:rPr lang="en-US"/><a:t>Log Likelihood</a:t></a:r></a:p></c:rich></c:tx>'
        '<c:overlay val="1"/></c:title>'
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/><c:crossAx val="10"/>'
    )
    assert new_block in xml, "y-axis styled block not found -- already reverted?"
    xml = xml.replace(new_block, old_block, 1)

    return xml


def sharpen(xml):
    before_caps = xml.count('cap="rnd"')
    before_joins = xml.count("<a:round/>")
    xml = xml.replace(' cap="rnd"', "")
    xml = xml.replace("<a:round/>", "")
    xml = xml.replace('<c:roundedCorners val="1"/>', '<c:roundedCorners val="0"/>')
    return xml, before_caps, before_joins


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    c4 = entries["xl/charts/chart4.xml"].decode("utf-8")
    c4 = revert_chart4_axes(c4)
    entries["xl/charts/chart4.xml"] = c4.encode("utf-8")

    totals = {"caps": 0, "joins": 0}
    for name in ("xl/charts/chart2.xml", "xl/charts/chart3.xml", "xl/charts/chart4.xml"):
        xml = entries[name].decode("utf-8")
        xml, n_caps, n_joins = sharpen(xml)
        totals["caps"] += n_caps
        totals["joins"] += n_joins
        entries[name] = xml.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print("Reverted chart4 axes to pre-restyle state (gridlines/title-position/label font/"
          "max=200,-600000 all restored).")
    print(f"Stripped {totals['caps']} cap=\"rnd\" attributes and {totals['joins']} <a:round/> "
          f"joins across chart2/3/4; set roundedCorners to 0 on all three.")


if __name__ == "__main__":
    main()
