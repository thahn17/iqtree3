#!/usr/bin/env python3
"""
One-off repair #3: style the new/added chart series to match the workbook's
existing hand-authored charts, and fix a leftover bug on the third chart's
axes.

Findings that motivate this:
- The two pre-existing charts ("r=10 vs r=1" / chart2.xml and
  "Standard/BLO/fullreopt" / chart3.xml) style their own native series with
  an explicit color, round line cap, and round line join
  (`<a:ln w="19050" cap="rnd"><a:solidFill>...</a:solidFill>
  <a:prstDash val="solid"/><a:round/></a:ln>`). Every series added by
  build_phylogeny_chart.py (in chart2, chart3, and the brand new chart4)
  instead has a bare `<a:ln w="19050"><a:prstDash val="solid"/></a:ln>` --
  no cap, no join, no color -- which renders with sharp/mitered corners at
  every direction change ("chamfered" look) and Excel's arbitrary auto
  color cycling.
- chart4.xml ("fast r10 regular vs shrink vs investigate") still has
  `<c:max val="200"/>` on its x-axis and `<c:max val="-600000"/>` on its
  y-axis -- leftover contamination from openpyxl's default Axis object that
  the earlier "fresh NumericAxis" fix evidently didn't actually eliminate.
  These clip the plotted lines against a nonsensical boundary, which is the
  real source of the "chamfered"/cut-off appearance.
- chart4.xml is a plain openpyxl-constructed chart, so it lacks the header
  removal / axis-title placement / fonts / chart-area border that the other
  two charts have (those were hand-authored in Excel, not built by us).

Fix:
1. Per-series: give every added series an explicit color (grouped by run
   "type", reusing the workbook's own colors where a type already existed)
   plus the round cap/join the native series use.
2. chart4 only: remove the bogus axis max values, remove the auto-title
   header (matching autoTitleDeleted=1 like the other charts), move the
   x-axis title into the vacated header space (same manual-layout position
   the other two charts use), style axis titles/labels/gridlines/legend
   fonts and the chart-area border to match, and stop the title/legend from
   overlaying the plot.

Only touches xl/charts/chart2.xml, chart3.xml, chart4.xml. No cell data,
sheet structure, or numeric caches are modified.
"""

import re
import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"

FILLS = {
    "accent1": '<a:schemeClr val="accent1"/>',
    "accent2": '<a:schemeClr val="accent2"/>',
    "accent3": '<a:schemeClr val="accent3"/>',
    "accent4": '<a:schemeClr val="accent4"/>',
    "accent5": '<a:schemeClr val="accent5"/>',
    "red": '<a:srgbClr val="FF0000"/>',
}

BARE_SPPR = '<c:spPr><a:ln w="19050"><a:prstDash val="solid"/></a:ln></c:spPr>'

TITLE_RE = re.compile(r"<c:tx><c:v>(.*?)</c:v></c:tx>")


def styled_sppr(color_key):
    return (f'<c:spPr><a:ln w="19050" cap="rnd"><a:solidFill>{FILLS[color_key]}</a:solidFill>'
            f'<a:prstDash val="solid"/><a:round/></a:ln></c:spPr>')


def classify_chart2(title):
    if title.startswith("r10 findopt") and title.endswith("regular"):
        return "accent1"
    if title.startswith("r10 findopt") and title.endswith("findopt"):
        return "red"
    if title.startswith("r1 findopt") and title.endswith("regular"):
        return "accent2"
    if title.startswith("r1 findopt") and title.endswith("findopt"):
        return "red"
    return None


def classify_chart3(title):
    if title.startswith("r10 findopt") and title.endswith("regular"):
        return "accent1"
    if title.startswith("r10 findopt") and "findopt (BLO)" in title:
        return "red"
    if title.startswith("r10 fullreopt") and title.endswith("regular"):
        return "accent5"
    if title.startswith("r10 fullreopt") and "findopt (BLO)" in title:
        return "red"
    return None


def classify_chart4(title):
    if title.startswith("r10 fast") and title.endswith("regular"):
        return "accent1"
    if title.startswith("r10 investigate3"):
        return "accent3"
    if title.startswith("r10 shrink500"):
        return "accent4"
    return None


def recolor_series(xml, classify):
    parts = []
    pos = 0
    n_fixed = 0
    while True:
        start = xml.find("<c:ser>", pos)
        if start == -1:
            parts.append(xml[pos:])
            break
        end = xml.find("</c:ser>", start)
        assert end != -1
        end += len("</c:ser>")
        block = xml[start:end]
        m = TITLE_RE.search(block)
        title = m.group(1) if m else ""
        color = classify(title)
        if color and BARE_SPPR in block:
            block = block.replace(BARE_SPPR, styled_sppr(color), 1)
            n_fixed += 1
        parts.append(xml[pos:start])
        parts.append(block)
        pos = end
    return "".join(parts), n_fixed


# --- templates copied verbatim (minus text/manualLayout numbers that must
# differ) from the workbook's own chart2.xml / chart3.xml ---

GRIDLINES = ('<c:majorGridlines><c:spPr><a:ln w="9525" cap="flat" cmpd="sng" algn="ctr">'
             '<a:solidFill><a:schemeClr val="tx1"><a:lumMod val="15000"/><a:lumOff val="85000"/>'
             '</a:schemeClr></a:solidFill><a:prstDash val="solid"/><a:round/></a:ln></c:spPr>'
             '</c:majorGridlines>')

AXIS_LINE_AND_LABEL_TXPR = (
    '<c:spPr><a:noFill/><a:ln w="9525" cap="flat" cmpd="sng" algn="ctr">'
    '<a:solidFill><a:schemeClr val="tx1"><a:lumMod val="25000"/><a:lumOff val="75000"/>'
    '</a:schemeClr></a:solidFill><a:prstDash val="solid"/><a:round/></a:ln></c:spPr>'
    '<c:txPr><a:bodyPr rot="-60000000" spcFirstLastPara="1" vertOverflow="ellipsis" vert="horz" '
    'wrap="square" anchor="ctr" anchorCtr="1"/><a:lstStyle/><a:p><a:pPr><a:defRPr sz="900" b="0" '
    'i="0" strike="noStrike" kern="1200" baseline="0"><a:solidFill><a:schemeClr val="tx1">'
    '<a:lumMod val="65000"/><a:lumOff val="35000"/></a:schemeClr></a:solidFill>'
    '<a:latin typeface="+mn-lt"/><a:ea typeface="+mn-ea"/><a:cs typeface="+mn-cs"/></a:defRPr>'
    '</a:pPr><a:endParaRPr lang="en-US"/></a:p></c:txPr>'
)

TITLE_TEXT_PPR = ('<a:pPr><a:defRPr sz="1500" b="0" i="0" strike="noStrike" kern="1200" '
                   'baseline="0"><a:solidFill><a:schemeClr val="tx1"><a:lumMod val="65000"/>'
                   '<a:lumOff val="35000"/></a:schemeClr></a:solidFill><a:latin typeface="+mn-lt"/>'
                   '<a:ea typeface="+mn-ea"/><a:cs typeface="+mn-cs"/></a:defRPr></a:pPr>')

TITLE_NOFILL_SPPR = '<c:spPr><a:noFill/><a:ln><a:noFill/><a:prstDash val="solid"/></a:ln></c:spPr>'

X_AXIS_TITLE = (
    '<c:title><c:tx><c:rich><a:bodyPr rot="0" spcFirstLastPara="1" vertOverflow="ellipsis" '
    'vert="horz" wrap="square" anchor="ctr" anchorCtr="1"/><a:lstStyle/><a:p>' + TITLE_TEXT_PPR +
    '<a:r><a:rPr lang="en-US" sz="1500" baseline="0"/><a:t>Time (s)</a:t></a:r></a:p></c:rich>'
    '</c:tx><c:layout><c:manualLayout><c:xMode val="edge"/><c:yMode val="edge"/>'
    '<c:x val="0.45910591184852878"/><c:y val="2.7411860872958809E-2"/></c:manualLayout></c:layout>'
    '<c:overlay val="0"/>' + TITLE_NOFILL_SPPR + '</c:title>'
)

Y_AXIS_TITLE = (
    '<c:title><c:tx><c:rich><a:bodyPr rot="-5400000" spcFirstLastPara="1" vertOverflow="ellipsis" '
    'vert="horz" wrap="square" anchor="ctr" anchorCtr="1"/><a:lstStyle/><a:p>' + TITLE_TEXT_PPR +
    '<a:r><a:rPr lang="en-US" sz="1500" baseline="0"/><a:t>Log Likelihood</a:t></a:r></a:p>'
    '</c:rich></c:tx><c:overlay val="0"/>' + TITLE_NOFILL_SPPR + '</c:title>'
)

LEGEND_TXPR = (
    '<c:txPr><a:bodyPr rot="0" spcFirstLastPara="1" vertOverflow="ellipsis" vert="horz" '
    'wrap="square" anchor="ctr" anchorCtr="1"/><a:lstStyle/><a:p><a:pPr><a:defRPr sz="800" b="0" '
    'i="0" strike="noStrike" kern="1200" baseline="0"><a:ln><a:noFill/><a:prstDash val="solid"/>'
    '</a:ln><a:solidFill><a:schemeClr val="tx1"><a:lumMod val="65000"/><a:lumOff val="35000"/>'
    '</a:schemeClr></a:solidFill><a:latin typeface="+mn-lt"/><a:ea typeface="+mn-ea"/>'
    '<a:cs typeface="+mn-cs"/></a:defRPr></a:pPr><a:endParaRPr lang="en-US"/></a:p></c:txPr>'
)

CHART_AREA_SPPR = (
    '<c:spPr><a:solidFill><a:schemeClr val="bg1"/></a:solidFill>'
    '<a:ln w="9525" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:schemeClr val="tx1">'
    '<a:lumMod val="15000"/><a:lumOff val="85000"/></a:schemeClr></a:solidFill>'
    '<a:prstDash val="solid"/><a:round/></a:ln></c:spPr>'
)


def fix_chart4(xml):
    # 1. remove the header/title box entirely
    m = re.search(r"<c:title>.*?</c:title><c:autoTitleDeleted val=\"0\"/>", xml, re.DOTALL)
    assert m, "chart4 title block not found"
    xml = xml[:m.start()] + '<c:autoTitleDeleted val="1"/>' + xml[m.end():]

    # 2. x-axis (axId=10): drop bogus max, add gridlines/title/line+label styling
    assert '<c:axId val="10"/><c:scaling><c:orientation val="minMax"/><c:max val="200"/></c:scaling>' in xml
    xml = xml.replace(
        '<c:axId val="10"/><c:scaling><c:orientation val="minMax"/><c:max val="200"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="b"/><c:majorGridlines/>'
        '<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr/></a:pPr>'
        '<a:r><a:rPr lang="en-US"/><a:t>Time (s)</a:t></a:r></a:p></c:rich></c:tx>'
        '<c:overlay val="1"/></c:title>'
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/><c:crossAx val="20"/>',
        '<c:axId val="10"/><c:scaling><c:orientation val="minMax"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="b"/>' + GRIDLINES + X_AXIS_TITLE +
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>' + AXIS_LINE_AND_LABEL_TXPR +
        '<c:crossAx val="20"/>',
        1,
    )

    # 3. y-axis (axId=20): drop bogus max, add gridlines/title/line+label styling
    assert '<c:axId val="20"/><c:scaling><c:orientation val="minMax"/><c:max val="-600000"/></c:scaling>' in xml
    xml = xml.replace(
        '<c:axId val="20"/><c:scaling><c:orientation val="minMax"/><c:max val="-600000"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="l"/><c:majorGridlines/>'
        '<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/><a:p><a:pPr><a:defRPr/></a:pPr>'
        '<a:r><a:rPr lang="en-US"/><a:t>Log Likelihood</a:t></a:r></a:p></c:rich></c:tx>'
        '<c:overlay val="1"/></c:title>'
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/><c:crossAx val="10"/>',
        '<c:axId val="20"/><c:scaling><c:orientation val="minMax"/></c:scaling>'
        '<c:delete val="0"/><c:axPos val="l"/>' + GRIDLINES + Y_AXIS_TITLE +
        '<c:numFmt formatCode="General" sourceLinked="1"/><c:majorTickMark val="none"/>'
        '<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>' + AXIS_LINE_AND_LABEL_TXPR +
        '<c:crossAx val="10"/>',
        1,
    )

    # 4. legend: stop overlaying the plot, match font
    assert '<c:legend><c:legendPos val="r"/><c:overlay val="1"/></c:legend>' in xml
    xml = xml.replace(
        '<c:legend><c:legendPos val="r"/><c:overlay val="1"/></c:legend>',
        '<c:legend><c:legendPos val="r"/><c:overlay val="0"/>' + LEGEND_TXPR + '</c:legend>',
        1,
    )

    # 5. chart-area border/fill, inserted right after </c:chart>
    assert "</c:chart><c:printSettings>" in xml
    xml = xml.replace("</c:chart><c:printSettings>",
                       "</c:chart>" + CHART_AREA_SPPR + "<c:printSettings>", 1)

    return xml


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    c2 = entries["xl/charts/chart2.xml"].decode("utf-8")
    c2, n2 = recolor_series(c2, classify_chart2)
    entries["xl/charts/chart2.xml"] = c2.encode("utf-8")

    c3 = entries["xl/charts/chart3.xml"].decode("utf-8")
    c3, n3 = recolor_series(c3, classify_chart3)
    entries["xl/charts/chart3.xml"] = c3.encode("utf-8")

    c4 = entries["xl/charts/chart4.xml"].decode("utf-8")
    c4, n4 = recolor_series(c4, classify_chart4)
    c4 = fix_chart4(c4)
    entries["xl/charts/chart4.xml"] = c4.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print(f"Recolored: chart2={n2} series, chart3={n3} series, chart4={n4} series.")
    print("chart4: removed header, fixed axis max clipping, restyled axes/legend/border.")


if __name__ == "__main__":
    main()
