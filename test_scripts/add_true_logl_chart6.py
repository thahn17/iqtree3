#!/usr/bin/env python3
"""
Adds one flat "true logL" reference line per dataset to chart6.xml only
(identified via its alt text: "The chart where true logL should be
added"). chart6 has 6 series: r=10/r10-findopt-regular and r=1/r1-findopt-
regular for 3 datasets (07a, 07c, 07d -- note: not 07b).

True logL (= logL + true_minus_current) is a fixed property of the
simulated dataset, not a real time series, so each reference line is just
2 points spanning that dataset's actual time range at a constant y value.
Uses <c:numLit> (literal inline values) instead of <c:numRef>, so no new
worksheet/sheet reference is needed for these tiny 2-point lines.

Only touches xl/charts/chart6.xml.
"""

import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"

# true_logL and (tmin, tmax) per dataset, from yrecord_GTR_FO_fast_findopt.csv
DATASETS = {
    "07a": (-623642.1826232, 3.33, 147.37),
    "07c": (-649071.173, 3.25, 163.14),
    "07d": (-605880.034, 3.30, 153.60),
}


def numlit(values):
    pts = "".join(f'<c:pt idx="{i}"><c:v>{v!r}</c:v></c:pt>' for i, v in enumerate(values))
    return f'<c:formatCode>General</c:formatCode><c:ptCount val="{len(values)}"/>{pts}'


def build_ref_line(idx, label, true_logl, tmin, tmax):
    return (
        f'<c:ser><c:idx val="{idx}"/><c:order val="{idx}"/>'
        f'<c:tx><c:v>true logL {label}</c:v></c:tx>'
        '<c:spPr><a:ln w="12700"><a:solidFill><a:schemeClr val="tx1"><a:lumMod val="50000"/>'
        '<a:lumOff val="50000"/></a:schemeClr></a:solidFill><a:prstDash val="dash"/></a:ln></c:spPr>'
        '<c:marker><c:symbol val="none"/></c:marker>'
        f'<c:xVal><c:numLit>{numlit([tmin, tmax])}</c:numLit></c:xVal>'
        f'<c:yVal><c:numLit>{numlit([true_logl, true_logl])}</c:numLit></c:yVal>'
        '<c:smooth val="0"/></c:ser>'
    )


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    chart6 = entries["xl/charts/chart6.xml"].decode("utf-8")
    assert 'descr="The chart where true logL should be added"' not in chart6  # sanity: this is the chart part, not the drawing
    n_before = chart6.count("<c:ser>")
    assert n_before == 6, f"expected 6 existing series, found {n_before}"

    marker = "<c:dLbls>"
    idx = chart6.find(marker)
    assert idx != -1

    new_series = []
    for i, (label, (true_logl, tmin, tmax)) in enumerate(DATASETS.items()):
        new_series.append(build_ref_line(n_before + i, label, true_logl, tmin, tmax))

    chart6 = chart6[:idx] + "".join(new_series) + chart6[idx:]
    entries["xl/charts/chart6.xml"] = chart6.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print(f"Added {len(DATASETS)} true-logL reference lines to chart6.xml "
          f"({', '.join(DATASETS)}). Total series now {n_before + len(DATASETS)}.")


if __name__ == "__main__":
    main()
