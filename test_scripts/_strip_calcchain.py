#!/usr/bin/env python3
"""
One-off repair #5: Excel's repair log flagged a stale calcChain.xml formula
record. calcChain.xml is purely a calculation-order cache (no real data --
Excel rebuilds it from scratch whenever it's missing), so the standard safe
fix is to remove the part entirely along with its relationship/content-type
entries, rather than hunt for the one stale reference inside it.
"""

import re
import zipfile

TARGET = "PhylogenyTreeLengths.xlsx"


def main():
    with zipfile.ZipFile(TARGET) as z:
        entries = {n: z.read(n) for n in z.namelist()}
        order = z.namelist()

    assert "xl/calcChain.xml" in entries
    del entries["xl/calcChain.xml"]
    order = [n for n in order if n != "xl/calcChain.xml"]

    wb_rels = entries["xl/_rels/workbook.xml.rels"].decode("utf-8")
    new_wb_rels = re.sub(r'<Relationship [^>]*Target="calcChain\.xml"/>', "", wb_rels)
    assert new_wb_rels != wb_rels
    entries["xl/_rels/workbook.xml.rels"] = new_wb_rels.encode("utf-8")

    ct = entries["[Content_Types].xml"].decode("utf-8")
    new_ct = re.sub(r'<Override PartName="/xl/calcChain\.xml"[^>]*/>', "", ct)
    assert new_ct != ct
    entries["[Content_Types].xml"] = new_ct.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, entries[name])

    print("Removed xl/calcChain.xml and its relationship/content-type entries.")


if __name__ == "__main__":
    main()
