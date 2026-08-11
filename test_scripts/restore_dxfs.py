#!/usr/bin/env python3
"""
One-off repair #2: openpyxl also empties the workbook's <dxfs> (differential
formatting) list down to <dxfs count="0"/> on every save. xl/tables/table3-5.xml
(restored from backup by restore_query_tables.py) reference dxfs by index via
dataDxfId="2"/"1"/"0" -- after openpyxl's save those are dangling references,
which is exactly what Excel's "Repaired Records: Table" notice was about.

Fix: splice the original <dxfs count="3">...</dxfs> block (same 3 entries, in
the same order, so the existing dataDxfId indices stay valid) back in, in the
same schema position openpyxl used for its empty placeholder.

Only touches xl/styles.xml. No cell data, chart data, or table definitions
are modified.
"""

import re
import zipfile

BACKUP = "PhylogenyTreeLengths.backup.xlsx"
TARGET = "PhylogenyTreeLengths.xlsx"


def read_zip(path):
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist()}, z.namelist()


def main():
    backup, _ = read_zip(BACKUP)
    target, order = read_zip(TARGET)

    orig_styles = backup["xl/styles.xml"].decode("utf-8")
    m = re.search(r"<dxfs\b[^>]*>.*?</dxfs>", orig_styles, re.DOTALL)
    assert m, "no <dxfs>...</dxfs> block found in backup styles.xml"
    dxfs_block = m.group(0)
    assert 'count="3"' in dxfs_block

    new_styles = target["xl/styles.xml"].decode("utf-8")
    assert "<dxfs count=\"0\"/>" in new_styles, "expected empty dxfs placeholder not found"
    new_styles = new_styles.replace('<dxfs count="0"/>', dxfs_block)
    target["xl/styles.xml"] = new_styles.encode("utf-8")

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, target[name])

    print("Restored <dxfs> block (3 entries) into xl/styles.xml.")


if __name__ == "__main__":
    main()
