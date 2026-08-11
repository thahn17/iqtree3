#!/usr/bin/env python3
"""
One-off repair: openpyxl silently drops legacy query-table / external-data-
connection parts (xl/connections.xml, xl/queryTables/*, xl/tables/*) on any
save, regardless of what else was changed. This splices those parts back
into PhylogenyTreeLengths.xlsx from PhylogenyTreeLengths.backup.xlsx (the
pre-openpyxl original), via raw zip surgery -- never round-tripped through
openpyxl again, since that would just drop them a second time.

Only touches: xl/connections.xml, xl/queryTables/*, xl/tables/*,
xl/tables/_rels/*, xl/_rels/workbook.xml.rels, [Content_Types].xml,
xl/worksheets/sheet{1,2,3,4,5}.xml(.rels). No cell data, chart data ranges,
or chart XML are modified.
"""

import re
import zipfile

BACKUP = "PhylogenyTreeLengths.backup.xlsx"
TARGET = "PhylogenyTreeLengths.xlsx"

# sheetN.xml -> table index (matches original workbook's table<->sheet mapping)
SHEET_TABLE = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}


def read_zip(path):
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist()}, z.namelist()


def main():
    backup, _ = read_zip(BACKUP)
    target, target_order = read_zip(TARGET)

    # 1. copy over the dropped parts verbatim
    to_add = ["xl/connections.xml"]
    for i in range(1, 6):
        to_add.append(f"xl/queryTables/queryTable{i}.xml")
        to_add.append(f"xl/tables/table{i}.xml")
        to_add.append(f"xl/tables/_rels/table{i}.xml.rels")
    for name in to_add:
        assert name in backup, f"missing in backup: {name}"
        assert name not in target, f"unexpectedly already present: {name}"
        target[name] = backup[name]

    # 2. workbook.xml.rels: add a connections relationship with a fresh rId
    wb_rels = target["xl/_rels/workbook.xml.rels"].decode("utf-8")
    used_ids = [int(m) for m in re.findall(r'Id="rId(\d+)"', wb_rels)]
    new_id = max(used_ids) + 1
    rel = (f'<Relationship Id="rId{new_id}" '
           f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/connections" '
           f'Target="connections.xml"/>')
    wb_rels = wb_rels.replace("</Relationships>", rel + "</Relationships>")
    target["xl/_rels/workbook.xml.rels"] = wb_rels.encode("utf-8")

    # 3. Content_Types: add overrides for connections/queryTables/tables
    ct = target["[Content_Types].xml"].decode("utf-8")
    overrides = ['<Override PartName="/xl/connections.xml" '
                 'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.connections+xml"/>']
    for i in range(1, 6):
        overrides.append(
            f'<Override PartName="/xl/queryTables/queryTable{i}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.queryTable+xml"/>')
        overrides.append(
            f'<Override PartName="/xl/tables/table{i}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.table+xml"/>')
    ct = ct.replace("</Types>", "".join(overrides) + "</Types>")
    target["[Content_Types].xml"] = ct.encode("utf-8")

    # 4. per-sheet: add table relationship + <tableParts> back
    for sheet_num, table_num in SHEET_TABLE.items():
        sheet_path = f"xl/worksheets/sheet{sheet_num}.xml"
        rels_path = f"xl/worksheets/_rels/sheet{sheet_num}.xml.rels"
        rel_type = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/table"
        rel_target = f"../tables/table{table_num}.xml"

        if rels_path in target:
            rels_xml = target[rels_path].decode("utf-8")
            used = [int(m) for m in re.findall(r'Id="rId(\d+)"', rels_xml)]
            rid = max(used) + 1
            new_rel = f'<Relationship Id="rId{rid}" Type="{rel_type}" Target="{rel_target}"/>'
            rels_xml = rels_xml.replace("</Relationships>", new_rel + "</Relationships>")
        else:
            rid = 1
            rels_xml = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                        f'<Relationship Id="rId{rid}" Type="{rel_type}" Target="{rel_target}"/>'
                        '</Relationships>')
            target_order.append(rels_path)
        target[rels_path] = rels_xml.encode("utf-8")

        sheet_xml = target[sheet_path].decode("utf-8")
        assert "<tableParts" not in sheet_xml, f"{sheet_path} already has tableParts"
        table_parts = f'<tableParts count="1"><tablePart r:id="rId{rid}"/></tableParts>'
        assert sheet_xml.rstrip().endswith("</worksheet>")
        sheet_xml = sheet_xml.rstrip()[: -len("</worksheet>")] + table_parts + "</worksheet>"
        target[sheet_path] = sheet_xml.encode("utf-8")

    # 5. write back out, preserving original member order and appending new parts
    order = list(target_order)
    for name in to_add:
        if name not in order:
            order.append(name)

    with zipfile.ZipFile(TARGET, "w", zipfile.ZIP_DEFLATED) as z:
        for name in order:
            z.writestr(name, target[name])

    print(f"Restored {len(to_add)} parts, patched {len(SHEET_TABLE)} worksheets, "
          f"workbook.xml.rels (+rId{new_id}), Content_Types (+{len(overrides)} overrides).")


if __name__ == "__main__":
    main()
