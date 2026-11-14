"""Display helpers; statistical results come only from a verified frozen payload.

Formatting (rounding, percentages, thousands separators) is presentation only.
This module does not estimate effects, pool results, or derive confidence limits.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from build_repaired_freeze import read_frozen_payload


class FrozenPublication:
    def __init__(self, frozen_path):
        self.path = Path(frozen_path).resolve()
        self.payload = read_frozen_payload(self.path)
        self.tables = self.payload["tables"]
        self.metadata = self.payload["metadata"]
        self.claims = []

    def rows(self, key, **where):
        return [r for r in self.tables[key]["records"]
                if all(str(r.get(k)) == str(v) for k, v in where.items())]

    def one(self, key, **where):
        rows = self.rows(key, **where)
        if len(rows) != 1:
            raise ValueError(f"Expected one frozen row: {key} {where}; found {len(rows)}")
        self.claims.append({"table": key, "selection": where,
                            "source_file": rows[0]["source_file"],
                            "source_row": rows[0]["source_row"]})
        return rows[0]

    def save_claims(self, path):
        Path(path).write_text(json.dumps(self.claims, ensure_ascii=False, indent=2), encoding="utf-8")

    def source_rows(self, suffix, **where):
        """Select registered source records, de-duplicating identical registration."""
        found = {}
        for key, table in self.tables.items():
            for row in table["records"]:
                if row["source_file"].endswith(suffix) and all(
                    str(row.get(k)) == str(v) for k, v in where.items()
                ):
                    identity = (row["source_file"], row["source_row"])
                    if identity in found and found[identity] != row:
                        raise ValueError(f"Conflicting frozen duplicate: {identity}")
                    found[identity] = row
        return list(found.values())

    def source_one(self, suffix, **where):
        found = self.source_rows(suffix, **where)
        if len(found) != 1:
            raise ValueError(f"Expected one frozen source row: {suffix} {where}; found {len(found)}")
        row = found[0]
        self.claims.append({"source_file": row["source_file"], "source_row": row["source_row"],
                            "selection": where})
        return row


def num(value, digits=3):
    if value is None or str(value) in ("", "NA", "NaN", "nan", "None"):
        return "NA"
    v = float(value)
    if not math.isfinite(v):
        return "NA"
    return f"{v:.{digits}f}"


def count(value):
    v = float(value)
    if v != int(v):
        raise ValueError(f"Count is not integer: {value}")
    return f"{int(v):,}"


def percent(value, digits=1, fraction=False):
    return num(float(value) * (100 if fraction else 1), digits) + "%"


def interval(row, point="RR", low="lo", high="hi", digits=3):
    return f"{num(row[point], digits)} ({num(row[low], digits)}–{num(row[high], digits)})"


def source_csv(rows, path):
    if not rows:
        raise ValueError("Cannot export empty display source")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, keys)
        writer.writeheader()
        writer.writerows(rows)


def md_table(headers, rows):
    def escape(x):
        return str(x).replace("|", "&#124;").replace("\n", " ")
    lines = ["| " + " | ".join(map(escape, headers)) + " |",
             "| " + " | ".join("---" for _ in headers) + " |"]
    lines.extend("| " + " | ".join(map(escape, row)) + " |" for row in rows)
    return "\n".join(lines)
