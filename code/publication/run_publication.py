#!/usr/bin/env python3
"""Render publication figures and tables from a checksum-verified freeze."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from compose_revision_tables import compose
from render_revision_figures import render


FIGURE_STEMS = ("fig1_design", "fig2_forest", "fig3_contrast", "fig4_curves")
FIGURE_FORMATS = ("png", "pdf", "svg")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(frozen: Path, outdir: Path) -> dict:
    frozen = frozen.expanduser().resolve(strict=True)
    outdir = outdir.expanduser().resolve()
    if frozen.name != "frozen_numbers.md":
        raise ValueError("--frozen must point to a file named frozen_numbers.md")
    outdir.mkdir(parents=True, exist_ok=True)

    figure_manifest = render(frozen, outdir)
    # Apply the compact submission layout to the same verified numerical rows.
    from render_readable_figures import render as render_readable
    readable = render_readable(frozen, outdir / "outputs")
    figure_manifest["submission_layout_figures"] = [p.relative_to(outdir).as_posix() for p in readable]
    table_manifest = compose(frozen, outdir)
    outputs = outdir / "outputs"
    main_figures = [outputs / f"{stem}.{ext}" for stem in FIGURE_STEMS for ext in FIGURE_FORMATS]
    supplementary_figure = [outputs / "supplementary_figures" / f"figS1_mrt_exploratory.{ext}"
                            for ext in FIGURE_FORMATS]
    main_tables = [outputs / "tables" / f"table{number}.csv" for number in range(1, 5)]
    missing = [path for path in (*main_figures, *supplementary_figure, *main_tables) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Expected publication artifacts were not generated: "
                                + ", ".join(str(path) for path in missing))

    tables_text = (outputs / "tables.md").read_text(encoding="utf-8")
    table_numbers = {int(value) for value in re.findall(r"^## Table ([1-4])\s*$", tables_text, re.MULTILINE)}
    if table_numbers != {1, 2, 3, 4}:
        raise ValueError(f"Expected four main tables; found section numbers {sorted(table_numbers)}")
    full_calendar_note = ("Admission totals are full-health-calendar records, not unique people; "
                          "they include all calendar unit-days rather than only exposure-eligible days.")
    if full_calendar_note not in tables_text:
        raise ValueError("Table 1 is missing the full-health-calendar admission-count note")

    supplementary_text = (outputs / "supplementary_tables.md").read_text(encoding="utf-8")
    supplement_numbers = {int(value) for value in re.findall(r"^## Table S(\d+)\b", supplementary_text, re.MULTILINE)}
    expected_supplements = set(range(1, 17))
    if supplement_numbers != expected_supplements:
        raise ValueError("Expected supplementary tables S1–S16; found "
                         + ", ".join(f"S{n}" for n in sorted(supplement_numbers)))

    result = {
        "freeze": frozen.name,
        "freeze_sha256": sha256(frozen),
        "main_figure_count": len(FIGURE_STEMS),
        "supplementary_figure_count": 1,
        "main_table_count": len(main_tables),
        "supplementary_tables": [f"S{number}" for number in sorted(supplement_numbers)],
        "figures": [path.relative_to(outdir).as_posix() for path in main_figures],
        "supplementary_figure": [path.relative_to(outdir).as_posix() for path in supplementary_figure],
        "tables": [path.relative_to(outdir).as_posix() for path in main_tables],
        "numeric_policy": "All displayed estimates are selected from the checksum-verified freeze; no statistical estimation is run.",
        "figure_manifest": figure_manifest,
        "table_manifest": table_manifest,
    }
    (outputs / "publication_runner_manifest.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True, type=Path,
                        help="checksum-verified frozen_numbers.md")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="publication output root; generated files go under outputs/")
    args = parser.parse_args()
    run(args.frozen, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
