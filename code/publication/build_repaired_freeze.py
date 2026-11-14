#!/usr/bin/env python3
"""Build a traceable, string-preserving frozen-numbers package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import tempfile
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any


REPAIR_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = Path(__file__).resolve()
BEGIN_MARKER = "<!-- BEGIN_FROZEN_PAYLOAD -->"
END_MARKER = "<!-- END_FROZEN_PAYLOAD -->"
SOURCE_COLUMNS = ("source_file", "source_row")
FORBIDDEN_PATH_PARTS = ("invalid", "quarantine", "smoke", "fixture")


class FreezeError(ValueError):
    """Raised when an input or existing freeze violates the freeze contract."""


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise FreezeError(f"duplicate JSON object key: {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise FreezeError(f"non-standard JSON constant is not allowed: {value}")


def _parse_json_bytes(data: bytes, label: str) -> Any:
    try:
        return json.loads(
            data.decode("utf-8-sig"),
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FreezeError(f"could not parse {label} as UTF-8 JSON: {exc}") from exc


def _read_json_input(path: Path, label: str) -> tuple[Any, str]:
    try:
        data = path.read_bytes()
    except OSError as exc:
        raise FreezeError(f"cannot read {label} {path}: {exc}") from exc
    return _parse_json_bytes(data, label), hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise FreezeError(f"cannot hash input {path}: {exc}") from exc
    return digest.hexdigest()


def _display_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root).as_posix()
    except ValueError:
        return str(path.resolve())


def _resolve_registered_csv(raw_path: Any, root: Path) -> tuple[Path, str]:
    if not isinstance(raw_path, str) or not raw_path:
        raise FreezeError("registry CSV paths must be non-empty strings")

    normalized = raw_path.replace("\\", "/")
    posix_path = PurePosixPath(normalized)
    windows_path = PureWindowsPath(raw_path)
    if posix_path.is_absolute() or windows_path.is_absolute() or windows_path.drive:
        raise FreezeError(f"registry CSV path must be relative to repair root: {raw_path!r}")
    if ".." in normalized.split("/"):
        raise FreezeError(f"registry CSV path traversal is not allowed: {raw_path!r}")
    lowered = normalized.casefold()
    if any(part in lowered for part in FORBIDDEN_PATH_PARTS):
        raise FreezeError(f"registry CSV path uses a disallowed invalid/quarantine/smoke/fixture path: {raw_path!r}")

    candidate = root.joinpath(*posix_path.parts)
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as exc:
        raise FreezeError(f"registered CSV does not exist: {raw_path!r}: {exc}") from exc
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise FreezeError(f"registered CSV resolves outside repair root: {raw_path!r}") from exc
    if not resolved.is_file():
        raise FreezeError(f"registered CSV is not a file: {raw_path!r}")
    return resolved, posix_path.as_posix()


def _read_csv(path: Path, source_file: str) -> tuple[list[str], list[dict[str, Any]], str]:
    before_hash = _sha256_file(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, strict=True)
            try:
                headers = next(reader)
            except StopIteration as exc:
                raise FreezeError(f"registered CSV is empty: {source_file}") from exc

            if not headers or any(not header.strip() for header in headers):
                raise FreezeError(f"CSV headers must be non-empty: {source_file}")
            if len(headers) != len(set(headers)):
                raise FreezeError(f"CSV has duplicate headers: {source_file}")
            overlap = set(headers).intersection(SOURCE_COLUMNS)
            if overlap:
                names = ", ".join(sorted(overlap))
                raise FreezeError(f"CSV uses reserved source metadata column(s) {names}: {source_file}")

            records: list[dict[str, Any]] = []
            for record_number, values in enumerate(reader, start=2):
                if len(values) != len(headers):
                    raise FreezeError(
                        f"CSV record {record_number} in {source_file} has {len(values)} fields; "
                        f"expected {len(headers)}"
                    )
                record = dict(zip(headers, values))
                record["source_file"] = source_file
                # CSV logical record numbers include the header as record 1.
                record["source_row"] = record_number
                records.append(record)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise FreezeError(f"cannot read CSV {source_file}: {exc}") from exc

    if not records:
        raise FreezeError(f"registered CSV has headers but no data records: {source_file}")
    after_hash = _sha256_file(path)
    if before_hash != after_hash:
        raise FreezeError(f"registered CSV changed while being read: {source_file}")
    return headers + list(SOURCE_COLUMNS), records, before_hash


def _markdown_cell(value: Any) -> str:
    # HTML entities keep source markup inert; pipes are encoded so they cannot
    # split Markdown columns. Newlines are made visible within a single cell.
    text = html.escape(str(value), quote=True)
    text = text.replace("|", "&#124;")
    text = text.replace("`", "&#96;")
    text = text.replace("*", "&#42;").replace("_", "&#95;")
    text = text.replace("[", "&#91;").replace("]", "&#93;")
    return text.replace("\r\n", "\n").replace("\r", "\n").replace("\n", "<br>")


def _json_fence(value: Any) -> list[str]:
    rendered = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)
    longest = max((len(match.group(0)) for match in re.finditer(r"`+", rendered)), default=0)
    fence = "`" * max(3, longest + 1)
    return [f"{fence}json", rendered, fence]


def _render_markdown(metadata: Any, tables: dict[str, Any], payload: dict[str, Any]) -> str:
    lines = [
        "# Frozen numbers",
        "",
        "## Metadata",
        "",
        *_json_fence(metadata),
        "",
        "## Tables",
        "",
    ]
    for table_key, table in tables.items():
        lines.extend((f"### {_markdown_cell(table_key)}", "", f"Records: {len(table['records'])}", ""))
        headers = table["headers"]
        lines.append("| " + " | ".join(_markdown_cell(header) for header in headers) + " |")
        lines.append("| " + " | ".join("---" for _ in headers) + " |")
        for record in table["records"]:
            lines.append("| " + " | ".join(_markdown_cell(record[header]) for header in headers) + " |")
        lines.append("")

    lines.extend(("## Machine-readable payload", "", BEGIN_MARKER,
                  *_json_fence(payload), END_MARKER, ""))
    return "\n".join(lines)


def _atomic_write(path: Path, data: bytes) -> Path:
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temp_path = Path(handle.name)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        return temp_path
    except OSError as exc:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise FreezeError(f"cannot stage output {path}: {exc}") from exc


def build_freeze(
    registry_path: Path,
    metadata_path: Path,
    outdir: Path,
    *,
    replace: bool = False,
    repair_root: Path = REPAIR_ROOT,
) -> tuple[Path, Path, dict[str, Any]]:
    """Create a freeze from an explicit registry and root-authored metadata."""
    root = repair_root.resolve(strict=True)
    registry_file = registry_path.expanduser().resolve(strict=True)
    metadata_file = metadata_path.expanduser().resolve(strict=True)
    if not registry_file.is_file() or not metadata_file.is_file():
        raise FreezeError("registry and metadata paths must be files")

    registry, registry_hash = _read_json_input(registry_file, "registry")
    metadata, metadata_hash = _read_json_input(metadata_file, "metadata")
    if not isinstance(registry, dict) or not registry:
        raise FreezeError("registry must be a non-empty JSON object mapping table keys to CSV path lists")
    if not isinstance(metadata, dict):
        raise FreezeError("metadata must be a JSON object")

    tables: dict[str, Any] = {}
    source_inputs: list[dict[str, Any]] = []
    record_counts: dict[str, Any] = {}
    for table_key, raw_paths in registry.items():
        if not isinstance(table_key, str) or not table_key.strip():
            raise FreezeError("registry table keys must be non-empty strings")
        if not isinstance(raw_paths, list) or not raw_paths:
            raise FreezeError(f"registry entry {table_key!r} must be a non-empty list of CSV paths")

        table_headers: list[str] | None = None
        table_records: list[dict[str, Any]] = []
        file_counts: list[dict[str, Any]] = []
        seen_files: set[Path] = set()
        for raw_path in raw_paths:
            csv_path, source_file = _resolve_registered_csv(raw_path, root)
            if csv_path in seen_files:
                raise FreezeError(f"registry repeats CSV {raw_path!r} in table {table_key!r}")
            seen_files.add(csv_path)
            headers, records, source_hash = _read_csv(csv_path, source_file)
            if table_headers is None:
                table_headers = headers
            elif headers != table_headers:
                raise FreezeError(
                    f"CSV headers do not match within table {table_key!r}: {source_file}"
                )
            table_records.extend(records)
            file_counts.append({"source_file": source_file, "record_count": len(records)})
            source_inputs.append(
                {
                    "role": "source_csv",
                    "table_key": table_key,
                    "path": source_file,
                    "sha256": source_hash,
                    "record_count": len(records),
                }
            )

        assert table_headers is not None
        tables[table_key] = {"headers": table_headers, "records": table_records}
        record_counts[table_key] = {
            "files": file_counts,
            "total_records": len(table_records),
        }

    try:
        script_hash = hashlib.sha256(SCRIPT_PATH.read_bytes()).hexdigest()
    except OSError as exc:
        raise FreezeError(f"cannot hash generator script {SCRIPT_PATH}: {exc}") from exc
    provenance = {
        "inputs": [
            {
                "role": "registry",
                "path": _display_path(registry_file, root),
                "sha256": registry_hash,
            },
            {
                "role": "metadata",
                "path": _display_path(metadata_file, root),
                "sha256": metadata_hash,
            },
            *source_inputs,
            {
                "role": "generator_script",
                "path": _display_path(SCRIPT_PATH, root),
                "sha256": script_hash,
            },
        ],
        "record_counts": record_counts,
    }
    payload = {"metadata": metadata, "tables": tables, "provenance": provenance}
    markdown_bytes = _render_markdown(metadata, tables, payload).encode("utf-8")
    checksum_bytes = f"{hashlib.sha256(markdown_bytes).hexdigest()}  frozen_numbers.md\n".encode("ascii")

    admin_dir = outdir.expanduser() / "00_admin"
    markdown_path = admin_dir / "frozen_numbers.md"
    checksum_path = admin_dir / "frozen_numbers.sha256"
    if not replace and (markdown_path.exists() or checksum_path.exists()):
        raise FreezeError(
            f"freeze output already exists under {admin_dir}; pass --replace to write a new version"
        )

    try:
        admin_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise FreezeError(f"cannot create output directory {admin_dir}: {exc}") from exc
    staged_markdown = _atomic_write(markdown_path, markdown_bytes)
    try:
        staged_checksum = _atomic_write(checksum_path, checksum_bytes)
        os.replace(staged_markdown, markdown_path)
        os.replace(staged_checksum, checksum_path)
    except Exception:
        staged_markdown.unlink(missing_ok=True)
        if "staged_checksum" in locals():
            staged_checksum.unlink(missing_ok=True)
        raise
    return markdown_path, checksum_path, payload


def read_frozen_payload(path: str | Path) -> dict[str, Any]:
    """Verify the adjacent checksum and return the unique embedded JSON payload."""
    markdown_path = Path(path)
    if markdown_path.name != "frozen_numbers.md":
        raise FreezeError("payload path must be named frozen_numbers.md")
    checksum_path = markdown_path.with_name("frozen_numbers.sha256")
    try:
        markdown_bytes = markdown_path.read_bytes()
        checksum_bytes = checksum_path.read_bytes()
    except OSError as exc:
        raise FreezeError(f"cannot read frozen file or adjacent checksum: {exc}") from exc

    expected_checksum = f"{hashlib.sha256(markdown_bytes).hexdigest()}  frozen_numbers.md\n".encode("ascii")
    if checksum_bytes != expected_checksum:
        raise FreezeError(f"checksum mismatch or invalid checksum format: {checksum_path}")
    try:
        markdown = markdown_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise FreezeError(f"frozen Markdown is not UTF-8: {exc}") from exc

    lines = markdown.splitlines()
    begin_positions = [index for index, line in enumerate(lines) if line == BEGIN_MARKER]
    end_positions = [index for index, line in enumerate(lines) if line == END_MARKER]
    if len(begin_positions) != 1 or len(end_positions) != 1:
        raise FreezeError("frozen Markdown must contain exactly one begin and one end payload marker")
    begin, end = begin_positions[0], end_positions[0]
    if end - begin < 4:
        raise FreezeError("frozen payload markers are out of order or incomplete")
    opening_fence = re.fullmatch(r"(`{3,})json", lines[begin + 1])
    if opening_fence is None or lines[end - 1] != opening_fence.group(1):
        raise FreezeError("frozen payload markers must enclose one fenced JSON block")

    payload_text = "\n".join(lines[begin + 2 : end - 1])
    try:
        payload = json.loads(
            payload_text,
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except (json.JSONDecodeError, FreezeError) as exc:
        raise FreezeError(f"embedded frozen payload is invalid JSON: {exc}") from exc
    if not isinstance(payload, dict) or set(payload) != {"metadata", "tables", "provenance"}:
        raise FreezeError("embedded frozen payload must contain metadata, tables, and provenance")
    return payload


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True, type=Path, help="JSON mapping table keys to CSV paths")
    parser.add_argument("--metadata", required=True, type=Path, help="root-authored metadata JSON")
    parser.add_argument("--outdir", required=True, type=Path, help="package root; files go under 00_admin")
    parser.add_argument("--replace", action="store_true", help="replace an existing frozen_numbers package")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        markdown_path, checksum_path, payload = build_freeze(
            args.registry,
            args.metadata,
            args.outdir,
            replace=args.replace,
        )
    except (FreezeError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")
    counts = {key: item["total_records"] for key, item in payload["provenance"]["record_counts"].items()}
    print(f"Wrote {markdown_path}")
    print(f"Wrote {checksum_path}")
    print(f"Record counts: {json.dumps(counts, ensure_ascii=False, sort_keys=True)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
