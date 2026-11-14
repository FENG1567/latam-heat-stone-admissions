#!/usr/bin/env python3
"""Build privacy-safe Mexico admission-origin grouped counts.

Only selected SAEH EGRESO columns are streamed. No identifiers or record-level
rows are written. Per-source-year SQLite aggregates make the extraction
bounded in memory and resumable against SHA-256 checked grouped caches.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import re
import sqlite3
import zipfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


REPAIR = Path(__file__).resolve().parents[1]
DATA_DIR = REPAIR / "data"
ADMIN_DIR = REPAIR / "admin"
EXISTING_EXTRACTOR = REPAIR / "scripts" / "extract_mexico_clinical.py"
SOURCE_ROOT = Path(r"D:\数据库集合\新增数据\Mexico_DGIS")
CATALOG_ZIP = SOURCE_ROOT / "clues" / "2019" / "ssa_Catalogos_Egresos_2019.zip"
BASELINE_STATE_DAY = DATA_DIR / "mx_complete_source_state_day_2015_2022.csv"
ALL_SOURCE_YEARS = (2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023)
EVENT_YEARS = tuple(range(2015, 2024))
MODEL_YEARS = tuple(range(2015, 2023))
STATE_CODES = {f"{i:02d}" for i in range(1, 33)}
TARGET_DX = {"K80", "K81", "N20", "N21", "N22", "N23"}
CHUNK_SIZE = 200_000
CACHE_VERSION = 1
SEPARATOR_CANDIDATES = ("|", ";", ",", "\t")
MISSING = "__MISSING__"
YEAR_COLUMNS = [
    "state", "date", "dx3", "origin_raw", "origin_group", "los_group", "type_service", "n"
]
COMBINED_COLUMNS = YEAR_COLUMNS
BASELINE_ENDPOINT_COLUMNS = {
    "K80": "k80",
    "K81": "k81",
    "N20": "n20_only",
    "N21": "n21_only",
    "N22": "n22_other_disease",
    "N23": "n23_unspecified_renal_colic",
}


def load_existing_extractor():
    spec = importlib.util.spec_from_file_location("mexico_clinical_extractor", EXISTING_EXTRACTOR)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load existing source parser: {EXISTING_EXTRACTOR}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXISTING = load_existing_extractor()


def canonical_header(value: str) -> str:
    value = value.replace("\ufeff", "").strip().upper()
    return re.sub(r"[^A-Z0-9]", "", value)


def clean_category(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.replace(r"\.0+$", "", regex=True)


def category_code(value: str) -> str | None:
    if not value or value == MISSING:
        return None
    value = re.sub(r"\.0+$", "", value.strip())
    if value.isdigit():
        return str(int(value))
    return None


def origin_group(value: str) -> str:
    code = category_code(value)
    return {
        "1": "outpatient_consultation_origin",
        "2": "ED_origin",
        "3": "referred_origin",
        "4": "other_origin",
        "6": "pathological_nursery_origin",
        "9": "unspecified_origin",
    }.get(code, "unknown_origin_missing" if value == MISSING else "unknown_origin_code")


def inspect_source(path: Path, year: int) -> dict[str, Any]:
    encoding = EXISTING.source_encoding(year)
    with path.open("r", encoding=encoding, errors="replace", newline="") as handle:
        header_line = handle.readline()
    delimiter = next((sep for sep in SEPARATOR_CANDIDATES if header_line.count(sep) > 5), "|")
    raw_headers = next(csv.reader([header_line], delimiter=delimiter))
    normalized = {canonical_header(name): name for name in raw_headers}
    fields = ("INGRE", "EGRESO", "ENTIDAD", "AFECPRIN", "PROCED", "TIPSERV", "DIAS_ESTA")
    selected = {field: normalized.get(canonical_header(field)) for field in fields}
    missing = [field for field, raw_name in selected.items() if raw_name is None]
    if missing:
        raise ValueError(f"{year}: missing required origin fields {missing}; headers={raw_headers!r}")
    return {
        "encoding": encoding,
        "delimiter": delimiter,
        "header_fields": raw_headers,
        "selected_fields": selected,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def discover_catalog_archives() -> list[str]:
    found = []
    for base in (SOURCE_ROOT / "clues", SOURCE_ROOT / "evidence"):
        if not base.exists():
            continue
        found.extend(
            str(path)
            for path in base.rglob("*.zip")
            if "catalog" in path.name.lower() and "egresos" in path.name.lower()
        )
    return sorted(set(found))


def read_catalogs() -> dict[str, Any]:
    if not CATALOG_ZIP.is_file():
        raise FileNotFoundError(CATALOG_ZIP)
    archive_hash = sha256_file(CATALOG_ZIP)
    expected_hash = "f631223e61b5dec47fd0651249877fa2da3ff25aacb6f40576571252c85dd3b7"
    if archive_hash != expected_hash:
        raise ValueError(f"2019 catalog archive hash differs from the verified copy: {archive_hash}")
    wanted = {
        "CATPROCEDENCIA": "Catalogos Egresos 2019/CATPROCEDENCIA.csv",
        "CATTIPOSERV": "Catalogos Egresos 2019/CATTIPOSERV.csv",
    }
    contents: dict[str, Any] = {}
    with zipfile.ZipFile(CATALOG_ZIP) as archive:
        corrupt_member = archive.testzip()
        if corrupt_member is not None:
            raise ValueError(f"Catalog ZIP CRC failure: {corrupt_member}")
        for label, member_name in wanted.items():
            member = archive.getinfo(member_name)
            text = archive.read(member_name).decode("utf-8-sig")
            rows = list(csv.DictReader(text.splitlines()))
            contents[label] = {
                "member": member_name,
                "member_bytes": member.file_size,
                "member_crc32": f"{member.CRC:08x}",
                "rows": rows,
            }
    expected_origin = {
        "1": "CONSULTA EXTERNA",
        "2": "URGENCIAS",
        "3": "REFERIDO",
        "4": "OTRO",
        "6": "CUNERO PATOLÓGICO",
        "9": "NO ESPECIFICADO",
    }
    observed_origin = {
        str(row["IDCatProcedencia"]): str(row["Descrip"]).strip()
        for row in contents["CATPROCEDENCIA"]["rows"]
    }
    if observed_origin != expected_origin:
        raise ValueError(f"Unexpected CATPROCEDENCIA contents: {observed_origin}")
    observed_service = {
        str(row["IDTipoServ"]): str(row["Descrip"]).strip()
        for row in contents["CATTIPOSERV"]["rows"]
    }
    if observed_service != {"1": "NORMAL", "2": "CORTA ESTANCIA"}:
        raise ValueError(f"Unexpected CATTIPOSERV contents: {observed_service}")
    return {
        "archive_path": str(CATALOG_ZIP),
        "archive_bytes": CATALOG_ZIP.stat().st_size,
        "archive_sha256": archive_hash,
        "zip_crc_check": "all members passed",
        "annual_catalog_archives_found_locally": discover_catalog_archives(),
        "catalogs": contents,
        "interpretation": {
            "PROCED": "admission origin/procedencia; not a surgical procedure code",
            "TIPSERV": "type of service at admission; not an urgency/elective classification",
            "mapping_scope": "Only a 2019 annual catalog is available locally. The 2019 field definition and code labels are extrapolated to 2015-2023 raw source fields, with year-specific raw code frequencies audited.",
        },
    }


def source_cache_paths(year: int) -> tuple[Path, Path]:
    return (
        DATA_DIR / f"mx_origin_source_year_{year}.csv.gz",
        ADMIN_DIR / f"mx_origin_source_year_{year}.json",
    )


def source_cache_current(year: int, source_hash: str) -> dict[str, Any] | None:
    data_path, audit_path = source_cache_paths(year)
    if not (data_path.is_file() and audit_path.is_file()):
        return None
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if audit.get("cache_version") != CACHE_VERSION or audit.get("source", {}).get("sha256") != source_hash:
        return None
    return audit


def clean_state(series: pd.Series) -> pd.Series:
    raw = clean_category(series)
    raw = raw.str.replace(r"\.0+$", "", regex=True)
    numeric = raw.str.fullmatch(r"\d+")
    result = raw.copy()
    result.loc[numeric] = result.loc[numeric].str.zfill(2)
    return result


def build_source_year(year: int, path: Path, source_hash: str, info: dict[str, Any]) -> dict[str, Any]:
    data_path, audit_path = source_cache_paths(year)
    temp_db = DATA_DIR / f"mx_origin_work_{year}.sqlite"
    temp_csv = DATA_DIR / f"mx_origin_source_year_{year}.csv.gz.tmp"
    if temp_db.exists():
        temp_db.unlink()
    if temp_csv.exists():
        temp_csv.unlink()
    conn = sqlite3.connect(temp_db)
    conn.execute(
        """CREATE TABLE agg (
               state TEXT NOT NULL,
               date TEXT NOT NULL,
               dx3 TEXT NOT NULL,
               origin_raw TEXT NOT NULL,
               origin_group TEXT NOT NULL,
               los_group TEXT NOT NULL,
               type_service TEXT NOT NULL,
               n INTEGER NOT NULL,
               PRIMARY KEY(state,date,dx3,origin_raw,origin_group,los_group,type_service)
           )"""
    )
    upsert = (
        "INSERT INTO agg(state,date,dx3,origin_raw,origin_group,los_group,type_service,n) VALUES(?,?,?,?,?,?,?,?) "
        "ON CONFLICT(state,date,dx3,origin_raw,origin_group,los_group,type_service) DO UPDATE SET n=n+excluded.n"
    )
    selected_raw = list(info["selected_fields"].values())
    raw_to_canonical = {raw: canonical for canonical, raw in info["selected_fields"].items()}
    raw_origin_counts: Counter[str] = Counter()
    unknown_origin_counts: Counter[str] = Counter()
    raw_service_counts: Counter[str] = Counter()
    raw_state_counts: Counter[str] = Counter()
    diagnosis_counts: Counter[str] = Counter()
    event_year_counts: Counter[str] = Counter()
    all_rows = 0
    target_rows = 0
    target_rows_valid_ingre = 0
    target_rows_invalid_ingre = 0
    target_rows_in_event_window = 0
    target_rows_outside_event_window = 0
    target_rows_matched_state = 0
    target_rows_unmatched_state = 0
    target_rows_invalid_or_missing_egreso = 0
    target_rows_ingre_after_egreso = 0
    grouped_output_rows = 0

    try:
        chunks = pd.read_csv(
            path,
            sep=info["delimiter"],
            usecols=selected_raw,
            dtype=str,
            chunksize=CHUNK_SIZE,
            encoding=info["encoding"],
            on_bad_lines="skip",
            low_memory=False,
            keep_default_na=False,
        )
        for chunk in chunks:
            chunk = chunk.rename(columns=raw_to_canonical)
            all_rows += len(chunk)
            diagnosis = chunk["AFECPRIN"].fillna("").astype(str).str.upper().str.strip()
            dx3 = diagnosis.str.replace(r"[^A-Z0-9]", "", regex=True).str[:3]
            target = dx3.isin(TARGET_DX)
            if not target.any():
                continue
            target_rows += int(target.sum())

            origin = clean_category(chunk["PROCED"])
            origin = origin.mask(origin.eq(""), MISSING)
            service = clean_category(chunk["TIPSERV"])
            service = service.mask(service.eq(""), MISSING)
            raw_origin_counts.update(origin.loc[target].value_counts().to_dict())
            raw_service_counts.update(service.loc[target].value_counts().to_dict())
            groups = origin.map(origin_group)
            unknown_mask = groups.isin(["unknown_origin_missing", "unknown_origin_code"])
            unknown_origin_counts.update(origin.loc[target & unknown_mask].value_counts().to_dict())
            diagnosis_counts.update(dx3.loc[target].value_counts().to_dict())

            ingre = pd.to_datetime(chunk["INGRE"], errors="coerce")
            egreso = pd.to_datetime(chunk["EGRESO"], errors="coerce")
            ingre_ok = ingre.notna()
            egreso_ok = egreso.notna()
            target_rows_valid_ingre += int((target & ingre_ok).sum())
            target_rows_invalid_ingre += int((target & ~ingre_ok).sum())
            target_rows_invalid_or_missing_egreso += int((target & ~egreso_ok).sum())
            target_rows_ingre_after_egreso += int((target & ingre_ok & egreso_ok & ingre.gt(egreso)).sum())
            event_window = ingre.dt.year.between(EVENT_YEARS[0], EVENT_YEARS[-1]).fillna(False)
            in_window_target = target & ingre_ok & event_window
            target_rows_in_event_window += int(in_window_target.sum())
            target_rows_outside_event_window += int((target & ingre_ok & ~event_window).sum())
            if in_window_target.any():
                event_by_dx = pd.DataFrame({"year": ingre.dt.year.loc[in_window_target], "dx3": dx3.loc[in_window_target]})
                for row in event_by_dx.groupby(["year", "dx3"], sort=False).size().items():
                    (event_year, endpoint), count = row
                    event_year_counts[f"{int(event_year)}|{endpoint}"] += int(count)

            state = clean_state(chunk["ENTIDAD"])
            raw_state_counts.update(
                chunk.loc[target, "ENTIDAD"].fillna("").astype(str).str.strip().replace("", MISSING).value_counts().to_dict()
            )
            state_ok = state.isin(STATE_CODES)
            target_rows_matched_state += int((in_window_target & state_ok).sum())
            target_rows_unmatched_state += int((in_window_target & ~state_ok).sum())
            keep = in_window_target & state_ok
            if not keep.any():
                continue

            raw_los = clean_category(chunk["DIAS_ESTA"])
            los_value = pd.to_numeric(raw_los, errors="coerce")
            los_valid = los_value.notna() & los_value.ge(0) & los_value.eq(los_value.round())
            los_group = pd.Series("unknown", index=chunk.index, dtype="object")
            los_group.loc[los_valid & los_value.eq(0)] = "0"
            los_group.loc[los_valid & los_value.between(1, 30)] = "1to30"
            los_group.loc[los_valid & los_value.gt(30)] = "over30"

            selected = pd.DataFrame(
                {
                    "state": state.loc[keep].to_numpy(),
                    "date": ingre.loc[keep].dt.strftime("%Y-%m-%d").to_numpy(),
                    "dx3": dx3.loc[keep].to_numpy(),
                    "origin_raw": origin.loc[keep].to_numpy(),
                    "origin_group": groups.loc[keep].to_numpy(),
                    "los_group": los_group.loc[keep].to_numpy(),
                    "type_service": service.loc[keep].to_numpy(),
                }
            )
            grouped = selected.groupby(YEAR_COLUMNS[:-1], sort=False, observed=True).size().reset_index(name="n")
            conn.executemany(
                upsert,
                [
                    (r.state, r.date, r.dx3, r.origin_raw, r.origin_group, r.los_group, r.type_service, int(r.n))
                    for r in grouped.itertuples(index=False)
                ],
            )
            conn.commit()
            grouped_output_rows += len(selected)

        source_output_rows = int(conn.execute("SELECT COALESCE(SUM(n),0) FROM agg").fetchone()[0])
        with gzip.open(temp_csv, "wt", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(["source_year", *YEAR_COLUMNS])
            cursor = conn.execute(
                "SELECT state,date,dx3,origin_raw,origin_group,los_group,type_service,n FROM agg "
                "ORDER BY state,date,dx3,origin_raw,los_group,type_service"
            )
            while True:
                rows = cursor.fetchmany(100_000)
                if not rows:
                    break
                writer.writerows([[year, *row] for row in rows])
    finally:
        conn.close()
        if temp_db.exists():
            temp_db.unlink()

    os.replace(temp_csv, data_path)
    audit = {
        "cache_version": CACHE_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_year": year,
        "source": {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": source_hash,
            "encoding": info["encoding"],
            "delimiter": info["delimiter"],
            "field_count": len(info["header_fields"]),
            "field_names": info["header_fields"],
            "selected_fields": info["selected_fields"],
        },
        "ingestion": {
            "rows_parsed": all_rows,
            "chunk_rows": CHUNK_SIZE,
            "max_simultaneously_open_source_years": 3,
            "bad_line_policy": "skip (pandas C parser does not expose skipped-line count)",
            "target_rows_by_dx3_all_dates": dict(sorted(diagnosis_counts.items())),
            "target_rows_valid_ingre": target_rows_valid_ingre,
            "target_rows_invalid_or_missing_ingre": target_rows_invalid_ingre,
            "target_rows_valid_egreso": target_rows - target_rows_invalid_or_missing_egreso,
            "target_rows_invalid_or_missing_egreso": target_rows_invalid_or_missing_egreso,
            "target_rows_ingre_after_egreso": target_rows_ingre_after_egreso,
            "target_rows_in_event_window_2015_2023": target_rows_in_event_window,
            "target_rows_outside_event_window": target_rows_outside_event_window,
            "target_rows_in_window_state_01_32": target_rows_matched_state,
            "target_rows_in_window_unmatched_state": target_rows_unmatched_state,
            "grouped_output_rows_before_aggregation": grouped_output_rows,
            "grouped_output_count_after_aggregation": source_output_rows,
            "valid_ingre_event_year_dx3_counts": dict(sorted(event_year_counts.items())),
        },
        "field_quality": {
            "PROCED_raw_value_counts_all_target_diagnoses": dict(sorted(raw_origin_counts.items())),
            "PROCED_unrecognized_or_missing_raw_counts_all_target_diagnoses": dict(sorted(unknown_origin_counts.items())),
            "TIPSERV_raw_value_counts_all_target_diagnoses": dict(sorted(raw_service_counts.items())),
            "ENTIDAD_raw_value_counts_all_target_diagnoses": dict(sorted(raw_state_counts.items())),
        },
        "aggregation": {
            "grain": "source_year x ENTIDAD-residence-state x INGRE-admission-date x diagnosis-three-character-category x raw-PROCED x PROCED-origin-group x recorded-DIAS_ESTA-group x raw-TIPSERV",
            "target_diagnoses": sorted(TARGET_DX),
            "origin_group_rule": "PROCED code 2=ED origin; 1=outpatient-consultation origin; 3,4,6,9 each retained separately; all unrecognized/missing values retained separately by raw value and marked unknown.",
            "los_group_rule": "DIAS_ESTA integer 0 -> 0; integer 1-30 -> 1to30; integer >30 -> over30; missing, negative, or noninteger -> unknown.",
            "TIPSERV_rule": "Raw category only; not interpreted as urgency/elective.",
            "patient_level_rows_emitted": False,
        },
        "output": str(data_path),
    }
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return audit


def merge_source_caches(years: list[int], catalog: dict[str, Any]) -> dict[str, Any]:
    master_path = DATA_DIR / "mx_origin_merge_work.sqlite"
    if master_path.exists():
        master_path.unlink()
    conn = sqlite3.connect(master_path)
    conn.execute(
        """CREATE TABLE agg (
               state TEXT NOT NULL, date TEXT NOT NULL, dx3 TEXT NOT NULL,
               origin_raw TEXT NOT NULL, origin_group TEXT NOT NULL,
               los_group TEXT NOT NULL, type_service TEXT NOT NULL, n INTEGER NOT NULL,
               PRIMARY KEY(state,date,dx3,origin_raw,origin_group,los_group,type_service)
           )"""
    )
    upsert = (
        "INSERT INTO agg(state,date,dx3,origin_raw,origin_group,los_group,type_service,n) VALUES(?,?,?,?,?,?,?,?) "
        "ON CONFLICT(state,date,dx3,origin_raw,origin_group,los_group,type_service) DO UPDATE SET n=n+excluded.n"
    )
    try:
        for year in years:
            year_path, _ = source_cache_paths(year)
            for chunk in pd.read_csv(year_path, usecols=YEAR_COLUMNS, dtype=str, chunksize=100_000):
                if chunk.empty:
                    continue
                conn.executemany(
                    upsert,
                    [
                        (r.state, r.date, r.dx3, r.origin_raw, r.origin_group, r.los_group, r.type_service, int(r.n))
                        for r in chunk.itertuples(index=False)
                    ],
                )
                conn.commit()

        combined_path = DATA_DIR / "mx_origin_grouped_2015_2023.csv.gz"
        temp_combined = DATA_DIR / "mx_origin_grouped_2015_2023.csv.gz.tmp"
        with gzip.open(temp_combined, "wt", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(COMBINED_COLUMNS)
            cursor = conn.execute(
                "SELECT state,date,dx3,origin_raw,origin_group,los_group,type_service,n FROM agg "
                "ORDER BY state,date,dx3,origin_raw,los_group,type_service"
            )
            while True:
                rows = cursor.fetchmany(100_000)
                if not rows:
                    break
                writer.writerows(rows)
        os.replace(temp_combined, combined_path)

        state_day_path = DATA_DIR / "mx_origin_all_origin_state_day_2015_2022.csv.gz"
        temp_state_day = DATA_DIR / "mx_origin_all_origin_state_day_2015_2022.csv.gz.tmp"
        state_days = conn.execute(
            """SELECT state,date,
                      SUM(CASE WHEN dx3='K80' THEN n ELSE 0 END) AS k80,
                      SUM(CASE WHEN dx3='K81' THEN n ELSE 0 END) AS k81,
                      SUM(CASE WHEN dx3='N20' THEN n ELSE 0 END) AS n20_only,
                      SUM(CASE WHEN dx3='N21' THEN n ELSE 0 END) AS n21_only,
                      SUM(CASE WHEN dx3='N22' THEN n ELSE 0 END) AS n22_other_disease,
                      SUM(CASE WHEN dx3='N23' THEN n ELSE 0 END) AS n23_unspecified_renal_colic
                 FROM agg WHERE date >= '2015-01-01' AND date < '2023-01-01'
                 GROUP BY state,date ORDER BY state,date"""
        )
        with gzip.open(temp_state_day, "wt", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(["state", "date", *BASELINE_ENDPOINT_COLUMNS.values()])
            while True:
                rows = state_days.fetchmany(100_000)
                if not rows:
                    break
                writer.writerows(rows)
        os.replace(temp_state_day, state_day_path)

        reconciliation = reconcile_state_day(conn)
        return {
            "source_years_included": years,
            "source_years_excluded": sorted(set(ALL_SOURCE_YEARS) - set(years)),
            "catalog": catalog,
            "outputs": {
                "grouped_all_origins_2015_2023": str(combined_path),
                "all_origin_state_day_2015_2022": str(state_day_path),
            },
            "grouped_rows_2015_2023": int(conn.execute("SELECT COUNT(*) FROM agg").fetchone()[0]),
            "events_all_sources_by_event_year_dx3_2015_2023": event_year_counts(conn),
            "all_origin_reconciliation_vs_existing_complete_source": reconciliation,
        }
    finally:
        conn.close()
        if master_path.exists():
            master_path.unlink()


def event_year_counts(conn: sqlite3.Connection) -> dict[str, int]:
    rows = conn.execute(
        "SELECT substr(date,1,4),dx3,SUM(n) FROM agg GROUP BY substr(date,1,4),dx3 ORDER BY 1,2"
    )
    return {f"{year}|{dx3}": int(n) for year, dx3, n in rows}


def reconcile_state_day(conn: sqlite3.Connection) -> dict[str, Any]:
    if not BASELINE_STATE_DAY.is_file():
        return {"status": "baseline_missing", "path": str(BASELINE_STATE_DAY)}
    actual = pd.read_sql_query(
        """SELECT state,date,
                  SUM(CASE WHEN dx3='K80' THEN n ELSE 0 END) AS k80,
                  SUM(CASE WHEN dx3='K81' THEN n ELSE 0 END) AS k81,
                  SUM(CASE WHEN dx3='N20' THEN n ELSE 0 END) AS n20_only,
                  SUM(CASE WHEN dx3='N21' THEN n ELSE 0 END) AS n21_only,
                  SUM(CASE WHEN dx3='N22' THEN n ELSE 0 END) AS n22_other_disease,
                  SUM(CASE WHEN dx3='N23' THEN n ELSE 0 END) AS n23_unspecified_renal_colic
             FROM agg WHERE date >= '2015-01-01' AND date < '2023-01-01'
             GROUP BY state,date""",
        conn,
    )
    baseline_cols = ["state", "date", *BASELINE_ENDPOINT_COLUMNS.values()]
    baseline = pd.read_csv(
        BASELINE_STATE_DAY,
        usecols=baseline_cols,
        dtype={"state": str, "date": str},
        low_memory=False,
    )
    joined = baseline.merge(actual, on=["state", "date"], how="outer", suffixes=("_baseline", "_origin"), indicator=True)
    checks = {}
    example_differences = []
    for dx3, column in BASELINE_ENDPOINT_COLUMNS.items():
        left = joined[f"{column}_baseline"].fillna(0).astype("int64")
        right = joined[f"{column}_origin"].fillna(0).astype("int64")
        different = left.ne(right)
        checks[dx3] = {
            "baseline_total": int(left.sum()),
            "origin_total": int(right.sum()),
            "difference": int((right - left).sum()),
            "mismatched_state_days": int(different.sum()),
        }
        if different.any():
            for index in joined.index[different][:10]:
                example_differences.append(
                    {
                        "state": str(joined.at[index, "state"]),
                        "date": str(joined.at[index, "date"]),
                        "dx3": dx3,
                        "baseline": int(left.iloc[index]),
                        "origin": int(right.iloc[index]),
                    }
                )
    left_only = joined["_merge"].eq("left_only")
    right_only = joined["_merge"].eq("right_only")
    baseline_only_nonzero = left_only & joined[[f"{column}_baseline" for column in BASELINE_ENDPOINT_COLUMNS.values()]].fillna(0).sum(axis=1).gt(0)
    origin_only_nonzero = right_only & joined[[f"{column}_origin" for column in BASELINE_ENDPOINT_COLUMNS.values()]].fillna(0).sum(axis=1).gt(0)
    missing_nonzero_state_days = int(baseline_only_nonzero.sum() + origin_only_nonzero.sum())
    return {
        "status": "exact_match" if all(v["mismatched_state_days"] == 0 for v in checks.values()) and missing_nonzero_state_days == 0 else "differences_found",
        "baseline_path": str(BASELINE_STATE_DAY),
        "event_date_window": "2015-01-01 through 2022-12-31",
        "baseline_and_origin_state_day_union_rows": int(len(joined)),
        "baseline_only_zero_event_state_days": int(left_only.sum() - baseline_only_nonzero.sum()),
        "origin_only_nonzero_state_days": int(origin_only_nonzero.sum()),
        "baseline_only_nonzero_state_days": int(baseline_only_nonzero.sum()),
        "counts_by_dx3": checks,
        "first_differences_do_not_drive_or_tune_source_selection": example_differences,
    }


def parse_years(raw: str | None) -> list[int]:
    if not raw:
        return list(ALL_SOURCE_YEARS)
    try:
        years = sorted({int(part.strip()) for part in raw.split(",") if part.strip()})
    except ValueError as error:
        raise argparse.ArgumentTypeError("--years must be a comma-separated year list") from error
    invalid = set(years) - set(ALL_SOURCE_YEARS)
    if not years or invalid:
        raise argparse.ArgumentTypeError(f"unsupported source years: {sorted(invalid)}")
    return years


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--years", help="optional comma-separated source-year smoke subset")
    parser.add_argument("--workers", type=int, default=3, help="concurrent source-year readers (maximum 3)")
    parser.add_argument("--no-combine", action="store_true", help="write/update only per-source-year grouped caches")
    args = parser.parse_args()
    if not 1 <= args.workers <= 3:
        parser.error("--workers must be between 1 and 3 to bound concurrently open source files")
    years = parse_years(args.years)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ADMIN_DIR.mkdir(parents=True, exist_ok=True)

    catalog = read_catalogs()
    specification_path = ADMIN_DIR / "mx_origin_specification.json"
    if not specification_path.is_file():
        raise FileNotFoundError(f"Lock the scientific interpretation before extraction: {specification_path}")
    specification = json.loads(specification_path.read_text(encoding="utf-8"))
    if specification.get("specification_status") != "locked_before_origin_outcome_extraction_and_effect_estimation":
        raise ValueError("Mexico origin specification is not locked")
    locked_catalog_hash = specification.get("interpretation", {}).get("codebook_source_sha256")
    if locked_catalog_hash != catalog["archive_sha256"]:
        raise ValueError("Locked specification catalog hash does not match the locally verified 2019 catalog archive")

    inputs: dict[int, tuple[Path, str, dict[str, Any]]] = {}
    for year in years:
        path = EXISTING.source_year_path(year)
        info = inspect_source(path, year)
        digest = sha256_file(path)
        inputs[year] = (path, digest, info)

    results: dict[int, dict[str, Any]] = {}
    to_build = []
    for year in years:
        path, digest, info = inputs[year]
        cached = source_cache_current(year, digest)
        if cached is not None:
            results[year] = cached
            print(f"{year}: grouped cache reused ({cached['ingestion']['grouped_output_count_after_aggregation']:,} events)", flush=True)
        else:
            to_build.append((year, path, digest, info))

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(build_source_year, year, path, digest, info): year
            for year, path, digest, info in to_build
        }
        for future in as_completed(futures):
            year = futures[future]
            result = future.result()
            results[year] = result
            print(
                f"{year}: built rows={result['ingestion']['rows_parsed']:,} "
                f"target={result['ingestion']['target_rows_in_event_window_2015_2023']:,} "
                f"output={result['ingestion']['grouped_output_count_after_aggregation']:,}",
                flush=True,
            )

    if args.no_combine:
        print("Per-source-year extraction complete; no combined all-source output was written.", flush=True)
        return
    if years != list(ALL_SOURCE_YEARS):
        raise ValueError("Combined outputs require all source years 2015-2023; use --no-combine for subset smoke runs")

    combined = merge_source_caches(years, catalog)
    all_sources = [results[year] for year in years]
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "specification_path": str(specification_path),
        "specification_sha256": sha256_file(specification_path),
        "specification_locked_before_effects": True,
        "field_role_and_origin_catalog": catalog,
        "source_year_audits": all_sources,
        "combined_source_audit": combined,
        "privacy": {
            "patient_level_rows_emitted": False,
            "identifiers_selected": False,
            "only_streamed_fields": ["INGRE", "EGRESO", "ENTIDAD", "AFECPRIN", "PROCED", "TIPSERV", "DIAS_ESTA"],
        },
    }
    audit_path = ADMIN_DIR / "mx_origin_extraction_audit.json"
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"combined state-day reconciliation: {combined['all_origin_reconciliation_vs_existing_complete_source']['status']}", flush=True)
    print(f"wrote {audit_path}", flush=True)


if __name__ == "__main__":
    main()
