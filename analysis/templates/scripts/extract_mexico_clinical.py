#!/usr/bin/env python3
"""Build grouped Mexico SAEH clinical aggregates for stone outcomes.

This script reads only EGRESOS source tables and writes grouped counts. It
never emits record identifiers or patient-level rows. Completed year caches
are reused when the source SHA-256 and cache version still match.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
REPAIR = ROOT / "work" / "repair_20261112"
DATA_DIR = REPAIR / "data"
ADMIN_DIR = REPAIR / "admin"
BASE = Path(r"D:\数据库集合\新增数据\Mexico_DGIS\data_full")
RECOVERED_2019 = Path(
    r"D:\数据库集合\新增数据\Mexico_DGIS\data_test\2019\Egresos2019\EGRESOS_2019.txt"
)
RECOVERED_2019_FIELD_AUDIT = Path(
    r"D:\数据库集合\新增数据\Mexico_DGIS\field_audit_2019.md"
)
EXTENSION_2023_FIELD_AUDIT = BASE / "2023" / "field_audit_2023.md"
PANEL = REPAIR / "data" / "source" / "mex_state_day_stone.csv"

PRIMARY_YEARS = [2015, 2016, 2017, 2018, 2020, 2021, 2022]
RECOVERED_YEARS = [2019]
EXTENSION_YEARS = [2023]
YEARS = sorted(PRIMARY_YEARS + RECOVERED_YEARS + EXTENSION_YEARS)
STATE_CODES = [f"{i:02d}" for i in range(1, 33)]
CHUNK_SIZE = 200_000
CACHE_VERSION = 2
SEPARATOR_CANDIDATES = ("|", ";", ",", "\t")

CANONICAL_FIELDS = {
    "INGRE": {"INGRE"},
    "ENTIDAD": {"ENTIDAD"},
    "AFECPRIN": {"AFECPRIN"},
    "EDAD": {"EDAD"},
    "CVEEDAD": {"CVEEDAD", "CVE_EDAD"},
    "SEXO": {"SEXO"},
    "CLUES": {"CLUES"},
    "CEDOCVE": {"CEDOCVE", "CEDOC_VE"},
}
REQUIRED_FIELDS = ("INGRE", "ENTIDAD", "AFECPRIN", "EDAD", "CVEEDAD", "SEXO")
ENDPOINTS = {"K80", "K81", "N20", "N21", "N22", "N23"}
AGE_GROUPS = ("<18", "18-44", "45-64", "65plus", "unknown")
SEX_GROUPS = ("male", "female", "unknown")
LONG_COLUMNS = ["state", "date", "endpoint", "fullcode", "age_group", "sex", "n"]
WIDE_COUNT_COLUMNS = [
    "k80",
    "k81",
    "n20_only",
    "n21_only",
    "n22_other_disease",
    "n23_unspecified_renal_colic",
    "n20_n21_confirmed",
    "k80_k800",
    "k80_k803",
    "k80_acute_k800_k803",
]
WIDE_COLUMNS = ["state", "date", *WIDE_COUNT_COLUMNS]


def canonical_header(value: str) -> str:
    value = value.replace("\ufeff", "").strip().upper()
    return re.sub(r"[^A-Z0-9]", "", value)


def source_year_path(year: int) -> Path:
    if year == 2019:
        if not RECOVERED_2019.is_file():
            raise FileNotFoundError(RECOVERED_2019)
        return RECOVERED_2019
    year_dir = BASE / str(year)
    candidates = []
    for root, _dirs, files in os.walk(year_dir):
        for filename in files:
            low = filename.lower()
            if "egreso" in low and low.endswith((".csv", ".txt")) and not low.startswith("defunc"):
                candidates.append(Path(root) / filename)
    if not candidates:
        raise FileNotFoundError(f"No EGRESO CSV/TXT source found under {year_dir}")
    candidates.sort(key=lambda item: item.stat().st_size, reverse=True)
    return candidates[0]


def source_encoding(year: int) -> str:
    if year == 2019:
        return "utf-8-sig"
    if year == 2023:
        return "utf-8"
    return "latin-1"


def inspect_source(path: Path, year: int) -> dict[str, Any]:
    encoding = source_encoding(year)
    with path.open("r", encoding=encoding, errors="replace", newline="") as handle:
        header_line = handle.readline()
    sep = next((s for s in SEPARATOR_CANDIDATES if header_line.count(s) > 5), "|")
    raw_headers = next(csv.reader([header_line], delimiter=sep))
    chosen: dict[str, str] = {}
    normalized = {canonical_header(name): name for name in raw_headers}
    for canonical, aliases in CANONICAL_FIELDS.items():
        for alias in aliases:
            matched = normalized.get(canonical_header(alias))
            if matched is not None:
                chosen[canonical] = matched
                break
    missing = [field for field in REQUIRED_FIELDS if field not in chosen]
    if missing:
        raise ValueError(f"{year}: missing required source fields {missing}; headers={raw_headers!r}")
    return {
        "encoding": encoding,
        "delimiter": sep,
        "header_fields": raw_headers,
        "selected_fields": chosen,
        "missing_optional_fields": [field for field in ("CLUES", "CEDOCVE") if field not in chosen],
        "raw_age_unit_field": chosen["CVEEDAD"],
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


def normalize_category(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.replace(r"\.0+$", "", regex=True)


def assign_age_group(age: pd.Series, age_unit: pd.Series) -> pd.Series:
    age_num = pd.to_numeric(age, errors="coerce")
    unit = normalize_category(age_unit)
    valid_age = age_num.notna() & (age_num >= 0) & (age_num != 999)
    result = pd.Series("unknown", index=age.index, dtype="object")
    result.loc[valid_age & unit.isin(["0", "1", "2"])] = "<18"
    years = valid_age & unit.eq("3")
    result.loc[years & age_num.lt(18)] = "<18"
    result.loc[years & age_num.ge(18) & age_num.le(44)] = "18-44"
    result.loc[years & age_num.ge(45) & age_num.le(64)] = "45-64"
    result.loc[years & age_num.ge(65)] = "65plus"
    return result


def assign_sex(series: pd.Series) -> pd.Series:
    raw = normalize_category(series)
    return raw.map({"1": "male", "2": "female"}).fillna("unknown")


def empty_counter() -> Counter[str]:
    return Counter()


def year_cache_paths(year: int) -> tuple[Path, Path, Path]:
    long_cache = DATA_DIR / f"mx_clinical_year_{year}.csv.gz"
    wide_cache = DATA_DIR / f"mx_clinical_state_day_{year}.csv.gz"
    audit_cache = ADMIN_DIR / f"mx_clinical_year_{year}.json"
    return long_cache, wide_cache, audit_cache


def cache_is_current(source_hash: str, paths: tuple[Path, Path, Path]) -> dict[str, Any] | None:
    long_cache, wide_cache, audit_cache = paths
    if not (long_cache.is_file() and wide_cache.is_file() and audit_cache.is_file()):
        return None
    try:
        audit = json.loads(audit_cache.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if audit.get("cache_version") != CACHE_VERSION or audit.get("source", {}).get("sha256") != source_hash:
        return None
    return audit
    return audit


def build_year(year: int, path: Path, source_hash: str, source_info: dict[str, Any]) -> dict[str, Any]:
    long_cache, wide_cache, audit_cache = year_cache_paths(year)
    temp_db = DATA_DIR / f"mx_clinical_work_{year}.sqlite"
    if temp_db.exists():
        temp_db.unlink()

    chosen = source_info["selected_fields"]
    selected_raw = set(chosen.values())
    canonical_by_raw = {v: k for k, v in chosen.items()}
    totals: dict[int, Counter[str]] = {}
    age_unit_counts: Counter[str] = empty_counter()
    sex_code_counts: Counter[str] = empty_counter()
    source_state_counts: Counter[str] = empty_counter()
    unmatched_state_counts: Counter[str] = empty_counter()
    unmatched_diagnosis_state_counts: Counter[str] = empty_counter()
    all_rows = 0
    invalid_dates = 0
    invalid_stone_dates = 0
    selected_diagnosis_rows = 0
    valid_date_stone_rows = 0
    valid_output_rows = 0
    outside_analysis_calendar_rows = 0
    age999_rows = 0
    clues_nonempty_rows = 0

    conn = sqlite3.connect(temp_db)
    conn.execute(
        """CREATE TABLE agg (
               state TEXT NOT NULL,
               date TEXT NOT NULL,
               endpoint TEXT NOT NULL,
               fullcode TEXT NOT NULL,
               age_group TEXT NOT NULL,
               sex TEXT NOT NULL,
               n INTEGER NOT NULL,
               PRIMARY KEY(state,date,endpoint,fullcode,age_group,sex)
           )"""
    )
    upsert = (
        "INSERT INTO agg(state,date,endpoint,fullcode,age_group,sex,n) VALUES(?,?,?,?,?,?,?) "
        "ON CONFLICT(state,date,endpoint,fullcode,age_group,sex) DO UPDATE SET n=n+excluded.n"
    )
    try:
        chunks = pd.read_csv(
            path,
            sep=source_info["delimiter"],
            usecols=lambda col: col in selected_raw,
            dtype=str,
            chunksize=CHUNK_SIZE,
            encoding=source_info["encoding"],
            on_bad_lines="skip",
            low_memory=False,
            keep_default_na=False,
        )
        for chunk in chunks:
            chunk = chunk.rename(columns=canonical_by_raw)
            all_rows += len(chunk)
            age_unit_counts.update(normalize_category(chunk["CVEEDAD"]).replace("", "__blank__").value_counts().to_dict())
            sex_code_counts.update(normalize_category(chunk["SEXO"]).replace("", "__blank__").value_counts().to_dict())
            age999_rows += int(pd.to_numeric(chunk["EDAD"], errors="coerce").eq(999).sum())
            if "CLUES" in chunk.columns:
                clues_nonempty_rows += int(chunk["CLUES"].fillna("").astype(str).str.strip().ne("").sum())

            date = pd.to_datetime(chunk["INGRE"], errors="coerce")
            date_ok = date.notna()
            invalid_dates += int((~date_ok).sum())
            date_iso = date.dt.strftime("%Y-%m-%d")

            original_dx = chunk["AFECPRIN"].fillna("").astype(str).str.upper().str.strip()
            legacy_prefix = original_dx.str[:3]
            fullcode = original_dx.str.replace(r"[^A-Z0-9]", "", regex=True)
            endpoint = fullcode.str[:3]
            stone = endpoint.isin(ENDPOINTS)
            selected_diagnosis_rows += int(stone.sum())
            invalid_stone_dates += int((stone & ~date_ok).sum())
            valid_date_stone_rows += int((stone & date_ok).sum())

            raw_state = chunk["ENTIDAD"].fillna("").astype(str).str.strip()
            state = raw_state.str.zfill(2)
            state_ok = state.isin(STATE_CODES)
            source_state_counts.update(raw_state.replace("", "__blank__").value_counts().to_dict())
            if (~state_ok).any():
                unmatched_state_counts.update(
                    raw_state.loc[~state_ok].replace("", "__blank__").value_counts().to_dict()
                )
                unmatched_stone = ~state_ok & stone
                if unmatched_stone.any():
                    unmatched_diagnosis_state_counts.update(
                        raw_state.loc[unmatched_stone].replace("", "__blank__").value_counts().to_dict()
                    )

            # Match the historical panel's date intersection. Primary source
            # files can contain admissions dated in another primary year, so
            # retain all primary calendar years and combine the source caches
            # later. Recovered 2019 remains a separate calendar-year arm.
            if year == 2019:
                allowed_years = RECOVERED_YEARS
            elif year == 2023:
                allowed_years = PRIMARY_YEARS + EXTENSION_YEARS
            else:
                allowed_years = PRIMARY_YEARS
            in_analysis_calendar = date.dt.year.isin(allowed_years)
            outside_analysis_calendar_rows += int((stone & date_ok & ~in_analysis_calendar).sum())
            keep = stone & date_ok & state_ok & in_analysis_calendar
            legacy_valid = keep
            if legacy_valid.any():
                by_calendar_year = date.dt.year.loc[legacy_valid]
                for calendar_year, idx in by_calendar_year.groupby(by_calendar_year).groups.items():
                    prefix = legacy_prefix.loc[idx]
                    year_totals = totals.setdefault(int(calendar_year), Counter())
                    year_totals["k80"] += int(prefix.eq("K80").sum())
                    year_totals["k81"] += int(prefix.eq("K81").sum())
                    year_totals["n20_collapsed_n20_to_n23"] += int(
                        prefix.isin(["N20", "N21", "N22", "N23"]).sum()
                    )
            if not keep.any():
                continue
            selected = pd.DataFrame(
                {
                    "state": state.loc[keep].to_numpy(),
                    "date": date_iso.loc[keep].to_numpy(),
                    "endpoint": endpoint.loc[keep].to_numpy(),
                    "fullcode": fullcode.loc[keep].to_numpy(),
                    "age_group": assign_age_group(chunk.loc[keep, "EDAD"], chunk.loc[keep, "CVEEDAD"]).to_numpy(),
                    "sex": assign_sex(chunk.loc[keep, "SEXO"]).to_numpy(),
                }
            )
            grouped = selected.groupby(LONG_COLUMNS[:-1], sort=False, observed=True).size().reset_index(name="n")
            conn.executemany(
                upsert,
                [
                    (row.state, row.date, row.endpoint, row.fullcode, row.age_group, row.sex, int(row.n))
                    for row in grouped.itertuples(index=False)
                ],
            )
            conn.commit()
            valid_output_rows += len(selected)

        long_cache.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(long_cache, "wt", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(LONG_COLUMNS)
            cursor = conn.execute(
                "SELECT state,date,endpoint,fullcode,age_group,sex,n FROM agg "
                "ORDER BY state,date,endpoint,fullcode,age_group,sex"
            )
            while True:
                rows = cursor.fetchmany(100_000)
                if not rows:
                    break
                writer.writerows(rows)

        wide_observed = pd.read_sql_query(
            """SELECT state,date,
                      SUM(CASE WHEN endpoint='K80' THEN n ELSE 0 END) AS k80,
                      SUM(CASE WHEN endpoint='K81' THEN n ELSE 0 END) AS k81,
                      SUM(CASE WHEN endpoint='N20' THEN n ELSE 0 END) AS n20_only,
                      SUM(CASE WHEN endpoint='N21' THEN n ELSE 0 END) AS n21_only,
                      SUM(CASE WHEN endpoint='N22' THEN n ELSE 0 END) AS n22_other_disease,
                      SUM(CASE WHEN endpoint='N23' THEN n ELSE 0 END) AS n23_unspecified_renal_colic,
                      SUM(CASE WHEN endpoint IN ('N20','N21') THEN n ELSE 0 END) AS n20_n21_confirmed,
                      SUM(CASE WHEN endpoint='K80' AND substr(fullcode,1,4)='K800' THEN n ELSE 0 END) AS k80_k800,
                      SUM(CASE WHEN endpoint='K80' AND substr(fullcode,1,4)='K803' THEN n ELSE 0 END) AS k80_k803,
                      SUM(CASE WHEN endpoint='K80' AND substr(fullcode,1,4) IN ('K800','K803') THEN n ELSE 0 END) AS k80_acute_k800_k803
                 FROM agg GROUP BY state,date ORDER BY state,date""",
            conn,
        )
        if wide_observed.empty:
            wide_observed = pd.DataFrame(columns=WIDE_COLUMNS)
        wide_observed.to_csv(wide_cache, index=False, compression="gzip")

        actual = {column: int(wide_observed[column].sum()) for column in WIDE_COUNT_COLUMNS}
        k80_subcodes = {}
        for row in conn.execute(
            "SELECT substr(fullcode,1,4),SUM(n) FROM agg WHERE endpoint='K80' GROUP BY substr(fullcode,1,4) ORDER BY 1"
        ):
            k80_subcodes[row[0]] = int(row[1])
        grouped_unknown = {
            "age_unknown_rows": int(
                conn.execute("SELECT COALESCE(SUM(n),0) FROM agg WHERE age_group='unknown'").fetchone()[0]
            ),
            "sex_unknown_rows": int(conn.execute("SELECT COALESCE(SUM(n),0) FROM agg WHERE sex='unknown'").fetchone()[0]),
            "both_age_and_sex_unknown_rows": int(
                conn.execute("SELECT COALESCE(SUM(n),0) FROM agg WHERE age_group='unknown' AND sex='unknown'").fetchone()[0]
            ),
        }
    finally:
        conn.close()
        if temp_db.exists():
            temp_db.unlink()

    if year == 2019:
        source_role = "recovered_additional_2019"
    elif year == 2023:
        source_role = "extension_2023_source"
    else:
        source_role = "primary_panel_year"
    audit = {
        "cache_version": CACHE_VERSION,
        "year": year,
        "source_role": source_role,
        "source": {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": source_hash,
            "encoding": source_info["encoding"],
            "delimiter": source_info["delimiter"],
            "field_count": len(source_info["header_fields"]),
            "field_names": source_info["header_fields"],
            "selected_fields": source_info["selected_fields"],
            "cveedad_header_present": "CVEEDAD" in source_info["selected_fields"],
            "raw_age_unit_field": source_info["raw_age_unit_field"],
            "clues_present": "CLUES" in source_info["selected_fields"],
            "cedocve_present": "CEDOCVE" in source_info["selected_fields"],
            "clues_nonempty_rows": clues_nonempty_rows,
        },
        "ingestion": {
            "rows_parsed": all_rows,
            "chunk_rows": CHUNK_SIZE,
            "bad_line_policy": "skip (same policy as historical S02 source builder; pandas C parser does not expose skipped-line count)",
            "date_field": "INGRE (admission date)",
            "invalid_or_missing_ingre_rows": invalid_dates,
            "selected_diagnosis_rows_K80_K81_N20_N23": selected_diagnosis_rows,
            "selected_diagnosis_rows_with_invalid_ingre": invalid_stone_dates,
            "selected_rows_with_valid_ingre": valid_date_stone_rows,
            "selected_rows_with_valid_ingre_and_matched_state_and_analysis_calendar": valid_output_rows,
            "selected_rows_with_valid_ingre_outside_analysis_calendar": outside_analysis_calendar_rows,
            "invalid_date_rule": "Invalid/missing INGRE rows are excluded from clinical aggregates.",
        },
        "field_quality": {
            "cveedad_code_counts_all_rows": dict(sorted(age_unit_counts.items())),
            "sex_code_counts_all_rows": dict(sorted(sex_code_counts.items())),
            "edad_999_rows_all_rows": age999_rows,
            "entidad_raw_counts_all_rows": dict(sorted(source_state_counts.items())),
            "unmatched_entidad_raw_counts_all_rows": dict(sorted(unmatched_state_counts.items())),
            "unmatched_entidad_raw_counts_selected_diagnoses": dict(sorted(unmatched_diagnosis_state_counts.items())),
        },
        "aggregation": {
            "grain": "state x admission-date x principal fullcode x age_group x sex",
            "principal_diagnosis_field": "AFECPRIN",
            "geography_field": "ENTIDAD (habitual residence state)",
            "state_codes_in_outputs": STATE_CODES,
            "unmatched_state_rule": "Only zero-padded numeric ENTIDAD codes 01-32 are emitted; all other raw values are excluded and retained in this audit count.",
            "diagnosis_normalization": "Uppercase AFECPRIN with non-alphanumeric separators removed; fullcode preserved at normalized source granularity.",
            "age_group_rule": "CVEEDAD 0/1/2 with EDAD 0-998 maps to <18; CVEEDAD 3 groups completed EDAD years into <18, 18-44, 45-64, 65plus; CVEEDAD 9, missing/other units, invalid/negative EDAD, and EDAD=999 map to unknown.",
            "sex_rule": "SEXO 1=male, 2=female; 9, blank, and all other values=unknown.",
            "age_sex_unknown_endpoint_rows": grouped_unknown,
            "k80_four_character_subcode_counts": k80_subcodes,
            "k80_acute_rule": "K800 (gallbladder calculus with acute cholecystitis) plus K803 (bile duct calculus with cholangitis), based on coding_evidence.json; kept separate and combined.",
            "n20_n21_confirmed_rule": "N20 plus N21 only; N22 and N23 remain separate sensitivity endpoints.",
            "TIPSERV_used": False,
            "zero_group_rule": "Per-source year caches contain observed positive grouped cells; final wide output combines source caches by admission date, crosses all 32 states with every analysis-calendar day, and fills missing endpoint counts with integer zero.",
            "patient_level_rows_emitted": False,
        },
        "legacy_panel_reproduction": {
            "historical_logic": "S02 uses AFECPRIN first three characters; it sums K80, K81, and maps each of N20/N21/N22/N23 to one historical N20 column.",
            "reconstructed_annual_totals_by_ingre_calendar_year": {
                str(calendar_year): {key: int(value) for key, value in sorted(values.items())}
                for calendar_year, values in sorted(totals.items())
            },
            "baseline_panel_comparison": "not_applicable_recovered_2019_excluded_from_primary_panel"
            if year == 2019
            else "pending_final_reconciliation",
        },
        "canonical_endpoint_totals": actual,
    }
    if year == 2019:
        audit["recovery_provenance"] = {
            "field_audit_path": str(RECOVERED_2019_FIELD_AUDIT),
            "field_audit_sha256": sha256_file(RECOVERED_2019_FIELD_AUDIT),
            "archive_zip_bytes": 85_756_134,
            "archive_zip_sha256_from_field_audit": "74a449ca69620c5e01adce4d719d03234575501184078512a4741b2398d46827",
            "raw_egresos_rows_expected_from_field_audit": 2_629_434,
            "raw_egresos_bytes_expected_from_field_audit": 562_533_021,
            "source_status": "Recovered additional SAEH-SSA 2019 file. This is not part of the original primary panel and is not used to infer missing exposure data.",
        }
    if year == 2023:
        audit["extension_provenance"] = {
            "field_audit_path": str(EXTENSION_2023_FIELD_AUDIT),
            "field_audit_sha256": sha256_file(EXTENSION_2023_FIELD_AUDIT),
            "archive_zip_bytes": 111_507_585,
            "archive_zip_sha256_from_field_audit": "8ad3e2e8b9c4e532e3f9c0ba2624f42f28864a7d4c278025658901987274b2c2",
            "raw_egresos_rows_expected_from_field_audit": 2_591_278,
            "raw_egresos_bytes_expected_from_field_audit": 467_960_465,
            "source_status": "Additional 2023 SAEH-SSA source. INGRE dates in primary-panel years contribute only to exact legacy S02 reconciliation; INGRE=2023 is emitted under the separate extension_2023 role.",
        }
    audit_cache.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return audit


def read_cached_audit(year: int, source_hash: str, source_info: dict[str, Any], path: Path) -> dict[str, Any]:
    paths = year_cache_paths(year)
    current = cache_is_current(source_hash, paths)
    if current is not None:
        return current
    return build_year(year, path, source_hash, source_info)


def baseline_panel_totals() -> dict[int, dict[str, int]]:
    panel = pd.read_csv(PANEL, dtype={"state": str})
    panel["date"] = pd.to_datetime(panel["date"], errors="coerce")
    panel["year"] = panel["date"].dt.year
    out: dict[int, dict[str, int]] = {}
    for year, frame in panel.groupby("year"):
        out[int(year)] = {
            "k80": int(frame["k80"].sum()),
            "k81": int(frame["k81"].sum()),
            "n20_collapsed_n20_to_n23": int(frame["n20"].sum()),
        }
    return out


def combine_long_caches(year_paths: list[tuple[int, Path]], output_path: Path) -> None:
    """Sum per-source-year aggregate caches to the requested state-day grain."""
    temp_db = DATA_DIR / "mx_clinical_work_combine.sqlite"
    if temp_db.exists():
        temp_db.unlink()
    conn = sqlite3.connect(temp_db)
    conn.execute(
        """CREATE TABLE agg (
               source_role TEXT NOT NULL,
               state TEXT NOT NULL,
               date TEXT NOT NULL,
               endpoint TEXT NOT NULL,
               fullcode TEXT NOT NULL,
               age_group TEXT NOT NULL,
               sex TEXT NOT NULL,
               n INTEGER NOT NULL,
               PRIMARY KEY(source_role,state,date,endpoint,fullcode,age_group,sex)
           )"""
    )
    upsert = (
        "INSERT INTO agg VALUES(?,?,?,?,?,?,?,?) "
        "ON CONFLICT(source_role,state,date,endpoint,fullcode,age_group,sex) "
        "DO UPDATE SET n=n+excluded.n"
    )
    try:
        for year, cache_path in year_paths:
            for frame in pd.read_csv(
                cache_path,
                dtype={"state": str, "date": str, "endpoint": str, "fullcode": str, "age_group": str, "sex": str},
                chunksize=100_000,
                compression="gzip",
            ):
                event_year = pd.to_numeric(frame["date"].str[:4], errors="coerce")
                frame = frame.loc[event_year.between(2015, 2023)].copy()
                if frame.empty:
                    continue
                if year == 2019:
                    frame.insert(0, "source_role", "recovered_additional_2019")
                elif year == 2023:
                    frame.insert(
                        0,
                        "source_role",
                        frame["date"].str.startswith("2023-").map(
                            {True: "extension_2023", False: "primary_panel_year"}
                        ),
                    )
                else:
                    frame.insert(0, "source_role", "primary_panel_year")
                conn.executemany(
                    upsert,
                    [
                        (row.source_role, row.state, row.date, row.endpoint, row.fullcode, row.age_group, row.sex, int(row.n))
                        for row in frame.itertuples(index=False)
                    ],
                )
                conn.commit()
        columns = ["source_role", *LONG_COLUMNS]
        with gzip.open(output_path, "wt", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            writer.writerow(columns)
            cursor = conn.execute(
                "SELECT source_role,state,date,endpoint,fullcode,age_group,sex,n FROM agg "
                "ORDER BY source_role,state,endpoint,fullcode,date,age_group,sex"
            )
            while True:
                rows = cursor.fetchmany(100_000)
                if not rows:
                    break
                writer.writerows(rows)
    finally:
        conn.close()
        if temp_db.exists():
            temp_db.unlink()


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ADMIN_DIR.mkdir(parents=True, exist_ok=True)
    audits = []
    for year in YEARS:
        path = source_year_path(year)
        source_info = inspect_source(path, year)
        source_hash = sha256_file(path)
        cache_paths = year_cache_paths(year)
        was_cached = cache_is_current(source_hash, cache_paths) is not None
        audit = read_cached_audit(year, source_hash, source_info, path)
        audit["cache_reused"] = was_cached
        if year == 2019:
            audit["source_crossyear_probe"] = {
                "method": "Independent streaming count of AFECPRIN prefix and INGRE in the 2019 source; state match uses ENTIDAD 01-32.",
                "valid_date_selected_stone_rows_outside_ingre_2019_by_year": {
                    "2014": 1,
                    "2016": 1,
                    "2017": 1,
                    "2018": 327,
                },
                "same_rows_with_matched_state_01_32_by_year": {
                    "2014": 1,
                    "2016": 1,
                    "2017": 1,
                    "2018": 325,
                },
                "same_rows_with_unmatched_state_by_year": {"2018": 2},
                "diagnosis_prefix_counts_across_these_rows": {
                    "K80": 181,
                    "K81": 120,
                    "N20": 21,
                    "N21": 3,
                    "N22": 1,
                    "N23": 4,
                },
                "included_in_current_2019_recovered_cache": False,
                "note": "These cross-year rows remain distinguishable in this audit. The primary cache holds INGRE=2019 from the recovered file only; an all-source 2019 fill requires a separate source-cohort reconstruction.",
            }
            cache_paths[2].write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        audits.append(audit)
        print(
            f"{year}: {'cached' if was_cached else 'built'} rows={audit['ingestion']['rows_parsed']:,} "
            f"selected={audit['ingestion'].get('selected_rows_with_valid_ingre_and_matched_state_and_analysis_calendar', audit['ingestion'].get('selected_rows_with_valid_ingre_matched_state_and_analysis_calendar', 0)):,} "
            f"source={path}",
            flush=True,
        )

    long_paths = [(year, year_cache_paths(year)[0]) for year in YEARS]
    wide_paths = [(year, year_cache_paths(year)[1]) for year in YEARS]
    long_output = DATA_DIR / "mx_clinical_long.csv.gz"
    combine_long_caches(long_paths, long_output)

    wide_parts = []
    for year, path in wide_paths:
        frame = pd.read_csv(path, dtype={"state": str, "date": str})
        if year == 2019:
            frame.insert(0, "source_role", "recovered_additional_2019")
        elif year == 2023:
            frame.insert(
                0,
                "source_role",
                frame["date"].str.startswith("2023-").map(
                    {True: "extension_2023", False: "primary_panel_year"}
                ),
            )
        else:
            frame.insert(0, "source_role", "primary_panel_year")
        wide_parts.append(frame)
    observed = pd.concat(wide_parts, ignore_index=True)
    observed = observed.groupby(["source_role", "state", "date"], as_index=False)[WIDE_COUNT_COLUMNS].sum()
    wide_grids = []
    for year in YEARS:
        if year == 2019:
            role = "recovered_additional_2019"
        elif year == 2023:
            role = "extension_2023"
        else:
            role = "primary_panel_year"
        days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D").strftime("%Y-%m-%d")
        grid = pd.MultiIndex.from_product([STATE_CODES, days], names=["state", "date"]).to_frame(index=False)
        part = observed.loc[(observed["source_role"] == role) & observed["date"].str.startswith(f"{year}-")]
        grid.insert(0, "source_role", role)
        wide = grid.merge(part, on=["source_role", "state", "date"], how="left")
        for column in WIDE_COUNT_COLUMNS:
            wide[column] = wide[column].fillna(0).astype("int64")
        wide_grids.append(wide[["source_role", *WIDE_COLUMNS]])
    wide_all = pd.concat(wide_grids, ignore_index=True)
    wide_all = wide_all.sort_values(["state", "date"], kind="stable").reset_index(drop=True)
    wide_output = DATA_DIR / "mx_clinical_state_day.csv"
    wide_all.to_csv(wide_output, index=False)

    panel = baseline_panel_totals()
    reconstructed_by_calendar_year: dict[str, Counter[str]] = {}
    for audit in audits:
        source_totals = audit["legacy_panel_reproduction"]["reconstructed_annual_totals_by_ingre_calendar_year"]
        for calendar_year, year_totals in source_totals.items():
            accumulated = reconstructed_by_calendar_year.setdefault(calendar_year, Counter())
            accumulated.update({key: int(value) for key, value in year_totals.items()})
    reconciliation = {}
    for year in PRIMARY_YEARS:
        reconstructed = dict(reconstructed_by_calendar_year.get(str(year), Counter()))
        baseline = panel.get(year, {})
        comparisons = {
            key: {
                "reconstructed": reconstructed.get(key, 0),
                "baseline_panel": baseline.get(key),
                "difference": None if baseline.get(key) is None else reconstructed.get(key, 0) - baseline[key],
                "matches": None if baseline.get(key) is None else reconstructed.get(key, 0) == baseline[key],
            }
            for key in ("k80", "k81", "n20_collapsed_n20_to_n23")
        }
        reconciliation[str(year)] = comparisons
    reconciliation["2019"] = {
        "status": "recovered_additional_2019_excluded_from_primary_panel_comparison",
        "reconstructed_clinical_totals": dict(reconstructed_by_calendar_year.get("2019", Counter())),
    }
    audit_2023 = next(audit for audit in audits if audit["year"] == 2023)
    source_2023_totals = audit_2023["legacy_panel_reproduction"]["reconstructed_annual_totals_by_ingre_calendar_year"]
    reconciliation["2023"] = {
        "status": "extension_2023_actual_admission_dates_not_in_primary_panel",
        "reconstructed_clinical_totals": dict(reconstructed_by_calendar_year.get("2023", Counter())),
        "source_2023_contributions_to_primary_calendar_years": {
            year: source_2023_totals.get(year, {}) for year in map(str, PRIMARY_YEARS)
        },
    }
    for audit in audits:
        audit["legacy_panel_reproduction"]["baseline_panel_comparison"] = (
            "not_applicable_recovered_2019_excluded_from_primary_panel"
            if audit["year"] == 2019
            else "see_top_level_annual_reconciliation_combined_across_primary_source_files"
        )

    audit_output = ADMIN_DIR / "mx_clinical_audit.json"
    audit_document = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "primary_panel_years": PRIMARY_YEARS,
        "additional_recovered_years": RECOVERED_YEARS,
        "extension_years": EXTENSION_YEARS,
        "years_in_long_and_wide_outputs": YEARS,
        "years_excluded": {
            "2019_primary": "not in original primary panel; recovered separately from historical SAEH archive",
            "2023_primary": "actual INGRE=2023 records are a separate extension; the 2023 source file's earlier INGRE dates in primary years are included only to reproduce the historical S02 source selection",
        },
        "source_builder_reference": r"F:\分析专用\climate_attribution_20261106\work\newproj\S02_mexico_rebuild.py",
        "coding_reference": str(REPAIR / "references" / "coding_evidence.json"),
        "baseline_panel_path": str(PANEL),
        "outputs": {
            "long_grouped_csv_gzip": str(long_output),
            "state_day_wide_csv": str(wide_output),
            "per_source_year_long_caches": {
                str(year): str(year_cache_paths(year)[0]) for year in YEARS
            },
            "per_source_year_state_day_caches": {
                str(year): str(year_cache_paths(year)[1]) for year in YEARS
            },
        },
        "zero_group_rules": {
            "long": "Only observed positive grouped cells are emitted. Unknown age and sex are explicit categories and are not dropped.",
            "state_day_wide": "Every output year has all state codes 01-32 crossed with all calendar days; no observed group is zero-filled to integer 0 after a left join.",
            "unmatched_states": "ENTIDAD outside 01-32 is excluded from outcomes and retained in yearly audit counts.",
            "invalid_dates": "Missing/unparseable INGRE records are excluded from outcomes and counted in yearly audits.",
        },
        "year_audits": audits,
        "annual_reconciliation": reconciliation,
        "reconciliation_summary": {
            "all_primary_year_metrics_match": all(
                value["matches"] is True
                for by_year in reconciliation.values()
                if isinstance(by_year, dict) and "status" not in by_year
                for value in by_year.values()
            ),
            "interpretation": "The old panel's n20 field includes N20+N21+N22+N23; new outputs keep these endpoints separate. 2019 is not compared because it was absent from the primary panel. INGRE=2023 outcomes are an extension; 2023-source records with earlier primary-panel INGRE dates are included in the legacy reconciliation to mirror S02.",
        },
    }
    audit_output.write_text(json.dumps(audit_document, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {long_output}", flush=True)
    print(f"wrote {wide_output} ({len(wide_all):,} state-days)", flush=True)
    print(f"wrote {audit_output}", flush=True)
    for year, result in reconciliation.items():
        if "status" in result:
            print(f"{year}: {result['status']}", flush=True)
        else:
            matches = all(cell["matches"] is True for cell in result.values())
            print(f"{year}: legacy panel reconciliation {'MATCH' if matches else 'MISMATCH'} {result}", flush=True)


if __name__ == "__main__":
    main()
