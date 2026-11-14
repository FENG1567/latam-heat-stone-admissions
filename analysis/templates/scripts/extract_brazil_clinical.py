#!/usr/bin/env python3
"""Build deidentified regional SIH-RD clinical aggregates for the repair run.

Only selected aggregate files and aggregate audit JSON are durable outputs.
Temporary HMACs of N_AIH and full source rows are used for duplicate auditing,
never written to final outputs, and removed after a successful run.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import hmac
import json
import os
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Keep each worker single-threaded inside numerical libraries.
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_name] = "1"

import pandas as pd  # noqa: E402


SOURCE_ROOT = Path("/data/admin/fengshuo/readmit_pred_20261104/work/sih_csv")
CROSSWALK = Path(
    "/data/admin/fengshuo/brazil_biliary_continuation_20261009/original15role_main_v1/"
    "efbe05b9767d421b99d9d6acc6281ca1/candidate_inputs/primary_fixed_support_v1/"
    "reference_step03_fixed_region_crosswalk.csv"
)
REMOTE_ROOT = Path("/data/admin/fengshuo/latam_stone_repair_20261112")
YEARS = tuple(range(2014, 2022))
TARGET_DX3 = {"K80", "K81", "N20", "N21", "N22", "N23"}
CLINICAL_DIMS = [
    "hrid",
    "admission_date",
    "diag_princ4",
    "sexo_raw",
    "age_group",
    "car_int_raw",
    "car_int_group",
]
AGG_DIMS = ["source_year", *CLINICAL_DIMS]
DATE_FORMATS = (
    ("yyyy-mm-dd", "%Y-%m-%d"),
    ("dd/mm/yyyy", "%d/%m/%Y"),
    ("yyyy/mm/dd", "%Y/%m/%d"),
    ("dd-mm-yyyy", "%d-%m-%Y"),
)
_CROSSWALK: dict[str, str] = {}
_HASH_KEY = b""


def _json_write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp-{os.getpid()}")
    with tmp.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, path)


def _csv_gz_write(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp-{os.getpid()}.gz")
    frame.to_csv(tmp, index=False, encoding="utf-8", compression="gzip")
    os.replace(tmp, path)


def _clean_code(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)


def _normalize_municipality(series: pd.Series) -> pd.Series:
    digits = series.astype(str).str.strip().str.replace(r"\.0+$", "", regex=True)
    digits = digits.str.replace(r"\D", "", regex=True)
    return digits.where(digits.eq(""), digits.str[-6:].str.zfill(6))


def _load_crosswalk(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    required = {"municipality_sih6", "health_region_2023"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Crosswalk is missing required columns: {sorted(missing)}")
    frame["municipality_key"] = _normalize_municipality(frame["municipality_sih6"])
    frame["hrid_value"] = frame["health_region_2023"].astype(str).str.strip()
    frame = frame[frame["municipality_key"].ne("")]
    conflicts = frame.groupby("municipality_key", dropna=False)["hrid_value"].nunique()
    conflicts = conflicts[conflicts.gt(1)]
    if not conflicts.empty:
        raise ValueError(f"Crosswalk has {len(conflicts)} conflicting municipality mappings")
    return dict(zip(frame["municipality_key"], frame["hrid_value"]))


def _age_group(age: pd.Series, unit: pd.Series) -> pd.Series:
    age_n = pd.to_numeric(age.astype(str).str.strip(), errors="coerce")
    unit_n = pd.to_numeric(unit.astype(str).str.strip(), errors="coerce")
    result = pd.Series("unknown", index=age.index, dtype="object")

    # COD_IDADE 2=days and 3=months. Valid values are necessarily under 18.
    result.loc[unit_n.eq(2) & age_n.between(0, 30)] = "<18"
    result.loc[unit_n.eq(3) & age_n.between(0, 11)] = "<18"

    valid_years = unit_n.eq(4) & age_n.between(0, 99)
    result.loc[valid_years & age_n.lt(18)] = "<18"
    result.loc[valid_years & age_n.between(18, 44)] = "18-44"
    result.loc[valid_years & age_n.between(45, 64)] = "45-64"
    result.loc[valid_years & age_n.ge(65)] = "65plus"

    # Dictionary code 5 denotes 100+ years; do not interpret it as five years.
    result.loc[unit_n.eq(5)] = "65plus"
    return result


def _car_group(raw: pd.Series) -> pd.Series:
    clean = raw.astype(str).str.strip().str.upper()
    missing = clean.isin({"", "0", "00", "NA", "N/A", "NAN", "NULL"})
    normalized = clean.where(~clean.str.fullmatch(r"\d{1,2}"), clean.str.zfill(2))
    out = pd.Series("other", index=raw.index, dtype="object")
    out.loc[normalized.eq("01")] = "01"
    out.loc[normalized.eq("02")] = "02"
    out.loc[missing] = "unknown"
    return out


def _parse_dates(series: pd.Series) -> tuple[pd.Series, pd.Series]:
    raw = series.astype(str).str.strip()
    parsed = pd.Series(pd.NaT, index=series.index, dtype="datetime64[ns]")
    fmt = pd.Series("invalid", index=series.index, dtype="object")
    empty = raw.eq("") | raw.str.upper().isin({"NA", "N/A", "NAN", "NULL"})
    fmt.loc[empty] = "empty"

    eight = raw.str.fullmatch(r"\d{8}") & ~empty
    if eight.any():
        ymd = raw.str.slice(0, 4).between("1900", "2100")
        ymd_mask = eight & ymd
        dmy_mask = eight & ~ymd
        parsed.loc[ymd_mask] = pd.to_datetime(raw.loc[ymd_mask], format="%Y%m%d", errors="coerce")
        parsed.loc[dmy_mask] = pd.to_datetime(raw.loc[dmy_mask], format="%d%m%Y", errors="coerce")
        fmt.loc[ymd_mask] = "yyyymmdd"
        fmt.loc[dmy_mask] = "ddmmyyyy"

    for label, date_format in DATE_FORMATS:
        mask = parsed.isna() & ~empty & raw.str.contains(r"[-/]", regex=True)
        if mask.any():
            candidate = pd.to_datetime(raw.loc[mask], format=date_format, errors="coerce")
            matched = candidate.notna()
            if matched.any():
                idx = candidate.index[matched]
                parsed.loc[idx] = candidate.loc[idx]
                fmt.loc[idx] = label

    fmt.loc[parsed.isna() & ~empty] = "invalid"
    return parsed, fmt


def _hmac(text: str) -> str:
    return hmac.new(_HASH_KEY, text.encode("utf-8", errors="replace"), hashlib.sha256).hexdigest()


def _worker_init(crosswalk: dict[str, str], hash_key: bytes) -> None:
    global _CROSSWALK, _HASH_KEY
    _CROSSWALK = crosswalk
    _HASH_KEY = hash_key


def _paths_for_source(source: Path, output_root: Path) -> tuple[Path, Path, Path]:
    source_year = source.parent.parent.name
    uf = source.parent.name
    stem = source.name.removesuffix(".csv.gz")
    base = output_root / "intermediate" / "per_file" / source_year / uf / stem
    agg_path = base.with_suffix(".agg.csv.gz")
    audit_path = base.with_suffix(".audit.json")
    id_path = output_root / "intermediate" / "private_id_hash" / source_year / uf / (stem + ".csv.gz")
    return agg_path, audit_path, id_path


def _process_file(source: Path, output_root: Path, chunksize: int) -> dict[str, Any]:
    source_year = int(source.parent.parent.name)
    agg_path, audit_path, id_path = _paths_for_source(source, output_root)
    size = source.stat().st_size
    mtime_ns = source.stat().st_mtime_ns

    if agg_path.exists() and audit_path.exists() and id_path.exists():
        try:
            old = json.loads(audit_path.read_text(encoding="utf-8"))
            if old.get("source_size_bytes") == size and old.get("source_mtime_ns") == mtime_ns:
                return {"status": "skipped", "source_year": source_year, "source_file": source.name}
        except (OSError, json.JSONDecodeError):
            pass

    started = time.time()
    total_rows = 0
    selected_rows = 0
    selected_missing_naih = 0
    unknown_muni = 0
    selected_unknown_muni = 0
    admission_year_counts: Counter[str] = Counter()
    selected_admission_year_counts: Counter[str] = Counter()
    date_format_counts: Counter[str] = Counter()
    age_code_counts: Counter[str] = Counter()
    age_value_counts: Counter[str] = Counter()
    sex_counts: Counter[str] = Counter()
    car_counts: Counter[str] = Counter()
    principal_base_counts: Counter[str] = Counter()
    selected_diag4_counts: Counter[str] = Counter()
    aggregate_parts: list[pd.DataFrame] = []
    id_parts: list[pd.DataFrame] = []
    columns: list[str] | None = None

    for chunk in pd.read_csv(
        source,
        dtype=str,
        keep_default_na=False,
        encoding="latin1",
        compression="gzip",
        chunksize=chunksize,
        on_bad_lines="error",
    ):
        if columns is None:
            columns = list(chunk.columns)
            required = {
                "N_AIH",
                "DT_INTER",
                "DIAG_PRINC",
                "IDADE",
                "COD_IDADE",
                "SEXO",
                "MUNIC_RES",
                "CAR_INT",
            }
            missing = required - set(columns)
            if missing:
                raise ValueError(f"{source.name} missing required columns: {sorted(missing)}")
        total_rows += len(chunk)

        # Aggregate raw dictionary/value code distributions without writing row-level data.
        for col, counter in (("COD_IDADE", age_code_counts), ("SEXO", sex_counts), ("CAR_INT", car_counts)):
            vals = chunk[col].astype(str).str.strip().replace("", "unknown")
            counter.update(vals.value_counts(dropna=False).to_dict())
        raw_age = chunk["IDADE"].astype(str).str.strip().replace("", "unknown")
        raw_unit = chunk["COD_IDADE"].astype(str).str.strip().replace("", "unknown")
        age_pairs = raw_unit + "|" + raw_age
        age_value_counts.update(age_pairs.value_counts(dropna=False).to_dict())

        dx_clean = _clean_code(chunk["DIAG_PRINC"])
        dx_base = dx_clean.str.slice(0, 3)
        dx_valid = dx_clean.str.fullmatch(r"[A-Z]\d{2}\d?").fillna(False)
        dx_base_audit = dx_base.where(dx_valid, "invalid_or_unknown")
        principal_base_counts.update(dx_base_audit.value_counts(dropna=False).to_dict())

        dates, formats = _parse_dates(chunk["DT_INTER"])
        date_format_counts.update(formats.value_counts(dropna=False).to_dict())
        admission_year = dates.dt.year.astype("Int64").astype(str).replace("<NA>", "unknown")
        admission_year_counts.update(admission_year.value_counts(dropna=False).to_dict())

        muni = _normalize_municipality(chunk["MUNIC_RES"])
        mapped = muni.map(_CROSSWALK)
        unknown_mask = mapped.isna() | mapped.astype(str).str.strip().eq("")
        unknown_muni += int(unknown_mask.sum())

        select_mask = dx_valid & dx_base.isin(TARGET_DX3)
        if not select_mask.any():
            continue

        selected_rows += int(select_mask.sum())
        selected = chunk.loc[select_mask].copy()
        selected_dates = dates.loc[select_mask]
        selected_years = admission_year.loc[select_mask]
        selected_muni_unknown = unknown_mask.loc[select_mask]
        selected_unknown_muni += int(selected_muni_unknown.sum())
        selected_admission_year_counts.update(selected_years.value_counts(dropna=False).to_dict())

        diag4 = dx_clean.loc[select_mask].str.slice(0, 4)
        selected_diag4_counts.update(diag4.value_counts(dropna=False).to_dict())
        selected["source_year"] = source_year
        selected["hrid"] = mapped.loc[select_mask].fillna("unknown")
        selected.loc[selected["hrid"].astype(str).str.strip().eq(""), "hrid"] = "unknown"
        selected["admission_date"] = selected_dates.dt.strftime("%Y-%m-%d").fillna("unknown")
        selected["diag_princ4"] = diag4
        selected["sexo_raw"] = selected["SEXO"].astype(str).str.strip().replace("", "unknown")
        selected["age_group"] = _age_group(selected["IDADE"], selected["COD_IDADE"])
        selected["car_int_raw"] = selected["CAR_INT"].astype(str).str.strip().replace("", "unknown")
        selected["car_int_group"] = _car_group(selected["CAR_INT"])

        grouped = (
            selected.groupby(AGG_DIMS, dropna=False, sort=False)
            .size()
            .rename("admission_records")
            .reset_index()
        )
        aggregate_parts.append(grouped)

        # Temporary HMAC sidecar supports a distinct-N_AIH sensitivity and repeat audit.
        # It contains no raw admission ID and is deleted after the final audit is written.
        id_values = selected["N_AIH"].astype(str).str.strip()
        selected_missing_naih += int(id_values.eq("").sum())
        nonmissing_id = id_values.ne("")
        row_hmacs: list[str] = []
        id_hmacs: list[str] = []
        selected_raw = selected.loc[nonmissing_id]
        selected_ids = id_values.loc[nonmissing_id]
        for value in selected_ids.tolist():
            id_hmacs.append(_hmac(value))
        # Row hashes cover every original CSV field, so exact-repeat detection is strict.
        for values in selected_raw[columns].itertuples(index=False, name=None):
            row_hmacs.append(_hmac("\x1f".join(str(v) for v in values)))
        if len(selected_raw):
            side = selected_raw[AGG_DIMS].copy()
            side["naih_hmac"] = id_hmacs
            side["row_hmac"] = row_hmacs
            side["source_file"] = source.name
            side = (
                side.groupby(
                    ["naih_hmac", "row_hmac", "source_file", *AGG_DIMS],
                    dropna=False,
                    sort=False,
                )
                .size()
                .rename("record_rows")
                .reset_index()
            )
            id_parts.append(side)

    if columns is None:
        raise ValueError(f"No CSV rows/headers found in {source}")

    if aggregate_parts:
        aggregate = pd.concat(aggregate_parts, ignore_index=True)
        aggregate = (
            aggregate.groupby(AGG_DIMS, dropna=False, sort=False)["admission_records"]
            .sum()
            .reset_index()
        )
    else:
        aggregate = pd.DataFrame(columns=[*AGG_DIMS, "admission_records"])
    _csv_gz_write(aggregate, agg_path)

    if id_parts:
        sidecar = pd.concat(id_parts, ignore_index=True)
        sidecar = (
            sidecar.groupby(
                ["naih_hmac", "row_hmac", "source_file", *AGG_DIMS],
                dropna=False,
                sort=False,
            )["record_rows"]
            .sum()
            .reset_index()
        )
    else:
        sidecar = pd.DataFrame(columns=["naih_hmac", "row_hmac", "source_file", *AGG_DIMS, "record_rows"])
    _csv_gz_write(sidecar, id_path)
    try:
        id_path.chmod(0o600)
    except OSError:
        pass

    audit = {
        "source_year": source_year,
        "source_file": source.name,
        "source_size_bytes": size,
        "source_mtime_ns": mtime_ns,
        "source_rows": total_rows,
        "selected_target_diagnosis_rows": selected_rows,
        "selected_missing_n_aih_rows": selected_missing_naih,
        "unknown_municipality_rows": unknown_muni,
        "selected_unknown_municipality_rows": selected_unknown_muni,
        "admission_year_counts": dict(sorted(admission_year_counts.items())),
        "selected_admission_year_counts": dict(sorted(selected_admission_year_counts.items())),
        "date_format_counts": dict(sorted(date_format_counts.items())),
        "cod_idade_raw_counts": dict(sorted(age_code_counts.items())),
        "cod_idade_idade_raw_pair_counts": dict(sorted(age_value_counts.items())),
        "sexo_raw_counts": dict(sorted(sex_counts.items())),
        "car_int_raw_counts": dict(sorted(car_counts.items())),
        "principal_dx3_counts": dict(sorted(principal_base_counts.items())),
        "selected_diag_princ4_counts": dict(sorted(selected_diag4_counts.items())),
        "elapsed_seconds": round(time.time() - started, 3),
    }
    _json_write(audit_path, audit)
    return {
        "status": "processed",
        "source_year": source_year,
        "source_file": source.name,
        "source_rows": total_rows,
        "selected_rows": selected_rows,
        "elapsed_seconds": audit["elapsed_seconds"],
    }


def _sum_dict_counts(audits: list[dict[str, Any]], key: str) -> dict[str, int]:
    result: Counter[str] = Counter()
    for audit in audits:
        result.update({str(k): int(v) for k, v in audit.get(key, {}).items()})
    return dict(sorted(result.items()))


def _sidecar_paths(output_root: Path) -> list[Path]:
    return sorted((output_root / "intermediate" / "private_id_hash").rglob("*.csv.gz"))


def _build_final_outputs(output_root: Path, source_files: list[Path], crosswalk: dict[str, str]) -> dict[str, Any]:
    per_file_audits: list[dict[str, Any]] = []
    aggregate_paths: list[Path] = []
    for source in source_files:
        agg_path, audit_path, _ = _paths_for_source(source, output_root)
        if not agg_path.exists() or not audit_path.exists():
            raise RuntimeError(f"Missing per-file result for {source.name}")
        per_file_audits.append(json.loads(audit_path.read_text(encoding="utf-8")))
        aggregate_paths.append(agg_path)

    all_aggregates = [pd.read_csv(p, dtype={d: str for d in AGG_DIMS}, keep_default_na=False) for p in aggregate_paths]
    aggregate = pd.concat(all_aggregates, ignore_index=True)
    if aggregate.empty:
        aggregate = pd.DataFrame(columns=[*AGG_DIMS, "admission_records"])
    aggregate["admission_records"] = pd.to_numeric(aggregate["admission_records"], errors="coerce").fillna(0).astype("int64")
    aggregate = (
        aggregate.groupby(AGG_DIMS, dropna=False, sort=False)["admission_records"]
        .sum()
        .reset_index()
    )

    sidecar_paths = _sidecar_paths(output_root)
    sidecars = []
    for path in sidecar_paths:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if not frame.empty:
            frame["record_rows"] = pd.to_numeric(frame["record_rows"], errors="coerce").fillna(0).astype("int64")
            sidecars.append(frame)
    id_audit: dict[str, Any]
    if sidecars:
        ids = pd.concat(sidecars, ignore_index=True)
        id_by_cell = ids.drop_duplicates(subset=["naih_hmac", *AGG_DIMS]).groupby(AGG_DIMS, dropna=False).size().rename("unique_n_aih").reset_index()
        aggregate = aggregate.merge(id_by_cell, on=AGG_DIMS, how="left")
        aggregate["unique_n_aih"] = aggregate["unique_n_aih"].fillna(0).astype("int64")

        all_id_rows = ids.groupby("naih_hmac", sort=False).agg(
            record_rows=("record_rows", "sum"),
            distinct_full_rows=("row_hmac", "nunique"),
            distinct_admission_dates=("admission_date", "nunique"),
            distinct_diag4=("diag_princ4", "nunique"),
            distinct_source_files=("source_file", "nunique"),
            distinct_source_years=("source_year", "nunique"),
        )
        repeated = all_id_rows[all_id_rows["record_rows"].gt(1)]
        exact_rows = ids.groupby("row_hmac", sort=False)["record_rows"].sum()
        exact_duplicate_row_groups = exact_rows[exact_rows.gt(1)]
        exact_duplicate_rows = int((exact_rows[exact_rows.gt(1)] - 1).sum())
        duplicate_id_patterns = {
            "unique_nonmissing_n_aih": int(len(all_id_rows)),
            "n_aih_with_multiple_selected_rows": int(len(repeated)),
            "extra_selected_rows_beyond_one_per_n_aih": int((repeated["record_rows"] - 1).sum()),
            "repeated_n_aih_with_multiple_admission_dates": int(repeated["distinct_admission_dates"].gt(1).sum()),
            "repeated_n_aih_with_multiple_diag4": int(repeated["distinct_diag4"].gt(1).sum()),
            "repeated_n_aih_across_source_files": int(repeated["distinct_source_files"].gt(1).sum()),
            "repeated_n_aih_across_source_years": int(repeated["distinct_source_years"].gt(1).sum()),
            "repeated_n_aih_same_date_same_diag4_only": int(
                (repeated["distinct_admission_dates"].eq(1) & repeated["distinct_diag4"].eq(1)).sum()
            ),
            "repeated_n_aih_with_multiple_full_row_signatures": int(repeated["distinct_full_rows"].gt(1).sum()),
            "exact_duplicate_full_row_signature_groups": int(len(exact_duplicate_row_groups)),
            "extra_selected_rows_that_are_exact_full_row_duplicates": exact_duplicate_rows,
        }
        id_audit = {
            **duplicate_id_patterns,
            "n_aih_type_column_available": False,
            "n_aih_type_column_name": None,
            "interpretation": (
                "The source CSV header has no IDENT field. Repeated-ID patterns are reported, "
                "but the available data cannot identify type-5 continuation versus type-1 initial "
                "AIH or establish whether a repeated authorization is a resubmission."
            ),
            "unique_n_aih_sensitivity": (
                "The aggregate files include unique_n_aih by cell as a sensitivity measure. "
                "Primary admission_records retain every selected SIH-RD row."
            ),
        }
    else:
        aggregate["unique_n_aih"] = 0
        id_audit = {
            "unique_nonmissing_n_aih": 0,
            "n_aih_type_column_available": False,
            "interpretation": "No nonmissing selected N_AIH sidecar records were available.",
        }

    aggregate = aggregate[[*AGG_DIMS, "admission_records", "unique_n_aih"]]
    aggregate.sort_values(AGG_DIMS, kind="mergesort", inplace=True, ignore_index=True)
    results = output_root / "results"
    results.mkdir(parents=True, exist_ok=True)
    all_path = results / "brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz"
    final_path = results / "brazil_clinical_aggregate_adm2014_2020_by_source_year.csv.gz"
    _csv_gz_write(aggregate, all_path)

    # Keep only parsed admission dates in the requested 2014–2020 analytic window.
    final_aggregate = aggregate.copy()
    final_aggregate["admission_year"] = pd.to_numeric(final_aggregate["admission_date"].str.slice(0, 4), errors="coerce")
    final_aggregate = final_aggregate[final_aggregate["admission_year"].between(2014, 2020, inclusive="both")]
    final_aggregate = final_aggregate.drop(columns=["admission_year"])
    _csv_gz_write(final_aggregate, final_path)

    total_rows_by_source_year: Counter[str] = Counter()
    selected_by_source_year: Counter[str] = Counter()
    selected_by_admission_year: Counter[str] = Counter()
    selected_source_admission_year: dict[str, dict[str, int]] = {}
    for audit in per_file_audits:
        sy = str(audit["source_year"])
        total_rows_by_source_year[sy] += int(audit["source_rows"])
        selected_by_source_year[sy] += int(audit["selected_target_diagnosis_rows"])
        selected_by_admission_year.update(
            {str(k): int(v) for k, v in audit.get("selected_admission_year_counts", {}).items()}
        )
        row = selected_source_admission_year.setdefault(sy, {})
        for ay, count in audit.get("selected_admission_year_counts", {}).items():
            row[str(ay)] = row.get(str(ay), 0) + int(count)

    source_admission_crosswalk = {
        sy: dict(sorted(values.items())) for sy, values in sorted(selected_source_admission_year.items())
    }
    source_file_counts = Counter(str(audit["source_year"]) for audit in per_file_audits)
    audit_json = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_root": str(SOURCE_ROOT),
        "crosswalk_path": str(CROSSWALK),
        "crosswalk_rows": len(crosswalk),
        "source_files_processed": len(per_file_audits),
        "source_file_counts_by_year": dict(sorted(source_file_counts.items())),
        "source_rows_by_source_year": dict(sorted(total_rows_by_source_year.items())),
        "target_diagnosis_rows_by_source_year": dict(sorted(selected_by_source_year.items())),
        "target_diagnosis_rows_by_admission_year": dict(sorted(selected_by_admission_year.items())),
        "target_diagnosis_source_year_by_admission_year": source_admission_crosswalk,
        "target_diag3": sorted(TARGET_DX3),
        "target_diag4_counts": _sum_dict_counts(per_file_audits, "selected_diag_princ4_counts"),
        "diagnosis_definition": {
            "source_field": "DIAG_PRINC",
            "normalization": (
                "Uppercase, remove punctuation/non-alphanumeric characters, then retain up to the first "
                "four characters. Source values that contain only the three-character category remain "
                "three characters; no subcode is imputed."
            ),
            "selection": "Normalized first three characters in K80, K81, N20, N21, N22, N23.",
        },
        "region_definition": {
            "output_field": "hrid",
            "crosswalk_field": "health_region_2023",
            "unmapped_policy": "Retain selected rows with hrid=unknown; do not silently drop them.",
        },
        "aggregate_measures": {
            "admission_records": "Primary count of every selected SIH-RD source row.",
            "unique_n_aih": (
                "Number of distinct nonmissing N_AIH values within each source-year and full subgroup cell. "
                "An N_AIH appearing in more than one cell contributes once to each affected cell."
            ),
            "sexo_raw": "Trimmed source SEXO value, with blank mapped to unknown; no sex recoding applied.",
        },
        "principal_dx3_counts_all_source_rows": _sum_dict_counts(per_file_audits, "principal_dx3_counts"),
        "cod_idade_raw_counts_all_source_rows": _sum_dict_counts(per_file_audits, "cod_idade_raw_counts"),
        "cod_idade_idade_raw_pair_counts_all_source_rows": _sum_dict_counts(
            per_file_audits, "cod_idade_idade_raw_pair_counts"
        ),
        "sexo_raw_counts_all_source_rows": _sum_dict_counts(per_file_audits, "sexo_raw_counts"),
        "car_int_raw_counts_all_source_rows": _sum_dict_counts(per_file_audits, "car_int_raw_counts"),
        "date_format_counts_all_source_rows": _sum_dict_counts(per_file_audits, "date_format_counts"),
        "admission_year_counts_all_source_rows_by_source_year": {
            str(year): _sum_dict_counts(
                [a for a in per_file_audits if str(a["source_year"]) == str(year)],
                "admission_year_counts",
            )
            for year in sorted(source_file_counts)
        },
        "unknown_municipality_rows_all_source_rows": int(sum(a["unknown_municipality_rows"] for a in per_file_audits)),
        "unknown_municipality_target_diagnosis_rows": int(
            sum(a["selected_unknown_municipality_rows"] for a in per_file_audits)
        ),
        "selected_target_diagnosis_rows_missing_n_aih": int(
            sum(a["selected_missing_n_aih_rows"] for a in per_file_audits)
        ),
        "aggregate_all_years_path": str(all_path),
        "aggregate_admission_2014_2020_path": str(final_path),
        "all_years_aggregate_cells": int(len(aggregate)),
        "admission_2014_2020_aggregate_cells": int(len(final_aggregate)),
        "admission_2014_2020_records": int(final_aggregate["admission_records"].sum()),
        "admission_2014_2020_unique_n_aih_cell_counts_sum": int(final_aggregate["unique_n_aih"].sum()),
        "date_window_definition": "DT_INTER admission date parsed from SIH-RD; 2014-01-01 through 2020-12-31 inclusive.",
        "sourceyear_interpretation": (
            "All available source-year folders 2014–2021 were processed. Source year is retained "
            "in aggregate dimensions and audit tables; no completeness is inferred from 2021 files."
        ),
        "age_definition": {
            "source_fields": ["IDADE", "COD_IDADE"],
            "dictionary_units": {"0": "ignored/unknown", "2": "days", "3": "months", "4": "years", "5": "100+ years"},
            "bins": ["<18", "18-44", "45-64", "65plus", "unknown"],
            "code_2_valid_raw_range": "0-30 days",
            "code_3_valid_raw_range": "0-11 months",
            "code_4_valid_raw_range": "0-99 years",
        },
        "car_int_definition": {
            "car_int_raw": "Trimmed source value retained as a separate output dimension; blank mapped to unknown.",
            "car_int_group": "01 and 02 separate; other nonmissing codes pooled as other; 0/blank/null pooled as unknown.",
        },
        "duplicate_n_aih_audit": id_audit,
        "per_file_elapsed_seconds_sum": round(sum(float(a.get("elapsed_seconds", 0)) for a in per_file_audits), 3),
    }
    _json_write(results / "brazil_clinical_extraction_audit.json", audit_json)
    _json_write(results / "brazil_clinical_duplicate_naih_audit.json", id_audit)

    # Remove only this run's private HMAC sidecars and key, after the audit and
    # unique-N_AIH cell counts have been materialized into aggregate outputs.
    private_root = output_root / "intermediate" / "private_id_hash"
    if private_root.exists():
        for private_file in private_root.rglob("*.csv.gz"):
            private_file.unlink()
        for folder in sorted(private_root.rglob("*"), reverse=True):
            if folder.is_dir():
                try:
                    folder.rmdir()
                except OSError:
                    pass
        try:
            private_root.rmdir()
        except OSError:
            pass
    key_path = output_root / "intermediate" / "private" / "id_hmac_key.bin"
    if key_path.exists():
        key_path.unlink()
    try:
        key_path.parent.rmdir()
    except OSError:
        pass
    return audit_json


def _get_hash_key(output_root: Path) -> bytes:
    private_dir = output_root / "intermediate" / "private"
    private_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        private_dir.chmod(0o700)
    except OSError:
        pass
    key_path = private_dir / "id_hmac_key.bin"
    if key_path.exists():
        key = key_path.read_bytes()
    else:
        key = os.urandom(32)
        fd = os.open(key_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "wb") as handle:
            handle.write(key)
    if len(key) != 32:
        raise ValueError(f"Invalid private HMAC key length at {key_path}")
    return key


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--crosswalk", type=Path, default=CROSSWALK)
    parser.add_argument("--output-root", type=Path, default=REMOTE_ROOT)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunksize", type=int, default=100_000)
    parser.add_argument("--years", nargs="+", type=int, default=list(YEARS))
    args = parser.parse_args()

    if args.workers < 1 or args.workers > 8:
        parser.error("--workers must be between 1 and 8")
    if args.chunksize < 10_000:
        parser.error("--chunksize must be at least 10000")
    if not args.crosswalk.is_file():
        raise FileNotFoundError(args.crosswalk)
    crosswalk = _load_crosswalk(args.crosswalk)
    source_files: list[Path] = []
    for year in args.years:
        year_dir = args.source_root / str(year)
        source_files.extend(sorted(year_dir.glob("*/*.csv.gz")))
    if not source_files:
        raise FileNotFoundError(f"No source files found under {args.source_root} for {args.years}")

    output_root = args.output_root
    for folder in (output_root, output_root / "intermediate", output_root / "results", output_root / "scripts"):
        folder.mkdir(parents=True, exist_ok=True)
    key = _get_hash_key(output_root)
    private_hash_root = output_root / "intermediate" / "private_id_hash"
    private_hash_root.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        private_hash_root.chmod(0o700)
    except OSError:
        pass
    print(
        json.dumps(
            {
                "event": "starting",
                "source_files": len(source_files),
                "source_years": sorted(set(int(p.parent.parent.name) for p in source_files)),
                "workers": args.workers,
                "chunksize": args.chunksize,
                "crosswalk_municipality_count": len(crosswalk),
                "output_root": str(output_root),
            }
        ),
        flush=True,
    )

    completed = 0
    selected_rows = 0
    source_rows = 0
    failures: list[dict[str, str]] = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_worker_init,
        initargs=(crosswalk, key),
    ) as pool:
        futures = {
            pool.submit(_process_file, source, output_root, args.chunksize): source for source in source_files
        }
        for future in concurrent.futures.as_completed(futures):
            source = futures[future]
            try:
                result = future.result()
                completed += 1
                selected_rows += int(result.get("selected_rows", 0))
                source_rows += int(result.get("source_rows", 0))
                if completed % 100 == 0 or completed == len(source_files):
                    print(
                        json.dumps(
                            {
                                "event": "progress",
                                "completed_files": completed,
                                "total_files": len(source_files),
                                "source_rows_processed": source_rows,
                                "selected_rows_processed_this_run": selected_rows,
                                "last_status": result["status"],
                                "last_source_year": result["source_year"],
                            }
                        ),
                        flush=True,
                    )
            except Exception as exc:  # keep processing independent files, then fail closed
                failures.append({"source_file": source.name, "error": f"{type(exc).__name__}: {exc}"})
                print(json.dumps({"event": "file_error", **failures[-1]}), flush=True)

    if failures:
        _json_write(output_root / "results" / "brazil_clinical_extraction_failures.json", failures)
        raise RuntimeError(f"{len(failures)} source files failed; see aggregate failure audit")

    audit = _build_final_outputs(output_root, source_files, crosswalk)
    print(
        json.dumps(
            {
                "event": "complete",
                "source_files": len(source_files),
                "source_rows": sum(int(v) for v in audit["source_rows_by_source_year"].values()),
                "target_diagnosis_rows_by_source_year": audit["target_diagnosis_rows_by_source_year"],
                "admission_2014_2020_records": audit["admission_2014_2020_records"],
                "aggregate_cells_2014_2020": audit["admission_2014_2020_aggregate_cells"],
                "duplicate_n_aih_audit": audit["duplicate_n_aih_audit"],
            }
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
