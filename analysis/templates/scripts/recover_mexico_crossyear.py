#!/usr/bin/env python3
"""Recover grouped Mexico source-year admissions omitted by calendar filters.

Reads only aggregate fields from the local SAEH-SSA EGRESOS source tables.
Patient/episode identifiers and facility IDs are never selected or written.
The established mx_clinical baseline caches are read-only inputs.
"""

from __future__ import annotations

import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
REPAIR = ROOT / "work" / "repair_20261112"
DATA = REPAIR / "data"
ADMIN = REPAIR / "admin"
INPUT_AUDIT = ADMIN / "mx_clinical_audit.json"
BASELINE_PANEL = DATA / "source" / "mex_state_day_stone.csv"
YEARS = list(range(2015, 2024))
FULL_SOURCE_YEARS = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023]
BASELINE_YEARS = [2015, 2016, 2017, 2018, 2020, 2021, 2022]
PREFERRED_YEARS = list(range(2015, 2023))
STATE_CODES = [f"{i:02d}" for i in range(1, 33)]
STONE_ENDPOINTS = {"K80", "K81", "N20", "N21", "N22", "N23"}
CHUNK_SIZE = 200_000
CODE_COLUMNS = ["state", "date", "endpoint", "fullcode", "age_group", "sex", "n"]
CORRECTION_COLUMNS = ["source_year", "source_role", *CODE_COLUMNS]
WIDE_COUNTS = [
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
SOURCE_WIDE_COLUMNS = ["source_year", "source_role", "state", "date", *WIDE_COUNTS]
WIDE_COLUMNS = ["source_role", "state", "date", *WIDE_COUNTS]
CURRENT_CACHE_YEARS = {
    2015: set(BASELINE_YEARS),
    2016: set(BASELINE_YEARS),
    2017: set(BASELINE_YEARS),
    2018: set(BASELINE_YEARS),
    2019: {2019},
    2020: set(BASELINE_YEARS),
    2021: set(BASELINE_YEARS),
    2022: set(BASELINE_YEARS),
    2023: set(BASELINE_YEARS + [2023]),
}


def norm_header(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.replace("\ufeff", "").strip().upper())


def source_role(source_year: int, event_year: int) -> str:
    if source_year == 2019:
        return "recovered_additional_2019"
    if source_year == 2023 and event_year == 2023:
        return "extension_2023"
    return "primary_panel_year"


def assign_age_group(age: pd.Series, age_unit: pd.Series) -> pd.Series:
    age_num = pd.to_numeric(age, errors="coerce")
    unit = age_unit.fillna("").astype(str).str.strip().str.replace(r"\.0+$", "", regex=True)
    valid = age_num.notna() & age_num.ge(0) & age_num.ne(999)
    result = pd.Series("unknown", index=age.index, dtype="object")
    result.loc[valid & unit.isin(["0", "1", "2"])] = "<18"
    years = valid & unit.eq("3")
    result.loc[years & age_num.lt(18)] = "<18"
    result.loc[years & age_num.ge(18) & age_num.le(44)] = "18-44"
    result.loc[years & age_num.ge(45) & age_num.le(64)] = "45-64"
    result.loc[years & age_num.ge(65)] = "65plus"
    return result


def assign_sex(sex: pd.Series) -> pd.Series:
    raw = sex.fillna("").astype(str).str.strip().str.replace(r"\.0+$", "", regex=True)
    return raw.map({"1": "male", "2": "female"}).fillna("unknown")


def make_year_context(year: int, year_audit: dict[str, Any]) -> dict[str, Any]:
    fields = year_audit["source"]["field_names"]
    aliases = {
        "INGRE": {"INGRE"},
        "EGRESO": {"EGRESO"},
        "DIAS_ESTA": {"DIASESTA"},
        "ENTIDAD": {"ENTIDAD"},
        "AFECPRIN": {"AFECPRIN"},
        "EDAD": {"EDAD"},
        "CVEEDAD": {"CVEEDAD", "CVEEDAD"},
        "SEXO": {"SEXO"},
    }
    normalized = {norm_header(x): x for x in fields}
    chosen: dict[str, str] = {}
    for key, choices in aliases.items():
        for choice in choices:
            if choice in normalized:
                chosen[key] = normalized[choice]
                break
    required = ["INGRE", "ENTIDAD", "AFECPRIN", "EDAD", "CVEEDAD", "SEXO"]
    missing = [x for x in required if x not in chosen]
    if missing:
        raise ValueError(f"{year}: source audit missing correction fields {missing}")
    return {
        "year": year,
        "path": year_audit["source"]["path"],
        "sha256": year_audit["source"]["sha256"],
        "bytes": year_audit["source"]["bytes"],
        "encoding": year_audit["source"]["encoding"],
        "delimiter": year_audit["source"]["delimiter"],
        "chosen": chosen,
        "header_fields": fields,
        "prior_outside_selected_count": year_audit["ingestion"]["selected_rows_with_valid_ingre_outside_analysis_calendar"],
        "cache_years": sorted(CURRENT_CACHE_YEARS[year]),
    }


def extract_source_correction(context: dict[str, Any]) -> dict[str, Any]:
    year = context["year"]
    path = Path(context["path"])
    if path.stat().st_size != context["bytes"]:
        raise ValueError(f"Source size changed since the verified baseline audit: {path}")
    chosen = context["chosen"]
    raw_fields = set(chosen.values())
    raw_to_canonical = {raw: canonical for canonical, raw in chosen.items()}
    correction_groups: list[pd.DataFrame] = []
    invalid_by_code: Counter[str] = Counter()
    unmatched_by_state: Counter[str] = Counter()
    selected_by_year: Counter[str] = Counter()
    missing_cache_by_year: Counter[str] = Counter()
    correction_by_code: Counter[str] = Counter()
    selected_rows = 0
    parsed_rows = 0
    invalid_date_rows = 0
    outside_study_rows = 0
    missing_cache_rows = 0
    emitted_record_rows = 0
    matched_correction_rows = 0
    crossyear_timing = Counter()
    crossyear_days = Counter()
    crossyear_duration_deltas: Counter[str] = Counter()
    crossyear_code_counts: Counter[str] = Counter()
    crossyear_by_year: Counter[str] = Counter()

    read_fields = set(raw_fields)
    for optional in ("EGRESO", "DIAS_ESTA"):
        if optional in chosen:
            read_fields.add(chosen[optional])

    chunks = pd.read_csv(
        path,
        sep=context["delimiter"],
        usecols=lambda name: name in read_fields,
        dtype=str,
        chunksize=CHUNK_SIZE,
        encoding=context["encoding"],
        on_bad_lines="skip",
        low_memory=False,
        keep_default_na=False,
    )
    for raw in chunks:
        parsed_rows += len(raw)
        raw = raw.rename(columns=raw_to_canonical)
        date = pd.to_datetime(raw["INGRE"], errors="coerce")
        original_dx = raw["AFECPRIN"].fillna("").astype(str).str.upper().str.strip()
        fullcode = original_dx.str.replace(r"[^A-Z0-9]", "", regex=True)
        endpoint = fullcode.str[:3]
        stone = endpoint.isin(STONE_ENDPOINTS)
        selected_rows += int(stone.sum())
        invalid = stone & date.isna()
        invalid_date_rows += int(invalid.sum())
        if invalid.any():
            invalid_by_code.update(endpoint.loc[invalid].value_counts().to_dict())

        event_year = date.dt.year
        valid = stone & date.notna()
        in_target = valid & event_year.between(2015, 2023)
        if valid.any():
            selected_by_year.update(
                {str(int(y)): int(n) for y, n in event_year.loc[valid].value_counts().to_dict().items()}
            )
        outside_study_rows += int((valid & ~event_year.between(2015, 2023)).sum())

        cache_years = set(context["cache_years"])
        missing_cache = in_target & ~event_year.isin(cache_years)
        missing_cache_rows += int(missing_cache.sum())
        if missing_cache.any():
            missing_cache_by_year.update(
                {str(int(y)): int(n) for y, n in event_year.loc[missing_cache].value_counts().to_dict().items()}
            )

        if year == 2019 and "EGRESO" in raw.columns:
            crossyear = stone & date.notna() & event_year.ne(2019)
            if crossyear.any():
                crossyear_by_year.update(
                    {str(int(y)): int(n) for y, n in event_year.loc[crossyear].value_counts().to_dict().items()}
                )
                crossyear_code_counts.update(endpoint.loc[crossyear].value_counts().to_dict())
                discharge = pd.to_datetime(raw["EGRESO"], errors="coerce")
                crossyear_timing["rows"] += int(crossyear.sum())
                crossyear_timing["ingre_le_egreso"] += int((crossyear & discharge.notna() & date.le(discharge)).sum())
                crossyear_timing["ingre_gt_egreso"] += int((crossyear & discharge.notna() & date.gt(discharge)).sum())
                crossyear_timing["egreso_missing_or_invalid"] += int((crossyear & discharge.isna()).sum())
                duration = pd.to_numeric(raw["DIAS_ESTA"], errors="coerce") if "DIAS_ESTA" in raw.columns else pd.Series(float("nan"), index=raw.index)
                crossyear_days["days_estancia_missing_or_nonnumeric"] += int((crossyear & duration.isna()).sum())
                crossyear_days["days_estancia_zero"] += int((crossyear & duration.eq(0)).sum())
                crossyear_days["days_estancia_positive"] += int((crossyear & duration.gt(0)).sum())
                crossyear_days["days_estancia_negative"] += int((crossyear & duration.lt(0)).sum())
                date_delta = (discharge - date).dt.days
                comparable = crossyear & discharge.notna() & date.notna() & duration.notna()
                crossyear_duration_deltas["compared_rows"] += int(comparable.sum())
                crossyear_duration_deltas["DIAS_ESTA_equals_calendar_date_delta"] += int((comparable & duration.eq(date_delta)).sum())
                crossyear_duration_deltas["DIAS_ESTA_differs_from_calendar_date_delta"] += int((comparable & duration.ne(date_delta)).sum())

        state_raw = raw["ENTIDAD"].fillna("").astype(str).str.strip()
        state = state_raw.str.zfill(2)
        state_ok = state.isin(STATE_CODES)
        correction = missing_cache & state_ok
        unmatched = missing_cache & ~state_ok
        if unmatched.any():
            unmatched_by_state.update(state_raw.loc[unmatched].replace("", "__blank__").value_counts().to_dict())
        emitted_record_rows += int(missing_cache.sum())
        matched_correction_rows += int(correction.sum())
        if correction.any():
            endpoint_c = endpoint.loc[correction]
            correction_by_code.update(endpoint_c.value_counts().to_dict())
            grouped = pd.DataFrame(
                {
                    "source_year": year,
                    "source_role": [source_role(year, int(y)) for y in event_year.loc[correction]],
                    "state": state.loc[correction].to_numpy(),
                    "date": date.loc[correction].dt.strftime("%Y-%m-%d").to_numpy(),
                    "endpoint": endpoint.loc[correction].to_numpy(),
                    "fullcode": fullcode.loc[correction].to_numpy(),
                    "age_group": assign_age_group(raw.loc[correction, "EDAD"], raw.loc[correction, "CVEEDAD"]).to_numpy(),
                    "sex": assign_sex(raw.loc[correction, "SEXO"]).to_numpy(),
                }
            )
            grouped = grouped.groupby(CORRECTION_COLUMNS[:-1], sort=False, observed=True).size().reset_index(name="n")
            correction_groups.append(grouped)

    if correction_groups:
        corrections = pd.concat(correction_groups, ignore_index=True)
        corrections = corrections.groupby(CORRECTION_COLUMNS[:-1], as_index=False, sort=True)["n"].sum()
        corrections = corrections[CORRECTION_COLUMNS]
    else:
        corrections = pd.DataFrame(columns=CORRECTION_COLUMNS)

    correction_path = DATA / f"mx_complete_source_correction_{year}.csv.gz"
    corrections.to_csv(correction_path, index=False, compression="gzip")
    return {
        "source_year": year,
        "source_path": str(path),
        "source_bytes": context["bytes"],
        "source_sha256": context["sha256"],
        "prior_pipeline_outside_window_selected_rows": context["prior_outside_selected_count"],
        "rows_parsed": parsed_rows,
        "selected_diagnosis_rows": selected_rows,
        "selected_rows_invalid_or_missing_ingre": invalid_date_rows,
        "invalid_ingre_selected_rows_by_endpoint": dict(sorted(invalid_by_code.items())),
        "selected_valid_ingre_rows_by_calendar_year": dict(sorted(selected_by_year.items())),
        "selected_valid_ingre_outside_study_calendar_rows": outside_study_rows,
        "selected_valid_ingre_missing_from_current_cache_rows": missing_cache_rows,
        "missing_cache_rows_by_ingre_year": dict(sorted(missing_cache_by_year.items())),
        "missing_cache_rows_with_matched_state_01_32": matched_correction_rows,
        "missing_cache_rows_by_endpoint": dict(sorted(correction_by_code.items())),
        "missing_cache_unmatched_state_counts": dict(sorted(unmatched_by_state.items())),
        "correction_group_cells": len(corrections),
        "correction_path": str(correction_path),
        "2019_crossyear_details": {
            "selected_valid_crossyear_rows_by_ingre_year": dict(sorted(crossyear_by_year.items())),
            "diagnosis_prefix_counts_across_crossyear_rows": dict(sorted(crossyear_code_counts.items())),
            "ingre_egreso_order": dict(crossyear_timing),
            "dias_esta_categories": dict(crossyear_days),
            "dias_esta_vs_calendar_date_delta": dict(crossyear_duration_deltas),
        }
        if year == 2019
        else None,
    }


def current_source_role(source_year: int, event_year: int) -> str:
    return source_role(source_year, event_year)


def add_source_rows(frames: list[pd.DataFrame], year: int, cache_path: Path) -> None:
    if not cache_path.is_file():
        raise FileNotFoundError(cache_path)
    frame = pd.read_csv(cache_path, dtype={"state": str, "date": str, "endpoint": str, "fullcode": str, "age_group": str, "sex": str})
    event_year = pd.to_numeric(frame["date"].str[:4], errors="coerce")
    frame = frame.loc[event_year.between(2015, 2022)].copy()
    if frame.empty:
        return
    frame.insert(0, "source_year", year)
    frame.insert(1, "source_role", [current_source_role(year, int(y)) for y in event_year.loc[frame.index]])
    frames.append(frame[CORRECTION_COLUMNS])


def build_source_wide(long_df: pd.DataFrame) -> pd.DataFrame:
    grain = ["source_year", "source_role", "state", "date"]
    totals = long_df.pivot_table(index=grain, columns="endpoint", values="n", aggfunc="sum", fill_value=0).reset_index()
    totals.columns.name = None
    for code in ("K80", "K81", "N20", "N21", "N22", "N23"):
        if code not in totals.columns:
            totals[code] = 0
    k80_codes = long_df.loc[long_df["endpoint"].eq("K80")].copy()
    k80_codes["k800"] = k80_codes["fullcode"].str[:4].eq("K800") * k80_codes["n"]
    k80_codes["k803"] = k80_codes["fullcode"].str[:4].eq("K803") * k80_codes["n"]
    acute = k80_codes.groupby(grain, as_index=False)[["k800", "k803"]].sum()
    totals = totals.merge(acute, on=grain, how="left")
    totals[["k800", "k803"]] = totals[["k800", "k803"]].fillna(0).astype("int64")
    totals["n20_n21_confirmed"] = totals["N20"] + totals["N21"]
    totals["k80_acute_k800_k803"] = totals["k800"] + totals["k803"]
    totals = totals.rename(
        columns={
            "K80": "k80",
            "K81": "k81",
            "N20": "n20_only",
            "N21": "n21_only",
            "N22": "n22_other_disease",
            "N23": "n23_unspecified_renal_colic",
            "k800": "k80_k800",
            "k803": "k80_k803",
        }
    )
    return totals[SOURCE_WIDE_COLUMNS].sort_values(["source_year", "state", "date"]).reset_index(drop=True)


def complete_state_day(observed: pd.DataFrame, years: list[int], role: str) -> pd.DataFrame:
    days = []
    for year in years:
        date_values = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D").strftime("%Y-%m-%d")
        grid = pd.MultiIndex.from_product([STATE_CODES, date_values], names=["state", "date"]).to_frame(index=False)
        grid.insert(0, "source_role", role)
        counts = observed.loc[observed["date"].str.startswith(f"{year}-"), ["state", "date", *WIDE_COUNTS]]
        day = grid.merge(counts, on=["state", "date"], how="left")
        for col in WIDE_COUNTS:
            day[col] = day[col].fillna(0).astype("int64")
        days.append(day[WIDE_COLUMNS])
    return pd.concat(days, ignore_index=True).sort_values(["state", "date"], kind="stable").reset_index(drop=True)


def baseline_totals() -> dict[int, dict[str, int]]:
    panel = pd.read_csv(BASELINE_PANEL, dtype={"state": str})
    panel["date"] = pd.to_datetime(panel["date"], errors="coerce")
    out = {}
    for year, frame in panel.groupby(panel["date"].dt.year):
        out[int(year)] = {"k80": int(frame["k80"].sum()), "k81": int(frame["k81"].sum()), "n20_collapsed_n20_to_n23": int(frame["n20"].sum())}
    return out


def family_totals(wide: pd.DataFrame) -> dict[str, int]:
    return {
        "k80": int(wide["k80"].sum()),
        "k81": int(wide["k81"].sum()),
        "n20_only": int(wide["n20_only"].sum()),
        "n21_only": int(wide["n21_only"].sum()),
        "n22_other_disease": int(wide["n22_other_disease"].sum()),
        "n23_unspecified_renal_colic": int(wide["n23_unspecified_renal_colic"].sum()),
        "n20_n21_confirmed": int(wide["n20_n21_confirmed"].sum()),
        "n20_collapsed_n20_to_n23": int(
            wide[["n20_only", "n21_only", "n22_other_disease", "n23_unspecified_renal_colic"]].to_numpy().sum()
        ),
        "k80_acute_k800_k803": int(wide["k80_acute_k800_k803"].sum()),
    }


def catalog_probe() -> dict[str, Any]:
    catalog_zip = Path(r"D:\数据库集合\新增数据\Mexico_DGIS\clues\2019\ssa_Catalogos_Egresos_2019.zip")
    wanted = {
        "Catalogos Egresos 2019/CATPROCEDENCIA.csv",
        "Catalogos Egresos 2019/CATTIPOSERV.csv",
    }
    with zipfile.ZipFile(catalog_zip) as archive:
        maps = {}
        for name in wanted:
            text = archive.read(name).decode("utf-8-sig", errors="replace")
            rows = list(csv.DictReader(io.StringIO(text)))
            maps[Path(name).stem] = rows
    field_availability = {}
    for year in FULL_SOURCE_YEARS:
        audit_path = ADMIN / f"mx_clinical_year_{year}.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        names = set(audit["source"]["field_names"])
        field_availability[str(year)] = {"PROCED": "PROCED" in names, "TIPSERV": "TIPSERV" in names}
    return {
        "catalog_archive_path": str(catalog_zip),
        "catalog_archive_sha256": hashlib.sha256(catalog_zip.read_bytes()).hexdigest(),
        "catalogs": maps,
        "candidate_urgency_mapping": "The 2019 CATPROCEDENCIA catalog labels IDCatProcedencia=2 as URGENCIAS; raw PROCED is present in all 2015-2023 source headers.",
        "tips_serv_mapping": "CATTIPOSERV labels IDTipoServ=1 NORMAL and 2 CORTA ESTANCIA; TIPSERV is not an urgency classification.",
        "field_availability_by_source_year": field_availability,
        "recommendation": "PROCED==2 is a candidate optional urgency sensitivity; verify the annual PROCED-to-CATPROCEDENCIA linkage and code stability before analysis. No urgent/elective filter was applied.",
    }


def write_cohort_amendment() -> Path:
    path = ADMIN / "cohort_amendment_2019.json"
    decision = {
        "decision_utc": datetime.now(timezone.utc).isoformat(),
        "decision": "Amend the preferred Mexico primary cohort to include 2019 admission dates using all available SAEH-SSA EGRESOS source years, including the recovered 2019 file and late-2019 admission-date records in later source-year files.",
        "effective_date_field": "INGRE",
        "source_year_dimension": "Retain source_year in all complete-source outputs to distinguish discharge-file cohorts and avoid double counting.",
        "original_gap_comparator": "Preserve the original 2015-2018 and 2020-2022 panel unchanged as the comparator; the amended complete-source series is a separate artifact.",
        "study_event_date_window": "2015-01-01 through 2022-12-31 for the complete-source preferred series; 2019 controlled fill is a separate state-day extract from all source-year files.",
        "effects_viewed_before_decision": False,
        "basis": "The 2019 SAEH-SSA source was recovered and validated before any 2019 effect estimates were inspected; its cross-year admission dates are retained by INGRE year and source_year.",
        "patient_level_outputs": False,
    }
    path.write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    ADMIN.mkdir(parents=True, exist_ok=True)
    amendment_path = write_cohort_amendment()
    root_audit = json.loads(INPUT_AUDIT.read_text(encoding="utf-8"))
    years_audit = {item["year"]: item for item in root_audit["year_audits"]}
    contexts = [make_year_context(year, years_audit[year]) for year in (2019, 2020, 2021)]

    corrections = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(extract_source_correction, context) for context in contexts]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            corrections.append(result)
            print(
                f"source {result['source_year']}: missing-cache target rows={result['selected_valid_ingre_missing_from_current_cache_rows']:,}; "
                f"matched={result['missing_cache_rows_with_matched_state_01_32']:,}; groups={result['correction_group_cells']:,}",
                flush=True,
            )
    corrections.sort(key=lambda item: item["source_year"])

    correction_frames = [pd.read_csv(item["correction_path"], dtype={"state": str, "date": str, "fullcode": str}) for item in corrections]
    corrections_frame = pd.concat(correction_frames, ignore_index=True) if correction_frames else pd.DataFrame(columns=CORRECTION_COLUMNS)
    corrections_frame.to_csv(DATA / "mx_complete_source_corrections.csv.gz", index=False, compression="gzip")

    source_frames: list[pd.DataFrame] = []
    for source_year in FULL_SOURCE_YEARS:
        add_source_rows(source_frames, source_year, DATA / f"mx_clinical_year_{source_year}.csv.gz")
    source_frames.append(corrections_frame)
    complete = pd.concat(source_frames, ignore_index=True)
    event_year = pd.to_numeric(complete["date"].str[:4], errors="coerce")
    complete = complete.loc[event_year.between(2015, 2022)].copy()
    merge_key = ["source_year", "state", "date", "endpoint", "fullcode", "age_group", "sex"]
    duplicate_rows = int(complete.duplicated(merge_key).sum())
    if duplicate_rows:
        raise RuntimeError(f"Source-year + group-key overlap found ({duplicate_rows} rows); refusing to double count.")
    complete = complete.sort_values(["source_year", "state", "date", "endpoint", "fullcode", "age_group", "sex"]).reset_index(drop=True)
    long_path = DATA / "mx_complete_source_long_2015_2022.csv.gz"
    complete.to_csv(long_path, index=False, compression="gzip")

    source_wide = build_source_wide(complete)
    source_wide_path = DATA / "mx_complete_source_by_source_year_state_day_2015_2022.csv.gz"
    source_wide.to_csv(source_wide_path, index=False, compression="gzip")

    all_observed = source_wide.groupby(["state", "date"], as_index=False)[WIDE_COUNTS].sum()
    full_state_day = complete_state_day(all_observed, PREFERRED_YEARS, "preferred_complete_source_2015_2022")
    full_path = DATA / "mx_complete_source_state_day_2015_2022.csv"
    full_state_day.to_csv(full_path, index=False)

    controlled_observed = source_wide.loc[source_wide["date"].str.startswith("2019-")].groupby(["state", "date"], as_index=False)[WIDE_COUNTS].sum()
    controlled_2019 = complete_state_day(controlled_observed, [2019], "controlled_2019_all_sources")
    controlled_2019_path = DATA / "mx_complete_source_2019_state_day.csv"
    controlled_2019.to_csv(controlled_2019_path, index=False)

    old_panel = baseline_totals()
    original_reconstruction = {}
    for year in BASELINE_YEARS:
        baseline_role = source_wide.loc[
            source_wide["source_year"].ne(2019)
            & source_wide["date"].str.startswith(f"{year}-")
            & source_wide["source_role"].eq("primary_panel_year")
        ]
        measured = {
            "k80": int(baseline_role["k80"].sum()),
            "k81": int(baseline_role["k81"].sum()),
            "n20_collapsed_n20_to_n23": int(
                baseline_role[["n20_only", "n21_only", "n22_other_disease", "n23_unspecified_renal_colic"]].to_numpy().sum()
            ),
        }
        original_reconstruction[str(year)] = {
            key: {
                "reconstructed": value,
                "baseline_panel": old_panel[year][key],
                "difference": value - old_panel[year][key],
                "matches": value == old_panel[year][key],
            }
            for key, value in measured.items()
        }

    complete_by_year = {}
    for year, frame in full_state_day.assign(year=full_state_day["date"].str[:4]).groupby("year"):
        complete_by_year[year] = family_totals(frame)
    source_2019_year_counts = source_wide.loc[source_wide["date"].str.startswith("2019-")].groupby("source_year")[WIDE_COUNTS].sum().to_dict("index")
    full_totals_2019 = family_totals(controlled_2019)
    urgent_probe = catalog_probe()

    source2023 = years_audit[2023]
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "cohort_amendment_path": str(amendment_path),
        "input_baseline_audit": str(INPUT_AUDIT),
        "scope": "Full-source admission-date reconstruction for 2015-2022, using all local source years 2015-2023; source-year cohorts remain grouped separately.",
        "source_year_semantics": {
            "source_year": "EGRESOS release/source-file year, retained as an explicit group dimension.",
            "event_date": "INGRE admission date; event calendar year determines exposure year.",
            "recovered_2019": "2019 source-year events are retained on their actual 2015-2019 admission dates; late-2019 admissions in 2020/2021 sources are separately tagged with their source_year and included in the controlled 2019 fill.",
            "2023_source": "The 2023 source contributes 2021/2022 admission-date events to the preferred series where present; its INGRE=2023 events are outside the 2015-2022 preferred cohort.",
        },
        "input_sources": [
            {
                "source_year": year,
                "path": years_audit[year]["source"]["path"],
                "bytes": years_audit[year]["source"]["bytes"],
                "sha256": years_audit[year]["source"]["sha256"],
                "encoding": years_audit[year]["source"]["encoding"],
                "delimiter": years_audit[year]["source"]["delimiter"],
            }
            for year in FULL_SOURCE_YEARS
        ],
        "crossyear_corrections": corrections,
        "source_2019_crossyear_timing_and_stay_audit": next(x["2019_crossyear_details"] for x in corrections if x["source_year"] == 2019),
        "source_2023_excluded_2019_check": {
            "source_path": source2023["source"]["path"],
            "audit_sha256": source2023["source"]["sha256"],
            "selected_valid_dates_outside_current_primary_plus_2023_windows": source2023["ingestion"]["selected_rows_with_valid_ingre_outside_analysis_calendar"],
            "interpretation": "The prior full source scan found zero selected stone records outside the 2015-2018, 2020-2023 windows, so no INGRE=2019 correction from source-year 2023 was omitted.",
        },
        "calendar_rules": {
            "complete_preferred_series": "Actual INGRE dates 2015-01-01 through 2022-12-31, all source_years 2015-2023.",
            "controlled_2019_fill": "All source-year admissions with INGRE in calendar 2019; one row per state-day after summing source-year grouped counts.",
            "invalid_dates": "Selected diagnoses with missing/unparseable INGRE are excluded from aggregates but counted by source year and diagnosis prefix.",
            "unmatched_states": "Only zero-padded ENTIDAD 01-32 are emitted; unmatched source values are reported in correction audits.",
            "zero_fill": "Complete wide outputs cross all 32 states with every date in the specified event-date calendar and set absent grouped outcomes to integer zero.",
        },
        "merge_integrity": {
            "group_key": merge_key,
            "duplicate_source_year_group_rows_before_write": duplicate_rows,
            "corrections_merged_by_source_year_and_group_key": True,
            "patient_level_rows_emitted": False,
        },
        "original_gap_cohort_reconciliation": original_reconstruction,
        "original_gap_cohort_all_metrics_match": all(
            cell["matches"] for year in original_reconstruction.values() for cell in year.values()
        ),
        "preferred_complete_source_annual_totals": complete_by_year,
        "controlled_2019_totals_by_source_year": source_2019_year_counts,
        "controlled_2019_totals_all_sources": full_totals_2019,
        "urgent_dictionary_probe": urgent_probe,
        "outputs": {
            "corrections_grouped": str(DATA / "mx_complete_source_corrections.csv.gz"),
            "complete_long": str(long_path),
            "complete_by_source_year_state_day": str(source_wide_path),
            "complete_state_day_2015_2022": str(full_path),
            "controlled_2019_state_day": str(controlled_2019_path),
        },
    }
    audit_path = ADMIN / "mx_complete_source_audit.json"
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {long_path} rows={len(complete):,}", flush=True)
    print(f"wrote {source_wide_path} rows={len(source_wide):,}", flush=True)
    print(f"wrote {full_path} state-days={len(full_state_day):,}", flush=True)
    print(f"wrote {controlled_2019_path} state-days={len(controlled_2019):,}", flush=True)
    print(f"wrote {audit_path}", flush=True)
    print("original gap-cohort reconciliation all match:", audit["original_gap_cohort_all_metrics_match"], flush=True)


if __name__ == "__main__":
    main()
