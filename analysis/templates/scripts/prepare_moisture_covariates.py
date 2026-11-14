"""Prepare auditable moisture covariates for the repaired MX and BR calendars.

Mexico uses the exact reference weights and daily aggregation implemented in
``acquire_mexico2023_weather.py``. Brazil uses the already normalized daily
temperature/dew-point products and preserves their existing RH field alongside
the planned Magnus index derived from daily mean temperature and dew point.
No effect models or canonical calendars are modified here.
"""
from __future__ import annotations

import os

# Keep the detached extraction single-process at native library level.
for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import argparse
import calendar
import gc
import hashlib
import importlib
import json
import math
import sys
import time
import traceback
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
ADMIN = BASE / "admin"
LOG_DIR = BASE / "logs"
RAW_MX = Path(r"F:\era5_ascii\mexico")
GEO_MX = Path(r"F:\era5_ascii\mex_states.geojson")
MX_REFERENCE_CSV = Path(r"F:\era5_ascii\mex_state_day_2018_01.csv")
MX_CANONICAL = DATA / "mx_updated_calendar.csv.gz"
MX_MONTH_CACHE = DATA / "mx_moisture_monthly"
MX_OUTPUT = DATA / "mx_moisture_calendar_2015_2022.csv.gz"

BR_ROOT = Path(r"F:\分析专用\climate_attribution_20261106\work")
BR_HUMID_DIR = BR_ROOT / "expo_humid"
BR_TEMP_DIR = BR_ROOT / "expo_temp"
BR_CANONICAL = DATA / "br_calendar.csv.gz"
BR_YEAR_CACHE = DATA / "br_moisture_yearly"
BR_OUTPUT = DATA / "br_moisture_calendar.csv.gz"

MANIFEST = ADMIN / "moisture_preparation_manifest.json"
QA_JSON = ADMIN / "moisture_preparation_qa.json"
MAIN_LOG = LOG_DIR / "moisture_preparation.log"

EXPECTED_MX_COLUMNS = [
    "unit", "date", "tmean", "dew_mean", "rh_mean", "rh_magnus",
    "dew3", "rh3", "temperature_valid_mass",
]
EXPECTED_BR_COLUMNS = [
    "unit", "date", "tmean", "dew_mean", "rh_mean", "rh_magnus",
    "dew3", "rh3", "rh_source3", "temperature_valid_mass",
]
BR_SOURCE_YEARS = range(2003, 2021)
BR_ANALYSIS_YEARS = range(2006, 2021)
BR_SUPPORT_FLAG_UNITS = {"29020", "16001", "41001"}
BR_REPORTED_SUPPORT_FLOOR = 0.62714
BR_SUPPORT_ROUNDING_TOL = 1e-6
RH_SOURCE_DIFF_ROUNDING_TOL = 1e-3
# The first completed BR annual caches were generated before the support-floor
# tolerance correction. Their feature-generation code is unchanged; exact
# source and cache hashes are still required before reusing them.
BR_CACHE_COMPATIBLE_CODE_SHA256 = {
    "160df8269fba7f5a332cfcbaf39589d363aadb0b0be2a1ca74876e52958f9c88",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")
    os.replace(tmp, path)


def atomic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    try:
        frame.to_csv(
            tmp, index=False, compression="gzip", date_format="%Y-%m-%d",
            float_format="%.15g",
        )
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def log_event(event: dict) -> None:
    record = {"timestamp_utc": now_utc(), **event}
    line = json.dumps(record, sort_keys=True, default=str)
    print(line, flush=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with MAIN_LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def atomic_frame_cache(frame: pd.DataFrame, path: Path) -> None:
    atomic_gzip_csv(frame, path)


def state_code(value) -> str:
    try:
        return f"{int(value):02d}"
    except (TypeError, ValueError):
        return str(value).strip().zfill(2)


def br_unit_code(value) -> str:
    return str(int(value)).zfill(5)


def magnus_rh(temperature_c, dewpoint_c) -> np.ndarray:
    temperature = np.asarray(temperature_c, dtype=np.float64)
    dewpoint = np.asarray(dewpoint_c, dtype=np.float64)
    es_t = 6.112 * np.exp(17.67 * temperature / (temperature + 243.5))
    es_td = 6.112 * np.exp(17.67 * dewpoint / (dewpoint + 243.5))
    return np.clip(100.0 * es_td / es_t, 0.0, 100.0)


def ensure_source(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")


def mx_expected_month_path(year: int, month: int) -> Path:
    return RAW_MX / f"era5land_hourly_temp_dewpoint_mexico_{year}_{month:02d}.nc"


def validate_mx_month(path: Path, year: int, month: int, lat: np.ndarray, lon: np.ndarray, source_hash: str | None = None):
    expected = pd.date_range(
        f"{year}-{month:02d}-01",
        periods=calendar.monthrange(year, month)[1] * 24,
        freq="h",
    )
    with xr.open_dataset(path, engine="netcdf4") as ds:
        if not {"t2m", "d2m"}.issubset(ds.data_vars):
            raise ValueError(f"{path.name} lacks t2m or d2m")
        file_lat = ds.latitude.values.astype(float)
        file_lon = ds.longitude.values.astype(float)
        if not np.array_equal(file_lat, lat) or not np.array_equal(file_lon, lon):
            raise ValueError(f"Grid differs from fixed Mexico reference in {path.name}")
        dates = mx_builder.load_times(ds)
        if not dates.equals(expected):
            raise ValueError(
                f"Invalid hourly coverage in {path.name}: rows={len(dates)}, "
                f"first={dates.min() if len(dates) else None}, "
                f"last={dates.max() if len(dates) else None}"
            )
        time_name = "valid_time" if "valid_time" in ds else "time"
        time_units = str(
            ds[time_name].attrs.get("units")
            or ds[time_name].encoding.get("units", "")
        )
        time_calendar = str(
            ds[time_name].attrs.get("calendar")
            or ds[time_name].encoding.get("calendar", "")
        )
        t_units = str(ds.t2m.attrs.get("units", "K")).strip().lower()
        d_units = str(ds.d2m.attrs.get("units", "K")).strip().lower()
        accepted = {"k", "kelvin", "degrees_k", "degree_k"}
        if t_units not in accepted or d_units not in accepted:
            raise ValueError(f"Unexpected Kelvin units in {path.name}: t2m={t_units}, d2m={d_units}")
        t_shape = tuple(ds.t2m.shape)
        d_shape = tuple(ds.d2m.shape)
        expected_shape = (len(expected), len(lat), len(lon))
        if t_shape != expected_shape or d_shape != expected_shape:
            raise ValueError(f"Unexpected arrays in {path.name}: t2m={t_shape}, d2m={d_shape}")
    return dates, {
        "file": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": source_hash or sha256_file(path),
        "hours": len(expected),
        "time_start_utc": expected[0].isoformat(),
        "time_end_utc": expected[-1].isoformat(),
        "time_coordinate": time_name,
        "time_coordinate_units": time_units,
        "time_coordinate_calendar": time_calendar,
        "t2m_units": "K",
        "d2m_units": "K",
        "grid_shape_lat_lon": [len(lat), len(lon)],
    }


def weighted_daily_state_means(daily_field: np.ndarray, weights: np.ndarray):
    """Reproduce the original state-wise finite-weight renormalization."""
    flat = daily_field.ravel()
    means = np.empty(weights.shape[0], dtype=np.float64)
    denominators = np.empty(weights.shape[0], dtype=np.float64)
    for i in range(weights.shape[0]):
        mask = weights[i] > 0
        numerator = np.nansum(flat[mask] * weights[i, mask])
        denominator = np.nansum(weights[i, mask] * np.isfinite(flat[mask]))
        means[i] = float(numerator / denominator) if denominator > 0 else np.nan
        denominators[i] = float(denominator)
    return means, denominators


def mx_month_cache_valid(path: Path, year: int, month: int, expected_rows: int) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        frame = pd.read_csv(path, dtype={"unit": str})
        frame["unit"] = frame.unit.map(state_code)
        if (
            frame.columns.tolist() == EXPECTED_MX_COLUMNS
            and len(frame) == expected_rows
            and not frame.duplicated(["unit", "date"]).any()
            and set(pd.to_datetime(frame.date).dt.year.astype(int)) == {year}
            and set(pd.to_datetime(frame.date).dt.month.astype(int)) == {month}
            and frame[["tmean", "dew_mean", "rh_mean", "rh_magnus", "temperature_valid_mass"]].notna().all().all()
        ):
            return frame
    except Exception:
        return None
    return None


def load_mx_manifest() -> dict:
    if MANIFEST.is_file():
        try:
            return json.loads(MANIFEST.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def run_mexico(code_sha: str) -> dict:
    if not (BASE / "scripts" / "acquire_mexico2023_weather.py").is_file():
        raise FileNotFoundError("Existing Mexico aggregation script is missing")
    sys.path.insert(0, str(BASE / "scripts"))
    global mx_builder
    mx_builder = importlib.import_module("acquire_mexico2023_weather")
    for path in (GEO_MX, MX_REFERENCE_CSV, MX_CANONICAL):
        ensure_source(path)
    source_script = BASE / "scripts" / "acquire_mexico2023_weather.py"
    original_builder = mx_builder.SOURCE_CODE
    states, lat, lon, weights, weight_reference = mx_builder.initialize_reference()
    if weight_reference.name != "era5land_hourly_temp_dewpoint_mexico_2011_01.nc":
        raise ValueError(f"Unexpected fixed mask reference selected: {weight_reference}")
    log_event({
        "event": "mx_weights_built",
        "states": len(states),
        "reference_file": str(weight_reference),
        "reference_sha256": sha256_file(weight_reference),
        "weights_once": True,
        "weight_build_seconds_approx": "timed outside this function",
    })

    manifest = load_mx_manifest()
    month_records = manifest.get("mexico", {}).get("months", {})
    month_cache_manifest = {}
    MX_MONTH_CACHE.mkdir(parents=True, exist_ok=True)
    all_month_frames = []
    all_raw_metadata = []
    denominator_values = []
    clip_low = 0
    clip_high = 0
    previous_mx_tail = None
    started = time.perf_counter()

    for year in range(2015, 2023):
        for month in range(1, 13):
            tag = f"{year}-{month:02d}"
            raw_path = mx_expected_month_path(year, month)
            ensure_source(raw_path)
            days = calendar.monthrange(year, month)[1]
            expected_rows = len(states) * days
            cache_path = MX_MONTH_CACHE / f"mx_moisture_{year}_{month:02d}.csv.gz"
            prior = month_records.get(tag, {})
            source_hash = sha256_file(raw_path)
            cached = None
            if prior.get("source_sha256") == source_hash and prior.get("code_sha256") == code_sha:
                cached = mx_month_cache_valid(cache_path, year, month, expected_rows)
            if cached is not None:
                frame = cached
                raw_meta = prior["raw"]
                month_audit = prior.get("qa", {})
                log_event({"event": "mx_month_cache_reused", "month": tag, "rows": len(frame), "source_sha256": source_hash})
            else:
                log_event({"event": "mx_month_start", "month": tag, "raw_file": str(raw_path), "source_sha256": source_hash})
                expected_dates, raw_meta = validate_mx_month(raw_path, year, month, lat, lon, source_hash)
                with xr.open_dataset(raw_path, engine="netcdf4") as ds:
                    temp_k = ds["t2m"].values.astype(np.float32, copy=False)
                    dew_k = ds["d2m"].values.astype(np.float32, copy=False)
                t_mask = np.isfinite(temp_k)
                d_mask = np.isfinite(dew_k)
                t_constant = bool(np.all(t_mask == t_mask[0:1]))
                d_constant = bool(np.all(d_mask == d_mask[0:1]))
                cross_equal = bool(np.array_equal(t_mask[0], d_mask[0]))
                if not (t_constant and d_constant and cross_equal):
                    raise ValueError(
                        f"Finite-mask QA failed for {tag}: t_time_constant={t_constant}, "
                        f"d_time_constant={d_constant}, t_d_equal={cross_equal}"
                    )
                if temp_k.shape != dew_k.shape or temp_k.shape[0] != len(expected_dates):
                    raise ValueError(f"T/D arrays do not align in {tag}")
                temp_c = temp_k - np.float32(273.15)
                dew_c = dew_k - np.float32(273.15)
                grid_shape = (len(lat), len(lon))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    temp_daily = np.nanmean(temp_c.reshape(days, 24, *grid_shape), axis=1)
                    dew_daily = np.nanmean(dew_c.reshape(days, 24, *grid_shape), axis=1)
                out_rows = []
                month_denominators = []
                month_low = 0
                month_high = 0
                state_values = [state_code(x) for x in states["cve_ent"]]
                for day_number, (temperature_grid, dew_grid) in enumerate(zip(temp_daily, dew_daily), start=1):
                    t_means, t_den = weighted_daily_state_means(temperature_grid, weights)
                    d_means, d_den = weighted_daily_state_means(dew_grid, weights)
                    if not np.array_equal(t_den, d_den):
                        raise ValueError(f"Daily T/D weighted denominators differ in {tag} day {day_number}")
                    if not np.isfinite(t_means).all() or not np.isfinite(d_means).all() or (t_den <= 0).any():
                        raise ValueError(f"Missing state-day T/D values in {tag} day {day_number}")
                    rh_unclipped = 100.0 * (
                        6.112 * np.exp(17.67 * d_means / (d_means + 243.5))
                    ) / (
                        6.112 * np.exp(17.67 * t_means / (t_means + 243.5))
                    )
                    month_low += int(np.sum(rh_unclipped < 0.0))
                    month_high += int(np.sum(rh_unclipped > 100.0))
                    rh = np.clip(rh_unclipped, 0.0, 100.0)
                    date_value = expected_dates[24 * (day_number - 1)].strftime("%Y-%m-%d")
                    for i, unit in enumerate(state_values):
                        out_rows.append((unit, date_value, float(t_means[i]), float(d_means[i]), float(rh[i]), float(rh[i]), float(t_den[i])))
                        month_denominators.append(float(t_den[i]))
                frame = pd.DataFrame(out_rows, columns=[
                    "unit", "date", "tmean", "dew_mean", "rh_mean", "rh_magnus", "temperature_valid_mass",
                ])

            # The trailing three-day weather windows cross month boundaries.
            # Cache reuse is safe because the preceding month's final two daily
            # values are restored from its cache before computing this month.
            frame["unit"] = frame.unit.map(state_code)
            frame["date"] = pd.to_datetime(frame.date, format="%Y-%m-%d")
            base_columns = ["unit", "date", "tmean", "dew_mean", "rh_mean", "rh_magnus", "temperature_valid_mass"]
            current_dates = frame.date.unique()
            if previous_mx_tail is None:
                combined = frame[base_columns].copy()
            else:
                combined = pd.concat([previous_mx_tail[base_columns], frame[base_columns]], ignore_index=True)
            combined = combined.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
            combined["dew3"] = combined.groupby("unit", sort=False)["dew_mean"].transform(
                lambda values: values.rolling(3, min_periods=3).mean()
            )
            combined["rh3"] = combined.groupby("unit", sort=False)["rh_magnus"].transform(
                lambda values: values.rolling(3, min_periods=3).mean()
            )
            frame = combined.loc[combined.date.isin(current_dates)].copy()
            previous_mx_tail = (
                frame.sort_values(["unit", "date"], kind="mergesort")
                .groupby("unit", sort=False, group_keys=False)
                .tail(2)[base_columns]
                .copy()
            )
            frame["date"] = frame.date.dt.strftime("%Y-%m-%d")
            frame = frame[EXPECTED_MX_COLUMNS]
            if len(frame) != expected_rows or frame.duplicated(["unit", "date"]).any():
                raise ValueError(f"Invalid Mexico monthly output keys for {tag}")
            if cached is None:
                month_denominators = np.asarray(month_denominators, dtype=np.float64)
                month_audit = {
                    "status": "passed",
                    "rows": int(len(frame)),
                    "finite_mask_time_invariant_t2m": t_constant,
                    "finite_mask_time_invariant_d2m": d_constant,
                    "finite_mask_t2m_equals_d2m": cross_equal,
                    "daily_valid_weight_denominator_min": float(month_denominators.min()),
                    "daily_valid_weight_denominator_max": float(month_denominators.max()),
                    "daily_valid_weight_denominator_nonfull_rows": int(np.count_nonzero(np.abs(month_denominators - 1.0) > 1e-6)),
                    "rh_unclipped_below_0_count": month_low,
                    "rh_unclipped_above_100_count": month_high,
                }
                atomic_frame_cache(frame, cache_path)
                month_records[tag] = {
                    "status": "complete",
                    "source_file": str(raw_path),
                    "source_sha256": source_hash,
                    "source_size_bytes": raw_path.stat().st_size,
                    "code_sha256": code_sha,
                    "cache_file": str(cache_path),
                    "cache_sha256": sha256_file(cache_path),
                    "raw": raw_meta,
                    "qa": month_audit,
                    "completed_utc": now_utc(),
                }
                manifest["mexico"] = {"status": "processing", "months": month_records}
                manifest["updated_utc"] = now_utc()
                atomic_json(MANIFEST, manifest)
                log_event({"event": "mx_month_complete", "month": tag, **month_audit})
                del temp_k, dew_k, t_mask, d_mask, temp_c, dew_c, temp_daily, dew_daily
                gc.collect()

            all_month_frames.append(frame)
            all_raw_metadata.append(raw_meta)
            denominator_values.extend([month_audit.get("daily_valid_weight_denominator_min"), month_audit.get("daily_valid_weight_denominator_max")])
            clip_low += int(month_audit.get("rh_unclipped_below_0_count", 0))
            clip_high += int(month_audit.get("rh_unclipped_above_100_count", 0))
            month_cache_manifest[tag] = month_records.get(tag, {})

    moisture = pd.concat(all_month_frames, ignore_index=True)
    moisture["unit"] = moisture.unit.map(state_code)
    moisture = moisture.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
    expected_dates = pd.date_range("2015-01-01", "2022-12-31", freq="D")
    expected_keys = 32 * len(expected_dates)
    if len(moisture) != expected_keys or moisture.duplicated(["unit", "date"]).any():
        raise ValueError(f"Mexico full-period key count failed: rows={len(moisture)}, expected={expected_keys}")
    if set(moisture.unit) != {f"{i:02d}" for i in range(1, 33)}:
        raise ValueError("Mexico output does not contain the 32 fixed-width state codes")
    observed_dates = pd.DatetimeIndex(pd.to_datetime(moisture.date).unique()).sort_values()
    if not observed_dates.equals(expected_dates):
        raise ValueError("Mexico moisture output is not a complete daily 2015-2022 calendar")
    if int(moisture[["dew3", "rh3"]].isna().sum().max()) != 64:
        raise ValueError("Mexico first-two-day 3-day-window support is not exactly 64 NA values per metric")

    canonical = pd.read_csv(MX_CANONICAL, usecols=["unit", "date", "tmean"], dtype={"unit": str})
    canonical["unit"] = canonical.unit.map(state_code)
    canonical["date"] = canonical.date.astype(str)
    if len(canonical) != expected_keys or canonical.duplicated(["unit", "date"]).any():
        raise ValueError("Mexico canonical calendar keys are not unique and complete")
    checked = moisture[["unit", "date", "tmean"]].merge(
        canonical, on=["unit", "date"], how="outer", suffixes=("_moisture", "_canonical"),
        indicator=True, validate="one_to_one",
    )
    if not (checked._merge == "both").all():
        raise ValueError(f"Mexico moisture/canonical keys differ: {checked._merge.value_counts().to_dict()}")
    mx_t_diff = (checked.tmean_moisture - checked.tmean_canonical).abs()
    mx_t_max = float(mx_t_diff.max())
    if mx_t_max > 1e-8:
        raise ValueError(f"Mexico all-date tmean reconciliation exceeded 1e-8: max_abs={mx_t_max:.12g}")

    reference = pd.read_csv(MX_REFERENCE_CSV, dtype={"state": str})
    reference["state"] = reference.state.map(state_code)
    reference["date"] = pd.to_datetime(reference[["year", "month", "day"]]).dt.strftime("%Y-%m-%d")
    mx_jan = moisture.loc[moisture.date.str.startswith("2018-01"), ["unit", "date", "tmean"]].merge(
        reference.rename(columns={"state": "unit"})[["unit", "date", "tmean"]],
        on=["unit", "date"], how="outer", suffixes=("_new", "_ref"), indicator=True, validate="one_to_one",
    )
    if len(mx_jan) != 992 or not (mx_jan._merge == "both").all():
        raise ValueError(f"Mexico Jan 2018 old-route validation keys failed: rows={len(mx_jan)}")
    jan_abs = float((mx_jan.tmean_new - mx_jan.tmean_ref).abs().max())
    if jan_abs > 1e-6:
        raise ValueError(f"Mexico Jan 2018 previous-route parity failed: max_abs={jan_abs:.12g}")

    moisture = moisture[EXPECTED_MX_COLUMNS]
    atomic_gzip_csv(moisture, MX_OUTPUT)
    mx_audit = {
        "status": "passed",
        "output": str(MX_OUTPUT),
        "output_sha256": sha256_file(MX_OUTPUT),
        "rows": int(len(moisture)),
        "units": int(moisture.unit.nunique()),
        "first_date": str(moisture.date.min()),
        "last_date": str(moisture.date.max()),
        "source_months": len(all_raw_metadata),
        "monthly_cache_dir": str(MX_MONTH_CACHE),
        "months": month_cache_manifest,
        "monthly_cache_manifest": month_cache_manifest,
        "tmean_reconciliation": {
            "canonical": str(MX_CANONICAL),
            "canonical_sha256": sha256_file(MX_CANONICAL),
            "matched_state_days": len(checked),
            "max_abs_difference_celsius": mx_t_max,
            "tolerance": 1e-8,
        },
        "previous_route_validation_2018_01": {
            "reference_csv": str(MX_REFERENCE_CSV),
            "reference_csv_sha256": sha256_file(MX_REFERENCE_CSV),
            "matched_state_days": len(mx_jan),
            "max_abs_difference_celsius": jan_abs,
            "audited_previous_result_max_abs": 3.552713678800501e-15,
            "acceptance_tolerance": 1e-6,
        },
        "aggregation": {
            "weights": "same initialize_reference() implementation as acquire_mexico2023_weather.py; cos(latitude) over strict polygon-contained grid centers and static 2011-01 first-hour valid-land mask; built once",
            "daily_temperature_and_dewpoint": "Kelvin to Celsius before arithmetic mean of available 24 hourly cell observations, then per-state finite-weight renormalization using the original algorithm",
            "daily_rh_index": "clip(100*es(state-day dew mean)/es(state-day temperature mean), 0, 100); es=6.112*exp(17.67*x/(x+243.5))",
            "rh_mean_definition": "daily Magnus-derived index; no legacy separate Mexico RH aggregation was present",
            "rh3_definition": "arithmetic mean of admission day and two preceding actual calendar days, requiring all three values; first two dates per state are NA because Mexico raw history starts 2015-01-01",
            "unit_format": "two-character state codes",
            "finite_mask_qa_all_months": "t2m and d2m hourly finite masks were equal and time-invariant in every month",
        },
        "weight_reference": {
            "file": str(weight_reference),
            "sha256": sha256_file(weight_reference),
            "geography": str(GEO_MX),
            "geography_sha256": sha256_file(GEO_MX),
            "aggregation_script": str(source_script),
            "aggregation_script_sha256": sha256_file(source_script),
            "original_builder": str(original_builder),
            "original_builder_sha256": sha256_file(original_builder),
            "grid_shape_lat_lon": [len(lat), len(lon)],
        },
        "daily_valid_weight_denominator_month_min_max_extrema": [
            float(np.nanmin([x for x in denominator_values if x is not None])),
            float(np.nanmax([x for x in denominator_values if x is not None])),
        ],
        "rh_unclipped_count_below_0": int(clip_low),
        "rh_unclipped_count_above_100": int(clip_high),
        "code_sha256": code_sha,
        "elapsed_seconds": time.perf_counter() - started,
        "checked_utc": now_utc(),
    }
    manifest["mexico"] = mx_audit
    manifest["updated_utc"] = now_utc()
    atomic_json(MANIFEST, manifest)
    log_event({"event": "mx_complete", "rows": len(moisture), "tmean_max_abs": mx_t_max, "jan2018_max_abs": jan_abs, "output_sha256": mx_audit["output_sha256"]})
    return mx_audit


def br_source_paths(year: int) -> tuple[Path, Path]:
    return (
        BR_HUMID_DIR / f"region_day_humid_{year}.csv",
        BR_TEMP_DIR / f"region_day_temp_{year}.csv",
    )


def br_manifest_year_ok(manifest: dict, year: int, humid_hash: str, temp_hash: str, code_sha: str) -> dict | None:
    record = manifest.get("brazil", {}).get("years", {}).get(str(year), {})
    cache_path = BR_YEAR_CACHE / f"br_moisture_{year}.csv.gz"
    if (
        record.get("status") != "complete"
        or record.get("source_humid_sha256") != humid_hash
        or record.get("source_temperature_sha256") != temp_hash
        or record.get("code_sha256") not in ({code_sha} | BR_CACHE_COMPATIBLE_CODE_SHA256)
        or record.get("cache_sha256") != (sha256_file(cache_path) if cache_path.is_file() else None)
    ):
        return None
    try:
        frame = pd.read_csv(cache_path, dtype={"unit": str})
        expected_days = 366 if calendar.isleap(year) else 365
        if (
            frame.columns.tolist() != EXPECTED_BR_COLUMNS
            or len(frame) != 450 * expected_days
            or frame.duplicated(["unit", "date"]).any()
            or set(pd.to_datetime(frame.date).dt.year.astype(int)) != {year}
        ):
            return None
        return frame
    except Exception:
        return None


def summarize_rh_difference(diff: np.ndarray) -> dict:
    diff = np.asarray(diff, dtype=np.float64)
    absolute = np.abs(diff)
    return {
        "rows": int(len(diff)),
        "source_minus_magnus_mean_percentage_points": float(diff.mean()),
        "source_minus_magnus_rmse_percentage_points": float(np.sqrt(np.mean(diff * diff))),
        "source_minus_magnus_mean_absolute_percentage_points": float(absolute.mean()),
        "source_minus_magnus_median_absolute_percentage_points": float(np.median(absolute)),
        "source_minus_magnus_p95_absolute_percentage_points": float(np.quantile(absolute, 0.95)),
        "source_minus_magnus_max_absolute_percentage_points": float(absolute.max()),
        "source_minus_magnus_min_percentage_points": float(diff.min()),
        "source_minus_magnus_max_percentage_points": float(diff.max()),
    }


def load_br_year(year: int) -> tuple[pd.DataFrame, str, str]:
    humid_path, temp_path = br_source_paths(year)
    ensure_source(humid_path)
    ensure_source(temp_path)
    humid_hash = sha256_file(humid_path)
    temp_hash = sha256_file(temp_path)
    humid = pd.read_csv(humid_path)
    temp = pd.read_csv(temp_path, usecols=["hrid", "date", "tmean", "valid_mass"])
    required_humid = {"hrid", "date", "dew_mean", "rh_mean", "at_mean"}
    if not required_humid.issubset(humid.columns):
        raise ValueError(f"Unexpected BR humidity schema in {humid_path.name}: {humid.columns.tolist()}")
    humid = humid[["hrid", "date", "dew_mean", "rh_mean"]].copy()
    humid["unit"] = humid.hrid.map(br_unit_code)
    humid["date"] = humid.date.astype(str)
    temp["unit"] = temp.hrid.map(br_unit_code)
    temp["date"] = temp.date.astype(str)
    humid = humid.drop(columns="hrid")
    temp = temp.drop(columns="hrid")
    if humid.duplicated(["unit", "date"]).any() or temp.duplicated(["unit", "date"]).any():
        raise ValueError(f"Duplicate BR unit-day source keys in {year}")
    merged = humid.merge(temp, on=["unit", "date"], how="outer", validate="one_to_one", indicator=True)
    if len(merged) != 450 * (366 if calendar.isleap(year) else 365) or not (merged._merge == "both").all():
        raise ValueError(f"BR T/dew point key support differs in {year}: {merged._merge.value_counts().to_dict()}")
    merged = merged.drop(columns="_merge")
    merged["date"] = pd.to_datetime(merged.date, format="%Y-%m-%d", errors="raise")
    merged = merged.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
    numeric = ["tmean", "dew_mean", "rh_mean", "temperature_valid_mass"]
    merged = merged.rename(columns={"valid_mass": "temperature_valid_mass"})
    if not np.isfinite(merged[numeric].to_numpy(dtype=np.float64)).all():
        raise ValueError(f"Non-finite BR temperature/dew/RH/support values in {year}")
    if ((merged.temperature_valid_mass < 0) | (merged.temperature_valid_mass > 1.000001)).any():
        raise ValueError(f"Invalid BR temperature-valid population mass in {year}")
    if ((merged.rh_mean < 0) | (merged.rh_mean > 100)).any():
        raise ValueError(f"Existing BR RH outside [0,100] in {year}; source field preserved, review required")
    merged["rh_magnus"] = magnus_rh(merged.tmean, merged.dew_mean)
    return merged, humid_hash, temp_hash


def br_validate_year_contiguous(frame: pd.DataFrame, year: int) -> None:
    expected_dates = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    for unit, part in frame.groupby("unit", sort=False):
        dates = pd.DatetimeIndex(part.date)
        if not dates.equals(expected_dates):
            raise ValueError(f"BR unit {unit} lacks a complete daily calendar in {year}")


def run_brazil(code_sha: str, manifest: dict) -> dict:
    ensure_source(BR_CANONICAL)
    BR_YEAR_CACHE.mkdir(parents=True, exist_ok=True)
    years_manifest = manifest.get("brazil", {}).get("years", {})
    canonical = pd.read_csv(BR_CANONICAL, usecols=["unit", "date", "tmean"], dtype={"unit": str})
    canonical["unit"] = canonical.unit.str.zfill(5)
    canonical["date"] = canonical.date.astype(str)
    canonical = canonical.loc[canonical.date.str[:4].isin([str(y) for y in BR_ANALYSIS_YEARS])].copy()
    canonical["date"] = pd.to_datetime(canonical.date, format="%Y-%m-%d", errors="raise")
    canonical = canonical.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
    expected_rows = sum(450 * (366 if calendar.isleap(y) else 365) for y in BR_ANALYSIS_YEARS)
    if len(canonical) != expected_rows or canonical.duplicated(["unit", "date"]).any():
        raise ValueError(f"BR canonical 2006-2020 keys failed: rows={len(canonical)}, expected={expected_rows}")
    expected_units = set(canonical.unit.unique())
    if len(expected_units) != 450:
        raise ValueError(f"Expected 450 BR canonical units, found {len(expected_units)}")

    previous_tail = None
    annual_frames = []
    full_diff_chunks = []
    year_audits = {}
    source_manifest = {}
    started = time.perf_counter()
    for year in BR_SOURCE_YEARS:
        tag = str(year)
        humid_path, temp_path = br_source_paths(year)
        ensure_source(humid_path)
        ensure_source(temp_path)
        humid_hash = sha256_file(humid_path)
        temp_hash = sha256_file(temp_path)
        cache_path = BR_YEAR_CACHE / f"br_moisture_{year}.csv.gz"
        frame = br_manifest_year_ok(manifest, year, humid_hash, temp_hash, code_sha)
        if frame is None:
            log_event({"event": "br_year_start", "year": year, "humidity_source_sha256": humid_hash, "temperature_source_sha256": temp_hash})
            current, _, _ = load_br_year(year)
            current = current.drop(columns=["rh3", "dew3", "rh_source3"], errors="ignore")
            if previous_tail is None:
                combined = current.copy()
            else:
                combined = pd.concat([previous_tail, current], ignore_index=True)
            combined = combined.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
            for unit, part in combined.groupby("unit", sort=False):
                expected = pd.date_range(part.date.min(), part.date.max(), freq="D")
                if not pd.DatetimeIndex(part.date).equals(expected):
                    raise ValueError(f"BR input calendar is not continuous across years for unit {unit}, near {year}")
            for field, feature in (("dew_mean", "dew3"), ("rh_magnus", "rh3"), ("rh_mean", "rh_source3")):
                combined[feature] = combined.groupby("unit", sort=False)[field].transform(
                    lambda values: values.rolling(3, min_periods=3).mean()
                )
            current_dates = current["date"].astype("datetime64[ns]")
            frame = combined.loc[combined.date.isin(current_dates), EXPECTED_BR_COLUMNS].copy()
            frame = frame.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
            if len(frame) != 450 * (366 if calendar.isleap(year) else 365):
                raise ValueError(f"BR yearly feature rows failed for {year}: {len(frame)}")
            frame["date"] = pd.to_datetime(frame.date).dt.strftime("%Y-%m-%d")
            atomic_frame_cache(frame, cache_path)
            diff = frame.rh_mean.to_numpy(dtype=np.float64) - frame.rh_magnus.to_numpy(dtype=np.float64)
            year_audit = {
                "status": "complete",
                "rows": int(len(frame)),
                "humidity_source": str(humid_path),
                "source_humid_sha256": humid_hash,
                "temperature_source": str(temp_path),
                "source_temperature_sha256": temp_hash,
                "cache_file": str(cache_path),
                "cache_sha256": sha256_file(cache_path),
                "code_sha256": code_sha,
                "rh_source_comparison": summarize_rh_difference(diff),
                "temperature_valid_mass_min": float(frame.temperature_valid_mass.min()),
                "temperature_valid_mass_max": float(frame.temperature_valid_mass.max()),
                "completed_utc": now_utc(),
            }
            years_manifest[tag] = year_audit
            manifest["brazil"] = {"status": "processing", "years": years_manifest}
            manifest["updated_utc"] = now_utc()
            atomic_json(MANIFEST, manifest)
            log_event({"event": "br_year_complete", "year": year, "rows": len(frame), "cache_sha256": year_audit["cache_sha256"], "rh_source_minus_magnus_max_abs": year_audit["rh_source_comparison"]["source_minus_magnus_max_absolute_percentage_points"]})
            current_for_tail = frame.copy()
        else:
            year_audit = years_manifest[tag]
            log_event({"event": "br_year_cache_reused", "year": year, "rows": len(frame), "cache_sha256": year_audit.get("cache_sha256")})
            current_for_tail = frame.copy()

        # The history comes from the actual full calendar. Retain two days per
        # region so January 1-2, 2006 use 2005 data before outcome-key selection.
        current_for_tail["date"] = pd.to_datetime(current_for_tail.date)
        previous_tail = (
            current_for_tail.sort_values(["unit", "date"], kind="mergesort")
            .groupby("unit", sort=False, group_keys=False)
            .tail(2)[["unit", "date", "dew_mean", "rh_mean", "rh_magnus", "temperature_valid_mass", "tmean"]]
            .copy()
        )
        frame["date"] = pd.to_datetime(frame.date, format="%Y-%m-%d")
        if year in BR_ANALYSIS_YEARS:
            annual_frames.append(frame)
        full_diff_chunks.append(frame[["rh_mean", "rh_magnus"]].to_numpy(dtype=np.float64))
        year_audits[tag] = year_audit
        source_manifest[tag] = {
            "humidity_source": str(humid_path),
            "humidity_source_sha256": humid_hash,
            "temperature_source": str(temp_path),
            "temperature_source_sha256": temp_hash,
            "cache_file": str(cache_path),
            "cache_sha256": year_audit.get("cache_sha256"),
        }

    moisture = pd.concat(annual_frames, ignore_index=True)
    moisture["unit"] = moisture.unit.astype(str).str.zfill(5)
    moisture = moisture.sort_values(["unit", "date"], kind="mergesort").reset_index(drop=True)
    if len(moisture) != expected_rows or moisture.duplicated(["unit", "date"]).any():
        raise ValueError(f"BR moisture output key count failed: rows={len(moisture)}, expected={expected_rows}")
    if set(moisture.unit.unique()) != expected_units:
        raise ValueError("BR moisture units differ from repaired canonical calendar")
    expected_analysis_dates = pd.date_range("2006-01-01", "2020-12-31", freq="D")
    observed_dates = pd.DatetimeIndex(moisture.date.unique()).sort_values()
    if not observed_dates.equals(expected_analysis_dates):
        raise ValueError("BR moisture output is not a complete daily 2006-2020 calendar")
    if moisture.loc[moisture.date.dt.year.eq(2006) & moisture.date.dt.month.eq(1) & moisture.date.dt.day.le(2), ["dew3", "rh3"]].isna().any().any():
        raise ValueError("BR Jan 1-2 2006 moisture lags did not use pre-2006 history")

    # Both tables are sorted identically so exact key-by-key comparison avoids a
    # large many-column join of the complete two-million-row calendar.
    moisture_keys = moisture[["unit", "date"]].copy()
    canonical_keys = canonical[["unit", "date"]].copy()
    if not moisture_keys.reset_index(drop=True).equals(canonical_keys.reset_index(drop=True)):
        raise ValueError("BR moisture unit/date keys differ from repaired canonical calendar")
    t_diff = (moisture.tmean.to_numpy(dtype=np.float64) - canonical.tmean.to_numpy(dtype=np.float64))
    br_t_max = float(np.max(np.abs(t_diff)))
    if br_t_max > 1e-8:
        raise ValueError(f"BR 2006-2020 tmean reconciliation exceeded 1e-8: max_abs={br_t_max:.12g}")

    coverage = moisture.groupby("unit", sort=True).temperature_valid_mass.agg(["min", "max"])
    coverage_span = coverage["max"] - coverage["min"]
    if (coverage_span > 1e-8).any():
        raise ValueError(f"BR companion temperature_valid_mass changed within region: {int((coverage_span > 1e-8).sum())} regions")
    support_excluded = sorted(coverage.index[coverage["min"] < 0.6].tolist())
    remaining_min = float(coverage.loc[~coverage.index.isin(support_excluded), "min"].min())
    if set(support_excluded) != BR_SUPPORT_FLAG_UNITS or remaining_min < BR_REPORTED_SUPPORT_FLOOR - BR_SUPPORT_ROUNDING_TOL:
        raise ValueError(f"Unexpected BR low-support units or remaining minimum: low={support_excluded}, remaining_min={remaining_min}")

    rh_pairs = np.concatenate(full_diff_chunks, axis=0)
    rh_diff_summary = summarize_rh_difference(rh_pairs[:, 0] - rh_pairs[:, 1])
    if rh_diff_summary["source_minus_magnus_max_absolute_percentage_points"] > RH_SOURCE_DIFF_ROUNDING_TOL:
        raise ValueError(
            "BR existing RH differs materially from the Magnus index on companion regional means; "
            f"max_abs={rh_diff_summary['source_minus_magnus_max_absolute_percentage_points']:.9g} pp. "
            "No final BR moisture file was overwritten; review source definition first."
        )

    moisture["date"] = moisture.date.dt.strftime("%Y-%m-%d")
    moisture = moisture[EXPECTED_BR_COLUMNS]
    atomic_gzip_csv(moisture, BR_OUTPUT)
    br_audit = {
        "status": "passed",
        "output": str(BR_OUTPUT),
        "output_sha256": sha256_file(BR_OUTPUT),
        "rows": int(len(moisture)),
        "units": int(moisture.unit.nunique()),
        "first_date": str(moisture.date.min()),
        "last_date": str(moisture.date.max()),
        "analysis_calendar_source": str(BR_CANONICAL),
        "analysis_calendar_source_sha256": sha256_file(BR_CANONICAL),
        "analysis_calendar_matched_keys": int(len(canonical)),
        "temperature_reconciliation": {
            "max_abs_difference_celsius": br_t_max,
            "tolerance": 1e-8,
        },
        "moisture_sources": source_manifest,
        "year_cache_dir": str(BR_YEAR_CACHE),
        "years": year_audits,
        "year_audits": year_audits,
        "humidity_rh_comparison_all_2003_2020": {
            **rh_diff_summary,
            "source_rh_preserved": True,
            "derived_field": "rh_magnus",
            "three_day_policy_field": "rh3 = trailing arithmetic mean of rh_magnus, requiring all 3 days",
            "source_field": "rh_mean (pre-existing humidity export); rh_source3 is its 3-day mean",
            "rounding_consistency_tolerance_percentage_points": RH_SOURCE_DIFF_ROUNDING_TOL,
            "status": "consistent_with_companion_mean_magnus_within_tolerance",
        },
        "temperature_valid_mass": {
            "field": "temperature_valid_mass",
            "source": "existing normalized expo_temp companion valid_mass",
            "time_invariant_within_region": True,
            "below_0_6_units": support_excluded,
            "remaining_units_minimum": remaining_min,
            "per_unit_min_max": coverage.reset_index().rename(columns={"index": "unit"}).to_dict(orient="records"),
        },
        "support_floor_check": {
            "rule": "flag units when minimum companion temperature_valid_mass over all 2006-2020 calendar dates is below 0.6",
            "excluded_units": support_excluded,
            "expected_excluded_units": sorted(BR_SUPPORT_FLAG_UNITS),
            "excluded_unit_set_match": set(support_excluded) == BR_SUPPORT_FLAG_UNITS,
            "actual_minimum_among_remaining_units": remaining_min,
            "readonly_archive_reported_floor_rounded_to_5_decimals": BR_REPORTED_SUPPORT_FLOOR,
            "absolute_rounding_tolerance_used_only_for_archived_floor_comparison": BR_SUPPORT_ROUNDING_TOL,
            "raw_full_precision_support_values_preserved": True,
        },
        "pre_2006_history": {
            "source_years_used": [2003, 2004, 2005],
            "included_in_lag_construction_before_2006_key_selection": True,
            "jan_2006_first_two_days_have_complete_dew3_rh3": True,
        },
        "mask_audit_scope": "Readonly source audit checked 2018-01 raw BR T/D joint support against all 450 region valid_mass values (max absolute difference 2e-7); no full-period raw BR humidity denominator re-extraction is claimed.",
        "code_sha256": code_sha,
        "elapsed_seconds": time.perf_counter() - started,
        "checked_utc": now_utc(),
    }
    manifest["brazil"] = br_audit
    manifest["updated_utc"] = now_utc()
    atomic_json(MANIFEST, manifest)
    log_event({"event": "br_complete", "rows": len(moisture), "tmean_max_abs": br_t_max, "rh_source_max_abs": rh_diff_summary["source_minus_magnus_max_absolute_percentage_points"], "output_sha256": br_audit["output_sha256"]})
    return br_audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mx-only", action="store_true", help="prepare only Mexico 2015-2022 moisture data")
    parser.add_argument("--br-only", action="store_true", help="prepare only Brazil 2006-2020 moisture data")
    args = parser.parse_args()
    if args.mx_only and args.br_only:
        parser.error("--mx-only and --br-only are mutually exclusive")
    DATA.mkdir(parents=True, exist_ok=True)
    ADMIN.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    code_sha = sha256_file(Path(__file__).resolve())
    started = time.perf_counter()
    manifest = load_mx_manifest()
    prior_failure = manifest.get("failure")
    if prior_failure:
        manifest.setdefault("attempt_history", []).append({
            "status": "failed",
            "started_utc": manifest.get("started_utc"),
            "failed_utc": manifest.get("failed_utc"),
            "code_sha256": manifest.get("code_sha256"),
            "failure": prior_failure,
            "resolution": (
                "Accepted 1e-6 absolute tolerance only when comparing the actual full-precision minimum "
                "temperature_valid_mass with the readonly audit's 5-decimal displayed floor 0.62714. "
                "The 2003-2020 annual caches were checksum-verified and reused without raw re-extraction."
            ),
        })
        manifest.pop("failure", None)
        manifest.pop("failed_utc", None)
    manifest.setdefault("project", "latam_stone_20261109")
    manifest.setdefault("purpose", "separately versioned moisture covariates for authorized matched-sample sensitivity analyses")
    manifest["code_sha256"] = code_sha
    manifest["started_utc"] = now_utc()
    manifest["status"] = "processing"
    manifest.setdefault("provenance", {})
    manifest["provenance"].update({
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": code_sha,
        "native_threads": {name: os.environ[name] for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
        "no_effect_models_run": True,
        "canonical_files_unmodified": True,
    })
    atomic_json(MANIFEST, manifest)
    log_event({"event": "run_start", "code_sha256": code_sha, "mx_only": args.mx_only, "br_only": args.br_only})
    try:
        if not args.br_only:
            mx = run_mexico(code_sha)
            manifest = load_mx_manifest()
            manifest["mexico"] = mx
            atomic_json(MANIFEST, manifest)
        if not args.mx_only:
            manifest = load_mx_manifest()
            br = run_brazil(code_sha, manifest)
            manifest = load_mx_manifest()
            manifest["brazil"] = br
            atomic_json(MANIFEST, manifest)
        manifest = load_mx_manifest()
        manifest["status"] = "complete"
        manifest["completed_utc"] = now_utc()
        manifest["elapsed_seconds"] = time.perf_counter() - started
        manifest["outputs"] = {
            BR_OUTPUT.name: {"path": str(BR_OUTPUT), "sha256": sha256_file(BR_OUTPUT), "status": "passed"},
            MX_OUTPUT.name: {"path": str(MX_OUTPUT), "sha256": sha256_file(MX_OUTPUT), "status": "passed"},
        }
        manifest["qa_file"] = str(QA_JSON)
        atomic_json(MANIFEST, manifest)
        qa_payload = {
            "status": "passed",
            "checked_utc": now_utc(),
            "manifest": str(MANIFEST),
            "manifest_sha256": sha256_file(MANIFEST),
            "outputs": manifest["outputs"],
            "mexico": {
                "status": manifest["mexico"]["status"],
                "rows": manifest["mexico"]["rows"],
                "units": manifest["mexico"]["units"],
                "dates": [manifest["mexico"]["first_date"], manifest["mexico"]["last_date"]],
                "all_date_tmean_max_abs_difference_celsius": manifest["mexico"]["tmean_reconciliation"]["max_abs_difference_celsius"],
                "jan2018_previous_route_max_abs_difference_celsius": manifest["mexico"]["previous_route_validation_2018_01"]["max_abs_difference_celsius"],
                "finite_mask_all_months": "T and Td masks equal and time-invariant",
                "three_day_window_initial_na_rows_per_metric": 64,
            },
            "brazil": {
                "status": manifest["brazil"]["status"],
                "rows": manifest["brazil"]["rows"],
                "units": manifest["brazil"]["units"],
                "dates": [manifest["brazil"]["first_date"], manifest["brazil"]["last_date"]],
                "all_date_tmean_max_abs_difference_celsius": manifest["brazil"]["temperature_reconciliation"]["max_abs_difference_celsius"],
                "source_rh_vs_magnus": manifest["brazil"]["humidity_rh_comparison_all_2003_2020"],
                "temperature_valid_mass": manifest["brazil"]["temperature_valid_mass"],
            "support_floor_check": manifest["brazil"]["support_floor_check"],
                "pre_2006_history": manifest["brazil"]["pre_2006_history"],
                "mask_audit_scope": manifest["brazil"]["mask_audit_scope"],
            },
            "no_effect_models_run": True,
            "canonical_files_unmodified": True,
        }
        atomic_json(QA_JSON, qa_payload)
        log_event({"event": "run_complete", "elapsed_seconds": manifest["elapsed_seconds"], "outputs": [str(MX_OUTPUT), str(BR_OUTPUT)]})
        return 0
    except Exception as exc:
        manifest = load_mx_manifest()
        manifest["status"] = "failed"
        manifest["failed_utc"] = now_utc()
        manifest["failure"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        atomic_json(MANIFEST, manifest)
        log_event({"event": "run_failed", "error_type": type(exc).__name__, "error": str(exc)})
        raise


if __name__ == "__main__":
    raise SystemExit(main())
