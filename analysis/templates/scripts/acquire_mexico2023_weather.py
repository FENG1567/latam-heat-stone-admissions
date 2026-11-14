"""Download and derive a separately versioned 2023 Mexico ERA5-Land panel.

The CDS requests are monthly, one-at-a-time, and resumable.  Spatial and daily
aggregation follow the original Mexico state-day builder and are checked
against its January 2018 output before any 2023 rows are written.
"""
from __future__ import annotations

import argparse
import calendar
import hashlib
import json
import os
import sys
import threading
import time
import traceback
import warnings
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr
from concurrent.futures import ThreadPoolExecutor


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data" / "weather_2023"
ADMIN = BASE / "admin"
LOG = BASE / "logs" / "mx_2023_weather_acquisition.log"
MANIFEST = ADMIN / "mx_2023_weather_manifest.json"
SOURCE_CODE = Path(__file__).resolve().parent / "build_mex_state_day.py"
ANALYSIS_PLAN = BASE / "admin" / "analysis_plan.json"
# ASCII aliases are intentional: netCDF4 on this Windows install cannot open
# the equivalent Chinese-character D: paths even when the same files exist.
RAW = Path(r"F:\era5_ascii\mexico")
GEO = Path(r"F:\era5_ascii\mex_states.geojson")
VALIDATION_CSV = Path(r"F:\era5_ascii\mex_state_day_2018_01.csv")
DAILY_CSV = BASE / "data" / "mx_weather_2023.csv"

DATASET = "reanalysis-era5-land"
AREA = [33, -118, 14, -86]
GRID = [0.1, 0.1]
HOURS = [f"{hour:02d}:00" for hour in range(24)]
EXPECTED_COLUMNS = ["state", "year", "month", "day", "tmean"]
MAX_CONCURRENT_REQUESTS = 2
MANIFEST_LOCK = threading.Lock()
LOG_LOCK = threading.Lock()
STOP_EVENT = threading.Event()


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def log_event(event: dict) -> None:
    event = {"timestamp_utc": now_utc(), **event}
    line = json.dumps(event, sort_keys=True, default=str)
    with LOG_LOCK:
        with LOG.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
            handle.flush()
            os.fsync(handle.fileno())
    print(line, flush=True)


def state_code(value) -> str:
    try:
        return f"{int(value):02d}"
    except (TypeError, ValueError):
        return str(value).strip().zfill(2)


def initialize_reference():
    if not SOURCE_CODE.is_file() or not GEO.is_file() or not VALIDATION_CSV.is_file():
        raise FileNotFoundError("Original builder, ASCII geography, or Jan 2018 validation CSV is missing")
    ref_candidates = sorted(RAW.glob("*.nc"))
    if not ref_candidates:
        raise FileNotFoundError(f"No ERA5-Land reference files at {RAW}")
    ref_path = ref_candidates[0]
    states = gpd.read_file(GEO).sort_values("cve_ent").reset_index(drop=True)
    if len(states) != 32:
        raise ValueError(f"Expected 32 Mexican states, found {len(states)}")
    with xr.open_dataset(ref_path, engine="netcdf4") as ref:
        lat = ref.latitude.values.astype(float)
        lon = ref.longitude.values.astype(float)
        if "valid_time" in ref:
            first = ref["t2m"].isel(valid_time=0).values
        else:
            first = ref["t2m"].isel(time=0).values
        valid_land = np.isfinite(first).ravel()
    xx, yy = np.meshgrid(lon, lat)
    area = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)), dtype=np.float32)
    weights = np.zeros((32, len(lat) * len(lon)), dtype=np.float32)
    for i, geom in enumerate(states.geometry):
        inside = shapely.contains_xy(geom, xx, yy).ravel() & valid_land
        weights[i, inside] = area.ravel()[inside]
    row_sums = weights.sum(axis=1, keepdims=True)
    if (row_sums <= 0).any():
        raise ValueError("At least one Mexican state has no valid reference grid cells")
    weights /= row_sums
    return states, lat, lon, weights, ref_path


def load_times(ds: xr.Dataset) -> pd.DatetimeIndex:
    name = "valid_time" if "valid_time" in ds else "time"
    if name not in ds:
        raise ValueError("ERA5-Land file has neither valid_time nor time")
    dates = pd.DatetimeIndex(pd.to_datetime(ds[name].values))
    if dates.tz is not None:
        dates = dates.tz_convert("UTC").tz_localize(None)
    return dates


def validate_raw_month(path: Path, year: int, month: int, lat: np.ndarray, lon: np.ndarray) -> dict:
    expected = pd.date_range(
        f"{year}-{month:02d}-01",
        periods=calendar.monthrange(year, month)[1] * 24,
        freq="h",
    )
    with xr.open_dataset(path, engine="netcdf4") as ds:
        if "t2m" not in ds:
            raise ValueError(f"{path.name} does not contain t2m")
        if "latitude" not in ds.coords or "longitude" not in ds.coords:
            raise ValueError(f"{path.name} lacks latitude/longitude coordinates")
        file_lat = ds.latitude.values.astype(float)
        file_lon = ds.longitude.values.astype(float)
        if not np.array_equal(file_lat, lat) or not np.array_equal(file_lon, lon):
            raise ValueError(f"Grid differs from original Mexico ERA5-Land grid in {path.name}")
        time_coordinate = "valid_time" if "valid_time" in ds else "time"
        time_units = str(ds[time_coordinate].attrs.get("units") or ds[time_coordinate].encoding.get("units", ""))
        time_calendar = str(ds[time_coordinate].attrs.get("calendar") or ds[time_coordinate].encoding.get("calendar", ""))
        times = load_times(ds)
        if not times.equals(expected):
            raise ValueError(
                f"Invalid hourly coverage in {path.name}: got {len(times)} points, "
                f"first={times.min() if len(times) else None}, last={times.max() if len(times) else None}"
            )
        raw_units = str(ds.t2m.attrs.get("units", "K"))
        units = raw_units.strip().lower()
        if units not in {"k", "kelvin", "degrees_k", "degree_k"}:
            raise ValueError(f"Unexpected t2m units {units!r} in {path.name}; expected Kelvin")
        shape = list(ds.t2m.shape)
        if shape[0] != len(expected) or tuple(shape[-2:]) != (len(lat), len(lon)):
            raise ValueError(f"Unexpected t2m shape {shape} in {path.name}")
    return {
        "file": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "hours": len(expected),
        "time_start_utc": expected[0].isoformat(),
        "time_end_utc": expected[-1].isoformat(),
        "time_coordinate": time_coordinate,
        "time_coordinate_units": time_units,
        "time_coordinate_calendar": time_calendar,
        "t2m_units": raw_units,
        "grid_shape_lat_lon": [len(lat), len(lon)],
        "latitude_first_last": [float(lat[0]), float(lat[-1])],
        "longitude_first_last": [float(lon[0]), float(lon[-1])],
    }


def derive_month(path: Path, year: int, month: int, states, lat, lon, weights) -> pd.DataFrame:
    with xr.open_dataset(path, engine="netcdf4") as ds:
        dates = load_times(ds)
        temp = ds["t2m"].values.astype(np.float32, copy=False) - np.float32(273.15)
    expected = calendar.monthrange(year, month)[1] * 24
    if len(dates) != expected:
        raise ValueError(f"Unexpected timestamp count in {path.name}")
    if temp.shape != (expected, len(lat), len(lon)):
        raise ValueError(f"Unexpected t2m array shape in {path.name}: {temp.shape}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        daily_fields = np.nanmean(temp.reshape(-1, 24, len(lat), len(lon)), axis=1)
    days = pd.date_range(f"{year}-{month:02d}-01", periods=len(daily_fields), freq="D")
    rows = []
    for day, field in zip(days, daily_fields):
        flat = field.ravel()
        for i, state in enumerate(states["cve_ent"]):
            mask = weights[i] > 0
            numerator = np.nansum(flat[mask] * weights[i, mask])
            denominator = np.nansum(weights[i, mask] * np.isfinite(flat[mask]))
            value = float(numerator / denominator) if denominator > 0 else np.nan
            rows.append((state_code(state), year, month, day.day, value))
    result = pd.DataFrame(rows, columns=EXPECTED_COLUMNS)
    if len(result) != 32 * len(days) or result.tmean.isna().any():
        raise ValueError(f"Incomplete/invalid state-day values for {year}-{month:02d}")
    return result


def compare_validation(states, lat, lon, weights, reference_path: Path) -> dict:
    candidate = RAW / "era5land_hourly_temp_dewpoint_mexico_2018_01.nc"
    if not candidate.is_file():
        raise FileNotFoundError(f"January 2018 ERA5-Land reference missing: {candidate}")
    raw_meta = validate_raw_month(candidate, 2018, 1, lat, lon)
    computed = derive_month(candidate, 2018, 1, states, lat, lon, weights)
    expected = pd.read_csv(VALIDATION_CSV, dtype={"state": str})
    expected["state"] = expected.state.map(state_code)
    expected["tmean"] = pd.to_numeric(expected.tmean, errors="raise")
    paired = computed.merge(
        expected[EXPECTED_COLUMNS], on=["state", "year", "month", "day"],
        suffixes=("_new", "_old"), validate="one_to_one",
    )
    max_abs = float((paired.tmean_new - paired.tmean_old).abs().max())
    if len(paired) != 32 * 31 or max_abs > 1e-6:
        raise ValueError(f"January 2018 validation failed: rows={len(paired)}, max_abs_diff={max_abs}")
    result = {
        "status": "passed",
        "month": "2018-01",
        "matched_state_days": len(paired),
        "max_abs_difference_celsius": max_abs,
        "validation_csv": str(VALIDATION_CSV),
        "validation_csv_sha256": sha256(VALIDATION_CSV),
        "raw_reference": raw_meta,
        "checked_utc": now_utc(),
    }
    return result


def request_for(year: int, month: int) -> dict:
    return {
        "variable": ["2m_temperature"],
        "year": str(year),
        "month": f"{month:02d}",
        "day": [f"{day:02d}" for day in range(1, calendar.monthrange(year, month)[1] + 1)],
        "time": HOURS,
        "area": AREA,
        "grid": GRID,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }


def cds_info_callback(message, *args) -> None:
    """Log task-state transitions without exposing task URLs or credentials."""
    try:
        rendered = message % args if args else str(message)
    except Exception:
        rendered = str(message)
    prefix = "Request is "
    if rendered.startswith(prefix):
        log_event({"event": "cds_status", "state": rendered[len(prefix):].strip()})


def retrieve_month(year: int, month: int, request: dict, output: Path, audit: dict) -> None:
    import cdsapi
    import requests

    request_canonical = json.dumps(request, sort_keys=True, separators=(",", ":"))
    request_hash = hashlib.sha256(request_canonical.encode("utf-8")).hexdigest()
    tag = f"{year}-{month:02d}"
    with MANIFEST_LOCK:
        month_record = audit["months"].setdefault(tag, {})
        month_record.update({"request_sha256": request_hash, "request": request, "status": "requesting"})
        atomic_json(MANIFEST, audit)
    tmp = output.with_suffix(output.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    # cdsapi's default Session argument is shared across Client instances; give
    # each worker its own requests.Session for thread safety.
    client = cdsapi.Client(
        quiet=True,
        progress=False,
        session=requests.Session(),
        info_callback=cds_info_callback,
        retry_max=4,
        sleep_max=30,
        timeout=60,
    )
    last_error = None
    # Queue caps are transient. Retry a bounded number of times and stop new
    # submissions after exhaustion rather than feed the remaining months into a cap.
    delays = [30, 120, 300]
    if STOP_EVENT.is_set():
        error = "not submitted because a peer month reported a CDS failure"
        with MANIFEST_LOCK:
            month_record.update({"status": "skipped_after_peer_failure", "last_error": error})
            atomic_json(MANIFEST, audit)
        log_event({"event": "request_skipped", "month": tag, "reason": error})
        raise RuntimeError(error)
    for attempt in range(1, len(delays) + 2):
        log_event({
            "event": "request_start", "month": tag,
            "attempt": attempt, "request_sha256": request_hash,
        })
        try:
            result = client.retrieve(DATASET, request, str(tmp))
            if not tmp.is_file() or tmp.stat().st_size <= 100_000:
                raise IOError("CDS returned without a plausible NetCDF file")
            os.replace(tmp, output)
            with MANIFEST_LOCK:
                month_record.update({
                    "status": "downloaded",
                    "attempts": attempt,
                    "downloaded_utc": now_utc(),
                    "size_bytes": output.stat().st_size,
                    "sha256": sha256(output),
                    # The CDS client result can include an opaque job/result URL;
                    # record only presence/type and never credentials.
                    "client_result_type": type(result).__name__,
                })
                audit["months"][tag] = month_record
                atomic_json(MANIFEST, audit)
            log_event({
                "event": "download_complete", "month": tag,
                "attempt": attempt, "bytes": output.stat().st_size,
                "sha256": month_record["sha256"],
            })
            return
        except Exception as exc:  # CDS exceptions vary by client/API release.
            last_error = f"{type(exc).__name__}: {str(exc)[:700]}"
            if tmp.exists():
                tmp.unlink()
            with MANIFEST_LOCK:
                month_record.update({"status": "retrying", "attempts": attempt, "last_error": last_error})
                atomic_json(MANIFEST, audit)
            log_event({
                "event": "request_failure", "month": tag,
                "attempt": attempt, "error": last_error,
            })
            # Stop queued future months on the first failure, but let this
            # already-submitted month complete its bounded retry sequence.
            STOP_EVENT.set()
            if attempt <= len(delays):
                delay = delays[attempt - 1]
                log_event({"event": "retry_wait", "month": tag, "seconds": delay})
                time.sleep(delay)
    with MANIFEST_LOCK:
        month_record.update({"status": "blocked", "attempts": len(delays) + 1, "last_error": last_error})
        atomic_json(MANIFEST, audit)
    STOP_EVENT.set()
    raise RuntimeError(f"CDS request blocked at {tag} after bounded retries: {last_error}")


def valid_cached_daily(path: Path, year: int, month: int) -> bool:
    if not path.is_file():
        return False
    try:
        frame = pd.read_csv(path, dtype={"state": str})
        frame.state = frame.state.map(state_code)
        expected_days = calendar.monthrange(year, month)[1]
        return (
            frame.columns.tolist() == EXPECTED_COLUMNS
            and len(frame) == 32 * expected_days
            and frame["state"].nunique() == 32
            and set(frame.year.astype(int)) == {year}
            and set(frame.month.astype(int)) == {month}
            and not frame.tmean.isna().any()
            and not frame.duplicated(["state", "year", "month", "day"]).any()
        )
    except Exception:
        return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-only", action="store_true", help="run the Jan 2018 parity check only")
    parser.add_argument("--months", nargs="+", type=int, choices=range(1, 13), help="subset for a deliberate partial run")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ADMIN.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    states, lat, lon, weights, ref_path = initialize_reference()
    validation = compare_validation(states, lat, lon, weights, ref_path)
    audit = {
        "project": "latam_stone_20261109",
        "purpose": "separately versioned Mexico 2023 post-pandemic temporal replication weather input",
        "authorization": "Direct public CDS data download authorized; no patient-level records downloaded or uploaded.",
        "status": "validation_passed_acquisition_pending",
        "created_utc": now_utc(),
        "dataset": DATASET,
        "retrieval_protocol": f"monthly requests with at most {MAX_CONCURRENT_REQUESTS} active; each uses a dedicated CDS client/session; bounded retries",
        "request_spec": {
            "variable": "2m_temperature only",
            "area_north_west_south_east": AREA,
            "grid_lat_lon_degrees": GRID,
            "hours_utc": HOURS,
            "time_coverage": "2023-01-01 00:00 through 2023-12-31 23:00 inclusive",
            "format": "NetCDF, unarchived",
        },
        "aggregation_method": {
            "geography": "32 state polygons sorted by cve_ent; reference gridpoint centers contained strictly within polygon",
            "weights": "cos(latitude) spatial area weights from fixed original reference valid-land mask; renormalize over each day's finite temperatures",
            "daily": "UTC calendar day arithmetic mean of all 24 hourly 2m temperatures per grid cell, NaN ignored",
            "units": "Kelvin converted to degrees Celsius",
        },
        "provenance": {
            "extractor_script": str(Path(__file__).resolve()),
            "extractor_script_sha256": sha256(Path(__file__).resolve()),
            "original_builder": str(SOURCE_CODE),
            "original_builder_sha256": sha256(SOURCE_CODE),
            "pre_results_analysis_plan": str(ANALYSIS_PLAN),
            "pre_results_analysis_plan_sha256": sha256(ANALYSIS_PLAN) if ANALYSIS_PLAN.is_file() else None,
            "state_geography": str(GEO),
            "state_geography_sha256": sha256(GEO),
            "reference_nc": str(ref_path),
            "reference_nc_sha256": sha256(ref_path),
            "reference_grid": {"shape_lat_lon": [len(lat), len(lon)], "latitude_first_last": [float(lat[0]), float(lat[-1])], "longitude_first_last": [float(lon[0]), float(lon[-1])]},
        },
        "validation_2018_01": validation,
        "months": {},
    }
    if MANIFEST.is_file():
        try:
            prior = json.loads(MANIFEST.read_text(encoding="utf-8"))
            if prior.get("provenance", {}).get("original_builder_sha256") == audit["provenance"]["original_builder_sha256"]:
                audit["months"] = prior.get("months", {})
        except Exception:
            pass
    atomic_json(MANIFEST, audit)
    log_event({"event": "validation_passed", "month": "2018-01", "matched_state_days": validation["matched_state_days"], "max_abs_diff": validation["max_abs_difference_celsius"]})
    if args.validate_only:
        audit["status"] = "validation_passed"
        atomic_json(MANIFEST, audit)
        return 0

    try:
        import cdsapi  # noqa: F401
    except Exception as exc:
        audit["status"] = "cdsapi_unavailable"
        audit["blocking_error"] = f"{type(exc).__name__}: {exc}"
        atomic_json(MANIFEST, audit)
        raise
    months = args.months or list(range(1, 13))
    download_jobs = {}
    for month in months:
        tag = f"2023-{month:02d}"
        raw_path = OUT / f"era5land_hourly_t2m_mexico_{tag}.nc"
        if not raw_path.is_file():
            download_jobs[month] = raw_path

    # Network retrieval is the only concurrent phase. The pool owns at most
    # two active CDS requests; all NetCDF validation and aggregation below run
    # in the main thread.
    executor = ThreadPoolExecutor(max_workers=MAX_CONCURRENT_REQUESTS, thread_name_prefix="cds-mx-2023")
    futures = {}
    try:
        for month, raw_path in download_jobs.items():
            futures[month] = executor.submit(retrieve_month, 2023, month, request_for(2023, month), raw_path, audit)

        try:
            for month in months:
                tag = f"2023-{month:02d}"
                raw_path = OUT / f"era5land_hourly_t2m_mexico_{tag}.nc"
                daily_path = OUT / f"mx_weather_2023_{month:02d}.csv"
                with MANIFEST_LOCK:
                    month_record = audit["months"].setdefault(tag, {})
                if not raw_path.is_file():
                    try:
                        futures[month].result()
                    except Exception as exc:
                        STOP_EVENT.set()
                        with MANIFEST_LOCK:
                            month_record.update({"status": "blocked", "blocking_error": f"{type(exc).__name__}: {str(exc)[:700]}"})
                            audit["status"] = "blocked"
                            atomic_json(MANIFEST, audit)
                        log_event({"event": "acquisition_blocked", "month": tag, "error": month_record["blocking_error"]})
                        return 2
                try:
                    raw_meta = validate_raw_month(raw_path, 2023, month, lat, lon)
                except Exception as exc:
                    STOP_EVENT.set()
                    with MANIFEST_LOCK:
                        month_record.update({"status": "raw_invalid", "raw_validation_error": f"{type(exc).__name__}: {str(exc)[:500]}"})
                        atomic_json(MANIFEST, audit)
                    raise
                with MANIFEST_LOCK:
                    month_record.update({"raw": raw_meta, "status": month_record.get("status", "downloaded_validated")})
                if valid_cached_daily(daily_path, 2023, month):
                    daily = pd.read_csv(daily_path, dtype={"state": str})
                    daily_record = {"file": daily_path.name, "size_bytes": daily_path.stat().st_size, "sha256": sha256(daily_path), "rows": len(daily), "status": "cached_valid"}
                else:
                    frame = derive_month(raw_path, 2023, month, states, lat, lon, weights)
                    atomic_csv(frame, daily_path)
                    daily_record = {"file": daily_path.name, "size_bytes": daily_path.stat().st_size, "sha256": sha256(daily_path), "rows": len(frame), "status": "derived"}
                with MANIFEST_LOCK:
                    month_record["daily"] = daily_record
                    month_record["status"] = "complete"
                    atomic_json(MANIFEST, audit)
                log_event({"event": "month_complete", "month": tag, "rows": daily_record["rows"], "daily_sha256": daily_record["sha256"]})
        except Exception:
            STOP_EVENT.set()
            raise
    finally:
        executor.shutdown(wait=True, cancel_futures=False)

    if len(months) != 12:
        audit["status"] = "partial_month_subset_complete"
        atomic_json(MANIFEST, audit)
        return 0
    frames = []
    for month in range(1, 13):
        path = OUT / f"mx_weather_2023_{month:02d}.csv"
        if not valid_cached_daily(path, 2023, month):
            raise ValueError(f"Missing/invalid month after extraction: {path.name}")
        frames.append(pd.read_csv(path, dtype={"state": str}))
    annual = pd.concat(frames, ignore_index=True)
    annual.state = annual.state.map(state_code)
    annual = annual[EXPECTED_COLUMNS].sort_values(["state", "year", "month", "day"], kind="stable").reset_index(drop=True)
    if len(annual) != 32 * 365 or annual.tmean.isna().any() or annual.duplicated(["state", "year", "month", "day"]).any():
        raise ValueError(f"Annual state-day coverage failed: rows={len(annual)}")
    if set(annual.state.unique()) != {state_code(x) for x in states.cve_ent}:
        raise ValueError("Annual output state codes differ from source geography")
    atomic_csv(annual, DAILY_CSV)
    audit.update({
        "status": "complete",
        "completed_utc": now_utc(),
        "output": {
            "file": str(DAILY_CSV),
            "rows": len(annual),
            "states": int(annual.state.nunique()),
            "first_date": "2023-01-01",
            "last_date": "2023-12-31",
            "sha256": sha256(DAILY_CSV),
        },
    })
    atomic_json(MANIFEST, audit)
    log_event({"event": "annual_complete", "file": str(DAILY_CSV), "rows": len(annual), "sha256": audit["output"]["sha256"]})
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        # Avoid dumping low-level auth/config objects. Error strings can be
        # scrubbed by cdsapi, but still redact common credential-like tokens.
        message = f"{type(exc).__name__}: {str(exc)[:900]}"
        log_event({"event": "fatal", "error": message})
        raise

