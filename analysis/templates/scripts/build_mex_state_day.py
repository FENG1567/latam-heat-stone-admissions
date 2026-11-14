"""Compute Mexico state-day mean temperature from ERA5-Land hourly nc files.

Usage:
  python build_mex_state_day.py --years 2018 --months 1 --out test.csv   # validation
  python build_mex_state_day.py --years 2015 2016 2017 2023 --out new.csv
"""
import argparse
import glob
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

BASE = r"F:\era5_ascii\mexico"
STATES = r"F:\era5_ascii\mex_states.geojson"


def grid_weights(states, lat, lon, valid=None):
    from shapely.prepared import prep
    ny, nx = len(lat), len(lon)
    W = np.zeros((len(states), ny * nx), dtype=np.float32)
    area = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nx), dtype=np.float32)
    for i, geom in enumerate(states.geometry):
        pg = prep(geom)
        for iy in range(ny):
            for ix in range(nx):
                if pg.contains(__import__("shapely").geometry.Point(lon[ix], lat[iy])):
                    W[i, iy * nx + ix] = area[iy, ix] if (valid is None or valid[iy * nx + ix]) else 0.0
    s = W.sum(axis=1, keepdims=True)
    s[s == 0] = 1
    return W / s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", nargs="+", type=int, required=True)
    ap.add_argument("--months", nargs="+", type=int, default=list(range(1, 13)))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    states = gpd.read_file(STATES).sort_values("cve_ent").reset_index(drop=True)
    ref = xr.open_dataset(sorted(glob.glob(os.path.join(BASE, "*.nc")))[0])
    lat = ref.latitude.values.astype(float)
    lon = ref.longitude.values.astype(float)
    valid = ~np.isnan(ref["t2m"].values[0]).ravel()
    W = grid_weights(states, lat, lon, valid)
    print("weights built", W.shape, "states", len(states), flush=True)

    rows = []
    for y in args.years:
        for m in args.months:
            p = os.path.join(BASE, f"era5land_hourly_temp_dewpoint_mexico_{y}_{m:02d}.nc")
            if not os.path.exists(p):
                print("MISSING", p, flush=True)
                continue
            ds = xr.open_dataset(p, engine="netcdf4")
            t = ds["t2m"].values - 273.15  # (time, lat, lon)
            times = pd.to_datetime(ds["valid_time"].values)
            ds.close()
            df = pd.DataFrame({"date": times.normalize(), "idx": np.arange(len(times))})
            day_mean = df.groupby("date")["idx"].apply(lambda ix: np.nanmean(t[ix.values], axis=0))
            for d, field in day_mean.items():
                flat = field.ravel()
                for i, st in enumerate(states["cve_ent"].values):
                    w = W[i]
                    msk = w > 0
                    num = np.nansum(flat[msk] * w[msk])
                    den = np.nansum(w[msk] * (~np.isnan(flat[msk])))
                    val = num / den if den > 0 else np.nan
                    rows.append((st, d.year, d.month, d.day, float(val)))
            print(f"{y}-{m:02d} days={len(day_mean)}", flush=True)
    out = pd.DataFrame(rows, columns=["state", "year", "month", "day", "tmean"])
    out.to_csv(args.out, index=False)
    print("wrote", args.out, out.shape, flush=True)


if __name__ == "__main__":
    main()



