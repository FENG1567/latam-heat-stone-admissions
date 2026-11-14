"""Audit fixed temperature targets against each independent refit period.

Uses existing estimates without refitting. Support requires each target within
the period's P2–P98 and at least 30 eligible days within one degree Celsius.
Reports the prespecified-style support rule as an additional sensitivity audit.
"""
from pathlib import Path
import json
import re

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / "temporal_support"
OUT.mkdir(parents=True, exist_ok=True)


def dl_pool(frame):
    y = frame.estimate.to_numpy()
    v = frame.variance.to_numpy()
    if len(y) < 2:
        return {"k": len(y), "status": "insufficient_supported_units"}
    if frame.unit.duplicated().any():
        raise ValueError("Duplicate units in temporal pool")
    w = 1 / v
    mean = np.sum(w * y) / np.sum(w)
    q = np.sum(w * (y - mean) ** 2)
    c = np.sum(w) - np.sum(w * w) / np.sum(w)
    tau = max(0.0, (q - len(y) + 1) / c)
    w = 1 / (v + tau)
    mean = np.sum(w * y) / np.sum(w)
    se = np.sqrt(1 / np.sum(w))
    z = norm.ppf(0.975)
    return {"k": len(y), "log_estimate": mean, "log_se": se,
            "RRR": np.exp(mean), "lo": np.exp(mean - z * se),
            "hi": np.exp(mean + z * se), "tau2": tau,
            "I2": 100 * max(0.0, (q - len(y) + 1) / q) if q else 0.0,
            "concordance": np.mean(y > 0), "status": "ok"}


def read_calendar(name, country):
    cols = ["unit", "date", "E3", "y_K80", "y_urinary"]
    if country == "BR":
        cols.append("population")
    d = pd.read_csv(BASE / "data" / name, usecols=cols, dtype={"unit": str})
    eligible = np.isfinite(d.E3) & np.isfinite(d.y_K80) & np.isfinite(d.y_urinary)
    eligible &= (d.y_K80 >= 0) & (d.y_urinary >= 0)
    if country == "BR":
        eligible &= np.isfinite(d.population) & (d.population > 0)
    return d.loc[eligible]


def scenario_range(path):
    match = re.search(r"temporal_(\d{8})_(\d{8})", path.name)
    if match:
        return tuple(pd.Timestamp(x).strftime("%Y-%m-%d") for x in match.groups())
    if "2015_2018" in path.name:
        return "2015-01-01", "2018-12-31"
    if "2020_2022" in path.name:
        return "2020-01-01", "2022-12-31"
    if "2023_jan_nov" in path.name:
        return "2023-01-01", "2023-11-30"
    if "2023" in path.name:
        return "2023-01-01", "2023-12-31"
    raise ValueError(f"Unrecognized temporal period {path}")


def audit_scenario(path, calendar, country):
    estimates = pd.read_csv(path, dtype={"unit": str})
    estimates = estimates.loc[(estimates.contrast == "y_urinary_vs_y_K80") &
                              (estimates.covariance == "hac14")].copy()
    if estimates.unit.duplicated().any():
        raise ValueError(f"Duplicate units in {path}")
    start, end = scenario_range(path)
    period = calendar.loc[calendar.date.between(start, end)]
    grouped = dict(tuple(period.groupby("unit")))
    audit = []
    for row in estimates.itertuples():
        e = grouped[row.unit].E3.to_numpy()
        if len(e) != row.n_complete:
            raise ValueError(f'Eligible-day count differs from fitted sample for {row.unit} in {row.scenario}')
        p2, p98 = np.quantile(e, [0.02, 0.98])
        n_high = np.count_nonzero(np.abs(e - row.p95) <= 1)
        n_ref = np.count_nonzero(np.abs(e - row.median) <= 1)
        supported = p2 <= row.median <= p98 and p2 <= row.p95 <= p98
        supported = supported and n_high >= 30 and n_ref >= 30
        audit.append({"country": country, "scenario": row.scenario, "unit": row.unit,
                      "period_start": start, "period_end": end, "eligible_days": len(e),
                      "target_P95": row.p95, "target_median": row.median,
                      "period_P2": p2, "period_P98": p98,
                      "days_near_P95": n_high, "days_near_median": n_ref,
                      "supported": bool(supported)})
    checked = pd.DataFrame(audit)
    supported_units = set(checked.loc[checked.supported, "unit"])
    return estimates, checked, supported_units


records = []
pools = []
families = [
    ("BR_original", "BR", "br_calendar.csv", "core", "*temporal*fixed_full_ref_BR_per_unit.csv"),
    ("MX_original_gap", "MX", "mx_calendar.csv.gz", "core", "*temporal*fixed_full_ref_MX_per_unit.csv"),
    ("MX_complete", "MX", "mx_updated_calendar.csv.gz", "mx_updated", "temporal_*fixed_full_reference_MX_per_unit.csv"),
    ("MX_2023", "MX", "mx_updated_2023_calendar.csv.gz", "mx_updated", "*2023*fixed*reference_MX_per_unit.csv"),
]
for family, country, calendar_name, folder, pattern in families:
    calendar = read_calendar(calendar_name, country)
    scenarios = []
    for path in sorted((BASE / "results" / folder).glob(pattern)):
        estimates, checked, supported = audit_scenario(path, calendar, country)
        checked["family"] = family
        records.append(checked)
        scenarios.append((estimates, supported))
    if not scenarios:
        raise ValueError(f"No temporal scenarios found for {family}")
    common = set.intersection(*(units for _, units in scenarios))
    for estimates, supported in scenarios:
        for scope, units in [("all_fitted", set(estimates.unit)),
                             ("period_supported", supported),
                             ("common_supported_across_family", common)]:
            pooled = dl_pool(estimates.loc[estimates.unit.isin(units)])
            pools.append({"country": country, "family": family,
                          "scenario": estimates.scenario.iloc[0], "scope": scope,
                          "covariance": "hac14", "method": "DL", **pooled})
    print(f"{family}: {len(scenarios)} scenarios, {len(common)} common-supported units", flush=True)

pd.concat(records, ignore_index=True).to_csv(OUT / "per_unit_support.csv", index=False)
pd.DataFrame(pools).to_csv(OUT / "support_restricted_pools.csv", index=False)
(OUT / "audit.json").write_text(json.dumps({
    "rule": "Both fixed targets inside period P2–P98 with at least 30 eligible days within 1 C of each target",
    "model_refits_performed": False, "source": "Existing per-unit HAC14 estimates",
    "period_unit_rows": sum(len(x) for x in records), "pools": len(pools),
    "common_support_grouping": "Within calendar family, separately for original gaps, complete cohort, and 2023",
    "decision_timing": "Post-result support audit, not a prospective protocol amendment"
}, indent=2), encoding="utf8")
