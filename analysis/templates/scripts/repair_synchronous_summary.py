"""Repair six HAC90 summary rows from saved national daily influence series.

No model fitting or fit-cache reading is performed. Original outputs and run
manifests retain their original provenance in the invalid-export archive.
"""
from pathlib import Path
import csv
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "results/synchronous_meta"
ARCHIVE = RESULTS / "invalid_pre_hac90_label_fix"
ARCHIVE.mkdir(parents=True, exist_ok=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf8", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def preserve(path):
    target = ARCHIVE / path.relative_to(RESULTS)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise RuntimeError(f"Refusing to replace previous repair archive {target}")
    shutil.copy2(path, target)
    return sha(target)


paths = [RESULTS / "BR_primary_threeoutcome_E3_separate_summary.csv",
         RESULTS / "MX_primary_threeoutcome_E3_separate_summary.csv",
         RESULTS / "mx_updated_primary/MX_primary_complete_source_2015_2022_summary.csv"]
records = []
for path in paths:
    fields, rows = read_csv(path)
    original_hash = preserve(path)
    daily_path = path.with_name(path.name.replace("_summary.csv", "_country_daily_influence.csv"))
    d = pd.read_csv(daily_path, parse_dates=["date"])
    if d.date.duplicated().any():
        raise ValueError("Duplicate dates in national daily influence")
    dates = pd.date_range(d.date.min(), d.date.max())
    psi = d.set_index("date").weighted_country_influence.reindex(dates, fill_value=0).to_numpy()
    variance = float(psi @ psi)
    for lag in range(1, 91):
        variance += 2 * (1 - lag / 91) * float(psi[lag:] @ psi[:-lag])
    changed = 0
    for row in rows:
        if int(row["synchronized_Bartlett_bandwidth_days"]) != 90:
            continue
        total = variance + float(row["V_tau_component"])
        se = np.sqrt(total)
        mean = float(row["estimate_logRR"])
        halfwidth = norm.ppf(0.975) * se
        updates = {"V_sync": variance, "variance": total, "se": se,
                   "lo_logRR": mean - halfwidth, "hi_logRR": mean + halfwidth,
                   "lo_RR": np.exp(mean - halfwidth), "hi_RR": np.exp(mean + halfwidth)}
        row.update({key: repr(float(value)) for key, value in updates.items()})
        changed += 1
    if changed != 2:
        raise ValueError(f"Expected two HAC90 rows in {path}")
    write_csv(path, fields, rows)
    records.append({"path": str(path.relative_to(BASE)), "old_sha256": original_hash,
                    "new_sha256": sha(path), "corrected_rows": changed,
                    "daily_influence_sha256": sha(daily_path), "true_HAC90_variance": variance})

for directory, members in [(RESULTS, paths[:2]), (RESULTS / "mx_updated_primary", paths[2:])]:
    combined = directory / "synchronous_meta_summary.csv"
    original_hash = preserve(combined)
    all_rows = []
    for path in members:
        fields, rows = read_csv(path)
        all_rows.extend(rows)
    write_csv(combined, fields, all_rows)
    records.append({"path": str(combined.relative_to(BASE)), "old_sha256": original_hash,
                    "new_sha256": sha(combined), "rebuilt_combined_export": True})

audit = RESULTS / "synchronous_meta_independent_audit.csv"
if audit.exists():
    preserve(audit)
manifest = {"reason": "Original summary loop reused HAC28 variance for HAC90 label",
            "procedure": "Recompute actual-calendar Bartlett HAC90 from saved national daily influences",
            "model_refits": 0, "RDS_caches_reread": 0,
            "only_HAC90_numeric_fields_changed": True,
            "source_run_manifests_preserved": True,
            "original_script_sha256": sha(ARCHIVE / "run_synchronous_meta.R"),
            "corrected_script_sha256": sha(BASE / "scripts/run_synchronous_meta.R"),
            "artifacts": records}
(RESULTS / "hac90_repair_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf8")
print(json.dumps({"repaired_country_summaries": 3, "corrected_rows": 6,
                  "combined_exports_rebuilt": 2, "model_refits": 0}, indent=2))
