#!/usr/bin/env python3
"""Select and audit existing Mexico follow-up fit outputs; no refitting."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
INPUTS = HERE / "inputs"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path.name}")
    if fields is None:
        fields = list(rows[0])
    if any(name in {"source_file", "source_row"} for name in fields):
        raise ValueError("reserved source_file/source_row column name")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def verify_result_manifest(name: str) -> dict:
    manifest = read_json(RESULTS / name)
    if manifest.get("status") != "complete":
        raise ValueError(f"{name} is not complete")
    for rel, expected in manifest.get("output_sha256", {}).items():
        path = RESULTS / rel
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"manifest hash mismatch or missing output: {rel}")
    return manifest


def numerics(row: dict, keys: list[str]) -> dict[str, float]:
    return {key: float(row[key]) for key in keys}


def select_unique(rows: list[dict[str, str]], **criteria: str) -> dict[str, str]:
    match = [row for row in rows if all(row.get(k) == v for k, v in criteria.items())]
    if len(match) != 1:
        raise ValueError(f"expected one row for {criteria}; found {len(match)}")
    return match[0]


def pooled_rows(path: Path) -> list[dict[str, str]]:
    rows = read_csv(path)
    return [r for r in rows if r.get("contrast") == "y_urinary_vs_y_K80"
            and r.get("covariance") in {"hac14", "hac28"}
            and r.get("method") in {"DL", "REML_HK"}]


def per_unit_rows(path: Path, covariance: str) -> list[dict[str, str]]:
    return [r for r in read_csv(path) if r.get("contrast") == "y_urinary_vs_y_K80"
            and r.get("covariance") == covariance]


def same_sign_fraction(sensitivity: list[dict[str, str]], reference: list[dict[str, str]]) -> tuple[str, str]:
    ref = {r["unit"]: float(r["estimate"]) for r in reference}
    paired = [(float(r["estimate"]), ref[r["unit"]]) for r in sensitivity if r["unit"] in ref]
    if not paired:
        return "", "0"
    same = sum((a > 0) == (b > 0) if a != 0 and b != 0 else a == b for a, b in paired)
    return f"{100.0 * same / len(paired):.6f}", str(len(paired))


def attempts_summary(attempt_path: Path, per_path: Path, failure_path: Path | None = None) -> dict[str, int]:
    attempts = read_csv(attempt_path) if attempt_path and attempt_path.is_file() else []
    fitted = {r["unit"] for r in read_csv(per_path)}
    if attempts:
        attempted = len({r["unit"] for r in attempts})
        eligible = sum(r.get("eligible", "").lower() == "true" for r in attempts)
    else:
        attempted = 32
        eligible = len(fitted)
    fit_failures = 0
    eligibility_failures = 0
    if failure_path and failure_path.is_file():
        failures = read_csv(failure_path)
        fit_failures = sum(r.get("phase") == "fit" for r in failures)
        eligibility_failures = sum(r.get("phase") == "eligibility" for r in failures)
    return {"states_attempted": attempted, "states_eligible": eligible,
            "states_fit": len(fitted), "states_ineligible": eligibility_failures,
            "fit_failures": fit_failures}


def main() -> None:
    a_manifest = verify_result_manifest("a_manifest.json")
    b_manifest = verify_result_manifest("b_manifest.json")

    full_events: dict[str, dict[str, float]] = {}
    for row in read_csv(RESULTS / "a_cohort_event_counts.csv"):
        full_events.setdefault(row["cohort"], {})[row["outcome"]] = float(row["events"])
    for row in read_csv(RESULTS / "b_candidate_repeat_subtraction_counts.csv"):
        if row["sensitivity"] == "primary_2015_2022_nonID_repeat_upper_bound":
            full_events.setdefault("primary_registered", {})[row["outcome"]] = float(row["original_events"])
            full_events.setdefault("primary_nonID_repeat_bound", {})[row["outcome"]] = float(row["bounded_events"])
        elif row["sensitivity"] == "extended_2023_nonID_repeat_upper_bound":
            full_events.setdefault("extended_nonID_repeat_bound", {})[row["outcome"]] = float(row["bounded_events"])
    for row in read_csv(RESULTS / "b_sourceyear_concordance_event_counts.csv"):
        if row["sensitivity"] == "primary_2015_2022_source_year_equals_discharge_year":
            full_events.setdefault("primary_registered", {})[row["outcome"]] = float(row["registered_events"])
            full_events.setdefault("primary_sourceyear_concordant", {})[row["outcome"]] = float(row["concordant_events"])
        elif row["sensitivity"] == "extended_2023_source_year_equals_discharge_year":
            full_events.setdefault("extended_2023_2024", {})[row["outcome"]] = float(row["registered_events"])
            full_events.setdefault("extended_sourceyear_concordant", {})[row["outcome"]] = float(row["concordant_events"])

    scenario_event_key = {
        "mexico_primary_registered": "primary_registered",
        "mexico_primary_nonID_repeat_upper_bound": "primary_nonID_repeat_bound",
        "mexico_primary_sourceyear_concordant": "primary_sourceyear_concordant",
        "mexico2023_source2023_baseline": "source2023_baseline",
        "mexico2023_followup": "source2023_plus_2024_extended_2023",
        "mexico2023_followup_nonID_repeat_upper_bound": "extended_nonID_repeat_bound",
        "mexico2023_followup_sourceyear_concordant": "extended_sourceyear_concordant",
    }

    scenarios = [
        {
            "scenario_tag": "mexico_primary_registered",
            "cohort_label": "Registered complete-source Mexico; admission dates 2015-2022",
            "source_selection": "Existing registered cohort from annual source files 2015-2023",
            "per": INPUTS / "primary_complete_source_2015_2022_MX_per_unit.csv",
            "pool": INPUTS / "primary_complete_source_2015_2022_MX_pooled.csv",
            "reference_per": None, "reference_label": "registered primary; baseline",
            "attempt": None, "failure": None,
        },
        {
            "scenario_tag": "mexico_primary_nonID_repeat_upper_bound",
            "cohort_label": "Registered complete-source Mexico; admission dates 2015-2022",
            "source_selection": "Subtract exact non-ID row-repeat candidate extras from source years 2015-2023; conservative upper bound",
            "per": RESULTS / "b_primary_repeat_bound_per_unit.csv",
            "pool": RESULTS / "b_primary_repeat_bound_pooled.csv",
            "reference_per": INPUTS / "primary_complete_source_2015_2022_MX_per_unit.csv",
            "reference_label": "registered primary",
            "attempt": RESULTS / "b_primary_repeat_bound_attempts.csv",
            "failure": RESULTS / "b_primary_repeat_bound_failure_audit.csv",
        },
        {
            "scenario_tag": "mexico_primary_sourceyear_concordant",
            "cohort_label": "Registered complete-source Mexico; admission dates 2015-2022",
            "source_selection": "Retain only rows where source_year equals discharge_year; source-selection sensitivity",
            "per": RESULTS / "b_primary_sourceyear_concordant_per_unit.csv",
            "pool": RESULTS / "b_primary_sourceyear_concordant_pooled.csv",
            "reference_per": INPUTS / "primary_complete_source_2015_2022_MX_per_unit.csv",
            "reference_label": "registered primary",
            "attempt": RESULTS / "b_primary_sourceyear_concordant_attempts.csv",
            "failure": RESULTS / "b_primary_sourceyear_concordant_failure_audit.csv",
        },
        {
            "scenario_tag": "mexico2023_source2023_baseline",
            "cohort_label": "Mexico admission dates in 2023; discharge/source year 2023 only",
            "source_selection": "Rebuilt from parser-verified source-year 2023 daily aggregate",
            "per": RESULTS / "a_baseline_2023_rebuilt_per_unit.csv",
            "pool": RESULTS / "a_baseline_2023_rebuilt_pooled.csv",
            "reference_per": None, "reference_label": "2023 source-only baseline",
            "attempt": RESULTS / "a_baseline_2023_rebuilt_attempts.csv",
            "failure": RESULTS / "a_baseline_2023_rebuilt_failure_audit.csv",
        },
        {
            "scenario_tag": "mexico2023_followup",
            "cohort_label": "Mexico admission dates in 2023; discharge/source years 2023 and 2024",
            "source_selection": "Add source-year 2024 discharges admitted in 2023",
            "per": RESULTS / "a_extended_2023_2024_per_unit.csv",
            "pool": RESULTS / "a_extended_2023_2024_pooled.csv",
            "reference_per": RESULTS / "a_baseline_2023_rebuilt_per_unit.csv",
            "reference_label": "2023 source-only baseline",
            "attempt": RESULTS / "a_extended_2023_2024_attempts.csv",
            "failure": RESULTS / "a_extended_2023_2024_failure_audit.csv",
        },
        {
            "scenario_tag": "mexico2023_followup_nonID_repeat_upper_bound",
            "cohort_label": "Mexico admission dates in 2023; discharge/source years 2023 and 2024",
            "source_selection": "Subtract exact non-ID row-repeat candidate extras from source years 2023-2024; conservative upper bound",
            "per": RESULTS / "b_extended_2023_repeat_bound_per_unit.csv",
            "pool": RESULTS / "b_extended_2023_repeat_bound_pooled.csv",
            "reference_per": RESULTS / "a_extended_2023_2024_per_unit.csv",
            "reference_label": "2023-2024 extended fit",
            "attempt": RESULTS / "b_extended_2023_repeat_bound_attempts.csv",
            "failure": RESULTS / "b_extended_2023_repeat_bound_failure_audit.csv",
        },
        {
            "scenario_tag": "mexico2023_followup_sourceyear_concordant",
            "cohort_label": "Mexico admission dates in 2023; discharge/source years 2023 and 2024",
            "source_selection": "Retain only source_year=discharge_year; no counts excluded in this 2023 follow-up",
            "per": RESULTS / "b_extended_2023_sourceyear_concordant_per_unit.csv",
            "pool": RESULTS / "b_extended_2023_sourceyear_concordant_pooled.csv",
            "reference_per": RESULTS / "a_extended_2023_2024_per_unit.csv",
            "reference_label": "2023-2024 extended fit",
            "attempt": RESULTS / "b_extended_2023_sourceyear_concordant_attempts.csv",
            "failure": RESULTS / "b_extended_2023_sourceyear_concordant_failure_audit.csv",
        },
    ]

    pooled_summary: list[dict] = []
    attempt_audit: list[dict] = []
    fit_audit_summary: list[dict] = []
    for scenario in scenarios:
        if not scenario["per"].is_file() or not scenario["pool"].is_file():
            raise FileNotFoundError(scenario["scenario_tag"])
        units_for_scenario = {r["unit"] for r in read_csv(scenario["per"])}
        attempt_info = attempts_summary(scenario["attempt"], scenario["per"], scenario["failure"])
        cohort_counts = full_events[scenario_event_key[scenario["scenario_tag"]]]
        totals = {
            "full_calendar_events_K80": cohort_counts["y_K80"],
            "full_calendar_events_K81": cohort_counts["y_K81"],
            "full_calendar_events_urinary": cohort_counts["y_urinary"],
        }
        attempt_audit.append({
            "scenario_tag": scenario["scenario_tag"], "cohort_label": scenario["cohort_label"],
            "source_selection": scenario["source_selection"], **attempt_info, **totals,
            "minimum_complete_dates": 250, "minimum_events_each_endpoint": 40,
        })

        fit_audit_path = RESULTS / (scenario["per"].name.replace("_per_unit.csv", "_fit_audit.csv"))
        if fit_audit_path.is_file():
            fit_rows = read_csv(fit_audit_path)
            fit_audit_summary.append({
                "scenario_tag": scenario["scenario_tag"], "fit_units_audited": len(fit_rows),
                "not_converged": sum(r.get("converged", "").lower() != "true" for r in fit_rows),
                "rank_deficient": sum(r.get("rank", "") != r.get("columns", "") for r in fit_rows),
                "boundary_fits": sum(r.get("boundary", "").lower() == "true" for r in fit_rows),
                "warning_rows": sum(bool(r.get("warnings", "").strip()) for r in fit_rows),
            })

        refs_by_cov = {}
        if scenario["reference_per"]:
            for cov in ("hac14", "hac28"):
                refs_by_cov[cov] = per_unit_rows(scenario["reference_per"], cov)
        unit_counts_by_cov = {}
        for cov in ("hac14", "hac28"):
            sub = per_unit_rows(scenario["per"], cov)
            vals = [float(r["estimate"]) for r in sub]
            unit_counts_by_cov[cov] = {
                "positive": sum(v > 0 for v in vals), "negative": sum(v < 0 for v in vals),
                "zero": sum(v == 0 for v in vals), "n": len(vals),
            }

        for pooled in pooled_rows(scenario["pool"]):
            cov = pooled["covariance"]
            counts = unit_counts_by_cov[cov]
            positive_fraction = counts["positive"] / counts["n"] if counts["n"] else 0.0
            concordance = float(pooled["concordance"])
            if abs(positive_fraction - concordance) > 1e-10:
                raise ValueError(f"unit sign proportion disagrees with pool concordance: {scenario['scenario_tag']} {cov}")
            same_pct, paired_n = ("", "0")
            if scenario["reference_per"]:
                same_pct, paired_n = same_sign_fraction(per_unit_rows(scenario["per"], cov), refs_by_cov[cov])
            pooled_summary.append({
                "scenario_tag": scenario["scenario_tag"], "cohort_label": scenario["cohort_label"],
                "source_selection": scenario["source_selection"], "contrast": pooled["contrast"],
                "covariance": cov.upper(), "pooling_method": pooled["method"],
                "RR": pooled["RR"], "lo": pooled["lo"], "hi": pooled["hi"], "k": pooled["k"],
                "I2_percent": pooled["I2"], "tau2": pooled["tau2"],
                "unit_estimates_positive": counts["positive"], "unit_estimates_negative": counts["negative"],
                "unit_estimates_zero": counts["zero"], "unit_positive_fraction": f"{positive_fraction:.8f}",
                "pool_concordance_positive_fraction": pooled["concordance"],
                "same_direction_vs_reference_percent": same_pct,
                "paired_units_vs_reference": paired_n,
                "direction_reference": scenario["reference_label"],
            })

    write_csv(RESULTS / "mexico_followup_pooled_sensitivity.csv", pooled_summary)
    write_csv(RESULTS / "mexico_followup_count_and_attempt_audit.csv", attempt_audit)
    write_csv(RESULTS / "mexico_followup_fit_audit_summary.csv", fit_audit_summary)

    followup_pool_path = RESULTS / "a_extended_2023_2024_pooled.csv"
    followup_pool = read_csv(followup_pool_path)
    chosen = select_unique(followup_pool, contrast="y_urinary_vs_y_K80", covariance="hac14", method="DL")
    display = [{
        "tag": "mexico2023_followup", "country": "MX", "cohort": "admission dates 2023; discharge/source years 2023 and 2024",
        "contrast": "y_urinary_vs_y_K80", "covariance": "HAC14", "method": "DL",
        "RR": chosen["RR"], "lo": chosen["lo"], "hi": chosen["hi"], "k": chosen["k"],
        "I2": chosen["I2"], "concordance": chosen["concordance"],
        "reference": "per-state 2015-2022 primary E3 basis/P95/shared median fixed",
        "upstream_scenario": chosen["scenario"], "upstream_source_file": followup_pool_path.name,
    }]
    display_path = RESULTS / "mexico2023_followup_display.csv"
    write_csv(display_path, display)

    validation: list[dict] = []
    def add_validation(check: str, passed: bool, evidence: str, value: str = "") -> None:
        validation.append({"check": check, "passed": str(bool(passed)).lower(), "value": value, "evidence": evidence})

    a_source = read_csv(RESULTS / "a_source_2023_calendar_count_check.csv")
    a_per = read_csv(RESULTS / "a_baseline_reproduction_per_unit_check.csv")
    a_pool = read_csv(RESULTS / "a_baseline_reproduction_pooled_check.csv")
    b_primary = read_csv(RESULTS / "b_primary_aggregate_reproduction_check.csv")
    b_extension = read_csv(RESULTS / "b_extension_aggregate_reproduction_check.csv")
    targets = read_csv(RESULTS / "b_registered_primary_integrity_check.csv")
    add_validation("source_2023_rebuild_exact", all(r["exact_reproduction"].lower() == "true" for r in a_source),
                   "a_source_2023_calendar_count_check.csv", str(max(float(r["max_abs_difference"]) for r in a_source)))
    add_validation("old_2023_fit_per_unit_reproduced", max(float(r["max_abs_difference"]) for r in a_per) <= 1e-10,
                   "a_baseline_reproduction_per_unit_check.csv", str(max(float(r["max_abs_difference"]) for r in a_per)))
    add_validation("old_2023_fit_pooled_reproduced", max(float(r["max_abs_difference"]) for r in a_pool) <= 1e-10,
                   "a_baseline_reproduction_pooled_check.csv", str(max(float(r["max_abs_difference"]) for r in a_pool)))
    add_validation("registered_primary_daily_identity_exact", all(r["exact_reproduction"].lower() == "true" for r in b_primary),
                   "b_primary_aggregate_reproduction_check.csv", str(max(float(r["max_abs_difference"]) for r in b_primary)))
    add_validation("extended_2023_all_source_identity_exact", all(r["exact_reproduction"].lower() == "true" for r in b_extension),
                   "b_extension_aggregate_reproduction_check.csv", str(max(float(r["max_abs_difference"]) for r in b_extension)))
    add_validation("registered_primary_fixed_targets_exact", max(float(r["max_abs_difference"]) for r in targets) <= 1e-10,
                   "b_registered_primary_integrity_check.csv", str(max(float(r["max_abs_difference"]) for r in targets)))
    all_fit_audits = [r for r in fit_audit_summary if r["scenario_tag"] != "mexico_primary_registered"]
    add_validation("all_followup_fit_audits_converged_full_rank",
                   all(r["not_converged"] == 0 and r["rank_deficient"] == 0 for r in all_fit_audits),
                   "*_fit_audit.csv", f"fit audit scenarios={len(all_fit_audits)}")
    add_validation("A_manifest_hashes", True, "a_manifest.json", a_manifest["created_utc"])
    add_validation("B_manifest_hashes", True, "b_manifest.json", b_manifest["created_utc"])
    write_csv(RESULTS / "mexico_followup_validation.csv", validation)

    input_paths = [
        RESULTS / "a_manifest.json", RESULTS / "b_manifest.json",
        RESULTS / "a_extended_2023_2024_pooled.csv", RESULTS / "b_sensitivity_deltas_pooled.csv",
        RESULTS / "b_candidate_repeat_subtraction_counts.csv", RESULTS / "b_sourceyear_concordance_event_counts.csv",
        INPUTS / "all_daily.csv.gz", INPUTS / "all_duplicates.csv", INPUTS / "complete.json",
        INPUTS / "count_identity_verified.json", INPUTS / "primary_complete_source_2015_2022_MX_per_unit.csv",
        INPUTS / "primary_complete_source_2015_2022_MX_pooled.csv",
    ]
    summary_manifest = {
        "status": "complete", "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "reporting-only selection from existing per-unit/pooled fits; no model refit",
        "summary_script_sha256": digest(Path(__file__)),
        "input_sha256": {p.name: digest(p) for p in input_paths},
        "output_sha256": {},
        "display_row_source": "a_extended_2023_2024_pooled.csv where contrast=y_urinary_vs_y_K80, covariance=hac14, method=DL",
    }
    for name in ["mexico2023_followup_display.csv", "mexico_followup_pooled_sensitivity.csv",
                 "mexico_followup_count_and_attempt_audit.csv", "mexico_followup_fit_audit_summary.csv",
                 "mexico_followup_validation.csv"]:
        summary_manifest["output_sha256"][name] = digest(RESULTS / name)
    manifest_path = RESULTS / "mexico_followup_summary_manifest.json"
    manifest_path.write_text(json.dumps(summary_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    accepted_handoff = {
        "status": "accepted",
        "csv_sources": [
            {
                "path": str(display_path.resolve()), "role": "followup_results",
                "display_selections": {
                    "mexico2023_followup": {
                        "where": {"tag": "mexico2023_followup"},
                        "fields": {"RR": "RR", "lo": "lo", "hi": "hi"},
                    }
                },
            },
            {"path": str((RESULTS / "mexico_followup_pooled_sensitivity.csv").resolve()), "role": "followup_results"},
            {"path": str((RESULTS / "mexico_followup_count_and_attempt_audit.csv").resolve()), "role": "followup_results"},
            {"path": str((RESULTS / "mexico_followup_validation.csv").resolve()), "role": "followup_results"},
            {"path": str((RESULTS / "b_candidate_repeat_subtraction_counts.csv").resolve()), "role": "followup_results"},
            {"path": str((RESULTS / "b_sourceyear_concordance_event_counts.csv").resolve()), "role": "followup_results"},
        ],
    }
    (HERE / "accepted_handoff.json").write_text(json.dumps(accepted_handoff, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

