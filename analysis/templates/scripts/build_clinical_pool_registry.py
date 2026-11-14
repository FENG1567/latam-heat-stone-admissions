#!/usr/bin/env python3
"""Build compact selector and coverage registries from accepted clinical pools."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]
BR = BASE / "results" / "clinical_full_BR"
MX_BASELINE = BASE / "results" / "clinical_full_MX_baseline"
MX_COMPLETE = BASE / "results" / "clinical_full_MX_complete"


SUMMARY_FIELDS = [
    "country", "result_type", "selection_role", "preferred_selector", "scenario",
    "source_selector", "source_years", "admission_years", "nuisance", "calendar_mode",
    "endpoint", "endpoint_label", "mode", "group", "contrast", "weight_set", "covariance", "pooling",
    "units_total", "units_passing_event_count", "units_below_event_count",
    "units_event_count_unclassified", "units_fit", "units_fit_failed",
    "units_count_eligible_but_fit_failed", "units_joint_fit_successful", "units_joint_fit_failed",
    "k", "logRRR", "se", "RRR", "RR_lo", "RR_hi",
    "I2", "tau2", "concordance", "median_joint_wald_chisq3", "joint_wald_p_lt_0_05_units",
    "status", "failure",
]

COVERAGE_FIELDS = [
    "country", "scenario", "cohort_role", "source_selector", "source_years", "admission_years",
    "units_in_calendar_window", "units_threshold_ineligible", "units_eligible",
    "units_fit_successfully", "units_nonthreshold_failed", "units_numeric_failure",
    "units_other_failure", "pooled_result_available", "primary_weight_set", "primary_pooling",
    "primary_k", "primary_RRR", "primary_RR_lo", "primary_RR_hi", "availability_source",
    "failure_source", "estimate_source",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_selector(scenario: str, label: str) -> tuple[str, str, str]:
    if scenario == "latebilling_reconstruction_2014_2021":
        return label, "2014-2021 (source-year 2021 adds late 2020 billing)", "2014-2020"
    if scenario == "matched_source_2014_2020":
        return label, "2014-2020", "2014-2020"
    if scenario == "mx_primary_panel_year_exact_baseline":
        return label, "2015-2018 and 2020-2023; primary_panel_year only", "2015-2018 and 2020-2022"
    if scenario == "mx_complete_source_2015_2022":
        return label, "2015-2022; all source roles", "2015-2022"
    return label, "", ""


def row_base(result_type: str, row: dict[str, str], country: str | None = None) -> dict[str, object]:
    scenario = row.get("scenario", "")
    label = row.get("scenario_label", "")
    selector, source_years, admission_years = source_selector(scenario, label)
    return {
        "country": country or row.get("country", ""),
        "result_type": result_type,
        "scenario": scenario,
        "source_selector": selector,
        "source_years": source_years,
        "admission_years": admission_years,
        "nuisance": row.get("nuisance", ""),
        "calendar_mode": {"separate": "outcome_specific", "shared": "shared"}.get(row.get("nuisance", ""), ""),
        "endpoint": row.get("endpoint", ""),
        "endpoint_label": row.get("endpoint_label", ""),
        "mode": row.get("mode", ""),
        "group": row.get("group", ""),
        "weight_set": row.get("weight_set", ""),
        "covariance": row.get("covariance", ""),
        "pooling": row.get("pooling", ""),
        "units_total": row.get("units_total", ""),
        "units_passing_event_count": row.get("units_passing_event_count", ""),
        "units_below_event_count": row.get("units_below_event_count", ""),
        "units_event_count_unclassified": row.get("units_event_count_unclassified", ""),
        "units_fit": row.get("units_fit", ""),
        "units_fit_failed": row.get("units_fit_failed", row.get("units_joint_fit_failed", "")),
        "units_count_eligible_but_fit_failed": row.get("units_count_eligible_but_fit_failed", ""),
        "units_joint_fit_successful": row.get("units_joint_fit_successful", ""),
        "units_joint_fit_failed": row.get("units_joint_fit_failed", ""),
        "k": row.get("k", ""),
        "logRRR": row.get("logRRR", row.get("logRRR_difference", "")),
        "se": row.get("se", ""),
        "RRR": row.get("RRR", row.get("RR_ratio", "")),
        "RR_lo": row.get("RR_lo", ""),
        "RR_hi": row.get("RR_hi", ""),
        "I2": row.get("I2", ""),
        "tau2": row.get("tau2", ""),
        "concordance": row.get("concordance", ""),
        "median_joint_wald_chisq3": row.get("median_joint_wald_chisq3", ""),
        "joint_wald_p_lt_0_05_units": row.get("joint_wald_p_lt_0_05_units", ""),
        "status": row.get("status", ""),
        "failure": row.get("failure", ""),
    }


def role_for(result_type: str, row: dict[str, object]) -> tuple[str, bool]:
    scenario = str(row.get("scenario", ""))
    nuisance = str(row.get("nuisance", ""))
    covariance = str(row.get("covariance", ""))
    pooling = str(row.get("pooling", ""))
    weight_set = str(row.get("weight_set", ""))
    primary_scenario = scenario == "latebilling_reconstruction_2014_2021"
    matched = scenario == "matched_source_2014_2020"
    if result_type == "endpoint_min30_events":
        return "min30_event_count_sensitivity", False
    if result_type in {"endpoint", "interaction", "admission_group"}:
        if primary_scenario and nuisance == "separate" and covariance == "hac14" and pooling == "DL":
            if result_type != "admission_group" or (row.get("mode") == "admission_type" and row.get("group") in {"elective", "urgent"}):
                return "preferred_primary_selector", True
        if matched and nuisance == "separate" and covariance == "hac14" and pooling == "DL":
            return "matched_source_sensitivity", False
        if primary_scenario and nuisance == "shared" and covariance == "hac14" and pooling == "DL":
            return "shared_calendar_sensitivity", False
        return "covariance_pooling_or_calendar_sensitivity", False
    if result_type in {"case_mix", "case_mix_ageband"}:
        if primary_scenario and weight_set == "pooled_observed_case_mix" and pooling == "DL":
            return "preferred_primary_selector", True
        if weight_set == "uniform_1_over_6":
            return "uniform_weight_sensitivity", False
        return "cohort_or_ageband_sensitivity", False
    return "other", False


def add_pooled_rows(summary: list[dict[str, object]]) -> None:
    endpoint_file_map = [
        ("endpoint", BR / "clinical_brazil_endpoint_pooled.csv"),
        ("endpoint_min30_events", BR / "clinical_brazil_endpoint_pooled_min30_events.csv"),
        ("interaction", BR / "clinical_brazil_interaction_coverage.csv"),
        ("admission_group", BR / "clinical_brazil_interaction_group_pooled.csv"),
        ("case_mix", BR / "clinical_brazil_case_mix_pooled.csv"),
        ("case_mix_ageband", BR / "clinical_brazil_case_mix_ageband_pooled.csv"),
    ]
    for result_type, path in endpoint_file_map:
        for raw in read_csv(path):
            row = row_base(result_type, raw)
            if result_type == "endpoint_min30_events":
                row["status"] = raw.get("status", "")
            if result_type == "interaction":
                row["contrast"] = {
                    "admission_type": "urgent minus elective log-RRR; joint three-coefficient Wald test",
                    "sex": "male minus female log-RRR; joint three-coefficient Wald test",
                }.get(raw.get("mode", ""), "")
                row["covariance"] = "hac14"
            elif result_type == "admission_group":
                row["contrast"] = raw.get("contrast", "")
                row["covariance"] = "hac14"
            elif result_type == "case_mix":
                row["mode"] = "adult_common_case_mix_standardization"
                row["nuisance"] = "not_applicable"
                row["calendar_mode"] = "stratum_specific_date_aligned"
                row["covariance"] = "hac14"
            elif result_type == "case_mix_ageband":
                row["mode"] = "age_band_standardization"
                row["nuisance"] = "not_applicable"
                row["calendar_mode"] = "stratum_specific_date_aligned"
                row["covariance"] = "hac14"
                row["weight_set"] = raw.get("age_group", "")
            role, preferred = role_for(result_type, row)
            row["selection_role"] = role
            row["preferred_selector"] = str(preferred).lower()
            summary.append(row)


def case_mix_coverage_row(country: str, scenario: str, label: str, total: int, success: int,
                          failures: list[dict[str, str]], availability_source: Path,
                          failure_source: Path, pooled_rows: list[dict[str, str]], primary_weight_set: str) -> dict[str, object]:
    threshold_failures = [
        row for row in failures
        if row.get("failure", "").lower().startswith("fewer than 30 cases")
        or "fewer than 400 eligible days" in row.get("failure", "").lower()
    ]
    numeric_failures = [
        row for row in failures
        if any(token in row.get("failure", "").lower() for token in ("quasi-poisson", "hessian", "singular"))
    ]
    other_failures = [row for row in failures if row not in threshold_failures and row not in numeric_failures]
    selector, source_years, admission_years = source_selector(scenario, label)
    preferred_pool = next((r for r in pooled_rows
                           if r.get("scenario") == scenario
                           and r.get("weight_set") == primary_weight_set
                           and r.get("pooling") == "DL"), None)
    if country == "brazil":
        role = "preferred_primary_cohort" if scenario == "latebilling_reconstruction_2014_2021" else "matched_source_sensitivity"
    else:
        role = "preferred_complete_source_cohort" if scenario == "mx_complete_source_2015_2022" else "original_gap_comparator"
    return {
        "country": country,
        "scenario": scenario,
        "cohort_role": role,
        "source_selector": selector,
        "source_years": source_years,
        "admission_years": admission_years,
        "units_in_calendar_window": total,
        "units_threshold_ineligible": len(threshold_failures),
        "units_eligible": total - len(threshold_failures),
        "units_fit_successfully": success,
        "units_nonthreshold_failed": len(failures) - len(threshold_failures),
        "units_numeric_failure": len(numeric_failures),
        "units_other_failure": len(other_failures),
        "pooled_result_available": str(preferred_pool is not None).lower(),
        "primary_weight_set": primary_weight_set,
        "primary_pooling": "DL",
        "primary_k": preferred_pool.get("k", "") if preferred_pool else "",
        "primary_RRR": preferred_pool.get("RRR", "") if preferred_pool else "",
        "primary_RR_lo": preferred_pool.get("RR_lo", "") if preferred_pool else "",
        "primary_RR_hi": preferred_pool.get("RR_hi", "") if preferred_pool else "",
        "availability_source": str(availability_source.relative_to(BASE)),
        "failure_source": str(failure_source.relative_to(BASE)),
        "estimate_source": str((availability_source.parent / "clinical_mexico_case_mix_pooled.csv" if country == "mexico"
                                 else BR / "clinical_brazil_case_mix_pooled.csv").relative_to(BASE)),
    }


def build_coverage_registry() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for directory, country in ((BR, "brazil"), (MX_BASELINE, "mexico"), (MX_COMPLETE, "mexico")):
        availability = read_csv(directory / f"clinical_{country}_case_mix_availability.csv")
        failures = read_csv(directory / f"clinical_{country}_case_mix_failures.csv")
        pooled_path = directory / f"clinical_{country}_case_mix_pooled.csv"
        pooled = read_csv(pooled_path)
        elig_rows = read_csv(directory / f"clinical_{country}_case_mix_eligibility.csv")
        total_by_scenario = {r["scenario"]: int(r["units_in_calendar_window"]) for r in availability}
        success_by_scenario = {r["scenario"]: int(r["units_fit_successfully"]) for r in availability}
        labels = {r["scenario"]: r["scenario_label"] for r in availability}
        if country == "brazil":
            primary_weight_set = "pooled_observed_case_mix"
        else:
            primary_weight_set = "pooled_observed_case_mix"
        for scenario in sorted(total_by_scenario):
            scenario_failures = [r for r in failures if r.get("scenario") == scenario]
            scenario_pooled = [r for r in pooled if r.get("scenario") == scenario]
            row = case_mix_coverage_row(
                country, scenario, labels[scenario], total_by_scenario[scenario],
                success_by_scenario[scenario], scenario_failures, directory / f"clinical_{country}_case_mix_availability.csv",
                directory / f"clinical_{country}_case_mix_failures.csv", scenario_pooled, primary_weight_set,
            )
            row["adult_known_sex_records_used_for_weights"] = elig_rows[0].get("adult_known_sex_records_used_for_weights", "") if elig_rows else ""
            rows.append(row)
    return rows


def build_case_mix_directional_summary() -> list[dict[str, object]]:
    sources = [
        (BR, "brazil"),
        (MX_COMPLETE, "mexico"),
    ]
    output: list[dict[str, object]] = []
    for directory, country in sources:
        per_unit = read_csv(directory / f"clinical_{country}_case_mix_standardized_per_unit.csv")
        pooled = read_csv(directory / f"clinical_{country}_case_mix_pooled.csv")
        availability = read_csv(directory / f"clinical_{country}_case_mix_availability.csv")
        availability_by_scenario = {r["scenario"]: r for r in availability}
        pooled_dl = {(r["scenario"], r["weight_set"]): r for r in pooled if r.get("pooling") == "DL"}
        groups: dict[tuple[str, str], list[dict[str, str]]] = {}
        for row in per_unit:
            groups.setdefault((row["scenario"], row["weight_set"]), []).append(row)
        for (scenario, weight_set), rows in sorted(groups.items()):
            valid = []
            for row in rows:
                try:
                    estimate = float(row.get("logRRR", ""))
                    variance = float(row.get("variance", ""))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(estimate) and math.isfinite(variance) and variance > 0:
                    valid.append((row, estimate))
            pool = pooled_dl.get((scenario, weight_set))
            if pool is None:
                raise ValueError(f"Missing DL case-mix pool for {(country, scenario, weight_set)}")
            k = len(valid)
            n_positive = sum(estimate > 0 for _, estimate in valid)
            if int(pool["k"]) != k:
                raise ValueError(f"Case-mix directional k mismatch for {(country, scenario, weight_set)}: {k} != {pool['k']}")
            avail = availability_by_scenario[scenario]
            label = rows[0].get("scenario_label", "")
            selector, source_years, admission_years = source_selector(scenario, label)
            output.append({
                "country": country, "scenario": scenario, "scenario_label": label,
                "source_selector": selector, "source_years": source_years, "admission_years": admission_years,
                "weight_set": weight_set, "covariance": "HAC14", "pooling": "DL",
                "units_total": avail.get("units_in_calendar_window", ""),
                "units_fit_successfully": avail.get("units_fit_successfully", ""),
                "units_fit_failed": avail.get("units_failed_or_unavailable", ""),
                "k": k, "n_positive": n_positive, "n_nonpositive": k - n_positive,
                "concordance": n_positive / k if k else "",
                "RRR": pool.get("RRR", ""), "RR_lo": pool.get("RR_lo", ""), "RR_hi": pool.get("RR_hi", ""),
            })
    return output


def main() -> None:
    summary: list[dict[str, object]] = []
    add_pooled_rows(summary)
    summary.sort(key=lambda r: (
        str(r.get("result_type", "")), str(r.get("scenario", "")), str(r.get("nuisance", "")),
        str(r.get("endpoint", "")), str(r.get("mode", "")), str(r.get("weight_set", "")),
        str(r.get("covariance", "")), str(r.get("pooling", "")),
    ))
    summary_path = BR / "clinical_brazil_national_pooled_summary.csv"
    write_csv(summary_path, summary, SUMMARY_FIELDS)

    coverage = build_coverage_registry()
    coverage_path = BASE / "results" / "clinical_case_mix_coverage_summary.csv"
    write_csv(coverage_path, coverage, COVERAGE_FIELDS)
    directional = build_case_mix_directional_summary()
    directional_path = BASE / "results" / "clinical_case_mix_directional_summary.csv"
    directional_fields = [
        "country", "scenario", "scenario_label", "source_selector", "source_years", "admission_years",
        "weight_set", "covariance", "pooling", "units_total", "units_fit_successfully", "units_fit_failed",
        "k", "n_positive", "n_nonpositive", "concordance", "RRR", "RR_lo", "RR_hi",
    ]
    write_csv(directional_path, directional, directional_fields)

    baseline_pool = read_csv(MX_BASELINE / "clinical_mexico_case_mix_pooled.csv")
    complete_pool = read_csv(MX_COMPLETE / "clinical_mexico_case_mix_pooled.csv")
    if len(summary) != 308:
        raise SystemExit(f"Expected 308 merged Brazil pool rows, observed {len(summary)}")
    mx = {(r["country"], r["scenario"]): r for r in coverage}
    expected = {
        ("mexico", "mx_primary_panel_year_exact_baseline"): (32, 0, 13, 19),
        ("mexico", "mx_complete_source_2015_2022"): (32, 16, 10, 6),
    }
    for key, counts in expected.items():
        row = mx[key]
        observed = (
            int(row["units_in_calendar_window"]), int(row["units_fit_successfully"]),
            int(row["units_threshold_ineligible"]), int(row["units_numeric_failure"]),
        )
        if observed != counts:
            raise SystemExit(f"Mexico case-mix coverage mismatch for {key}: {observed} != {counts}")
    if baseline_pool:
        raise SystemExit("Old-gap Mexico case-mix pool is expected to be empty and must not be registered as an estimate")
    if not any(r.get("weight_set") == "pooled_observed_case_mix" and r.get("pooling") == "DL"
               and r.get("scenario") == "mx_complete_source_2015_2022" for r in complete_pool):
        raise SystemExit("Preferred complete-source Mexico observed-case-mix DL pool is missing")

    manifest = {
        "summary_rows": len(summary),
        "summary_path": str(summary_path.relative_to(BASE)),
        "summary_sha256": sha256(summary_path),
        "coverage_rows": len(coverage),
        "coverage_path": str(coverage_path.relative_to(BASE)),
        "coverage_sha256": sha256(coverage_path),
        "directional_rows": len(directional),
        "directional_path": str(directional_path.relative_to(BASE)),
        "directional_sha256": sha256(directional_path),
        "source_files": {
            str(path.relative_to(BASE)): sha256(path)
            for path in [
                BR / "clinical_brazil_endpoint_pooled.csv",
                BR / "clinical_brazil_endpoint_pooled_min30_events.csv",
                BR / "clinical_brazil_interaction_coverage.csv",
                BR / "clinical_brazil_interaction_group_pooled.csv",
                BR / "clinical_brazil_case_mix_pooled.csv",
                BR / "clinical_brazil_case_mix_ageband_pooled.csv",
                MX_BASELINE / "clinical_mexico_case_mix_availability.csv",
                MX_BASELINE / "clinical_mexico_case_mix_failures.csv",
                MX_BASELINE / "clinical_mexico_case_mix_pooled.csv",
                MX_COMPLETE / "clinical_mexico_case_mix_availability.csv",
                MX_COMPLETE / "clinical_mexico_case_mix_failures.csv",
                MX_COMPLETE / "clinical_mexico_case_mix_pooled.csv",
            ]
        },
        "mexico_case_mix_expected_coverage": {
            "mx_primary_panel_year_exact_baseline": {"units_total": 32, "fit_success": 0, "threshold_ineligible": 13, "numeric_failure": 19, "estimate_registered": False},
            "mx_complete_source_2015_2022": {"units_total": 32, "fit_success": 16, "threshold_ineligible": 10, "numeric_failure": 6, "observed_case_mix_DL_k": 16},
        },
    }
    manifest_path = BASE / "admin" / "clinical_pool_registry_manifest_20261113.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary_rows": len(summary), "summary_sha256": manifest["summary_sha256"],
                      "coverage_rows": len(coverage), "coverage_sha256": manifest["coverage_sha256"],
                      "directional_rows": len(directional), "directional_sha256": manifest["directional_sha256"],
                      "mexico_baseline_pool_rows": len(baseline_pool),
                      "mexico_complete_pool_rows": len(complete_pool)}, indent=2))


if __name__ == "__main__":
    main()
