#!/usr/bin/env python3
"""Render four main tables and accepted-result supplements from a verified freeze.

This script formats frozen observations only. It does not estimate effects,
re-pool geographic units, compute confidence limits, or read source CSVs. Run it
only after the root declares the final checksum-bound freeze ready.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from publication_common import FrozenPublication, count, md_table, num, percent, source_csv


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLE4_REQUIRED_TAGS = {
    "br_strict", "br_urgent", "br_elective", "br_admission_interaction",
    "br_later", "mx_later",
    "br_moisture_dew", "br_moisture_rh", "mx_moisture_dew",
    "mx_moisture_rh", "mx_moisture_rh_hk",
}
TABLE4_ADMISSION_TAGS = {"br_urgent", "br_elective", "br_admission_interaction"}
CURVE_SCOPE_BY_KEY = {
    "curves_supported": "support_restricted",
    "curves_all_units": "all_successful",
}
CURVE_METHODS = {"DL", "REML_profile_HK"}
CSV_ONLY_SUPPLEMENT_KEYS = {"main_mrt_selection_unit", "temporal_support_units"}
SUPPLEMENT_KEY_ROW_LIMIT = 500
SUPPLEMENT_SECTION_ROW_LIMIT = 900
SUPPLEMENT_COLUMN_LIMIT = 14

SUPPLEMENT_DIMENSION_FIELDS = (
    "country", "dataset", "cohort", "scenario", "scenario_label", "period",
    "outcome", "endpoint", "contrast", "estimand", "method", "pooling",
    "covariance", "scope", "curve_scope", "temperature", "mode", "group", "stratum",
    "category", "subgroup", "age_group", "age_band", "ageband", "weight_set",
    "admission_type", "setting", "origin", "sample_scope", "exposure", "calendar",
    "metric", "measure",
    "status", "omitted_unit", "year", "source_year", "nuisance", "adjustment",
    "weather_variable", "moisture_metric", "diagnosis", "admission_origin",
    "origin_class", "origin_group", "classification", "reference_period",
)
SUPPLEMENT_ESSENTIAL_FIELDS = (
    "country", "dataset", "cohort", "scenario", "scenario_label", "period",
    "outcome", "endpoint", "contrast", "method", "pooling", "covariance",
    "mode", "group", "stratum", "weight_set", "age_group", "age_band", "ageband",
    "scope", "curve_scope", "sample_scope", "temperature", "admission_type",
    "origin", "nuisance",
)
SUPPLEMENT_METRIC_FIELDS = (
    "RRR", "RR", "RR_ratio", "estimate", "point", "logRRR", "logRRR_difference",
    "lo", "hi", "RR_lo", "RR_hi", "RR_ratio_lo", "RR_ratio_hi", "se", "I2",
    "tau2", "k", "units_fit", "units_joint_fit_successful", "units_total",
    "units_failed", "units_fit_failed", "units_joint_fit_failed", "concordance",
    "p", "p_value", "joint_wald_p_lt_0_05_units", "n", "n_days", "eligible_rows",
    "row_count", "distinct_date_count", "total_events", "events", "event_count",
    "mean_C", "median_C", "p05_C", "p95_C", "main_k", "main_country_k_min",
    "main_country_k_max", "successful_regions_main_unit", "failed_regions_main_effect",
    "attempted_regions_success_plus_failure", "observed_mrt_boundary_P10_regions",
    "observed_mrt_boundary_P60_regions", "observed_mrt_boundary_interior_regions",
    "observed_mrt_boundary_missing_or_other_regions", "mean_simulated_boundary_frequency",
    "median_simulated_boundary_frequency", "mean_simulated_P10_boundary_frequency",
    "median_simulated_P10_boundary_frequency", "pooled_simulated_P10_boundary_frequency",
    "mean_simulated_P60_boundary_frequency", "median_simulated_P60_boundary_frequency",
    "pooled_simulated_P60_boundary_frequency", "simulation_draws_per_region",
)
SUPPLEMENT_METRIC_TOKENS = (
    "count", "events", "units", "regions", "days", "rows", "mean", "median",
    "frequency", "boundary", "eligible", "failed", "successful", "support", "p_value",
)
SUPPLEMENT_BLOCKED_FIELD_TOKENS = (
    "source_file", "source_row", "_json", "json", "checksum", "sha256", "path",
    "long_unit", "unit_ids", "case_id", "raw_code", "artifact_id",
)


def _record_claims(pub: FrozenPublication, key: str, rows: list[dict], where: dict,
                   artifact: str) -> list[dict]:
    if not rows:
        raise ValueError(f"Frozen selection returned no rows: {key} {where}")
    for row in rows:
        pub.claims.append({
            "artifact": artifact,
            "table": key,
            "selection": where,
            "source_file": row["source_file"],
            "source_row": row["source_row"],
        })
    return rows


def _rows(pub: FrozenPublication, key: str, artifact: str, **where) -> list[dict]:
    if key not in pub.tables:
        raise KeyError(f"Required frozen table key is missing: {key}")
    return _record_claims(pub, key, pub.rows(key, **where), where, artifact)


def _one(pub: FrozenPublication, key: str, artifact: str, **where) -> dict:
    rows = _rows(pub, key, artifact, **where)
    if len(rows) != 1:
        raise ValueError(f"Expected one frozen row: {key} {where}; found {len(rows)}")
    return rows[0]


def _select_from_publication_metadata(pub: FrozenPublication, tag: str, artifact: str) -> dict:
    selections = pub.metadata.get("publication_selections", {})
    if tag not in selections:
        raise KeyError(f"Final metadata is missing publication selector {tag!r}")
    spec = selections[tag]
    suffix = spec["source_suffix"]
    where = spec["where"]
    rows = pub.source_rows(suffix, **where)
    if len(rows) != 1:
        raise ValueError(f"Expected one selected source for {tag}: found {len(rows)}")
    row = rows[0]
    pub.claims.append({"artifact": artifact, "selector_tag": tag,
                       "selection": where, "source_file": row["source_file"],
                       "source_row": row["source_row"]})
    return row


def _source_ref(row: dict) -> str:
    return "{}:{}".format(row["source_file"], row["source_row"])


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    source_csv(rows, path)


def _write_table(path: Path, number: int, caption: str, headers: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"## Table {number}\n\n" + caption + "\n\n" + md_table(headers, rows) + "\n",
        encoding="utf-8",
    )


def _supplement_rows(key: str, records: list[dict]) -> list[dict]:
    if key not in CURVE_SCOPE_BY_KEY:
        return records
    scope = CURVE_SCOPE_BY_KEY[key]
    methods = {str(row.get("method", "")) for row in records}
    unexpected_methods = methods - CURVE_METHODS
    if unexpected_methods:
        raise ValueError(
            f"Verified curve table {key} has unaccepted methods: {sorted(unexpected_methods)}"
        )
    selected = [row for row in records if str(row.get("curve_scope", "")) == scope]
    if not selected:
        raise ValueError(f"Verified curve table {key} has no {scope!r} rows")
    unique = {}
    for row in selected:
        identity = (row.get("source_file"), row.get("source_row"))
        if identity in unique and unique[identity] != row:
            raise ValueError(f"Conflicting verified curve row identity: {identity}")
        unique[identity] = row
    return list(unique.values())


def _supplement_field_blocked(field: str, rows: list[dict]) -> bool:
    normalized = field.casefold()
    if any(token in normalized for token in SUPPLEMENT_BLOCKED_FIELD_TOKENS):
        return True
    return any(len(str(row.get(field, ""))) > 120 for row in rows if row.get(field) is not None)


def _supplement_columns(key: str, rows: list[dict]) -> list[str]:
    fields = list(dict.fromkeys(field for row in rows for field in row))
    usable = [field for field in fields if not _supplement_field_blocked(field, rows)]
    if key == "pub_main_models":
        preferred = (
            "country", "endpoint", "main_k", "failed_regions_main_effect",
            "observed_mrt_boundary_P10_regions", "observed_mrt_boundary_P60_regions",
            "observed_mrt_boundary_interior_regions",
            "observed_mrt_boundary_missing_or_other_regions",
            "mean_simulated_P10_boundary_frequency", "median_simulated_P10_boundary_frequency",
            "pooled_simulated_P10_boundary_frequency",
            "mean_simulated_P60_boundary_frequency", "median_simulated_P60_boundary_frequency",
            "pooled_simulated_P60_boundary_frequency",
        )
        return [field for field in preferred if field in usable]
    dimensions = [field for field in SUPPLEMENT_ESSENTIAL_FIELDS if field in usable]
    metric_candidates = [field for field in SUPPLEMENT_METRIC_FIELDS if field in usable]
    metric_tokens = tuple(token.casefold() for token in SUPPLEMENT_METRIC_TOKENS)
    additional_metrics = []
    for field in usable:
        normalized = field.casefold()
        if field in dimensions or field in metric_candidates:
            continue
        if (any(token in normalized for token in metric_tokens)
                and any(isinstance(row.get(field), (int, float))
                        or (isinstance(row.get(field), str) and _is_number(row[field]))
                        for row in rows if row.get(field) not in (None, ""))):
            additional_metrics.append(field)
    metrics_available = metric_candidates + additional_metrics
    ordered_metrics = []
    point_field = next(
        (field for field in ("RRR", "RR", "RR_ratio", "estimate", "point", "logRRR")
         if field in metrics_available),
        None,
    )
    if point_field:
        ordered_metrics.append(point_field)
    low_field = next(
        (field for field in ("lo", "RR_lo", "RR_ratio_lo", "lower") if field in metrics_available),
        None,
    )
    high_field = next(
        (field for field in ("hi", "RR_hi", "RR_ratio_hi", "upper") if field in metrics_available),
        None,
    )
    for field in (low_field, high_field, "k", "units_fit", "units_joint_fit_successful",
                  "I2", "concordance", "p", "p_value", "n", "row_count",
                  "eligible_rows", "total_events", "events"):
        if field and field in metrics_available and field not in ordered_metrics:
            ordered_metrics.append(field)
    ordered_metrics.extend(field for field in metrics_available if field not in ordered_metrics)
    metric_budget = max(1, SUPPLEMENT_COLUMN_LIMIT - len(dimensions))
    metrics = ordered_metrics[:metric_budget]
    optional_dimensions = [
        field for field in SUPPLEMENT_DIMENSION_FIELDS
        if field in usable and field not in dimensions
    ]
    optional_budget = max(0, SUPPLEMENT_COLUMN_LIMIT - len(dimensions) - len(metrics))
    return dimensions + optional_dimensions[:optional_budget] + metrics


def _is_number(value) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def _supplement_cell(field: str, value) -> str:
    if value is None or str(value).strip() == "":
        return "NA"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if _is_number(value):
        numeric = float(value)
        normalized = field.casefold()
        if "year" in normalized and numeric.is_integer():
            return str(int(numeric))
        if "temperature" in normalized:
            return num(numeric, 1)
        if "concordance" in normalized or "frequency" in normalized:
            return percent(numeric * 100 if 0 <= numeric <= 1 else numeric, 1)
        count_tokens = ("count", "events", "units", "regions", "days", "rows", "draws")
        if ((any(token in normalized for token in count_tokens)
                or normalized in {"k", "n", "main_k", "eligible_rows", "row_count"})
                and numeric.is_integer()):
            return f"{int(numeric):,}"
        if normalized in {"i2"}:
            return num(numeric, 1)
        return num(numeric, 3)
    return str(value)


def _supplement_header(field: str) -> str:
    labels = {
        "I2": "I² (%)",
        "RRR": "RRR",
        "RR": "RR",
        "RR_ratio": "RR ratio",
        "RR_lo": "Lower 95% CI",
        "RR_hi": "Upper 95% CI",
        "concordance": "Direction >0 (%)",
    }
    if field in labels:
        return labels[field]
    title = field.replace("_", " ").replace("I2", "I²").title()
    return title + " (%)" if "frequency" in field.casefold() else title


def _supplement_key_markdown(key: str, rows: list[dict]) -> tuple[list[str], int]:
    csv_link = f"[Registered source rows (CSV)](source_data/{key}.csv)"
    if key in CSV_ONLY_SUPPLEMENT_KEYS:
        return [
            f"#### {key}", "", f"{csv_link}. Unit-level detail is provided in the CSV.", "",
        ], 0
    if len(rows) > SUPPLEMENT_KEY_ROW_LIMIT:
        return [
            f"#### {key}", "",
            f"{len(rows):,} registered rows are retained in the CSV because this table exceeds the compact Markdown row limit. {csv_link}.",
            "",
        ], 0
    columns = _supplement_columns(key, rows)
    if not columns:
        return [
            f"#### {key}", "", f"{csv_link}. No compact numeric fields were available for Markdown display.", "",
        ], 0
    display_rows = [
        [_supplement_cell(field, row.get(field)) for field in columns]
        for row in rows
    ]
    lines = [
        f"#### {key}", "",
        f"{csv_link}; {len(rows):,} registered rows. Exact source filenames and row numbers remain in the CSV.",
        "",
        md_table([_supplement_header(field) for field in columns], display_rows),
        "",
    ]
    return lines, len(rows)


def _table1(pub: FrozenPublication, tables_dir: Path) -> tuple[list[dict], str]:
    output = []
    display_rows = []
    for dataset, country, label in [
        ("BR_calendar", "BR", "Brazil, primary 2006–2020"),
        ("MX_2015_2022", "MX", "Mexico, preferred complete-source 2015–2022"),
    ]:
        calendar = _one(pub, "pub_calendar", "table1", dataset=dataset, country=country)
        exposure = _one(pub, "pub_exposure", "table1", dataset=dataset, country=country)
        outcome_pool = pub.rows("pub_outcomes", dataset=dataset, country=country,
                                row_scope="all_health")
        outcome_rows = [r for r in outcome_pool if r["outcome"] in {"K80", "K81", "urinary"}]
        outcome_rows = _record_claims(
            pub, "pub_outcomes", outcome_rows,
            {"dataset": dataset, "country": country, "row_scope": "all_health",
             "outcome_in": ["K80", "K81", "urinary"]}, "table1")
        by_outcome = {r["outcome"]: r for r in outcome_rows}
        if set(by_outcome) != {"K80", "K81", "urinary"}:
            raise ValueError(f"Table 1 outcome totals incomplete for {dataset}")
        source_rows = [calendar, exposure, *outcome_rows]
        refs = "; ".join(_source_ref(r) for r in source_rows)
        row = {
            "country_period": label,
            "geographic_units": calendar["unit_count"],
            "calendar_unit_days": calendar["row_count"],
            "calendar_dates": calendar["distinct_date_count"],
            "eligible_unit_days": calendar["eligible_rows"],
            "zero_admission_eligible_days": calendar["all3_zero_admission_days_among_eligible"],
            "mean_E3_C": exposure["mean_C"],
            "median_E3_C": exposure["median_C"],
            "P05_E3_C": exposure["p05_C"],
            "P95_E3_C": exposure["p95_C"],
            "K80_admission_records": by_outcome["K80"]["total_events"],
            "K81_admission_records": by_outcome["K81"]["total_events"],
            "urinary_N20_N23_admission_records": by_outcome["urinary"]["total_events"],
            "source_rows": refs,
        }
        output.append(row)
        display_rows.append([
            label, count(calendar["unit_count"]), count(calendar["row_count"]),
            count(calendar["eligible_rows"]),
            count(calendar["all3_zero_admission_days_among_eligible"]),
            num(exposure["mean_C"], 2), num(exposure["median_C"], 2),
            num(exposure["p05_C"], 2), num(exposure["p95_C"], 2),
            count(by_outcome["K80"]["total_events"]),
            count(by_outcome["K81"]["total_events"]),
            count(by_outcome["urinary"]["total_events"]),
            refs,
        ])
    _write_csv(tables_dir / "table1.csv", output)
    caption = (
        "Table 1 | Coverage, admission totals and three-day temperature distribution. "
        "Calendar unit-days and exposure-eligible unit-days are separate. Admission totals "
        "cover the full primary calendar period and are records, not unique people; zero-admission "
        "days are eligible calendar unit-days. "
        "Exposure percentiles summarize the pooled eligible unit-day distribution; they are not "
        "means of unit-specific percentiles and are distinct from the unit-specific P95 targets "
        "used in association models."
    )
    _write_table(tables_dir / "table1.md", 1, caption,
                 ["Country and period", "Units", "Calendar unit-days", "Eligible unit-days",
                  "Zero-admission eligible days", "Mean E3 °C", "Median E3 °C", "P05 E3 °C",
                  "P95 E3 °C", "K80 records", "K81 records", "N20–N23 records", "Source rows"],
                 display_rows)
    return output, caption


def _table2(pub: FrozenPublication, tables_dir: Path) -> tuple[list[dict], str]:
    primary_covariance = pub.metadata.get("definitions", {}).get("primary_covariance", "")
    covariance_label = "Primary calendar-aligned Bartlett HAC14"
    output = []
    display_rows = []
    main_rows = pub.rows("main_effects", country="BR", contrast="p95_vs_mrt")
    main_rows = [r for r in main_rows if r["method"] in {"DL", "REML_HK"}]
    expected = {"y_K80", "y_K81", "y_urinary", "y_sum3"}
    if {r["endpoint"] for r in main_rows} != expected:
        raise ValueError("Table 2 must include all four Brazilian main-effect endpoints")
    main_rows = _record_claims(pub, "main_effects", main_rows,
                               {"country": "BR", "contrast": "p95_vs_mrt",
                                "method_in": ["DL", "REML_HK"]}, "table2")
    for row in main_rows:
        item = {
            "country": "BR", "cohort": "Brazil primary 2006–2020",
            "endpoint": row["endpoint"], "estimand": "RR, P95 versus constrained MRT",
            "point": row["RR"], "lo": row["lo"], "hi": row["hi"],
            "method": row["method"], "covariance": covariance_label,
            "k": row["k"], "I2": row["I2"], "concordance": row["concordance"],
            "source_file": row["source_file"], "source_row": row["source_row"],
        }
        output.append(item)
        display_rows.append(["BR", row["endpoint"], "RR: P95 vs constrained MRT",
                             item["method"], covariance_label, num(row["RR"], 3),
                             "{}–{}".format(num(row["lo"], 3), num(row["hi"], 3)),
                             count(row["k"]), num(row["I2"], 1),
                             percent(row["concordance"], 1, fraction=True),
                             _source_ref(row)])

    for country, table_key, scenario, cohort in [
        ("BR", "core_primary_contrast", "primary_threeoutcome_E3_separate", "Brazil primary 2006–2020"),
        ("MX", "mx_updated_primary", "primary_complete_source_2015_2022", "Mexico complete-source 2015–2022"),
    ]:
        rows = _rows(pub, table_key, "table2", country=country, scenario=scenario,
                     contrast="y_urinary_vs_y_K80", covariance="hac14")
        for row in rows:
            if row["method"] not in {"DL", "REML_HK"}:
                continue
            item = {
                "country": country, "cohort": cohort,
                "endpoint": "urinary versus K80", "estimand": "Direct RRR: P95 versus shared median",
                "point": row["RR"], "lo": row["lo"], "hi": row["hi"],
                "method": row["method"], "covariance": row["covariance"],
                "k": row["k"], "I2": row["I2"], "concordance": row["concordance"],
                "source_file": row["source_file"], "source_row": row["source_row"],
            }
            output.append(item)
            display_rows.append([country, cohort, "Direct RRR: urinary vs K80, P95 vs shared median",
                                 row["method"], row["covariance"], num(row["RR"], 3),
                                 "{}–{}".format(num(row["lo"], 3), num(row["hi"], 3)),
                                 count(row["k"]), num(row["I2"], 1),
                                 percent(row["concordance"], 1, fraction=True),
                                 _source_ref(row)])
    if primary_covariance and "HAC14" not in primary_covariance.upper():
        raise ValueError("Final metadata primary covariance description does not identify HAC14")
    _write_csv(tables_dir / "table2.csv", output)
    caption = (
        "Table 2 | Main effects and direct shared-reference country contrasts. "
        "All intervals are 95% confidence intervals. Main effects are endpoint-specific "
        "P95-versus-constrained-MRT rate ratios; direct RRRs compare diagnoses at the same "
        "unit-specific P95 and shared median. k denotes geographic units."
    )
    _write_table(tables_dir / "table2.md", 2, caption,
                 ["Country", "Endpoint/cohort", "Estimand", "Pooling", "Covariance", "Estimate",
                  "95% CI", "k", "I² (%)", "Direction >0 (%)", "Source row"], display_rows)
    return output, caption


def _original8_label(scenario: str) -> str:
    parts = scenario.split("_")
    exposure = parts[2]
    adjustment = "shared calendar" if parts[3] == "shared" else "outcome-specific year-month FE"
    period = "2006–2020" if parts[4] == "full" else "2008–2015"
    return "{}; {}; {}".format(exposure, adjustment, period)


def _table3(pub: FrozenPublication, tables_dir: Path) -> tuple[list[dict], str]:
    output = []
    display_rows = []
    base = {"contrast": "y_urinary_vs_y_K80", "covariance": "hac14", "method": "DL"}
    candidates = pub.rows("core_robustness", country="BR", **base)
    grid = [r for r in candidates if r["scenario"].startswith("sensitivity_original8_")]
    if len(grid) != 8:
        raise ValueError(f"Table 3 needs eight original-grid rows; found {len(grid)}")
    grid = _record_claims(pub, "core_robustness", grid,
                          {"country": "BR", "scenario_prefix": "sensitivity_original8_", **base},
                          "table3")
    for row in grid:
        label = _original8_label(row["scenario"])
        item = {
            "family": "Original eight-cell grid", "cohort": "Brazil primary 2006–2020",
            "specification": label, "country": row["country"], "scenario": row["scenario"],
            "RRR": row["RR"], "lo": row["lo"], "hi": row["hi"], "k": row["k"],
            "I2": row["I2"], "concordance": row["concordance"],
            "method": row["method"], "covariance": row["covariance"],
            "source_file": row["source_file"], "source_row": row["source_row"],
        }
        output.append(item)
        display_rows.append(["Original eight-cell grid", "BR primary 2006–2020", label,
                             num(row["RR"], 3), "{}–{}".format(num(row["lo"], 3), num(row["hi"], 3)),
                             count(row["k"]), num(row["I2"], 1),
                             percent(row["concordance"], 1, fraction=True),
                             row["method"], row["covariance"], _source_ref(row)])

    rolling_rows = []
    for country in ("BR", "MX"):
        for days in (3, 7, 14, 21):
            scenario = f"sensitivity_lag_E{days}"
            table_rows = pub.rows("core_robustness", country=country, scenario=scenario, **base)
            if len(table_rows) != 1:
                raise ValueError(f"Expected one rolling-window summary for {country}/{scenario}; found {len(table_rows)}")
            rolling_rows.extend(table_rows)
    rolling_rows = _record_claims(pub, "core_robustness", rolling_rows,
                                  {"country_in": ["BR", "MX"],
                                   "scenario_in": [f"sensitivity_lag_E{i}" for i in (3, 7, 14, 21)],
                                   **base}, "table3")
    for row in rolling_rows:
        cohort = "Brazil primary 2006–2020" if row["country"] == "BR" else "Mexico original-gap comparator"
        item = {
            "family": "Rolling exposure window", "cohort": cohort,
            "specification": row["scenario"].replace("sensitivity_lag_", ""),
            "country": row["country"], "scenario": row["scenario"],
            "RRR": row["RR"], "lo": row["lo"], "hi": row["hi"], "k": row["k"],
            "I2": row["I2"], "concordance": row["concordance"],
            "method": row["method"], "covariance": row["covariance"],
            "source_file": row["source_file"], "source_row": row["source_row"],
        }
        output.append(item)
        display_rows.append(["Rolling exposure window", cohort, item["specification"],
                             num(row["RR"], 3), "{}–{}".format(num(row["lo"], 3), num(row["hi"], 3)),
                             count(row["k"]), num(row["I2"], 1),
                             percent(row["concordance"], 1, fraction=True),
                             row["method"], row["covariance"], _source_ref(row)])

    _write_csv(tables_dir / "table3.csv", output)
    caption = (
        "Table 3 | Original eight-cell grid and rolling exposure windows. "
        "Every point is a direct urinary-versus-K80 ratio of rate ratios, not a risk difference. "
        "The Mexico rolling-window rows use the repaired original-gap comparator cohort; the "
        "preferred complete-source cohort is shown separately in Table 2. The distributed-lag "
        "cumulative contrast is a different estimand."
    )
    _write_table(tables_dir / "table3.md", 3, caption,
                 ["Analysis family", "Country/cohort", "Specification", "RRR", "95% CI", "k",
                  "I² (%)", "Direction >0 (%)", "Pooling", "Covariance", "Source row"],
                 display_rows)
    return output, caption


def _table4(pub: FrozenPublication, tables_dir: Path) -> tuple[list[dict], str]:
    specs = pub.metadata.get("publication_table4_rows")
    if not isinstance(specs, list) or not specs:
        raise ValueError("Final metadata must define non-empty publication_table4_rows")
    tags = []
    for spec in specs:
        if not isinstance(spec, dict):
            raise ValueError("Every publication_table4_rows item must be an object")
        for field in ("tag", "label", "category", "estimand"):
            if not str(spec.get(field, "")).strip():
                raise ValueError(f"Every Table 4 selector requires a non-empty {field}")
        tags.append(str(spec["tag"]))
    if len(tags) != len(set(tags)):
        raise ValueError("Table 4 selector tags must be unique")
    missing_tags = TABLE4_REQUIRED_TAGS - set(tags)
    if missing_tags:
        raise ValueError(f"Table 4 selectors missing required tags: {sorted(missing_tags)}")
    selections = pub.metadata.get("publication_selections", {})
    missing_selectors = sorted(set(tags) - set(selections))
    if missing_selectors:
        raise KeyError(f"Final metadata is missing Table 4 publication selectors: {missing_selectors}")

    constants = pub.metadata.get("presentation_constants", {})
    clinical_source_period = constants.get("Brazil_clinical_source_period_primary")
    clinical_admission_period = constants.get("Brazil_clinical_admission_period")
    for name, value in (("clinical source period", clinical_source_period),
                        ("clinical admission period", clinical_admission_period)):
        if not isinstance(value, list) or len(value) != 2:
            raise ValueError(f"Final metadata must include the Brazil clinical {name}")
    source_years = "{}–{}".format(*clinical_source_period)
    admission_years = "{}–{}".format(*clinical_admission_period)
    clinical_period_label = (
        f"clinical source {source_years}; admission dates {admission_years}"
    )
    mexico_origin_period = constants.get("Mexico_origin_dictionary_period")
    if not isinstance(mexico_origin_period, list) or len(mexico_origin_period) != 2:
        raise ValueError("Final metadata must include the Mexico origin dictionary period")
    mexico_origin_period_label = "{}–{}".format(*mexico_origin_period)
    output = []
    display_rows = []
    for spec in specs:
        tag = spec["tag"]
        row = _select_from_publication_metadata(pub, tag, "table4")
        selection = selections[tag]
        point_field = spec.get("point", selection.get("point", "RR" if "RR" in row else "RRR"))
        low_field = spec.get("low", selection.get("low", "lo"))
        high_field = spec.get("high", selection.get("high", "hi"))
        point = row.get(point_field)
        low = row.get(low_field)
        high = row.get(high_field)
        if point is None:
            raise ValueError(f"Table 4 selector {tag} point field {point_field} is missing")
        k = row.get("k", row.get("units_joint_fit_successful", row.get("units_fit", "")))
        comparison = spec["label"]
        estimand = spec["estimand"]
        if tag == "br_strict":
            where = selection.get("where", {})
            if where.get("endpoint") != "strict_N20_only_vs_K80":
                raise ValueError("The Brazil strict clinical row must select strict N20-only versus K80")
            comparison = f"Brazil strict N20 vs K80 | {clinical_period_label}"
            estimand = "Strict N20/K80 RRR"
        elif tag in {"br_urgent", "br_elective"}:
            expected_group = "urgent" if tag == "br_urgent" else "elective"
            where = selection.get("where", {})
            suffix = str(selection.get("source_suffix", ""))
            contrast = str(row.get("contrast", "")).casefold()
            if (where.get("mode") != "admission_type" or where.get("group") != expected_group
                    or "interaction_group_pooled" not in suffix):
                raise ValueError(f"Table 4 selector {tag} must select the {expected_group} admission stratum")
            if "n20-only" not in contrast or "k80" not in contrast or "urinary" in contrast:
                raise ValueError(f"Table 4 {expected_group} estimate is not the strict N20-versus-K80 contrast")
            if count(k) != "411":
                raise ValueError(f"Table 4 {expected_group} stratum must use the shared joint-fit sample k=411")
            comparison = (
                f"Brazil {expected_group} | strict N20 vs K80 | {clinical_period_label}; "
                "same joint-fit sample k=411"
            )
            estimand = "Strict N20/K80 RRR; joint admission-type model"
        elif tag == "br_admission_interaction":
            where = selection.get("where", {})
            suffix = str(selection.get("source_suffix", ""))
            if (where.get("mode") != "admission_type" or point_field != "RR_ratio"
                    or "interaction_coverage" not in suffix):
                raise ValueError("The admission interaction must be the direct strict N20/K80 ratio from the joint fit")
            if count(k) != "411":
                raise ValueError("Table 4 admission interaction must use the joint-fit sample k=411")
            comparison = (
                f"Brazil urgent vs elective | strict N20/K80 | {clinical_period_label}; "
                "joint-fit sample k=411"
            )
            estimand = (
                "Ratio of strict N20/K80 RRRs from the joint fit; not a quotient of separately pooled margins"
            )
        elif tag == "mx_origin_ed":
            comparison = f"Mexico ED origin | dictionary-supported {mexico_origin_period_label}"
        elif tag == "mx_origin_interaction":
            comparison = (
                f"Mexico ED vs outpatient origin | dictionary-supported {mexico_origin_period_label}"
            )
        i2 = row.get("I2", "")
        concordance = row.get("concordance", row.get("concordance_positive", ""))
        method = row.get("method", row.get("pooling", ""))
        covariance = row.get("covariance", row.get("unit_variance", ""))
        item = {
            "category": spec["category"], "comparison": comparison,
            "estimand": estimand, "point_field": point_field,
            "point": point, "lo": "" if low is None else low,
            "hi": "" if high is None else high, "method": method,
            "covariance": covariance, "k": k, "I2": i2,
            "concordance": concordance, "source_file": row["source_file"],
            "source_row": row["source_row"],
        }
        output.append(item)
        interval_text = "{}–{}".format(num(low, 3), num(high, 3))
        display_rows.append([spec["category"], comparison, estimand,
                             num(point, 3), interval_text, method, covariance,
                             count(k) if k not in (None, "") else "NA",
                             num(i2, 1),
                             percent(concordance, 1, fraction=True)
                             if concordance not in (None, "") else "NA",
                             _source_ref(row)])
    _write_csv(tables_dir / "table4.csv", output)
    caption = (
        "Table 4 | Compact clinical, cohort, model-sensitivity and moisture checks. "
        "Rows are selected before display through the checksum-verified publication selector metadata. "
        "Intervals are the source 95% confidence limits; a blank source field remains unavailable. "
        "Brazil urgent and elective estimates use strict N20 versus K80, clinical source years "
        f"{source_years}, and admission dates {admission_years}; both use the same joint-fit sample "
        "of 411 geographic units. The formal urgent-versus-elective ratio is taken from the joint "
        "model and is not reconstructed from separately pooled margins. Case-mix comparisons use "
        "broad N20–N23 codes on a distinct eligible subset; Mexico-origin and moisture comparisons "
        "retain their own populations and definitions. Clinical subgroup reference temperatures are "
        "derived from the clinical-window calendars, separate from the fixed full-primary targets "
        "used in the independent period refits."
    )
    _write_table(tables_dir / "table4.md", 4, caption,
                 ["Category", "Comparison", "Estimand", "Estimate", "95% CI", "Pooling",
                  "Covariance/model", "k", "I² (%)", "Direction >0 (%)", "Source row"],
                 display_rows)
    return output, caption


def _write_supplement(pub: FrozenPublication, output_root: Path) -> list[Path]:
    outputs_dir = output_root / "outputs"
    source_dir = outputs_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    records_by_key = {}
    for key, table in pub.tables.items():
        rows = _supplement_rows(key, table["records"])
        selection = {}
        if key in CURVE_SCOPE_BY_KEY:
            selection = {
                "curve_scope": CURVE_SCOPE_BY_KEY[key],
                "method_in": sorted(CURVE_METHODS),
            }
        rows = _record_claims(
            pub, key, rows, selection, f"outputs/source_data/{key}.csv"
        )
        _write_csv(source_dir / f"{key}.csv", rows)
        records_by_key[key] = rows

    if {"curves_supported", "curves_all_units"} <= set(records_by_key):
        support_ids = {
            (row["source_file"], row["source_row"])
            for row in records_by_key["curves_supported"]
        }
        all_successful_ids = {
            (row["source_file"], row["source_row"])
            for row in records_by_key["curves_all_units"]
        }
        duplicates = support_ids & all_successful_ids
        if duplicates:
            raise ValueError(
                f"Supplementary curve scopes contain duplicate source rows: {sorted(duplicates)[:5]}"
            )

    clinical_keys = sorted(k for k in pub.tables if k.startswith("clinical_"))
    moisture_keys = sorted(k for k in pub.tables if k.startswith("moisture_"))
    sections = [
        ("Table S1. Full core robustness summary", ["core_robustness", "core_leave_one_state_out"]),
        ("Table S2. Mexico cohort restoration and year/period sensitivities",
         ["mx_updated_original_gap", "mx_updated_sensitivity", "mx_updated_leave_one_state_out"]),
        ("Table S3. Synchronized uncertainty", ["sync_meta", "sync_hac90"]),
        ("Table S4. Distributed-lag results", ["dlnm_original", "dlnm_mx_complete"]),
        ("Table S5. Common-supported temporal refits", ["temporal_support_pools", "temporal_support_units"]),
        ("Table S6. Supported and all-successful curve vertices", ["curves_supported", "curves_all_units"]),
        ("Table S7. Mexico admission-origin evidence", ["mx_origin_dictionary", "mx_origin_original_gap",
                                                           "mx_origin_dictionary_codes",
                                                           "mx_origin_dictionary_source_coverage"]),
        ("Table S8. Clinical endpoints, admission interactions and case-mix coverage",
         clinical_keys),
        ("Table S9. Matched-sample moisture sensitivity and support exclusions",
         moisture_keys),
        ("Table S10. Main-model boundary frequencies and MRT selection",
         ["pub_main_models", "main_mrt_selection_unit"]),
    ]
    lines = ["# Supplementary tables", "",
             "Exact registered CSV rows for every accepted result-table key are in `source_data/`; "
             "source row numbers and source filenames are retained. Markdown cells use presentation formatting only; "
             "no estimates, intervals, tests or pooled results are recomputed.", ""]
    for title, keys in sections:
        present = [key for key in keys if key in pub.tables]
        lines.extend([f"## {title}", ""])
        if present:
            if title.startswith("Table S8."):
                constants = pub.metadata.get("presentation_constants", {})
                source_period = constants.get("Brazil_clinical_source_period_primary", [])
                admission_period = constants.get("Brazil_clinical_admission_period", [])
                source_years = "{}–{}".format(*source_period) if len(source_period) == 2 else "registered years"
                admission_years = "{}–{}".format(*admission_period) if len(admission_period) == 2 else "registered dates"
                lines.extend([
                    f"Admission-type estimates and their formal interaction use strict N20 versus K80 for clinical source years {source_years} and admission dates {admission_years}; urgent and elective estimates share the joint-fit sample (k=411). The interaction is the direct joint-model ratio. Case-mix results use broad N20–N23 codes and a different eligible subset.",
                    "",
                ])
            elif title.startswith("Table S10."):
                lines.extend([
                    "The four pooled main-model rows summarize observed and simulated MRT-boundary frequencies. Unit-level MRT-selection simulation details are retained in the linked CSV only.",
                    "",
                ])
            section_rows = 0
            for key in present:
                if key == "pub_main_models" and len(records_by_key[key]) != 4:
                    raise ValueError(
                        f"Table S10 expects four pooled main-model rows; found {len(records_by_key[key])}"
                    )
                key_rows = records_by_key[key]
                if (key not in CSV_ONLY_SUPPLEMENT_KEYS
                        and len(key_rows) <= SUPPLEMENT_KEY_ROW_LIMIT
                        and section_rows + len(key_rows) > SUPPLEMENT_SECTION_ROW_LIMIT):
                    lines.extend([
                        f"#### {key}", "",
                        f"{len(key_rows):,} rows are retained in [the registered source CSV](source_data/{key}.csv) to keep this Markdown section below the compact row limit.",
                        "",
                    ])
                    continue
                key_lines, rendered_rows = _supplement_key_markdown(key, key_rows)
                lines.extend(key_lines)
                section_rows += rendered_rows
        else:
            lines.append("No accepted source key is present in the final frozen payload.")
        lines.append("")
    supplement_md = outputs_dir / "supplementary_tables.md"
    supplement_md.write_text("\n".join(lines), encoding="utf-8")
    manifest = source_dir / "source_rows_manifest.json"
    pub.save_claims(manifest)
    return [supplement_md, manifest, *sorted(source_dir.glob("*.csv"))]


def render(frozen: Path, output_root: Path) -> list[Path]:
    pub = FrozenPublication(frozen)
    required = {
        "pub_calendar", "pub_exposure", "pub_outcomes", "main_effects",
        "core_primary_contrast", "core_robustness", "mx_updated_primary",
        "mx_updated_sensitivity", "support_32c",
    }
    missing = sorted(required - set(pub.tables))
    if missing:
        raise KeyError(f"Final freeze is missing required table keys: {missing}")
    outputs_dir = output_root / "outputs"
    tables_dir = outputs_dir / "tables"
    source_data_dir = outputs_dir / "source_data"
    tables_dir.mkdir(parents=True, exist_ok=True)
    table_sets = [_table1(pub, tables_dir), _table2(pub, tables_dir),
                  _table3(pub, tables_dir), _table4(pub, tables_dir)]
    combined = ["# Main tables", ""]
    for ix, (_, caption) in enumerate(table_sets, 1):
        table_md = tables_dir / f"table{ix}.md"
        combined.extend([table_md.read_text(encoding="utf-8").rstrip(), ""])
    combined_path = outputs_dir / "tables.md"
    combined_path.write_text("\n".join(combined), encoding="utf-8")
    supplemental_paths = _write_supplement(pub, output_root)
    main_paths = [tables_dir / f"table{i}.{ext}" for i in range(1, 5)
                  for ext in ("md", "csv")]
    source_manifest = source_data_dir / "table_source_rows.json"
    main_paths.extend([combined_path, source_manifest])
    pub.save_claims(source_manifest)
    return main_paths + supplemental_paths


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True, type=Path,
                        help="checksum-verified frozen_numbers.md")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="repaired publication package root")
    return parser


def main() -> int:
    args = _parser().parse_args()
    outputs = render(args.frozen, args.outdir)
    print("Rendered {} table/supplement artifacts under {}".format(len(outputs), args.outdir))
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
