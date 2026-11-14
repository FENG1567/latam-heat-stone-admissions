"""Compose the revision tables from the verified revision freeze.

This is a presentation-only generator. Every displayed result is selected from
the frozen payload; it does not estimate, pool, derive confidence limits, or
reconstruct contrasts from separate marginal estimates.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
OLD = HERE
sys.path.insert(0, str(HERE))

from publication_common import FrozenPublication, count, md_table, num, percent, source_csv


_prior_spec = importlib.util.spec_from_file_location(
    "prior_revision_table_helpers", OLD / "render_repaired_tables.py"
)
if _prior_spec is None or _prior_spec.loader is None:
    raise ImportError("Could not load inherited supplementary table helpers")
prior_tables = importlib.util.module_from_spec(_prior_spec)
_prior_spec.loader.exec_module(prior_tables)


PRIMARY_ENDPOINTS = (
    ("y_K80", "K80"),
    ("y_urinary", "Urinary N20–N23"),
    ("urinary_vs_K80", "Direct urinary-versus-K80 RRR"),
)
TABLE4_ROWS = (
    ("br_flexible", "Calendar", "Brazil, flexible calendar and holiday", "Broad urinary/K80 RRR"),
    ("mx_yearmonth", "Calendar", "Mexico, year-month effects; preferred cohort", "Broad urinary/K80 RRR"),
    ("mx_flexible", "Calendar", "Mexico, flexible calendar and holiday; preferred cohort", "Broad urinary/K80 RRR"),
    ("mexico2023_followup", "Independent period", "Mexico 2023 refit with next-year discharge follow-up", "Broad urinary/K80 RRR; refit, not prediction"),
    ("br_urgent", "Urgency", "Brazil urgent admissions; strict N20 versus K80", "Strict N20/K80 RRR"),
    ("br_elective", "Urgency", "Brazil elective admissions; strict N20 versus K80", "Strict N20/K80 RRR"),
    ("br_admission_interaction", "Urgency interaction", "Brazil urgent versus elective", "Direct ratio of strict-code RRRs from joint fit"),
    ("br_urgency_flexible", "Urgency sensitivity", "Brazil urgency interaction with flexible calendar", "Direct ratio of strict-code RRRs"),
    ("br_urgency_no2020", "Urgency sensitivity", "Brazil urgency interaction excluding 2020; fixed clinical-period targets", "Direct ratio of strict-code RRRs"),
    ("mx_strict_N20", "Phenotype", "Mexico strict N20 versus K80", "Strict-code RRR"),
    ("mx_confirmed_N20_N21", "Phenotype", "Mexico N20–N21 versus K80", "N20–N21 coded RRR"),
    ("mx_N23", "Phenotype", "Mexico N23 versus K80", "Separate N23 RRR"),
    ("br_adult_matched", "Adult bridge", "Brazil adult broad-code, unstandardized; matched sample", "Unstandardized broad-code RRR"),
    ("br_standardized_matched", "Adult bridge", "Brazil adult broad-code, common-case-weight standardized; same sample", "Common-case-weight standardized RRR"),
    ("mx_adult_matched", "Adult bridge", "Mexico adult broad-code, unstandardized; matched sample", "Unstandardized broad-code RRR"),
    ("mx_standardized_matched", "Adult bridge", "Mexico adult broad-code, common-case-weight standardized; same sample", "Common-case-weight standardized RRR"),
    ("br_moisture_dew", "Moisture", "Brazil, dew-point adjustment", "Conditional broad urinary/K80 RRR"),
    ("br_moisture_rh", "Moisture", "Brazil, relative-humidity adjustment", "Conditional broad urinary/K80 RRR"),
    ("mx_moisture_dew", "Moisture", "Mexico, dew-point adjustment", "Conditional broad urinary/K80 RRR"),
    ("mx_moisture_rh", "Moisture", "Mexico, relative-humidity adjustment; DL", "Conditional broad urinary/K80 RRR"),
    ("mx_moisture_rh_hk", "Moisture", "Mexico, relative-humidity adjustment; REML–HK", "Conditional broad urinary/K80 RRR"),
)

SUPPLEMENT_GROUPS = (
    (11, "Source records, discharge-year audit and 2023 follow-up", ("source_audit", "mexico_followup"),
     ("source_audit", "brazil_duplicate_audit", "parser_summary", "discharge_years",
      "diagnosis_source_audit", "source_2023_2024", "source_followup_counts",
      "source_crossyear_identity")),
    (12, "Temporal models and independent-period refits", ("temporal",),
     ("temporal_revision", "temporal_models")),
    (13, "Adult bridge and strict-code phenotype controls", ("clinical_bridge",),
     ("clinical_bridge", "adult_bridge", "strict_phenotype", "code_composition")),
    (14, "Synchronized marginal estimates and uncertainty", (),
     ("marginal_sync", "marginal_independent")),
    (15, "Supported curves, pooling weights and estimator diagnostics", (),
     ("curves", "curve_weight_diagnostics", "curve_identical_example",
      "cache_sample_acceptance", "curve_psd")),
    (16, "Unit-specific type-7 target-temperature distributions", (), ("target_type7",)),
)


def _role_key(pub: FrozenPublication, name: str, *, required: bool = True) -> str | None:
    roles = pub.metadata.get("revision", {}).get("tables_by_role", {})
    value = roles.get(name)
    if value is None:
        if required:
            raise KeyError(f"Revision freeze has no required tables_by_role alias {name!r}")
        return None
    if isinstance(value, dict):
        value = value.get("table_key", value.get("key"))
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            raise ValueError(f"Role {name!r} must name exactly one frozen table, got {value!r}")
        value = value[0]
    if not isinstance(value, str) or not value:
        raise ValueError(f"Invalid frozen table alias for role {name!r}: {value!r}")
    if value not in pub.tables:
        raise KeyError(f"Role {name!r} points to absent frozen table {value!r}")
    return value


def _rows(pub: FrozenPublication, key: str, artifact: str, **where) -> list[dict]:
    if key not in pub.tables:
        raise KeyError(f"Required frozen table is missing: {key}")
    rows = pub.rows(key, **where)
    if not rows:
        raise ValueError(f"Frozen selection returned no rows: {key} {where}")
    return prior_tables._record_claims(pub, key, rows, where, artifact)


def _one(pub: FrozenPublication, key: str, artifact: str, **where) -> dict:
    rows = _rows(pub, key, artifact, **where)
    if len(rows) != 1:
        raise ValueError(f"Expected one frozen row: {key} {where}; found {len(rows)}")
    return rows[0]


def _select(pub: FrozenPublication, tag: str) -> tuple[dict, dict]:
    revision = pub.metadata.get("revision", {})
    spec = revision.get("display_selections", {}).get(tag)
    if spec is not None:
        key = spec["table_key"]
        where = spec.get("where", {})
        row = pub.one(key, **where)
        aliases = {name: row[source] for name, source in spec.get("fields", {}).items()}
        return {**row, **aliases}, spec

    spec = pub.metadata.get("publication_selections", {}).get(tag)
    if spec is None:
        raise KeyError(f"No revision display selector or inherited selector for {tag!r}")
    rows = pub.source_rows(spec["source_suffix"], **spec["where"])
    if len(rows) != 1:
        raise ValueError(f"Expected one frozen source row for {tag}: found {len(rows)}")
    row = rows[0]
    pub.claims.append({"selector_tag": tag, "selection": spec["where"],
                       "source_file": row["source_file"], "source_row": row["source_row"]})
    return row, spec


def _source_ref(row: dict) -> str:
    return f"{row['source_file']}:{row['source_row']}"


def _shown_count(value) -> str:
    if value in (None, ""):
        return "NA"
    try:
        return count(value)
    except (TypeError, ValueError):
        return num(value, 2)


def _first(row: dict, fields: tuple[str, ...]):
    for field in fields:
        value = row.get(field)
        if value not in (None, ""):
            return value
    return None


def _table1(pub: FrozenPublication) -> tuple[str, list[dict]]:
    calendar_specs = (
        ("BR_calendar", "BR", "calendar_br", "Brazil, primary 2006–2020"),
        ("MX_2015_2022", "MX", "calendar_mx", "Mexico, preferred complete-source 2015–2022"),
    )
    rows_out = []
    display = []
    outcome_key = "pub_outcomes"
    for dataset, country, tag, label in calendar_specs:
        calendar, _ = _select(pub, tag)
        exposure = _one(pub, "pub_exposure", "table1", dataset=dataset, country=country)
        outcome_rows = pub.rows(outcome_key, dataset=dataset, country=country, row_scope="all_health")
        outcome_rows = [row for row in outcome_rows if row.get("outcome") in {"K80", "K81", "urinary"}]
        outcome_rows = prior_tables._record_claims(
            pub, outcome_key, outcome_rows,
            {"dataset": dataset, "country": country, "row_scope": "all_health",
             "outcome_in": ["K80", "K81", "urinary"]}, "table1")
        outcome_by_name = {row["outcome"]: row for row in outcome_rows}
        if set(outcome_by_name) != {"K80", "K81", "urinary"}:
            raise ValueError(f"Table 1 outcome totals incomplete for {dataset}")
        target = _one(pub, _role_key(pub, "target_type7"), "table1 target",
                      country=country, variable="delta_p95_minus_median")
        source_rows = [calendar, exposure, *outcome_rows, target]
        refs = "; ".join(_source_ref(row) for row in source_rows)
        item = {
            "country_period": label,
            "geographic_units": calendar.get("unit_count", "NA"),
            "calendar_unit_days": calendar.get("row_count", "NA"),
            "calendar_dates": calendar.get("distinct_date_count", "NA"),
            "eligible_unit_days": calendar.get("eligible_rows", "NA"),
            "missing_health_unit_days": calendar.get("missing_health_rows", "NA"),
            "zero_admission_calendar_days": calendar.get("all3_zero_admission_days_among_health", "NA"),
            "zero_admission_eligible_days": calendar.get("all3_zero_admission_days_among_eligible", "NA"),
            "mean_E3_C": exposure.get("mean_C", "NA"),
            "median_E3_C": exposure.get("median_C", "NA"),
            "P05_E3_C": exposure.get("p05_C", "NA"),
            "P95_E3_C": exposure.get("p95_C", "NA"),
            "target_increment_median_C": target["median"],
            "target_increment_q25_C": target["q25"],
            "target_increment_q75_C": target["q75"],
            "K80_admission_records": outcome_by_name["K80"]["total_events"],
            "K81_admission_records": outcome_by_name["K81"]["total_events"],
            "urinary_N20_N23_admission_records": outcome_by_name["urinary"]["total_events"],
            "source_rows": refs,
        }
        rows_out.append(item)
    by_country = {item["country_period"].split(",", 1)[0]: item for item in rows_out}
    brazil = by_country["Brazil"]
    mexico = by_country["Mexico"]
    metrics = (
        ("Geographic units", "geographic_units", "count"),
        ("Calendar unit-days", "calendar_unit_days", "count"),
        ("Calendar dates", "calendar_dates", "count"),
        ("Exposure-eligible unit-days", "eligible_unit_days", "count"),
        ("Missing-health unit-days", "missing_health_unit_days", "count"),
        ("Zero-admission unit-days (health calendar)", "zero_admission_calendar_days", "count"),
        ("Zero-admission eligible unit-days", "zero_admission_eligible_days", "count"),
        ("Mean E3 temperature (°C)", "mean_E3_C", "temperature"),
        ("Median E3 temperature (°C)", "median_E3_C", "temperature"),
        ("P05 E3 temperature (°C)", "P05_E3_C", "temperature"),
        ("P95 E3 temperature (°C)", "P95_E3_C", "temperature"),
        ("Unit-specific P95-minus-median increment, median (IQR; °C)", "target_increment_median_C", "target"),
        ("K80 admission records", "K80_admission_records", "count"),
        ("K81 admission records", "K81_admission_records", "count"),
        ("Urinary N20–N23 admission records", "urinary_N20_N23_admission_records", "count"),
    )
    display = []
    for label, field, value_type in metrics:
        row_values = []
        for item in (brazil, mexico):
            if value_type == "count":
                value = _shown_count(item[field])
            elif value_type == "temperature":
                value = num(item[field], 2)
            else:
                value = (f"{num(item['target_increment_median_C'], 2)} "
                         f"({num(item['target_increment_q25_C'], 2)}–{num(item['target_increment_q75_C'], 2)})")
            row_values.append(value)
        display.append([label, *row_values])
    caption = (
        "Coverage, diagnostic counts and exposure descriptors for the two primary cohorts. "
        "Calendar and exposure-eligible unit-days are distinct; missing health unit-days and zero-admission "
        "eligible days are shown separately. Admission totals are full-health-calendar records, not unique people; "
        "they include all calendar unit-days rather than only exposure-eligible days. E3 percentiles "
        "describe pooled eligible unit-days. The target increment summarizes the already-defined unit-specific "
        "P95-minus-median contrasts; its across-unit median and quartiles use empirical type-7 quantiles. "
        "Full source-row provenance is in the companion CSV and source-claims file."
    )
    table = "## Table 1\n\n" + caption + "\n\n" + md_table(
        ["Characteristic", "Brazil, primary 2006–2020", "Mexico, preferred 2015–2022"], display) + "\n"
    return table, rows_out


def _table2(pub: FrozenPublication) -> tuple[str, list[dict]]:
    key = _role_key(pub, "marginal_sync")
    display = []
    source_rows = []
    main_rows = []
    for country in ("BR", "MX"):
        for endpoint, label in PRIMARY_ENDPOINTS:
            row = _one(pub, key, "table2 primary",
                       country=country, endpoint=endpoint, bandwidth_days="14",
                       inference="calendar_synchronous_plus_DL_tau_component")
            item = {"panel": "Aligned primary marginal and direct results", "country": country,
                    "endpoint": label, "estimate": row["RR"], "lo": row["lo"], "hi": row["hi"],
                    "pooling": row.get("method", "DL"),
                    "inference": row.get("inference", "calendar_synchronous_plus_DL_tau_component"),
                    "k": row.get("k", "NA"), "I2": row.get("I2", "NA"),
                    "concordance": row.get("concordance", "NA"),
                    "source_file": row["source_file"], "source_row": row["source_row"]}
            main_rows.append(item)
            display.append([f"{country} | {label}",
                            f"{num(row['RR'], 3)} ({num(row['lo'], 3)}–{num(row['hi'], 3)})",
                            _shown_count(item["k"]),
                            num(item["I2"], 1) if item["I2"] != "NA" else "NA",
                            percent(item["concordance"], 1, fraction=True)
                            if item["concordance"] != "NA" else "NA"])
            source_rows.append(row)

    burden_key = _role_key(pub, "burden_country", required=False) or "burden_country"
    absolute_display = []
    absolute_rows = []
    for endpoint, label in (("y_K80", "K80"), ("y_urinary", "Urinary N20–N23"),
                            ("urinary_minus_K80", "Urinary minus K80")):
        row = _one(pub, burden_key, "table2 absolute scenario", country="BR", endpoint=endpoint)
        estimate = row["national_point_sum_of_unit_means"]
        item = {"panel": "Additional Brazilian absolute daily scenario", "country": "BR",
                "endpoint": label, "estimate": estimate, "lo": row["lo"], "hi": row["hi"],
                "k": row.get("included_region_count", "NA"), "source_file": row["source_file"],
                "source_row": row["source_row"]}
        absolute_rows.append(item)
        absolute_display.append([label,
                                f"{num(estimate, 1)} ({num(row['lo'], 1)}–{num(row['hi'], 1)})",
                                _shown_count(row.get("included_region_count", "NA"))])
        source_rows.append(row)

    caption = (
        "Aligned country-specific marginal RRs and directly estimated urinary-versus-K80 RRRs for the "
        "same primary joint-fit unit sets, using each unit’s P95 versus its shared median. Marginal RRs "
        "retain endpoint-specific pooling weights; dividing them does not recover the directly pooled RRR. "
        "Intervals are the frozen calendar-synchronized Bartlett HAC14 inference with the modeled DL "
        "between-unit component, conditional on observed weights and heterogeneity. k, I² and positive-unit "
        "direction are frozen unit-level summaries. Full inference fields and source-row provenance are in the "
        "companion CSV and source-claims file."
    )
    table = "## Table 2\n\n" + caption + "\n\n" + md_table(
        ["Country / endpoint", "RR or direct RRR (95% CI)", "k", "I² (%)", "Direction >0 (%)"], display)
    table += ("\n\n**Additional panel: standardized Brazilian daily absolute-count scenario.** "
              "These frozen panel summaries compare simultaneous local P95 settings with local medians; "
              "they are not observed national heat-event totals.\n\n")
    table += md_table(["Endpoint", "Daily count difference (95% CI)", "Included regions"], absolute_display) + "\n"
    return table, main_rows + absolute_rows


def _grid_label(scenario: str) -> str:
    parts = scenario.split("_")
    if len(parts) < 5:
        return scenario
    exposure = parts[2]
    adjustment = "shared calendar" if parts[3] == "shared" else "outcome-specific year-month fixed effects"
    period = "2006–2020" if parts[4] == "full" else "2008–2015"
    return f"{exposure}; {adjustment}; {period}"


def _table3(pub: FrozenPublication) -> tuple[str, list[dict]]:
    key = "core_robustness"
    base = {"contrast": "y_urinary_vs_y_K80", "covariance": "hac14", "method": "DL"}
    grid = [row for row in pub.rows(key, country="BR", **base)
            if str(row.get("scenario", "")).startswith("sensitivity_original8_")]
    if len(grid) != 8:
        raise ValueError(f"Table 3 requires the retained eight-cell grid; found {len(grid)}")
    grid = prior_tables._record_claims(
        pub, key, grid, {"country": "BR", "scenario_prefix": "sensitivity_original8_", **base}, "table3")
    selected = [(row, "Original eight-cell grid", "Brazil 2008–2015" if row['scenario'].endswith('_sub') else "Brazil 2006–2020", _grid_label(row["scenario"]))
                for row in grid]
    rolling = []
    for country in ("BR", "MX"):
        for days in (3, 7, 14, 21):
            scenario = f"sensitivity_lag_E{days}"
            rows = pub.rows(key, country=country, scenario=scenario, **base)
            if len(rows) != 1:
                raise ValueError(f"Expected one rolling-window row for {country}/{scenario}; found {len(rows)}")
            rolling.extend(rows)
    rolling = prior_tables._record_claims(
        pub, key, rolling,
        {"country_in": ["BR", "MX"], "scenario_in": [f"sensitivity_lag_E{x}" for x in (3, 7, 14, 21)], **base},
        "table3")
    for row in rolling:
        cohort = ("Brazil primary 2006–2020" if row["country"] == "BR"
                  else "Mexico original-gap cohort, historical comparator")
        selected.append((row, "Rolling exposure-window sensitivity", cohort,
                         str(row["scenario"]).replace("sensitivity_lag_", "")))

    display = []
    rows_out = []
    for row, family, cohort, specification in selected:
        item = {"analysis_family": family, "country_cohort": cohort,
                "specification": specification, "RRR": row["RR"], "lo": row["lo"], "hi": row["hi"],
                "k": row["k"], "I2": row["I2"], "concordance": row["concordance"],
                "pooling": row["method"], "covariance": row["covariance"],
                "source_file": row["source_file"], "source_row": row["source_row"]}
        rows_out.append(item)
        if family == "Original eight-cell grid":
            scenario = f"Eight-cell grid | BR primary | {specification}"
        else:
            scenario = f"Rolling {specification} | {cohort}"
        display.append([scenario,
                        f"{num(row['RR'], 3)} ({num(row['lo'], 3)}–{num(row['hi'], 3)})",
                        _shown_count(row["k"]), num(row["I2"], 1),
                        percent(row["concordance"], 1, fraction=True)])
    caption = (
        "The original eight-cell sensitivity grid and rolling exposure windows. Every estimate is a direct "
        "urinary-versus-K80 ratio of rate ratios. Mexican rolling-window estimates retain the historical "
        "original-gap cohort label; the preferred complete-source cohort is reported in Table 2. These are "
        "archived DL/HAC14 results from the verified revision freeze; the distributed-lag cumulative contrast "
        "is a different estimand. Full inference metadata and source-row provenance are in the companion CSV "
        "and source-claims file."
    )
    table = "## Table 3\n\n" + caption + "\n\n" + md_table(
        ["Scenario / cohort", "Direct RRR (95% CI)", "k", "I² (%)", "Direction >0 (%)"], display) + "\n"
    return table, rows_out


def _coverage_status(row: dict) -> str:
    attempted = _first(row, ("units_attempted", "attempted_units", "attempted_regions",
                             "attempted_regions_success_plus_failure", "units_total", "n_attempted", "states_attempted"))
    successful = _first(row, ("units_joint_fit_successful", "units_fit", "successful_units",
                              "successful_regions", "k", "n_success", "states_fit"))
    failed = _first(row, ("units_failed", "units_fit_failed", "failed_units", "failed_regions", "n_failed", "fit_failures"))
    ineligible = _first(row, ("states_ineligible", "units_ineligible"))
    same_status = _first(row, ("same_sample_status", "same_sample", "common_success_status"))
    common_count = _first(row, ("common_success_units", "units_common_success", "common_success_k"))
    values = []
    if row.get("selected_subset_full_panel_units") is not None:
        values.append(f"selected subset {_shown_count(row.get('units_total',row.get('k')))} of {_shown_count(row['selected_subset_full_panel_units'])} panel units")
    if successful is not None:
        values.append(f"fit { _shown_count(successful) }")
    if attempted is not None:
        values.append(f"attempted { _shown_count(attempted) }")
    if failed is not None:
        values.append(f"failed { _shown_count(failed) }")
    if ineligible is not None:
        values.append(f"ineligible {_shown_count(ineligible)}")
    if same_status is not None:
        if isinstance(same_status, bool) or str(same_status).casefold() in {"true", "false"}:
            same_status = "yes" if str(same_status).casefold() == "true" else "no"
        values.append(f"same-sample status {same_status}")
    if common_count is not None:
        values.append(f"common-success units {_shown_count(common_count)}")
    return "; ".join(values) if values else "See revision audit supplement"


def _table4(pub: FrozenPublication) -> tuple[str, list[dict]]:
    display = []
    rows_out = []
    for tag, category, comparison, estimand in TABLE4_ROWS:
        row, spec = _select(pub, tag)
        if tag in ("br_flexible","mx_yearmonth","mx_flexible","br_urgency_flexible","br_urgency_no2020"):
            calendar,_=_select(pub,"calendar_mx" if tag.startswith("mx_") else "calendar_br")
            row={**row,"units_attempted":calendar["unit_count"]}
        if tag == "mexico2023_followup":
            coverage=pub.source_one("mexico_followup_count_and_attempt_audit.csv",scenario_tag="mexico2023_followup")
            row={**row,**{key:coverage[key] for key in ("states_attempted","states_fit","states_ineligible","fit_failures")}}
        if tag in ("br_adult_matched","br_standardized_matched","mx_adult_matched","mx_standardized_matched"):
            calendar,_=_select(pub,"calendar_br" if tag.startswith("br_") else "calendar_mx")
            row={**row,"selected_subset_full_panel_units":calendar["unit_count"]}
        point_field = spec.get("point") or ("RR" if "RR" in row else "RRR" if "RRR" in row else "RR_ratio")
        low_field = spec.get("low", "lo")
        high_field = spec.get("high", "hi")
        if point_field not in row:
            raise KeyError(f"Table 4 selector {tag!r} has no point field {point_field!r}")
        estimate, low, high = row[point_field], row.get(low_field), row.get(high_field)
        item = {
            "tag": tag, "category": category, "comparison": comparison, "estimand": estimand,
            "point_field": point_field, "estimate": estimate,
            "lo": "" if low is None else low, "hi": "" if high is None else high,
            "pooling": row.get("method", row.get("pooling", "")),
            "covariance_model": row.get("covariance", row.get("inference", row.get("model", ""))),
            "k": _first(row, ("k", "units_joint_fit_successful", "units_fit")),
            "I2": row.get("I2", ""),
            "concordance": row.get("concordance", row.get("concordance_positive", "")),
            "coverage_and_fit_status": _coverage_status(row),
            "source_file": row["source_file"], "source_row": row["source_row"],
        }
        rows_out.append(item)
        display.append([f"{category}: {comparison}",
                        f"{num(estimate, 3)} ({num(low, 3)}–{num(high, 3)})",
                        _shown_count(item["k"]),
                        num(item["I2"], 1) if item["I2"] not in (None, "") else "NA",
                        percent(item["concordance"], 1, fraction=True)
                        if item["concordance"] not in (None, "") else "NA",
                        item["coverage_and_fit_status"]])
    caption = (
        "Targeted calendar, phenotype, adult-bridge and urgency comparisons, with concise moisture checks. "
        "Calendar, phenotype, adult and urgency rows use unit-HAC14 contrasts with independent-unit DL pooling; moisture rows retain their labeled source inference. These are distinct from the synchronized primary intervals in Table 2. "
        "Rows use the frozen selectors and source intervals. The urgency interaction is selected from the "
        "joint model and is not obtained by dividing independently pooled strata. Fit attempts, exclusions, "
        "failures and common-success denominators are shown when present in the selected record; detailed "
        "attempt and same-sample audits are retained in the linked supplementary revision tables. The 2023 "
        "estimate is a period refit, not prediction validation. Full scenario definitions, inference fields "
        "and source-row provenance are in the companion CSV and source-claims file."
    )
    table = "## Table 4\n\n" + caption + "\n\n" + md_table(
        ["Scenario / contrast", "Effect (95% CI)", "k", "I² (%)", "Direction >0 (%)",
         "Fit / sample status"], display) + "\n"
    return table, rows_out


def _revision_supplement(pub: FrozenPublication) -> tuple[list[str], list[str], list[str]]:
    lines = []
    skipped = []
    included = []
    revision = pub.metadata.get("revision", {})
    new_keys = set(revision.get("new_table_keys", []))
    module_keys = {
        module: {key for key in new_keys if key.startswith(f"revision_{module}_")}
        for _, _, modules, _ in SUPPLEMENT_GROUPS for module in modules
    }
    for number, title, modules, roles in SUPPLEMENT_GROUPS:
        keys = set().union(*(module_keys.get(module, set()) for module in modules)) if modules else set()
        for role in roles:
            key = _role_key(pub, role, required=False)
            if key is not None:
                keys.add(key)
        keys = sorted(key for key in keys if key in new_keys)
        if not keys:
            skipped.append(f"S{number}: no accepted revision table keys for {title}")
            continue
        included.append(f"S{number}")
        lines.extend([f"## Table S{number}. {title}", ""])
        if number == 11:
            lines.extend(["The 2023 admission-period result uses available next-year discharge follow-up. Source-year and discharge-year overlaps are cohort diagnostics, not proof of unique patients or clinical episodes.", ""])
        elif number == 12:
            lines.extend(["Calendar controls and period refits retain each model's registered cohort, target, common-success set and failure accounting. Independent period estimates are refits, not prediction validation.", ""])
        elif number == 13:
            lines.extend(["Adult unstandardized and common-case-weight summaries retain their matched successful units and coverage. Strict-code controls are shown separately; these analyses do not estimate population-standardized incidence.", ""])
        elif number == 14:
            lines.extend(["Synchronized intervals are conditional on observed pooling weights and the modeled between-unit component. Direct RRRs are preserved separately from the marginal RRs.", ""])
        elif number == 15:
            lines.extend(["Pointwise supported curves and pooling diagnostics use frozen estimates; they do not imply a biological threshold. Numerical-stability warnings and estimator sensitivity are retained.", ""])
        elif number == 16:
            lines.extend(["These rows are the registered unit-specific type-7 target distributions used by the revision analyses.", ""])
        for key in keys:
            rows = pub.tables[key]["records"]
            lines.extend([f"#### {key}", ""])
            lines.extend([f"[Complete labeled source table](source_data/{key}.csv). All source rows and provenance are retained in this CSV.", ""])
            if len(rows)>200:
                lines.extend(["The detailed unit-level or temperature-level rows are supplied in the complete labeled CSV above; no rows are excluded from that file.", ""])
                continue
            # New schemas must retain their actual identifiers. The inherited
            # formatter recognized old field names and silently dropped variant,
            # baseline role, metric and variable labels from revision tables.
            columns=list(dict.fromkeys(c for r in rows for c in r if c not in ('source_file','source_row')
                and not c.startswith('upstream_') and c not in ('formula','basis_json','basis_spec_json','B_json','coef_json','coefficients_json')))
            def shown(value):
                if value is None or str(value) in ('','NA','NaN','nan'):return 'NA'
                try:
                    x=float(value)
                    if x==int(x):return str(int(x))
                    return format(x,'.5g')
                except (ValueError,TypeError,OverflowError):return str(value)
            lines.extend([md_table(columns,[[shown(r.get(c,'')) for c in columns] for r in rows]),''])
    return lines, skipped, included


def compose(frozen, output_root) -> dict:
    output_root = Path(output_root)
    outputs = output_root / "outputs"
    pub = FrozenPublication(frozen)
    revision = pub.metadata.get("revision")
    if not isinstance(revision, dict):
        raise ValueError("Input freeze is not the accepted revision freeze")
    _role_key(pub, "marginal_sync")
    _role_key(pub, "target_type7")
    _role_key(pub, "curves")

    outputs.mkdir(parents=True, exist_ok=True)
    table_dir = outputs / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    # Keep S1–S10 bound to inherited freeze tables. Accepted revision tables
    # are exported and appended only under the six revision supplement headings.
    all_tables = pub.tables
    new_keys = set(revision.get("new_table_keys", []))
    if not new_keys <= set(all_tables):
        raise KeyError(f"Revision freeze omits accepted table keys: {sorted(new_keys - set(all_tables))}")
    pub.tables = {key: table for key, table in all_tables.items() if key not in new_keys}
    prior_tables._write_supplement(pub, output_root)
    pub.tables = all_tables
    supplement_path = outputs / "supplementary_tables.md"
    inherited_supplement = supplement_path.read_text(encoding="utf-8")
    inherited_supplement = inherited_supplement.replace(
        "The four pooled main-model rows summarize observed and simulated MRT-boundary frequencies. "
        "Unit-level MRT-selection simulation details are retained in the linked CSV only.",
        "The MRT summaries are secondary exploratory selected-reference results. Their intervals "
        "condition on the selected P10–P60 minimum; unit-coefficient simulations are not national "
        "selection-aware inference. Unit-level selection-simulation details are retained in the linked CSV.")
    if "secondary exploratory selected-reference results" not in inherited_supplement:
        raise ValueError("Inherited Table S10 must explicitly qualify selected-MRT results as exploratory")
    supplement_path.write_text(inherited_supplement, encoding="utf-8")
    source_dir = outputs / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    for key in sorted(new_keys):
        rows = pub.tables[key]["records"]
        source_csv(rows, source_dir / f"{key}.csv")
        prior_tables._record_claims(pub, key, rows, {}, f"outputs/source_data/{key}.csv")

    table_results = [_table1(pub), _table2(pub), _table3(pub), _table4(pub)]
    combined = ["# Main tables", ""]
    csv_sources = []
    for index, (markdown, rows) in enumerate(table_results, start=1):
        combined.extend([markdown.rstrip(), ""])
        csv_path = table_dir / f"table{index}.csv"
        source_csv(rows, csv_path)
        csv_sources.append(str(csv_path))
    tables_path = outputs / "tables.md"
    tables_path.write_text("\n".join(combined).rstrip() + "\n", encoding="utf-8")

    supplement_lines, skipped, included = _revision_supplement(pub)
    previous = supplement_path.read_text(encoding="utf-8")
    supplement_path.write_text(previous.rstrip() + "\n\n" + "\n".join(supplement_lines).rstrip() + "\n",
                                encoding="utf-8")

    claims_path = outputs / "source_data" / "table_source_rows_revision.json"
    claims_path.parent.mkdir(parents=True, exist_ok=True)
    pub.save_claims(claims_path)
    manifest = {
        "freeze": pub.path.name,
        "tables": 4,
        "main_markdown": "outputs/tables.md",
        "supplementary_markdown": "outputs/supplementary_tables.md",
        "main_csv": [f"outputs/tables/table{i}.csv" for i in range(1, 5)],
        "required_roles": {name: _role_key(pub, name) for name in ("marginal_sync", "target_type7", "curves")},
        "supplement_sections_included": included,
        "optional_roles_skipped": skipped,
        "source_claim_count": len(pub.claims),
        "numbers": "Selected from the verified revision freeze; no effects or intervals recomputed.",
    }
    manifest_path = outputs / "revision_tables_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))
    return manifest


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    compose(args.frozen, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
