#!/usr/bin/env python3
"""Render the four repaired publication figures from a verified frozen payload.

No source CSVs are read here. Every numerical observation is selected from the
checksum-verified freeze by table key and row identity. Run this only after the
root declares the final freeze ready.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from publication_common import FrozenPublication, count, num, percent, source_csv


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BLUE = "#1f4e79"
RED = "#c0392b"
NEUTRAL = "#5b6570"
GRID = "#e8ebee"


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "axes.axisbelow": True,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "savefig.facecolor": "white",
    })


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


def _save_source_rows(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    source_csv(rows, path)


def _save_figure(fig, path_without_extension: Path) -> list[Path]:
    path_without_extension.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_without_extension.with_suffix(".png"), dpi=300,
                bbox_inches="tight", facecolor="white")
    fig.savefig(path_without_extension.with_suffix(".svg"),
                bbox_inches="tight", facecolor="white")
    fig.savefig(path_without_extension.with_suffix(".pdf"),
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [path_without_extension.with_suffix(ext) for ext in (".png", ".svg", ".pdf")]


def _label(ax, label: str) -> None:
    ax.text(-0.075, 1.02, label, transform=ax.transAxes,
            fontsize=10, fontweight="bold", ha="left", va="bottom")


def _year_position(value: str) -> float:
    current = date.fromisoformat(value)
    first = date(current.year, 1, 1)
    following = date(current.year + 1, 1, 1)
    return current.year + (current - first).days / (following - first).days


def _inclusive_period(start_year: int, end_year: int) -> tuple[float, float]:
    return float(start_year), float(end_year + 1)


def _figure1(pub: FrozenPublication, out: Path, supplement: Path) -> tuple[list[Path], str]:
    artifact = "figure1"
    cal_br = _one(pub, "pub_calendar", artifact, dataset="BR_calendar", country="BR")
    cal_mx = _one(pub, "pub_calendar", artifact, dataset="MX_2015_2022", country="MX")
    cal_mx23 = _one(pub, "pub_calendar", artifact, dataset="MX_2023_extension", country="MX")
    _save_source_rows([cal_br, cal_mx, cal_mx23], supplement / "figure1_calendar_source.csv")

    br_all = pub.rows("pub_outcomes", dataset="BR_calendar", country="BR",
                      row_scope="all_health")
    mx_all = pub.rows("pub_outcomes", dataset="MX_2015_2022", country="MX",
                      row_scope="all_health")
    diagnostic_rows = [r for r in br_all + mx_all
                       if r["outcome"] in {"K80", "K81", "urinary", "sum3"}]
    br_outcomes = [r for r in diagnostic_rows if r["dataset"] == "BR_calendar"]
    mx_outcomes = [r for r in diagnostic_rows if r["dataset"] == "MX_2015_2022"]
    _record_claims(pub, "pub_outcomes", br_outcomes,
                   {"dataset": "BR_calendar", "country": "BR", "row_scope": "all_health",
                    "outcome_in": ["K80", "K81", "urinary", "sum3"]}, artifact)
    _record_claims(pub, "pub_outcomes", mx_outcomes,
                   {"dataset": "MX_2015_2022", "country": "MX", "row_scope": "all_health",
                    "outcome_in": ["K80", "K81", "urinary", "sum3"]}, artifact)
    _save_source_rows(diagnostic_rows, supplement / "figure1_diagnostic_admission_totals_source.csv")

    constants = pub.metadata.get("presentation_constants", {})
    clinical_period = constants.get("Brazil_clinical_admission_period")
    if not isinstance(clinical_period, list) or len(clinical_period) != 2:
        raise ValueError("Final metadata must include the Brazil clinical admission period")
    clinical_source_period = constants.get("Brazil_clinical_source_period_primary")
    if not isinstance(clinical_source_period, list) or len(clinical_source_period) != 2:
        raise ValueError("Final metadata must include the Brazil clinical source-year period")

    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.2, 1.6), hspace=0.34, wspace=0.25)
    timeline = fig.add_subplot(gs[0, :])
    axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    timeline.set_axisbelow(True)
    timeline.set_ylim(-0.65, 3.75)
    timeline.set_xlim(2005.5, 2024.1)
    timeline.set_yticks([3, 2.25, 1, 0])
    timeline.set_yticklabels([
        "Brazil primary", "Brazil clinical admissions", "Mexico preferred", "Mexico extension",
    ])
    timeline.set_xticks([2006, 2010, 2014, 2018, 2022, 2023])
    timeline.set_xticklabels(["2006", "2010", "2014", "2018", "2022", "2023"])
    timeline.set_xlabel("Admission year")
    timeline.grid(axis="x", color=GRID, linewidth=0.7)
    _label(timeline, "a")

    def period_bar(row: dict, y: float, color: str, label: str, linestyle="-") -> None:
        x0 = _year_position(row["first_date"])
        x1 = _year_position(row["last_date"]) + 1 / 366
        timeline.plot([x0, x1], [y, y], color=color, linewidth=7,
                      solid_capstyle="butt", linestyle=linestyle)
        timeline.text(x0, y + 0.18, label, fontsize=7.2, color=color, ha="left", va="bottom")

    period_bar(cal_br, 3, BLUE, "Brazil 2006–2020")
    period_bar(cal_mx, 1, BLUE, "Mexico 2015–2022")
    period_bar(cal_mx23, 0, RED, "Independent 2023 extension")
    # The 2019 band denotes reconstructed source coverage within the preferred Mexico period.
    timeline.fill_between([2019, 2020], 0.66, 1.34, color=RED, alpha=0.10, zorder=0)
    timeline.text(2019.5, 0.54, "2019 restored", ha="center", va="top", fontsize=7.2,
                  color=RED)
    c0, c1 = _inclusive_period(int(clinical_period[0]), int(clinical_period[1]))
    clinical_y = 2.25
    timeline.plot([c0, c1], [clinical_y, clinical_y], color=RED, linewidth=2.0,
                  linestyle=(0, (2, 2)))
    timeline.text(c0, clinical_y + 0.16,
                  "Admission dates {}–{}".format(*clinical_period),
                  fontsize=7.2, color=RED, ha="left", va="bottom")

    outcome_colors = {"K80": BLUE, "K81": BLUE, "urinary": RED, "sum3": NEUTRAL}
    outcome_labels = {
        "K80": "K80", "K81": "K81", "urinary": "N20–N23 urinary", "sum3": "sum3 (combined)",
    }
    for ax, rows, title, panel in zip(
        axes,
        [br_outcomes, mx_outcomes],
        ["Brazil | full primary calendar", "Mexico | preferred 2015–2022 cohort"],
        ["b", "c"],
    ):
        by_outcome = {r["outcome"]: r for r in rows}
        missing = {"K80", "K81", "urinary", "sum3"} - set(by_outcome)
        if missing:
            raise ValueError(f"Outcome total source rows missing in {title}: {sorted(missing)}")
        outcomes = ["K80", "K81", "urinary", "sum3"]
        vals = [int(by_outcome[o]["total_events"]) for o in outcomes]
        ys = np.arange(len(outcomes))
        bars = ax.barh(ys, vals, color=[outcome_colors[o] for o in outcomes], height=0.58)
        bars.patches[1].set_alpha(0.58)
        bars.patches[1].set_hatch("//")
        ax.set_yticks(ys)
        ax.set_yticklabels([outcome_labels[o] for o in outcomes])
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=8.6, fontweight="bold", pad=8)
        ax.set_xlabel("Admission records")
        ax.grid(axis="x", color=GRID, linewidth=0.7)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        _label(ax, panel)
        upper = max(vals) * 1.27 if vals else 1
        ax.set_xlim(0, upper)
        for bar, value in zip(bars, vals):
            ax.text(value + upper * 0.018, bar.get_y() + bar.get_height() / 2,
                    f"{value:,}", ha="left", va="center", fontsize=7.2)
        ax.text(0.0, -0.22, "Admissions, not unique people",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.0, color=NEUTRAL)

    fig.suptitle("Coverage and diagnostic admission counts", fontsize=12, fontweight="bold")
    paths = _save_figure(fig, out / "fig1_design")
    legend = (
        "## Figure 1\n\n"
        "**Coverage and diagnostic admission counts.** a, Brazilian primary admission coverage (2006–2020), with the shorter clinical admission-date period ({}–{}; source years {}–{} including late billings), and the preferred complete-source Mexican period (2015–2022) plus the independent 2023 extension. The highlighted Mexican 2019 segment was restored in the preferred cohort. The Brazilian and Mexican primary periods include k={} and k={} geographic units, respectively. b,c, K80, K81, N20–N23 and combined sum3 admission records for each primary period, each country shown on its own count scale. K81 does not establish calculous disease; the broad urinary group includes unspecified renal colic. Counts are admissions, not unique people.\n\n"
        "The accompanying source data retain the registered row-level provenance for the calendar descriptors and diagnosis totals."
    ).format(clinical_period[0], clinical_period[1], clinical_source_period[0],
             clinical_source_period[1], count(cal_br["unit_count"]), count(cal_mx["unit_count"]))
    legend += ("\n\nMexico's full calendar contains " + count(cal_mx["row_count"]) +
               " state-days, including " + count(cal_mx["all3_zero_admission_days_among_health"]) +
               " zero-admission days; " + count(cal_mx["eligible_rows"]) +
               " state-days have complete exposure windows. I² and directional concordance are not applicable to these descriptive counts.")
    return paths, legend


def _figure2(pub: FrozenPublication, out: Path, supplement: Path) -> tuple[list[Path], str]:
    artifact = "figure2"
    rows = _rows(pub, "main_effects", artifact, country="BR", contrast="p95_vs_mrt", method="DL")
    endpoints = ["y_K80", "y_K81", "y_urinary", "y_sum3"]
    by_endpoint = {r["endpoint"]: r for r in rows}
    if set(by_endpoint) != set(endpoints):
        raise ValueError(f"Expected four Brazilian endpoints in main_effects; found {sorted(by_endpoint)}")
    ordered = [by_endpoint[e] for e in endpoints]
    _save_source_rows(ordered, supplement / "figure2_brazil_main_effects_source.csv")

    colors = {"y_K80": BLUE, "y_K81": BLUE, "y_urinary": RED, "y_sum3": NEUTRAL}
    markers = {"y_K80": "o", "y_K81": "s", "y_urinary": "o", "y_sum3": "D"}
    labels = {"y_K80": "K80", "y_K81": "K81", "y_urinary": "N20–N23 urinary", "y_sum3": "Combined sum3"}
    lows = [float(r["lo"]) for r in ordered]
    highs = [float(r["hi"]) for r in ordered]
    left = min(lows)
    right = max(highs)
    pad = max((right - left) * 0.08, 0.012)

    fig, (ax, stats) = plt.subplots(1, 2, figsize=(7.2, 4.6),
                                   gridspec_kw={"width_ratios": [2.7, 2.0]},
                                   constrained_layout=True)
    y = np.arange(len(ordered))[::-1]
    for yi, row in zip(y, ordered):
        est, lo, hi = float(row["RR"]), float(row["lo"]), float(row["hi"])
        ax.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt=markers[row["endpoint"]],
                    color=colors[row["endpoint"]], ecolor=colors[row["endpoint"]],
                    markersize=5.3, elinewidth=1.25, capsize=2.5, zorder=3)
    ax.axvline(1, color=NEUTRAL, linewidth=1.0, linestyle="--")
    ax.set_xlim(min(left - pad, 0.99), max(right + pad, 1.01))
    xticks = list(ax.get_xticks())
    if 1.0 not in xticks:
        xticks.append(1.0)
        ax.set_xticks(sorted(xticks))
    ax.set_ylim(-0.7, len(ordered) - 0.3)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[r["endpoint"]] for r in ordered])
    ax.set_xlabel("Rate ratio | exact P95 versus constrained MRT")
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_title("Brazil | pooled unit-specific main effects", loc="left", fontsize=9,
                 fontweight="bold")
    _label(ax, "a")
    stats.axis("off")
    stats.set_xlim(0, 1)
    stats.set_ylim(ax.get_ylim())
    column_x = {"rr": 0.0, "k": 0.55, "i2": 0.66, "direction": 0.79}
    heading_y = len(ordered) - 0.45
    stats.text(column_x["rr"], heading_y, "RR (95% CI)", fontsize=7.0,
               fontweight="bold", va="bottom")
    stats.text(column_x["k"], heading_y, "k", fontsize=7.0,
               fontweight="bold", va="bottom")
    stats.text(column_x["i2"], heading_y, "I² (%)", fontsize=7.0,
               fontweight="bold", va="bottom")
    stats.text(column_x["direction"], heading_y, "Direction >0", fontsize=7.0,
               fontweight="bold", va="bottom")
    for yi, row in zip(y, ordered):
        stats.text(column_x["rr"], yi,
                   "{} ({}–{})".format(num(row["RR"], 3), num(row["lo"], 3), num(row["hi"], 3)),
                   fontsize=7.0, va="center")
        stats.text(column_x["k"], yi, count(row["k"]), fontsize=7.0, va="center")
        stats.text(column_x["i2"], yi, num(row["I2"], 1), fontsize=7.0, va="center")
        stats.text(column_x["direction"], yi,
                   percent(row["concordance"], 1, fraction=True), fontsize=7.0, va="center")
    stats.text(0.0, -0.55, "DL pooling; 95% CI", fontsize=7.0, color=NEUTRAL, va="bottom")
    fig.suptitle("Brazilian P95-versus-MRT associations", fontsize=12, fontweight="bold")
    paths = _save_figure(fig, out / "fig2_forest")
    legend = (
        "## Figure 2\n\n"
        "**Brazilian P95-versus-MRT associations.**\n"
        "a, Pooled relative rates for K80, K81, urinary N20–N23 and the combined sum3 endpoint. Points and horizontal bars show the saved DerSimonian–Laird estimates and 95% confidence intervals. The vertical line marks a rate ratio of 1. The companion column gives the fitted geographic-unit count k, I² and the share of unit estimates above zero. K81 does not establish calculous disease; the broad urinary group includes unspecified renal colic. Pooled marginal estimates are distinct endpoints and are not used to reconstruct the direct urinary-versus-K80 contrast.\n"
        "The accompanying source data retain row-level provenance for every displayed estimate."
    )
    legend += "\n\nEndpoint-specific statistics: " + "; ".join(
        labels[r["endpoint"]] + ": k=" + count(r["k"]) + ", I²=" + percent(r["I2"], 1) +
        ", direction >0=" + percent(r["concordance"], 1, fraction=True) for r in ordered) + "."
    return paths, legend


def _core_display_label(scenario: str) -> str:
    parts = scenario.split("_")
    if len(parts) < 5:
        return scenario
    exposure = "mean" if parts[2] == "E3" else "max 3 d"
    calendar = "shared" if parts[3] == "shared" else "year-month"
    period = "full" if parts[4] == "full" else "sub"
    return f"{exposure} | {calendar} | {period}"


def _figure3(pub: FrozenPublication, out: Path, supplement: Path) -> tuple[list[Path], str]:
    artifact = "figure3"
    keywhere = {"contrast": "y_urinary_vs_y_K80", "covariance": "hac14", "method": "DL"}
    all_core = pub.rows("core_robustness", country="BR", **keywhere)
    grid = [r for r in all_core if r["scenario"].startswith("sensitivity_original8_")]
    if len(grid) != 8:
        raise ValueError(f"Expected the repaired original eight-cell grid; found {len(grid)} rows")
    grid.sort(key=lambda r: r["scenario"])
    _record_claims(pub, "core_robustness", grid,
                   {"country": "BR", "scenario_prefix": "sensitivity_original8_", **keywhere},
                   artifact)

    primary_br = _one(pub, "core_primary_contrast", artifact, country="BR",
                      scenario="primary_threeoutcome_E3_separate", **keywhere)
    primary_mx = _one(pub, "mx_updated_primary", artifact, country="MX",
                      scenario="primary_complete_source_2015_2022", **keywhere)

    time_rows = []
    time_specs = [
        ("core_robustness", "sensitivity_temporal_20060101_20161231_fixed_full_ref", "BR",
         "Brazil 2006–2016"),
        ("core_robustness", "sensitivity_temporal_20170101_20201231_fixed_full_ref", "BR",
         "Brazil 2017–2020"),
        ("mx_updated_sensitivity", "temporal_2015_2018_fixed_full_reference", "MX",
         "Mexico 2015–2018"),
        ("mx_updated_sensitivity", "temporal_2020_2022_fixed_full_reference", "MX",
         "Mexico 2020–2022"),
        ("mx_updated_sensitivity", "temporal_replication_2023_fixed_2015_2022_reference", "MX",
         "Mexico 2023"),
    ]
    for key, scenario, country, label in time_specs:
        row = _one(pub, key, artifact, country=country, scenario=scenario, **keywhere)
        time_rows.append((row, label))
    flat_rows = grid + [primary_br, primary_mx] + [r for r, _ in time_rows]
    _save_source_rows(flat_rows, supplement / "figure3_RRR_source_rows.csv")

    grid_labels = [_core_display_label(r["scenario"]) for r in grid]
    preferred = [
        (primary_br, "Brazil"),
        (primary_mx, "Mexico"),
    ]
    all_ci_lows = [float(r["lo"]) for r in flat_rows]
    all_ci_highs = [float(r["hi"]) for r in flat_rows]
    xmin, xmax = min(all_ci_lows), max(all_ci_highs)
    padding = max((xmax - xmin) * 0.06, 0.025)
    xmin, xmax = min(xmin - padding, 0.99), max(xmax + padding, 1.01)

    fig = plt.figure(figsize=(10.8, 5.7), constrained_layout=False)
    gs = fig.add_gridspec(1, 3, width_ratios=(1.34, 1.0, 1.10),
                          left=0.08, right=0.99, top=0.80, bottom=0.16, wspace=0.55)
    spec_data = [(grid, grid_labels), ([r for r, _ in preferred], [l for _, l in preferred]),
                 ([r for r, _ in time_rows], [l for _, l in time_rows])]
    panel_titles = ["a  Brazil: eight-cell grid", "b  Country estimates", "c  Period refits"]
    label_ratios = [(0.46, 0.54), (0.42, 0.58), (0.43, 0.57)]
    axes = []
    outer_titles = []
    for ix, ((rows, labels), title, label_ratio) in enumerate(
            zip(spec_data, panel_titles, label_ratios)):
        panel_gs = gs[0, ix].subgridspec(1, 2, width_ratios=label_ratio, wspace=0.02)
        label_ax = fig.add_subplot(panel_gs[0, 0])
        ax = fig.add_subplot(panel_gs[0, 1], sharex=axes[0] if axes else None)
        axes.append(ax)
        outer_titles.append((gs[0, ix], title))
        yy = np.arange(len(rows))[::-1]
        ylimits = (-0.65, len(rows) - 0.35)
        label_ax.set_xlim(0, 1)
        label_ax.set_ylim(ylimits)
        label_ax.axis("off")
        for yv, label in zip(yy, labels):
            label_ax.text(0.99, yv, label, ha="right", va="center", fontsize=10.5)
        ax.set_yticks(yy)
        ax.tick_params(axis="y", left=False, labelleft=False)
        for yv, row in zip(yy, rows):
            color = BLUE if row["country"] == "BR" else RED
            est, low, high = float(row["RR"]), float(row["lo"]), float(row["hi"])
            ax.errorbar(est, yv, xerr=[[est - low], [high - est]], fmt="o",
                        color=color, ecolor=color, markersize=5.2, elinewidth=1.25,
                        capsize=2.5, zorder=3)
        ax.axvline(1, color=NEUTRAL, linestyle="--", linewidth=1.0)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ylimits)
        ax.grid(axis="x", color=GRID, linewidth=0.75)
        ax.tick_params(axis="x", labelsize=10.5)
        ax.spines["left"].set_visible(False)

    for panel_spec, title in outer_titles:
        box = panel_spec.get_position(fig)
        fig.text(box.x0, box.y1 + 0.012, title, ha="left", va="bottom",
                 fontsize=11.0, fontweight="bold")
    for ax in axes:
        ax.set_xticks([1.00, 1.10, 1.20])
    fig.legend(handles=[
        Line2D([], [], color=BLUE, marker="o", linestyle="none", label="Brazil"),
        Line2D([], [], color=RED, marker="o", linestyle="none", label="Mexico"),
    ], loc="upper center", bbox_to_anchor=(0.5, 0.91), ncol=2, frameon=False,
        fontsize=10.5, handletextpad=0.4, columnspacing=1.4)
    fig.text(0.5, 0.055, "Direct ratio of rate ratios (RRR) | 95% CI",
             ha="center", va="center", fontsize=11.0)
    fig.suptitle("Direct diagnosis contrasts across specifications and periods",
                 fontsize=14.0, fontweight="bold", y=0.985)
    paths = _save_figure(fig, out / "fig3_contrast")
    stat_rows = list(zip(grid_labels, grid))
    stat_rows.extend((label, row) for row, label in preferred)
    stat_rows.extend((label, row) for row, label in time_rows)
    exact_stats = "; ".join(
        "{}: RRR {} ({}–{}), k={}, I²={}%, direction >0={}%".format(
            label, num(row["RR"], 3), num(row["lo"], 3), num(row["hi"], 3),
            count(row["k"]), num(row["I2"], 1),
            num(float(row["concordance"]) * 100, 1))
        for label, row in stat_rows
    )
    legend = (
        "## Figure 3\n\n"
        "**Direct diagnosis contrasts across specifications and periods.** a, The eight archived Brazil combinations selected before the repaired results were generated: mean (E3) or maximum 3-day exposure (Emax3), shared-calendar or year-month fixed effects, and full (2006–2020) or subperiod (2008–2015) analysis. These archived settings are not a claim of prospective registration. b, Preferred Brazil 2006–2020 and complete-source Mexico 2015–2022 estimates. c, Independent period refits: Brazil 2006–2016 and 2017–2020; Mexico 2015–2018 and 2020–2022; and Mexico 2023 with the 2015–2022 reference targets. Points and bars show DerSimonian–Laird estimates and 95% confidence intervals; all three panels share one horizontal RRR scale that contains every interval, including the 2023 interval, which spans 1. Brazil and Mexico are blue and red, respectively. k counts geographic units. Every estimate is a direct ratio of rate ratios, not a risk difference. The primary contrast uses its joint covariance and is not reconstructed by dividing separately pooled effects.\n\n"
        "Frozen estimates and interval/statistic values for all 15 plotted points: {}.\n\n"
        "The accompanying source data retain the exact source rows and selection provenance for all plotted estimates."
    ).format(exact_stats)
    return paths, legend


def _figure4(pub: FrozenPublication, out: Path, supplement: Path) -> tuple[list[Path], str]:
    artifact = "figure4"
    max_temp = float(pub.metadata.get("presentation_constants", {}).get("curve_display_upper_temperature_C", 30))
    if max_temp != 30:
        raise ValueError("The repaired main curve display contract requires a 30 C upper bound")
    selection = {"country": "BR", "method": "DL", "curve_scope": "support_restricted"}
    rows = pub.rows("curves_supported", **selection)
    rows = [r for r in rows if float(r["temperature"]) <= max_temp]
    _record_claims(pub, "curves_supported", rows,
                   {**selection, "temperature_lte_C": max_temp}, artifact)
    specs = [
        ("y_K80", "RR_vs_shared_median:support_restricted", "K80", BLUE),
        ("y_urinary", "RR_vs_shared_median:support_restricted", "Urinary N20–N23", RED),
        ("urinary_vs_K80", "direct_RRR:support_restricted", "Direct urinary / K80 RRR", NEUTRAL),
    ]
    curves = []
    for endpoint, contrast, label, color in specs:
        current = [r for r in rows if r["endpoint"] == endpoint and r["contrast"] == contrast]
        current.sort(key=lambda r: float(r["temperature"]))
        temps = [float(r["temperature"]) for r in current]
        if len(current) < 2 or len(temps) != len(set(temps)):
            raise ValueError(f"Invalid supported curve rows for {endpoint}/{contrast}")
        if max(temps) > max_temp:
            raise ValueError(f"Curve source exceeds the contracted {max_temp} C display bound")
        curves.append((current, label, color))
    _save_source_rows([r for curve, _, _ in curves for r in curve],
                      supplement / "figure4_supported_curve_vertices_source.csv")
    reference_temps = [float(r["temperature"]) for r in curves[0][0]]
    reference_k = [str(r["k"]) for r in curves[0][0]]
    for curve, label, _ in curves[1:]:
        curve_temps = [float(r["temperature"]) for r in curve]
        curve_k = [str(r["k"]) for r in curve]
        if curve_temps != reference_temps or curve_k != reference_k:
            raise ValueError(
                f"The shared contributor track requires identical temperature/k arrays; {label} differs"
            )
    support32 = _one(pub, "support_32c", artifact, country="BR")
    _save_source_rows([support32], supplement / "figure4_32C_support_source.csv")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=False)
    gs = fig.add_gridspec(2, 3, height_ratios=(4.8, 1.25), left=0.07, right=0.98,
                          top=0.91, bottom=0.12, wspace=0.30, hspace=0.34)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    support_ax = fig.add_subplot(gs[1, :])
    at_max_summaries = []
    for ax, (curve, label, color), panel in zip(axes, curves, ["a", "b", "c"]):
        temp = np.asarray([float(r["temperature"]) for r in curve])
        estimate = np.asarray([float(r["RR"]) for r in curve])
        low = np.asarray([float(r["lo"]) for r in curve])
        high = np.asarray([float(r["hi"]) for r in curve])
        ax.fill_between(temp, low, high, color=color, alpha=0.18, linewidth=0)
        ax.plot(temp, estimate, color=color, linewidth=1.8)
        ax.axhline(1, color=NEUTRAL, linestyle="--", linewidth=0.8)
        ax.set_xlim(min(temp), max_temp)
        ax.set_ylabel("RR" if panel != "c" else "Direct RRR")
        ax.set_title(label, loc="left", fontsize=8.4, fontweight="bold")
        ax.grid(color=GRID, linewidth=0.55)
        _label(ax, panel)
        maxrow = next(r for r in curve if float(r["temperature"]) == max_temp)
        at_max_summaries.append(
            "{}: k={}, I²={}%, direction >0={}%".format(
                label, count(maxrow["k"]), num(maxrow["I2"], 1),
                num(float(maxrow["concordance"]) * 100, 1)))
        ax.text(0.02, 0.02,
                "30 °C: k={} | I²={}%\ndirection >0={}%".format(
                    count(maxrow["k"]), num(maxrow["I2"], 1),
                    num(float(maxrow["concordance"]) * 100, 1)),
                transform=ax.transAxes, fontsize=7.0, color=NEUTRAL,
                ha="left", va="bottom", linespacing=1.15)

    # The contributor counts are identical across the three frozen curve series.
    common_temp = np.asarray(reference_temps)
    common_k = np.asarray([int(value) for value in reference_k])
    support_ax.step(common_temp, common_k, where="mid", color=NEUTRAL, linewidth=1.45)
    support_ax.set_xlim(min(common_temp), max_temp)
    support_ax.set_ylim(bottom=0)
    support_ax.set_xlabel("Three-day mean temperature (°C)")
    support_ax.set_ylabel("Contributing units, k")
    support_ax.grid(axis="y", color=GRID, linewidth=0.55)
    _label(support_ax, "d")
    fig.suptitle("Locally supported Brazilian temperature-response curves",
                 fontsize=12, fontweight="bold", y=0.975)
    paths = _save_figure(fig, out / "fig4_curves")
    legend = (
        "## Figure 4\n\n"
        "**Locally supported Brazilian temperature-response curves.** a,b, Pooled K80 and urinary N20–N23 relative rates versus each unit's shared median. c, Directly estimated joint urinary-versus-K80 ratio of rate ratios, with its own covariance. Lines and ribbons give DerSimonian–Laird estimates and pointwise 95% confidence intervals. Only regions with the displayed temperature inside the unit-specific E3 P2–P98 range and at least 30 eligible days within ±1 °C contribute; the lower track shows the shared changing k for panels a–c. This is one common contributor count, since all three frozen k arrays match at every displayed temperature. The display ends at 30 °C. Movement along a pooled curve combines changes in contributing geography and meta-analytic weights; it is not a fixed-national-population response. At 30 °C, the frozen summaries are {}. Other temperatures have their own support statistics in the source data. At 32 °C, {} of {} regions and {} region-days occur within ±1 °C, but only {} region meets the stricter curve-support flag. No curve crossing implies a biological threshold. All-successful-unit extrapolation-inclusive curve vertices are provided in the source-data supplement.\n\n".format(
            "; ".join(at_max_summaries),
            count(support32["regions_with_at_least_1_E3_day_within_1C"]),
            count(support32["regions_assessed"]),
            count(support32["total_region_days_within_1C_of_32C"]),
            count(support32["regions_marked_strictly_supported_at_32C"]))
        + "The accompanying source data retain row-level provenance, curve support counts and 32 °C support diagnostics."
    )
    return paths, legend


def render(frozen: Path, output_root: Path) -> list[Path]:
    _style()
    pub = FrozenPublication(frozen)
    required = {
        "pub_calendar", "pub_outcomes", "pub_exposure", "core_primary_contrast",
        "core_robustness", "mx_updated_primary", "mx_updated_sensitivity",
        "main_effects", "curves_supported", "curves_all_units", "support_32c",
    }
    missing = sorted(required - set(pub.tables))
    if missing:
        raise KeyError(f"Final freeze is missing required figure table keys: {missing}")
    outputs_dir = output_root / "outputs"
    source_data_dir = outputs_dir / "source_data"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    source_data_dir.mkdir(parents=True, exist_ok=True)
    legends = []
    paths = []
    for maker in (_figure1, _figure2, _figure3, _figure4):
        outputs, legend = maker(pub, outputs_dir, source_data_dir)
        paths.extend(outputs)
        legends.append(legend)
    legend_path = outputs_dir / "figure_legends.md"
    legend_path.write_text("\n\n".join(legends) + "\n", encoding="utf-8")
    manifest_path = source_data_dir / "figure_source_rows.json"
    pub.save_claims(manifest_path)
    paths.extend([legend_path, manifest_path])
    return paths


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
    print("Rendered {} figure/legend/source artifacts under {}".format(len(outputs), args.outdir))
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
