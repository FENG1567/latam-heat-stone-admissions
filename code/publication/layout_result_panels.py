"""Readable layouts for frozen Fig2, FigS1 and Fig4 results.

The parent renderer supplies publication helpers through shared so this module
can be imported without depending on a mutable rendering driver.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def _get(shared: Any, name: str, default: Any = None) -> Any:
    if isinstance(shared, Mapping):
        return shared.get(name, default)
    return getattr(shared, name, default)


def _font(shared: Any) -> float:
    return float(_get(shared, "BASE_FONT", 10.5))


def _bold_ticks(ax, size: float) -> None:
    ax.tick_params(axis="both", labelsize=size)
    for axis in (ax.xaxis, ax.yaxis):
        for tick in axis.get_major_ticks():
            for label in (tick.label1, tick.label2):
                label.set_fontfamily("Arial")
                label.set_fontweight("bold")
                label.set_fontsize(size)


def _text_kwargs(size: float) -> dict[str, Any]:
    return {"fontsize": size, "fontfamily": "Arial", "fontweight": "bold"}


def forest4(pub, editor, rows, labels, title, stem, *, shared):
    """Render four frozen estimates with separate labels and statistic columns."""
    if len(rows) != 4 or len(labels) != 4:
        raise ValueError(f"Expected four rows and labels for {stem}")
    font = _font(shared)
    BLUE, RED, GRAY, GRID = (_get(shared, key) for key in ("BLUE", "RED", "GRAY", "GRID"))
    num = _get(shared, "num")
    count = _get(shared, "count")
    percent = _get(shared, "percent")
    save_all = _get(shared, "save_all")
    if any(value is None for value in (num, count, percent, save_all)):
        raise ValueError("shared must provide num, count, percent and save_all")

    yy = np.arange(4, dtype=float)[::-1]
    lows = np.asarray([float(row["lo"]) for row in rows])
    highs = np.asarray([float(row["hi"]) for row in rows])
    xmin = float(min(lows.min(), 1.0))
    xmax = float(max(highs.max(), 1.0))
    pad = max((xmax - xmin) * 0.055, 0.015)
    xmin -= pad
    xmax += pad

    fig = plt.figure(figsize=(8.7, 3.3))
    grid = fig.add_gridspec(
        1, 6, width_ratios=(1.86, 3.05, 2.03, 0.56, 0.78, 1.27),
        left=0.035, right=0.99, top=0.835, bottom=0.205, wspace=0.045,
    )
    label_ax = fig.add_subplot(grid[0, 0])
    plot_ax = fig.add_subplot(grid[0, 1])
    rr_ax = fig.add_subplot(grid[0, 2], sharey=plot_ax)
    k_ax = fig.add_subplot(grid[0, 3], sharey=plot_ax)
    i2_ax = fig.add_subplot(grid[0, 4], sharey=plot_ax)
    positive_ax = fig.add_subplot(grid[0, 5], sharey=plot_ax)
    ylimits = (-0.64, 3.82)

    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(ylimits)
    label_ax.axis("off")
    plot_ax.set_xlim(xmin, xmax)
    plot_ax.set_ylim(ylimits)
    plot_ax.set_yticks(yy)
    plot_ax.tick_params(axis="y", left=False, labelleft=False)
    plot_ax.axvline(1.0, color=GRAY, linestyle="--", linewidth=1.0, zorder=0)
    plot_ax.grid(axis="x", color=GRID, linewidth=0.7)
    plot_ax.set_axisbelow(True)
    plot_ax.set_xlabel("Relative rate (RR)", **_text_kwargs(font), labelpad=5)
    _bold_ticks(plot_ax, font)

    marker = {"y_K80": "s", "y_K81": "D", "y_urinary": "o", "y_sum3": "P"}
    for y, row, label in zip(yy, rows, labels):
        color = BLUE if row.get("country") == "BR" else RED
        rr, lo, hi = (float(row[key]) for key in ("RR", "lo", "hi"))
        plot_ax.errorbar(
            rr, y, xerr=[[rr - lo], [hi - rr]],
            fmt=marker.get(row.get("endpoint"), "o"), color=color, ecolor=color,
            markersize=5.4, elinewidth=1.45, capsize=2.7, zorder=3,
        )
        visible_label = str(label).replace(" | ", "\n")
        label_ax.text(0.035, y, visible_label, ha="left", va="center",
                      clip_on=True, **_text_kwargs(font))

    column_specs = [
        (rr_ax, "RR (95% CI)", 0.025, "left"),
        (k_ax, "k", 0.5, "center"),
        (i2_ax, "I² (%)", 0.5, "center"),
        (positive_ax, "Positive\nunits (%)", 0.5, "center"),
    ]
    for ax, heading, x, align in column_specs:
        ax.set_xlim(0, 1)
        ax.set_ylim(ylimits)
        ax.axis("off")
        ax.text(x, 3.53, heading, ha=align, va="bottom", **_text_kwargs(font))
    for y, row in zip(yy, rows):
        values = [
            (rr_ax, f"{num(row['RR'])} ({num(row['lo'])}–{num(row['hi'])})", 0.025, "left"),
            (k_ax, count(row["k"]), 0.5, "center"),
            (i2_ax, num(row["I2"], 1), 0.5, "center"),
            (positive_ax, percent(row["concordance"], 1, fraction=True), 0.5, "center"),
        ]
        for ax, value, x, align in values:
            ax.text(x, y, value, ha=align, va="center", clip_on=True,
                    **_text_kwargs(font))

    fig.text(0.035, 0.966, title, ha="left", va="top", **_text_kwargs(12.0))
    outstem = Path(editor) / stem
    save_all(fig, outstem)
    return outstem


def fig4(pub, editor, *, shared):
    """Render the original six frozen curve series and shared support counts."""
    font = _font(shared)
    BLUE, RED, GRAY, GRID = (_get(shared, key) for key in ("BLUE", "RED", "GRAY", "GRID"))
    role_table = _get(shared, "role_table")
    record = _get(shared, "record")
    save_all = _get(shared, "save_all")
    if any(value is None for value in (role_table, record, save_all)):
        raise ValueError("shared must provide role_table, record and save_all")

    table = role_table(pub, "curves")
    selection = {"country": "BR", "curve_scope": "support_restricted"}
    rows = pub.rows(table, **selection)
    rows = [row for row in rows if 14 <= float(row["temperature"]) <= 30]
    if not rows:
        raise ValueError("No locally supported temperature-curve rows in 14–30 °C")
    for row in rows:
        record(
            pub, "fig4_curves", "a-c", table, row,
            {**selection, "temperature": row["temperature"],
             "endpoint": row["endpoint"], "method": row["method"]},
            ("endpoint", "temperature", "method", "RR", "lo", "hi",
             "k", "I2", "concordance"),
        )

    endpoints = [
        ("y_K80", "K80", BLUE),
        ("y_urinary", "Urinary N20–N23", RED),
        ("urinary_vs_K80", "Direct urinary/K80 RRR", GRAY),
    ]
    methods = [("DL", "-", 0.13), ("REML_profile_HK", "--", 0.055)]
    curves = {}
    for endpoint, label, color in endpoints:
        for method, linestyle, alpha in methods:
            current = sorted(
                [row for row in rows if row["endpoint"] == endpoint and row["method"] == method],
                key=lambda row: float(row["temperature"]),
            )
            if len(current) < 2:
                raise ValueError(f"Missing frozen curve vertices: {endpoint}/{method}")
            curves[(endpoint, method)] = current

    ref = curves[("y_K80", "DL")]
    ref_temp = [float(row["temperature"]) for row in ref]
    ref_k = [str(row["k"]) for row in ref]
    for current in curves.values():
        if [float(row["temperature"]) for row in current] != ref_temp:
            raise ValueError("Curve series do not share the frozen 14–30 °C vertices")
        if [str(row["k"]) for row in current] != ref_k:
            raise ValueError("Curve series do not share common local-support k")

    fig = plt.figure(figsize=(8.8, 4.9))
    fig.text(0.10, 0.97, "A", ha="left", va="top", **_text_kwargs(18.0))
    fig.text(0.145, 0.97, "K80 and urinary relative rates",
             ha="left", va="top", **_text_kwargs(12.0))
    fig.text(0.64, 0.97, "B", ha="left", va="top", **_text_kwargs(18.0))
    fig.text(0.685, 0.97, "Direct urinary/K80 contrast",
             ha="left", va="top", **_text_kwargs(12.0))

    handles = [
        Line2D([], [], color=BLUE, lw=1.8, label="K80"),
        Line2D([], [], color=RED, lw=1.8, label="Urinary N20–N23"),
        Line2D([], [], color=GRAY, lw=1.8, label="DL (solid)"),
        Line2D([], [], color=GRAY, lw=1.8, linestyle="--", label="REML–HK (dashed)"),
    ]
    legend = fig.legend(
        handles=handles, loc="center", bbox_to_anchor=(0.5, 0.855),
        ncol=4, frameon=False, columnspacing=1.4, handlelength=1.8,
        handletextpad=0.5,
        prop={"family": "Arial", "weight": "bold", "size": font},
    )
    for text in legend.get_texts():
        text.set_fontfamily("Arial")
        text.set_fontweight("bold")
        text.set_fontsize(font)

    ax_a = fig.add_axes([0.10, 0.40, 0.44, 0.39])
    ax_b = fig.add_axes([0.64, 0.40, 0.345, 0.39], sharex=ax_a)
    for endpoint, label, color in endpoints[:2]:
        for method, linestyle, alpha in methods:
            current = curves[(endpoint, method)]
            x = np.asarray([float(row["temperature"]) for row in current])
            point = np.asarray([float(row["RR"]) for row in current])
            low = np.asarray([float(row["lo"]) for row in current])
            high = np.asarray([float(row["hi"]) for row in current])
            ax_a.fill_between(x, low, high, color=color, alpha=alpha, linewidth=0)
            ax_a.plot(x, point, color=color, linestyle=linestyle, linewidth=1.8,
                      zorder=3 if method == "DL" else 2)

    endpoint = "urinary_vs_K80"
    for method, linestyle, alpha in methods:
        current = curves[(endpoint, method)]
        x = np.asarray([float(row["temperature"]) for row in current])
        point = np.asarray([float(row["RR"]) for row in current])
        low = np.asarray([float(row["lo"]) for row in current])
        high = np.asarray([float(row["hi"]) for row in current])
        ax_b.fill_between(x, low, high, color=GRAY, alpha=alpha, linewidth=0)
        ax_b.plot(x, point, color=GRAY, linestyle=linestyle, linewidth=1.9)

    for ax, ylabel in ((ax_a, "RR"), (ax_b, "RRR")):
        ax.axhline(1.0, color=GRAY, linestyle=":", linewidth=1.0)
        ax.set_xlim(14, 30)
        ax.set_ylabel(ylabel, **_text_kwargs(font), labelpad=5)
        ax.grid(color=GRID, linewidth=0.7)
        ax.set_xticks([14, 18, 22, 26, 30])
        ax.set_xlabel("")
        _bold_ticks(ax, font)
        ax.tick_params(labelbottom=True)

    fig.text(0.10, 0.32, "C", ha="left", va="top", **_text_kwargs(18.0))
    fig.text(0.145, 0.32, "Contributing regions",
             ha="left", va="top", **_text_kwargs(12.0))

    ax_c = fig.add_axes([0.10, 0.15, 0.885, 0.10], sharex=ax_a)
    temp = np.asarray(ref_temp)
    kvals = np.asarray([int(row["k"]) for row in ref])
    ax_c.step(temp, kvals, where="mid", color=GRAY, linewidth=1.8)
    ax_c.set_xlim(14, 30)
    ax_c.set_ylim(bottom=0)
    ax_c.set_xticks([14, 18, 22, 26, 30])
    ax_c.set_xlabel("Three-day mean temperature (°C)",
                    **_text_kwargs(font), labelpad=2)
    ax_c.set_ylabel("Regions, k", **_text_kwargs(font), labelpad=2)
    ax_c.grid(axis="y", color=GRID, linewidth=0.7)
    _bold_ticks(ax_c, font)

    outstem = Path(editor) / "fig4_curves"
    save_all(fig, outstem)
    return outstem
