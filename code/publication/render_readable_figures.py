#!/usr/bin/env python3
"""Render compact, source-traceable publication figures from FrozenPublication only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import render_repaired_figures as base
from publication_common import FrozenPublication, count, num, percent

BLUE, RED, GRAY, GRID = base.BLUE, base.RED, base.NEUTRAL, base.GRID
WIDTH_IN = 8.4
BASE_FONT = 10.5
AUDIT: list[dict] = []


def style():
    import os
    from matplotlib import font_manager
    base._style()
    font_path = os.environ.get("LATAM_ARIAL_BOLD")
    if font_path:
        font_manager.fontManager.addfont(font_path)
        family = font_manager.FontProperties(fname=font_path).get_name()
    else:
        try:
            font_manager.findfont(font_manager.FontProperties(family="Arial", weight="bold"), fallback_to_default=False)
            family = "Arial"
        except ValueError:
            family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family, "font.sans-serif": [family, "Arial", "DejaVu Sans"],
        "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
        "font.size": 10.5, "axes.labelsize": 10.5, "axes.titlesize": 12,
        "xtick.labelsize": 10, "ytick.labelsize": 10.5,
        "lines.linewidth": 1.5, "axes.linewidth": 0.9,
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
        "text.color": "#202020", "axes.labelcolor": "#202020",
        "xtick.color": "#202020", "ytick.color": "#202020",
    })



def record(pub, figure, panel, table, row, selection, visible=()):
    pub.claims.append({
        "artifact": figure, "table": table, "selection": dict(selection),
        "source_file": row["source_file"], "source_row": row["source_row"],
    })
    item = {
        "figure": figure, "panel": panel, "table": table,
        "source_file": row["source_file"], "source_row": row["source_row"],
        "selection": dict(selection),
    }
    if visible:
        item["frozen_display_fields"] = {key: row.get(key) for key in visible}
    AUDIT.append(item)


def selected_rows(pub, figure, panel, table, **where):
    rows = pub.rows(table, **where)
    if not rows:
        raise ValueError(f"Empty frozen selection: {table} {where}")
    for row in rows:
        record(pub, figure, panel, table, row, where,
               ("endpoint", "scenario", "temperature", "method", "RR", "lo", "hi",
                "k", "I2", "concordance", "total_events", "first_date", "last_date"))
    return rows


def selected_one(pub, figure, panel, table, **where):
    rows = selected_rows(pub, figure, panel, table, **where)
    if len(rows) != 1:
        raise ValueError(f"Expected one frozen row: {table} {where}; got {len(rows)}")
    return rows[0]


def role_table(pub, role):
    return pub.metadata["revision"]["tables_by_role"][role]


def primary(pub, figure, panel, country, endpoint):
    return selected_one(
        pub, figure, panel, role_table(pub, "marginal_sync"),
        country=country, endpoint=endpoint, bandwidth_days="14",
        inference="calendar_synchronous_plus_DL_tau_component",
    )


def save_all(fig, stem: Path):
    from matplotlib.text import Text
    import os
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = []
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        axis = text.axes
        if axis is not None and not axis.axison and text not in axis.texts:
            continue
        box = text.get_window_extent(renderer)
        visible.append({"text": text.get_text(), "font": text.get_fontproperties().get_name(),
                        "weight": str(text.get_fontweight()), "size_pt": text.get_fontsize(),
                        "bbox_px": list(box.extents)})
    outside = [v for v in visible if v["bbox_px"][0] < -1 or v["bbox_px"][1] < -1 or
               v["bbox_px"][2] > fig.bbox.width+1 or v["bbox_px"][3] > fig.bbox.height+1]
    collisions = []
    for i, one in enumerate(visible):
        a = one["bbox_px"]
        for two in visible[i+1:]:
            b = two["bbox_px"]
            iw = min(a[2],b[2])-max(a[0],b[0]); ih = min(a[3],b[3])-max(a[1],b[1])
            if iw > 1.0 and ih > 1.0:
                collisions.append({"a":one["text"],"b":two["text"],"overlap_px":[iw,ih]})
    info = {"figure":stem.name,"native_size_inches":list(fig.get_size_inches()),
            "text":visible,"outside_canvas":outside,"text_collisions":collisions}
    evidence = os.environ.get("LATAM_LAYOUT_EVIDENCE")
    if evidence:
        Path(evidence).mkdir(parents=True,exist_ok=True)
        (Path(evidence)/(stem.name+"_layout.json")).write_text(json.dumps(info,ensure_ascii=False,indent=2))
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)



def validate_exports(stems):
    for stem in stems:
        png = stem.with_suffix(".png")
        pdf = stem.with_suffix(".pdf")
        svg = stem.with_suffix(".svg")
        tiff = stem.with_suffix(".tiff")
        for path in (png, pdf, svg, tiff):
            if not path.is_file() or path.stat().st_size < 1000:
                raise ValueError(f"Missing or undersized export: {path}")
        with Image.open(png) as im:
            if im.width < 1600:
                raise ValueError(f"PNG width below 1600 px: {png} ({im.width})")
        with Image.open(tiff) as im:
            dpi = im.info.get("dpi", (0, 0))
            compression = im.tag_v2.get(259)
            if min(dpi) < 599 or compression != 5:
                raise ValueError(
                    f"TIFF export contract failed for {tiff}: dpi={dpi}, compression={compression}"
                )
        svg_text = svg.read_text(encoding="utf-8")
        if "<text" not in svg_text:
            raise ValueError(f"SVG contains no editable text nodes: {svg}")
        if b"/TrueType" not in pdf.read_bytes() and b"/FontFile2" not in pdf.read_bytes():
            raise ValueError(f"PDF does not appear to embed editable TrueType fonts: {pdf}")
    print("EXPORT_QA " + json.dumps([
        {
            "stem": str(stem),
            "png_px": list(Image.open(stem.with_suffix(".png")).size),
            "tiff_dpi": list(map(float, Image.open(stem.with_suffix(".tiff")).info.get("dpi", (0, 0)))),
            "tiff_compression_tag_259": Image.open(stem.with_suffix(".tiff")).tag_v2.get(259),
            "svg_text_nodes": stem.with_suffix(".svg").read_text(encoding="utf-8").count("<text"),
        } for stem in stems
    ], ensure_ascii=False), flush=True)


def year_position(value):
    return base._year_position(value)


def inclusive_period(start, end):
    return base._inclusive_period(start, end)


def panel_heading(fig, x, y, letter, title):
    gap = 24.0 / (fig.get_figwidth() * 72.0)
    fig.text(x, y, letter, ha="left", va="top", fontsize=18, fontweight="bold")
    fig.text(x+gap, y-0.005, title, ha="left", va="top", fontsize=12, fontweight="bold")

def fig1(pub, editor):
    table = "pub_calendar"
    cal_br = selected_one(pub, "fig1_design", "a", table, dataset="BR_calendar", country="BR")
    cal_mx = selected_one(pub, "fig1_design", "a", table, dataset="MX_2015_2022", country="MX")
    cal_mx23 = selected_one(pub, "fig1_design", "a", table, dataset="MX_2023_extension", country="MX")
    all_br = pub.rows("pub_outcomes", dataset="BR_calendar", country="BR", row_scope="all_health")
    all_mx = pub.rows("pub_outcomes", dataset="MX_2015_2022", country="MX", row_scope="all_health")
    chosen = {"K80", "K81", "urinary", "sum3"}
    br = [r for r in all_br if r["outcome"] in chosen]
    mx = [r for r in all_mx if r["outcome"] in chosen]
    for dataset, country, rows in [
        ("BR_calendar", "BR", br), ("MX_2015_2022", "MX", mx)
    ]:
        if {r["outcome"] for r in rows} != chosen or len(rows) != 4:
            raise ValueError(f"Figure 1 diagnosis totals incomplete for {dataset}")
        for row in rows:
            record(pub, "fig1_design", "b" if country == "BR" else "c",
                   "pub_outcomes", row,
                   {"dataset": dataset, "country": country, "row_scope": "all_health",
                    "outcome": row["outcome"]},
                   ("outcome", "total_events", "dataset"))
    constants = pub.metadata.get("presentation_constants", {})
    clinical = constants.get("Brazil_clinical_admission_period")
    clinical_source = constants.get("Brazil_clinical_source_period_primary")
    if not isinstance(clinical, list) or len(clinical) != 2:
        raise ValueError("Missing frozen Brazilian clinical admission period")
    if not isinstance(clinical_source, list) or len(clinical_source) != 2:
        raise ValueError("Missing frozen Brazilian clinical source period")

    fig = plt.figure(figsize=(8.4,5.2))
    panel_heading(fig,0.020,0.967,"A","Study periods")
    fig.text(0.98,0.982,"Coverage and admission counts",ha="right",va="top",fontsize=10.5,fontweight="bold")
    timeline=fig.add_axes([0.205,0.57,0.765,0.330])
    timeline.set_ylim(-0.50,3.50);timeline.set_xlim(2005.5,2024.4)
    timeline.set_yticks([3,2,1,0],["Brazil primary","Brazil clinical","Mexico primary","Mexico extension"])
    timeline.set_xticks([2006,2010,2014,2018,2022]);timeline.set_xlabel("Admission year",labelpad=3)
    timeline.tick_params(axis="y",length=0,pad=7)
    timeline.grid(axis="x",color=GRID,linewidth=0.7);timeline.set_axisbelow(True)
    def period(row,y,color,label,linestyle="-",lw=6.5):
        x0=year_position(row["first_date"]);x1=year_position(row["last_date"])+1/366
        timeline.plot([x0,x1],[y,y],color=color,linewidth=lw,solid_capstyle="butt",linestyle=linestyle,zorder=2)
        timeline.text((x0+x1)/2,y+0.19,label,ha="center",va="bottom",fontsize=10.5,fontweight="bold",color=color)
    period(cal_br,3,BLUE,"2006–2020")
    period(cal_mx,1,BLUE,"2015–2022")
    period(cal_mx23,0,RED,"2023")
    timeline.fill_between([2019,2020],0.74,1.27,color=RED,alpha=0.10,zorder=0)
    timeline.text(2019.5,0.66,"2019 included",ha="center",va="top",fontsize=10,fontweight="bold",color=RED)
    c0,c1=inclusive_period(int(clinical[0]),int(clinical[1]))
    timeline.plot([c0,c1],[2,2],color=RED,linewidth=2,linestyle=(0,(2,2)),zorder=3)
    timeline.text((c0+c1)/2,2.18,"{}–{}".format(*clinical),fontsize=10.5,fontweight="bold",color=RED,ha="center",va="bottom")
    colors={"K80":BLUE,"K81":BLUE,"urinary":RED,"sum3":GRAY}
    names={"K80":"K80","K81":"K81","urinary":"Urinary\nN20–N23","sum3":"Combined"}
    for rows,country,period_label,letter,rect,hx in [
        (br,"Brazil","2006–2020","B",[0.130,0.16,0.335,0.265],0.020),
        (mx,"Mexico","2015–2022","C",[0.635,0.16,0.335,0.265],0.525),
    ]:
        panel_heading(fig,hx,0.465,letter,f"{country} | {period_label}")
        ax=fig.add_axes(rect);by={r["outcome"]:r for r in rows}
        vals=[int(by[o]["total_events"]) for o in ("K80","K81","urinary","sum3")]
        largest=max(vals)
        scale,unit=(1000000,"millions") if largest>=1000000 else (1000,"thousands") if largest>=10000 else (1,"records")
        ys=np.arange(4);bars=ax.barh(ys,np.asarray(vals)/scale,color=[colors[o] for o in ("K80","K81","urinary","sum3")],height=0.59)
        bars.patches[1].set_alpha(0.55);bars.patches[1].set_hatch("//")
        ax.set_yticks(ys,[names[o] for o in ("K80","K81","urinary","sum3")]);ax.invert_yaxis()
        upper=max(vals)/scale*1.48;ax.set_xlim(0,upper)
        ax.set_xlabel(f"Admission records ({unit})",labelpad=4);ax.tick_params(axis="y",length=0,pad=5)
        ax.grid(axis="x",color=GRID,linewidth=0.65);ax.set_axisbelow(True)
        for bar,value in zip(bars,vals):
            ax.text(value/scale+upper*0.017,bar.get_y()+bar.get_height()/2,f"{value:,}",ha="left",va="center",fontsize=10.5,fontweight="bold")
    fig.text(0.5,0.025,"Counts: administrative admission records",ha="center",va="bottom",fontsize=10.5,fontweight="bold",color=GRAY)
    stem=editor/"fig1_design";save_all(fig,stem);return stem



def forest4(pub,editor,rows,labels,title,stem):
    import layout_result_panels
    return layout_result_panels.forest4(pub,editor,rows,labels,title,stem,shared=sys.modules[__name__])



def fig2(pub, editor):
    rows = []
    labels = []
    for country, country_label in [("BR", "Brazil"), ("MX", "Mexico")]:
        for endpoint, endpoint_label in [("y_K80", "K80"),
                                         ("y_urinary", "Urinary N20–N23")]:
            rows.append(primary(pub, "fig2_forest", "a", country, endpoint))
            labels.append(f"{country_label} | {endpoint_label}")
    return forest4(pub, editor, rows, labels,
                   "Diagnosis-specific heat associations", "fig2_forest")


def fig_s1(pub, editor):
    table = "main_effects"
    rows = selected_rows(pub, "figS1_mrt_exploratory", "a", table,
                         country="BR", contrast="p95_vs_mrt", method="DL")
    expected = ["y_K80", "y_K81", "y_urinary", "y_sum3"]
    by = {row["endpoint"]: row for row in rows}
    if len(rows) != 4 or set(by) != set(expected):
        raise ValueError(f"Unexpected selected-MRT rows: {sorted(by)}")
    rows = [by[key] for key in expected]
    labels = ["K80", "K81", "Urinary N20–N23", "Combined"]
    return forest4(
        pub, editor / "supplementary_figures", rows, labels,
        "P95 versus minimum-risk temperature", "figS1_mrt_exploratory",
    )


def grid_label(scenario):
    parts = scenario.split("_")
    if len(parts) < 5 or parts[0:2] != ["sensitivity", "original8"]:
        raise ValueError(f"Unexpected original-grid scenario: {scenario}")
    exposure = parts[2]
    calendar = "shared" if parts[3] == "shared" else "year–month"
    if parts[4] == "full":
        period = "2006–2020"
    elif parts[4] == "sub":
        period = "2008–2015"
    else:
        raise ValueError(f"Unknown archived period: {scenario}")
    return f"{exposure} · {calendar}\n{period}"


def fig3(pub, editor):
    keywhere = {"contrast": "y_urinary_vs_y_K80", "covariance": "hac14", "method": "DL"}
    core_table = "core_robustness"
    all_core = pub.rows(core_table, country="BR", **keywhere)
    grid = [r for r in all_core if r["scenario"].startswith("sensitivity_original8_")]
    if len(grid) != 8:
        raise ValueError(f"Expected eight frozen grid rows, got {len(grid)}")
    grid.sort(key=lambda row: row["scenario"])
    for row in grid:
        record(pub, "fig3_contrast", "a", core_table, row,
               {"country": "BR", "scenario": row["scenario"], **keywhere},
               ("scenario", "RR", "lo", "hi", "k", "I2", "concordance"))
    br = primary(pub, "fig3_contrast", "b", "BR", "urinary_vs_K80")
    mx = primary(pub, "fig3_contrast", "b", "MX", "urinary_vs_K80")
    time_specs = [
        ("core_robustness", "sensitivity_temporal_20060101_20161231_fixed_full_ref", "BR", "BR 2006–2016"),
        ("core_robustness", "sensitivity_temporal_20170101_20201231_fixed_full_ref", "BR", "BR 2017–2020"),
        ("mx_updated_sensitivity", "temporal_2015_2018_fixed_full_reference", "MX", "MX 2015–2018"),
        ("mx_updated_sensitivity", "temporal_2020_2022_fixed_full_reference", "MX", "MX 2020–2022"),
        ("mx_updated_sensitivity", "temporal_replication_2023_fixed_2015_2022_reference", "MX", "MX 2023"),
    ]
    periods = [
        (selected_one(pub, "fig3_contrast", "c", table,
                      country=country, scenario=scenario, **keywhere), label)
        for table, scenario, country, label in time_specs
    ]
    all_rows = grid + [br, mx] + [row for row, _ in periods]
    lo_all = min(float(row["lo"]) for row in all_rows)
    hi_all = max(float(row["hi"]) for row in all_rows)
    pad = max((hi_all - lo_all) * 0.04, 0.004)
    xmin = min(0.94, lo_all - pad)
    xmax = max(1.25, hi_all + pad)
    fig=plt.figure(figsize=(9.0,5.4))
    panel_heading(fig,0.018,0.965,"A","Brazil sensitivity grid")
    panel_heading(fig,0.500,0.965,"B","Primary estimates")
    panel_heading(fig,0.500,0.560,"C","Independent period estimates")
    specifications=[
        ([0.018,0.16,0.218,0.71],[0.247,0.16,0.210,0.71],grid,[grid_label(r["scenario"]) for r in grid],list(np.arange(7,-1,-1,dtype=float)),(-0.65,7.65)),
        ([0.500,0.66,0.078,0.21],[0.588,0.66,0.215,0.21],[br,mx],["Brazil","Mexico"],[1.,0.],(-0.7,1.7)),
        ([0.500,0.16,0.166,0.335],[0.677,0.16,0.310,0.335],[r for r,_ in periods],[label.replace("BR ","Brazil ").replace("MX ","Mexico ") for _,label in periods],list(np.arange(4,-1,-1,dtype=float)),(-0.65,4.65)),
    ]
    plot_axes=[]
    for index,(label_rect,plot_rect,rows,labels,positions,ylims) in enumerate(specifications):
        la=fig.add_axes(label_rect);ax=fig.add_axes(plot_rect);plot_axes.append(ax)
        la.set_xlim(0,1);la.set_ylim(*ylims);la.axis("off")
        for yy,label,row in zip(positions,labels,rows):
            color=BLUE if row["country"]=="BR" else RED
            la.text(0.98,yy,label,ha="right",va="center",fontsize=10.5,fontweight="bold",color="#202020" if index==0 else color)
        ax.set_xlim(xmin,xmax);ax.set_ylim(*ylims);ax.set_yticks([])
        ax.axvline(1,color=GRAY,linestyle="--",linewidth=1,zorder=0)
        ax.grid(axis="x",color=GRID,linewidth=0.7);ax.set_xticks([1.0,1.1,1.2])
        for yy,row in zip(positions,rows):
            color=BLUE if row["country"]=="BR" else RED
            est,low,high=(float(row[k]) for k in ("RR","lo","hi"))
            ax.errorbar(est,yy,xerr=[[est-low],[high-est]],fmt="o",color=color,ecolor=color,markersize=5.4,elinewidth=1.4,capsize=2.5,zorder=3)
        ax.spines["left"].set_visible(False)
        if index==1:
            st=fig.add_axes([0.817,0.66,0.170,0.21]);st.set_xlim(0,1);st.set_ylim(*ylims);st.axis("off")
            st.text(0,1.52,"RRR (95% CI)",ha="left",va="bottom",fontsize=10.5,fontweight="bold")
            for yy,row in zip(positions,rows):
                st.text(0,yy,f"{num(row['RR'])} ({num(row['lo'])}–{num(row['hi'])})",ha="left",va="center",fontsize=10.5,fontweight="bold")
    fig.text(0.50,0.075,"Direct ratio of rate ratios (RRR), 95% CI",ha="center",va="center",fontsize=11,fontweight="bold")
    fig.text(0.50,0.023,"Brazil: blue    |    Mexico: red",ha="center",va="bottom",fontsize=10.5,fontweight="bold",color=GRAY)
    stem=editor/"fig3_contrast";save_all(fig,stem)
    print("FIG3_X_RANGE "+json.dumps({"all_frozen_ci_min":lo_all,"all_frozen_ci_max":hi_all,"display_xlim":[xmin,xmax],"ci_coverage_pass":xmin<=lo_all and xmax>=hi_all}),flush=True)
    return stem



def fig4(pub,editor):
    import layout_result_panels
    return layout_result_panels.fig4(pub,editor,shared=sys.modules[__name__])



def render(frozen, editor):
    editor = Path(editor)
    editor.mkdir(parents=True, exist_ok=True)
    (editor / "supplementary_figures").mkdir(parents=True, exist_ok=True)
    style()
    pub = FrozenPublication(frozen)
    stems = [
        fig1(pub, editor),
        fig2(pub, editor),
        fig3(pub, editor),
        fig4(pub, editor),
        fig_s1(pub, editor),
    ]
    validate_exports(stems)
    print("SOURCE_ROW_AUDIT " + json.dumps(AUDIT, ensure_ascii=False, separators=(",", ":")),
          flush=True)
    print("FIGURE_OUTPUTS " + json.dumps([str(stem) for stem in stems], ensure_ascii=False),
          flush=True)
    print("FROZEN_VALIDATION " + json.dumps({
        "path": str(pub.path), "checksum_verified_by": "FrozenPublication/read_frozen_payload",
        "selected_source_rows": len(AUDIT),
        "point_estimates_ci_k_I2_positive_share": "directly selected from frozen rows; no estimation",
    }, ensure_ascii=False), flush=True)
    return stems


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    stems = render(args.frozen, args.outdir)
    print(f"Rendered {len(stems)} figure groups and {len(stems) * 4} files.", flush=True)


if __name__ == "__main__":
    main()
