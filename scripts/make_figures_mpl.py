#!/usr/bin/env python3
"""Matplotlib-based manuscript figures 1-7 for Nature Communications submission.

Run:  python scripts/make_figures_mpl.py
All PDFs land in manuscript/figures/ (same paths as the ReportLab versions).
"""

from __future__ import annotations

import math
import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

ROOT     = Path(__file__).resolve().parents[1]
FINAL    = ROOT / "outputs" / "final_framework_outputs"
ANALYSIS = ROOT / "outputs" / "analysis_outputs"
FIG_DIR  = ROOT / "manuscript" / "figures"
SOURCE_DIR = ROOT / "manuscript" / "source_data"

# Palette
BLUE    = "#2B6CB0"
ORANGE  = "#C56B45"
GREEN   = "#2F855A"
GRAY    = "#6B7280"
L_GRAY  = "#E5E7EB"
M_GRAY  = "#9CA3AF"
DARK    = "#111827"
PALE    = "#F8FAFC"
RED     = "#B64B4B"

matplotlib.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size":           8,
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
    "axes.linewidth":      0.8,
    "axes.edgecolor":      M_GRAY,
    "xtick.major.width":   0.8,
    "ytick.major.width":   0.8,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "axes.grid":           False,
})

SITE_ORDER   = {"A-2018": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3}
DRUG_SHORT   = {
    "Ciprofloxacin": "Cipro",
    "Norfloxacin": "Norflox",
    "Amoxicillin-Clavulanic acid": "Amox-Clav",
    "Ceftriaxone": "CRO",
    "Ceftazidime": "CAZ",
    "Cefepime": "FEP",
    "Oxacillin": "Oxa",
}
DRUG_ORDER   = ["Cipro", "Norflox", "Amox-Clav", "CRO", "CAZ", "FEP"]


# Shared utilities

def _parse_ci(text) -> tuple[float, float, float]:
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return math.nan, math.nan, math.nan
    m = re.match(r"\s*([0-9.]+)\s+\(([0-9.]+)-([0-9.]+)\)", str(text))
    if not m:
        try:
            return float(text), math.nan, math.nan
        except ValueError:
            return math.nan, math.nan, math.nan
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def _prep_ci(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in ["raw_auc_95ci", "stratum_centered_auc_95ci"]:
        parsed = df[col].map(_parse_ci)
        df[col + "_mid"]  = parsed.map(lambda t: t[0])
        df[col + "_low"]  = parsed.map(lambda t: t[1])
        df[col + "_high"] = parsed.map(lambda t: t[2])
    return df


def _errbars(ax, y, lo, hi, color, lw=0.45):
    if math.isnan(lo) or math.isnan(hi):
        return
    ax.plot([lo, hi], [y, y], color=color, linewidth=lw, zorder=2)
    for xv in [lo, hi]:
        ax.plot([xv, xv], [y - 0.10, y + 0.10], color=color, linewidth=lw, zorder=2)


def _auc_cell_color(val: float) -> str:
    if math.isnan(val):
        return "white"
    v = max(0.45, min(0.90, float(val)))
    if v < 0.58:
        t = (v - 0.45) / 0.13
        r = int(250 + (197 - 250) * t)
        g = int(245 + (107 - 245) * t)
        b = int(241 + (69  - 241) * t)
    elif v < 0.68:
        t = (v - 0.58) / 0.10
        r = int(239 + (169 - 239) * t)
        g = int(246 + (200 - 246) * t)
        b = int(255 + (232 - 255) * t)
    else:
        t = (v - 0.68) / 0.22
        r = int(219 + (43  - 219) * t)
        g = int(235 + (108 - 235) * t)
        b = int(248 + (176 - 248) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def _phi_color(phi: float) -> str:
    if math.isnan(phi):
        return "white"
    t = max(0.0, min(1.0, float(phi)))
    r = int(247 + (43  - 247) * t)
    g = int(250 + (108 - 250) * t)
    b = int(252 + (176 - 252) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def _fig_header(fig, title: str, subtitle: str) -> None:
    """Write figure title and subtitle at the top of the figure."""
    fig.text(0.01, 0.985, title,
             fontsize=12, fontweight="bold", color=DARK,
             va="top", ha="left", transform=fig.transFigure)
    fig.text(0.01, 0.946, subtitle,
             fontsize=8, color=GRAY,
             va="top", ha="left", transform=fig.transFigure)


def _panel_label_above(fig, ax, text_str: str, subtitle_str: str = "") -> None:
    """Write a bold panel label and optional subtitle just above axes.

    Uses ax.set_title for the label (automatically fits within figure canvas)
    and ax.text in axes-transform coordinates for the subtitle (no overflow).
    Keep subtitle_str <= ~50 chars to stay within panel width.
    """
    ax.set_title(text_str, loc="left", pad=7, fontsize=9.5,
                 fontweight="bold", color=DARK)
    if subtitle_str:
        # y=1.22 sits above the ax.set_title text without overlap
        ax.text(0.0, 1.22, subtitle_str,
                transform=ax.transAxes, fontsize=7.5, color=GRAY,
                va="bottom", ha="left", clip_on=False)


def _style_data_ax(ax, xmin: float, xmax: float, n_rows: int,
                   xtick_step: float = 0.1) -> None:
    """Apply uniform axis style to a dumbbell-plot axis."""
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ticks = [round(t, 2) for t in np.arange(
        round(math.ceil(xmin * 10) / 10, 1),
        xmax + 0.001,
        xtick_step,
    ) if t <= xmax + 0.001]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.1f}" for t in ticks], fontsize=7)
    for t in ticks:
        ax.axvline(t, color=L_GRAY, linewidth=0.7, zorder=0)
    ax.axvline(0.5, color=L_GRAY, linewidth=0.7, linestyle="--", zorder=1)
    ax.set_xlabel("AUC", fontsize=8, color=GRAY, labelpad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(M_GRAY)
    ax.tick_params(axis="x", length=3, color=M_GRAY, labelcolor=GRAY)
    ax.tick_params(axis="y", length=0, pad=4)


def _style_forest_ax(ax, xmin: float, xmax: float,
                     y_positions: list[float], ylabels: list[str],
                     xtick_step: float = 0.1) -> None:
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(min(y_positions) - 0.45, max(y_positions) + 0.45)
    ticks = [round(t, 2) for t in np.arange(
        round(math.ceil(xmin * 10) / 10, 1),
        xmax + 0.001,
        xtick_step,
    ) if t <= xmax + 0.001]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.1f}" for t in ticks], fontsize=7)
    for t in ticks:
        ax.axvline(t, color=L_GRAY, linewidth=0.7, zorder=0)
    ax.axvline(0.5, color=L_GRAY, linewidth=0.7, linestyle="--", zorder=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(ylabels, fontsize=8, color=DARK)
    ax.set_xlabel("AUC", fontsize=8, color=GRAY, labelpad=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(M_GRAY)
    ax.tick_params(axis="x", length=3, color=M_GRAY, labelcolor=GRAY)
    ax.tick_params(axis="y", length=0, pad=4)



# Dumbbell helper

def _dumbbell_panel(ax, sub: pd.DataFrame, color: str,
                    xmin: float = 0.40, xmax: float = 1.00,
                    delta_col: bool = True) -> None:
    sub = sub.copy()
    sub["_ord"] = sub["site"].map(SITE_ORDER).fillna(99)
    sub = sub.sort_values("_ord").reset_index(drop=True)
    n   = len(sub)
    _style_data_ax(ax, xmin, xmax, n)
    blend = ax.get_yaxis_transform()   # x in axes [0-1], y in data coords

    for i, row in sub.iterrows():
        y    = n - 1 - i
        caut = "caution" in str(getattr(row, "adequacy", ""))
        raw_c = GRAY
        cen_c = M_GRAY if caut else color

        raw = float(row["raw_auc_95ci_mid"])
        cen = float(row["stratum_centered_auc_95ci_mid"])

        if not math.isnan(raw):
            ax.scatter(raw, y, color=raw_c, s=22, marker="o", zorder=5, linewidths=0)
        if not math.isnan(cen):
            if caut:
                ax.scatter(cen, y, facecolors="white", edgecolors=cen_c,
                           s=26, marker="s", linewidths=1.1, zorder=6)
            else:
                ax.scatter(cen, y, color=cen_c, s=26, marker="s", zorder=6, linewidths=0)
        _errbars(ax, y, float(row["raw_auc_95ci_low"]),
                 float(row["raw_auc_95ci_high"]), raw_c, lw=0.40)
        _errbars(ax, y, float(row["stratum_centered_auc_95ci_low"]),
                 float(row["stratum_centered_auc_95ci_high"]), cen_c, lw=0.55)

        if delta_col:
            delta = float(row.get("raw_to_centered_delta", math.nan))
            ret   = float(row.get("matched_retention_pct", math.nan))
            if not math.isnan(delta):
                label = f"Delta={delta:+.3f};  ret={ret:.1f}%"
                ax.text(1.03, y, label, fontsize=6.0, color=GRAY,
                        ha="left", va="center", transform=blend, clip_on=False)

    ax.set_yticks(list(range(n)))
    ylabels = list(reversed(sub["site"].tolist()))
    ax.set_yticklabels(ylabels, fontsize=8, color=DARK)




def _primary_forest_panel(ax, sub: pd.DataFrame, color: str,
                          xmin: float = 0.40, xmax: float = 1.00) -> None:
    sub = sub.copy()
    sub["_ord"] = sub["site"].map(SITE_ORDER).fillna(99)
    sub = sub.sort_values("_ord").reset_index(drop=True)
    y_positions = list(reversed(range(len(sub))))
    ylabels = list(reversed(sub["site"].tolist()))
    _style_forest_ax(ax, xmin, xmax, y_positions, ylabels)

    for i, row in sub.iterrows():
        y_base = len(sub) - 1 - i
        raw_y = y_base + 0.15
        cen_y = y_base - 0.15
        caut = "caution" in str(getattr(row, "adequacy", ""))
        cen_c = M_GRAY if caut else color
        raw = float(row["raw_auc_95ci_mid"])
        cen = float(row["stratum_centered_auc_95ci_mid"])

        ax.axhspan(y_base - 0.31, y_base + 0.31, color=PALE, zorder=-2)

        if not math.isnan(raw):
            _errbars(ax, raw_y, float(row["raw_auc_95ci_low"]),
                     float(row["raw_auc_95ci_high"]), GRAY, lw=0.45)
            ax.scatter(raw, raw_y, color=GRAY, s=22, marker="o", zorder=5, linewidths=0)
        if not math.isnan(cen):
            _errbars(ax, cen_y, float(row["stratum_centered_auc_95ci_low"]),
                     float(row["stratum_centered_auc_95ci_high"]), cen_c, lw=0.55)
            if caut:
                ax.scatter(cen, cen_y, facecolors="white", edgecolors=cen_c,
                           s=28, marker="s", linewidths=1.1, zorder=6)
            else:
                ax.scatter(cen, cen_y, color=cen_c, s=28, marker="s", zorder=6, linewidths=0)


# Figure 3: Model-class matrix

def _model_class_col(mc: str, mv: str) -> str:
    if mc == "CNN/Mega":
        return "CNN/Mega"
    if mc == "LightGBM" and mv == "multi-task":
        return "LGBM multi"
    if mc == "LightGBM" and mv == "single-task":
        return "LGBM single"
    if mc == "Weis LR":
        return "Weis LR"
    return f"{mc} {mv}".strip()


def _figure_3_matrix_rows(matrix_df: pd.DataFrame) -> pd.DataFrame:
    matrix = matrix_df[matrix_df["status"].eq("complete")].copy()
    targets = [
        ("E. coli", "Ciprofloxacin", "A-2018"),
        ("E. coli", "Ciprofloxacin", "DRIAMS-C"),
        ("E. coli", "Ciprofloxacin", "DRIAMS-D"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "A-2018"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "DRIAMS-C"),
        ("E. coli", "Amoxicillin-Clavulanic acid", "DRIAMS-D"),
        ("S. aureus", "Oxacillin", "A-2018"),
        ("S. aureus", "Oxacillin", "DRIAMS-B"),
        ("S. aureus", "Oxacillin", "DRIAMS-C"),
    ]
    target_index = {target: index for index, target in enumerate(targets)}

    def organism_short(organism: str) -> str:
        if organism == "Escherichia coli":
            return "E. coli"
        if organism == "Staphylococcus aureus":
            return "S. aureus"
        return organism

    matrix["organism_short"] = matrix["organism"].map(organism_short)
    matrix["target_key"] = list(zip(matrix["organism_short"], matrix["drug"], matrix["site"]))
    matrix = matrix[matrix["target_key"].isin(target_index)].copy()
    matrix["target_order"] = matrix["target_key"].map(target_index)
    matrix["display_target"] = matrix.apply(
        lambda row: f"{row['organism_short']} / {DRUG_SHORT.get(str(row['drug']), str(row['drug']))} / {row['site']}",
        axis=1,
    )
    matrix["model_column"] = matrix.apply(
        lambda row: _model_class_col(str(row["model_class"]), str(row["model_variant"])),
        axis=1,
    )
    matrix = matrix[matrix["model_column"].isin(["CNN/Mega", "LGBM multi", "LGBM single", "Weis LR"])].copy()
    for col in ["raw_auc", "centered_auc", "matched_retention"]:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")

    atlas_path = ANALYSIS / "audit_atlas_current" / "atlas_matrix.csv"
    if atlas_path.exists():
        atlas = pd.read_csv(atlas_path)
        weis_temporal = atlas[
            atlas["atlas_id"].eq("weis_lr_ecoli_temporal_a2018")
            & atlas["drug"].isin(["Ciprofloxacin", "Amoxicillin-Clavulanic acid"])
            & atlas["site"].eq("DRIAMS-A")
        ].copy()
        weis_ecoli = atlas[
            atlas["atlas_id"].eq("weis_lr_ecoli")
            & atlas["drug"].isin(["Ciprofloxacin", "Amoxicillin-Clavulanic acid"])
            & atlas["site"].isin(["DRIAMS-C", "DRIAMS-D"])
        ].copy()
        weis = pd.concat([weis_temporal, weis_ecoli], ignore_index=True)
        if not weis.empty:
            weis.loc[weis["site"].eq("DRIAMS-A"), "site"] = "A-2018"
            weis["organism_short"] = "E. coli"
            weis["target_key"] = list(zip(weis["organism_short"], weis["drug"], weis["site"]))
            weis["target_order"] = weis["target_key"].map(target_index)
            weis = weis[weis["target_order"].notna()].copy()
            weis["display_target"] = weis.apply(
                lambda row: f"E. coli / {DRUG_SHORT.get(str(row['drug']), str(row['drug']))} / A-2018",
                axis=1,
            )
            weis.loc[~weis["site"].eq("A-2018"), "display_target"] = weis.loc[
                ~weis["site"].eq("A-2018")
            ].apply(
                lambda row: f"E. coli / {DRUG_SHORT.get(str(row['drug']), str(row['drug']))} / {row['site']}",
                axis=1,
            )
            weis["model_column"] = "Weis LR"
            weis["model_class"] = "Weis LR"
            weis["model_variant"] = weis["atlas_id"].map(
                {
                    "weis_lr_ecoli_temporal_a2018": "temporal A-2018",
                    "weis_lr_ecoli": "Weis-style panel",
                }
            )
            weis["centered_auc"] = pd.to_numeric(weis["stratum_centered_auc"], errors="coerce")
            weis["valid_strata"] = weis["n_valid_strata"]
            weis["adequacy_label"] = weis["audit_adequacy_label"]
            weis["interpretation"] = weis["audit_interpretation_category"]
            weis["source_path"] = "outputs/analysis_outputs/audit_atlas_current/atlas_matrix.csv"
            weis["status"] = "complete"
            keep_cols = [
                "status",
                "model_class",
                "model_variant",
                "scope",
                "organism",
                "drug",
                "site",
                "raw_auc",
                "centered_auc",
                "pairwise_accuracy",
                "matched_retention",
                "n_total",
                "n_r",
                "n_matched",
                "n_matched_r",
                "valid_strata",
                "adequacy_label",
                "interpretation",
                "source_path",
                "notes",
                "organism_short",
                "target_key",
                "target_order",
                "display_target",
                "model_column",
            ]
            matrix = pd.concat([matrix, weis[keep_cols]], ignore_index=True)
    return matrix.sort_values(["target_order", "model_column"])


def _normalize_primary_source(primary: pd.DataFrame) -> pd.DataFrame:
    source = primary.copy()
    mask = source["pair"].eq("E. coli / Amox-Clav") & source["site"].eq("A-2018")
    source.loc[mask, "interpretation"] = "Weak/borderline residual signal"
    return source


def _normalize_figure_3_source(fig3: pd.DataFrame) -> pd.DataFrame:
    source = fig3.copy()
    caution = source["adequacy_label"].astype(str).str.contains("caution", na=False)
    source.loc[caution, "interpretation"] = "Cautionary; low matched support"
    return source


def _write_source_data(primary: pd.DataFrame, matrix_df: pd.DataFrame) -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    _normalize_primary_source(primary).to_csv(
        SOURCE_DIR / "source_data_fig2_primary_background_audit.csv",
        index=False,
    )
    if not matrix_df.empty:
        fig3 = _figure_3_matrix_rows(matrix_df)
        fig3 = fig3[
            [
                "display_target",
                "model_column",
                "model_class",
                "model_variant",
                "scope",
                "organism",
                "drug",
                "site",
                "raw_auc",
                "centered_auc",
                "pairwise_accuracy",
                "matched_retention",
                "valid_strata",
                "adequacy_label",
                "interpretation",
                "source_path",
            ]
        ]
        _normalize_figure_3_source(fig3).to_csv(
            SOURCE_DIR / "source_data_fig3_model_family_replication.csv",
            index=False,
        )
    atlas_summary = ANALYSIS / "audit_atlas_current" / "atlas_summary.csv"
    if atlas_summary.exists():
        pd.read_csv(atlas_summary).to_csv(
            SOURCE_DIR / "source_data_fig8_audit_atlas_summary.csv",
            index=False,
        )






# Figure 6: Three-way decomposition

def _build_three_way(primary_df, sa_summary_df, ecoli_bg_df, saureus_bg_df):
    data = _prep_ci(primary_df)
    rows = []
    for pair_label, drug_long in [
        ("E. coli / Cipro",     "Ciprofloxacin"),
        ("E. coli / Amox-Clav", "Amoxicillin-Clavulanic acid"),
    ]:
        sub    = data[data["pair"].eq(pair_label)].copy()
        bg_sub = ecoli_bg_df[ecoli_bg_df["drug"].eq(drug_long)]
        for _, row in sub.iterrows():
            site   = row["site"]
            bg_row = bg_sub[bg_sub["site"].eq(site)]
            bg_auc = float(bg_row["exact_background_auc"].iloc[0]) if len(bg_row) else math.nan
            rows.append({
                "pair": pair_label, "site": site,
                "raw_auc": row["raw_auc_95ci_mid"],
                "raw_lo":  row["raw_auc_95ci_low"],  "raw_hi": row["raw_auc_95ci_high"],
                "bg_auc":  bg_auc,
                "cen_auc": row["stratum_centered_auc_95ci_mid"],
                "cen_lo":  row["stratum_centered_auc_95ci_low"],
                "cen_hi":  row["stratum_centered_auc_95ci_high"],
                "caution": "caution" in str(row["adequacy"]),
            })
    sa_oxa = sa_summary_df[sa_summary_df["drug"].eq("Oxacillin")]
    bg_sa  = saureus_bg_df[saureus_bg_df["drug"].eq("Oxacillin")]
    for _, row in sa_oxa.iterrows():
        site   = row["site"]
        bg_row = bg_sa[bg_sa["site"].eq(site)]
        bg_auc = float(bg_row["exact_background_auc"].iloc[0]) if len(bg_row) else math.nan
        n_str  = int(row["n_valid_strata"])
        has_s  = n_str > 0
        rows.append({
            "pair": "S. aureus / Oxacillin", "site": site,
            "raw_auc": float(row["raw_auc"]),
            "raw_lo":  float(row["raw_auc_ci_low"]), "raw_hi": float(row["raw_auc_ci_high"]),
            "bg_auc":  bg_auc,
            "cen_auc": float(row["stratum_centered_auc"])        if has_s else math.nan,
            "cen_lo":  float(row["stratum_centered_auc_ci_low"]) if has_s else math.nan,
            "cen_hi":  float(row["stratum_centered_auc_ci_high"]) if has_s else math.nan,
            "caution": n_str <= 1 or "caution" in str(row["adequacy_label"]),
        })
    atlas_path = ANALYSIS / "audit_atlas_current" / "atlas_matrix.csv"
    if atlas_path.exists():
        atlas = pd.read_csv(atlas_path)
        kleb = atlas[
            atlas["atlas_id"].eq("klebsiella_cnn")
            & atlas["drug"].eq("Ceftriaxone")
        ].copy()
        for _, row in kleb.iterrows():
            centered = float(row["stratum_centered_auc"]) if pd.notna(row.get("stratum_centered_auc")) else math.nan
            n_str = int(float(row.get("n_valid_strata", 0) or 0))
            rows.append({
                "pair": "K. pneumoniae / Ceftriaxone",
                "site": row["site"],
                "raw_auc": float(row["raw_auc"]) if pd.notna(row.get("raw_auc")) else math.nan,
                "raw_lo": math.nan,
                "raw_hi": math.nan,
                "bg_auc": float(row["exact_background_auc"]) if pd.notna(row.get("exact_background_auc")) else math.nan,
                "cen_auc": centered,
                "cen_lo": math.nan,
                "cen_hi": math.nan,
                "caution": n_str <= 1 or "caution" in str(row.get("audit_adequacy_label", "")),
            })
    return pd.DataFrame(rows)



# Figure 7: Framework comparison table

COMPARISON_ROWS = [
    ("Reports external-site AUC",              "Yes",     "Required",      "Yes",  ""),
    ("Controls for co-resistance background",  "No",      "Not required",  "Yes",  "Core audit metric"),
    ("Quantifies signal collapse vs. retention","No",     "No",            "Yes",  "raw - centered delta"),
    ("Site-specific evaluation",               "Yes",     "Recommended",   "Yes",  ""),
    ("Multi-drug panel (> 2 pairs)",           "Some",    "N/A",           "Yes",  "saureus_panel, ecoli_mechanism6"),
    ("Second organism validation",             "No",      "N/A",           "Yes",  "E. coli + S. aureus"),
    ("Sensitivity analysis (stratum threshold)","No",     "No",            "Yes",  "2, 3, 5, 10 isolates per stratum"),
    ("Model-agnostic (accepts any CSV)",       "No",      "N/A",           "Yes",  "No DRIAMS/PyTorch dependency"),
    ("Open code + reproducible example data",  "Partial", "N/A",           "Yes",  "github + SCHEMA.md"),
]




def _clean_header(fig, title: str, subtitle: str = "") -> None:
    fig.text(0.045, 0.93, title, fontsize=12.5, fontweight="bold",
             color=DARK, va="top", ha="left")
    if subtitle:
        fig.text(0.045, 0.875, subtitle, fontsize=8.0, color=GRAY,
                 va="top", ha="left")


def _clean_auc_axis(ax, xmin=0.45, xmax=0.95, rows=4) -> None:
    _style_data_ax(ax, xmin, xmax, rows)
    ax.grid(False)
    ax.set_axisbelow(True)


def figure_1_framework() -> Path:
    path = FIG_DIR / "figure_1_framework.pdf"
    fig, ax = plt.subplots(figsize=(8.6, 3.05))
    fig.subplots_adjust(left=0.05, right=0.98, top=0.72, bottom=0.12)
    _clean_header(
        fig,
        "Fig. 1 | Background-matched MALDI-AMR audit",
        "Model-agnostic test of whether focal-drug prediction survives co-resistance background control.",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    boxes = [
        ("Model\npredictions", "isolate, site/year,\norganism, drug,\nlabel, score"),
        ("Raw\ntransfer", "external AUC/AUPR\nbefore background\ncontrol"),
        ("Background\nmatch", "stratify by labels\nof the other\nantibiotics"),
        ("Centered\naudit", "subtract each\nstratum mean\nmodel score"),
        ("Interpretation", "retained signal,\npartial retention,\nor collapse"),
    ]
    n = len(boxes)
    gap = 0.026
    bw = (1 - gap * (n - 1)) / n
    by, bh = 0.04, 0.80
    for i, (title, body) in enumerate(boxes):
        x = i * (bw + gap)
        ax.add_patch(FancyBboxPatch(
            (x, by), bw, bh, boxstyle="round,pad=0.012",
            facecolor=PALE, edgecolor=L_GRAY, linewidth=0.9,
        ))
        ax.text(x + bw / 2, by + bh - 0.16, title, ha="center", va="top",
                fontsize=7.5, fontweight="bold", color=DARK, linespacing=1.15)
        ax.text(x + bw / 2, by + 0.25, body, ha="center", va="center",
                fontsize=6.3, color=GRAY, linespacing=1.18)
        if i < n - 1:
            ax.annotate("", xy=(x + bw + gap * 0.82, by + bh / 2),
                        xytext=(x + bw + gap * 0.18, by + bh / 2),
                        arrowprops=dict(arrowstyle="->", color=M_GRAY, lw=1.0))
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_2_primary_audit(primary: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_2_primary_background_audit.pdf"
    data = _prep_ci(primary)
    panels = [("E. coli / Cipro", BLUE), ("E. coli / Amox-Clav", ORANGE)]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.15))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.78, bottom=0.21, wspace=0.36)
    _clean_header(
        fig,
        "Fig. 2 | Focal-drug signal after background matching",
        "Ciprofloxacin retained partial within-background signal; amoxicillin-clavulanic acid did not.",
    )
    for ax, (pair_label, color) in zip(axes, panels):
        sub = data[data["pair"].eq(pair_label)].copy()
        _primary_forest_panel(ax, sub, color)
        ax.set_title(pair_label, loc="left", fontsize=10, fontweight="bold",
                     color=DARK, pad=8)
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GRAY,
               markersize=6, label="Raw AUC"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=DARK,
               markersize=6, label="Background-centered AUC"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="white",
               markeredgecolor=M_GRAY, markeredgewidth=1.1,
               markersize=6, label="Sparse matched support"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.52, 0.055))
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_3_model_replication(model_df: pd.DataFrame,
                               matrix_df: pd.DataFrame | None = None) -> Path:
    path = FIG_DIR / "figure_3_model_family_replication.pdf"
    df = model_df.copy()
    if matrix_df is not None and not matrix_df.empty:
        df = _figure_3_matrix_rows(matrix_df)
    elif "display_target" not in df.columns:
        source_path = SOURCE_DIR / "source_data_fig3_model_family_replication.csv"
        if source_path.exists():
            df = pd.read_csv(source_path)
    target_order = [
        "E. coli / Cipro / A-2018", "E. coli / Cipro / DRIAMS-C",
        "E. coli / Cipro / DRIAMS-D", "E. coli / Amox-Clav / A-2018",
        "E. coli / Amox-Clav / DRIAMS-C", "E. coli / Amox-Clav / DRIAMS-D",
        "S. aureus / Oxa / A-2018", "S. aureus / Oxa / DRIAMS-B",
        "S. aureus / Oxa / DRIAMS-C",
    ]
    columns = ["CNN/Mega", "LGBM multi", "LGBM single", "Weis LR"]
    df = df[df["display_target"].isin(target_order) & df["model_column"].isin(columns)].copy()
    pivot = df.set_index(["display_target", "model_column"])

    fig, ax = plt.subplots(figsize=(10.6, 4.9))
    fig.subplots_adjust(left=0.03, right=0.99, top=0.82, bottom=0.08)
    _clean_header(
        fig,
        "Fig. 3 | Completed model-class background-audit matrix",
        "Background sensitivity is consistent across model architectures.",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    x_label = 0.24
    x0 = 0.29
    col_w = 0.17
    row_h = 0.076
    y_top = 0.88
    ax.text(0.04, y_top, "Target / site", fontsize=8, fontweight="bold",
            color=DARK, va="center")
    for j, col in enumerate(columns):
        ax.text(x0 + j * col_w + col_w / 2, y_top, col, fontsize=8,
                fontweight="bold", color=DARK, ha="center", va="center")
    for i, target in enumerate(target_order):
        y = y_top - 0.075 - i * row_h
        if i in {3, 6}:
            ax.plot([0.03, 0.98], [y + row_h * 0.52, y + row_h * 0.52],
                    color=L_GRAY, lw=0.8)
        ax.text(0.04, y, target, fontsize=7.2, color=DARK, va="center")
        for j, col in enumerate(columns):
            x = x0 + j * col_w
            rec = pivot.loc[(target, col)] if (target, col) in pivot.index else None
            if rec is None:
                face = PALE
                label = "not run"
                text_color = M_GRAY
                edge_color = "white"
            else:
                raw = float(rec["raw_auc"])
                cen = float(rec["centered_auc"])
                caution = "caution" in str(rec.get("adequacy_label", ""))
                if math.isnan(raw) or math.isnan(cen):
                    face = PALE
                    label = "n/a"
                    text_color = M_GRAY
                else:
                    face = _auc_cell_color(cen)
                    label = f"{raw:.2f} -> {cen:.2f}"
                    text_color = DARK
                if caution:
                    label += "*"
                edge_color = M_GRAY if caution else "white"
            ax.add_patch(mpatches.Rectangle((x, y - row_h * 0.42), col_w - 0.004,
                                            row_h * 0.84, facecolor=face,
                                            edgecolor=edge_color,
                                            lw=1.0))
            ax.text(x + col_w / 2, y, label, fontsize=7.0, color=text_color,
                    ha="center", va="center")
    caption = (
        "Cells show raw AUC followed by background-centered AUC for CNN/Mega, multi-task LightGBM, "
        "single-task LightGBM and Weis/Borgwardt LR rows. CNN and LightGBM use the expanded study panels. "
        "Weis/Borgwardt rows are shown separately: exact upstream parity was verified only for the official LR "
        "subset, while other Weis-format rows test audit portability rather than serve as pooled benchmarks. "
        "Asterisks mark low matched support; retention and valid-stratum counts are reported in the source data."
    )
    ax.text(
        0.04,
        0.02,
        "\n".join(textwrap.wrap(caption, 170)),
        fontsize=6.2,
        color=GRAY,
        va="bottom",
    )
    note = (
        "* low matched support; N/A = no usable S/R oxacillin-labeled S. aureus prediction rows in "
        "DRIAMS-D export; 'not run' denotes no corresponding Weis/Borgwardt LR audit row for that target/site."
    )
    ax.text(
        0.70,
        0.02,
        "\n".join(textwrap.wrap(note, 72)),
        fontsize=5.9,
        color=GRAY,
        va="bottom",
        ha="left",
    )
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_4_cross_resistance() -> Path:
    path = FIG_DIR / "figure_4_cross_resistance_network.pdf"
    edges_df = pd.read_csv(ANALYSIS / "cross_resistance_network" / "cross_resistance_edges.csv")
    all_e = edges_df[edges_df["site"].eq("ALL")]
    mat = pd.DataFrame(index=DRUG_ORDER, columns=DRUG_ORDER, dtype=float)
    for d in DRUG_ORDER:
        mat.loc[d, d] = 1.0
    for _, row in all_e.iterrows():
        a = DRUG_SHORT.get(str(row["drug_a"]), str(row["drug_a"]))
        b = DRUG_SHORT.get(str(row["drug_b"]), str(row["drug_b"]))
        if a in mat.index and b in mat.columns:
            mat.loc[a, b] = row["phi"]
            mat.loc[b, a] = row["phi"]

    fig = plt.figure(figsize=(10.2, 4.7))
    _clean_header(
        fig,
        "Fig. 4 | Co-resistance blocks define exploitable background",
        "Phi correlations across E. coli resistance labels show structured resistance ecology.",
    )
    ax = fig.add_axes([0.07, 0.14, 0.54, 0.66])
    im = ax.imshow(mat.values.astype(float), vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks(range(len(DRUG_ORDER)))
    ax.set_yticks(range(len(DRUG_ORDER)))
    ax.set_xticklabels(DRUG_ORDER, fontsize=8)
    ax.set_yticklabels(DRUG_ORDER, fontsize=8)
    ax.tick_params(length=0)
    ax.set_title("Phi co-resistance correlation", loc="left", fontsize=9,
                 fontweight="bold", pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(len(DRUG_ORDER)):
        for j in range(len(DRUG_ORDER)):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=7, color="white" if val > 0.55 else DARK)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("phi", fontsize=7, color=GRAY)

    ax2 = fig.add_axes([0.69, 0.18, 0.27, 0.58])
    ax2.axis("off")
    ax2.text(0, 1.0, "Strongest blocks", fontsize=9, fontweight="bold",
             color=DARK, va="top")
    top = all_e.sort_values("phi", ascending=False).head(5).reset_index(drop=True)
    for i, row in top.iterrows():
        y = 0.84 - i * 0.15
        label = f"{DRUG_SHORT.get(str(row['drug_a']), row['drug_a'])} / {DRUG_SHORT.get(str(row['drug_b']), row['drug_b'])}"
        ax2.text(0, y, label, fontsize=8, fontweight="bold", color=DARK, va="top")
        ax2.text(0, y - 0.065, f"phi={row['phi']:.3f}   RR lift={row['rr_lift']:.2f}",
                 fontsize=7.5, color=GRAY, va="top")
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_5_public_support(wgs: pd.DataFrame, enrichment: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_5_public_wgs_proteomic_support.pdf"
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 3.95))
    fig.subplots_adjust(left=0.09, right=0.95, top=0.70, bottom=0.21, wspace=0.52)
    _clean_header(
        fig,
        "Fig. 5 | Public WGS-linked MALDI data support lineage encoding",
        "Public Basel UPEC spectra link MALDI peak features, WGS-derived lineage labels, and resistance phenotypes.",
    )
    ax = axes[0]
    wgs_sorted = wgs.sort_values("auc", ascending=True)
    labels = [str(x).replace("_", " ") for x in wgs_sorted["target"]]
    colors = [BLUE if x == "ST131" else ORANGE for x in wgs_sorted["target"]]
    y = np.arange(len(wgs_sorted))
    ax.barh(y, wgs_sorted["auc"] - 0.5, left=0.5, color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(0.5, 1.0)
    ax.set_xlabel("AUC", fontsize=8, color=GRAY)
    ax.set_title("a  MALDI peak-feature AUC", loc="left", fontsize=9,
                 fontweight="bold", color=DARK, pad=8)
    ax.grid(axis="x", color=L_GRAY, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for yi, val in zip(y, wgs_sorted["auc"]):
        ax.text(val + 0.004, yi, f"{val:.3f}", fontsize=7.5,
                fontweight="bold", va="center", color=DARK)

    ax2 = axes[1]
    show = enrichment[enrichment["target"].isin(
        ["ST131", "Ciprofloxacin_R", "Ceftriaxone_R"])].copy()
    show = show.sort_values("fold_enrichment", ascending=True)
    labels = [str(x).replace("_", " ") for x in show["target"]]
    colors = [BLUE if x == "ST131" else ORANGE for x in show["target"]]
    y = np.arange(len(show))
    ax2.barh(y, show["fold_enrichment"], color=colors, height=0.55)
    ax2.axvline(1, color=M_GRAY, linestyle="--", lw=0.9)
    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, fontsize=8)
    ax2.set_xlim(0, max(3.9, show["fold_enrichment"].max() + 0.75))
    ax2.set_xlabel("Fold enrichment", fontsize=8, color=GRAY)
    ax2.set_title("b  Overlap with published ST131 biomarkers", loc="left",
                  fontsize=9, fontweight="bold", color=DARK, pad=8)
    ax2.grid(axis="x", color=L_GRAY, linewidth=0.7)
    ax2.set_axisbelow(True)
    ax2.spines[["top", "right", "left"]].set_visible(False)
    for yi, (_, row) in zip(y, show.iterrows()):
        ax2.text(row["fold_enrichment"] + 0.03, yi,
                 f"{row['fold_enrichment']:.2f}x; p={row['empirical_p_ge_observed']:.4f}",
                 fontsize=7.4, va="center", color=DARK)
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_6_three_way_decomposition(primary_df, sa_summary_df,
                                      ecoli_bg_df, saureus_bg_df) -> Path:
    path = FIG_DIR / "figure_8_three_way_decomposition.pdf"
    df = _build_three_way(primary_df, sa_summary_df, ecoli_bg_df, saureus_bg_df)
    df.rename(columns={"bg_auc": "bg_only_auc", "cen_auc": "centered_auc"}).to_csv(
        SOURCE_DIR / "source_data_fig6_three_way_decomposition.csv",
        index=False,
    )
    panels_cfg = [
        ("E. coli / Cipro", "Focal signal retained"),
        ("E. coli / Amox-Clav", "Background-sensitive"),
        ("S. aureus / Oxacillin", "Signal retained"),
        ("K. pneumoniae / Ceftriaxone", "Third-organism sparse control"),
    ]
    metrics = [
        ("raw_auc", "Raw\nMALDI"),
        ("bg_auc", "Co-resistance\nonly"),
        ("cen_auc", "Background-\ncentered"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13.8, 4.7))
    fig.subplots_adjust(left=0.07, right=0.92, top=0.72, bottom=0.16, wspace=0.34)
    _clean_header(
        fig,
        "Fig. 6 | Four-panel three-score AUC audit",
        "Each panel reports AUC for the same organism-drug pair under raw MALDI, co-resistance-only, and background-centered scoring.",
    )
    im = None
    for ax, (pair_label, subtitle) in zip(axes, panels_cfg):
        sub = df[df["pair"].eq(pair_label)].copy()
        sub["_ord"] = sub["site"].map(SITE_ORDER).fillna(99)
        sub = sub.sort_values("_ord").reset_index(drop=True)
        vals = sub[[m[0] for m in metrics]].to_numpy(float)
        im = ax.imshow(vals, vmin=0.45, vmax=1.0, cmap="Blues", aspect="auto")
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels([m[1] for m in metrics], fontsize=7.4)
        site_labels = [
            f"{site}{'*' if caut else ''}"
            for site, caut in zip(sub["site"].tolist(), sub["caution"].tolist())
        ]
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(site_labels, fontsize=8)
        ax.tick_params(length=0)
        ax.set_title(pair_label, loc="left", fontsize=10, fontweight="bold",
                     color=DARK, pad=12)
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=7.2,
                color=GRAY, va="bottom", ha="left")
        for spine in ax.spines.values():
            spine.set_visible(False)
        for y in range(vals.shape[0]):
            for x in range(vals.shape[1]):
                v = vals[y, x]
                label = "--" if math.isnan(v) else f"{v:.2f}"
                ax.text(x, y, label, ha="center", va="center",
                        fontsize=8, fontweight="bold" if x == 2 else "normal",
                        color="white" if (not math.isnan(v) and v >= 0.74) else DARK)
        ax.set_xticks(np.arange(-0.5, len(metrics), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(sub), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.4)
        ax.tick_params(which="minor", bottom=False, left=False)
    if im is not None:
        cax = fig.add_axes([0.935, 0.22, 0.015, 0.46])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("AUC", fontsize=7, color=GRAY)
        cb.ax.tick_params(labelsize=7, length=2)
    fig.text(0.08, 0.07, "* sparse matched support", fontsize=6.8, color=GRAY)
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_6_falsification_controls(df: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_6_falsification_controls.pdf"
    df = df.copy()
    drug_order = ["Cipro", "Norflox", "Amox-Clav", "CRO", "CAZ", "FEP"]
    sites = ["A-2018", "DRIAMS-B", "DRIAMS-C", "DRIAMS-D"]
    df["drug_short"] = df["drug"].map(lambda x: DRUG_SHORT.get(str(x), str(x)))
    mat = np.full((len(sites), len(drug_order)), np.nan)
    pvals = np.full_like(mat, np.nan)
    for _, row in df.iterrows():
        if row["site"] in sites and row["drug_short"] in drug_order:
            mat[sites.index(row["site"]), drug_order.index(row["drug_short"])] = row["observed_minus_burden_auc"]
            pvals[sites.index(row["site"]), drug_order.index(row["drug_short"])] = row["shuffle_empirical_p_ge_observed"]
    fig = plt.figure(figsize=(9.2, 4.25))
    _clean_header(
        fig,
        "Fig. 6 supplement | Falsification controls",
        "Observed AUC compared with background-burden and shuffled-label controls.",
    )
    ax = fig.add_axes([0.08, 0.16, 0.58, 0.62])
    masked = np.ma.masked_invalid(mat)
    cmap = plt.get_cmap("RdBu").copy()
    cmap.set_bad("#F3F4F6")
    im = ax.imshow(masked, vmin=-0.35, vmax=0.15, cmap=cmap)
    ax.set_xticks(range(len(drug_order)))
    ax.set_xticklabels(drug_order, fontsize=8)
    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites, fontsize=8)
    ax.set_title("Observed AUC minus background-burden AUC", loc="left",
                 fontsize=9, fontweight="bold", color=DARK, pad=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                star = "*" if np.isfinite(pvals[i, j]) and pvals[i, j] < 0.05 else ""
                ax.text(j, i, f"{mat[i, j]:+.2f}{star}", ha="center", va="center",
                        fontsize=6.8, color=DARK)
            else:
                ax.text(j, i, "n/a", ha="center", va="center",
                        fontsize=6.8, color=GRAY)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cb.ax.tick_params(labelsize=7, length=2)
    ax_note = fig.add_axes([0.72, 0.25, 0.23, 0.45])
    ax_note.axis("off")
    ax_note.text(0, 1, "Interpretation", fontsize=9, fontweight="bold",
                 color=DARK, va="top")
    note = ("Negative values mean the background-burden score is competitive "
            "with the MALDI model. Asterisks mark p < 0.05 vs shuffled labels.")
    ax_note.text(0, 0.82, "\n".join(textwrap.wrap(note, 38)), fontsize=7.3,
                 color=GRAY, va="top", linespacing=1.35)
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_6_saureus_oxacillin_audit(df: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_7_saureus_oxacillin_audit.pdf"
    df = df[df["pair"].astype(str).str.contains("S. aureus", regex=False)].copy()
    df = df.sort_values("cnn_minus_bg")
    colors = [GREEN if v > 0 else M_GRAY for v in df["cnn_minus_bg"]]
    fig, ax = plt.subplots(figsize=(7.9, 3.8))
    fig.subplots_adjust(left=0.20, right=0.93, top=0.70, bottom=0.20)
    _clean_header(
        fig,
        "Fig. 6 supplement | S. aureus oxacillin audit",
        "MALDI/CNN AUC compared with exact-background-only prediction.",
    )
    y = np.arange(len(df))
    ax.barh(y, df["cnn_minus_bg"], color=colors, edgecolor="none", height=0.58)
    ax.axvline(0, color=M_GRAY, linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(df["site"].astype(str), fontsize=8)
    ax.set_xlabel("CNN AUC - background-only AUC", fontsize=8, color=GRAY)
    ax.set_title("Residual MALDI advantage by site", loc="left", fontsize=9,
                 fontweight="bold", color=DARK, pad=8)
    xmax = max(0.11, float(df["cnn_minus_bg"].max()) + 0.02)
    xmin = min(-0.035, float(df["cnn_minus_bg"].min()) - 0.02)
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="x", color=L_GRAY, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for yi, val in zip(y, df["cnn_minus_bg"]):
        ha = "left" if val >= 0 else "right"
        dx = 0.004 if val >= 0 else -0.004
        ax.text(val + dx, yi, f"{val:+.2f}", ha=ha, va="center",
                fontsize=7.5, color=DARK)
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_7_deployment_decision_flow(df: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_9_deployment_decision_flow.pdf"
    fig = plt.figure(figsize=(10.8, 5.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    _clean_header(
        fig,
        "Fig. 7 supplement | Deployment decision flow",
        "Audit outcomes translated into conservative validation and deployment actions.",
    )
    xs = [0.05, 0.37, 0.62]
    headers = ["Audit pattern", "Decision category", "Recommended action"]
    widths = [34, 24, 54]
    for x, label in zip(xs, headers):
        ax.text(x, 0.795, label, fontsize=8, fontweight="bold",
                color=DARK, transform=ax.transAxes)
    status_colors = [GREEN, BLUE, ORANGE, M_GRAY, RED]
    row_h = 0.105
    for i, row in df.reset_index(drop=True).iterrows():
        y = 0.665 - i * 0.125
        ax.add_patch(mpatches.Rectangle((0.04, y), 0.92, row_h,
                                        facecolor=PALE if i % 2 == 0 else "white",
                                        edgecolor="none", transform=ax.transAxes))
        col = status_colors[i % len(status_colors)]
        ax.add_patch(mpatches.Rectangle((0.36, y + 0.02), 0.012, row_h - 0.04,
                                        facecolor=col, edgecolor="none",
                                        transform=ax.transAxes))
        texts = [
            str(row["scenario"]),
            str(row["decision_category"]).replace("_", " "),
            str(row["recommended_action"]),
        ]
        for x, txt, wrap_width, color, weight in [
            (xs[0], texts[0], widths[0], DARK, "normal"),
            (xs[1] + 0.02, texts[1], widths[1], col, "bold"),
            (xs[2], texts[2], widths[2], GRAY, "normal"),
        ]:
            ax.text(x, y + row_h - 0.03, "\n".join(textwrap.wrap(txt, wrap_width)),
                    fontsize=7.0, color=color, fontweight=weight, va="top",
                    transform=ax.transAxes, linespacing=1.22)
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


def figure_8_audit_atlas_summary(summary: pd.DataFrame) -> Path:
    path = FIG_DIR / "figure_11_audit_atlas_summary.pdf"
    df = summary.copy()
    preferred = [
        "ecoli_cnn",
        "ecoli_lgbm_multi",
        "ecoli_lgbm_single",
        "weis_lr_ecoli_temporal_a2018",
        "saureus_cnn",
        "saureus_lgbm_multi",
        "saureus_lgbm_single",
        "klebsiella_cnn",
        "weis_lr_klebsiella_ceftriaxone",
        "weis_lr_ecoli",
    ]
    rank = {name: i for i, name in enumerate(preferred)}
    df["_rank"] = df["atlas_id"].map(rank).fillna(len(rank))
    df = df.sort_values(["_rank", "atlas_id"]).reset_index(drop=True)
    numeric_cols = [
        "n_matrix_rows",
        "n_interpretable",
        "n_retained_or_partial",
        "n_background_sensitive",
        "n_not_interpretable",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    retained = df["n_retained_or_partial"].to_numpy(float)
    background = df["n_background_sensitive"].to_numpy(float)
    sparse = df["n_not_interpretable"].to_numpy(float)
    other = np.maximum(
        df["n_matrix_rows"].to_numpy(float) - retained - background - sparse,
        0,
    )
    y = np.arange(len(df))
    labels = [
        f"{model}\n{scope}"
        for model, scope in zip(df["model_name"].astype(str), df["scope"].astype(str))
    ]

    fig, ax = plt.subplots(figsize=(11.2, 5.8))
    fig.subplots_adjust(left=0.31, right=0.96, top=0.78, bottom=0.16)
    _clean_header(
        fig,
        "Fig. 8 | Audit atlas coverage and outcomes",
        "Each bar is one prediction set; segments summarize retained/partial signal, background-sensitive rows, sparse rows, and other interpretable outcomes.",
    )
    ax.barh(y, retained, color=GREEN, edgecolor="none", height=0.58, label="retained/partial")
    ax.barh(y, background, left=retained, color=ORANGE, edgecolor="none", height=0.58, label="background-sensitive")
    ax.barh(y, sparse, left=retained + background, color=L_GRAY, edgecolor="none", height=0.58, label="sparse/not interpretable")
    ax.barh(y, other, left=retained + background + sparse, color=BLUE, alpha=0.55, edgecolor="none", height=0.58, label="other")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlabel("Number of atlas rows", fontsize=8, color=GRAY)
    ax.grid(axis="x", color=L_GRAY, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=7.5)
    xmax = max(8, int(df["n_matrix_rows"].max()) + 3)
    ax.set_xlim(0, xmax)
    for yi, total, interp in zip(y, df["n_matrix_rows"], df["n_interpretable"]):
        ax.text(total + 0.35, yi, f"{int(total)} rows / {int(interp)} interp.",
                fontsize=6.9, color=GRAY, va="center")
    handles, legend_labels = ax.get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=4,
               frameon=False, fontsize=7.3, bbox_to_anchor=(0.57, 0.045))
    fig.savefig(path, format="pdf", dpi=300)
    plt.close(fig)
    return path


# Main

def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    primary    = pd.read_csv(FINAL / "table_1_primary_background_matched_audit.csv")
    model_df   = pd.read_csv(FINAL / "table_2_cnn_vs_lgbm_multi_background_audit.csv")
    wgs        = pd.read_csv(FINAL / "table_7_public_wgs_maldi_auc.csv")
    enrichment = pd.read_csv(FINAL / "table_9_published_st131_biomarker_enrichment.csv")
    falsification = pd.read_csv(SOURCE_DIR / "source_data_fig6_falsification_controls.csv")
    saureus_oxa = pd.read_csv(SOURCE_DIR / "source_data_fig6_saureus_oxa_audit.csv")
    deployment_flow = pd.read_csv(SOURCE_DIR / "source_data_fig7_deployment_decision_flow.csv")
    atlas_summary_path = ANALYSIS / "audit_atlas_current" / "atlas_summary.csv"
    atlas_summary = pd.read_csv(atlas_summary_path) if atlas_summary_path.exists() else pd.DataFrame()
    matrix_path = ANALYSIS / "model_class_matrix" / "model_class_matrix.csv"
    matrix_df   = pd.read_csv(matrix_path) if matrix_path.exists() else pd.DataFrame()
    _write_source_data(primary, matrix_df)

    sa_path  = ANALYSIS / "saureus_panel_oxa_background_audit" / "background_matched_audit_summary.csv"
    ebg_path = ANALYSIS / "co_resistance_only_baseline_ecoli"   / "co_resistance_only_baseline.csv"
    sbg_path = ANALYSIS / "co_resistance_only_baseline_saureus" / "co_resistance_only_baseline.csv"

    created = []
    for name, fn in [
        ("figure_1", lambda: figure_1_framework()),
        ("figure_2", lambda: figure_2_primary_audit(primary)),
        ("figure_3", lambda: figure_3_model_replication(model_df, matrix_df)),
        ("figure_4", lambda: figure_4_cross_resistance()),
        ("figure_5", lambda: figure_5_public_support(wgs, enrichment)),
        ("figure_6_falsification", lambda: figure_6_falsification_controls(falsification)),
        ("figure_6_saureus_oxa", lambda: figure_6_saureus_oxacillin_audit(saureus_oxa)),
        ("figure_7_deployment", lambda: figure_7_deployment_decision_flow(deployment_flow)),
        ("figure_8_atlas", lambda: figure_8_audit_atlas_summary(atlas_summary)),
    ]:
        if name == "figure_8_atlas" and atlas_summary.empty:
            print("  [figure_8_atlas] atlas_summary.csv missing - skipped.")
            continue
        print(f"Generating {name} ...")
        p = fn()
        created.append(p)
        print(f"  -> {p.relative_to(ROOT)}")

    if all(p.exists() for p in [sa_path, ebg_path, sbg_path]):
        print("Generating figure_6 ...")
        sa  = pd.read_csv(sa_path)
        ebg = pd.read_csv(ebg_path)
        sbg = pd.read_csv(sbg_path)
        p6  = figure_6_three_way_decomposition(primary, sa, ebg, sbg)
        created.append(p6)
        print(f"  -> {p6.relative_to(ROOT)}")
    else:
        print("  [figure_6] input files missing - skipped.")

    print("\nDone. Created:")
    for p in created:
        data = open(p, "rb").read(3000).decode("latin-1", errors="ignore")
        import re as re2
        m = re2.search(r"/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]", data)
        dims = ""
        if m:
            w = (float(m.group(3)) - float(m.group(1))) / 72
            h = (float(m.group(4)) - float(m.group(2))) / 72
            dims = f"  ({w:.2f}\" x {h:.2f}\")"
        print(f"  {p.relative_to(ROOT)}{dims}")


if __name__ == "__main__":
    main()
