#!/usr/bin/env python3
"""Create manuscript tables from the current audit atlas outputs."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_DIR = ROOT / "outputs" / "analysis_outputs" / "audit_atlas_current"
TABLE_DIR = ROOT / "manuscript" / "tables"
SOURCE_DIR = ROOT / "manuscript" / "source_data"
WEIS_SAUREUS_OXA = (
    ROOT
    / "outputs"
    / "analysis_outputs"
    / "weis_lr_official_panel_parity"
    / "audit_summaries"
    / "saureus_oxacillin_background_matched_audit_summary.csv"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in rows])


def latex_escape(value: object) -> str:
    text = "" if value is None else str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def number(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def fmt_auc(value: str) -> str:
    parsed = number(value)
    return "--" if math.isnan(parsed) else f"{parsed:.3f}"


def fmt_pct(value: str) -> str:
    parsed = number(value)
    return "--" if math.isnan(parsed) else f"{100 * parsed:.1f}"


def short_drug(drug: str) -> str:
    return {
        "Amoxicillin-Clavulanic acid": "Amox-Clav",
        "Ciprofloxacin": "Cipro",
        "Norfloxacin": "Norflox",
        "Ceftriaxone": "CRO",
        "Ceftazidime": "CAZ",
        "Cefepime": "FEP",
        "Oxacillin": "Oxa",
    }.get(drug, drug)


def short_scope(scope: str) -> str:
    return (
        scope.replace("E. coli", r"\textit{E. coli}")
        .replace("S. aureus", r"\textit{S. aureus}")
        .replace("K. pneumoniae", r"\textit{K. pneumoniae}")
    )


def latex_table(
    path: Path,
    caption: str,
    label: str,
    headers: list[str],
    rows: list[list[str]],
    align: str,
    notes: str,
    font_size: str = r"\small",
) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        font_size,
        rf"\caption{{{latex_escape(caption)}}}",
        rf"\label{{{label}}}",
        r"\resizebox{\linewidth}{!}{%",
        rf"\begin{{tabular}}{{{align}}}",
        r"\toprule",
        " & ".join(latex_escape(header) for header in headers) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(" & ".join(latex_escape(cell) for cell in row) + r" \\")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}%",
            r"}",
            rf"\vspace{{2mm}}\parbox{{0.95\linewidth}}{{\footnotesize {latex_escape(notes)}}}",
            r"\end{table}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def model_label(atlas_id: str, model_name: str) -> str:
    if atlas_id == "ecoli_cnn":
        return "CNN/Mega"
    if atlas_id == "ecoli_lgbm_multi":
        return "LGBM multi"
    if atlas_id == "ecoli_lgbm_single":
        return "LGBM single"
    if atlas_id == "weis_lr_ecoli_temporal_a2018":
        return "Weis LR"
    return ""


def value_cell(row: dict[str, str] | None) -> str:
    if row is None:
        return "--"
    suffix = "*" if "caution" in row.get("audit_adequacy_label", "") else ""
    return f"{fmt_auc(row.get('raw_auc', ''))}->{fmt_auc(row.get('stratum_centered_auc', ''))}{suffix}"


def split_value(row: dict[str, str] | None, column: str) -> str:
    if row is None:
        return "--"
    suffix = "*" if column == "stratum_centered_auc" and "caution" in row.get("audit_adequacy_label", "") else ""
    return f"{fmt_auc(row.get(column, ''))}{suffix}"


def consensus_for(drug: str, site: str) -> str:
    if drug == "Ciprofloxacin":
        return {
            "A-2018": "retained",
            "DRIAMS-C": "retained",
            "DRIAMS-D": "partial",
        }.get(site, "partial")
    if drug == "Amoxicillin-Clavulanic acid":
        return {
            "A-2018": "weak",
            "DRIAMS-C": "collapsed",
            "DRIAMS-D": "collapsed",
        }.get(site, "weak/collapsed")
    return ""


def build_model_class_table(atlas: list[dict[str, str]]) -> list[dict[str, str]]:
    target_sites = ["A-2018", "DRIAMS-C", "DRIAMS-D"]
    target_drugs = ["Ciprofloxacin", "Amoxicillin-Clavulanic acid"]
    model_order = ["CNN/Mega", "LGBM multi", "LGBM single", "Weis LR"]
    lookup: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in atlas:
        label = model_label(row.get("atlas_id", ""), row.get("model_name", ""))
        if row.get("atlas_id") == "weis_lr_ecoli":
            label = "Weis LR"
        if label not in model_order:
            continue
        site = row.get("site", "")
        if row.get("atlas_id") == "weis_lr_ecoli_temporal_a2018" and site == "DRIAMS-A":
            site = "A-2018"
        key = (site, row.get("drug", ""), label)
        lookup[key] = row

    output = []
    for drug in target_drugs:
        for site in target_sites:
            out = {
                "Site": site,
                "Drug": short_drug(drug),
                "CNN/Mega": value_cell(lookup.get((site, drug, "CNN/Mega"))),
                "LGBM multi": value_cell(lookup.get((site, drug, "LGBM multi"))),
                "LGBM single": value_cell(lookup.get((site, drug, "LGBM single"))),
                "Weis LR": value_cell(lookup.get((site, drug, "Weis LR"))),
            }
            output.append(out)
    return output


def build_model_class_split_table(atlas: list[dict[str, str]]) -> list[dict[str, str]]:
    target_sites = ["A-2018", "DRIAMS-C", "DRIAMS-D"]
    target_drugs = ["Ciprofloxacin", "Amoxicillin-Clavulanic acid"]
    model_order = ["CNN/Mega", "LGBM multi", "LGBM single", "Weis LR"]
    lookup: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in atlas:
        label = model_label(row.get("atlas_id", ""), row.get("model_name", ""))
        if row.get("atlas_id") == "weis_lr_ecoli":
            label = "Weis LR"
        if label not in model_order:
            continue
        site = row.get("site", "")
        if row.get("atlas_id") == "weis_lr_ecoli_temporal_a2018" and site == "DRIAMS-A":
            site = "A-2018"
        lookup[(site, row.get("drug", ""), label)] = row

    output = []
    for drug in target_drugs:
        for site in target_sites:
            row = {"Site": site, "Drug": short_drug(drug)}
            for model in model_order:
                item = lookup.get((site, drug, model))
                row[f"{model} raw"] = split_value(item, "raw_auc")
                row[f"{model} centered"] = split_value(item, "stratum_centered_auc")
            row["Consensus"] = consensus_for(drug, site)
            output.append(row)
    return output


def build_klebsiella_table(atlas: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [
        row
        for row in atlas
        if row.get("atlas_id") in {"klebsiella_cnn", "weis_lr_klebsiella_ceftriaxone"}
    ]
    rows.sort(key=lambda row: (row.get("atlas_id", ""), row.get("drug", ""), row.get("site", "")))
    output = []
    for row in rows:
        output.append(
            {
                "Model": row.get("model_name", ""),
                "Site": row.get("site", ""),
                "Drug": short_drug(row.get("drug", "")),
                "Raw AUC": fmt_auc(row.get("raw_auc", "")),
                "Centered AUC": fmt_auc(row.get("stratum_centered_auc", "")),
                "Background AUC": "--"
                if row.get("atlas_id") == "weis_lr_klebsiella_ceftriaxone"
                else fmt_auc(row.get("exact_background_auc", "")),
                "Retention (%)": fmt_pct(row.get("matched_retention", "")),
                "Strata": row.get("n_valid_strata", ""),
                "Interpretation": row.get("atlas_interpretation", "").replace("_", " "),
            }
        )
    return output


def organism_pair_label(organism: str, drug: str) -> str:
    if organism == "Staphylococcus aureus":
        return r"S. aureus / Oxacillin"
    if organism == "Klebsiella pneumoniae":
        return r"K. pneumoniae / Ceftriaxone"
    return f"{organism} / {short_drug(drug)}"


def build_organism_extension_table(atlas: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in atlas:
        include = (
            row.get("atlas_id") == "saureus_cnn"
            and row.get("drug") == "Oxacillin"
        ) or (
            row.get("atlas_id") == "klebsiella_cnn"
            and row.get("drug") == "Ceftriaxone"
        ) or row.get("atlas_id") == "weis_lr_klebsiella_ceftriaxone"
        if not include:
            continue
        rows.append(
            {
                "Organism/drug": organism_pair_label(row.get("organism", ""), row.get("drug", "")),
                "Model": row.get("model_name", ""),
                "Site": row.get("site", ""),
                "Raw AUC": fmt_auc(row.get("raw_auc", "")),
                "Centered AUC": fmt_auc(row.get("stratum_centered_auc", "")),
                "Background AUC": "--"
                if row.get("model_name", "") == "Weis LR"
                else fmt_auc(row.get("exact_background_auc", "")),
                "Retention (%)": fmt_pct(row.get("matched_retention", "")),
                "Strata": row.get("n_valid_strata", ""),
                "Interpretation": row.get("atlas_interpretation", "").replace("_", " "),
            }
        )

    if WEIS_SAUREUS_OXA.exists():
        for row in read_csv(WEIS_SAUREUS_OXA):
            rows.append(
                {
                    "Organism/drug": r"S. aureus / Oxacillin",
                    "Model": "Weis LR",
                    "Site": row.get("site", ""),
                    "Raw AUC": fmt_auc(row.get("raw_auc", "")),
                    "Centered AUC": fmt_auc(row.get("stratum_centered_auc", "")),
                    "Background AUC": "--",
                    "Retention (%)": fmt_pct(row.get("matched_retention", "")),
                    "Strata": row.get("n_valid_strata", ""),
                    "Interpretation": row.get("interpretation_category", "").replace("_", " "),
                }
            )

    site_rank = {"A-2018": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3}
    model_rank = {"CNN/Mega": 0, "Weis LR": 1}
    pair_rank = {
        r"S. aureus / Oxacillin": 0,
        r"K. pneumoniae / Ceftriaxone": 1,
    }
    return sorted(
        rows,
        key=lambda row: (
            pair_rank.get(row["Organism/drug"], 99),
            model_rank.get(row["Model"], 99),
            site_rank.get(row["Site"], 99),
        ),
    )


def build_atlas_table(summary: list[dict[str, str]]) -> list[dict[str, str]]:
    order = [
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
    rank = {name: index for index, name in enumerate(order)}
    rows = sorted(summary, key=lambda row: rank.get(row.get("atlas_id", ""), len(rank)))
    output = []
    for row in rows:
        output.append(
            {
                "Atlas set": row.get("atlas_id", ""),
                "Model": row.get("model_name", ""),
                "Scope": row.get("scope", ""),
                "Rows": row.get("n_matrix_rows", ""),
                "Drugs": row.get("n_drugs", ""),
                "Sites": row.get("n_sites", ""),
                "Interpretable": row.get("n_interpretable", ""),
                "Retained/partial": row.get("n_retained_or_partial", ""),
                "Background-sensitive": row.get("n_background_sensitive", ""),
            }
        )
    return output


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    atlas = read_csv(ATLAS_DIR / "atlas_matrix.csv")
    summary = read_csv(ATLAS_DIR / "atlas_summary.csv")

    model_rows = build_model_class_table(atlas)
    model_fields = ["Site", "Drug", "CNN/Mega", "LGBM multi", "LGBM single", "Weis LR"]
    write_csv(SOURCE_DIR / "source_data_table3_model_class_with_weis.csv", model_rows, model_fields)
    latex_table(
        TABLE_DIR / "table_2_model_replication_compact.tex",
        "Compact model-family replication of raw-to-background-centered attenuation.",
        "tab:model-replication-compact",
        model_fields,
        [[row[field] for field in model_fields] for row in model_rows],
        "llcccc",
        "Cells show raw AUC -> background-centered AUC. Asterisks mark caution rows with low matched support. Weis LR cells combine the temporal A-2018 export with existing Weis-style DRIAMS-C/DRIAMS-D compatibility rows; they are compatibility evidence rather than pooled benchmarks with the CNN/LightGBM panel.",
    )
    latex_table(
        TABLE_DIR / "table_3_model_class_with_weis.tex",
        "Model-class audit comparison with Weis LR compatibility rows.",
        "tab:model-class-with-weis",
        model_fields,
        [[row[field] for field in model_fields] for row in model_rows],
        "llcccc",
        "Cells show raw AUC -> background-centered AUC. Asterisks mark caution rows with low matched support. Weis LR cells combine the temporal A-2018 export with existing Weis-style DRIAMS-C/DRIAMS-D compatibility rows; they are not pooled benchmarks with the CNN/LightGBM panel.",
    )

    split_rows = build_model_class_split_table(atlas)
    split_fields = [
        "Site",
        "Drug",
        "CNN/Mega raw",
        "CNN/Mega centered",
        "LGBM multi raw",
        "LGBM multi centered",
        "LGBM single raw",
        "LGBM single centered",
        "Weis LR raw",
        "Weis LR centered",
        "Consensus",
    ]
    write_csv(SOURCE_DIR / "source_data_table3_model_family_split_columns.csv", split_rows, split_fields)
    latex_table(
        TABLE_DIR / "table_2_model_replication.tex",
        "Raw and background-centered AUCs across model families for the primary E. coli contrast.",
        "tab:model-replication",
        split_fields,
        [[row[field] for field in split_fields] for row in split_rows],
        "llrrrrrrrrl",
        "Columns separate raw and background-centered AUCs for direct model-family comparison. Asterisks mark caution rows with low matched support. Weis LR cells combine the temporal A-2018 export with existing Weis-style DRIAMS-C/DRIAMS-D compatibility rows; they are compatibility evidence rather than pooled benchmarks with the CNN/LightGBM panel.",
        font_size=r"\scriptsize",
    )

    extension_rows = build_organism_extension_table(atlas)
    extension_fields = [
        "Organism/drug",
        "Model",
        "Site",
        "Raw AUC",
        "Centered AUC",
        "Background AUC",
        "Retention (%)",
        "Strata",
        "Interpretation",
    ]
    write_csv(SOURCE_DIR / "source_data_table_organism_extension_representative.csv", extension_rows, extension_fields)
    latex_table(
        TABLE_DIR / "table_organism_extension_representative.tex",
        "Representative organism-extension audit rows.",
        "tab:organism-extension",
        extension_fields,
        [[row[field] for field in extension_fields] for row in extension_rows],
        "lllrrrrll",
        "Representative rows test whether the audit is only an E. coli/ciprofloxacin story. S. aureus/oxacillin is a second-organism retained-signal case; K. pneumoniae/ceftriaxone is a third-organism background-sensitive or sparse-control case. Weis LR rows are compatibility exports; background-only AUC was not computed for those rows.",
        font_size=r"\scriptsize",
    )

    kleb_rows = build_klebsiella_table(atlas)
    kleb_fields = [
        "Model",
        "Site",
        "Drug",
        "Raw AUC",
        "Centered AUC",
        "Background AUC",
        "Retention (%)",
        "Strata",
        "Interpretation",
    ]
    write_csv(SOURCE_DIR / "source_data_table_klebsiella_panel_audit.csv", kleb_rows, kleb_fields)
    latex_table(
        TABLE_DIR / "table_klebsiella_panel_audit.tex",
        "Full K. pneumoniae six-drug audit results and Weis ceftriaxone extension.",
        "tab:klebsiella-panel-audit",
        kleb_fields,
        [[row[field] for field in kleb_fields] for row in kleb_rows],
        "lllrrrrll",
        "The CNN/Mega rows summarize the K. pneumoniae six-drug panel. The Weis LR ceftriaxone rows are reported separately because they come from the official-panel compatibility export; background-only AUC was not computed for those rows.",
        font_size=r"\scriptsize",
    )

    atlas_rows = build_atlas_table(summary)
    atlas_fields = [
        "Atlas set",
        "Model",
        "Scope",
        "Rows",
        "Drugs",
        "Sites",
        "Interpretable",
        "Retained/partial",
        "Background-sensitive",
    ]
    write_csv(SOURCE_DIR / "source_data_table_audit_atlas_summary.csv", atlas_rows, atlas_fields)
    latex_table(
        TABLE_DIR / "table_audit_atlas_summary.tex",
        "Audit atlas coverage and outcome summary.",
        "tab:audit-atlas-summary",
        atlas_fields,
        [[row[field] for field in atlas_fields] for row in atlas_rows],
        "lllrrrrrr",
        "The atlas combines background-matched audit rows, exact-background baselines, calibration checks and sensitivity summaries. Retained/partial counts include rows classified as retained or partial residual signal by the atlas rules.",
        font_size=r"\scriptsize",
    )

    manifest_path = SOURCE_DIR / "source_data_manifest.csv"
    manifest_fields = ["display_item", "file", "description"]
    manifest = read_csv(manifest_path) if manifest_path.exists() else []
    additions = [
        {
            "display_item": "Table 3",
            "file": "source_data_table3_model_class_with_weis.csv",
            "description": "Model-class audit cells comparing CNN/Mega, LightGBM multi-task, LightGBM single-task and Weis LR compatibility rows.",
        },
        {
            "display_item": "Klebsiella table",
            "file": "source_data_table_klebsiella_panel_audit.csv",
            "description": "Klebsiella pneumoniae CNN/Mega six-drug panel rows plus Weis LR ceftriaxone compatibility rows.",
        },
        {
            "display_item": "Atlas table",
            "file": "source_data_table_audit_atlas_summary.csv",
            "description": "Audit atlas coverage and outcome counts by prediction set.",
        },
        {
            "display_item": "Figure 8",
            "file": "source_data_fig8_audit_atlas_summary.csv",
            "description": "Audit atlas coverage and outcome counts used for the atlas summary graph.",
        },
    ]
    existing = {(row.get("display_item", ""), row.get("file", "")) for row in manifest}
    manifest.extend(row for row in additions if (row["display_item"], row["file"]) not in existing)
    write_csv(manifest_path, manifest, manifest_fields)

    print("Wrote manuscript atlas/model tables.")


if __name__ == "__main__":
    main()
