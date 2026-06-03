#!/usr/bin/env python3
"""Run a background-matched audit atlas over one or more prediction tables.

The atlas runner is intentionally model-agnostic. It consumes long-format
isolate-level prediction CSVs, runs the existing audit modules for each
prediction set, and writes one combined matrix for paper/reviewer reporting.

Example:
    python scripts/run_audit_atlas.py \
        --prediction-set ecoli_cnn=outputs/analysis_outputs/background_matched_predictions.csv \
        --prediction-set saureus_cnn=outputs/analysis_outputs/saureus_panel_oxa_background_audit/background_matched_predictions.csv \
        --output-dir outputs/analysis_outputs/audit_atlas

Manifest CSVs are also supported with columns:
    atlas_id,predictions_csv,model_name,scope,notes

Optional per-row column overrides are supported in the manifest:
    id_col,site_col,year_col,organism_col,drug_col,label_col,prob_col,background_signature_col
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
AUDIT_SCRIPT = ROOT / "run_background_audit_framework.py"
BASELINE_SCRIPT = ROOT / "scripts" / "co_resistance_only_baseline.py"
CALIBRATION_SCRIPT = ROOT / "scripts" / "calibration_analysis.py"
SENSITIVITY_SCRIPT = ROOT / "scripts" / "sensitivity_sweep.py"


AUDIT_SUMMARY = "background_matched_audit_summary.csv"
BASELINE_SUMMARY = "co_resistance_only_baseline.csv"
CALIBRATION_SUMMARY = "calibration_summary.csv"
SENSITIVITY_DETAIL = "sensitivity_detail.csv"


MATRIX_FIELDS = [
    "atlas_id",
    "model_name",
    "scope",
    "site",
    "organism",
    "drug",
    "raw_auc",
    "raw_auc_ci_low",
    "raw_auc_ci_high",
    "matched_auc",
    "matched_auc_ci_low",
    "matched_auc_ci_high",
    "stratum_centered_auc",
    "stratum_centered_auc_ci_low",
    "stratum_centered_auc_ci_high",
    "pairwise_accuracy",
    "pairwise_comparisons",
    "permutation_p",
    "matched_retention",
    "n_total",
    "n_r",
    "n_matched",
    "n_matched_r",
    "n_valid_strata",
    "audit_adequacy_label",
    "audit_interpretation_category",
    "exact_background_auc",
    "background_burden_auc",
    "observed_minus_exact_background_auc",
    "observed_minus_burden_auc",
    "co_resistance_adequacy_label",
    "co_resistance_interpretation",
    "calibration_auc",
    "brier",
    "expected_calibration_error",
    "balanced_accuracy",
    "calibration_label",
    "centered_auc_n3",
    "retention_n3",
    "valid_strata_n3",
    "centered_auc_n5",
    "retention_n5",
    "valid_strata_n5",
    "centered_auc_n10",
    "retention_n10",
    "valid_strata_n10",
    "atlas_interpretation",
    "pre_specified_outcome",
    "prediction_match",
    "clonal_conservation_score",
    "co_resistance_collinearity",
    "pre_specified_rationale",
    "locked_before_run",
    "notes",
]

SUMMARY_FIELDS = [
    "atlas_id",
    "model_name",
    "scope",
    "prediction_csv",
    "n_matrix_rows",
    "n_organisms",
    "n_drugs",
    "n_sites",
    "n_interpretable",
    "n_retained_or_partial",
    "n_background_sensitive",
    "n_not_interpretable",
    "n_pre_specified",
    "n_pre_specified_matched",
]

EXPECTATION_CHECK_FIELDS = [
    "atlas_id",
    "model_name",
    "scope",
    "site",
    "organism",
    "drug",
    "predicted_outcome",
    "observed_outcome",
    "prediction_match",
    "raw_auc",
    "stratum_centered_auc",
    "matched_retention",
    "n_valid_strata",
    "clonal_conservation_score",
    "co_resistance_collinearity",
    "rationale",
    "locked_before_run",
    "expectation_source",
]


@dataclass(frozen=True)
class PredictionSet:
    atlas_id: str
    predictions_csv: Path
    model_name: str = ""
    scope: str = ""
    notes: str = ""
    id_col: str = ""
    site_col: str = ""
    year_col: str = ""
    organism_col: str = ""
    drug_col: str = ""
    label_col: str = ""
    prob_col: str = ""
    background_signature_col: str = ""


@dataclass(frozen=True)
class AtlasConfig:
    output_dir: Path
    min_pos_per_stratum: int
    min_neg_per_stratum: int
    bootstrap_n: int
    permutation_n: int
    sensitivity_thresholds: str
    baseline_min_n: int
    baseline_min_pos: int
    baseline_min_neg: int
    calibration_bins: int
    calibration_threshold: float
    match_year: bool
    background_drugs: str
    background_signature_col: str
    id_col: str
    site_col: str
    year_col: str
    organism_col: str
    drug_col: str
    label_col: str
    prob_col: str
    dry_run: bool
    skip_existing: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the background-matched audit, co-resistance-only baseline, "
            "calibration analysis, and stratum-size sensitivity sweep across "
            "one or more long prediction CSVs."
        )
    )
    parser.add_argument(
        "--prediction-set",
        action="append",
        default=[],
        metavar="ID=CSV",
        help="Prediction table to include. Repeatable. Example: ecoli_cnn=predictions.csv",
    )
    parser.add_argument(
        "--manifest-csv",
        action="append",
        default=[],
        help=(
            "Optional CSV manifest with columns atlas_id,predictions_csv,model_name,scope,notes. "
            "Column mapping overrides are also supported per row."
        ),
    )
    parser.add_argument(
        "--expectation-csv",
        "--locked-expectations-csv",
        action="append",
        default=[],
        help=(
            "Optional locked expectation CSV. Recommended columns: organism,drug,"
            "predicted_outcome,clonal_conservation_score,co_resistance_collinearity,"
            "rationale,locked_before_run. Optional columns: atlas_id,model_name,site."
        ),
    )
    parser.add_argument("--output-dir", default="outputs/analysis_outputs/audit_atlas")
    parser.add_argument("--min-pos-per-stratum", type=int, default=3)
    parser.add_argument("--min-neg-per-stratum", type=int, default=3)
    parser.add_argument("--bootstrap-n", type=int, default=500)
    parser.add_argument("--permutation-n", type=int, default=500)
    parser.add_argument("--sensitivity-thresholds", default="3,5,10")
    parser.add_argument("--baseline-min-n", type=int, default=30)
    parser.add_argument("--baseline-min-pos", type=int, default=3)
    parser.add_argument("--baseline-min-neg", type=int, default=3)
    parser.add_argument("--calibration-bins", type=int, default=10)
    parser.add_argument("--calibration-threshold", type=float, default=0.5)
    parser.add_argument("--match-year", action="store_true")
    parser.add_argument("--background-drugs", default="")
    parser.add_argument(
        "--background-signature-col",
        default="background_signature",
        help=(
            "Use this precomputed background-signature column when present. "
            "Leave empty to always derive signatures from the long table."
        ),
    )
    parser.add_argument("--id-col", default="isolate_id")
    parser.add_argument("--site-col", default="site")
    parser.add_argument("--year-col", default="year")
    parser.add_argument("--organism-col", default="organism")
    parser.add_argument("--drug-col", default="drug")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--prob-col", default="prob")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them.")
    parser.add_argument("--skip-existing", action="store_true", help="Reuse existing per-set outputs when present.")
    return parser.parse_args()


def sanitize_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "prediction_set"


def resolve_path(path_text: str, base_dir: Path | None = None) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path
    candidates: list[Path] = []
    if base_dir is not None:
        candidates.append(base_dir / path)
    candidates.extend([ROOT / path, Path.cwd() / path, path])
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def portable_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_prediction_set(text: str) -> PredictionSet:
    if "=" not in text:
        raise ValueError(f"--prediction-set must be formatted as ID=CSV, got: {text}")
    atlas_id, csv_path = text.split("=", 1)
    return PredictionSet(atlas_id=sanitize_id(atlas_id), predictions_csv=resolve_path(csv_path))


def load_manifest(path: Path) -> list[PredictionSet]:
    path = path.expanduser()
    base_dir = path.parent
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    sets: list[PredictionSet] = []
    for idx, row in enumerate(rows, start=2):
        pred = (row.get("predictions_csv") or row.get("prediction_csv") or "").strip()
        if not pred:
            raise ValueError(f"{path}:{idx} is missing predictions_csv")
        atlas_id = row.get("atlas_id") or Path(pred).stem
        sets.append(
            PredictionSet(
                atlas_id=sanitize_id(atlas_id),
                predictions_csv=resolve_path(pred, base_dir=base_dir),
                model_name=(row.get("model_name") or "").strip(),
                scope=(row.get("scope") or "").strip(),
                notes=(row.get("notes") or "").strip(),
                id_col=(row.get("id_col") or "").strip(),
                site_col=(row.get("site_col") or "").strip(),
                year_col=(row.get("year_col") or "").strip(),
                organism_col=(row.get("organism_col") or "").strip(),
                drug_col=(row.get("drug_col") or "").strip(),
                label_col=(row.get("label_col") or "").strip(),
                prob_col=(row.get("prob_col") or "").strip(),
                background_signature_col=(row.get("background_signature_col") or "").strip(),
            )
        )
    return sets


def load_expectations(paths: Sequence[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        path = path.expanduser()
        with path.open(newline="") as f:
            for idx, row in enumerate(csv.DictReader(f), start=2):
                organism = (row.get("organism") or "").strip()
                drug = (row.get("drug") or "").strip()
                predicted = (row.get("predicted_outcome") or row.get("expected_outcome") or "").strip()
                if not organism or not drug or not predicted:
                    raise ValueError(
                        f"{path}:{idx} must include organism, drug, and predicted_outcome"
                    )
                rows.append(
                    {
                        "atlas_id": sanitize_id(row.get("atlas_id") or "") if row.get("atlas_id") else "",
                        "model_name": (row.get("model_name") or "").strip(),
                        "site": (row.get("site") or "").strip(),
                        "organism": organism,
                        "drug": drug,
                        "predicted_outcome": predicted,
                        "clonal_conservation_score": (row.get("clonal_conservation_score") or "").strip(),
                        "co_resistance_collinearity": (row.get("co_resistance_collinearity") or "").strip(),
                        "rationale": (row.get("rationale") or "").strip(),
                        "locked_before_run": (row.get("locked_before_run") or "").strip(),
                        "expectation_source": portable_path(path),
                    }
                )
    return rows


def csv_has_column(path: Path, column: str) -> bool:
    if not column:
        return False
    with path.open(newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return False
    return column in header


def unique_prediction_sets(sets: Sequence[PredictionSet]) -> list[PredictionSet]:
    seen: set[str] = set()
    output: list[PredictionSet] = []
    for item in sets:
        if item.atlas_id in seen:
            raise ValueError(f"Duplicate atlas_id: {item.atlas_id}")
        seen.add(item.atlas_id)
        output.append(item)
    return output


def command_str(command: Sequence[str]) -> str:
    return " ".join(shlex_quote(part) for part in command)


def shlex_quote(value: object) -> str:
    text = str(value)
    if re.fullmatch(r"[A-Za-z0-9_./:=,+-]+", text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def run_command(command: Sequence[str], *, expected: Path, dry_run: bool, skip_existing: bool) -> None:
    print(f"$ {command_str(command)}", flush=True)
    if dry_run:
        return
    if skip_existing and expected.exists():
        print(f"  using existing {expected}", flush=True)
        return
    subprocess.run(list(command), cwd=ROOT, check=True)


def column_name(prediction_set: PredictionSet, config: AtlasConfig, field: str) -> str:
    override = getattr(prediction_set, field)
    return override or getattr(config, field)


def audit_command(config: AtlasConfig, prediction_set: PredictionSet, out_dir: Path) -> list[str]:
    command = [
        sys.executable,
        str(AUDIT_SCRIPT),
        "--predictions-csv",
        str(prediction_set.predictions_csv),
        "--output-dir",
        str(out_dir),
        "--model-name",
        prediction_set.model_name or prediction_set.atlas_id,
        "--id-col",
        column_name(prediction_set, config, "id_col"),
        "--site-col",
        column_name(prediction_set, config, "site_col"),
        "--year-col",
        column_name(prediction_set, config, "year_col"),
        "--organism-col",
        column_name(prediction_set, config, "organism_col"),
        "--drug-col",
        column_name(prediction_set, config, "drug_col"),
        "--label-col",
        column_name(prediction_set, config, "label_col"),
        "--prob-col",
        column_name(prediction_set, config, "prob_col"),
        "--min-pos-per-stratum",
        str(config.min_pos_per_stratum),
        "--min-neg-per-stratum",
        str(config.min_neg_per_stratum),
        "--bootstrap-n",
        str(config.bootstrap_n),
        "--permutation-n",
        str(config.permutation_n),
    ]
    if config.background_drugs:
        command.extend(["--background-drugs", config.background_drugs])
    if config.match_year:
        command.append("--match-year")
    background_signature_col = column_name(prediction_set, config, "background_signature_col")
    if csv_has_column(prediction_set.predictions_csv, background_signature_col):
        command.extend(["--background-signature-col", background_signature_col])
    return command


def baseline_command(config: AtlasConfig, prediction_set: PredictionSet, out_dir: Path) -> list[str]:
    return [
        sys.executable,
        str(BASELINE_SCRIPT),
        "--predictions-csv",
        str(prediction_set.predictions_csv),
        "--output-dir",
        str(out_dir),
        "--min-n",
        str(config.baseline_min_n),
        "--min-pos",
        str(config.baseline_min_pos),
        "--min-neg",
        str(config.baseline_min_neg),
    ]


def calibration_command(config: AtlasConfig, prediction_set: PredictionSet, out_dir: Path) -> list[str]:
    return [
        sys.executable,
        str(CALIBRATION_SCRIPT),
        "--predictions-csv",
        str(prediction_set.predictions_csv),
        "--output-dir",
        str(out_dir),
        "--n-bins",
        str(config.calibration_bins),
        "--threshold",
        str(config.calibration_threshold),
    ]


def sensitivity_command(config: AtlasConfig, prediction_set: PredictionSet, out_dir: Path) -> list[str]:
    command = [
        sys.executable,
        str(SENSITIVITY_SCRIPT),
        "--predictions-csv",
        str(prediction_set.predictions_csv),
        "--output-dir",
        str(out_dir),
        "--thresholds",
        config.sensitivity_thresholds,
        "--id-col",
        column_name(prediction_set, config, "id_col"),
        "--site-col",
        column_name(prediction_set, config, "site_col"),
        "--year-col",
        column_name(prediction_set, config, "year_col"),
        "--organism-col",
        column_name(prediction_set, config, "organism_col"),
        "--drug-col",
        column_name(prediction_set, config, "drug_col"),
        "--label-col",
        column_name(prediction_set, config, "label_col"),
        "--prob-col",
        column_name(prediction_set, config, "prob_col"),
    ]
    if config.background_drugs:
        command.extend(["--background-drugs", config.background_drugs])
    if config.match_year:
        command.append("--match-year")
    return command


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def key_for(row: dict) -> tuple[str, str, str]:
    return (row.get("site", ""), row.get("organism", ""), row.get("drug", ""))


def index_rows(rows: Iterable[dict]) -> dict[tuple[str, str, str], dict]:
    return {key_for(row): row for row in rows}


def index_sensitivity(rows: Iterable[dict]) -> dict[tuple[str, str, str], dict[str, dict]]:
    indexed: dict[tuple[str, str, str], dict[str, dict]] = {}
    for row in rows:
        threshold = str(row.get("min_stratum", "")).strip()
        indexed.setdefault(key_for(row), {})[threshold] = row
    return indexed


def parse_float(value: object) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def parse_int(value: object) -> int | None:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def is_interpretable(row: dict) -> bool:
    adequacy = str(row.get("adequacy_label", "")).lower()
    if adequacy.startswith("not_interpretable"):
        return False
    valid = parse_int(row.get("n_valid_strata"))
    centered = parse_float(row.get("stratum_centered_auc"))
    return bool(valid and valid > 0 and centered is not None)


def classify_atlas_row(audit: dict, baseline: dict) -> str:
    if not is_interpretable(audit):
        return "not_interpretable_or_sparse"

    adequacy = str(audit.get("adequacy_label", "")).lower()
    caution_prefix = "cautionary_" if adequacy.startswith("caution") else ""
    raw = parse_float(audit.get("raw_auc"))
    centered = parse_float(audit.get("stratum_centered_auc"))
    retention = parse_float(audit.get("matched_retention"))
    exact = parse_float(baseline.get("exact_background_auc"))

    if centered is None:
        return "not_interpretable_or_sparse"
    if exact is not None and raw is not None and exact >= raw - 0.03:
        if centered < 0.60:
            return "background_only_competitive_and_weak_residual"
        return f"{caution_prefix}background_only_competitive_but_residual_signal"
    if centered >= 0.70:
        return f"{caution_prefix}strong_within_background_signal"
    if centered >= 0.60:
        return f"{caution_prefix}moderate_within_background_signal"
    if centered >= 0.55:
        return f"{caution_prefix}weak_or_partial_within_background_signal"
    if raw is not None and raw >= 0.65 and centered < 0.55:
        return "background_sensitive_collapse"
    if retention is not None and retention < 0.10:
        return "too_little_matched_support"
    return "near_chance_or_ambiguous"


def normalize_expected_outcome(value: str) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")
    if text in {"retain", "retained", "retention", "strong_retention", "strong_signal", "survive", "survives"}:
        return "retention"
    if text in {"partial", "partial_retention", "weak_retention", "weak", "borderline", "attenuated_retention"}:
        return "partial_retention"
    if text in {"collapse", "collapsed", "near_chance", "weak_collapse", "background_sensitive", "collapse_or_weak"}:
        return "collapse"
    if text in {"collapse_or_sparse", "collapse_or_uninterpretable", "collapse_or_not_interpretable"}:
        return "collapse_or_sparse"
    if text in {"sparse", "uninterpretable", "not_interpretable", "structurally_uninterpretable", "low_support"}:
        return "not_interpretable"
    return text or "unspecified"


def expectation_match(predicted_outcome: str, observed_outcome: str) -> str:
    predicted = normalize_expected_outcome(predicted_outcome)
    observed = str(observed_outcome or "")
    retained = {
        "strong_within_background_signal",
        "moderate_within_background_signal",
        "background_only_competitive_but_residual_signal",
        "cautionary_strong_within_background_signal",
        "cautionary_moderate_within_background_signal",
        "cautionary_background_only_competitive_but_residual_signal",
    }
    partial = {
        "weak_or_partial_within_background_signal",
        "background_only_competitive_but_residual_signal",
        "cautionary_weak_or_partial_within_background_signal",
        "cautionary_background_only_competitive_but_residual_signal",
    }
    collapsed = {
        "background_sensitive_collapse",
        "background_only_competitive_and_weak_residual",
        "near_chance_or_ambiguous",
    }
    sparse = {"not_interpretable_or_sparse", "too_little_matched_support"}

    if predicted == "retention":
        return "matched" if observed in retained or observed in partial else "mismatched"
    if predicted == "partial_retention":
        if observed in partial:
            return "matched"
        if observed in retained:
            return "stronger_than_predicted"
        return "mismatched"
    if predicted == "collapse":
        return "matched" if observed in collapsed else "mismatched"
    if predicted == "collapse_or_sparse":
        return "matched" if observed in collapsed or observed in sparse else "mismatched"
    if predicted == "not_interpretable":
        return "matched" if observed in sparse else "mismatched"
    return "not_assessed"


def expectation_for_row(row: dict, expectations: Sequence[dict]) -> dict | None:
    candidates = []
    for expectation in expectations:
        if expectation.get("organism") != row.get("organism"):
            continue
        if expectation.get("drug") != row.get("drug"):
            continue
        score = 0
        for field in ("atlas_id", "model_name", "site"):
            expected_value = expectation.get(field, "")
            if expected_value:
                if expected_value != row.get(field):
                    score = -1
                    break
                score += 1
        if score >= 0:
            candidates.append((score, expectation))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def attach_expectations(matrix_rows: Sequence[dict], expectations: Sequence[dict]) -> list[dict]:
    output = []
    for row in matrix_rows:
        enriched = dict(row)
        expectation = expectation_for_row(enriched, expectations)
        if expectation:
            observed = enriched.get("atlas_interpretation", "")
            match = expectation_match(expectation.get("predicted_outcome", ""), observed)
            enriched.update(
                {
                    "pre_specified_outcome": expectation.get("predicted_outcome", ""),
                    "prediction_match": match,
                    "clonal_conservation_score": expectation.get("clonal_conservation_score", ""),
                    "co_resistance_collinearity": expectation.get("co_resistance_collinearity", ""),
                    "pre_specified_rationale": expectation.get("rationale", ""),
                    "locked_before_run": expectation.get("locked_before_run", ""),
                }
            )
        else:
            enriched.setdefault("pre_specified_outcome", "")
            enriched.setdefault("prediction_match", "")
            enriched.setdefault("clonal_conservation_score", "")
            enriched.setdefault("co_resistance_collinearity", "")
            enriched.setdefault("pre_specified_rationale", "")
            enriched.setdefault("locked_before_run", "")
        output.append(enriched)
    return output


def build_expectation_check_rows(matrix_rows: Sequence[dict], expectations: Sequence[dict]) -> list[dict]:
    rows = []
    for row in matrix_rows:
        if not row.get("pre_specified_outcome"):
            continue
        expectation = expectation_for_row(row, expectations) or {}
        rows.append(
            {
                "atlas_id": row.get("atlas_id", ""),
                "model_name": row.get("model_name", ""),
                "scope": row.get("scope", ""),
                "site": row.get("site", ""),
                "organism": row.get("organism", ""),
                "drug": row.get("drug", ""),
                "predicted_outcome": row.get("pre_specified_outcome", ""),
                "observed_outcome": row.get("atlas_interpretation", ""),
                "prediction_match": row.get("prediction_match", ""),
                "raw_auc": row.get("raw_auc", ""),
                "stratum_centered_auc": row.get("stratum_centered_auc", ""),
                "matched_retention": row.get("matched_retention", ""),
                "n_valid_strata": row.get("n_valid_strata", ""),
                "clonal_conservation_score": row.get("clonal_conservation_score", ""),
                "co_resistance_collinearity": row.get("co_resistance_collinearity", ""),
                "rationale": row.get("pre_specified_rationale", ""),
                "locked_before_run": row.get("locked_before_run", ""),
                "expectation_source": expectation.get("expectation_source", ""),
            }
        )
    return rows


def build_matrix_for_set(prediction_set: PredictionSet, run_dir: Path) -> list[dict]:
    audit_rows = read_csv(run_dir / "audit" / AUDIT_SUMMARY)
    baseline_rows = index_rows(read_csv(run_dir / "co_resistance_only" / BASELINE_SUMMARY))
    calibration_rows = index_rows(read_csv(run_dir / "calibration" / CALIBRATION_SUMMARY))
    sensitivity_rows = index_sensitivity(read_csv(run_dir / "sensitivity" / SENSITIVITY_DETAIL))

    matrix_rows: list[dict] = []
    for audit in audit_rows:
        key = key_for(audit)
        baseline = baseline_rows.get(key, {})
        calibration = calibration_rows.get(key, {})
        sensitivity = sensitivity_rows.get(key, {})
        row = {
            "atlas_id": prediction_set.atlas_id,
            "model_name": prediction_set.model_name or prediction_set.atlas_id,
            "scope": prediction_set.scope,
            "site": audit.get("site", ""),
            "organism": audit.get("organism", ""),
            "drug": audit.get("drug", ""),
            "raw_auc": audit.get("raw_auc", ""),
            "raw_auc_ci_low": audit.get("raw_auc_ci_low", ""),
            "raw_auc_ci_high": audit.get("raw_auc_ci_high", ""),
            "matched_auc": audit.get("matched_auc", ""),
            "matched_auc_ci_low": audit.get("matched_auc_ci_low", ""),
            "matched_auc_ci_high": audit.get("matched_auc_ci_high", ""),
            "stratum_centered_auc": audit.get("stratum_centered_auc", ""),
            "stratum_centered_auc_ci_low": audit.get("stratum_centered_auc_ci_low", ""),
            "stratum_centered_auc_ci_high": audit.get("stratum_centered_auc_ci_high", ""),
            "pairwise_accuracy": audit.get("pairwise_accuracy", ""),
            "pairwise_comparisons": audit.get("pairwise_comparisons", ""),
            "permutation_p": audit.get("permutation_p", ""),
            "matched_retention": audit.get("matched_retention", ""),
            "n_total": audit.get("n_total", ""),
            "n_r": audit.get("n_r", ""),
            "n_matched": audit.get("n_matched", ""),
            "n_matched_r": audit.get("n_matched_r", ""),
            "n_valid_strata": audit.get("n_valid_strata", ""),
            "audit_adequacy_label": audit.get("adequacy_label", ""),
            "audit_interpretation_category": audit.get("interpretation_category", ""),
            "exact_background_auc": baseline.get("exact_background_auc", ""),
            "background_burden_auc": baseline.get("background_burden_auc", ""),
            "observed_minus_exact_background_auc": baseline.get("observed_minus_exact_background_auc", ""),
            "observed_minus_burden_auc": baseline.get("observed_minus_burden_auc", ""),
            "co_resistance_adequacy_label": baseline.get("adequacy_label", ""),
            "co_resistance_interpretation": baseline.get("baseline_interpretation", ""),
            "calibration_auc": calibration.get("auc", ""),
            "brier": calibration.get("brier", ""),
            "expected_calibration_error": calibration.get("expected_calibration_error", ""),
            "balanced_accuracy": calibration.get("balanced_accuracy", ""),
            "calibration_label": calibration.get("calibration_label", ""),
            "notes": prediction_set.notes,
        }
        for threshold in ("3", "5", "10"):
            sens = sensitivity.get(threshold, {})
            row[f"centered_auc_n{threshold}"] = sens.get("centered_auc", "")
            row[f"retention_n{threshold}"] = sens.get("matched_retention", "")
            row[f"valid_strata_n{threshold}"] = sens.get("n_valid_strata", "")
        row["atlas_interpretation"] = classify_atlas_row(audit, baseline)
        matrix_rows.append(row)
    return matrix_rows


def summarize_set(prediction_set: PredictionSet, matrix_rows: Sequence[dict]) -> dict:
    interpretations = [row.get("atlas_interpretation", "") for row in matrix_rows]
    return {
        "atlas_id": prediction_set.atlas_id,
        "model_name": prediction_set.model_name or prediction_set.atlas_id,
        "scope": prediction_set.scope,
        "prediction_csv": portable_path(prediction_set.predictions_csv),
        "n_matrix_rows": len(matrix_rows),
        "n_organisms": len({row.get("organism", "") for row in matrix_rows if row.get("organism", "")}),
        "n_drugs": len({row.get("drug", "") for row in matrix_rows if row.get("drug", "")}),
        "n_sites": len({row.get("site", "") for row in matrix_rows if row.get("site", "")}),
        "n_interpretable": sum(1 for row in matrix_rows if row.get("atlas_interpretation") != "not_interpretable_or_sparse"),
        "n_retained_or_partial": sum(
            1
            for value in interpretations
            if value
            in {
                "strong_within_background_signal",
                "moderate_within_background_signal",
                "weak_or_partial_within_background_signal",
                "background_only_competitive_but_residual_signal",
                "cautionary_strong_within_background_signal",
                "cautionary_moderate_within_background_signal",
                "cautionary_weak_or_partial_within_background_signal",
                "cautionary_background_only_competitive_but_residual_signal",
            }
        ),
        "n_background_sensitive": sum(
            1
            for value in interpretations
            if value
            in {
                "background_sensitive_collapse",
                "background_only_competitive_and_weak_residual",
            }
        ),
        "n_not_interpretable": sum(1 for value in interpretations if value == "not_interpretable_or_sparse"),
        "n_pre_specified": sum(1 for row in matrix_rows if row.get("pre_specified_outcome")),
        "n_pre_specified_matched": sum(1 for row in matrix_rows if row.get("prediction_match") == "matched"),
    }


def write_csv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def markdown_table(rows: Sequence[dict], fields: Sequence[str], limit: int = 40) -> str:
    if not rows:
        return "_No rows._"
    display = list(rows[:limit])
    widths = {
        field: max(len(field), *(len(str(row.get(field, ""))) for row in display))
        for field in fields
    }
    lines = [
        "| " + " | ".join(field.ljust(widths[field]) for field in fields) + " |",
        "| " + " | ".join("-" * widths[field] for field in fields) + " |",
    ]
    for row in display:
        lines.append("| " + " | ".join(str(row.get(field, "")).ljust(widths[field]) for field in fields) + " |")
    if len(rows) > limit:
        lines.extend(["", f"_Showing first {limit} of {len(rows)} rows. See CSV for the complete atlas._"])
    return "\n".join(lines)


def write_markdown_summary(output_dir: Path, summary_rows: Sequence[dict], matrix_rows: Sequence[dict]) -> None:
    spotlight_fields = [
        "atlas_id",
        "site",
        "organism",
        "drug",
        "raw_auc",
        "stratum_centered_auc",
        "exact_background_auc",
        "matched_retention",
        "n_valid_strata",
        "atlas_interpretation",
        "pre_specified_outcome",
        "prediction_match",
    ]
    md = [
        "# Audit Atlas Summary",
        "",
        "This atlas combines background-matched audit outputs, co-resistance-only baselines, calibration metrics, and stratum-size sensitivity checks across prediction tables.",
        "",
        "## Prediction Sets",
        "",
        markdown_table(summary_rows, SUMMARY_FIELDS),
        "",
        "## Atlas Matrix Preview",
        "",
        markdown_table(matrix_rows, spotlight_fields),
        "",
        "Primary CSV outputs:",
        "",
        "- `atlas_matrix.csv`: one row per atlas/model/site/organism/drug result",
        "- `atlas_summary.csv`: one row per prediction set",
        "- `runs/<atlas_id>/`: per-set audit, baseline, calibration, and sensitivity outputs",
        "",
    ]
    (output_dir / "atlas_summary.md").write_text("\n".join(md))


def config_from_args(args: argparse.Namespace) -> AtlasConfig:
    return AtlasConfig(
        output_dir=Path(args.output_dir),
        min_pos_per_stratum=args.min_pos_per_stratum,
        min_neg_per_stratum=args.min_neg_per_stratum,
        bootstrap_n=args.bootstrap_n,
        permutation_n=args.permutation_n,
        sensitivity_thresholds=args.sensitivity_thresholds,
        baseline_min_n=args.baseline_min_n,
        baseline_min_pos=args.baseline_min_pos,
        baseline_min_neg=args.baseline_min_neg,
        calibration_bins=args.calibration_bins,
        calibration_threshold=args.calibration_threshold,
        match_year=args.match_year,
        background_drugs=args.background_drugs,
        background_signature_col=args.background_signature_col,
        id_col=args.id_col,
        site_col=args.site_col,
        year_col=args.year_col,
        organism_col=args.organism_col,
        drug_col=args.drug_col,
        label_col=args.label_col,
        prob_col=args.prob_col,
        dry_run=args.dry_run,
        skip_existing=args.skip_existing,
    )


def main() -> None:
    args = parse_args()
    prediction_sets = [parse_prediction_set(text) for text in args.prediction_set]
    for manifest in args.manifest_csv:
        prediction_sets.extend(load_manifest(Path(manifest)))
    prediction_sets = unique_prediction_sets(prediction_sets)
    if not prediction_sets:
        raise SystemExit("No prediction sets provided. Use --prediction-set ID=CSV or --manifest-csv.")
    expectations = load_expectations([Path(path) for path in args.expectation_csv])

    config = config_from_args(args)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    for prediction_set in prediction_sets:
        if not prediction_set.predictions_csv.exists():
            raise FileNotFoundError(f"Prediction CSV not found: {prediction_set.predictions_csv}")

    all_matrix_rows: list[dict] = []
    summary_rows: list[dict] = []

    for prediction_set in prediction_sets:
        run_dir = config.output_dir / "runs" / prediction_set.atlas_id
        audit_dir = run_dir / "audit"
        baseline_dir = run_dir / "co_resistance_only"
        calibration_dir = run_dir / "calibration"
        sensitivity_dir = run_dir / "sensitivity"

        print(f"\n=== Atlas set: {prediction_set.atlas_id} ===", flush=True)
        run_command(
            audit_command(config, prediction_set, audit_dir),
            expected=audit_dir / AUDIT_SUMMARY,
            dry_run=config.dry_run,
            skip_existing=config.skip_existing,
        )
        run_command(
            baseline_command(config, prediction_set, baseline_dir),
            expected=baseline_dir / BASELINE_SUMMARY,
            dry_run=config.dry_run,
            skip_existing=config.skip_existing,
        )
        run_command(
            calibration_command(config, prediction_set, calibration_dir),
            expected=calibration_dir / CALIBRATION_SUMMARY,
            dry_run=config.dry_run,
            skip_existing=config.skip_existing,
        )
        run_command(
            sensitivity_command(config, prediction_set, sensitivity_dir),
            expected=sensitivity_dir / SENSITIVITY_DETAIL,
            dry_run=config.dry_run,
            skip_existing=config.skip_existing,
        )

        if config.dry_run:
            continue

        matrix_rows = build_matrix_for_set(prediction_set, run_dir)
        if expectations:
            matrix_rows = attach_expectations(matrix_rows, expectations)
        all_matrix_rows.extend(matrix_rows)
        summary_rows.append(summarize_set(prediction_set, matrix_rows))

    if config.dry_run:
        print("\nDry run complete. No atlas files were written.", flush=True)
        return

    write_csv(config.output_dir / "atlas_matrix.csv", all_matrix_rows, MATRIX_FIELDS)
    write_csv(config.output_dir / "atlas_summary.csv", summary_rows, SUMMARY_FIELDS)
    if expectations:
        write_csv(
            config.output_dir / "pre_specified_prediction_check.csv",
            build_expectation_check_rows(all_matrix_rows, expectations),
            EXPECTATION_CHECK_FIELDS,
        )
    write_markdown_summary(config.output_dir, summary_rows, all_matrix_rows)

    print(f"\nWrote atlas matrix rows: {len(all_matrix_rows)}", flush=True)
    print(f"Atlas output: {config.output_dir}", flush=True)


if __name__ == "__main__":
    main()
