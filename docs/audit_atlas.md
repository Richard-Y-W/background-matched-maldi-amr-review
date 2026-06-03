# Audit Atlas Runner

`scripts/run_audit_atlas.py` runs the background-matched audit workflow across one or more
long-format prediction CSVs and combines the outputs into a single atlas matrix.

The runner does not train models. It consumes isolate-level prediction tables from any model
class, then runs:

- background-matched audit
- co-resistance-only baseline
- calibration analysis
- minimum-stratum sensitivity sweep
- optional pre-specified prediction checks

## Minimal Command

```bash
python scripts/run_audit_atlas.py \
  --prediction-set ecoli_cnn=outputs/analysis_outputs/background_matched_predictions.csv \
  --output-dir outputs/analysis_outputs/audit_atlas \
  --bootstrap-n 500 \
  --permutation-n 500 \
  --sensitivity-thresholds 3,5,10
```

Multiple prediction sets can be supplied by repeating `--prediction-set`:

```bash
python scripts/run_audit_atlas.py \
  --prediction-set ecoli_cnn=path/to/ecoli_cnn_predictions.csv \
  --prediction-set ecoli_lgbm_single=path/to/ecoli_lgbm_single_predictions.csv \
  --prediction-set saureus_cnn=path/to/saureus_cnn_predictions.csv \
  --output-dir outputs/analysis_outputs/audit_atlas
```

## Manifest Mode

For a larger atlas, use a manifest CSV:

```csv
atlas_id,predictions_csv,model_name,scope,notes
ecoli_cnn,path/to/ecoli_cnn_predictions.csv,CNN/Mega,E. coli six-drug panel,primary model
saureus_cnn,path/to/saureus_cnn_predictions.csv,CNN/Mega,S. aureus oxacillin panel,generality check
```

Template: [`docs/atlas_manifest_template.csv`](atlas_manifest_template.csv).

Run:

```bash
python scripts/run_audit_atlas.py \
  --manifest-csv atlas_manifest.csv \
  --output-dir outputs/analysis_outputs/audit_atlas
```

Relative `predictions_csv` paths are resolved relative to the manifest file.

## Locked Expectations

To make the atlas more than descriptive, pre-specify expected retain/collapse outcomes before
running new organism-drug pairs. The expectation CSV should be written before looking at the
new audit outputs.

```csv
organism,drug,predicted_outcome,clonal_conservation_score,co_resistance_collinearity,rationale,locked_before_run
Escherichia coli,Ciprofloxacin,partial_retention,high,high,fluoroquinolone resistance expected to retain lineage-linked residual signal,yes
Escherichia coli,Amoxicillin-Clavulanic acid,collapse,low,moderate,beta-lactam phenotype expected to be heterogeneous and background-sensitive,yes
Staphylococcus aureus,Oxacillin,retention,high,moderate,MRSA/MSSA status expected to retain within-background signal,yes
```

Template: [`docs/locked_atlas_expectations_template.csv`](locked_atlas_expectations_template.csv).

Run:

```bash
python scripts/run_audit_atlas.py \
  --manifest-csv atlas_manifest.csv \
  --expectation-csv locked_atlas_expectations.csv \
  --output-dir outputs/analysis_outputs/audit_atlas
```

The runner writes `pre_specified_prediction_check.csv`, with one row per observed atlas row that
matches a locked expectation. This table reports the predicted outcome, observed atlas outcome,
and whether they matched.

Supported `predicted_outcome` labels include:

- `retention`
- `partial_retention`
- `collapse`
- `collapse_or_sparse`
- `not_interpretable`

Use `collapse_or_sparse` when the biological expectation is that a row should either collapse
toward background-sensitive/near-chance behavior or become unsuitable for interpretation because
matched strata are too sparse. This is useful for atlas extensions where co-resistance biology and
sample support both make a strong retained-signal prediction unlikely.

Optional expectation columns `atlas_id`, `model_name`, and `site` can be used to make predictions
specific to a model or site. If these columns are absent, expectations match by organism and drug.

## Main Outputs

- `atlas_matrix.csv`: one row per atlas/model/site/organism/drug result
- `atlas_summary.csv`: one row per prediction set
- `atlas_summary.md`: compact markdown summary
- `pre_specified_prediction_check.csv`: locked prediction check, when `--expectation-csv` is used
- `runs/<atlas_id>/audit/`: per-set background-matched audit outputs
- `runs/<atlas_id>/co_resistance_only/`: no-spectrum co-resistance baseline
- `runs/<atlas_id>/calibration/`: calibration and threshold metrics
- `runs/<atlas_id>/sensitivity/`: minimum-stratum sensitivity sweep

## Fast Smoke Test

For a quick structural check, use zero bootstrap/permutation iterations:

```bash
python scripts/run_audit_atlas.py \
  --prediction-set example=example_predictions.csv \
  --output-dir /tmp/audit_atlas_smoke \
  --bootstrap-n 0 \
  --permutation-n 0 \
  --sensitivity-thresholds 3,5,10
```
