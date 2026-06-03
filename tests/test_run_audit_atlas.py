import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_audit_atlas.py"


def write_toy_predictions(path: pathlib.Path) -> None:
    rows = []
    isolate_idx = 0
    for drug_a_label in [0, 1]:
        for drug_b_label in [0, 1]:
            for rep in range(13):
                isolate_idx += 1
                isolate_id = f"toy_{isolate_idx}"
                rows.append(
                    {
                        "isolate_id": isolate_id,
                        "site": "A-2018",
                        "year": "2020",
                        "organism": "Toy species",
                        "drug": "Drug A",
                        "label": str(drug_a_label),
                        "prob": str(0.88 if drug_a_label else 0.12),
                    }
                )
                rows.append(
                    {
                        "isolate_id": isolate_id,
                        "site": "A-2018",
                        "year": "2020",
                        "organism": "Toy species",
                        "drug": "Drug B",
                        "label": str(drug_b_label),
                        "prob": str(0.82 if drug_b_label else 0.18),
                    }
                )

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_expectations(path: pathlib.Path) -> None:
    rows = [
        {
            "organism": "Toy species",
            "drug": "Drug A",
            "predicted_outcome": "retention",
            "clonal_conservation_score": "high",
            "co_resistance_collinearity": "moderate",
            "rationale": "synthetic focal signal should survive matching",
            "locked_before_run": "yes",
        },
        {
            "organism": "Toy species",
            "drug": "Drug B",
            "predicted_outcome": "retention",
            "clonal_conservation_score": "high",
            "co_resistance_collinearity": "moderate",
            "rationale": "synthetic focal signal should survive matching",
            "locked_before_run": "yes",
        },
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


class RunAuditAtlasTests(unittest.TestCase):
    def test_collapse_or_sparse_expectation_matches_collapse_and_sparse_outcomes(self):
        sys.path.insert(0, str(ROOT))
        from scripts.run_audit_atlas import expectation_match

        self.assertEqual(
            expectation_match("collapse_or_sparse", "background_sensitive_collapse"),
            "matched",
        )
        self.assertEqual(
            expectation_match("collapse_or_sparse", "not_interpretable_or_sparse"),
            "matched",
        )
        self.assertEqual(
            expectation_match("collapse_or_sparse", "strong_within_background_signal"),
            "mismatched",
        )

    def test_cli_writes_atlas_outputs_and_pre_specified_prediction_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = pathlib.Path(tmp)
            predictions = tmp_path / "toy_predictions.csv"
            expectations = tmp_path / "locked_expectations.csv"
            output_dir = tmp_path / "atlas"
            write_toy_predictions(predictions)
            write_expectations(expectations)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--prediction-set",
                    f"toy={predictions}",
                    "--expectation-csv",
                    str(expectations),
                    "--output-dir",
                    str(output_dir),
                    "--min-pos-per-stratum",
                    "3",
                    "--min-neg-per-stratum",
                    "3",
                    "--baseline-min-n",
                    "52",
                    "--baseline-min-pos",
                    "3",
                    "--baseline-min-neg",
                    "3",
                    "--bootstrap-n",
                    "0",
                    "--permutation-n",
                    "0",
                    "--sensitivity-thresholds",
                    "3,5,10",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "atlas_matrix.csv").exists())
            self.assertTrue((output_dir / "atlas_summary.csv").exists())
            self.assertTrue((output_dir / "pre_specified_prediction_check.csv").exists())

            with (output_dir / "pre_specified_prediction_check.csv").open(newline="") as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 2)
            self.assertEqual({row["prediction_match"] for row in rows}, {"matched"})
            self.assertEqual({row["locked_before_run"] for row in rows}, {"yes"})


if __name__ == "__main__":
    unittest.main()
