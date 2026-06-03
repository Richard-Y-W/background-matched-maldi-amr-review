import importlib.util
import pathlib
import sys
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_weis_ecoli6_temporal_export.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("run_weis_ecoli6_temporal_export", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class RunWeisEcoli6TemporalExportTests(unittest.TestCase):
    def test_command_targets_a2018_temporal_ecoli_six_drug_export(self):
        runner = load_runner()
        args = runner.build_parser().parse_args(
            [
                "--weis-repo",
                "/kaggle/working/maldi_amr",
                "--driams-root",
                "/kaggle/input/datasets/drscarlat/driams",
                "--output-dir",
                "/kaggle/working/weis_lr_ecoli6_temporal_a2018",
                "--dry-run",
            ]
        )

        command = runner.build_command(args)
        joined = " ".join(command)

        self.assertIn("--species Escherichia coli", joined)
        self.assertIn("--test-sites DRIAMS-A", joined)
        self.assertIn("--train-years 2015,2016,2017", joined)
        self.assertIn("--test-years 2018", joined)
        self.assertIn("--train-row-policy all", joined)
        self.assertIn("--external-row-policy all", joined)
        self.assertIn("--model lr", joined)
        self.assertIn("Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,Ceftriaxone,Ceftazidime,Cefepime", joined)


if __name__ == "__main__":
    unittest.main()
