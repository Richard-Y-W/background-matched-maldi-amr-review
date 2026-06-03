#!/usr/bin/env python3
"""Run the Weis LR E. coli six-drug temporal A-2018 export.

This wrapper produces the missing temporal row needed to include Weis cleanly
in the model-class matrix:

    train: DRIAMS-A 2015, 2016, 2017
    test : DRIAMS-A 2018

It delegates to ``export_weis_predictions_for_audit.py`` so the prediction CSV
and background-matched audit outputs have the same schema as the external-site
Weis exports.
"""

from __future__ import annotations

import argparse
import pathlib
import shlex
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
EXPORTER = ROOT / "scripts" / "export_weis_predictions_for_audit.py"
ECOLI6_DRUGS = (
    "Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,"
    "Ceftriaxone,Ceftazidime,Cefepime"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weis-repo", type=pathlib.Path, required=True)
    parser.add_argument("--driams-root", type=pathlib.Path, required=True)
    parser.add_argument(
        "--audit-script",
        type=pathlib.Path,
        default=ROOT / "run_background_audit_framework.py",
        help="Path to the background-matched audit framework script.",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=pathlib.Path("/kaggle/working/weis_lr_ecoli6_temporal_a2018"),
    )
    parser.add_argument("--model", default="lr")
    parser.add_argument("--seed", type=int, default=35)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--bootstrap-n", type=int, default=500)
    parser.add_argument("--permutation-n", type=int, default=500)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def build_command(args: argparse.Namespace) -> list[str]:
    return [
        sys.executable,
        str(EXPORTER),
        "--weis-repo",
        str(args.weis_repo),
        "--driams-root",
        str(args.driams_root),
        "--audit-script",
        str(args.audit_script),
        "--output-dir",
        str(args.output_dir),
        "--panel",
        "custom",
        "--species",
        "Escherichia coli",
        "--drugs",
        ECOLI6_DRUGS,
        "--model",
        str(args.model),
        "--seed",
        str(args.seed),
        "--n-folds",
        str(args.n_folds),
        "--train-site",
        "DRIAMS-A",
        "--test-sites",
        "DRIAMS-A",
        "--train-years",
        "2015,2016,2017",
        "--test-years",
        "2018",
        "--train-row-policy",
        "all",
        "--external-row-policy",
        "all",
        "--bootstrap-n",
        str(args.bootstrap_n),
        "--permutation-n",
        str(args.permutation_n),
    ]


def main() -> None:
    args = build_parser().parse_args()
    command = build_command(args)
    print(" ".join(shlex.quote(part) for part in command), flush=True)
    if args.dry_run:
        return
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
