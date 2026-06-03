#!/usr/bin/env python3
"""Build manuscript figures and tables.

This wrapper regenerates both the legacy final-framework artifacts under
outputs/final_framework_outputs and the current manuscript vector PDFs under
manuscript/figures.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILDER = ROOT / "scripts" / "make_final_framework_tables_figures.py"
MANUSCRIPT_FIGURE_BUILDER = ROOT / "scripts" / "make_figures_mpl.py"
ATLAS_TABLE_BUILDER = ROOT / "scripts" / "make_audit_atlas_tables.py"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build final paper tables and figures.")
    p.add_argument("--builder", type=Path, default=DEFAULT_BUILDER)
    p.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory forwarded to the final artifact builder.",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    cmd = [sys.executable, str(args.builder)]
    if args.output_dir is not None:
        cmd.extend(["--output-dir", str(args.output_dir)])
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    atlas_table_cmd = [sys.executable, str(ATLAS_TABLE_BUILDER)]
    print("Running:", " ".join(atlas_table_cmd), flush=True)
    subprocess.run(atlas_table_cmd, check=True)

    manuscript_cmd = [sys.executable, str(MANUSCRIPT_FIGURE_BUILDER)]
    print("Running:", " ".join(manuscript_cmd), flush=True)
    subprocess.run(manuscript_cmd, check=True)


if __name__ == "__main__":
    main()
