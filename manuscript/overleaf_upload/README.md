# Manuscript Package

This folder contains an Overleaf-ready manuscript and vector PDF figures.

## Main files

- `main.tex` - Nature Communications-style manuscript with embedded references.
- `references.bib` - BibTeX convenience file for citation managers.
- `figures/*.pdf` - vector figures generated from repository CSV outputs.
- `tables/*.tex` - `booktabs` LaTeX tables included by `main.tex`.
- `source_data/*.csv` - source-data files for plotted figure panels.
- `supplementary/supplementary_information.tex` - supplementary notes and tables.

## Regenerate figures and tables

From the repository root:

```bash
python scripts/make_figures_mpl.py
```

The script reads committed source data and analysis outputs, then writes the
current vector PDFs to `manuscript/figures/`.

The public UPEC clone-control support analysis is regenerated separately:

```bash
python scripts/upec_clone_control_analysis.py
```

It writes WGS-lineage-controlled source tables to
`outputs/analysis_outputs/upec_clone_control_outputs/`, with manuscript copies
stored in `source_data/`.

## Overleaf

Upload the full `manuscript/` folder into Overleaf, preserving the `main.tex`,
`figures/`, `tables/`, `source_data/`, and `supplementary/` paths. The manuscript
uses standard LaTeX packages (`graphicx`, `booktabs`, `natbib`, `lineno`,
`hyperref`) and compiles without a custom Nature class file.
