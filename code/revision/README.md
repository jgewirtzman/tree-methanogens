# Revision analyses — NPH-MS-2026-56441

Self-contained set of scripts for the manuscript revision. Integrated into the main
repo (no separate top-level `revision/` folder): **code** lives here under `code/revision/`
and is git-tracked; **data** lives in the one consolidated `data/` tree (Zenodo drop-in);
**outputs** are written to `outputs/revision/`.

## Reproduce everything
```
Rscript code/revision/run_all.R
```
Runs every final generator (stats then figures), then assembles the numbered
manuscript figure set into `outputs/revision/figures/{main,SI}/`. Reads `data/` and
`code/`; writes only to `outputs/revision/`.

## Layout
```
code/revision/
├── run_all.R                  # one-command reproduce
├── rev_00_assemble_figures.R  # copies finals -> outputs/revision/figures/{main,SI}/ (numbered)
├── rev_prep_species_data.R    # shared, side-effect-free species-level data prep (sourced)
├── rev_fig01_*.R … rev_fig09_*.R, rev_figS02/04/11/12/15/17_*.R   # final figure generators
│                              #   (named by manuscript figure number; see MANIFEST)
├── rev_fig02a/07b/07c_*.R     # supporting panel scripts
├── rev_stat_*.R               # final analyses / tables (arcsinh redo, RMA, multigene, isotopes, …)
├── exploratory/               # superseded / provenance scripts (NOT run by run_all)
└── notes/                     # response-to-referees, captions, findings, planning docs
```

## Data
All inputs come from the consolidated `data/` archive (drop-in from Zenodo). Revision-
specific curated inputs were added there, alongside the existing files:
- `data/processed/molecular/methanotroph_definitions_revised.csv` (Known→Putative reclass)
- `data/processed/molecular/black_oak/` (extraction/core/soil masses for the felled-oak copies/g)

## Outputs
- `outputs/revision/` — all generated figures, reports, tables (CSV/TXT/PNG).
- `outputs/revision/exploratory/` — outputs of the exploratory scripts.
- `outputs/revision/figures/main/`, `.../SI/` — the numbered manuscript figure set.
- `outputs/revision/figures/MANIFEST.md` — number → content → source script → source PNG.

## Notes
`notes/` holds the prose deliverables and planning records: `response_to_referees_skeleton.md`,
`plant_traits_SI_and_discussion.md`, the isotope/microbiology findings, `SI_FIGURES_PLAN.md`,
`PENDING_CHANGES.md`, `STATUS.md`, and the verbatim reviews.
