# Revision analyses — NPH-MS-2026-56441

Self-contained set of scripts for the manuscript revision. Integrated into the main
repo (no separate top-level `revision/` folder): **code** lives here under `code/revision/`
and is git-tracked; **data** lives in the one consolidated `data/` tree (Zenodo drop-in);
**outputs** are written to `outputs/revision/`.

## Reproduce everything
```
Rscript code/run_all.R
```
Runs the original figure pipeline first (`generate_all_figures.R`, which populates the
**unchanged** figures like Fig 8 and several SI — two scripts fail there by design, fig5 &
S12, both replaced by revision versions), then every revision generator (stats then figures
via glob), then assembles the numbered set into `outputs/figures/{main,SI,photos}/`.
A full run yields **9 main + 23 SI + 1 photo, 0 missing**. Reads `data/` and `code/`; writes
to `outputs/figures/` (originals) and `outputs/revision/`.

**Master inventory of everything new/updated in the revision:** `notes/REVISION_INVENTORY.md`.

## Layout
```
code/revision/
├── run_all.R                  # one-command reproduce (globs rev_stat_*/rev_fig*)
├── 00_assemble_figures.R  # copies finals -> outputs/figures/{main,SI,photos}/ (numbered)
├── prep_species_data.R    # shared, side-effect-free species-level data prep (sourced)
├── rev_fig01_*.R … rev_fig09_*.R    # main-figure generators (by manuscript number)
│     rev_fig06_hydrogenotrophy      #   NEW synthesis main fig
│     rev_fig07_decay-methanogenesis #   NEW expanded Fig 7 (decay + fungi + felled-oak; replaces felled-oak-only)
├── rev_figS02/04/11/12/15/17/19/20/21_*.R   # SI-figure generators
│     rev_figS19 probe validation · rev_figS20 stem deterioration · rev_figS21 ddPCR-16S concordance  # NEW
├── figS_black-oak-cross-sections.R      # NEW felled-oak cross-section photo plate (-> photos/)
├── rev_fig02a/07b/07c_*.R     # supporting panel scripts
├── rev_stat_*.R               # final analyses / tables (campaign counts, RMA, multigene, isotopes, …)
├── bo-funguild-saprotroph.R  # exploratory FUNGuild on felled-oak ITS (not wired into a figure)
├── exploratory/               # superseded / provenance scripts (NOT run by run_all)
└── notes/                     # REVISION_INVENTORY.md, response-to-referees, captions, findings, planning
```

## Data
All inputs come from the consolidated `data/` archive (drop-in from Zenodo). Revision-
specific curated inputs were added there, alongside the existing files:
- `data/processed/molecular/methanotroph_definitions_revised.csv` (Known→Putative reclass)
- `data/processed/molecular/black_oak/` — extraction/core/soil masses (felled-oak copies/g) and
  `bo_its_load.csv` (felled-oak QUVE ITS load for Fig 7 panel f; makes Fig 7 standalone)

## Outputs
- `outputs/revision/` — all generated figures, reports, tables (CSV/TXT/PNG).
- `outputs/scratch/` — outputs of the exploratory scripts (incl. `explore_*.png`).
- `outputs/figures/original/main/`, `.../SI/`, `.../photos/` — the numbered manuscript figure set
  (photos = field plates, e.g. cross-sections + chamber photos; separate from SI data figures).
- Some assembled figures are original-pipeline figures (unchanged); run `generate_all_figures.R`
  first to populate them, else the assembler reports them MISSING.

## Notes
`notes/` holds the prose deliverables and planning records: **`REVISION_INVENTORY.md`** (master
stock-take), `response_to_referees_skeleton.md`, `plant_traits_SI_and_discussion.md`, the
isotope/microbiology findings, `SI_FIGURES_PLAN.md`, `PENDING_CHANGES.md`, `STATUS.md`, and the
verbatim reviews.
