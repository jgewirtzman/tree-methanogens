# Tree Microbiomes and Methane Emissions in Upland Forests

[![DOI](https://zenodo.org/badge/1058191477.svg)](https://doi.org/10.5281/zenodo.18779715)

Code for a study of CH₄ exchange by tree stems and soils at Yale-Myers Forest, Connecticut, and its relation to the methanogens and methanotrophs in wood and soil:

> Gewirtzman J, Arnold W, et al. *Methane emissions from upland tree stems reflect a balance of microbial production and consumption.* Preprint: [bioRxiv 10.1101/2025.09.30.679632](https://www.biorxiv.org/content/10.1101/2025.09.30.679632v2)

The data: 1,200 stem and 288 soil chamber deployments (2020–2023, goFlux fits of analyzer traces); ddPCR of *mcrA*, *pmoA* and *mmoX* in wood and soil; 16S amplicon sequencing with PICRUSt2 and FAPROTAX; internal stem gas and its isotopes; and random forests that scale stem and soil flux to the 8,194 stems of a 3.72-ha censused stand.

## Getting started

**Data.** Zenodo record [10.5281/zenodo.18779715](https://doi.org/10.5281/zenodo.18779715) (concept DOI: always the latest version). Unzip each file at the repository root:

| File | Contents |
|---|---|
| `tree-methanogens_data.zip` | `data/raw/` (field, instrument and lab records; `README_raw.md`, hand corrections in `CORRECTIONS.md`), `data/tables/` (the canonical tables), `data/processed/` (stage-A products), `data/compiled/` (`README_data.md`, `results/`) |
| `tree-methanogens_picrust.zip` | optional: `data/raw/picrust/`, the PICRUSt2 predictions (7.9 GB unzipped) |
| `tree-methanogens_code.zip` | this repository at the archived commit |

The PICRUSt2 download is needed only to re-derive `data/processed/molecular/picrust/picrust_pathways.csv` and `picrust_pathway_info.csv`, which are in the data zip. Without it, a run from raw uses those two files, checks them against checksums recorded in `code/01_import/06_picrust_pathways.R` and says so in `outputs/logs/06_picrust_pathways.txt`. The 16S reads are in the NCBI SRA (BioProject [PRJNA1124946](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1124946)); the pipeline starts from the ASV table. Git tracks only the folder layout of `data/` and `outputs/`.

**Environment.** R 4.4.3; every package is pinned in `renv.lock` (including goFlux 0.2.0 from GitHub). Start R at the repository root and run `renv::restore()` first. Run under a UTF-8 locale (some column names contain µ and °); if `locale` shows `LC_CTYPE="C"`, prefix commands with `LC_ALL=en_US.UTF-8`.

**Running.** From the repository root:

```bash
Rscript code/run_all.R                    # full run, stages B–G, from data/processed/ and data/raw/
Rscript code/run_all.R --from raw         # full run, stages A–G: refit every flux from the analyzer traces
Rscript code/run_all.R --from C           # reuse the fitted model in outputs/models/
Rscript code/run_all.R --only E           # one stage (figures; also Rscript code/make_figures.R)
Rscript code/run_all.R --lite --from raw  # every script in minutes, in a copy outside the repository
Rscript code/check_consistency.R          # the outputs agree with one another
```

The full run (about 4.5 hours from raw, mostly the model audits) makes the canonical outputs, the numbers in the paper. A `--lite` run uses small forests and few repeats (`code/lib/run_mode.R`) in a fresh copy of the repository with `data/raw/` linked read-only; it tests that the pipeline works, and its numbers are not results. Each script's log is in `outputs/logs/`.

## How it is organised

`code/pipeline.csv` lists every step in run order with its stage and purpose; `code/run_all.R` runs it.

- **A, raw → processed** (`code/01_import/`, `code/02_flux/`): flux fitting from the analyzer traces with goFlux, from the hand-picked closure windows and recorded model choices in `data/raw/flux_windows/`, then the closure screen, analyzer precision and detection limits; GC calibration of internal gas; ddPCR import; tree-ID crosswalk and wood cores; the reduction of the PICRUSt2 predictions to pathway abundances per library.
- **B, canonical tables** (`code/tables/`): `data/tables/`, one table per kind of observation, every row kept and flagged; then the random-forest training rows, fit and held-out skill (`code/05_model/`).
- **C–E, analyses**: drivers and scaling (`04_drivers/`, `06_upscale/`), audits, molecular analyses and statistics (`05_model/`, `07_molecular/`, `09_tables_stats/`), figures (`08_figures/`).
- **F–G, archive and checks**: `data/compiled/` (`code/zenodo/`), `validate_tables.R`, `check_consistency.R`.
- **Outputs**: `outputs/figures/{main,SI,tables}/` (numbered), `outputs/data/` (result tables), `outputs/models/`, `outputs/audit/`, `outputs/logs/`.

Every cleaning rule is listed in `code/DATA_HYGIENE.md`; hand-run tools are described in `code/tools/README.md`.

**Canonical tables** (`data/tables/`; columns, units and checks in `code/tables/dictionary*.csv`, rendered in `data/compiled/README_data.md`):

| Table | One row per |
|---|---|
| `species`, `site`, `collars` | species; the stand reference; soil collar |
| `trees`, `tree_tags`, `inventory_stems` | measured tree; spelling of a tree name; censused stem (`in_stand`: the 8,194 scaled) |
| `flux_stem`, `flux_soil` | chamber deployment (fit, detection, `qc_pass`, `in_rf_training`) |
| `soil_readings`, `weather`, `moisture_survey`, `stream_points` | soil reading at a tree or collar; station hour (true time); point of the December 2020 moisture survey; GPS flag along the stream channel |
| `samples`, `ddpcr`, `qpcr` | DNA sample; sample × gene × assay × plate; sample × target (facility 16S and ITS qPCR) |
| `samples_16S`, `asv_counts_16S*`, `asv_taxonomy_16S*` | 16S library; ASV (main and felled-oak runs) |
| `picrust_pathways`, `picrust_pathway_info` | library × MetaCyc pathway; pathway |
| `wood_cores` | core section |
| `internal_gas` | internal stem gas sample |
| `isotopes` | Picarro isotope run |
| `species_traits`, `species_trait_sources` | species × trait (TRY, GWDD); contributing record |

The methanotroph classification rules are in `code/lib/methanotroph_definitions_README.md`.

## Figures and tables

Scripts are in `code/08_figures/` unless a folder is given; the assembler `zz_assemble_figures.R` numbers them and writes `outputs/figures/MANIFEST.md`. Several supplementary generators keep an older number in their file name.

| Figure | Script | Figure | Script |
|---|---|---|---|
| 1 | `fig01_temporal-flux.R` | S13–S15 | `figS11-13_picrust-heatmaps.R` |
| 2 | `fig02_height-flux.R` | S16–S17 | `figS15-16_internal-gas.R` |
| 3 | `fig03_variance-partition.R` | S18 | `figS17_isotope-sources.R` |
| 4, S8 | `fig04_gene-abundance.R` | S19 | `figS18_stem-deterioration.R` |
| 5 | `fig05_methane-cycling.R` | S20 | `figS19_black-oak-cross-sections.R` |
| 6 | `fig06_hydrogenotrophy.R` | S21 | `figS20_black-oak-methanome.R` |
| 7 | `fig07_decay-methanogenesis.R` | S22 | `figS21_scale-dependent-genes.R` |
| 8, S24 | `fig08_radial-species.R` | S23 | `figS22_radial-sections.R` |
| 9 | `fig09_budget.R` | S25 | `figS28_plant-traits.R` |
| S1 | `figS01_moisture-overlay.R` | S26 | `figS24_rf-model-summary.R` |
| S2a, S2b | photographs (`data/raw/photos/`) | S27 | `figS25_rf-calibration.R` |
| S3 | `figS03_detection-limits.R` | S28 | `figS27_scaling-profiles.R` |
| S4 | `figS04_mcra-probe-validation.R` | S29 | `figS26_scaling-heatmap.R` |
| S5 | `figS05_16s-ddpcr-association.R` | Methods S3 | `methodsS3_extraction_recovery.R` |
| S6 | `figS06_height-slope-moisture.R` | Table 1 | `09_tables_stats/02_campaign_counts.R` |
| S7 | `09_tables_stats/09_tree-distribution.R` | Table S2 | `09_tables_stats/04_known-putative-table.R` |
| S9 | `figS08_pmoa-mmox.R` | Table S3 | `09_tables_stats/27_table-s3-mmo-screen.R` |
| S10 | `figS14_taxonomy-mcra.R` | Table S4 | `07_molecular/03_methanotroph_resolution.R` |
| S11 | `figS09_taxonomy-pmoa.R` | Table S5 | `09_tables_stats/05_table_ddpcr-16s-concordance.R` |
| S12 | `figS10_faprotax.R` | Table S6 | `09_tables_stats/08_dbh_by_species_campaign.R` |

Tables S2 and S6 also go to `outputs/tables/SI_Tables_S2_S6.xlsx` (`09_tables_stats/23_si-table-workbook.R`); Table S8 is `data/tables/species_trait_sources.csv`. Every number quoted in the text is computed by `09_tables_stats/20_manuscript_statistics.R` into `outputs/audit/manuscript_statistics.txt`.

## License and citation

Code: MIT (`LICENSE`). Data: CC-BY 4.0. Please cite the paper and the archive, [10.5281/zenodo.18779715](https://doi.org/10.5281/zenodo.18779715), as given in `CITATION.cff`.
