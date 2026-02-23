# Tree Microbiomes and Methane Emissions in Upland Forests

Analysis code and data processing pipeline for the study of methane fluxes from tree stems and soils at Yale-Myers Forest (YMF), Connecticut, and their relationship to methanogenic and methanotrophic microbial communities.

**Preprint:** [bioRxiv 10.1101/2025.09.30.679632v2](https://www.biorxiv.org/content/10.1101/2025.09.30.679632v2)

## Overview

This project integrates:
- **Flux measurements**: 1,148 stem CH4 flux measurements from 482 trees + 276 soil measurements (2020-2023) using semi-rigid and rigid chamber designs with an LGR3 off-axis ICOS analyzer
- **Molecular quantification**: ddPCR of mcrA (methanogens), pmoA/mmoX (methanotrophs) from 564 wood and soil samples
- **Community characterization**: 16S rRNA amplicon sequencing, PICRUSt2, FAPROTAX
- **Spatial scaling**: Random Forest models to upscale plot-level CH4 flux across the forest inventory

## Directory Structure

```
tree-methanogens/
├── data/
│   ├── raw/                    # Unprocessed instrument and field data
│   │   ├── lgr/                # LGR3 analyzer output (.txt files)
│   │   ├── inventory/          # ForestGEO census data
│   │   ├── dbh/                # Tree diameter measurements
│   │   ├── tree_cores/         # Increment core data
│   │   ├── ddpcr/              # Digital droplet PCR raw data
│   │   ├── ddpcr_raw/          # Collated ddPCR gene data
│   │   ├── internal_gas/       # Picarro & GC internal gas data
│   │   ├── weather/            # YMF met tower data
│   │   ├── 16s/                # 16S amplicon OTU/taxonomy tables
│   │   │   └── black_oak/      # Black Oak 16S amplicon data (QUVE)
│   │   ├── picrust/            # PICRUSt2 functional predictions & metadata
│   │   ├── environmental/      # Environmental covariates
│   │   ├── external/           # External comparison data (Eastern US)
│   │   ├── flux_chamber_corrections/ # Chamber dimension corrections
│   │   └── field_data/         # iPad field data, soil measurements
│   │       └── black_oak/      # Felled oak flux, GC, and standards data
│   └── processed/              # Cleaned/derived datasets
│       ├── flux/               # Calculated CH4/CO2 flux results
│       ├── molecular/          # Processed ddPCR data
│       ├── tree_data/          # DBH consensus, tree ID mapping
│       ├── tree_cores/         # Processed wood property data
│       ├── internal_gas/       # Processed GC/isotope data
│       ├── integrated/         # Multi-source merged datasets
│       ├── metadata/           # Phylogenetic data (PhytoPhylo)
│       └── soil/               # Soil flux data
├── code/                       # See "Pipeline" section below
├── outputs/
│   ├── figures/                # All generated figures
│   │   ├── main/               #   Main text figures
│   │   └── supplementary/      #   Supplementary figures
│   ├── tables/                 # Summary statistics & model tables
│   ├── models/                 # Saved RF model objects (.RData)
│   └── flux_predictions/       # Upscaled monthly/annual flux estimates
└── deprecated/                 # Superseded scripts and data (not tracked)
```

## Pipeline

Scripts are numbered within each folder to indicate execution order. Scripts prefixed `util_` are optional utilities/diagnostics.

### Stage 1: Raw Data Processing (parallel tracks)

These scripts process raw instrument data into analysis-ready datasets. The four tracks below can be run independently.

#### Track A — Tree metadata (`code/03_tree_data/`)
| # | Script | Produces |
|---|--------|----------|
| 01 | `01_tree_id_consensus.R` | `tree_id_comprehensive_mapping.csv`, `tree_dbh_consensus_comprehensive.csv` |
| 02 | `02_process_wood_cores.R` | `tree_core_filled_complete.csv` |
| 03 | `03_process_internal_gas.R` | `sample_data_only.csv` |

#### Track B — Molecular data (`code/02_ddpcr/`)
| # | Script | Produces |
|---|--------|----------|
| 01 | `01_import_and_align.R` | `processed_ddpcr_data.csv` |
| 02-03 | `02_check_completeness.R`, `03_confirm_import.R` | QC validation |

#### Track C — Semi-rigid flux (`code/01_flux_processing/semirigid/`)
| # | Script | Produces |
|---|--------|----------|
| 01 | `01_calc_chamber_dims.R` | Chamber geometry calculations |
| 02 | `02_join_flux_geometry.R` | `flux_with_geometry_fixed.csv` |
| 03 | `03_prep_tree_auxfile.R` / `03_prep_soil_auxfile.R` | goFlux auxfiles |
| 04 | `04_goflux_trees.R` / `04_goflux_soils.R` | `semirigid_tree_final_complete_dataset.csv`, `*_soil.csv` |
| 05 | `05_fix_december.R` | Updates flux files with Dec 2020 data |

#### Track D — Static chamber flux (`code/01_flux_processing/static/`)
| # | Script | Produces |
|---|--------|----------|
| 01 | `01_prep_auxfile.R` / `01_prep_auxfile_2023.R` | goFlux auxfiles |
| 02 | `02_goflux_trees_2021.R` / `02_goflux_trees_2023.R` | `methanogen_tree_flux_complete_dataset.csv` |
| 03 | `03_fix_failed.R` | Corrected flux identifications |

### Stage 2: Data Integration (`code/00_harmonization/`)

Run **after** all Stage 1 tracks complete.

| # | Script | Produces |
|---|--------|----------|
| 01 | `01_fix_soil_flux.R` | Corrected soil flux dataset |
| 02 | `02_harmonize_all_data.R` | **`merged_tree_dataset_final.csv`** (master dataset) |

### Stage 3: Analysis & Upscaling

#### Random Forest workflow (`code/04_scaling/`)
| # | Script | Produces |
|---|--------|----------|
| 01 | `01_load_and_prep_data.R` | `rf_workflow_input_data_with_2023.RData` |
| 02 | `02_rf_models.R` | `RF_MODELS.RData`, `TAXONOMY_PRIORS.RData`, flux tables |
| 03 | `03_predict_tree_flux.R` | `tree_monthly_predictions.RData` |
| 03 | `03_predict_soil_flux.R` | `soil_monthly_predictions.RData` |

#### Gene-flux analysis (`code/05_gene_flux_analysis/`)
These scripts can be run in any order after Stage 2. See `01_gene_flux_linear_models.R` for the main gene-flux analysis, plus `methanotrophs/` and `methanogens/` subfolders for organism-specific analyses.

### Stage 4: Figures & Maps

These scripts generate publication figures. Run after Stages 2-3.

| Folder | Key scripts |
|--------|-------------|
| `code/04_scaling/` | `04`-`09_*.R` (diagnostics, maps, publication plots) |
| `code/06_figures/` | `01`-`07_*.R` (correlation, radial, variance, timeseries) |
| `code/07_maps/` | `01`-`05_*.R` (spatial interpolation, seasonal maps) |
| `code/01_flux_processing/static/` | `04_height_effect_analysis.R` (Figure 2) |
| `code/02_ddpcr/` | `04_species_barplots.R` (Figure 4) |

See `code/04_scaling/RF_CH4_workflow_spec.md` for the detailed RF technical specification.

## Figure-Script Reference

| Figure | Description | Script |
|--------|-------------|--------|
| Fig 1 | Temporal flux across hydrological gradient | `04_scaling/08_rf_publication_plots.R`, `06_figures/06_soil_tree_timeseries.R` |
| Fig 2 | Height-dependent flux patterns | `01_flux_processing/static/04_height_effect_analysis.R` |
| Fig 3 | Variance partitioning | `06_figures/04_variance_partition.R`, `06_figures/01_correlation_plots.R` |
| Fig 4 | Methanogen/methanotroph abundance | `02_ddpcr/04_species_barplots.R` |
| Fig 5 | Methanogen/methanotroph 16S composition | `06_figures/08_methanogen_16s_composition.R`, `06_figures/08b_methanotroph_16s_composition.R`, `06_figures/08c_combined_methane_cycling_composition.R` |
| Fig 6 | Bacterial family heatmap | `05_gene_flux_analysis/methanotrophs/07_picrust_analysis.R` |
| Fig 7 | Felled oak vertical profiles | `06_figures/09_felled_oak_profiles.R` |
| Fig 8 | Radial cross-sections + gene-flux | `06_figures/02_radial_cross_sections.R`, `05_gene_flux_analysis/methanotrophs/05_radial_gene_plots.R` |
| Fig 9 | RF seasonal flux maps | `04_scaling/07_seasonal_maps.R`, `04_scaling/08_rf_publication_plots.R` |
| S1-S2 | Methanotroph composition | `05_gene_flux_analysis/methanotrophs/02_rf_models.R` |
| S3 | PICRUSt2 pathways | `05_gene_flux_analysis/methanotrophs/07_picrust_analysis.R` |
| S4-S5 | Internal CH4 vs mcrA | `06_figures/05_internal_gas_plots.R`, `05_gene_flux_analysis/methanotrophs/06_pmoa_mmox_analysis.R` |
| S6 | RF model predictions | `04_scaling/08_rf_publication_plots.R` |
| S7 | Individual gene-flux models | `05_gene_flux_analysis/methanotrophs/04_species_gene_flux.R` |
| S8 | Species gene-flux correlations | `05_gene_flux_analysis/methanotrophs/04_species_gene_flux.R` |
| S9 | Individual tree cross-sections | `06_figures/02_radial_cross_sections.R`, `06_figures/03_threshold_analysis.R` |
| S10 | pMMO:sMMO balance | `05_gene_flux_analysis/methanotrophs/06_pmoa_mmox_analysis.R` |
| S11 | Methanogen-methanotroph independence | `05_gene_flux_analysis/methanotrophs/04_species_gene_flux.R` |
| S12 | Black oak methanome heatmap | `06_figures/10_black_oak_methanome_heatmap.R` |
| Isotope (Dataset 1) | δ¹³CH₄ — single dataset | `06_figures/11a_isotope_d13ch4_single.R` |
| Isotope (Dataset 2) | δ¹³CH₄ — paired dataset | `06_figures/11b_isotope_d13ch4_paired.R` |
| Taxonomy heatmap | Family-level 16S × mcrA associations | `06_figures/12a_taxonomy_mcra_heatmap.R` |
| PICRUSt heatmap | MetaCyc pathway × mcrA associations | `06_figures/12b_picrust_pathway_heatmap.R` |

## Key Datasets

| File | Location | Description |
|------|----------|-------------|
| `merged_tree_dataset_final.csv` | `data/processed/integrated/` | Master dataset merging all tree-level measurements |
| `rf_workflow_input_data_with_2023.RData` | `data/processed/integrated/` | Integrated data ready for RF modeling |
| `processed_ddpcr_data.csv` | `data/processed/molecular/` | ddPCR gene quantification results |
| `methanogen_tree_flux_complete_dataset.csv` | `data/processed/flux/` | Tree flux + methanogen data combined |
| `tree_id_comprehensive_mapping.csv` | `data/processed/tree_data/` | Authoritative tree ID cross-reference |
| `RF_MODELS.RData` | `outputs/models/` | Trained tree and soil RF models |
| `MONTHLY_FLUXES.csv` | `outputs/flux_predictions/` | RF-predicted monthly plot-level fluxes |
| `ANNUAL_SUMMARY.csv` | `outputs/flux_predictions/` | Annual flux totals with uncertainty |

## Dependencies

Key R packages:
- `goFlux` - Flux calculation from continuous gas analyzer data
- `randomForest` - CH4 flux prediction models
- `vegan` - Variance partitioning (varpart)
- `tidyverse` - Data manipulation and visualization
- `sf`, `terra`, `akima` - Spatial analysis and interpolation
- `ape`, `phytools` - Phylogenetic analysis
- `lme4`, `lmerTest` - Mixed-effects models
- `pheatmap` - Heatmap visualization
- `ggridges`, `patchwork`, `cowplot` - Figure layout

## Data Availability

Raw data files are archived on Zenodo: **[DOI forthcoming]**

To reproduce analyses, download the data archive and extract its contents into `data/raw/`. The directory structure is preserved in git via `.gitkeep` files so all paths are pre-configured. Processed data in `data/processed/` can be regenerated by running the pipeline scripts in order.

## Notes

- Scripts use relative paths from the project root (requires opening `tree-methanogens.Rproj` in RStudio)
- All input data is stored within `data/` — no external path dependencies
- Data file *contents* are excluded from git tracking (archived on Zenodo); only `.gitkeep` files are tracked to preserve directory structure
- The `deprecated/` folder contains superseded script versions and old data caches (not tracked)
- Scripts prefixed `util_` are optional diagnostics/utilities, not part of the core pipeline
