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
These scripts can be run in any order after Stage 2. See `01_gene_flux_linear_models.R` for the main gene-flux analysis, `07_picrust_pathway_associations.R` for PICRUSt2 MetaCyc pathway–gene LMER associations, plus `methanotrophs/` and `methanogens/` subfolders for organism-specific analyses.

### Stage 4: Figures & Maps

These scripts generate publication figures. Run after Stages 2-3.

| Folder | Key scripts |
|--------|-------------|
| `code/04_scaling/` | `04`-`09_*.R` (diagnostics, maps, publication plots) |
| `code/06_figures/` | `01`-`12_*.R` (correlation, radial, variance, timeseries, heatmaps) |
| `code/07_maps/` | `01`-`05_*.R` (spatial interpolation, seasonal maps) |
| `code/01_flux_processing/static/` | `04_height_effect_analysis.R` (Figure 2) |
| `code/02_ddpcr/` | `04_species_barplots.R` (Figure 4) |

See `code/04_scaling/RF_CH4_workflow_spec.md` for the detailed RF technical specification.

## Figure-Script Reference

| Figure | Output File | Description | Script |
|--------|------------|-------------|--------|
| Fig 1 | `fig1_temporal_flux_timeseries.png` | Temporal flux across hydrological gradient | `06_figures/06_soil_tree_timeseries.R` |
| Fig 2 | `fig2_height_dependent_flux.png` | Height-dependent flux patterns | `01_flux_processing/static/04_height_effect_analysis.R` |
| Fig 3 | `fig3_variance_partitioning.png` | Variance partitioning | `06_figures/04_variance_partition.R` |
| Fig 4 | `fig4_methanogen_methanotroph_abundance.png` | Methanogen/methanotroph abundance | `02_ddpcr/util_combined_plot.R` |
| Fig 5 | `fig5_combined_methane_cycling_composition.png` | Combined methane-cycling 16S composition | `06_figures/08c_combined_methane_cycling_composition.R` |
| Fig 6 | `fig6_picrust_mcra_no_mcra_heatmap.png` | MetaCyc pathway × mcrA associations (no-mcrA OTU) | `06_figures/12b_picrust_pathway_heatmap.R` |
| Fig 7 | `fig7_felled_oak_profiles.png` | Felled oak vertical profiles | `06_figures/09_felled_oak_profiles.R` |
| Fig 8 | `fig8_radial_species_comparison.png` | Radial cross-sections + species comparison | `05_gene_flux_analysis/methanotrophs/04_species_gene_flux.R` |
| Fig 9 | `fig9_upscaled_flux_seasonal.png` | Upscaled seasonal flux overview | `04_scaling/09_upscale_publication_plots.R` |
| S1 | `figS1_moisture_overlay.png` | Moisture interpolation overlay | — |
| S2 | `figS2_faprotax_heatmaps.png` | FAPROTAX functional heatmaps | `06_figures/08d_faprotax_heatmaps.R` |
| S3 | `figS3_taxonomy_mcra_heatmap.png` | Family-level 16S × mcrA associations | `06_figures/12a_taxonomy_mcra_heatmap.R` |
| S4 | `figS4_picrust_mcra_all_heatmap.png` | MetaCyc pathway × mcrA (full FDR < 0.01 set) | `06_figures/12b_picrust_pathway_heatmap.R` |
| S5 | `figS5_picrust_pmoa_heatmap.png` | MetaCyc pathway × pmoA associations (no-pmoA OTU) | `06_figures/12b_picrust_pathway_heatmap.R` |
| S6 | `figS6_internal_gas_beeswarm.png` | Internal gas beeswarm by species | `06_figures/05_internal_gas_plots.R` |
| S7 | `figS7_internal_gas_profiles.png` | Internal gas multi-panel profiles | `06_figures/05_internal_gas_plots.R` |
| S8 | `figS8_d13ch4_rainfall.png` | δ¹³CH₄ vs rainfall | `06_figures/11a_isotope_d13ch4_single.R` |
| S9 | `figS9_rf_predictions.png` | RF model predictions (3-row layout) | `04_scaling/08_rf_publication_plots.R` |
| S10 | `figS10_scale_dependent_genes.png` | Scale-dependent gene–flux patterns | `05_gene_flux_analysis/methanotrophs/01_summary_stats.R` |
| S11 | `figS11_tree_radial_sections.pdf` | Tree radial mcrA cross-sections | `06_figures/02_radial_cross_sections.R` |
| S12 | `figS12_methanotroph_abundance_patterns.pdf` | pmoA vs mmoX + ratio analysis | `05_gene_flux_analysis/methanotrophs/06_pmoa_mmox_analysis.R` |
| S13 | `figS13_mcra_vs_methanotroph.png` | mcrA vs methanotroph independence | `05_gene_flux_analysis/methanotrophs/01_summary_stats.R` |
| S14 | `figS14_black_oak_methanome.png` | Black oak methanome heatmap | `06_figures/10_black_oak_methanome_heatmap.R` |
| S15 | `figS15_taxonomy_pmoa_heatmap.png` | Family-level 16S × pmoA associations | `06_figures/12c_taxonomy_pmoa_heatmap.R` |

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
- `data.table` - Fast I/O for large PICRUSt2 contribution tables
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
