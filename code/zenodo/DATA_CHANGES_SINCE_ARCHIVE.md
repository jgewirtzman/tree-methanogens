# Data changes not yet in the Zenodo archive

`data/` is git-ignored: it is the Zenodo drop-in. Anything below exists only on the
local/Drive copy until the next archive upload. Code to reproduce every derived file is
committed; the **raw** additions cannot be regenerated and must be uploaded.

Recorded 2026-09-30.

## 1. New raw inputs — must be uploaded (no script can recreate them)

| file | source |
|---|---|
| `data/raw/inventory/fg_2018_multiple_stems.csv` | sheet "multiple stems" of `deprecated/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW.xlsx` (md5 e981d571732545b4c7c04591dd751ae3). 40 stems on 17 multi-stem 2018 trees; without it they are lost. |
| `data/raw/inventory/fg_2018_codes.csv` | sheet "codes" (census code legend; no dead code) |
| `data/raw/inventory/fg_2018_problems.csv` | sheet "problems" |
| `data/raw/inventory/README_fg_2018_sheets.txt` | provenance note |

## 2. Renamed files — an older archive still has the old names

Run `Rscript code/02_flux/migrate_campaign_filenames.R` on a restored archive; it renames by content.

| old name | new name |
|---|---|
| `methanogen_tree_flux_complete_dataset.csv` | `tree_flux_2023_cross_species.csv` |
| `CH4_best_flux_lgr_results.csv` | `CH4_best_flux_lgr_results_2021_multiheight.csv` |
| `CO2_best_flux_lgr_results.csv` | `CO2_best_flux_lgr_results_2021_multiheight.csv` |
| `CH4_flux_lgr_results.csv` | `CH4_flux_lgr_results_2021_multiheight.csv` |
| `CO2_flux_lgr_results.csv` | `CO2_flux_lgr_results_2023_cross_species.csv` (it held 2023 data) |
| `lgr_manual_identification_results.csv` | `lgr_manual_identification_2021_multiheight.csv` |
| `lgr_manual_identification_results_final.csv` | `lgr_manual_identification_2021_multiheight_final.csv` |
| `lgr_manual_identification_summary_final.csv` | `lgr_manual_identification_summary_2021_multiheight_final.csv` |

## 3. Regenerated — reproducible from committed code, in this order

1. `Rscript code/02_flux/assemble_campaign_flux.R 2021_multiheight` → `data/processed/flux/tree_flux_2021_multiheight.csv`
2. `Rscript code/03_merge/compile_soil_env.R` → `data/processed/environmental/soil_env_by_collar.csv` (now includes `Soil_temp_moisture_2020.xlsx`: 291 records, 22 dates)
3. `(cd code/01_import && Rscript 03_process_internal_gas.R)` → `data/processed/internal_gas/{sample_data_only,processed_GC_data_internal_conc,internal_gas_calibration_check}.csv`
4. `Rscript code/01_import/03b_process_internal_gas_2024.R` → `data/processed/internal_gas/stem_gas_2024_calibrated.csv`
5. `(cd code/03_merge && Rscript 02_harmonize_all_data.R)` → `data/processed/integrated/merged_tree_dataset_final.csv` (needs a UTF-8 locale)
6. `(cd code/05_model && Rscript 01_load_and_prep_data.R && Rscript 02_rf_models.R)` → `data/processed/integrated/rf_workflow_input_data_with_2023.RData`, `outputs/models/*`
7. `Rscript code/run_all.R` → every output, figure and audit
8. `Rscript code/zenodo/compile_zenodo_datasets.R` → `data/compiled/*` (then upload `data/`)

## 4. Pending from another session

The analyzer-volume correction (GLA131 cell 28 cm³, not goFlux's example 70 cm³) is
uncommitted in worktree `.claude/worktrees/admiring-saha-012e11`. When it lands it changes
the auxfiles and every flux by 0.5–6.8 %, so steps 1, 5–8 rerun after it.
