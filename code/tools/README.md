# Tools — run by hand, not by the pipeline

| script | when to run it |
|---|---|
| `click_closure_windows.R` | Interactive. Choose chamber-closure windows by hand for flux records goFlux cannot window automatically (used for the untagged monthly stems). Its products are archived in `data/processed/flux/untagged_rescue/`; rerunning means repeating the clicking. |
| `click_untagged_windows.R` | Interactive, for the 45 untagged monthly stems. After `02_flux/rescue/01_untagged_auxfile.R`, pick each closure window; save the picker's result as `data/processed/flux/untagged_rescue/untagged_manID.rds` (stage A's `02_untagged_fluxes.R` re-fits from it). The current `.rds` is the hand-picked record. |
| `migrate_campaign_filenames.R` | Once, on a data archive older than 2026-09-30, to rename the flux files to the per-campaign names (it identifies files by content). |
| `mmo_capacity_screen.R` | Queries NCBI. Rebuilds `outputs/data/mmo_capacity_screen.csv` (Table S5); needs network access and is not deterministic over time as databases grow. |
| `harvest_dictionary.R` | Once, to seed `code/zenodo/column_dictionary.csv` from hand-written column notes. The dictionary is then edited directly. |
