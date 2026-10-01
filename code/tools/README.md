# Tools — run by hand, not by the pipeline

| script | when to run it |
|---|---|
| `click_closure_windows.R` | Interactive. Choose chamber-closure windows by hand for flux records goFlux cannot window automatically (used for the untagged monthly stems). Its products are archived in `data/processed/flux/untagged_rescue/`; rerunning means repeating the clicking. |
| `click_untagged_windows.R` | Interactive, for the 45 untagged monthly stems. After `02_flux/rescue/01_untagged_auxfile.R`, pick each closure window; save the picker's result as `data/processed/flux/untagged_rescue/untagged_manID.rds` (stage A's `02_untagged_fluxes.R` re-fits from it). The current `.rds` is the hand-picked record. |
| `rebuild_closure_windows.R` | Run once (2026-10-01). Recovered the hand-picked windows for the monthly stems and the 2023 campaign, which had never been saved, from the stored fits (every window to machine precision), and wrote `flux_model_choices.csv`. Not needed again unless those stored fits change. |
| `repick_failed_2021_windows.R` | Interactive continuation of a 2021 picking session: re-picks the five deployments whose windows failed and writes `lgr_manual_identification_2021_multiheight_final.csv`, the 2021 record. Needs the session's objects (`ow.lgr3.complete`, `manID_batches`). |
| `migrate_campaign_filenames.R` | Once, on a data archive older than 2026-09-30, to rename the flux files to the per-campaign names (it identifies files by content). |
| `mmo_capacity_screen.R` | Queries NCBI. Rebuilds `outputs/data/mmo_capacity_screen.csv` (Table S5); needs network access and is not deterministic over time as databases grow. |
| `vendor_traits.R` | Copies per-species plant traits from the sibling tree-gas-traits analysis into `data/processed/traits/ymf_species_traits.csv`. Needs that repository's `ymf_with_traits.csv`; the vendored file is an archived input. |
| `harvest_dictionary.R` | Once, to seed `code/zenodo/column_dictionary.csv` from hand-written column notes. The dictionary is then edited directly. |
