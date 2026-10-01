# Legacy upscaling chain (archived 2026-09-30)

These scripts predict and map stand-level flux from `tree_monthly_predictions.RData`
and `soil_monthly_predictions.RData` (September 2025). Those predictions predate the
model lock and were replaced by the canonical chain in `run_all.R`
(`predict_tree_flux_current.R`, `predict_soil_surface.R`, `budget_canonical.R`,
`fig09_budget.R`).

None of them was run by `run_all.R`, except `09_upscale_publication_plots.R`, which
`generate_all_figures.R` still ran. Its `fig9_upscaled_flux_seasonal.png` is not the
assembled Figure 9 (`Figure_9_ch4-budget.png`). Several of these scripts also wrote the
same figure files as each other.
