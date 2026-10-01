# ==============================================================================
# figure_scripts.R -- the one list of figure generators.
# Sourced by code/make_figures.R (which runs them) and code/run_all.R (whose
# reachability check must know they are run). Two copies of this list existed
# until 2026-09-30, one per runner, and they had drifted: the old runner still
# named the archived pre-revision Figure 3 script.
# ==============================================================================
# Original-pipeline generators: figures carried over unchanged, labelled with the
# figure they feed.
FIGURE_SCRIPTS_ORIGINAL <- c(
  "code/archive/superseded_2026-10-01/08_figures/06_soil_tree_timeseries.R"                  = "fig1",
  "code/09_tables_stats/01_height_effects.R"            = "fig2",
  "code/archive/superseded_2026-10-01/08_figures/util_combined_plot.R"                       = "fig4",
  "code/archive/superseded_2026-10-01/08_figures/08c_combined_methane_cycling_composition.R" = "fig5",
  "code/08_figures/figS11-13_picrust-heatmaps.R"              = "fig6",
  "code/archive/superseded_2026-10-01/08_figures/09_felled_oak_profiles.R"                   = "fig7",
  "code/08_figures/fig08_radial-species.R"                   = "fig8",
  "code/08_figures/figS09_taxonomy-pmoa.R"                = "S2",
  "code/08_figures/figS10_faprotax.R"                    = "S3",
  "code/08_figures/figS14_taxonomy-mcra.R"                = "S6",
  "code/08_figures/figS15-16_internal-gas.R"                    = "S7,S8",
  "code/archive/superseded_2026-10-01/08_figures/11a_isotope_d13ch4_single.R"                = "S9",
  "code/archive/superseded_2026-10-01/07_molecular/methanotrophs/03_pmoa_mmox_analysis.R"    = "S10",
  "code/archive/superseded_2026-10-01/08_figures/10_black_oak_methanome_heatmap.R"           = "S12",
  "code/08_figures/figS22_radial-sections.R"                 = "S13",
  "code/archive/superseded_2026-10-01/08_figures/08_rf_publication_plots.R"                  = "S15",
  "code/08_figures/figS01_moisture-overlay.R"                    = "S1")

# Revision generators. NAMED, NOT GLOBBED: a glob on a filename prefix empties
# silently when the prefix changes.
FIGURE_SCRIPTS_REVISION <- c(
  "code/08_figures/fig01_temporal-flux.R",
  "code/08_figures/fig02_height-flux.R",
  "code/archive/superseded_2026-10-01/08_figures/fig02a_axis-support.R",
  "code/08_figures/fig03_variance-partition.R",
  "code/08_figures/fig04_gene-abundance.R",
  "code/08_figures/fig05_methane-cycling.R",
  "code/08_figures/fig06_hydrogenotrophy.R",
  "code/08_figures/fig07_decay-methanogenesis.R",
  "code/archive/superseded_2026-10-01/08_figures/fig07b_copies-per-g.R",
  "code/archive/superseded_2026-10-01/08_figures/fig07c_flux-unit.R",
  "code/08_figures/fig09_budget.R",
  "code/08_figures/figS06_height-slope-moisture.R",
  "code/08_figures/figS08_pmoa-mmox.R",
  "code/08_figures/figS21_scale-dependent-genes.R",
  "code/08_figures/figS17_isotope-sources.R",
  "code/08_figures/figS20_black-oak-methanome.R",
  "code/08_figures/figS28_plant-traits.R",
  "code/08_figures/figS04_mcra-probe-validation.R",
  "code/08_figures/figS18_stem-deterioration.R",
  "code/08_figures/figS24_rf-model-summary.R",
  "code/08_figures/figS03_detection-limits.R",
  "code/08_figures/figS05_16s-ddpcr-association.R",
  "code/08_figures/figS19_black-oak-cross-sections.R",
  "code/08_figures/figS25_rf-calibration.R",
  "code/archive/superseded_2026-10-01/08_figures/fig_height_curves.R",
  "code/archive/superseded_2026-10-01/08_figures/fig_scaling_diagnostics.R",
  "code/08_figures/figS26_scaling-heatmap.R",
  "code/08_figures/figS27_scaling-profiles.R")
