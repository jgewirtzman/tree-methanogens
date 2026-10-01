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
  "code/08_figures/06_soil_tree_timeseries.R"                  = "fig1",
  "code/02_flux/static/04_height_effect_analysis.R"            = "fig2",
  "code/08_figures/util_combined_plot.R"                       = "fig4",
  "code/08_figures/08c_combined_methane_cycling_composition.R" = "fig5",
  "code/08_figures/12b_picrust_pathway_heatmap.R"              = "fig6",
  "code/08_figures/09_felled_oak_profiles.R"                   = "fig7",
  "code/07_molecular/04_species_gene_flux.R"                   = "fig8",
  "code/08_figures/12c_taxonomy_pmoa_heatmap.R"                = "S2",
  "code/08_figures/08d_faprotax_heatmaps.R"                    = "S3",
  "code/08_figures/12a_taxonomy_mcra_heatmap.R"                = "S6",
  "code/08_figures/05_internal_gas_plots.R"                    = "S7,S8",
  "code/08_figures/11a_isotope_d13ch4_single.R"                = "S9",
  "code/07_molecular/methanotrophs/03_pmoa_mmox_analysis.R"    = "S10",
  "code/08_figures/10_black_oak_methanome_heatmap.R"           = "S12",
  "code/08_figures/02_radial_cross_sections.R"                 = "S13",
  "code/08_figures/08_rf_publication_plots.R"                  = "S15",
  "code/08_figures/05_methods_figure_map.R"                    = "S1")

# Revision generators. NAMED, NOT GLOBBED: a glob on a filename prefix empties
# silently when the prefix changes.
FIGURE_SCRIPTS_REVISION <- c(
  "code/08_figures/fig01_temporal-flux.R",
  "code/08_figures/fig02_height-flux.R",
  "code/08_figures/fig02a_axis-support.R",
  "code/08_figures/fig03_variance-partition.R",
  "code/08_figures/fig04_gene-abundance.R",
  "code/08_figures/fig05_methane-cycling.R",
  "code/08_figures/fig06_hydrogenotrophy.R",
  "code/08_figures/fig07_decay-methanogenesis.R",
  "code/08_figures/fig07b_copies-per-g.R",
  "code/08_figures/fig07c_flux-unit.R",
  "code/08_figures/fig09_budget.R",
  "code/08_figures/figS02_height-slope-moisture.R",
  "code/08_figures/figS04_pmoa-mmox-coupling.R",
  "code/08_figures/figS11_scale-dependent.R",
  "code/08_figures/figS12_isotope-sources.R",
  "code/08_figures/figS15_black-oak-methanome.R",
  "code/08_figures/figS17_plant-traits.R",
  "code/08_figures/figS19_mcra-probe-validation.R",
  "code/08_figures/figS20_stem-deterioration.R",
  "code/08_figures/figS21_rf-model-summary.R",
  "code/08_figures/figSI_detection.R",
  "code/08_figures/figS_black-oak-cross-sections.R",
  "code/08_figures/figS_rf-calibration.R",
  "code/08_figures/fig_height_curves.R",
  "code/08_figures/fig_scaling_diagnostics.R",
  "code/08_figures/fig_scaling_heatmap.R",
  "code/08_figures/fig_scaling_profiles.R")
