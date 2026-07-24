# ==============================================================================
# REVISION — assemble the numbered manuscript figure set.
# Copies each final figure into outputs/revision/figures/{main,SI}/ with its
# manuscript number + slug. Sources: the rev_* generators write to outputs/revision/;
# a few unchanged figures come from the original pipeline's outputs/figures/.
# Run after the generators (see run_all.R). Idempotent.
# ==============================================================================
mainD <- "outputs/revision/figures/main"; siD <- "outputs/revision/figures/SI"
dir.create(mainD, showWarnings = FALSE, recursive = TRUE)
dir.create(siD,   showWarnings = FALSE, recursive = TRUE)

# dest basename  <-  source path
MAIN <- c(
  "Figure_1_temporal-flux"               = "outputs/revision/fig1_final.png",
  "Figure_2_height-flux"                 = "outputs/revision/fig2_final.png",
  "Figure_3_variance-partition"          = "outputs/revision/fig3_final.png",
  "Figure_4_gene-abundance"              = "outputs/revision/fig4_final.png",
  "Figure_5_methane-cycling-composition" = "outputs/revision/fig5_final.png",
  "Figure_6_hydrogenotrophy"             = "outputs/revision/fig_hydrogenotrophy.png",
  "Figure_7_felled-oak"                  = "outputs/revision/fig7_felled_oak_final.png",
  "Figure_8_radial-species"              = "outputs/figures/main/fig8_radial_species_comparison.png",  # unchanged (orig pipeline)
  "Figure_9_ch4-budget"                  = "outputs/revision/fig_budget_maps.png")
SI <- c(
  "Figure_S01_moisture-overlay"          = "outputs/figures/supplementary/figS1_moisture_overlay.png",
  "Figure_S02_height-slope-moisture"     = "outputs/revision/height_slope_vs_moisture.png",
  "Figure_S03_mcra-vs-methanotroph"      = "outputs/figures/supplementary/figS14_mcra_vs_methanotroph.png",
  "Figure_S04_pmoa-mmox-coupling"        = "outputs/revision/figS10_final.png",
  "Figure_S05_taxonomy-mcra"             = "outputs/figures/supplementary/figS6_taxonomy_mcra_heatmap.png",
  "Figure_S06_taxonomy-pmoa"             = "outputs/figures/supplementary/figS2_taxonomy_pmoa_heatmap.png",
  "Figure_S07_faprotax"                  = "outputs/figures/supplementary/figS3_faprotax_heatmaps.png",
  "Figure_S08_picrust-mcra-all"          = "outputs/figures/supplementary/figS4_picrust_mcra_all_heatmap.png",
  "Figure_S09_picrust-mcra-no-mcra"      = "outputs/figures/main/fig6_picrust_mcra_no_mcra_heatmap.png",  # demoted main Fig 6
  "Figure_S10_picrust-pmoa"              = "outputs/figures/supplementary/figS5_picrust_pmoa_heatmap.png",
  "Figure_S11_scale-dependent-genes"     = "outputs/revision/figS11_final.png",
  "Figure_S12_isotope-sources"           = "outputs/revision/SI_isotopes_source_composite.png",
  "Figure_S13_internal-gas-beeswarm"     = "outputs/figures/supplementary/figS7_internal_gas_beeswarm.png",
  "Figure_S14_internal-gas-profiles"     = "outputs/figures/supplementary/figS8_internal_gas_profiles.png",
  "Figure_S15_black-oak-methanome"       = "outputs/revision/black_oak_methanome_revised.png",  # revised putative defs
  "Figure_S16_radial-sections"           = "outputs/figures/supplementary/figS13_tree_radial_sections.png",
  "Figure_S17_plant-traits"              = "outputs/revision/traits_heatmap_robust.png",
  "Figure_S18_rf-flux-predictions"       = "outputs/figures/supplementary/figS15_rf_predictions.png",
  "Figure_S19_mcra-probe-validation"     = "outputs/revision/figS19_mcra_probe_validation.png",
  "Figure_S20_stem-deterioration"        = "outputs/revision/figS20_stem_deterioration.png")

copy_set <- function(map, dest) {
  miss <- 0
  for (nm in names(map)) {
    src <- map[[nm]]
    if (file.exists(src)) file.copy(src, file.path(dest, paste0(nm, ".png")), overwrite = TRUE)
    else { cat("  MISSING:", src, "->", nm, "\n"); miss <- miss + 1 }
  }
  miss
}
m1 <- copy_set(MAIN, mainD); m2 <- copy_set(SI, siD)
cat(sprintf("Assembled %d main + %d SI figures (%d missing).\n", length(MAIN)-m1, length(SI)-m2, m1+m2))
