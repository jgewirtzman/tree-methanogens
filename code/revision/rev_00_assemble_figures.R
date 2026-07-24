# ==============================================================================
# REVISION — assemble the numbered manuscript figure set.
# Copies each final figure into outputs/revision/figures/{main,SI,photos}/ with its
# manuscript number + slug. Sources: the rev_* generators write to outputs/revision/;
# a few unchanged figures come from the original pipeline's outputs/figures/
# (run generate_all_figures.R to (re)create those). Run after the generators
# (see run_all.R). Idempotent.
#
# SI order = reference order agreed with Jon (see notes/REVISION_INVENTORY.md).
# PROVISIONAL: two placements still open — S16 black-oak methanome and S21 radial
# sections — finalize against the text's citation order.
# ==============================================================================
mainD <- "outputs/revision/figures/main"; siD <- "outputs/revision/figures/SI"
photoD <- "outputs/revision/figures/photos"
for (d in c(mainD, siD, photoD)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  unlink(list.files(d, "\\.png$", full.names = TRUE))   # clear stale/renamed figures each run
}

# dest basename  <-  source path
MAIN <- c(
  "Figure_1_temporal-flux"               = "outputs/revision/fig1_final.png",
  "Figure_2_height-flux"                 = "outputs/revision/fig2_final.png",
  "Figure_3_variance-partition"          = "outputs/revision/fig3_final.png",
  "Figure_4_gene-abundance"              = "outputs/revision/fig4_final.png",
  "Figure_5_methane-cycling-composition" = "outputs/revision/fig5_final.png",
  "Figure_6_hydrogenotrophy"             = "outputs/revision/fig_hydrogenotrophy.png",
  "Figure_7_decay-methanogenesis"        = "outputs/revision/fig7_decay_methanogenesis.png",  # NEW expanded Fig 7 (folds old felled-oak)
  "Figure_8_radial-species"              = "outputs/figures/main/fig8_radial_species_comparison.png",  # unchanged (orig pipeline); central takeaway
  "Figure_9_ch4-budget"                  = "outputs/revision/fig_budget_maps.png")

# SI — reference order (each cited near its main-text section)
SI <- c(
  # Flux & environment (Figs 1-2)
  "Figure_S01_moisture-overlay"          = "outputs/figures/supplementary/figS1_moisture_overlay.png",
  "Figure_S02_height-slope-moisture"     = "outputs/revision/height_slope_vs_moisture.png",
  # Genes & abundance (Fig 4)
  "Figure_S03_ddpcr-mcra-probe-validation" = "outputs/revision/figS19_mcra_probe_validation.png",
  "Figure_S04_ddpcr-16s-concordance"     = "outputs/revision/figS21_ddpcr_16s_concordance.png",
  "Figure_S05_mcra-vs-methanotroph"      = "outputs/figures/supplementary/figS14_mcra_vs_methanotroph.png",
  "Figure_S06_pmoa-mmox-coupling"        = "outputs/revision/figS10_final.png",
  "Figure_S07_pmoa-mmox-by-compartment"  = "outputs/revision/SI_fig_pmoa_mmox_separate.png",
  # Composition / taxonomy (Fig 5)
  "Figure_S08_taxonomy-mcra"             = "outputs/figures/supplementary/figS6_taxonomy_mcra_heatmap.png",
  "Figure_S09_taxonomy-pmoa"             = "outputs/figures/supplementary/figS2_taxonomy_pmoa_heatmap.png",
  # Function / pathway (Fig 6)
  "Figure_S10_faprotax"                  = "outputs/figures/supplementary/figS3_faprotax_heatmaps.png",
  "Figure_S11_picrust-mcra-all"          = "outputs/figures/supplementary/figS4_picrust_mcra_all_heatmap.png",
  "Figure_S12_picrust-mcra-no-mcra"      = "outputs/figures/main/fig6_picrust_mcra_no_mcra_heatmap.png",  # demoted old main Fig 6
  "Figure_S13_picrust-pmoa"              = "outputs/figures/supplementary/figS5_picrust_pmoa_heatmap.png",
  # Within-tree / decay (Fig 7)
  "Figure_S14_internal-gas-beeswarm"     = "outputs/figures/supplementary/figS7_internal_gas_beeswarm.png",
  "Figure_S15_internal-gas-profiles"     = "outputs/figures/supplementary/figS8_internal_gas_profiles.png",
  "Figure_S16_black-oak-methanome"       = "outputs/revision/black_oak_methanome_revised.png",  # PLACEMENT TBC
  "Figure_S17_stem-deterioration"        = "outputs/revision/figS20_stem_deterioration.png",    # supports Fig 7a
  # Isotopes (right after internal gas)
  "Figure_S18_d13ch4-raincloud"          = "outputs/figures/supplementary/figS9_d13ch4_rainfall.png",
  "Figure_S19_isotope-sources"           = "outputs/revision/SI_isotopes_source_composite.png",
  # Gene-flux scaling (Fig 8)
  "Figure_S20_scale-dependent-genes"     = "outputs/revision/figS11_final.png",
  "Figure_S21_radial-sections"           = "outputs/figures/supplementary/figS13_tree_radial_sections.png",  # PLACEMENT TBC
  # Upscaling (Fig 9)
  "Figure_S22_rf-flux-predictions"       = "outputs/figures/supplementary/figS15_rf_predictions.png",
  # Plant traits (discussion)
  "Figure_S23_plant-traits"              = "outputs/revision/traits_heatmap_robust.png")

# Photo plates — separate section (NOT SI data figures); chamber photos to be added
PHOTOS <- c(
  "Plate_black-oak-cross-sections"       = "outputs/revision/figS_black_oak_cross_sections.png")

copy_set <- function(map, dest) {
  miss <- 0
  for (nm in names(map)) {
    src <- map[[nm]]
    if (file.exists(src)) file.copy(src, file.path(dest, paste0(nm, ".png")), overwrite = TRUE)
    else { cat("  MISSING:", src, "->", nm, "\n"); miss <- miss + 1 }
  }
  miss
}
m1 <- copy_set(MAIN, mainD); m2 <- copy_set(SI, siD); m3 <- copy_set(PHOTOS, photoD)
cat(sprintf("Assembled %d main + %d SI + %d photo (%d missing; missing = original-pipeline figs, run generate_all_figures.R).\n",
            length(MAIN)-m1, length(SI)-m2, length(PHOTOS)-m3, m1+m2+m3))
