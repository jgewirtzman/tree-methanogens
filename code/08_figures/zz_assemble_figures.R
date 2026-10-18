# ==============================================================================
# zz_assemble_figures.R -- number and copy every figure and table of the manuscript
# ------------------------------------------------------------------------------
# Stage E, last (after every figure generator; code/pipeline.csv). Run from the
# repository root. Idempotent.
# Copies each final figure into outputs/figures/{main,SI}/ and each table into
# outputs/figures/tables/, named by manuscript number + slug, and writes
# outputs/figures/MANIFEST.md from the maps below. Sources are the outputs of the
# stage E generators (outputs/figures/generated/ and outputs/figures/original/),
# the stage D tables and two static files (photographs; the NCBI screen table).
# SI figures are numbered in order of first citation in the manuscript text.
# A missing source is reported (and counted); a source older than the start of the
# current pipeline run is listed as stale.
# ==============================================================================
mainD <- "outputs/figures/main"; siD <- "outputs/figures/SI"
photoD <- "outputs/figures/photos"; tableD <- "outputs/figures/tables"
for (d in c(mainD, siD, photoD, tableD)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  unlink(list.files(d, "\\.(png|csv)$", full.names = TRUE))   # clear stale/renamed outputs each run
}

# dest basename  <-  source path
MAIN <- c(
  "Figure_1_temporal-flux"               = "outputs/figures/generated/fig1_final.png",
  "Figure_2_height-flux"                 = "outputs/figures/generated/fig2_final.png",
  "Figure_3_variance-partition"          = "outputs/figures/generated/fig3_final.png",
  "Figure_4_gene-abundance"              = "outputs/figures/generated/fig4_final.png",
  "Figure_5_methane-cycling-composition" = "outputs/figures/generated/fig5_final.png",
  "Figure_6_hydrogenotrophy"             = "outputs/figures/generated/fig_hydrogenotrophy.png",
  "Figure_7_decay-methanogenesis"        = "outputs/figures/generated/fig7_decay_methanogenesis.png",
  "Figure_8_radial-species"              = "outputs/figures/original/main/fig8_radial_species_comparison.png",
  "Figure_9_ch4-budget"                  = "outputs/figures/generated/fig_budget_maps.png")

# SI, numbered in order of first citation in the main text; grouped by the section
# that first cites them. Figure S2 (chamber photographs) is in SI_PHOTOS.
SI <- c(
  # Methods
  "Figure_S01_moisture-overlay"          = "outputs/figures/original/supplementary/figS1_moisture_overlay.png",
  "Figure_S03_detection-limits"          = "outputs/figures/generated/fig_SI_detection.png",
  "Figure_S04_ddpcr-mcra-probe-validation" = "outputs/figures/generated/figS19_mcra_probe_validation.png",
  "Figure_S05_16s-ddpcr-association"     = "outputs/figures/generated/figSI_16s_ddpcr_association.png",
  # Results 2-3: height, species
  "Figure_S06_height-slope-moisture"     = "outputs/figures/generated/height_slope_vs_moisture.png",
  "Figure_S07_species-moisture-niche"    = "outputs/figures/generated/tree_species_moisture_niche.png",
  # Results 4-5: genes, composition
  "Figure_S08_gene-abundance-species"   = "outputs/figures/generated/figS_gene-abundance-species.png",
  "Figure_S09_pmoa-mmox"                 = "outputs/figures/generated/figS10_final.png",
  "Figure_S11_taxonomy-pmoa"             = "outputs/figures/original/supplementary/figS2_taxonomy_pmoa_heatmap.png",
  # Results 6: inferred function
  "Figure_S12_faprotax"                  = "outputs/figures/original/supplementary/figS3_faprotax_heatmaps.png",
  "Figure_S13_picrust-mcra-no-methanogen" = "outputs/figures/original/main/fig6_picrust_mcra_no_mcra_heatmap.png",
  "Figure_S14_picrust-mcra-all"          = "outputs/figures/original/supplementary/figS4_picrust_mcra_all_heatmap.png",
  "Figure_S15_picrust-pmoa"              = "outputs/figures/original/supplementary/figS5_picrust_pmoa_heatmap.png",
  "Figure_S10_taxonomy-mcra"             = "outputs/figures/original/supplementary/figS6_taxonomy_mcra_heatmap.png",
  # Results 7: internal gas, isotopes
  "Figure_S16_internal-gas-beeswarm"     = "outputs/figures/original/supplementary/figS7_internal_gas_beeswarm.png",
  "Figure_S17_internal-gas-profiles"     = "outputs/figures/original/supplementary/figS8_internal_gas_profiles.png",
  "Figure_S18_isotope-sources"           = "outputs/figures/generated/SI_isotopes_source_composite.png",
  # Results 8: decay and the felled oak
  "Figure_S19_stem-deterioration"        = "outputs/figures/generated/figS20_stem_deterioration.png",
  "Figure_S21_black-oak-methanome"       = "outputs/figures/generated/black_oak_methanome_revised.png",
  # Results 9: gene-flux across scales
  "Figure_S22_scale-dependent-genes"     = "outputs/figures/generated/figS11_final.png",
  "Figure_S23_radial-sections"           = "outputs/figures/original/supplementary/figS13_tree_radial_sections.png",
  "Figure_S24_mcra-vs-methanotroph"      = "outputs/figures/original/supplementary/figS14_mcra_vs_methanotroph.png",
  # Results 10: stand-scale bounds
  "Figure_S26_rf-model-summary"          = "outputs/figures/generated/figS21_rf_model_summary.png",
  "Figure_S27_rf-calibration"            = "outputs/figures/generated/figS_rf_calibration.png",
  "Figure_S29_scaling-heatmap"           = "outputs/figures/generated/fig_scaling_heatmap.png",
  "Figure_S28_scaling-profiles"          = "outputs/figures/generated/fig_scaling_profiles.png",
  # Discussion
  "Figure_S25_plant-traits"              = "outputs/figures/generated/plant_traits.png")

# Photographs (static: exempt from the staleness check) and the photo plate of Figure S20
SI_PHOTOS <- c(
  "Figure_S02a_semirigid-chamber"        = "data/raw/photos/semirigid_chamber.jpg",
  "Figure_S02b_rigid-chamber"            = "data/raw/photos/rigid_chamber.jpg",
  "Figure_S20_black-oak-cross-sections"  = "outputs/figures/generated/figS_black_oak_cross_sections.png")

# Manuscript tables with one data source. Not assembled here: Table S1 (ddPCR assays)
# and Table S7 (literature identifications), written in the manuscript, and Table S4,
# composed in the manuscript from outputs/data/methanotroph_resolution*.csv
# (code/07_molecular/03_methanotroph_resolution.R).
TABLES <- c(
  "Table_1_campaign-summary"             = "outputs/data/campaign_counts.csv",
  "Table_S2_known-putative-taxa"         = "outputs/data/known_putative_taxa_table.csv",
  "Table_S5_ddpcr-16s-concordance"       = "outputs/data/tbl_ddpcr_16s_concordance.csv",
  "Table_S5_ddpcr-16s-concordance-view"  = "outputs/figures/generated/tbl_ddpcr_16s_concordance.png",
  "Table_S6_dbh-by-species-campaign"     = "outputs/data/dbh_by_species_campaign.csv",
  "Table_S8_species-trait-sources"       = "data/tables/species_trait_sources.csv")
# Not rebuilt by run_all: mmo_capacity_screen.R queries NCBI, so the run shown in the
# SI is committed (code/lib/mmo_capacity_screen_README.md) and exempt from the staleness check.
TABLES_STATIC <- c(
  "Table_S3_mmo-capacity-screen"         = "code/lib/mmo_capacity_screen_table_S3_2026-09-30.csv")

# Photo plates outside the numbered SI figures (none at present)
PHOTOS <- c()

# STALENESS. "Exists" is not "current": a generator that fails leaves its previous
# output in place. The reference is the marker run_all.R writes when it starts, so
# "stale" means "its generator did not produce it during the last pipeline run".
# Photographs and the committed NCBI screen table are exempt.
STALE_REF <- local({
  mk <- "outputs/.pipeline_run_started"
  if (file.exists(mk)) file.mtime(mk) else NA
})
stale_list <- character(0)

copy_set <- function(map, dest, check_stale = TRUE) {
  miss <- 0
  for (nm in names(map)) {
    src <- map[[nm]]
    if (file.exists(src)) {
      if (check_stale && !is.na(STALE_REF) && file.mtime(src) < STALE_REF)
        stale_list <<- c(stale_list, sprintf("%s  <-  %s  (%s)", nm, src,
                          format(file.mtime(src), "%Y-%m-%d %H:%M")))
      file.copy(src, file.path(dest, paste0(nm, ".", tools::file_ext(src))), overwrite = TRUE)
    }
    else { cat("  MISSING:", src, "->", nm, "\n"); miss <- miss + 1 }
  }
  miss
}
m1 <- copy_set(MAIN, mainD); m2 <- copy_set(SI, siD) + copy_set(SI_PHOTOS, siD, check_stale = FALSE)
m3 <- copy_set(PHOTOS, photoD, check_stale = FALSE); m4 <- copy_set(TABLES, tableD) + copy_set(TABLES_STATIC, tableD, check_stale = FALSE)
cat(sprintf("Assembled %d main + %d SI + %d photo + %d table (%d missing; regenerate with: Rscript code/make_figures.R).\n",
            length(MAIN)-m1, length(SI)+length(SI_PHOTOS)-m2, length(PHOTOS)-m3, length(TABLES)+length(TABLES_STATIC)-m4, m1+m2+m3+m4))
if (length(stale_list)) {
  cat(sprintf("\n  !! %d assembled file(s) PREDATE this pipeline run (started %s) -- their\n",
              length(stale_list), format(STALE_REF, "%Y-%m-%d %H:%M")))
  cat("     generator did not run or failed, so these are LEFTOVERS:\n")
  for (s in stale_list) cat("       ", s, "\n")
  cat("     Fix the generator; do not ship these.\n")
} else if (!is.na(STALE_REF)) {
  cat("  all assembled figures/tables were produced by this pipeline run\n")
} else {
  cat("  (staleness not checked: no outputs/.pipeline_run_started marker;\n   run via code/run_all.R to enable the check)\n")
}

# MANIFEST.md — regenerated from the maps each run so it can never drift
man <- c("# Revised manuscript figure set (auto-generated by zz_assemble_figures.R)",
         "", "SI order = order of first citation in the manuscript text.", "")
sect <- function(title, map, dest) c(paste0("## ", title), "", "| Figure | Source |", "|---|---|",
  vapply(names(map), function(nm) sprintf("| %s | `%s`%s |", nm, map[[nm]],
    if (!file.exists(map[[nm]])) " (missing — run code/make_figures.R)" else ""), character(1)), "")
man <- c(man, sect("Main", MAIN, mainD), sect("Supplementary", c(SI, SI_PHOTOS), siD),
         sect("Photo plates", PHOTOS, photoD), sect("Tables (data CSV + rendered)", c(TABLES, TABLES_STATIC), tableD))
writeLines(man, "outputs/figures/MANIFEST.md")
cat("Wrote outputs/figures/MANIFEST.md\n")
