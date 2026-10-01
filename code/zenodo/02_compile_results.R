# ==============================================================================
# 02_compile_results.R -- the paper's computed results, archived in data/compiled/results/
# ------------------------------------------------------------------------------
# Runs at the end of the pipeline (stage F). Copies, without recomputing, the result
# tables every stand-level and summary number in the paper is read from, so each can
# be checked from the archive without rerunning the model. check_consistency.R
# asserts their agreement with one another.
# ==============================================================================
out <- "data/compiled/results"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(out, full.names = TRUE))
RESULTS <- c(
  # stand budget and scaling (Results 10, Fig 9, Figs S26-S27)
  "canonical_budget.csv"         = "outputs/data/canonical_budget.csv",
  "canonical_monthly.csv"        = "outputs/data/canonical_monthly.csv",
  "scaling_full_grid.csv"        = "outputs/data/scaling_full_grid.csv",
  "scaling_headline.csv"         = "outputs/data/scaling_headline.csv",
  "scaling_slope_diagnostics.csv"= "outputs/data/scaling_slope_diagnostics.csv",
  "wai_bottomup.csv"             = "outputs/data/wai_bottomup.csv",
  # model skill (Results 10, Methods S4, S7)
  "rf_grouped_cv.csv"            = "outputs/data/rf_grouped_cv.csv",
  "gene_rf_cv.csv"               = "outputs/data/gene_rf_cv.csv",
  # isotopes (Results 7, Fig 6d)
  "isotopes_summary.csv"         = "outputs/data/ISOTOPES_summary.csv",
  # tables
  "table1_campaign_counts.csv"   = "outputs/data/campaign_counts.csv",
  "tableS2_known_putative_taxa.csv" = "outputs/data/known_putative_taxa_table.csv",
  "tableS3_ddpcr_16s_concordance.csv" = "outputs/data/tbl_ddpcr_16s_concordance.csv",
  "tableS4_dbh_by_species_campaign.csv" = "outputs/data/dbh_by_species_campaign.csv",
  "tableS5_mmo_capacity_screen.csv" = "outputs/data/mmo_capacity_screen.csv",
  # PICRUSt2 pathway associations (Fig 6c, Figs S11-S13)
  "picrust_pathway_associations.csv" = "data/processed/molecular/picrust/pathway_associations_combined.csv",
  "picrust_pathway_associations_pmoa.csv" = "data/processed/molecular/picrust/pathway_associations_pmoa.csv")
miss <- RESULTS[!file.exists(RESULTS)]
if (length(miss)) stop("missing result file(s) -- run the pipeline first:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)
for (nm in names(RESULTS)) file.copy(RESULTS[[nm]], file.path(out, nm), overwrite = TRUE)
cat(sprintf("%d result tables in %s/\n", length(RESULTS), out))
