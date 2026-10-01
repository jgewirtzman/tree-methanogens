# ==============================================================================
# migrate_campaign_filenames.R -- bring an older data/ drop-in up to the
# campaign-named flux files (2026-09-30)
# ------------------------------------------------------------------------------
# The three tree goFlux scripts used to write shared filenames; they now write
# per-campaign names (see assemble_campaign_flux.R). A data/ folder restored from
# an archive made before that change still holds the old names. This script renames
# each old file to the campaign its CONTENT belongs to -- read from the deployment
# dates in UniqueID, or the table's own columns -- never from the old name alone,
# because the old names held whichever campaign ran last.
#
# Idempotent: a file is renamed only if its target does not exist; otherwise it is
# reported and left alone. Run from the repo root:
#   Rscript code/tools/migrate_campaign_filenames.R
# then, if tree_flux_2021_multiheight.csv is absent:
#   Rscript code/02_flux/assemble_campaign_flux.R 2021_multiheight
# ==============================================================================
fd <- Sys.getenv("FLUX_DIR", "data/processed/flux")   # override for testing
campaign_of <- function(f) {
  x <- read.csv(file.path(fd, f), check.names = FALSE, nrows = 5000, stringsAsFactors = FALSE)
  if ("Tree Tag" %in% names(x)) return("2023_cross_species")
  if (all(c("tree_id", "measurement_height") %in% names(x))) return("2021_multiheight")
  id <- if ("UniqueID" %in% names(x)) x$UniqueID else character(0)
  ym <- unique(substr(id[grepl("^[0-9]{8}", id)], 1, 6))
  if (!length(ym)) return(NA_character_)
  if (all(ym %in% c("202107", "202108"))) return("2021_multiheight")
  if (all(substr(ym, 1, 4) == "2023")) return("2023_cross_species")
  if (all(ym >= "202006" & ym <= "202105")) return("semirigid_tree")
  NA_character_
}
OLD <- c("methanogen_tree_flux_complete_dataset.csv", "CH4_best_flux_lgr_results.csv",
         "CO2_best_flux_lgr_results.csv", "CH4_flux_lgr_results.csv", "CO2_flux_lgr_results.csv",
         "lgr_manual_identification_results.csv", "lgr_manual_identification_results_final.csv",
         "lgr_manual_identification_summary_final.csv")
target <- function(f, cmp) {
  if (f == "methanogen_tree_flux_complete_dataset.csv") return(paste0("tree_flux_", cmp, ".csv"))
  if (grepl("^lgr_manual_identification_results_final", f)) return(paste0("lgr_manual_identification_", cmp, "_final.csv"))
  if (grepl("^lgr_manual_identification_summary_final", f)) return(paste0("lgr_manual_identification_summary_", cmp, "_final.csv"))
  if (grepl("^lgr_manual_identification_results", f)) return(paste0("lgr_manual_identification_", cmp, ".csv"))
  sub("\\.csv$", paste0("_", cmp, ".csv"), f)
}
for (f in OLD[file.exists(file.path(fd, OLD))]) {
  cmp <- if (grepl("summary", f)) "2021_multiheight" else campaign_of(f)
  if (is.na(cmp)) { cat(sprintf("  ?? %-48s campaign not identifiable; left alone\n", f)); next }
  to <- target(f, cmp)
  if (file.exists(file.path(fd, to))) { cat(sprintf("  == %-48s %s already exists; left alone\n", f, to)); next }
  file.rename(file.path(fd, f), file.path(fd, to)); cat(sprintf("  -> %-48s %s\n", f, to))
}
cat("done\n")
