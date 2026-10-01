# ==============================================================================
# dead_stems.R -- which flux deployments were made on dead stems
# ------------------------------------------------------------------------------
# One definition, used by 02_rf_models.R (to keep dead stems out of training) and
# 03_export_canonical_tables.R (to say why they are not in it).
#
# WHY (2026-09-30). The stand inventory is live stems only: neither census records
# a dead status, so standing dead trees are outside the upscaling. Ten dead stems
# were measured -- 6 untagged snags in the monthly survey, 1 dead red oak in 2021,
# 3 trees marked dead in the 2023 survey (bark loss or wounding = "dead") -- an
# unplanned subset, half of the monthly ones in the wetland swamp, with
# inconsistent fluxes (monthly snags ~10x live on a typical day; the 2023 dead
# trees at or below zero). With no live/dead predictor possible for the inventory,
# training on them shifts predictions for live trees of the same species. Tested:
# a live/dead predictor changed nothing (the forest rarely splits on a 3.6% flag).
# They stay in every descriptive analysis.
#
# Deployments are identified by campaign and CH4 flux value, which is unique within
# a campaign (03_export_canonical_tables.R asserts this).
# ==============================================================================
dead_stem_flux <- function(dir = NULL) {
  if (is.null(dir)) for (d in c("data/processed/flux", "../../data/processed/flux"))
    if (dir.exists(d)) { dir <- d; break }
  rd <- function(f) read.csv(file.path(dir, f), check.names = FALSE, stringsAsFactors = FALSE)
  m <- rd("semirigid_tree_final_complete_dataset_with_untagged.csv")
  h <- rd("tree_flux_2021_multiheight.csv")
  y <- rd("tree_flux_2023_cross_species.csv")
  y_dead <- grepl("^dead$", trimws(y[["Bark Missing (1-3)"]] %||% y[[grep("^Bark", names(y))[1]]]), ignore.case = TRUE) |
            grepl("^dead$", trimws(y[[grep("^Wounding", names(y))[1]]]), ignore.case = TRUE)
  out <- rbind(
    data.frame(campaign = "monthly_2020_2021",  flux = m$CH4_best.flux.x[grepl("_dead$", m$`Plot Tag` %||% m$Plot.Tag)]),
    data.frame(campaign = "2021_multiheight",   flux = h$CH4_best.flux[grepl("dead", h$tree_id, ignore.case = TRUE)]),
    data.frame(campaign = "2023_cross_species", flux = y$CH4_best.flux[y_dead]))
  out[is.finite(out$flux), ]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# TRUE where a flux value (any campaign) belongs to a dead stem
is_dead_stem_flux <- function(flux, dead = dead_stem_flux())
  round(as.numeric(flux), 8) %in% round(dead$flux, 8)
