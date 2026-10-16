# ==============================================================================
# dead_stems.R -- which flux deployments were made on dead stems
# ------------------------------------------------------------------------------
# The flag lives in data/tables/flux_stem.csv (dead_stem, by measurement_id; built by
# code/tables/10_flux_stem.R from the rules below). dead_stem_ids() is what
# 02_rf_models.R reads to keep dead stems out of training; dead_stem_flux() applies
# the same rules to the processed flux files, for the tree-level audit in
# code/tables/04_trees.R.
#
# WHY (2026-09-30). The stand inventory is live stems only: neither census records
# a dead status, so standing dead trees are outside the upscaling. Twelve dead stems
# were measured -- 6 untagged snags and 2 tagged stems in the monthly survey, 1 dead red oak in 2021,
# 3 trees marked dead in the 2023 survey (bark loss or wounding = "dead") -- an
# unplanned subset, half of the monthly ones in the wetland swamp, with
# inconsistent fluxes (monthly snags ~10x live on a typical day; the 2023 dead
# trees at or below zero). With no live/dead predictor possible for the inventory,
# training on them shifts predictions for live trees of the same species. Tested:
# a live/dead predictor changed nothing (the forest rarely splits on a 3.6% flag).
# They stay in every descriptive analysis.
#
# 2026-10-15: monthly stems flagged in the sheet's `dead` column are included
# (tagged red oak 4156 and hemlock 5091, also noted dead in the field notes and,
# for 5091, in 2023): 12 stems, 59 deployments.
# ==============================================================================
dead_stem_flux <- function(dir = NULL) {
  source(Filter(file.exists, c("code/lib/monthly_cols.R", "../lib/monthly_cols.R", "../../lib/monthly_cols.R"))[1])
  if (is.null(dir)) for (d in c("data/processed/flux", "../../data/processed/flux"))
    if (dir.exists(d)) { dir <- d; break }
  rd <- function(f) read.csv(file.path(dir, f), check.names = FALSE, stringsAsFactors = FALSE)
  m <- monthly_plain(rd("semirigid_tree_final_complete_dataset_with_untagged.csv"))
  h <- rd("tree_flux_2021_multiheight.csv")
  y <- rd("tree_flux_2023_cross_species.csv")
  y_dead <- grepl("^dead$", trimws(y[["Bark Missing (1-3)"]] %||% y[[grep("^Bark", names(y))[1]]]), ignore.case = TRUE) |
            grepl("^dead$", trimws(y[[grep("^Wounding", names(y))[1]]]), ignore.case = TRUE)
  # monthly: untagged snags (label ends "_dead") and tagged stems the sheet records
  # as dead (`dead` column: 4156 and 5091)
  m_dead <- grepl("_dead$", m$`Plot Tag` %||% m$Plot.Tag) |
            (!is.null(m$dead) & toupper(as.character(m$dead)) %in% "TRUE")
  out <- rbind(
    data.frame(campaign = "monthly_2020_2021",  flux = m$CH4_best.flux[m_dead]),
    data.frame(campaign = "2021_multiheight",   flux = h$CH4_best.flux[grepl("dead", h$tree_id, ignore.case = TRUE)]),
    data.frame(campaign = "2023_cross_species", flux = y$CH4_best.flux[y_dead]))
  out[is.finite(out$flux), ]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# measurement_id of every dead-stem deployment (data/tables/flux_stem.csv)
dead_stem_ids <- function() {
  f <- Filter(file.exists, c("data/tables/flux_stem.csv", "../../data/tables/flux_stem.csv"))[1]
  if (is.na(f)) stop("data/tables/flux_stem.csv not found", call. = FALSE)
  d <- utils::read.csv(f, colClasses = "character", na.strings = "")
  stopifnot(all(d$dead_stem %in% c("TRUE", "FALSE")))
  d$measurement_id[d$dead_stem == "TRUE"]
}
