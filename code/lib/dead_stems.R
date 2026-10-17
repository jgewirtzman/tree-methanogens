# ==============================================================================
# dead_stems.R -- which flux deployments were made on dead stems
# ------------------------------------------------------------------------------
# The definition lives in code/tables/04_trees.R (trees.csv dead, dead_from,
# dead_evidence, dead_basis): a tree is dead from the first date any record marks it
# dead, if it is marked dead on more than one date or is not measured again after
# that date. code/tables/10_flux_stem.R turns it into the per-deployment flag
# flux_stem.csv dead_stem (dead tree, deployment on or after dead_from). This file
# only reads that flag; it applies no rule of its own.
#
# WHY dead stems are kept out of training (2026-09-30). The stand inventory is live
# stems only: neither census records a dead status, so standing dead trees are outside
# the upscaling. Dead stems were an unplanned subset (half of the monthly ones in the
# wetland swamp) with inconsistent fluxes; with no live/dead predictor possible for the
# inventory, training on them shifts predictions for live trees of the same species.
# They stay in every descriptive analysis.
# ==============================================================================

# measurement_id of every dead-stem deployment (data/tables/flux_stem.csv)
dead_stem_ids <- function() {
  f <- Filter(file.exists, c("data/tables/flux_stem.csv", "../../data/tables/flux_stem.csv"))[1]
  if (is.na(f)) stop("data/tables/flux_stem.csv not found", call. = FALSE)
  d <- utils::read.csv(f, colClasses = "character", na.strings = "")
  stopifnot(all(d$dead_stem %in% c("TRUE", "FALSE")))
  d$measurement_id[d$dead_stem == "TRUE"]
}
