# ==============================================================================
# species_levels.R   --   SOURCED, never run alone
# ------------------------------------------------------------------------------
# THE single definition of how a raw species name maps onto a TreeRF model level.
#
# WHY THIS FILE EXISTS. The training side builds a taxonomic ladder
# (02_rf_models.R, "SPECIES LEVELS"): a species with >= 10 records keeps its own
# level, otherwise it falls back to GEN_<genus> if that genus was measured, and
# only otherwise to SPECIES_OTHER. The prediction side did not replicate it. Four
# scripts each carried
#
#     ifelse(species %in% trained, species, "SPECIES_OTHER")
#
# and `trained` holds the POOLED level names -- "GEN_Kalmia", "GEN_Quercus", ...
# -- while the inventory holds raw binomials. So "Kalmia latifolia" matched
# nothing and fell through to SPECIES_OTHER, taking 1,905 mountain laurel stems
# and the 19 Quercus velutina with it. Only 6 stems (Acer pensylvanicum) plus 41
# unidentified genuinely belonged there. tree_flux_predictions.csv contained zero
# GEN_* rows, which is the proof: the ladder never reached the inventory.
#
# It was not a harmless relabelling. SPECIES_OTHER has the second-highest observed
# mean of the 18 levels (0.2494 nmol m-2 s-1, calibration ratio 1.228) against
# GEN_Kalmia's 0.0350 and ratio 0.468, so 1,924 stems were predicted from a bucket
# whose mean is 7x their own. It carried 11.81% of the tree flux on 3.81% of the
# band area, and the fix moves the measured band term down ~7%.
#
# It also made 10_rf_species_pooling.R inert: the scheme its cross-validation
# selected was not the scheme the budget used.
#
# The mapping is fully determined by `trained` -- no taxonomy table is needed,
# because a level named GEN_<genus> is itself the evidence that <genus> was
# measured. Genus is the first whitespace-delimited token, exactly as on the
# training side.
# ==============================================================================

# Map raw species names onto TreeRF levels, reproducing the training ladder.
#
#   species  character vector of raw names (binomials, possibly NA)
#   trained  the model's own level vocabulary, i.e.
#            sort(unique(as.character(tree_train_complete$species_clean)))
#
# Returns a character vector: an element of `trained`, or NO_LEVEL for a stem with
# no own or genus level (unidentified, or of a genus never measured). NO_LEVEL is not
# a model level; such stems are predicted with predict_species_marginal() below.
NO_LEVEL <- "SPECIES_OTHER"
species_to_model_level <- function(species, trained) {
  stopifnot(is.character(trained), length(trained) > 0)

  sp <- as.character(species)

  # Unidentified, exactly as 02_rf_models.R defines it: missing, "unknown*", or
  # a single token (a genus-only or garbage entry is not a binomial).
  unident <- is.na(sp) | grepl("^unknown", sp, ignore.case = TRUE) |
             !grepl(" ", trimws(sp))

  gen_level <- paste0("GEN_", sub(" .*", "", trimws(sp)))

  out <- ifelse(unident,               NO_LEVEL,
         ifelse(sp %in% trained,        sp,
         ifelse(gen_level %in% trained, gen_level,
                                        NO_LEVEL)))

  # An own-species level must win over its genus level; the nesting above already
  # guarantees that, but assert it rather than trust the reading.
  stopifnot(all(out %in% c(trained, NO_LEVEL)))
  out
}

# STEMS WITHOUT A LEVEL (2026-10-16). An inventory stem whose species has no trained
# level -- unidentified, or a species of an unmeasured genus (Hamamelis, Acer
# pensylvanicum, Amelanchier, Vaccinium, Liriodendron, Quercus coccinea in this
# inventory) -- is predicted under every trained level in turn, with its own DBH,
# moisture and position, and the predictions are averaged. The weights are each
# level's share of the 0-2 m band area (pi x DBH x band; Kalmia 0.75 m) of the
# in-stand inventory stems that map to it, i.e. the stand's own woody-surface
# composition. There is no genus fallback for these stems. (Until 2026-10-16 they
# took a SPECIES_OTHER level whose only training record was tree 7809, an
# unidentified 2021 stem that trees.csv identifies as Kalmia latifolia.)

# Band-area share of each trained level among inventory stems mapped to it.
#   level      species_to_model_level() of the in-stand inventory stems
#   band_area  their 0-2 m band area (m2)
species_level_weights <- function(level, band_area, trained) {
  ok <- level %in% trained & is.finite(band_area)
  w <- tapply(band_area[ok], factor(level[ok], levels = trained), sum)
  w[is.na(w)] <- 0
  stopifnot(sum(w) > 0)
  w / sum(w)
}

# Predictions for stems with a level, and the weighted average over levels for
# stems with NO_LEVEL.
#   predict_fn  function(newdata) -> numeric predictions; newdata has a `species`
#               factor with levels names(weights) plus the columns of `newdata`
#   newdata     the other predictors, one row per stem (no species column)
#   level       model level per row (species_to_model_level())
#   weights     species_level_weights()
#   mult        optional per-level multiplier (the per-species calibration), applied
#               inside the average for NO_LEVEL rows only; rows with a level are
#               returned as predicted
predict_species_marginal <- function(predict_fn, newdata, level, weights, mult = NULL) {
  lv <- names(weights)
  stopifnot(length(level) == nrow(newdata), all(level %in% c(lv, NO_LEVEL)))
  out <- rep(NA_real_, nrow(newdata))
  own <- level != NO_LEVEL
  if (any(own))
    out[own] <- predict_fn(data.frame(species = factor(level[own], levels = lv), newdata[own, , drop = FALSE]))
  if (any(!own)) {
    nd <- newdata[!own, , drop = FALSE]
    m <- if (is.null(mult)) setNames(rep(1, length(lv)), lv) else
           setNames(ifelse(is.na(mult[lv]), 1, mult[lv]), lv)
    acc <- 0
    for (l in lv[weights > 0])
      acc <- acc + weights[[l]] * m[[l]] *
        predict_fn(data.frame(species = factor(rep(l, nrow(nd)), levels = lv), nd))
    out[!own] <- acc
  }
  out
}

# One-line audit of what the mapping did, so a silent fallback can never recur.
# Print this wherever the mapping is applied.
report_species_levels <- function(species, level, trained) {
  fell <- level == NO_LEVEL & !is.na(species) &
          !grepl("^unknown", as.character(species), ignore.case = TRUE) &
          grepl(" ", trimws(as.character(species)))
  cat(sprintf("species levels: %d stems -> %d of %d model levels | own %d | genus %d | no level (averaged over levels) %d\n",
              length(level), length(unique(level[level != NO_LEVEL])), length(trained),
              sum(!grepl("^GEN_", level) & level != NO_LEVEL),
              sum(grepl("^GEN_", level)), sum(level == NO_LEVEL)))
  if (any(fell)) {
    tb <- sort(table(as.character(species)[fell]), decreasing = TRUE)
    cat("  identified species with no own or genus level (genuinely unmeasured):\n")
    for (i in seq_along(tb))
      cat(sprintf("    %-26s %d\n", names(tb)[i], as.integer(tb[i])))
  }
  invisible(NULL)
}

# Species codes the random-forest workflow names (05_model/01_load_and_prep_data.R
# TAXONOMY, in this order; names and taxonomy from data/tables/species.csv). Census
# stems of any other code are left out of the model's inventory, and
# 01_import/08_inventory_build.R names only these codes.
MODEL_SPECIES <- c("ACRU", "ACSA", "ACPE", "BEAL", "BELE", "BEPA", "BEPO", "FAGR", "FRAM",
                   "PIST", "QURU", "QUAL", "QUVE", "TSCA", "CAOV", "CALA", "KALA", "PRSE", "SAAL")
