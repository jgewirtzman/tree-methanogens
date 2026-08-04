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
# It also made rf_species_pooling.R inert: the scheme its cross-validation
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
# Returns a character vector, every element guaranteed to be an element of
# `trained`, so it is always safe to wrap in factor(levels = trained).
species_to_model_level <- function(species, trained) {
  stopifnot(is.character(trained), length(trained) > 0)
  if (!"SPECIES_OTHER" %in% trained)
    stop("`trained` has no SPECIES_OTHER level -- wrong vocabulary passed in")

  sp <- as.character(species)

  # Unidentified, exactly as 02_rf_models.R defines it: missing, "unknown*", or
  # a single token (a genus-only or garbage entry is not a binomial).
  unident <- is.na(sp) | grepl("^unknown", sp, ignore.case = TRUE) |
             !grepl(" ", trimws(sp))

  gen_level <- paste0("GEN_", sub(" .*", "", trimws(sp)))

  out <- ifelse(unident,               "SPECIES_OTHER",
         ifelse(sp %in% trained,        sp,
         ifelse(gen_level %in% trained, gen_level,
                                        "SPECIES_OTHER")))

  # An own-species level must win over its genus level; the nesting above already
  # guarantees that, but assert it rather than trust the reading.
  stopifnot(all(out %in% trained))
  out
}

# One-line audit of what the mapping did, so a silent fallback can never recur.
# Print this wherever the mapping is applied.
report_species_levels <- function(species, level, trained) {
  fell <- level == "SPECIES_OTHER" & !is.na(species) &
          !grepl("^unknown", as.character(species), ignore.case = TRUE) &
          grepl(" ", trimws(as.character(species)))
  cat(sprintf("species levels: %d stems -> %d of %d model levels | own %d | genus %d | SPECIES_OTHER %d\n",
              length(level), length(unique(level)), length(trained),
              sum(!grepl("^GEN_|^SPECIES_OTHER$", level)),
              sum(grepl("^GEN_", level)), sum(level == "SPECIES_OTHER")))
  if (any(fell)) {
    tb <- sort(table(as.character(species)[fell]), decreasing = TRUE)
    cat("  identified species with NO own or genus level (genuinely unmeasured):\n")
    for (i in seq_along(tb))
      cat(sprintf("    %-26s %d\n", names(tb)[i], as.integer(tb[i])))
  }
  invisible(NULL)
}
