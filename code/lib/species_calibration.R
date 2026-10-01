# ==============================================================================
# species_calibration.R -- per-species stand-total calibration of the TreeRF
# ------------------------------------------------------------------------------
# One definition, used by predict_tree_flux_current.R (applies it) and
# rf_calibration_sensitivity.R (removes and re-applies it), so the two cannot drift.
#
# ratio = mean observed / mean OUT-OF-BAG predicted flux, per species level; it
# corrects the forest's shrinkage toward the mean for a SUM (see the long note in
# predict_tree_flux_current.R).
#
# MINIMUM SAMPLE (2026-09-30): a level with fewer than MIN_CAL_N training
# measurements, or a non-positive ratio, gets ratio 1 (no correction). After the
# untagged stems were given their species, SPECIES_OTHER rested on ONE measurement
# with a slightly negative flux and a ratio of -0.30, which turned 47 inventory
# stems from emitting to absorbing. A mean of 1-4 measurements cannot calibrate a
# level; this is the "low n -> 1" case rf_calibration_sensitivity.R already reports.
# ==============================================================================
MIN_CAL_N <- 5

species_calibration <- function(sp, obs, oob, min_n = MIN_CAL_N) {
  data.frame(sp = as.character(sp), obs = obs, oob = oob) |>
    subset(is.finite(obs) & is.finite(oob)) |>
    (\(x) do.call(rbind, lapply(split(x, x$sp), function(g) {
      om <- mean(g$oob); raw <- if (om > 0) mean(g$obs) / om else 1
      data.frame(sp = g$sp[1], n = nrow(g), obs_mean = mean(g$obs), oob_mean = om,
                 ratio_raw = raw, ratio = if (nrow(g) >= min_n && raw > 0) raw else 1)
    })))() |>
    (\(x) { rownames(x) <- NULL; x })()
}
