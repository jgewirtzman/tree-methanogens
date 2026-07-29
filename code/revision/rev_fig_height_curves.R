#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_height_curves.R
# ------------------------------------------------------------------------------
# Species-specific height-flux curves from the RF, compared against the linear
# mixed-model slopes the manuscript currently reports (Fig 2).
#
# WHY THIS ADDS SOMETHING
#   The LMM fits ONE slope per species, so a response that declines then flattens --
#   or peaks mid-stem, as the felled oak does at 4-6 m -- is forced into a straight
#   line. With measurement height as an RF feature the shape is free to be non-linear
#   and to differ by species, and height x moisture interactions are tested without
#   being pre-specified.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   Extrapolate above 200 cm. Random forests cannot extrapolate: beyond the training
#   range they return a constant, so predictions above 2 m would flat-line as an
#   ARTEFACT of the method and could be misread as "flux is constant with height" --
#   precisely the claim the whole-woody-surface scenario turns on. The felled oak
#   profile (to 10 m) remains the only evidence above 2 m.
#
# Output: outputs/revision/fig_height_curves.png
#         outputs/revision/height_curves_rf_vs_lmm.csv
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(lme4)})
set.seed(42)

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

X <- as.data.frame(X_tree); d <- tree_train_complete
stopifnot("height_cm" %in% names(X))

# The height campaign is the only source of within-tree height contrast
hs <- !is.na(d$measurement_height_cm)
cat("height-campaign measurements:", sum(hs), "on", length(unique(d$tree_id[hs])), "trees\n")

sp_n <- sort(table(d$species_clean[hs]), decreasing = TRUE)
species <- names(sp_n)[sp_n >= 20 & names(sp_n) != "SPECIES_OTHER"]
cat("species with >= 20 height measurements:", paste(species, collapse = ", "), "\n\n")

# ---- RF partial dependence on height, within each species ------------------
grid_h <- seq(50, 200, length.out = 16)
pd <- bind_rows(lapply(species, function(sp) {
  rows <- which(d$species_clean == sp & hs)
  xs <- X[rows, , drop = FALSE]
  bind_rows(lapply(grid_h, function(h) {
    xs$height_cm <- h
    # TreeRF predicts on the MEASURED scale (nmol m-2 s-1); the asinh transform was
    # removed. This was mean(sinh(...)), a leftover back-transform that inflated the
    # RF curve while the LMM curve below is a genuine asinh fit -- so the panel's
    # whole point, comparing the two, was comparing different scales.
    data.frame(species = sp, height_cm = h,
               flux = mean(predict(TreeRF, xs)$predictions, na.rm = TRUE))
  }))
}))

# ---- LMM slope per species, for comparison ---------------------------------
lmm_rows <- bind_rows(lapply(species, function(sp) {
  dd <- d[d$species_clean == sp & hs, ]
  dd$z <- asinh(dd$stem_flux_corrected)
  m <- tryCatch(lmer(z ~ measurement_height_cm + (1 | tree_id), data = dd),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- fixef(m)
  data.frame(species = sp, intercept = b[1], slope_per_cm = b[2],
             slope_per_m = b[2] * 100, n = nrow(dd))
}))
cat("=== LMM height slopes (asinh flux per metre) ===\n")
print(lmm_rows[, c("species", "n", "slope_per_m")], row.names = FALSE, digits = 3)

lmm_pred <- bind_rows(lapply(seq_len(nrow(lmm_rows)), function(i)
  data.frame(species = lmm_rows$species[i], height_cm = grid_h,
             flux = sinh(lmm_rows$intercept[i] + lmm_rows$slope_per_cm[i] * grid_h))))

obs <- d[hs & d$species_clean %in% species,
         c("species_clean", "measurement_height_cm", "stem_flux_corrected")]
names(obs) <- c("species", "height_cm", "flux")

p <- ggplot() +
  geom_jitter(data = obs, aes(height_cm, flux), width = 6, alpha = 0.25, size = 0.8,
              colour = "grey45") +
  geom_line(data = lmm_pred, aes(height_cm, flux, colour = "linear mixed model"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(data = pd, aes(height_cm, flux, colour = "random forest"), linewidth = 0.9) +
  facet_wrap(~ species, scales = "free_y") +
  scale_colour_manual(values = c("random forest" = "#B2182B",
                                 "linear mixed model" = "#2166AC"), name = NULL) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Stem CH4 flux vs measurement height, by species",
       subtitle = paste("RF partial dependence (solid) vs LMM slope (dashed);",
                        "points are observations. Range limited to 50-200 cm --",
                        "\nRFs cannot extrapolate, so no curve is drawn above the measured range."),
       x = "measurement height (cm)",
       y = expression(CH[4]~flux~(nmol~m^-2~s^-1)~", pseudo-log")) +
  theme_minimal(base_size = 9) + theme(legend.position = "top")

ggsave(file.path(outdir, "fig_height_curves.png"), p, width = 10, height = 7, dpi = 200, bg = "white")
write.csv(pd %>% rename(rf_flux = flux) %>%
            left_join(lmm_pred %>% rename(lmm_flux = flux), by = c("species", "height_cm")),
          file.path(outdir, "height_curves_rf_vs_lmm.csv"), row.names = FALSE)

# ---- is the RF shape actually non-linear? ----------------------------------
cat("\n=== curvature: does the RF depart from the LMM straight line? ===\n")
cat(sprintf("  %-26s %10s %10s %10s\n", "species", "RF 50cm", "RF 200cm", "% change"))
for (sp in species) {
  q <- pd[pd$species == sp, ]
  a <- q$flux[which.min(q$height_cm)]; b <- q$flux[which.max(q$height_cm)]
  cat(sprintf("  %-26s %10.4f %10.4f %9.0f%%\n", sp, a, b, 100 * (b / a - 1)))
}
cat("\nWritten: outputs/revision/fig_height_curves.png\n")
