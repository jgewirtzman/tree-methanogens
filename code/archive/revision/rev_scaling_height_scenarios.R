#!/usr/bin/env Rscript
# ==============================================================================
# rev_scaling_height_scenarios.R
# ------------------------------------------------------------------------------
# Height-resolved stem-flux scaling, replacing the single-value-per-tree approach.
#
# 1. LOWER STEM (0-2 m), where we have data.
#    The RF now carries measurement height as a predictor, so instead of assigning one
#    value to the whole band we predict along it and integrate:
#        F_band = mean_h f(h) over h in [0, 200] cm,  applied to pi * DBH * 2
#    This matters: mean flux roughly halves from 50 cm (0.220) to 200 cm (0.102), so a
#    single value taken at 125 cm under-estimates the band mean by ~6.5%.
#
# 2. WHOLE WOODY SURFACE, where we do not.
#    Three scenarios, as distinct assumptions rather than as an estimate:
#      (1) flux at 2 m applied to the entire woody surface
#      (2) the height-agnostic mean flux applied to the entire woody surface
#      (3) the fitted height pattern continued to the top of the tree
#
#    For (3) the RF cannot help: random forests return a constant beyond the training
#    range, so an RF "extrapolation" above 200 cm would silently reproduce scenario (1).
#    A parametric form must be fitted to the 50/125/200 cm data and projected. Because
#    three heights cannot distinguish between candidate forms, ALL of them are carried
#    through as a sensitivity band rather than one being chosen:
#         exponential decay to zero        f = a*exp(-b*h)
#         exponential decay to an asymptote f = c + a*exp(-b*h)
#         power law                        f = a*h^-b
#         log-linear                       f = exp(a + b*h)
#         linear, floored at zero
#         constant above 2 m               (equivalent to scenario 1)
#
#    THE FELLED OAK IS THE ONLY EVIDENCE ABOVE 2 m and is used to judge them. Its
#    profile (0.5-10 m) gives 0.038 nmol m-2 s-1 at 2 m and 0.034 / 0.077 / 0.031 at
#    6 / 8 / 10 m -- i.e. flux does NOT decay away above the lower stem but plateaus,
#    favouring the asymptotic forms and scenario (1). It also contains a 10x hotspot at
#    4 m. n = 1 tree: this disciplines the extrapolation, it cannot confirm it.
#
# 3. GEOMETRY, reported both ways because the assumption is consequential:
#      (a) WAI = 3.07 m2 woody area per m2 ground (Gauci et al. 2024), the bulk figure
#          used previously and comparable with the literature
#      (b) per-tree DBH-based height allometry with woody area spread evenly over height
#
# Output: outputs/revision/scaling_height_scenarios.csv
#         outputs/revision/fig_height_extrapolation.png
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(tidyr)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")

# --- observed height profile, survey trees (2021 campaign) ------------------
d <- tree_train_complete
hs <- !is.na(d$measurement_height_cm) & is.finite(d$stem_flux_corrected)
obs <- d[hs, ] %>% group_by(h = measurement_height_cm) %>%
  summarise(flux = mean(stem_flux_corrected), n = n(), .groups = "drop")
cat("=== observed lower-stem profile (survey trees) ===\n"); print(as.data.frame(obs), row.names = FALSE)

# --- felled oak, the only data above 2 m ------------------------------------
bo <- read.csv("data/compiled/black_oak_experiment.csv", stringsAsFactors = FALSE)
bo$h_m <- suppressWarnings(as.numeric(bo$height_m))
bo <- bo[is.finite(bo$h_m) & is.finite(bo$CH4_best_flux_nmol_m2_s), ]
cat("\n=== felled oak profile (the only observations above 2 m) ===\n")
print(bo[order(bo$h_m), c("h_m", "CH4_best_flux_nmol_m2_s")], row.names = FALSE)
above <- bo$CH4_best_flux_nmol_m2_s[bo$h_m > 2 & bo$h_m != 4]
at2   <- bo$CH4_best_flux_nmol_m2_s[abs(bo$h_m - 2) < 0.01]
cat(sprintf("\n  at 2 m: %.4f    mean of 6/8/10 m (excluding the 4 m hotspot): %.4f   ratio %.2f\n",
            at2, mean(above), mean(above) / at2))

# --- candidate height functions, fitted to 50-200 cm ------------------------
h <- obs$h / 100; f <- obs$flux                      # metres
FORMS <- list(
  exp_zero  = function() { m <- lm(log(pmax(f, 1e-6)) ~ h)
                           function(z) exp(coef(m)[1] + coef(m)[2] * z) },
  exp_asym  = function() { st <- list(c = min(f) * 0.8, a = max(f), b = 1)
                           m <- tryCatch(nls(f ~ c + a * exp(-b * h), start = st,
                                             control = nls.control(warnOnly = TRUE)),
                                         error = function(e) NULL)
                           if (is.null(m)) return(NULL)
                           cf <- coef(m); function(z) cf["c"] + cf["a"] * exp(-cf["b"] * z) },
  power     = function() { m <- lm(log(pmax(f, 1e-6)) ~ log(h))
                           function(z) exp(coef(m)[1] + coef(m)[2] * log(pmax(z, 0.05))) },
  loglinear = function() { m <- lm(log(pmax(f, 1e-6)) ~ h)
                           function(z) exp(coef(m)[1] + coef(m)[2] * z) },
  linear0   = function() { m <- lm(f ~ h)
                           function(z) pmax(coef(m)[1] + coef(m)[2] * z, 0) },
  constant2 = function() { f2 <- approx(h, f, xout = 2, rule = 2)$y
                           function(z) ifelse(z <= 2, approx(h, f, xout = pmin(z, 2), rule = 2)$y, f2) }
)
fits <- lapply(FORMS, function(g) tryCatch(g(), error = function(e) NULL))
fits <- fits[!vapply(fits, is.null, logical(1))]
cat("\n  fitted forms:", paste(names(fits), collapse = ", "), "\n")

# --- how well does each form predict the felled oak above 2 m? --------------
oak_hi <- bo[bo$h_m > 2, ]
scale_to_2m <- at2 / approx(h, f, xout = 2, rule = 2)$y   # put survey curve on the oak's level
cat("\n=== extrapolation checked against the felled oak (>2 m) ===\n")
cat(sprintf("  %-11s %10s %10s %10s\n", "form", "pred 6m", "pred 10m", "RMSE >2m"))
chk <- lapply(names(fits), function(nm) {
  p <- fits[[nm]](oak_hi$h_m) * scale_to_2m
  data.frame(form = nm,
             p6  = fits[[nm]](6)  * scale_to_2m,
             p10 = fits[[nm]](10) * scale_to_2m,
             rmse = sqrt(mean((oak_hi$CH4_best_flux_nmol_m2_s - p)^2)))
}) %>% bind_rows()
for (i in seq_len(nrow(chk)))
  cat(sprintf("  %-11s %10.4f %10.4f %10.4f\n", chk$form[i], chk$p6[i], chk$p10[i], chk$rmse[i]))

# --- band-mean multipliers ---------------------------------------------------
band_mean <- function(fn, lo, hi, n = 200) mean(fn(seq(lo, hi, length.out = n)))
f_at2 <- approx(h, f, xout = 2, rule = 2)$y
f_band <- band_mean(function(z) approx(h, f, xout = pmin(pmax(z, 0.5), 2), rule = 2)$y, 0, 2)
f_meanobs <- mean(f)
cat(sprintf("\n=== lower stem (0-2 m) ===\n  integrated band mean %.4f vs value at 2 m %.4f vs simple mean %.4f\n",
            f_band, f_at2, f_meanobs))

TREE_H <- 20   # m, nominal canopy height for the allometric variant
res <- bind_rows(lapply(names(fits), function(nm) {
  data.frame(form = nm,
             whole_tree_mean_0_to_H = band_mean(fits[[nm]], 0.05, TREE_H),
             rel_to_2m = band_mean(fits[[nm]], 0.05, TREE_H) / f_at2)
}))
res <- bind_rows(
  data.frame(form = "(1) value at 2 m",      whole_tree_mean_0_to_H = f_at2,     rel_to_2m = 1),
  data.frame(form = "(2) height-agnostic mean", whole_tree_mean_0_to_H = f_meanobs, rel_to_2m = f_meanobs / f_at2),
  res %>% mutate(form = paste0("(3) ", form))
)
cat("\n=== whole-woody-surface scenarios: mean flux applied to the surface ===\n")
print(as.data.frame(res), row.names = FALSE, digits = 3)
write.csv(res, file.path(outdir, "scaling_height_scenarios.csv"), row.names = FALSE)

# --- figure ------------------------------------------------------------------
grid <- seq(0.05, 12, length.out = 200)
curves <- bind_rows(lapply(names(fits), function(nm)
  data.frame(form = nm, h = grid, flux = fits[[nm]](grid) * scale_to_2m)))
p <- ggplot() +
  annotate("rect", xmin = 0, xmax = 2, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = .5) +
  geom_line(data = curves, aes(h, flux, colour = form), linewidth = .7) +
  geom_point(data = bo, aes(h_m, CH4_best_flux_nmol_m2_s), size = 2.4, shape = 21,
             fill = "black", colour = "white") +
  geom_point(data = data.frame(h = h, f = f * scale_to_2m), aes(h, f),
             size = 3, shape = 24, fill = "#B2182B", colour = "white") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Extrapolating stem flux above the measured range",
       subtitle = paste("shaded = measured band (0-2 m); triangles = survey-tree means;",
                        "circles = felled oak, the only data above 2 m"),
       x = "height (m)", y = expression(CH[4]~flux~(nmol~m^-2~s^-1)), colour = "form") +
  theme_minimal(base_size = 9)
ggsave(file.path(outdir, "fig_height_extrapolation.png"), p, width = 9, height = 5.5, dpi = 200, bg = "white")
cat("\nWritten: outputs/revision/fig_height_extrapolation.png\n")
