#!/usr/bin/env Rscript
# ==============================================================================
# 06_moisture_elevation_check.R -- an independent check on the moisture surface
# ------------------------------------------------------------------------------
# Stage C (non-fatal check), after 04_moisture_surface.R.
#
# Elevation never enters the moisture interpolation, so topography is a free
# validation: water runs downhill, and a moisture surface that is capturing real
# structure rather than interpolation artifact should be wetter in low ground.
# If the two are uncorrelated, the surface is not describing the plot's hydrology
# and should be treated as decoration.
#
# Elevation: the GPS elevations recorded with the moisture survey points, the monthly-
# survey plots, the wells, the monthly-survey trees (GPS outliers, > 2 SD from the
# trees' median position, left out) and the stream flags
# (data/raw/inventory/spatial_data/).
#
# Two tests, in increasing dependence on interpolation:
#   1. AT THE SURVEY POINTS. Some rows of the elevation archive were recorded at
#      the moisture survey locations themselves, so moisture and elevation can be
#      compared with no interpolation at all. This is the strongest form.
#   2. SURFACE AGAINST SURFACE. Elevation is interpolated onto the same 2 m stand
#      grid with the same thin-plate spline, and the two z-scored fields are
#      compared cell by cell. Weaker, since both are smoothed, but it tests the
#      thing the budget actually consumes.
#
# Elevation is not a model predictor: elevation and distance to river are static per
# collar (28 values across 266 observations) and act as collar identity, raising OOB
# R2 while lowering grouped-CV R2 (05_model/02_rf_models.R, soil feature set). Using
# topography to VALIDATE an independently built surface is sound; using it to
# fit one here is the overfitting route already rejected.
#
# Output: outputs/audit/moisture_elevation_check.txt
#         outputs/figures/generated/fig_moisture_elevation.png
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(fields); library(ggplot2)
                                library(patchwork)})
set.seed(42)
source("code/lib/geometry.R")

SPATIAL <- "data/raw/inventory/spatial_data"
elev_of <- function(d, source) data.frame(Longitude = d$Longitude, Latitude = d$Latitude,
                                          Elevation = d$Elevation, Source = source, stringsAsFactors = FALSE)
TREES <- read.csv(file.path(SPATIAL, "YM_trees_measured.csv.csv"))
dc <- sqrt((TREES$Longitude - median(TREES$Longitude, na.rm = TRUE))^2 +
           (TREES$Latitude - median(TREES$Latitude, na.rm = TRUE))^2)
TREES <- TREES[dc <= mean(dc, na.rm = TRUE) + 2 * sd(dc, na.rm = TRUE), ]
E <- rbind(elev_of(read.csv(file.path(SPATIAL, "soil_moisture_20201216.csv")), "Soil_Moisture"),
           elev_of(read.csv(file.path(SPATIAL, "plots.csv")), "Research_Plots"),
           elev_of(read.csv(file.path(SPATIAL, "wells.csv")), "Wells"),
           elev_of(TREES, "Trees"),
           elev_of(readxl::read_excel(file.path(SPATIAL, "River.xlsx")), "River"))
E <- E[is.finite(E$Longitude) & is.finite(E$Latitude) & is.finite(E$Elevation), ]
M <- read.csv("data/tables/moisture_survey.csv", stringsAsFactors = FALSE)   # December 2020 survey
M <- data.frame(lon = M$longitude, lat = M$latitude, vwc = M$soil_vwc_pct)

cat(sprintf("elevation points %d (sources: %s)\n", nrow(E),
            paste(names(sort(table(E$Source), decreasing = TRUE)), collapse = ", ")))
cat(sprintf("elevation range %.1f - %.1f (units as archived)\n\n",
            min(E$Elevation), max(E$Elevation)))

# ---- 1. at the survey points, no interpolation ------------------------------
key <- function(a, b) paste(round(a, 6), round(b, 6))
E$k <- key(E$Longitude, E$Latitude); M$k <- key(M$lon, M$lat)
J <- inner_join(M, E %>% select(k, Elevation) %>% distinct(k, .keep_all = TRUE), "k")
cat("=== TEST 1: moisture vs elevation AT THE SURVEY POINTS (no interpolation) ===\n")
if (nrow(J) >= 10) {
  cat(sprintf("  matched points: %d\n", nrow(J)))
  cat(sprintf("  Pearson  r   = %+.3f (p = %.3g)\n",
              cor(J$vwc, J$Elevation), cor.test(J$vwc, J$Elevation)$p.value))
  cat(sprintf("  Spearman rho = %+.3f\n", cor(J$vwc, J$Elevation, method = "spearman")))
  cat("  expected sign: NEGATIVE (low ground is wet)\n")
} else cat(sprintf("  only %d matched points -- too few\n", nrow(J)))

# ---- 2. surface against surface ---------------------------------------------
tr <- geo_transforms()
pe <- tr$inv(E$Latitude, E$Longitude); E$PX <- pe$PX; E$PY <- pe$PY
G <- read.csv("outputs/tables/moisture_surface_grid.csv", stringsAsFactors = FALSE)
efit <- suppressWarnings(Tps(cbind(E$PX, E$PY), E$Elevation))
G$elev <- as.numeric(predict(efit, cbind(G$PX, G$PY)))
z <- function(v) (v - mean(v))/sd(v)
G$z_moist <- z(G$vwc); G$z_elev <- z(G$elev)

cat("\n=== TEST 2: interpolated surfaces, z-scored, on the 2 m stand grid ===\n")
cat(sprintf("  cells %d\n", nrow(G)))
cat(sprintf("  Pearson  r   = %+.3f\n", cor(G$z_moist, G$z_elev)))
cat(sprintf("  Spearman rho = %+.3f\n", cor(G$z_moist, G$z_elev, method = "spearman")))
cat(sprintf("  variance of moisture explained by elevation: %.1f%%\n",
            100*cor(G$z_moist, G$z_elev)^2))

cat("\n=== moisture by elevation quintile (the monotonicity check) ===\n")
G$q <- cut(G$elev, quantile(G$elev, seq(0, 1, 0.2)), include.lowest = TRUE, labels = 1:5)
print(as.data.frame(G %>% group_by(q) %>%
  summarise(elev_mean = round(mean(elev),1), vwc_mean = round(mean(vwc),1),
            vwc_rel = round(mean(vwc_rel),2), n = dplyr::n(), .groups = "drop")), row.names = FALSE)
cat("  (quintile 1 = lowest ground; VWC should fall from 1 to 5)\n")

p1 <- ggplot(G, aes(PX, PY, fill = vwc)) + geom_raster() +
  scale_fill_viridis_c(option = "mako", name = "VWC %") + coord_equal() +
  labs(title = "(a) Moisture, thin-plate spline", x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
p2 <- ggplot(G, aes(PX, PY, fill = elev)) + geom_raster() +
  scale_fill_viridis_c(option = "cividis", name = "elevation") + coord_equal() +
  labs(title = "(b) Elevation, same spline", x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
p3 <- ggplot(G, aes(z_elev, z_moist)) +
  geom_bin2d(bins = 60) + scale_fill_viridis_c(option = "magma", trans = "log10", name = "cells") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "white", linewidth = 0.6, se = FALSE) +
  labs(title = sprintf("(c) z-scores, r = %+.3f", cor(G$z_moist, G$z_elev)),
       x = "elevation (z)", y = "moisture (z)") + theme_bw(base_size = 9)
ggsave("outputs/figures/generated/fig_moisture_elevation.png", p1 | p2 | p3,
       width = 12, height = 4.4, dpi = 190, bg = "white")

con <- file("outputs/audit/moisture_elevation_check.txt", "w")
cat(sprintf("MOISTURE vs ELEVATION -- independent check\nbuilt %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = con)
cat(sprintf("at survey points (n=%d): Pearson %+.3f, Spearman %+.3f\n", nrow(J),
            ifelse(nrow(J) >= 10, cor(J$vwc, J$Elevation), NA),
            ifelse(nrow(J) >= 10, cor(J$vwc, J$Elevation, method = "spearman"), NA)), file = con)
cat(sprintf("interpolated surfaces (n=%d cells): Pearson %+.3f, Spearman %+.3f\n",
            nrow(G), cor(G$z_moist, G$z_elev),
            cor(G$z_moist, G$z_elev, method = "spearman")), file = con)
close(con)
cat("\nwritten: outputs/audit/moisture_elevation_check.txt and fig_moisture_elevation.png\n")
