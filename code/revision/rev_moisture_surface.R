#!/usr/bin/env Rscript
# ==============================================================================
# rev_moisture_surface.R
# ------------------------------------------------------------------------------
# THE soil-moisture spatial pattern, on a regular plot-local grid.
#
# Built with a thin-plate spline over the 71 December survey points, chosen by
# the comparison in rev_moisture_interpolation.R. Replaces akima::interp() at its
# default linear = TRUE, which fits flat triangles between points -- the source of
# the planar facets and creases in the Figure 9 soil panel -- and which is only
# defined inside the convex hull, leaving 29% of the plot unfillable.
#
# Why TPS rather than the method with the lowest leave-one-out error. IDW power 2
# scores best on LOO (RMSE 15.6 against 16.8) and is unusable: it places a
# bullseye at every survey point with a dark halo around it. That artifact is
# exactly what LOO rewards, because reproducing each point at its own location is
# what the criterion measures, and it is not what this surface is for. The
# surface only has to carry large-scale wet/dry structure -- the plot-wide LEVEL
# is imposed separately by rev_moisture_climatology.R -- so smoothness and full
# coverage matter more than pointwise error. TPS is 5x smoother than akima
# (roughness 0.183 vs 0.949), fills 100% of the plot, and is near-unbiased.
#
# What this surface is NOT. It is a relative pattern, normalised to mean 1 over
# the censused stand; it carries no absolute moisture level. And its absolute
# skill is poor -- LOO RMSE ~16 VWC% against a survey mean of 22.6 -- so 71 points
# at 13 m median spacing support large-scale structure and nothing finer.
#
# Caveat inherited from the input, and worth stating in the methods: river points
# are inserted into the survey with an ASSUMED VWC of 100%, which is not a
# measurement, and it anchors the wet end of the gradient.
#
# Output: outputs/tables/moisture_surface_grid.csv
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(fields)})
set.seed(42)
source("code/revision/rev_geometry.R")

CELL_M <- 2      # metres; 100 x 100 cells over the plot, matching the survey's limits

SRC <- "data/raw/inventory/spatial_data/soil_moisture_20201216.csv"
if (!file.exists(SRC)) stop("missing ", SRC)
D <- read.csv(SRC, stringsAsFactors = FALSE)
xc <- intersect(c("Longitude","lon","x","X"), names(D))[1]
yc <- intersect(c("Latitude","lat","y","Y"), names(D))[1]
D <- D[is.finite(D[[xc]]) & is.finite(D[[yc]]) & is.finite(D$VWC), ]
D <- data.frame(lon = D[[xc]], lat = D[[yc]], vwc = D$VWC)
D <- D[!duplicated(round(cbind(D$lon, D$lat), 8)), ]

tr <- geo_transforms()
p <- tr$inv(D$lat, D$lon); D$PX <- p$PX; D$PY <- p$PY
cat(sprintf("survey points %d | VWC %.1f-%.1f (mean %.1f)\n", nrow(D), min(D$vwc), max(D$vwc), mean(D$vwc)))
cat(sprintf("in plot-local metres: PX %.0f..%.0f  PY %.0f..%.0f\n",
            min(D$PX), max(D$PX), min(D$PY), max(D$PY)))
inside <- in_square(D$PX, D$PY)
cat(sprintf("survey points inside the 200 m square: %d of %d\n", sum(inside), nrow(D)))

# --- fit and predict on the stand grid ---------------------------------------
# FIT ON EVERY SURVEY POINT, CLIP AFTERWARDS. The spline is fitted to all 71
# points, including the 13 that lie outside the 200 m square: they carry real
# information about the gradient near the boundary, and dropping them first would
# leave the edge of the surface extrapolating from nothing. Only the PREDICTION
# is restricted to censused ground.
fit <- suppressWarnings(Tps(cbind(D$PX, D$PY), D$vwc))
cat(sprintf("spline fitted on all %d survey points (%d of them outside the square)\n",
            nrow(D), sum(!in_square(D$PX, D$PY))))
half <- CELL_M/2
G <- expand.grid(PX = seq(half, PLOT_SIDE_M - half, by = CELL_M),
                 PY = seq(half, PLOT_SIDE_M - half, by = CELL_M))
G$in_stand <- in_stand(G$PX, G$PY)
G <- G[G$in_stand, ]
G$vwc <- as.numeric(predict(fit, cbind(G$PX, G$PY)))
# a spline can undershoot below zero in sparsely sampled corners; moisture cannot
n_neg <- sum(G$vwc < 0)
G$vwc <- pmax(G$vwc, min(D$vwc))
cat(sprintf("grid %d cells at %d m (%.2f ha); %d cells floored at the survey minimum\n",
            nrow(G), CELL_M, nrow(G)*CELL_M^2/1e4, n_neg))

ll <- tr$fwd(G$PX, G$PY); G$x <- ll$lon; G$y <- ll$lat
G$vwc_rel <- G$vwc / mean(G$vwc)          # relative pattern, mean 1 over the stand

cat(sprintf("predicted VWC %.1f - %.1f (mean %.1f) | relative %.2f - %.2f\n",
            min(G$vwc), max(G$vwc), mean(G$vwc), min(G$vwc_rel), max(G$vwc_rel)))

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(G[, c("PX","PY","x","y","vwc","vwc_rel")],
          "outputs/tables/moisture_surface_grid.csv", row.names = FALSE)
cat("written: outputs/tables/moisture_surface_grid.csv\n")
