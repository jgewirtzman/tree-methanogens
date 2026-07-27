#!/usr/bin/env Rscript
# ==============================================================================
# rev_geometry.R  --  THE stand geometry. Source this; never redefine plot area.
# ------------------------------------------------------------------------------
# Every script that divides by a ground area, clips a surface, or draws the plot
# outline must get that geometry from here. Before this module existed, the plot
# area was hardcoded as 200^2 in the budget while the soil surface was averaged
# over its own (larger, ragged) extent, so the two terms of the budget described
# different ground.
#
# THE STAND IS NOT THE FULL SQUARE. The census covers the 200 x 200 m ForestGEO
# plot minus a rectangular block in the bottom-right corner that was never
# surveyed: PX 160-200, PY 0-60. On a 20 m quadrat grid that block holds five
# quadrats with zero stems and one with two stems, and it is the only such gap in
# the plot -- every other quadrat carries 11-994 stems. Esri World Imagery shows
# continuous closed canopy over the block, so it is unsurveyed forest, not absent
# forest: the stem census is deficient there, but soil flux genuinely applies.
# We therefore exclude the block from BOTH terms rather than pretending either is
# complete over it, and the stand is defined as the censused ground.
#
# Consequence, and it differs by term:
#   * the TREE term is a sum over discrete stems divided by ground area, so it
#     scales with the denominator (x 40000/37600 = 1.0638);
#   * the SOIL term is the spatial MEAN of a predicted surface, so clipping does
#     not renormalise it -- it changes only because the mean over the clipped
#     domain differs from the mean over the old one.
# The tree:soil ratio is therefore NOT invariant under this change.
#
# Coordinates: PX/PY are plot-local metres. The plot is rotated ~10 degrees from
# north, so a lon/lat bounding box is NOT the plot -- always clip in PX/PY.
# ==============================================================================

PLOT_SIDE_M <- 200
NOTCH_PX_MIN <- 160    # unsurveyed block: PX 160-200 ...
NOTCH_PY_MAX <- 60     #                   ... PY 0-60

PLOT_AREA_M2  <- PLOT_SIDE_M^2                                    # 40,000 nominal
NOTCH_AREA_M2 <- (PLOT_SIDE_M - NOTCH_PX_MIN) * NOTCH_PY_MAX      #  2,400 excluded
STAND_AREA_M2 <- PLOT_AREA_M2 - NOTCH_AREA_M2                     # 37,600 = 3.76 ha

# --- membership ---------------------------------------------------------------
#' Strict inequalities, so stems recorded exactly ON the notch edge belong to the
#' stand, not the gap. Three stems sit on that line (two at PX = 160.0, one at
#' PY = 60.0); the notch INTERIOR contains zero stems, which is what makes the
#' block identifiable as unsurveyed rather than merely sparse.
in_notch <- function(PX, PY) {
  !is.na(PX) & !is.na(PY) & PX > NOTCH_PX_MIN & PY < NOTCH_PY_MAX
}
in_square <- function(PX, PY) {
  !is.na(PX) & !is.na(PY) & PX >= 0 & PX <= PLOT_SIDE_M & PY >= 0 & PY <= PLOT_SIDE_M
}
#' TRUE for points on censused ground. Points with no coordinate return FALSE;
#' such stems are real and must still be carried in the budget (see
#' rev_inventory_build.R), so test `located` separately from `in_stand`.
in_stand <- function(PX, PY) in_square(PX, PY) & !in_notch(PX, PY)

# --- outline, plot-local ------------------------------------------------------
#' Closed ring of the stand boundary in PX/PY. Six corners: the square with the
#' bottom-right notch cut out.
stand_ring_local <- function() {
  data.frame(
    PX = c(0, NOTCH_PX_MIN, NOTCH_PX_MIN, PLOT_SIDE_M, PLOT_SIDE_M, 0, 0),
    PY = c(0, 0, NOTCH_PY_MAX, NOTCH_PY_MAX, PLOT_SIDE_M, PLOT_SIDE_M, 0))
}
square_ring_local <- function() {
  data.frame(PX = c(0, PLOT_SIDE_M, PLOT_SIDE_M, 0, 0),
             PY = c(0, 0, PLOT_SIDE_M, PLOT_SIDE_M, 0))
}

# --- lon/lat projection -------------------------------------------------------
#' Geodetic parameters live in the workflow RData alongside the transform itself.
#' Returns a list with `fwd(PX,PY)` and `inv(lat,lon)`, both vectorised.
geo_transforms <- function(rdata = "data/processed/integrated/rf_workflow_input_data_with_2023.RData") {
  e <- new.env(); load(rdata, envir = e); w <- e$rf_workflow_data
  lat0 <- w$ref_lat * pi/180; lon0 <- w$ref_lon * pi/180
  rot  <- w$rotation_angle;   er   <- w$earth_radius
  list(
    fwd = function(PX, PY) {
      r <- w$geodetic_transform(PX, PY, lat0, lon0, rot, er)
      data.frame(lon = r$longitude, lat = r$latitude)
    },
    inv = function(lat, lon) {
      r <- w$geodetic_inverse(lat, lon, lat0, lon0, rot, er)
      data.frame(PX = r$x, PY = r$y)
    })
}

#' Stand outline in lon/lat. Edges are densified because the plot is rotated:
#' a straight line in PX/PY is not a straight line in lon/lat.
stand_ring_lonlat <- function(tr = geo_transforms(), n_per_edge = 40) {
  r <- stand_ring_local()
  seg <- do.call(rbind, lapply(seq_len(nrow(r) - 1), function(i)
    data.frame(PX = seq(r$PX[i], r$PX[i+1], length.out = n_per_edge),
               PY = seq(r$PY[i], r$PY[i+1], length.out = n_per_edge))))
  tr$fwd(seg$PX, seg$PY)
}

if (identical(environment(), globalenv()) && !interactive()) {
  cat(sprintf("stand geometry: %d x %d m minus notch (PX %d-%d, PY 0-%d)\n",
              PLOT_SIDE_M, PLOT_SIDE_M, NOTCH_PX_MIN, PLOT_SIDE_M, NOTCH_PY_MAX))
  cat(sprintf("  nominal %d m2 | notch %d m2 | STAND %d m2 = %.2f ha\n",
              PLOT_AREA_M2, NOTCH_AREA_M2, STAND_AREA_M2, STAND_AREA_M2/1e4))
}
