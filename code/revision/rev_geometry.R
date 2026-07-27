#!/usr/bin/env Rscript
# ==============================================================================
# rev_geometry.R  --  THE stand geometry. Source this; never redefine plot area.
# ------------------------------------------------------------------------------
# Every script that divides by a ground area, clips a surface, or draws the plot
# outline must get that geometry from here. Before this module existed, plot area
# was hardcoded as 200^2 in the budget while the soil surface was averaged over
# its own larger, ragged extent, so the two terms of the budget described
# different ground.
#
# THE STAND IS NOT THE FULL SQUARE, AND THE GAP IS NOT A RECTANGLE.
# The census covers the 200 x 200 m ForestGEO plot minus an unsurveyed region in
# the bottom-right corner. On a 10 m grid that region is a STAIRCASE of 26 cells
# (2,600 m2), contiguous and reaching the plot edge -- not the single block an
# earlier version of this file assumed, which both over-cut (it removed occupied
# ground at PY 20-30) and under-cut (it missed PX 140-160 at PY 0-20).
#
# The region is derived, not drawn: take located stems, bin to 10 m, and
# flood-fill empty cells from the bottom-right corner. That distinguishes the
# unsurveyed block from the 13 ISOLATED empty cells elsewhere in the plot, which
# are ordinary forest gaps and stay in the stand. rev_inventory_build.R re-derives
# it from the data each run and warns if it no longer matches the constant below.
#
# Esri World Imagery shows continuous closed canopy over the region, so it is
# unsurveyed forest, not absent forest: the stem census is deficient there but
# soil flux genuinely applies. We exclude it from BOTH terms rather than pretend
# either is complete over it, and the stand is defined as the censused ground.
#
# Consequence, and it differs by term:
#   * the TREE term is a sum over discrete stems divided by ground area, so it
#     scales with the denominator;
#   * the SOIL term is the spatial MEAN of a predicted surface, so clipping does
#     not renormalise it -- it changes only because the mean over the clipped
#     domain differs from the mean over the old one.
# The tree:soil ratio is therefore NOT invariant under this change.
#
# Coordinates: PX/PY are plot-local metres. The plot is rotated ~10 degrees from
# north, so a lon/lat bounding box is NOT the plot -- always clip in PX/PY.
# ==============================================================================

PLOT_SIDE_M <- 200
GAP_CELL_M  <- 10          # resolution at which the unsurveyed region is defined

# Per 10 m PY band, the PX at which the unsurveyed region begins; NA = fully
# censused. Index i covers PY [(i-1)*10, i*10). Derived by flood-fill; see header.
# The step back out at PY 20-30 (170 rather than 160) is real: the cell
# PX 160-170 there holds two stems, so it was surveyed.
GAP_PX_MIN_BY_PY <- c(140, 150, 170, 160, 160, 160,
                      rep(NA_real_, PLOT_SIDE_M/GAP_CELL_M - 6))

#' PX at which the gap starts for the band containing PY; Inf where none.
.gap_px <- function(PY) {
  i <- pmin(pmax(floor(PY/GAP_CELL_M) + 1, 1), length(GAP_PX_MIN_BY_PY))
  v <- GAP_PX_MIN_BY_PY[i]
  ifelse(is.na(v), Inf, v)
}

PLOT_AREA_M2 <- PLOT_SIDE_M^2                                        # 40,000 nominal
GAP_AREA_M2  <- sum(ifelse(is.na(GAP_PX_MIN_BY_PY), 0,
                           PLOT_SIDE_M - GAP_PX_MIN_BY_PY) * GAP_CELL_M)  # 2,600
STAND_AREA_M2 <- PLOT_AREA_M2 - GAP_AREA_M2                          # 37,400 = 3.74 ha

# --- membership ---------------------------------------------------------------
#' Inside the unsurveyed staircase. Uses >= on PX so a stem recorded exactly on a
#' cell edge belongs to the stand rather than the gap; no located stem falls
#' inside the region by construction.
in_gap <- function(PX, PY) {
  !is.na(PX) & !is.na(PY) & PY >= 0 & PY < PLOT_SIDE_M & PX >= .gap_px(PY)
}
in_square <- function(PX, PY) {
  !is.na(PX) & !is.na(PY) & PX >= 0 & PX <= PLOT_SIDE_M & PY >= 0 & PY <= PLOT_SIDE_M
}
#' TRUE for points on censused ground. Points with no coordinate return FALSE;
#' such stems are real and must still be carried in the budget (see
#' rev_inventory_build.R), so test `located` separately from `in_stand`.
in_stand <- function(PX, PY) in_square(PX, PY) & !in_gap(PX, PY)

# Retained for readability at call sites that predate the staircase.
in_notch <- in_gap
NOTCH_AREA_M2 <- GAP_AREA_M2

# --- outline, plot-local ------------------------------------------------------
#' Closed rectilinear ring of the stand boundary in PX/PY, built from the same
#' per-band thresholds as the membership test so the two cannot drift apart.
stand_ring_local <- function() {
  n <- length(GAP_PX_MIN_BY_PY)
  edge <- ifelse(is.na(GAP_PX_MIN_BY_PY), PLOT_SIDE_M, GAP_PX_MIN_BY_PY)
  PX <- 0; PY <- 0                              # start at the origin
  for (i in seq_len(n)) {                       # walk the right edge upward
    PX <- c(PX, edge[i], edge[i])
    PY <- c(PY, (i-1)*GAP_CELL_M, i*GAP_CELL_M)
  }
  PX <- c(PX, 0, 0); PY <- c(PY, PLOT_SIDE_M, 0)
  # drop consecutive duplicates left by equal-threshold neighbours
  keep <- c(TRUE, (PX[-1] != head(PX, -1)) | (PY[-1] != head(PY, -1)))
  data.frame(PX = PX[keep], PY = PY[keep])
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
stand_ring_lonlat <- function(tr = geo_transforms(), n_per_edge = 12) {
  r <- stand_ring_local()
  seg <- do.call(rbind, lapply(seq_len(nrow(r) - 1), function(i)
    data.frame(PX = seq(r$PX[i], r$PX[i+1], length.out = n_per_edge),
               PY = seq(r$PY[i], r$PY[i+1], length.out = n_per_edge))))
  tr$fwd(seg$PX, seg$PY)
}

# --- self-check ---------------------------------------------------------------
# The ring and the membership test are independent constructions from the same
# constant; a shoelace area on the ring must reproduce STAND_AREA_M2.
local({
  r <- stand_ring_local()
  a <- abs(sum(r$PX[-nrow(r)]*r$PY[-1] - r$PX[-1]*r$PY[-nrow(r)]))/2
  if (abs(a - STAND_AREA_M2) > 1e-6)
    stop(sprintf("rev_geometry.R: ring area %.1f != STAND_AREA_M2 %.1f", a, STAND_AREA_M2))
})

if (identical(environment(), globalenv()) && !interactive()) {
  cat(sprintf("stand geometry: %d x %d m minus a %d-cell staircase gap\n",
              PLOT_SIDE_M, PLOT_SIDE_M, GAP_AREA_M2/GAP_CELL_M^2))
  cat(sprintf("  nominal %d m2 | gap %d m2 | STAND %d m2 = %.2f ha\n",
              PLOT_AREA_M2, GAP_AREA_M2, STAND_AREA_M2, STAND_AREA_M2/1e4))
}
