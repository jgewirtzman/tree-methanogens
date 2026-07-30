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
# THE STAND IS NOT THE FULL SQUARE. The census is organised in 20 m quadrats, and
# SEVEN OF THE 100 WERE NEVER CENSUSED -- they appear nowhere in the raw census
# tables. The stand is the 93 that were:
#
#     quadrat (row, col, 0-indexed)   PX        PY
#     q7  q8  q9                      140-200    0- 20
#     q108 q109                       160-200   20- 40
#     q208 q209                       160-200   40- 60
#
# giving a monotonic rectilinear staircase, 2,800 m2 excluded, stand 37,200 m2.
#
# THIS IS THE SURVEY DESIGN, NOT AN INFERENCE. An earlier version of this file
# derived the gap by flood-filling empty 10 m cells from the bottom-right corner.
# That was the wrong basis twice over: it produced a non-monotonic boundary with a
# re-entrant step (an artifact of one 10 m cell that happened to hold two stems),
# and in principle it could not distinguish a quadrat that was never walked from
# one that was walked and found empty. Masking on the quadrat identifiers answers
# that question directly.
#
# Three located stems fall inside the seven uncensused quadrats (one in q7, two in
# q108) and are excluded with them. They are almost certainly stems recorded from
# an adjoining quadrat, or carrying a few metres of coordinate error.
#
# Esri World Imagery shows continuous closed canopy over the block, so it is
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
GAP_CELL_M  <- 20          # ForestGEO quadrat size; the census's own survey unit

# The seven quadrats absent from the raw census, as codes row*100 + col with row
# and col 0-indexed from the plot origin. Verified against the raw tables on every
# run by rev_inventory_build.R.
UNCENSUSED_QUADRATS <- c(7, 8, 9, 108, 109, 208, 209)

# Per 20 m PY band, the PX at which uncensused ground begins; NA = fully censused.
# Derived from UNCENSUSED_QUADRATS so the two cannot disagree.
GAP_PX_MIN_BY_PY <- local({
  v <- rep(NA_real_, PLOT_SIDE_M/GAP_CELL_M)
  for (q in UNCENSUSED_QUADRATS) {
    r <- q %/% 100; cc <- q %% 100
    v[r + 1] <- min(v[r + 1], cc*GAP_CELL_M, na.rm = TRUE)
  }
  v
})

#' PX at which uncensused ground starts for the 20 m band containing PY; Inf where
#' the whole band was censused.
.gap_px <- function(PY) {
  i <- pmin(pmax(floor(PY/GAP_CELL_M) + 1, 1), length(GAP_PX_MIN_BY_PY))
  v <- GAP_PX_MIN_BY_PY[i]
  ifelse(is.na(v), Inf, v)
}

PLOT_AREA_M2 <- PLOT_SIDE_M^2                                        # 40,000 nominal
GAP_AREA_M2  <- length(UNCENSUSED_QUADRATS) * GAP_CELL_M^2           # 2,800
STAND_AREA_M2 <- PLOT_AREA_M2 - GAP_AREA_M2                          # 37,200 = 3.72 ha

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
  # Callable from the repo root or from code/04_scaling/, whose scripts use ../../
  if (!file.exists(rdata) && file.exists(file.path("../..", rdata)))
    rdata <- file.path("../..", rdata)
  if (!file.exists(rdata)) stop("geo_transforms(): cannot find ", rdata)
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

# ==============================================================================
# THE STAND'S STEMS, AND THEIR HEIGHTS. Source of both, so no script re-derives.
# ------------------------------------------------------------------------------
# WHY THIS IS HERE. Seven scripts loaded the stem list themselves, and five of them
# reached for `rf_workflow_data$PLACEHOLDER_INVENTORY` with `A_PLOT <- 200*200`.
# That object is fg19-only: 7,010 stems, with the ENTIRE bytag survey (961 stems)
# missing, over the nominal square rather than the censused stand. So those scripts
# reported a stand 12% short of stems on ground 7.5% too large -- 3,554 m2 of band
# area and a stem-area index of 0.0888 against the canonical 4,325.9 and 0.1163,
# i.e. 24% low. One of them, rev_wai_bottomup_and_rf_interactions.R, is the sole
# producer of the bottom-up woody area index, which came out 1.69/2.11/2.57 when the
# canonical inventory gives 2.23/2.82/3.40 -- 32% low in a quantity that multiplies
# roughly 70% of the headline scaling estimate.
#
# The height allometry was separately copy-pasted verbatim into five scripts. Its
# scale coefficient is anchored on the 95th percentile of the inventory's own
# diameters, so it silently differed by which stem list the script happened to load.
#
# Both now live here. `PLACEHOLDER_INVENTORY` and `PLACEHOLDER_PLOT_AREA` must not
# be read by anything again; prefer canonical_inventory() and STAND_AREA_M2.
# ==============================================================================

# Hulshof et al. 2015 Ecol Evol 5:1193-1204 exponents; site-specific scale.
ALLOM_B_ANGIO <- 0.53
ALLOM_B_GYMNO <- 0.60
CANOPY_H_M    <- 25        # midpoint of the observed 20-30 m closed-canopy range
GYMNOSPERMS   <- c("Pinus strobus", "Tsuga canadensis")

#' Per-stem height from diameter, H = 1.37 + a * DBH^b.
#' `a` is fitted so a stem at the 95th percentile of diameter (over stems > 10 cm)
#' reaches CANOPY_H_M. No shrub cap, by decision -- see scaling_parameters.md.
stem_height_m <- function(dbh_m, species, canopy_h = CANOPY_H_M) {
  gy <- species %in% GYMNOSPERMS
  d_anchor <- stats::quantile(dbh_m[dbh_m > 0.10], 0.95, na.rm = TRUE)
  a <- ifelse(gy, (canopy_h - 1.37)/d_anchor^ALLOM_B_GYMNO,
                  (canopy_h - 1.37)/d_anchor^ALLOM_B_ANGIO)
  1.37 + a * dbh_m^ifelse(gy, ALLOM_B_GYMNO, ALLOM_B_ANGIO)
}

#' THE stem list. Reads what rev_inventory_build.R wrote; adds dbh, H and band area.
#' in_stand_only = TRUE gives the 8,006 stems the budget is defined over.
canonical_inventory <- function(in_stand_only = TRUE,
                                file = "outputs/tables/inventory_stems.csv") {
  # Callable from the repo root or from code/04_scaling/, whose scripts use ../../
  if (!file.exists(file) && file.exists(file.path("../..", file)))
    file <- file.path("../..", file)
  if (!file.exists(file))
    stop("canonical_inventory(): missing ", file,
         " -- run: Rscript code/revision/rev_inventory_build.R")
  INV <- utils::read.csv(file, stringsAsFactors = FALSE)
  INV <- INV[is.finite(INV$dbh_m) & INV$dbh_m > 0, ]
  if (in_stand_only) INV <- INV[INV$in_stand, ]
  INV$dbh <- INV$dbh_m
  INV$H   <- stem_height_m(INV$dbh_m, INV$species)
  # Kalmia is a shrub measured over 0.75 m, not 2 m -- matches the band used by
  # rev_predict_tree_flux_current.R. %in% is NA-safe; 41 stems have species = NA.
  INV$band_m     <- ifelse(INV$species %in% "Kalmia latifolia", 0.75, 2.00)
  INV$A_band_m2  <- pi * INV$dbh_m * INV$band_m
  INV
}

if (identical(environment(), globalenv()) && !interactive()) {
  cat(sprintf("stand geometry: %d x %d m minus %d uncensused %d m quadrats\n",
              PLOT_SIDE_M, PLOT_SIDE_M, length(UNCENSUSED_QUADRATS), GAP_CELL_M))
  cat(sprintf("  nominal %d m2 | gap %d m2 | STAND %d m2 = %.2f ha\n",
              PLOT_AREA_M2, GAP_AREA_M2, STAND_AREA_M2, STAND_AREA_M2/1e4))
}
