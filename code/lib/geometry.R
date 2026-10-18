#!/usr/bin/env Rscript
# ==============================================================================
# geometry.R -- the stand geometry, the canonical stem list and stem heights
# ------------------------------------------------------------------------------
# Source this; never redefine plot area. Every script that divides by a ground area,
# clips a surface, or draws the plot outline takes the geometry from here.
#
# The stand is not the full 200 m square. The census is organised in 20 m quadrats,
# and seven of the 100 were never censused (they appear nowhere in the raw census
# tables). The stand is the 93 that were:
#
#     quadrat (row, col, 0-indexed)   PX        PY
#     q7  q8  q9                      140-200    0- 20
#     q108 q109                       160-200   20- 40
#     q208 q209                       160-200   40- 60
#
# a rectilinear staircase: 2,800 m2 excluded, stand 37,200 m2. Three located stems
# fall inside the uncensused quadrats (one in q7, two in q108) and are excluded with
# them. The block is closed-canopy forest (Esri World Imagery) whose stem census is
# missing, so it is excluded from both the tree and the soil term.
# The tree term is a sum over stems divided by ground area, so it scales with the
# denominator; the soil term is the mean of a predicted surface over the stand.
#
# Coordinates: PX/PY are plot-local metres. The plot is rotated ~10 degrees from
# north, so a lon/lat bounding box is not the plot: always clip in PX/PY.
# ==============================================================================

PLOT_SIDE_M <- 200
GAP_CELL_M  <- 20          # plot census quadrat size; the census's own survey unit

# The seven quadrats absent from the raw census, as codes row*100 + col with row
# and col 0-indexed from the plot origin. Verified against the raw tables on every
# run by 08_inventory_build.R.
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
#' 08_inventory_build.R), so test `located` separately from `in_stand`.
in_stand <- function(PX, PY) in_square(PX, PY) & !in_gap(PX, PY)

# the name the inventory build uses
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

# --- lon/lat projection -------------------------------------------------------
# Plot-local PX/PY (m) <-> WGS84: great-circle offset from the reference point after
# rotating the plot axes. The one definition of the reference point and the
# transform; data/tables/site.csv (02_site.R) and the model's .RData carry copies.
GEO_REF_LAT      <- 41.989211   # plot reference point, decimal degrees
GEO_REF_LON      <- -72.13092
GEO_ROTATION_RAD <- 0.174533    # plot axes vs true north (10 degrees)
GEO_EARTH_RADIUS <- 6371000     # m

# Forward transformation: plot coords (m) to GPS (decimal degrees)
geodetic_transform <- function(px, py, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  distance <- sqrt(px^2 + py^2)
  original_bearing <- atan2(px, py)
  rotated_bearing <- original_bearing + rotation_angle

  px_rotated <- distance * sin(rotated_bearing)
  py_rotated <- distance * cos(rotated_bearing)

  bearing <- atan2(px_rotated, py_rotated)
  distance_m <- sqrt(px_rotated^2 + py_rotated^2)

  lat2_rad <- asin(sin(lat1_rad) * cos(distance_m/earth_radius) +
                     cos(lat1_rad) * sin(distance_m/earth_radius) * cos(bearing))

  lon2_rad <- lon1_rad + atan2(sin(bearing) * sin(distance_m/earth_radius) * cos(lat1_rad),
                               cos(distance_m/earth_radius) - sin(lat1_rad) * sin(lat2_rad))

  return(list(longitude = lon2_rad * 180/pi, latitude = lat2_rad * 180/pi))
}

# Inverse transformation: GPS to plot coords
geodetic_inverse <- function(lat, lon, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  lat_rad <- lat * pi/180
  lon_rad <- lon * pi/180
  delta_lon <- lon_rad - lon1_rad

  y_comp <- sin(delta_lon) * cos(lat_rad)
  x_comp <- cos(lat1_rad) * sin(lat_rad) - sin(lat1_rad) * cos(lat_rad) * cos(delta_lon)
  bearing <- atan2(y_comp, x_comp)

  a <- sin((lat_rad - lat1_rad)/2)^2 + cos(lat1_rad) * cos(lat_rad) * sin(delta_lon/2)^2
  c <- 2 * asin(sqrt(a))
  distance_m <- earth_radius * c

  bearing_adjusted <- bearing - rotation_angle
  x <- distance_m * sin(bearing_adjusted)
  y <- distance_m * cos(bearing_adjusted)

  return(list(x = x, y = y))
}

#' Returns a list with `fwd(PX,PY)` and `inv(lat,lon)`, both vectorised.
geo_transforms <- function() {
  lat0 <- GEO_REF_LAT * pi/180; lon0 <- GEO_REF_LON * pi/180
  rot  <- GEO_ROTATION_RAD;     er   <- GEO_EARTH_RADIUS
  list(
    fwd = function(PX, PY) {
      r <- geodetic_transform(PX, PY, lat0, lon0, rot, er)
      data.frame(lon = r$longitude, lat = r$latitude)
    },
    inv = function(lat, lon) {
      r <- geodetic_inverse(lat, lon, lat0, lon0, rot, er)
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
    stop(sprintf("geometry.R: ring area %.1f != STAND_AREA_M2 %.1f", a, STAND_AREA_M2))
})

# ==============================================================================
# THE STAND'S STEMS: the canonical stem list, stem heights, the measured band and
# the bole:branch split of woody area. Scripts take these from here; none
# re-derives them.
# ==============================================================================

# Measured band: stems 0-2 m; Kalmia latifolia (a shrub) over 0.75 m.
STEM_BAND_M   <- 2.00
KALMIA_BAND_M <- 0.75
band_length_m <- function(species) ifelse(species %in% "Kalmia latifolia", KALMIA_BAND_M, STEM_BAND_M)

# Woody area above the band: branch:stem area 3.35 (Whittaker & Woodwell; midpoint of
# 2.7-4.0), so the bole carries 1/4.35 of it and branches 3.35/4.35.
BRANCH_STEM_RATIO <- 3.35
BOLE_FRACTION     <- 1 / (1 + BRANCH_STEM_RATIO)
BRANCH_FRACTION   <- BRANCH_STEM_RATIO / (1 + BRANCH_STEM_RATIO)

# Hulshof et al. 2015 Ecol Evol 5:1193-1204 exponents; site-specific scale.
ALLOM_B_ANGIO <- 0.53
ALLOM_B_GYMNO <- 0.60
CANOPY_H_M    <- 25        # midpoint of the observed 20-30 m closed-canopy range
GYMNOSPERMS   <- c("Pinus strobus", "Tsuga canadensis")

#' Per-stem height from diameter, H = 1.37 + a * DBH^b.
#' `a` is fitted so a stem at the 95th percentile of diameter (over stems > 10 cm)
#' reaches CANOPY_H_M. No shrub cap: the allometry is applied uniformly.
stem_height_m <- function(dbh_m, species, canopy_h = CANOPY_H_M) {
  gy <- species %in% GYMNOSPERMS
  d_anchor <- stats::quantile(dbh_m[dbh_m > 0.10], 0.95, na.rm = TRUE)
  a <- ifelse(gy, (canopy_h - 1.37)/d_anchor^ALLOM_B_GYMNO,
                  (canopy_h - 1.37)/d_anchor^ALLOM_B_ANGIO)
  1.37 + a * dbh_m^ifelse(gy, ALLOM_B_GYMNO, ALLOM_B_ANGIO)
}

#' The stem list. Reads what 08_inventory_build.R wrote; adds dbh, H and band area.
#' in_stand_only = TRUE gives the stems the budget is defined over.
canonical_inventory <- function(in_stand_only = TRUE,
                                file = "outputs/tables/inventory_stems.csv") {
  # Callable from the repo root or from code/05_model/, whose scripts use ../../
  if (!file.exists(file) && file.exists(file.path("../..", file)))
    file <- file.path("../..", file)
  if (!file.exists(file))
    stop("canonical_inventory(): missing ", file,
         " -- run: Rscript code/01_import/08_inventory_build.R")
  INV <- utils::read.csv(file, stringsAsFactors = FALSE)
  INV <- INV[is.finite(INV$dbh_m) & INV$dbh_m > 0, ]
  if (in_stand_only) INV <- INV[INV$in_stand, ]
  INV$dbh <- INV$dbh_m
  INV$H   <- stem_height_m(INV$dbh_m, INV$species)
  INV$band_m     <- band_length_m(INV$species)   # NA-safe: some stems have species NA
  INV$A_band_m2  <- pi * INV$dbh_m * INV$band_m
  INV
}

if (identical(environment(), globalenv()) && !interactive()) {
  cat(sprintf("stand geometry: %d x %d m minus %d uncensused %d m quadrats\n",
              PLOT_SIDE_M, PLOT_SIDE_M, length(UNCENSUSED_QUADRATS), GAP_CELL_M))
  cat(sprintf("  nominal %d m2 | gap %d m2 | STAND %d m2 = %.2f ha\n",
              PLOT_AREA_M2, GAP_AREA_M2, STAND_AREA_M2, STAND_AREA_M2/1e4))
}
