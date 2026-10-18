#!/usr/bin/env Rscript
# ==============================================================================
# 07_scaling_full_grid.R -- the whole-surface scaling grid and the named pair
# ------------------------------------------------------------------------------
# Stage C; run from the repository root. Every unconstrained assumption of the
# whole-surface estimate, crossed:
#
#   WAI (total woody area)   5 levels
#       1.50   Whittaker & Woodwell 1967, Brookhaven oak-pine (low bound)
#       3.07   Gauci et al. 2024, temperate forest (high bound)
#       low / mid / high bottom-up for this inventory with W&W relations (taper
#              1.20-1.35, branch:stem 2.7-4.0), read from wai_bottomup.csv
#              (04_wai_bottomup.R)
#   BOLE SHAPE               2   cone / paraboloid (W&W: true stems lie between them)
#   BRANCH SHAPE             4   uniform_all / uniform_top50 / gaussian_50 / gaussian_75
#   FLUX ABOVE 2 m           6   no stem may increase substantially with height
#       constant          each stem keeps its 2 m value
#       exponential       decay toward zero, fitted to the campaign means at 0.5/1.25/2 m
#       power             gradual decay, never reaching zero
#       linear_floored    linear decline truncated at zero
#       linear_bounded_median  linear decline through zero into uptake, bounded at the
#                         median detected stem uptake
#       exp_band_slope    per-stem decay fitted to the stem's modelled 0-2 m profile,
#                         capped at its own 2 m value
#   An asymptotic form is phi*constant + (1-phi)*exponential, so it is spanned by two
#   forms already present.
#   = length(FORMS) x length(WAIS) x length(BRANCH) x length(BOLE) combinations (240).
#
# Each estimate = MEASURED (woody area at z <= 2 m: the canonical band term of
# 03_budget_canonical.R) + EXTRAPOLATED (area above 2 m x each stem's calibrated 2 m
# flux x the form's shape ratio; nothing above 2 m is measured at this site).
# The named pair (decline = exponential, hold = constant, at the bottom-up "this
# stand" WAI, gaussian_75, cone) is written to scaling_headline.csv.
#
# Inputs: outputs/tables/{inventory_stems, tree_flux_predictions, tree_band_profile,
#   tree_band_profile_edges}.csv, outputs/data/{canonical_budget, wai_bottomup}.csv,
#   data/tables/flux_stem.csv, outputs/models/TRAINING_DATA.RData
# Outputs (outputs/data/): scaling_full_grid.csv, scaling_headline.csv,
#   scaling_slope_diagnostics.csv, and for figS27 scaling_flux_shapes.csv,
#   scaling_shape_kernels.csv, scaling_band_profile.csv, scaling_flux_anchors.csv,
#   scaling_area_profiles.csv; outputs/figures/generated/fig_scaling_full_grid.png
# ==============================================================================
source("code/lib/outputs.R")

suppressMessages({library(ranger);library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
set.seed(42)
source("code/lib/geometry.R")
source("code/lib/species_levels.R")
source("code/lib/flux_transform.R")
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
d <- tree_train_complete; trained <- sort(unique(as.character(d$species_clean)))
CONV <- 86400*365.25*16e-6

# the canonical products: the measured band and the soil term are read, not recomputed
need <- c(inv="outputs/tables/inventory_stems.csv",
          tree="outputs/tables/tree_flux_predictions.csv",
          prof="outputs/tables/tree_band_profile.csv",
          pedge="outputs/tables/tree_band_profile_edges.csv",
          bud="outputs/data/canonical_budget.csv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing inputs:\n  ", paste(miss, collapse="\n  "),
  "\nrun 08_inventory_build.R, 01_predict_tree_flux.R, 02_predict_soil_surface.R, 03_budget_canonical.R first")
B   <- read.csv(need[["bud"]], stringsAsFactors=FALSE)
val <- function(q){v <- B$value[B$quantity==q]; if(!length(v)) stop("canonical_budget.csv lacks ",q); v}
SOIL_ANN  <- val("soil_mg_m2_yr")
MEAS_BAND <- val("tree_measured_mg_m2_yr")   # invariant across the grid, by construction
A_PLOT    <- STAND_AREA_M2                   # censused ground, not the nominal square

TRP <- read.csv(need[["tree"]], stringsAsFactors=FALSE)
INV <- read.csv(need[["inv"]],  stringsAsFactors=FALSE)
INV <- INV[INV$in_stand, ]; TRP <- TRP[TRP$in_stand, ]
stopifnot(nrow(INV)==nrow(TRP))
INV$dbh <- INV$dbh_m
INV$A_band_m2 <- TRP$A_stem_m2                # measured-band area, Kalmia at 0.75 m
INV$sp <- species_to_model_level(INV$species, trained)
INV$H <- stem_height_m(INV$dbh, INV$species)   # the one allometry (geometry.R); no shrub cap
n <- nrow(INV)

# --- per-tree band profile: the calibrated profile of 01_predict_tree_flux.R ----
PE    <- read.csv(need[["pedge"]], stringsAsFactors=FALSE)
PRF   <- read.csv(need[["prof"]],  stringsAsFactors=FALSE)
PRF   <- PRF[PRF$in_stand, ]
stopifnot(nrow(PRF)==n)
edges <- c(PE$edge_lo_cm, tail(PE$edge_hi_cm,1)); K <- nrow(PE)
PROFI <- as.matrix(PRF[, sprintf("f_i%d", seq_len(K)), drop=FALSE])
dimnames(PROFI) <- NULL
stopifnot(all(is.finite(PROFI)))
zh <- (head(edges,-1) + pmin(tail(edges,-1),201) - 1)/2/100   # interval midpoints, m
f2 <- PROFI[,K]                             # 2 m anchor for the above-band forms
# the anchor must agree with the canonical per-stem export, exactly
stopifnot(max(abs(f2 - TRP$flux_2m_nmol_m2_s)) < 1e-9)
cat(sprintf("band profile read: %d stems x %d intervals (edges %s cm), calibrated\n",
            n, K, paste(edges, collapse=", ")))
#' Step lookup: the model's value at height z (metres), constant within an interval.
prof_at <- function(i, zm) PROFI[i, findInterval(pmin(pmax(zm*100,50),200), edges,
                                                 rightmost.closed=TRUE)]
rf_slope <- as.numeric(apply(PROFI,1,function(f){
  if(any(!is.finite(f))||all(f<=0)||sd(f)==0) return(0)
  as.numeric(coef(lm(log(pmax(f,1e-9))~zh))[2])}))

# Diagnostics of the per-stem slopes behind exp_band_slope (scaling_slope_diagnostics.csv)
slope_fit <- function(f, z) {
  if(any(!is.finite(f))||all(f<=0)||sd(f)==0) return(c(NA, NA))
  m <- lm(log(pmax(f,1e-9))~z); c(coef(m)[2], summary(m)$r.squared) }
fit_all   <- t(apply(PROFI, 1, slope_fit, z = zh))
fit_upper <- t(apply(PROFI[, -1, drop=FALSE], 1, slope_fit, z = zh[-1]))   # without the lowest interval
slope_diag <- data.frame(
  quantity = c("median_slope_all_intervals", "median_slope_without_lowest", "lowest_interval_top_cm",
               "median_r2_all_intervals", "fit_span_m", "stems"),
  value = c(median(fit_all[,1], na.rm=TRUE), median(fit_upper[,1], na.rm=TRUE), edges[2],
            median(fit_all[,2], na.rm=TRUE), diff(range(zh)), n))
write.csv(slope_diag, out_path("scaling_slope_diagnostics.csv"), row.names=FALSE)
cat(sprintf("exp_band_slope diagnostics: median slope %.3f per m (all), %.3f without the interval below %d cm; median R2 %.3f over %.2f m\n",
            slope_diag$value[1], slope_diag$value[2], edges[2], slope_diag$value[4], slope_diag$value[5]))

# --- shapes and forms --------------------------------------------------------
hq <- c(.5,1.25,2)
fbar <- sapply(hq,function(z) mean(d$stem_flux_corrected[!is.na(d$measurement_height_cm) &
        abs(d$measurement_height_cm/100-z)<1e-6],na.rm=TRUE))
b_log <- as.numeric(coef(lm(log(pmax(fbar,1e-9))~hq))[2])
b_pow <- as.numeric(coef(lm(log(pmax(fbar,1e-9))~log(hq)))[2])
b_lin <- as.numeric(coef(lm(fbar~hq)))
# cylinder dropped: real stems taper, and W&W show true stems lie BETWEEN the conic
# and parabolic forms, so those two bracket reality.
BOLE <- list(cone=function(u) pmax(1-u,0), paraboloid=function(u) sqrt(pmax(1-u,0)))
# BRANCH PLACEMENT. None is empirically supported here; they bracket plausible crown
# geometries: area low (uniform_all), a hard cutoff (uniform_top50), and Gaussian
# crowns at two depths (0.50H, 0.75H). (uniform_crown, all branch area above 0.60H,
# duplicated uniform_top50 to within 3% and is not carried.)
BRANCH <- list(uniform_all=function(u) rep(1,length(u)),
  uniform_top50=function(u) as.numeric(u>=.50),
  gaussian_50=function(u) dnorm(u,.50,.18),
  gaussian_75=function(u) dnorm(u,.75,.15))
# WAI SCENARIOS. The bottom-up trio is read from 04_wai_bottomup.R (built from W&W's
# taper and branch:stem ratios, so closeness to their published values corroborates
# the arithmetic only). The two literature values are citations and stay literal.
WAI_BU <- local({
  f <- out_path("wai_bottomup.csv")
  if (!file.exists(f))
    stop("missing ", f, " -- run: Rscript code/06_upscale/04_wai_bottomup.R")
  w <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(nrow(w) == 3, all(is.finite(w$wai)), !is.unsorted(w$wai))
  if (!isTRUE(all.equal(w$stand_area_m2[1], STAND_AREA_M2)))
    stop("wai_bottomup.csv was computed over ", w$stand_area_m2[1],
         " m2 but the stand is ", STAND_AREA_M2, " m2")
  setNames(w$wai, sprintf("%.2f %s", w$wai, w$label))
})
WAIS <- c(`1.50 W&W Brookhaven (low bound)` = 1.50,
          WAI_BU,
          `3.07 Gauci 2024 (high bound)` = 3.07)
WAIS <- WAIS[order(WAIS)]
cat(sprintf("WAI axis (%d levels): %s\n", length(WAIS), paste(names(WAIS), collapse=" | ")))
# NEGATIVE BOUND for the linear form: the median detected stem uptake (data/tables/
# flux_stem.csv, detection_class "uptake"; 02_flux/03_precision_and_detection.R). The median, not
# a low percentile: uptake prevalence tracks instrument precision (sigma per campaign,
# 02b_instrument_precision.R) rather than any biological covariate, so a deep tail
# value would come from the two noisier campaigns. The felled oak (all four
# measurements above 2 m positive) contradicts even this bound.
.fs <- read.csv("data/tables/flux_stem.csv", stringsAsFactors = FALSE)
NEG_MEDIAN <- round(median(.fs$CH4_flux_nmol_m2_s[.fs$detection_class %in% "uptake"]), 4)
stopifnot(is.finite(NEG_MEDIAN), NEG_MEDIAN < 0)
cat(sprintf("NEG_MEDIAN (median detected stem uptake, flux_stem.csv) = %.4f\n", NEG_MEDIAN))
rm(.fs)
FORMS <- c("constant","exponential","power","linear_floored",
           "linear_bounded_median","exp_band_slope")

NU <- 90; u <- (seq_len(NU)-.5)/NU
Z  <- outer(INV$H, u)                       # n x NU absolute heights
above <- Z > 2
# per-tree flux inside the band, floored below 0.5 m
Zb <- pmin(pmax(Z,0.5),2)
Fband <- matrix(NA_real_, n, NU)
for (i in seq_len(n)) Fband[i,] <- prof_at(i, Zb[i,])

# Shape ratios relative to each stem's own 2 m flux. linear_bounded_median is the only
# form that goes negative: a decline that continues through zero into net uptake aloft,
# reported in the literature but not observed in this dataset.
f2m_lin <- max(b_lin[1]+b_lin[2]*2, 1e-9)
# The linear forms are absolute-valued below, then converted to a ratio, so the floor
# is applied in flux units rather than as a fraction of f(2 m).
lin_abs <- function(zz) b_lin[1] + b_lin[2]*zz
shape_ratio <- function(form, zz) switch(form,
  constant              = rep(1,length(zz)),
  exponential           = exp(b_log*(zz-2)),
  power                 = exp(b_pow*(log(pmax(zz,.05))-log(2))),
  linear_floored        = pmax(lin_abs(zz), 0)/f2m_lin,
  linear_bounded_median = pmax(lin_abs(zz), NEG_MEDIAN)/f2m_lin)

SHARED <- setdiff(FORMS,"exp_band_slope")
Rmat <- lapply(SHARED, function(fm) matrix(shape_ratio(fm, as.vector(Z)), nrow=n))
names(Rmat) <- SHARED
# exp_band_slope: a log-linear decay whose rate is fitted per stem to the forest's
# band values (interval midpoints), projected above 2 m and capped at 1 so no stem
# exceeds its own 2 m flux. Same family as `exponential`, but the rate comes from
# each stem's modelled within-band profile, a model-derived shape: the fit spans
# ~1.1 m, the lowest interval supplies most of the decline (dropping it changes the
# median slope; scaling_slope_diagnostics.csv), and it is projected to canopy height.
# A lower bound on emission rather than a central case.
Rmat$exp_band_slope <- pmin(exp(rf_slope * (Z-2)), 1)

res <- list()
for (bn in names(BOLE)) for (rn in names(BRANCH)) {
  wb <- BOLE[[bn]](u); wb <- wb/sum(wb)
  wr <- BRANCH[[rn]](u); wr <- wr/sum(wr)
  # Shape of area with height: bole and branch in W&W proportion (geometry.R)
  Ashape <- outer(rep(1,n), wb)*BOLE_FRACTION + outer(rep(1,n), wr)*BRANCH_FRACTION
  Ashape <- Ashape * (pi*INV$dbh*INV$H/2)   # per-tree weighting by conic stem area
  # The 0-2 m band area is MEASURED (pi*DBH*2 per stem) and does NOT depend on which
  # literature WAI is adopted; it is held fixed and spread uniformly over the band
  # (a 2 m cylinder, taper negligible). Only the area ABOVE 2 m absorbs the WAI choice.
  nb <- rowSums(!above)
  Aband <- (!above) * INV$A_band_m2 / pmax(nb,1)
  shp_ab <- Ashape*above; shp_ab <- shp_ab/pmax(rowSums(shp_ab),1e-12)
  for (wn in names(WAIS)) {
    A_above_tot <- WAIS[wn]*A_PLOT - sum(INV$A_band_m2)
    w_tree <- pi*INV$dbh*INV$H/2; w_tree <- w_tree/sum(w_tree)
    Aslab <- Aband + shp_ab * (w_tree * A_above_tot)
    meas <- MEAS_BAND   # the canonical band term
    for (fm in FORMS) {
      extr <- sum((f2*Rmat[[fm]])[above]*Aslab[above])/A_PLOT*CONV
      res[[length(res)+1]] <- data.frame(WAI=wn, bole=bn, branch=rn, flux=fm,
        measured_mg=meas, extrapolated_mg=extr, total_mg=meas+extr,
        pct_extrapolated=100*extr/(meas+extr),
        pct_of_soil=100*(meas+extr)/abs(SOIL_ANN), net_mg=SOIL_ANN+meas+extr)
    }
  }
}
R <- bind_rows(res)
write.csv(R, out_path("scaling_full_grid.csv"), row.names=FALSE)

# ---- THE NAMED PAIR ---------------------------------------------------------------
# The whole-surface estimate is quoted as two cases at one geometry:
#   decline  flux `exponential`: the decay fitted to the observed campaign means at
#            0.5 / 1.25 / 2 m, carried upward. This is the pooled measured decline,
#            and it is mostly yellow birch: without that species the pooled slope
#            is ~0. It is not a validated extrapolator (exponential decay scored last
#            on upward leave-one-height-out CV).
#   hold     flux `constant`: each stem keeps its own 2 m value. Best upward
#            extrapolator in the same CV (13_height_form_crossvalidation.R) and
#            consistent with the climbed oak (0.82x the 2 m flux at 10 m).
# Both shapes come from measurements (exp_band_slope, a model-derived shape, is not
# named). Shared geometry: WAI = the bottom-up estimate for this stand (selected by
# label); branch gaussian_75; bole cone (bole shape is worth <1 mg across the grid).
# The range across the grid remains the result.
HEAD_WAI <- grep("bottom-up, this stand", names(WAIS), value = TRUE)[1]
stopifnot(length(HEAD_WAI) == 1, !is.na(HEAD_WAI))
CASES <- data.frame(case = c("decline", "hold"), flux = c("exponential", "constant"),
                    branch = "gaussian_75", bole = "cone", WAI = HEAD_WAI,
                    stringsAsFactors = FALSE)
HL <- merge(CASES, R, by = c("flux", "branch", "bole", "WAI"), sort = FALSE)
stopifnot(nrow(HL) == nrow(CASES))
HL <- HL[match(CASES$case, HL$case), c("case", "flux", "branch", "bole", "WAI",
         "measured_mg", "extrapolated_mg", "total_mg", "pct_extrapolated",
         "pct_of_soil", "net_mg")]
write.csv(HL, out_path("scaling_headline.csv"), row.names = FALSE)
cat(sprintf("\n=== NAMED CASES === (%s | %s | %s)\n", HEAD_WAI, CASES$branch[1], CASES$bole[1]))
for (i in seq_len(nrow(HL)))
  cat(sprintf("  %-8s %-12s total %.1f mg CH4 m-2 yr-1 (%.1f%% of soil uptake, %.0f%% extrapolated), net %.1f\n",
              HL$case[i], HL$flux[i], HL$total_mg[i], HL$pct_of_soil[i],
              HL$pct_extrapolated[i], HL$net_mg[i]))

# ---- export the shapes for figS27_scaling-profiles.R (one implementation) ------
ZB <- seq(0, 26, by = 0.25)                       # absolute-height bins, metres
zmid <- head(ZB,-1) + diff(ZB)/2

# flux shape ratios, relative to each stem's own 2 m value. All forms except
# exp_band_slope are the same function of height for every stem; that one is a
# per-stem fitted slope, so its stand curve is the area-weighted mean over stems.
fx_prof <- bind_rows(lapply(FORMS, function(fm) {
  if (fm == "exp_band_slope") {
    r <- sapply(zmid, function(z) {
      v <- pmin(exp(rf_slope*(z-2)), 1); sum(v*f2*INV$dbh)/sum(f2*INV$dbh) })
  } else r <- shape_ratio(fm, zmid)
  data.frame(flux = fm, z = zmid, ratio = as.numeric(r))
}))
write.csv(fx_prof, out_path("scaling_flux_shapes.csv"), row.names=FALSE)

# Raw shape kernels on a single tree (the stand profiles below also carry the
# inventory's mix of stem heights)
kern <- bind_rows(
  bind_rows(lapply(names(BOLE), function(b) {
    w <- BOLE[[b]](u); data.frame(part = "bole", shape = b, u = u, w = w/sum(w)) })),
  bind_rows(lapply(names(BRANCH), function(b) {
    w <- BRANCH[[b]](u); data.frame(part = "branch", shape = b, u = u, w = w/sum(w)) })))
kern$part_fraction <- ifelse(kern$part == "bole", BOLE_FRACTION, BRANCH_FRACTION)
write.csv(kern, out_path("scaling_shape_kernels.csv"), row.names=FALSE)

# The measured band profile, 0-2 m, as a ratio to each stem's own 2 m value: the
# forest's height response (a step function), not one of the extrapolation forms.
zb <- seq(0, 2, by = 0.02)
kk <- findInterval(pmin(pmax(zb*100, 50), 200), edges, rightmost.closed = TRUE)
wgt <- INV$A_band_m2
band_prof <- data.frame(series = "all stems (area-weighted)", z = zb,
  ratio = sapply(kk, function(k) sum(PROFI[, k]*wgt)/sum(f2*wgt)))
# also per species (the aggregate describes no individual species)
sp_big <- names(sort(table(INV$sp), decreasing = TRUE))
sp_big <- sp_big[sp_big != "SPECIES_OTHER"][1:6]
band_sp <- bind_rows(lapply(sp_big, function(sp) {
  i <- INV$sp == sp
  v <- colMeans(PROFI[i, , drop = FALSE]); f2i <- mean(f2[i])
  data.frame(series = sp, z = zb, ratio = sapply(kk, function(k) v[k]/f2i)) }))
band_prof <- bind_rows(band_prof, band_sp)
# absolute flux too, and a 1.25 m normalisation, so the figure can present all
# three framings rather than committing to "relative to 2 m" by default
k125 <- findInterval(125, edges, rightmost.closed = TRUE)
abs_ag <- sapply(kk, function(k) sum(PROFI[, k]*wgt)/sum(wgt))
band_prof$absolute <- NA_real_; band_prof$rel125 <- NA_real_
ia <- band_prof$series == "all stems (area-weighted)"
band_prof$absolute[ia] <- abs_ag
band_prof$rel125[ia]   <- abs_ag / (sum(PROFI[, k125]*wgt)/sum(wgt))
for (sp in sp_big) {
  i <- INV$sp == sp; v <- colMeans(PROFI[i, , drop = FALSE])
  j <- band_prof$series == sp
  band_prof$absolute[j] <- sapply(kk, function(k) v[k])
  band_prof$rel125[j]   <- sapply(kk, function(k) v[k]/v[k125])
}
write.csv(band_prof, out_path("scaling_band_profile.csv"), row.names=FALSE)

# the stand-mean ABSOLUTE 2 m flux, so the extrapolation forms can be shown in
# absolute units as well as as ratios
write.csv(data.frame(f2_stand_mean = sum(f2*wgt)/sum(wgt),
                     f125_stand_mean = sum(PROFI[, k125]*wgt)/sum(wgt)),
          out_path("scaling_flux_anchors.csv"), row.names=FALSE)
ag <- band_prof[band_prof$series == "all stems (area-weighted)", ]
cat(sprintf("  band profile 0-2 m: aggregate %.2f at base -> 1.00 at 2 m (steps at %s cm)\n",
            ag$ratio[1], paste(head(edges[-1], -1), collapse=", ")))
cat(sprintf("  per species at the stem base: %s\n", paste(sprintf("%s %.2f",
      sub(" .*","",sp_big), sapply(sp_big, function(sp)
        band_sp$ratio[band_sp$series==sp][1])), collapse="  ")))

# stand woody-area density against absolute height, per bole x branch x WAI
ar_prof <- list()
for (bn in names(BOLE)) for (rn in names(BRANCH)) {
  wb <- BOLE[[bn]](u); wb <- wb/sum(wb)
  wr <- BRANCH[[rn]](u); wr <- wr/sum(wr)
  Ashape <- outer(rep(1,n), wb)*BOLE_FRACTION + outer(rep(1,n), wr)*BRANCH_FRACTION
  Ashape <- Ashape * (pi*INV$dbh*INV$H/2)
  nb <- rowSums(!above); Aband <- (!above) * INV$A_band_m2 / pmax(nb,1)
  shp_ab <- Ashape*above; shp_ab <- shp_ab/pmax(rowSums(shp_ab),1e-12)
  for (wn in names(WAIS)) {
    A_above_tot <- WAIS[wn]*A_PLOT - sum(INV$A_band_m2)
    w_tree <- pi*INV$dbh*INV$H/2; w_tree <- w_tree/sum(w_tree)
    Aslab <- Aband + shp_ab * (w_tree * A_above_tot)
    k <- findInterval(as.vector(Z), ZB, rightmost.closed = TRUE)
    ok <- k >= 1 & k <= length(zmid)
    a <- tapply(as.vector(Aslab)[ok], k[ok], sum)
    v <- rep(0, length(zmid)); v[as.integer(names(a))] <- as.numeric(a)
    # also the flux-weighted area, sum_i f2_i * Aslab_i(z): stem area and 2 m flux are
    # correlated across stems, so this (not area x mean flux) integrates to the grid cell
    fa <- tapply((as.vector(Aslab) * rep(f2, times = ncol(Aslab)))[ok], k[ok], sum)
    w2 <- rep(0, length(zmid)); w2[as.integer(names(fa))] <- as.numeric(fa)
    ar_prof[[length(ar_prof)+1]] <- data.frame(
      bole = bn, branch = rn, WAI = wn, z = zmid,
      area_m2 = v, area_m2_per_m = v/diff(ZB)[1],
      fluxarea = w2, fluxarea_per_m = w2/diff(ZB)[1])
  }
}
ar_prof <- bind_rows(ar_prof)
write.csv(ar_prof, out_path("scaling_area_profiles.csv"), row.names=FALSE)
cat(sprintf("  exported shape profiles: %d flux rows, %d area rows (bins 0-%g m)\n",
            nrow(fx_prof), nrow(ar_prof), max(ZB)))

cat(sprintf("\n=== FULL GRID: %d combinations ===\n", nrow(R)))
cat(sprintf("  %d WAI x %d bole x %d branch x %d flux forms\n\n",
            length(WAIS), length(BOLE), length(BRANCH), length(FORMS)))

# The grid spans zero (linear_bounded_median), so ranges are reported as intervals and
# spreads in mg, not ratios. The uptake form is also reported separately: the only
# measurements above 2 m (the felled oak) are all positive.
CONTRA <- "linear_bounded_median"
Rp <- R[R$flux != CONTRA, ]
rng <- function(v) sprintf("%.1f to %.1f", min(v), max(v))
cat(sprintf("  whole-surface total, all %d      %s mg CH4 m-2 yr-1\n", nrow(R), rng(R$total_mg)))
cat(sprintf("  excluding the uptake form (%d)   %s\n", nrow(Rp), rng(Rp$total_mg)))
cat(sprintf("  measured band (invariant)        %.2f\n", unique(round(R$measured_mg,2))[1]))
cat(sprintf("  as %% of soil uptake, all         %s%%\n", rng(R$pct_of_soil)))
cat(sprintf("    excluding the uptake form      %s%%\n", rng(Rp$pct_of_soil)))
# share of the total that is extrapolated is only interpretable where the total is
# positive and larger than the measured band
ok <- R$total_mg > R$measured_mg
cat(sprintf("  extrapolated share (%d of %d combinations where it is defined)  %s%%\n",
            sum(ok), nrow(R), rng(R$pct_extrapolated[ok])))
cat(sprintf("  net budget                       %s   -- SINK in %d of %d\n",
            rng(R$net_mg), sum(R$net_mg < 0), nrow(R)))
cat(sprintf("  median total %.1f (IQR %.1f-%.1f); excluding uptake form %.1f (IQR %.1f-%.1f)\n",
            median(R$total_mg), quantile(R$total_mg,.25), quantile(R$total_mg,.75),
            median(Rp$total_mg), quantile(Rp$total_mg,.25), quantile(Rp$total_mg,.75)))

cat("\n--- which assumption moves the answer most? ---\n")
cat("    spread = max - min of the median total across that assumption's levels,\n")
cat("    in mg CH4 m-2 yr-1. A ratio would be undefined across the sign change.\n")
att <- lapply(c("WAI","bole","branch","flux"), function(v) {
  g <- R %>% group_by(.data[[v]]) %>% summarise(m = median(total_mg), .groups = "drop")
  data.frame(assumption = v, levels = nrow(g), lo = min(g$m), hi = max(g$m),
             spread = max(g$m) - min(g$m))
}) %>% bind_rows() %>% arrange(-spread)
for (i in seq_len(nrow(att))) cat(sprintf("  %-8s (%d levels)  median total %7.1f to %7.1f   spread %6.1f\n",
  att$assumption[i], att$levels[i], att$lo[i], att$hi[i], att$spread[i]))
cat(sprintf("\n  the flux form dominates the geometry by %.1fx (%.1f vs %.1f mg)\n",
            att$spread[1]/att$spread[2], att$spread[1], att$spread[2]))

cat("\n--- the hold case (constant flux), by WAI ---\n")
hc <- R %>% filter(flux=="constant") %>% group_by(WAI) %>%
  summarise(total=mean(total_mg), measured=mean(measured_mg),
            pct_extrap=mean(pct_extrapolated), pct_soil=mean(pct_of_soil), .groups="drop")
print(as.data.frame(hc), row.names=FALSE, digits=4)
cat("\n  (constant flux is identical across all 8 area shapes -- the vertical\n   distribution cancels algebraically, so only WAI matters in this row)\n")

# --- figure -------------------------------------------------------------------
# panel order from FORMS (a typed list could silently drop renamed forms)
stopifnot(setequal(unique(R$flux), FORMS))
R$flux <- factor(R$flux, levels=FORMS[order(-tapply(R$total_mg, R$flux, median)[FORMS])])
# signed log axis (asinh, log10 units, cofactor 2): keeps the negative totals
slog <- flux_asinh_trans(cofactor = 2, base10 = TRUE)
sbrk <- c(-10, -3, 0, 3, 10, 30, 100)
p1 <- ggplot(R, aes(flux, total_mg, fill=WAI)) +
  geom_boxplot(outlier.size=.5, linewidth=.3) +
  geom_hline(yintercept=0, linewidth=.25, colour="grey40") +
  scale_y_continuous(trans=slog, breaks=sbrk, labels=sbrk) +
  labs(title=sprintf("a  whole-surface budget across all %d combinations", nrow(R)),
       subtitle=sprintf("box spans the %d area-distribution shapes; colour = WAI source",
                        length(BRANCH)*length(BOLE)),
       x=NULL, y=expression(mg~CH[4]~m^-2~yr^-1), fill="WAI") +
  theme_bw(base_size=8) + theme(legend.position="bottom",
    legend.text=element_text(size=6), axis.text.x=element_text(angle=20,hjust=1))
p2 <- ggplot(R, aes(total_mg, pct_extrapolated, colour=flux, shape=WAI)) +
  geom_point(size=1.6, alpha=.8) +
  geom_vline(xintercept=0, linewidth=.25, colour="grey40") +
  scale_x_continuous(trans=slog, breaks=sbrk, labels=sbrk) +
  labs(title="b  how much of each estimate is extrapolated above 2 m",
       x=expression(mg~CH[4]~m^-2~yr^-1), y="% from above 2 m", colour=NULL, shape=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6),
                                legend.box="vertical", legend.key.height=unit(7,"pt"))
ggsave(out_path("fig_scaling_full_grid.png"), p1/p2 +
  plot_annotation(title="Full stem-CH4 scaling sensitivity",
    subtitle=sprintf("%d combinations; soil uptake %.0f mg CH4 m-2 yr-1 for reference",
                     nrow(R), SOIL_ANN),
    theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=9,height=8,dpi=200,bg="white")
cat("\nWritten: outputs/data/scaling_full_grid.csv and fig_scaling_full_grid.png\n")
