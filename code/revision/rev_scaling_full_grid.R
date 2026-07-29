#!/usr/bin/env Rscript
# ==============================================================================
# rev_scaling_full_grid.R
# ------------------------------------------------------------------------------
# THE FULL SCALING SENSITIVITY GRID.
#
# Every unconstrained assumption, crossed:
#
#   WAI (total woody area)   3 sourced values
#       1.50   Whittaker & Woodwell 1967, Brookhaven oak-pine (open, mixed conifer)
#       2.11   bottom-up from THIS inventory + W&W relations (taper 1.275, br:stem 3.35);
#              the full bottom-up range is 1.69-2.57
#       3.07   Gauci et al. 2024, temperate forest -- the value used in the submitted ms
#
#   BOLE SHAPE               3   cylinder / cone / paraboloid   (W&W: true stems lie
#                                between the conic and parabolic forms)
#   BRANCH SHAPE             4   uniform_all / uniform_top50 / uniform_crown / gaussian_75
#   FLUX ABOVE 2 m           6   Each has a stated mechanism and is identifiable from
#                                three heights. One constraint is imposed: no stem may
#                                increase substantially with height (no mechanism, no
#                                precedent). Everything else is allowed.
#       constant          no decline aloft; equals the asymptote-at-f(2m) case
#       exponential       exponential decay toward zero; the asymptote-at-0 case
#       power             gradual decay, different curvature, never reaches zero
#       linear_floored    linear decline truncated at zero
#       linear_bounded_median  linear decline continuing THROUGH zero into net
#                         uptake, bounded at the median detected stem uptake
#       exp_band_slope per-stem slope learned by the RF from species, moisture,
#                         temperature and DBH, capped at that stem's own 2 m value
#
#   An asymptotic ("decline then plateau") form is NOT carried separately: it is
#   algebraically phi*constant + (1-phi)*exponential, so the family is spanned by the
#   two endpoints already present. See rev_asymptote_estimation.R -- the asymptote is
#   estimable only after aggregating to species or moisture bins, only converges for a
#   minority of groups, and where estimable sits at 0.83-0.97 of the 2 m flux, i.e.
#   against the 'constant' end of that continuum.
#
#   = 6 flux x 5 WAI x 4 branch x 2 bole = 240 combinations. Do not hardcode this
#   count anywhere: it is length(FORMS)*length(WAIS)*length(BRANCH)*length(BOLE),
#   and the stale "216" that stood here outlived three changes to the axes.
#
# STRUCTURE OF EACH ESTIMATE. The budget is split into two parts that are reported
# separately, because they have completely different evidentiary status:
#
#   MEASURED    woody area at z <= 2 m, flux from the RF predicted at that height for
#               that tree (species, DBH, moisture, temperature). Interpolation inside
#               the training range. Below 0.5 m the 0.5 m value is held (conservative;
#               tested earlier and worth <5%).
#
#   EXTRAPOLATED  woody area at z > 2 m, flux = each tree's RF value at 2 m times the
#               chosen form's shape ratio. NOTHING here is measured at this site.
#
# The extrapolated fraction is the honest headline uncertainty: it is typically ~85-90%
# of the whole-surface number.
#
# Output: outputs/revision/scaling_full_grid.csv
#         outputs/revision/fig_scaling_full_grid.png
# ==============================================================================

suppressMessages({library(ranger);library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
source("code/revision/rev_geometry.R")
source("code/revision/rev_species_levels.R")
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR <- rf_workflow_data$PLACEHOLDER_DRIVERS
d <- tree_train_complete; trained <- sort(unique(as.character(d$species_clean)))
CONV <- 86400*365.25*16e-6

# Inputs are the canonical products, not re-derivations. The measured band and the
# soil term were previously recomputed here and drifted from the budget: the grid
# reported a 5.7874 mg band while canonical_budget.csv said 5.41, and Figure 9
# stacked one on the other in a single waterfall. The soil term was the hardcoded
# constant -543.88. Both are now read.
need <- c(inv="outputs/tables/inventory_stems.csv",
          tree="outputs/tables/tree_flux_predictions.csv",
          prof="outputs/tables/tree_band_profile.csv",
          pedge="outputs/tables/tree_band_profile_edges.csv",
          bud="outputs/revision/canonical_budget.csv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing inputs:\n  ", paste(miss, collapse="\n  "),
  "\nrun rev_inventory_build.R, rev_predict_tree_flux_current.R, rev_predict_soil_surface.R, rev_budget_canonical.R first")
B   <- read.csv(need[["bud"]], stringsAsFactors=FALSE)
val <- function(q){v <- B$value[B$quantity==q]; if(!length(v)) stop("canonical_budget.csv lacks ",q); v}
SOIL_ANN  <- val("soil_mg_m2_yr")
MEAS_BAND <- val("tree_measured_mg_m2_yr")   # invariant across the grid, by construction
A_PLOT    <- STAND_AREA_M2                   # censused ground, not the nominal square

TRP <- read.csv(need[["tree"]], stringsAsFactors=FALSE)
INV <- read.csv(need[["inv"]],  stringsAsFactors=FALSE)
INV <- INV[INV$in_stand, ]; TRP <- TRP[TRP$in_stand, ]
stopifnot(nrow(INV)==nrow(TRP))
INV$dbh <- INV$dbh_m                          # already unit-checked and typo-repaired
INV$A_band_m2 <- TRP$A_stem_m2                # measured-band area, Kalmia at 0.75 m
INV <- INV %>% group_by(species) %>%
  ungroup() %>% as.data.frame()
INV$sp <- species_to_model_level(INV$species, trained)
GY <- c("Pinus strobus","Tsuga canadensis")
DA <- quantile(INV$dbh[INV$dbh>0.10],.95,na.rm=TRUE)
INV$H <- 1.37 + ifelse(INV$species %in% GY,(25-1.37)/DA^0.60,(25-1.37)/DA^0.53) *
         INV$dbh^ifelse(INV$species %in% GY,0.60,0.53)     # no shrub cap, by decision
n <- nrow(INV)

# --- per-tree band profile: READ, not re-derived ------------------------------
# This block used to re-run the forest here to rebuild the profile and the 2 m
# anchor. That re-derivation drifted from the canonical band in two ways at once:
# it drove the model from the CAMPAIGN monthly means (DRIVERS$soil_moisture_at_tree
# with one plot-mean value per month, December substituted from the soil collars,
# January and March absent and filled by approx(rule = 2), and campaign soil
# temperature rather than the climatology), and it applied NO per-species
# calibration. Since the extrapolated term hangs entirely off the 2 m anchor, that
# put 74% of the headline scenario on a different basis from the 26% read from
# canonical_budget.csv above -- inside a single reported total.
# rev_predict_tree_flux_current.R now exports the calibrated profile it already
# computes, on the climatology drivers and the per-stem moisture surface, and it is
# consumed here. One producer per quantity: the grid re-derives nothing.
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
# BRANCH PLACEMENT. None of these is empirically supported here; they bracket
# plausible crown geometries. gaussian_50 was added after the profile figure made
# the set's shape visible: uniform_crown puts every branch above 0.60H with a hard
# cutoff and nothing below, which no real crown does, and without a second Gaussian
# the smooth options were represented only by one centred high at 0.75H. The pair
# of Gaussians now brackets a deep crown and a shallow one.
# uniform_crown (all branch area above 0.60H) is DROPPED. It was a near-duplicate
# of uniform_top50: cutoffs at 0.60H and 0.50H gave median totals of 16.4 and 16.9
# mg, a 3% difference, so the two together consumed two of five levels while
# spanning almost none of the range. The remaining four bracket properly -- area
# low (uniform_all), a hard cutoff (uniform_top50), and two Gaussian crowns at
# different depths (0.50H, 0.75H).
BRANCH <- list(uniform_all=function(u) rep(1,length(u)),
  uniform_top50=function(u) as.numeric(u>=.50),
  gaussian_50=function(u) dnorm(u,.50,.18),
  gaussian_75=function(u) dnorm(u,.75,.15))
# WAI SCENARIOS. The bottom-up value is OURS and is labelled as such. It must not
# be presented as jointly supported by Whittaker & Woodwell's climax-forest value
# of 2.0-2.2: our estimate is built FROM their taper correction and branch:stem
# ratios applied to this inventory, so landing near their published figure is
# partly circular and is corroboration only of the arithmetic, not of the value.
# That belongs in the Discussion, not in a scenario label. The two literature
# values bound the range from either side.
# The three bottom-up values are READ from their producer, not retyped. They were
# hardcoded here as 1.69 / 2.11 / 2.57, which were computed by
# rev_wai_bottomup_and_rf_interactions.R on the PRE-REBUILD inventory (7,010 stems,
# the whole bytag survey missing) divided by the nominal 40,000 m2 square instead of
# the 37,200 m2 stand. Both errors pushed WAI down; the correct trio is 2.23 / 2.82 /
# 3.40. Because the literals lived here, fixing the producer alone would not have
# reached the headline -- and WAI multiplies the extrapolated term, which is ~70% of it.
# The two literature bounds stay literal: they are citations, not derived quantities.
WAI_BU <- local({
  f <- file.path(outdir, "wai_bottomup.csv")
  if (!file.exists(f))
    stop("missing ", f, " -- run: Rscript code/revision/rev_wai_bottomup_and_rf_interactions.R")
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
# NEGATIVE BOUND for the linear form, from the detected stem uptakes surviving Tier 1
# QC (MDF at 90% confidence, precision per field period; n = 89).
#
#   USED:  -0.0262, the median detected stem uptake. This is 2.5% of the median
#          detected soil uptake (-1.057): any stem sink implied here is ~40x weaker
#          than the soil sink.
#
#   The 5th percentile (-0.0984) was tested and REJECTED. Uptake prevalence tracks
#   instrument precision, not any biological covariate:
#       Height+molecular  sigma 1.200 ppb   0.0% uptake  (none detected)
#       Cross-species     sigma 1.725 ppb  10.4% uptake   deepest -0.3912   median r2 0.38
#       Monthly survey    sigma 2.181 ppb   7.0% uptake   deepest -0.1757   median r2 0.21
#   The 5th percentile is 3x deeper than anything seen under the best instrument
#   conditions and is drawn entirely from the two noisy campaigns' tails, where fits
#   are poor. Because the floor applies to ~90% of woody surface, adopting it would
#   make an extreme tail value the modal assumption for the whole canopy rather than a
#   mild outer bound.
#
#   The campaign with the lowest instrument noise detects NO stem uptake at all, so the
#   bound rests entirely on the two noisier periods. It is further contradicted by the
#   climbed tree (Felled black oak), where all four measurements above 2 m are positive
#   and each exceeds this bound in magnitude.
NEG_MEDIAN <- -0.0262
FORMS <- c("constant","exponential","power","linear_floored",
           "linear_bounded_median","exp_band_slope")

NU <- 90; u <- (seq_len(NU)-.5)/NU
Z  <- outer(INV$H, u)                       # n x NU absolute heights
above <- Z > 2
# per-tree flux inside the band, floored below 0.5 m
Zb <- pmin(pmax(Z,0.5),2)
Fband <- matrix(NA_real_, n, NU)
for (i in seq_len(n)) Fband[i,] <- prof_at(i, Zb[i,])

# Shape ratios relative to each stem's own 2 m flux. linear_unfloored is the only form
# permitted to go negative, representing a decline that continues through zero into net
# uptake aloft -- a behaviour reported in the literature and implied by ~15% of the
# converging stems here, but never directly observed in our dataset.
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
# exp_band_slope. Named for what it is rather than for the model that produced it:
# a log-linear (exponential) decay whose rate is fitted per stem to the four RF
# band values, then projected above 2 m and capped at 1 so no stem exceeds its own
# 2 m flux. It is the same functional FAMILY as `exponential`; what differs is that
# the rate comes from each stem's own within-band profile rather than from one fit
# to the campaign means.
#
# It is a deliberately conservative scenario and its limits should be stated. The
# fit spans 1.13 m (interval midpoints 0.685-1.815 m) and is projected to 25 m;
# median R2 of the four-point fit is 0.559, because a straight line in log space
# fits a non-monotonic step poorly; and the 50 cm value supplies essentially all
# of the decline -- refitting above 88 cm gives a median slope of +0.057 against
# -0.672 with it. The consequence is that 97% of stems fall below 1% of their 2 m
# flux by 25 m, against 0.82x measured on the one climbed tree at 10 m. It is
# retained because a single climbed stem cannot exclude decline in other species
# or conditions, but it is a lower bound on emission rather than a central case.
# Uncapped, the 31 stems (0.4%) with a slightly positive fitted slope reach up to 100x
# their 2 m value at canopy height -- the one behaviour excluded a priori.
#
# NOTE, and it matters for reading the grid: this rate is fitted to the RF's own
# within-band profile, so it is a MODEL-DERIVED shape, not a data-derived one.
# Dropping dbh_within_z from the forest moved it: the median slope became shallower
# (-0.914 to -0.672) while the area-weighted fraction surviving to 25 m FELL, from
# 0.108 to 0.030, because the shallow-sloped stems that carried most of that
# survival lost their shallow slopes. That is why the extrapolated term falls
# further (-53%) than the measured band (-32%) between the two model versions.
Rmat$exp_band_slope <- pmin(exp(rf_slope * (Z-2)), 1)

res <- list()
for (bn in names(BOLE)) for (rn in names(BRANCH)) {
  wb <- BOLE[[bn]](u); wb <- wb/sum(wb)
  wr <- BRANCH[[rn]](u); wr <- wr/sum(wr)
  # Shape of area with height: bole and branch in W&W proportion (branch:stem 3.35).
  Ashape <- outer(rep(1,n), wb)*(1/4.35) + outer(rep(1,n), wr)*(3.35/4.35)
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
    # the measured band is the canonical value, not a re-derivation (see header)
    meas <- MEAS_BAND
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
write.csv(R, file.path(outdir,"scaling_full_grid.csv"), row.names=FALSE)

# ---- THE HEADLINE SCENARIO, named once ---------------------------------------
# One combination is quoted in the text and drawn in the figures, so it is defined
# here rather than chosen ad hoc in each script:
#   flux    exp_band_slope  a decay rate taken from the model's own within-band
#                           profile rather than from an assumed shape.
#                           IT IS NOT "the lowest of the positive forms" -- that
#                           claim stood here and is false. At this geometry it is
#                           third of five positive forms (constant 55.9, power 25.5,
#                           exp_band_slope 15.7, exponential 7.7, linear_floored
#                           5.1 mg), so it cannot be defended as the conservative
#                           choice. It is also the form the project's own evidence
#                           argues against quoting: leave-one-height-out CV ranks
#                           this family LAST on upward extrapolation (exp_zero RMSE
#                           6.899 vs constant 0.425), rev_rf_height_extrapolation.R
#                           concludes in its own text that it should be reported
#                           "as a low-end bounding scenario, not as the best
#                           estimate", and the only direct measurements above 2 m
#                           (the climbed oak, 0.82x at 10 m) contradict it -- this
#                           form puts the median stem at 0.005x by 10 m.
#                           RETAINED AS THE NAMED CELL ONLY so that figures and text
#                           point at one place. THE RANGE IS THE RESULT. Whether
#                           this stays the quoted cell is an open decision; see
#                           code/revision/notes/AUDIT_2026-07-29.md, blocker 4.
#   WAI     2.11            our bottom-up estimate for THIS stand, not a
#                           literature value imported from another forest
#   branch  gaussian_75     a crown centred at 0.75H, which is what a closed-canopy
#                           dominant looks like after self-pruning; the hard-cutoff
#                           and low-centred alternatives are kept as scenarios
#   bole    cone            irrelevant (0.3 mg across the whole grid), so the
#                           simpler of the two W&W forms
# The RANGE, not this cell, is the result; this is the cell to point at.
# Select the headline WAI by its MEANING, not by its digits. This was
# grep("2.11", names(WAIS)) -- a numeric string match, so the moment the bottom-up
# estimate was corrected from 2.11 to 2.82 the selector silently returned NA.
HEADLINE <- list(flux = "exp_band_slope", branch = "gaussian_75", bole = "cone",
                 WAI = grep("bottom-up, this stand", names(WAIS), value = TRUE)[1])
stopifnot(length(HEADLINE$WAI) == 1, !is.na(HEADLINE$WAI))
hr <- R[R$flux == HEADLINE$flux & R$branch == HEADLINE$branch &
        R$bole == HEADLINE$bole & R$WAI == HEADLINE$WAI, ]
write.csv(cbind(as.data.frame(HEADLINE), hr[, c("measured_mg","extrapolated_mg",
          "total_mg","pct_extrapolated","pct_of_soil","net_mg")]),
          file.path(outdir,"scaling_headline.csv"), row.names=FALSE)
cat(sprintf("\n=== HEADLINE SCENARIO ===\n  %s | %s | %s | %s\n",
            HEADLINE$flux, HEADLINE$WAI, HEADLINE$branch, HEADLINE$bole))
cat(sprintf("  total %.1f mg CH4 m-2 yr-1 (%.1f%% of soil uptake, %.0f%% extrapolated)\n",
            hr$total_mg, hr$pct_of_soil, hr$pct_extrapolated))

# ---- export the SHAPES themselves, for rev_fig_scaling_profiles.R ------------
# The profile figure must show exactly what this grid computed, so the curves are
# written out here rather than re-derived from the same formulae elsewhere: a
# second implementation is a second thing that can drift.
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
write.csv(fx_prof, file.path(outdir,"scaling_flux_shapes.csv"), row.names=FALSE)

# RAW SHAPE KERNELS, before any application to the inventory. The stand profiles
# below are the sum over thousands of stems of differing height, so trees drop out
# of successive slabs and the curves acquire steps and bumps that belong to the
# INVENTORY, not to the shapes. Exporting the kernels lets the profile figure show
# the forms themselves, on a single tree, which is what a reader needs to
# understand them.
kern <- bind_rows(
  bind_rows(lapply(names(BOLE), function(b) {
    w <- BOLE[[b]](u); data.frame(part = "bole", shape = b, u = u, w = w/sum(w)) })),
  bind_rows(lapply(names(BRANCH), function(b) {
    w <- BRANCH[[b]](u); data.frame(part = "branch", shape = b, u = u, w = w/sum(w)) })))
# W&W branch:stem 3.35, so the bole carries 1/4.35 of woody area and branches 3.35/4.35
kern$part_fraction <- ifelse(kern$part == "bole", 1/4.35, 3.35/4.35)
write.csv(kern, file.path(outdir,"scaling_shape_kernels.csv"), row.names=FALSE)

# THE MEASURED BAND PROFILE, 0-2 m, as a ratio to each stem's own 2 m value. This
# is what the model actually uses below 2 m -- the RF's own height response, a step
# function with breakpoints where the forest splits -- and it is NOT one of the six
# extrapolation forms. The profile figure needs it to show the full 0-26 m picture.
zb <- seq(0, 2, by = 0.02)
kk <- findInterval(pmin(pmax(zb*100, 50), 200), edges, rightmost.closed = TRUE)
wgt <- INV$A_band_m2
band_prof <- data.frame(series = "all stems (area-weighted)", z = zb,
  ratio = sapply(kk, function(k) sum(PROFI[, k]*wgt)/sum(f2*wgt)))
# ALSO per species. The aggregate is what the budget uses, but it describes no
# individual species: the decline from 50 cm to 2 m runs from 1.6x in Quercus
# rubra to 3.5x in Pinus strobus. Reporting only the pooled curve would present a
# stand-level average as if it were a biological profile.
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
write.csv(band_prof, file.path(outdir,"scaling_band_profile.csv"), row.names=FALSE)

# the stand-mean ABSOLUTE 2 m flux, so the extrapolation forms can be shown in
# absolute units as well as as ratios
write.csv(data.frame(f2_stand_mean = sum(f2*wgt)/sum(wgt),
                     f125_stand_mean = sum(PROFI[, k125]*wgt)/sum(wgt)),
          file.path(outdir,"scaling_flux_anchors.csv"), row.names=FALSE)
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
  Ashape <- outer(rep(1,n), wb)*(1/4.35) + outer(rep(1,n), wr)*(3.35/4.35)
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
    # ALSO the flux-weighted area, sum_i f2_i * Aslab_i(z). The plain area profile
    # times a mean 2 m flux does NOT reproduce the grid, because stem area and stem
    # flux are correlated across the inventory -- large stems are not average ones.
    # Exporting this makes the profile figure integrate to exactly the grid cell.
    fa <- tapply((as.vector(Aslab) * rep(f2, times = ncol(Aslab)))[ok], k[ok], sum)
    w2 <- rep(0, length(zmid)); w2[as.integer(names(fa))] <- as.numeric(fa)
    ar_prof[[length(ar_prof)+1]] <- data.frame(
      bole = bn, branch = rn, WAI = wn, z = zmid,
      area_m2 = v, area_m2_per_m = v/diff(ZB)[1],
      fluxarea = w2, fluxarea_per_m = w2/diff(ZB)[1])
  }
}
ar_prof <- bind_rows(ar_prof)
write.csv(ar_prof, file.path(outdir,"scaling_area_profiles.csv"), row.names=FALSE)
cat(sprintf("  exported shape profiles: %d flux rows, %d area rows (bins 0-%g m)\n",
            nrow(fx_prof), nrow(ar_prof), max(ZB)))

cat(sprintf("\n=== FULL GRID: %d combinations ===\n", nrow(R)))
cat(sprintf("  %d WAI x %d bole x %d branch x %d flux forms\n\n",
            length(WAIS), length(BOLE), length(BRANCH), length(FORMS)))

# The grid now spans zero, because linear_bounded_median lets flux decline through
# zero into uptake aloft. Fold-changes and shares of a total are undefined across a
# sign change, so ranges are reported as intervals and spreads in mg, never as
# ratios. The uptake form is also reported separately: it is CONTRADICTED by the
# only direct evidence above 2 m (all four climbed-tree measurements there are
# positive, each larger in magnitude than the sink it assumes), so it bounds what
# the data cannot exclude rather than describing a likely outcome.
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

cat("\n--- headline candidate: constant flux, by WAI ---\n")
hc <- R %>% filter(flux=="constant") %>% group_by(WAI) %>%
  summarise(total=mean(total_mg), measured=mean(measured_mg),
            pct_extrap=mean(pct_extrapolated), pct_soil=mean(pct_of_soil), .groups="drop")
print(as.data.frame(hc), row.names=FALSE, digits=4)
cat("\n  (constant flux is identical across all 8 area shapes -- the vertical\n   distribution cancels algebraically, so only WAI matters in this row)\n")

# --- figure -------------------------------------------------------------------
# Order the panel by the forms that ACTUALLY exist. This was a hardcoded list of
# six names, three of which (asymptote50, rf_learned, log_linear) were retired or
# renamed: factor() silently mapped them to NA, so 120 of the 240 rows -- including
# exp_band_slope, the headline form -- vanished from the figure without a warning.
# Derive the order from FORMS so it can never fall out of step again.
stopifnot(setequal(unique(R$flux), FORMS))
R$flux <- factor(R$flux, levels=FORMS[order(-tapply(R$total_mg, R$flux, median)[FORMS])])
# A log10 axis cannot render the uptake form: linear_bounded_median goes negative, so
# log10 sent 40 of the 240 rows to NaN and ggplot dropped them with only a warning --
# the same silent-loss failure as the factor-level bug above, and it hid precisely the
# scenario the reviewers asked about. Use a signed log (asinh, base 10, unit cofactor)
# so sign is preserved and the whole grid is visible.
slog <- scales::trans_new("signed_log10",
          transform = function(x) asinh(x/2)/log(10),
          inverse   = function(x) 2*sinh(x*log(10)))
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
ggsave(file.path(outdir,"fig_scaling_full_grid.png"), p1/p2 +
  plot_annotation(title="Full stem-CH4 scaling sensitivity",
    subtitle=sprintf("%d combinations; soil uptake %.0f mg CH4 m-2 yr-1 for reference",
                     nrow(R), SOIL_ANN),
    theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=9,height=8,dpi=200,bg="white")
cat("\nWritten: outputs/revision/scaling_full_grid.csv and fig_scaling_full_grid.png\n")
