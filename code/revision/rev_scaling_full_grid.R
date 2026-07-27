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
#       linear_unfloored  linear decline continuing THROUGH zero into net uptake
#       rf_learned_capped per-stem slope learned by the RF from species, moisture,
#                         temperature and DBH, capped at that stem's own 2 m value
#
#   An asymptotic ("decline then plateau") form is NOT carried separately: it is
#   algebraically phi*constant + (1-phi)*exponential, so the family is spanned by the
#   two endpoints already present. See rev_asymptote_estimation.R -- the asymptote is
#   estimable only after aggregating to species or moisture bins, only converges for a
#   minority of groups, and where estimable sits at 0.83-0.97 of the 2 m flux, i.e.
#   against the 'constant' end of that continuum.
#
#   = 216 combinations.
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
  mutate(dbh_within_z=if(n()>1&&sd(dbh,na.rm=TRUE)>0) as.numeric(scale(dbh)) else 0) %>%
  ungroup() %>% as.data.frame()
INV$sp <- ifelse(!is.na(INV$species) & INV$species %in% trained, INV$species, "SPECIES_OTHER")
GY <- c("Pinus strobus","Tsuga canadensis")
DA <- quantile(INV$dbh[INV$dbh>0.10],.95,na.rm=TRUE)
INV$H <- 1.37 + ifelse(INV$species %in% GY,(25-1.37)/DA^0.60,(25-1.37)/DA^0.53) *
         INV$dbh^ifelse(INV$species %in% GY,0.60,0.53)     # no shrub cap, by decision
n <- nrow(INV)

# --- per-tree RF profile through the measured band ---------------------------
mm <- d %>% group_by(month) %>% summarise(m=mean(soil_moisture_at_tree,na.rm=TRUE),.groups="drop")
sm <- soil_train_complete %>% group_by(month) %>% summarise(ms=mean(soil_moisture_at_site,na.rm=TRUE),.groups="drop")
DR <- DR %>% left_join(mm,"month") %>% left_join(sm,"month") %>% mutate(m=ifelse(is.finite(m),m,ms))
fl <- function(v,mo) approx(mo[is.finite(v)],v[is.finite(v)],xout=mo,rule=2)$y
DR$m <- fl(DR$m,DR$month); DR$soil_temp_C_mean <- fl(DR$soil_temp_C_mean,DR$month)
# The RF is a STEP function in height: trained at 50/125/200 cm only, so it takes
# 3-4 distinct values across the band and the transitions sit at shared thresholds.
# Detect them, then evaluate once per interval -- exact, and ~40x fewer model calls
# than the 12.5 cm grid this used before (which also landed mid-step).
pred_at <- function(h, mo, idx=seq_len(n)) predict(TreeRF, data.frame(
    species=factor(INV$sp[idx],levels=trained), dbh_m=INV$dbh[idx],
    dbh_within_z=INV$dbh_within_z[idx], soil_moisture_at_tree=DR$m[mo],
    soil_temp_C_mean=DR$soil_temp_C_mean[mo], air_temp_C_mean=DR$air_temp_C_mean[mo],
    height_cm=h), num.threads=1)$predictions
samp <- sample(n, min(400,n)); Hs <- 50:200
Ps <- vapply(Hs, function(h) pred_at(h,7,samp), numeric(length(samp)))
brk <- Hs[-1][colSums(abs(Ps[,-1,drop=FALSE]-Ps[,-ncol(Ps),drop=FALSE]))>0]
edges <- c(50, brk, 201); K <- length(edges)-1
cat(sprintf("height steps at %s cm -> %d intervals; %d stems x %d x 12 model calls\n",
            paste(brk,collapse=", "), K, n, K))
FI <- array(NA_real_,c(n,K,12))
for (mo in 1:12) for (k in seq_len(K)) FI[,k,mo] <- pred_at(edges[k], mo)
PROFI <- apply(FI,c(1,2),mean)              # n x K, annual mean per interval
zh <- (head(edges,-1) + pmin(tail(edges,-1),201) - 1)/2/100   # interval midpoints, m
f2 <- PROFI[,K]                             # 2 m anchor for the above-band forms
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
BRANCH <- list(uniform_all=function(u) rep(1,length(u)),
  uniform_top50=function(u) as.numeric(u>=.50), uniform_crown=function(u) as.numeric(u>=.60),
  gaussian_75=function(u) dnorm(u,.75,.15))
WAIS <- c(`1.50 W&W Brookhaven`=1.50, `1.69 bottom-up low`=1.69,
          `2.11 bottom-up mid / W&W climax`=2.11, `2.57 bottom-up high`=2.57,
          `3.07 Gauci 2024`=3.07)
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
           "linear_bounded_median","rf_learned_capped")

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

SHARED <- setdiff(FORMS,"rf_learned_capped")
Rmat <- lapply(SHARED, function(fm) matrix(shape_ratio(fm, as.vector(Z)), nrow=n))
names(Rmat) <- SHARED
# rf_learned_capped: per-stem slope, CAPPED at 1 so no stem exceeds its own 2 m flux.
# Uncapped, the 98 stems (1.4%) with a slightly positive fitted slope reach up to 303x
# their 2 m value at canopy height -- the one behaviour excluded a priori.
Rmat$rf_learned_capped <- pmin(exp(rf_slope * (Z-2)), 1)

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

cat(sprintf("\n=== FULL GRID: %d combinations ===\n", nrow(R)))
cat(sprintf("  %d WAI x %d bole x %d branch x %d flux forms\n\n",
            length(WAIS), length(BOLE), length(BRANCH), length(FORMS)))
cat(sprintf("  whole-surface total   %8.1f  to %8.1f  mg CH4 m-2 yr-1  (%.1f-fold)\n",
            min(R$total_mg), max(R$total_mg), max(R$total_mg)/min(R$total_mg)))
cat(sprintf("  measured band only    %8.2f  to %8.2f\n", min(R$measured_mg), max(R$measured_mg)))
cat(sprintf("  extrapolated share    %8.1f%% to %8.1f%%\n",
            min(R$pct_extrapolated), max(R$pct_extrapolated)))
cat(sprintf("  as %% of soil uptake   %8.1f%% to %8.1f%%\n", min(R$pct_of_soil), max(R$pct_of_soil)))
cat(sprintf("  net budget            %8.1f  to %8.1f   -- SINK in %d of %d\n",
            min(R$net_mg), max(R$net_mg), sum(R$net_mg<0), nrow(R)))
cat(sprintf("  median total          %8.1f   IQR %.1f - %.1f\n", median(R$total_mg),
            quantile(R$total_mg,.25), quantile(R$total_mg,.75)))

cat("\n--- variance attribution: which assumption moves the answer most? ---\n")
for (v in c("WAI","bole","branch","flux")) {
  g <- R %>% group_by(.data[[v]]) %>% summarise(m=median(total_mg),.groups="drop")
  cat(sprintf("  %-8s median total ranges %7.1f to %7.1f   (%.2fx)\n",
              v, min(g$m), max(g$m), max(g$m)/min(g$m)))
}

cat("\n--- headline candidate: constant flux, by WAI ---\n")
hc <- R %>% filter(flux=="constant") %>% group_by(WAI) %>%
  summarise(total=mean(total_mg), measured=mean(measured_mg),
            pct_extrap=mean(pct_extrapolated), pct_soil=mean(pct_of_soil), .groups="drop")
print(as.data.frame(hc), row.names=FALSE, digits=4)
cat("\n  (constant flux is identical across all 12 area shapes -- the vertical\n   distribution cancels algebraically, so only WAI matters in this row)\n")

# --- figure -------------------------------------------------------------------
R$flux <- factor(R$flux, levels=c("constant","asymptote50","power","rf_learned",
                                  "log_linear","linear_floored"))
p1 <- ggplot(R, aes(flux, total_mg, fill=WAI)) +
  geom_boxplot(outlier.size=.5, linewidth=.3) +
  scale_y_continuous(trans="log10") +
  labs(title="a  whole-surface budget across all 216 combinations",
       subtitle="box spans the 12 area-distribution shapes; colour = WAI source",
       x=NULL, y=expression(mg~CH[4]~m^-2~yr^-1), fill="WAI") +
  theme_bw(base_size=8) + theme(legend.position="bottom",
    legend.text=element_text(size=6), axis.text.x=element_text(angle=20,hjust=1))
p2 <- ggplot(R, aes(total_mg, pct_extrapolated, colour=flux, shape=WAI)) +
  geom_point(size=1.6, alpha=.8) + scale_x_continuous(trans="log10") +
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
