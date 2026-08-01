source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# rev_wai_bottomup_and_rf_interactions.R
# ------------------------------------------------------------------------------
# Four questions, with sources.
#
#   PART 1  Height allometry: per-tree heights from DBH, from a published exponent.
#   PART 2  Bottom-up WAI from our own inventory using Whittaker & Woodwell (1967)
#           surface-area relations, compared with their published forest values and
#           with the Gauci et al. (2024) WAI of 3.07.
#   PART 3  Height-flux forms: what did power / exponential / linear / asymptotic
#           actually buy us, and can three heights distinguish them at all?
#   PART 4  Does the RF need explicit interactions? Tested, not asserted.
#
# ------------------------------------------------------------------------------
# SOURCES
# ------------------------------------------------------------------------------
# Whittaker RH, Woodwell GM (1967) Surface area relations of woody plants and forest
#   communities. Am J Bot 54:931-939.
#     - conic stem surface  S_c = pi * d * h / 2
#     - true stem surface / conic estimate = 1.20-1.35  (Table 1:6)
#     - Brookhaven oak-pine (open):  stem/ground 0.30, branch/ground ~1.2, total ~1.5
#     - Climax temperate deciduous (Great Smoky Mtns): stem/ground 0.5-0.6,
#       branch/ground 1.5-1.6, TOTAL BARK/GROUND = 2.0-2.2
#     - leaf/ground 4.0-6.0 in closed-canopy deciduous forest
#
# Hulshof CM, Swenson NG, Weiser MD (2015) Tree height-diameter allometry across the
#   United States. Ecol Evol 5:1193-1204.  (FIA, 2.98 million stems, 293 species)
#     - power-law scaling H ~ D^b, b = 0.53 angiosperms, 0.60 gymnosperms
#
# Yale-Myers Forest stand structure: mature stands 80-100 yr, closed canopy 20-30 m,
#   LAI 4.7-6.0. NOTE the LAI agrees with Whittaker & Woodwell's 4.0-6.0 for closed
#   deciduous forest, an independent indication the site is comparable to theirs.
#
# Gauci V et al. (2024) Nature 631:796-800.  WAI = 3.07 m2 m-2 for temperate forest.
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(ranger)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
source("code/lib/rev_geometry.R")
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
# CANONICAL stems and ground. This read rf_workflow_data$PLACEHOLDER_INVENTORY over
# A_PLOT <- 200*200: fg19 only, so the whole bytag survey (961 stems) was missing,
# 7,010 stems instead of 8,006, divided by the nominal 40,000 m2 square instead of
# the 37,200 m2 censused stand. Both errors push the WAI DOWN, and this script is the
# sole producer of the bottom-up WAI that rev_scaling_full_grid.R consumes, where it
# multiplies ~70% of the headline estimate. The trio it emitted (1.69/2.11/2.57) is
# 32% low; on the canonical inventory the same recipe gives 2.23/2.82/3.40.
INV <- canonical_inventory()
A_PLOT <- STAND_AREA_M2
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# fix_dbh() removed: it was the over-correcting species-specific DBH repair stack
# that rev_inventory_build.R replaced with one documented rule. canonical_inventory()
# returns diameters already unit-checked and typo-repaired, and H from the single
# allometry in rev_geometry.R.

# ==============================================================================
rule("PART 1  HEIGHT ALLOMETRY (per stem, not one canopy height)")
# ==============================================================================
# H = 1.37 + a * D^b, with b from Hulshof et al. 2015 (FIA, US-wide). 'a' is anchored
# so that a large canopy stem reaches the observed Yale-Myers canopy height: this uses
# the published EXPONENT (shape) with a site-specific SCALE, rather than importing
# both from a region that may differ.
GYMNO <- c("Pinus strobus","Tsuga canadensis")
INV$b_exp   <- ifelse(INV$species %in% GYMNO, 0.60, 0.53)
D_ANCHOR <- quantile(INV$dbh[INV$dbh > 0.10], 0.95, na.rm=TRUE)   # large canopy stem
cat(sprintf("\n  anchor DBH (95th pct of stems >10 cm) = %.1f cm\n", 100*D_ANCHOR))
cat(sprintf("  %-22s %8s %8s %8s %8s\n","canopy anchor H_c","a(ang)","median H","mean H","max H"))
HEIGHT_SCENARIOS <- c(20,25,30)
for (Hc in HEIGHT_SCENARIOS) {
  a_ang <- (Hc-1.37)/D_ANCHOR^0.53
  H <- 1.37 + ifelse(INV$species %in% GYMNO,
                     (Hc-1.37)/D_ANCHOR^0.60, a_ang) * INV$dbh^INV$b_exp
  cat(sprintf("  %-22s %8.2f %8.1f %8.1f %8.1f\n",
              sprintf("H_c = %d m", Hc), a_ang, median(H), mean(H), max(H)))
}
# One allometry, defined in rev_geometry.R. It was copy-pasted verbatim into five
# scripts, and because its scale coefficient is anchored on the 95th percentile of
# whichever stem list the script loaded, the duplicates silently disagreed.
H_CANOPY <- CANOPY_H_M
INV$H <- stem_height_m(INV$dbh_m, INV$species, canopy_h = H_CANOPY)
cat(sprintf("
  Using H_c = %d m. Because H scales as D^0.53, SUBCANOPY STEMS GET SUBCANOPY HEIGHTS:
  median stem height %.1f m against a %d m canopy -- the point you raised. A single
  canopy height applied to all %d stems would overstate total woody area by %.0f%%.
", H_CANOPY, median(INV$H), H_CANOPY, nrow(INV),
   100*(sum(pi*INV$dbh*H_CANOPY)/sum(pi*INV$dbh*INV$H) - 1)))

# ==============================================================================
rule("PART 2  BOTTOM-UP WAI, Whittaker & Woodwell (1967) relations")
# ==============================================================================
CONIC_CORR  <- c(1.20, 1.35)   # true stem surface / conic estimate, W&W Table 1:6
BRANCH_STEM <- c(2.7, 4.0)     # branch:stem ratio; 2.7-3.0 Smokies climax, 4.0 Brookhaven

S_conic <- sum(pi*INV$dbh*INV$H/2)                 # W&W conic stem surface
cat(sprintf("
  Conic stem surface  sum(pi*D*H/2)        %8.0f m2   index %.3f
  x true/conic correction %.2f-%.2f        %8.0f - %.0f m2   index %.3f - %.3f
", S_conic, S_conic/A_PLOT, CONIC_CORR[1], CONIC_CORR[2],
   S_conic*CONIC_CORR[1], S_conic*CONIC_CORR[2],
   S_conic*CONIC_CORR[1]/A_PLOT, S_conic*CONIC_CORR[2]/A_PLOT))

cat("\n  Total woody (bark) area index, stem x (1 + branch:stem):\n\n")
cat(sprintf("    %-28s %10s %10s\n","","branch:stem 2.7","branch:stem 4.0"))
grid <- expand.grid(cc=CONIC_CORR, bs=BRANCH_STEM)
for (cc in CONIC_CORR) {
  v <- sapply(BRANCH_STEM, function(bs) S_conic*cc*(1+bs)/A_PLOT)
  cat(sprintf("    conic corr %.2f %-15s %10.2f %10.2f\n", cc, "", v[1], v[2]))
}
wai_lo <- S_conic*CONIC_CORR[1]*(1+BRANCH_STEM[1])/A_PLOT
wai_hi <- S_conic*CONIC_CORR[2]*(1+BRANCH_STEM[2])/A_PLOT

# EXPORT, so rev_scaling_full_grid.R reads these instead of carrying its own copy.
# The grid used to hardcode `1.69 bottom-up, low` / `2.11 bottom-up, this stand` /
# `2.57 bottom-up, high` as literals, so correcting this producer would not have
# reached the headline at all. One producer per quantity.
WAI_OUT <- data.frame(
  label = c("bottom-up, low", "bottom-up, this stand", "bottom-up, high"),
  wai   = c(wai_lo, (wai_lo + wai_hi)/2, wai_hi),
  basis = sprintf("S_conic %.0f m2 over %.0f m2, conic corr %.2f-%.2f, branch:stem %.1f-%.1f",
                  S_conic, A_PLOT, CONIC_CORR[1], CONIC_CORR[2],
                  BRANCH_STEM[1], BRANCH_STEM[2]),
  n_stems = nrow(INV), stand_area_m2 = A_PLOT, stringsAsFactors = FALSE)
write.csv(WAI_OUT, out_path("wai_bottomup.csv"), row.names = FALSE)
cat(sprintf("\n  written: %s/wai_bottomup.csv  (%.2f / %.2f / %.2f from %d stems)\n",
            outdir, wai_lo, (wai_lo+wai_hi)/2, wai_hi, nrow(INV)))

cat(sprintf("
  ---------------------------------------------------------------------------
  WAI ESTIMATES, THREE INDEPENDENT ROUTES
  ---------------------------------------------------------------------------
    bottom-up, this inventory + W&W relations      %.2f - %.2f
    W&W published, climax temperate deciduous      2.00 - 2.20
    W&W published, Brookhaven oak-pine (open)      ~1.50
    Gauci et al. 2024, temperate forest            3.07
  ---------------------------------------------------------------------------
  READ THIS BEFORE QUOTING THE COMPARISON. Both statements that used to stand here
  are now false, because the bottom-up trio was recomputed on the canonical inventory
  (8,006 stems over the 37,200 m2 censused stand, not 7,010 over the nominal 40,000)
  and rose about 32%%:

    * that our bottom-up estimate brackets the W&W published temperate-deciduous
      value (2.00-2.20) -- it no longer brackets it, the range now sits ENTIRELY
      ABOVE it;
    * that the Gauci value sits at or above the top of that range, so the budget uses
      the HIGH end -- 3.07 now sits INSIDE our own range, between the mid and the
      high. It is no longer an upper bound and the budget is not using the high end.

  Note also that agreement with W&W's published value was never independent evidence:
  our estimate is built FROM their taper correction and branch:stem ratios, so any
  resemblance corroborates the arithmetic, not the value (scaling_parameters.md 4).
  The stand is 65%% conifer by basal area while those ratios come from deciduous and
  pitch-pine stands, which is the standing limitation. See decision E in
  manuscript/OPEN_DECISIONS.md.
", wai_lo, wai_hi))

# ==============================================================================
rule("PART 3  HEIGHT-FLUX FORMS: what did they buy us?")
# ==============================================================================
pair <- tree_train_complete %>% filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected))
hq <- c(50,125,200)/100
fb <- sapply(hq, function(z) mean(pair$stem_flux_corrected[abs(pair$measurement_height_cm/100-z)<1e-6]))
cat(sprintf("\n  stand-mean profile: %.4f (0.5 m)  %.4f (1.25 m)  %.4f (2 m)\n", fb[1],fb[2],fb[3]))
cat("
  THREE POINTS. A 2-parameter form has 1 residual df; a 3-parameter form has 0 and
  fits exactly. No model-selection statistic (AIC, F, cross-validation) can separate
  them -- any apparent ranking is an artefact. The forms are therefore not competing
  hypotheses to be selected among; they are BOUNDING SCENARIOS.
")
bo <- read.csv("data/compiled/black_oak_experiment.csv", stringsAsFactors=FALSE)
bo$h_m <- suppressWarnings(as.numeric(bo$height_m))
bo <- bo[is.finite(bo$h_m)&is.finite(bo$CH4_best_flux_nmol_m2_s),]
oak_hi <- bo[bo$h_m>2,]; at2 <- bo$CH4_best_flux_nmol_m2_s[abs(bo$h_m-2)<.01][1]
b_log <- coef(lm(log(pmax(fb,1e-9))~hq))[2]; b_lin <- coef(lm(fb~hq))
b_pow <- coef(lm(log(pmax(fb,1e-9))~log(hq)))[2]
FORMS <- list(
  constant  = function(z) rep(1,length(z)),
  exp_asym  = function(z) 0.5+0.5*exp(b_log*(z-2)),
  power     = function(z) exp(b_pow*(log(pmax(z,.05))-log(2))),
  exp_zero  = function(z) exp(b_log*(z-2)),
  linear0   = function(z) pmax(b_lin[1]+b_lin[2]*z,0)/max(b_lin[1]+b_lin[2]*2,1e-9))
cat(sprintf("  Validation against the felled oak (the ONLY data >2 m), n = %d:\n\n", nrow(oak_hi)))
cat(sprintf("    %-10s %10s %10s %10s\n","form","pred 6 m","pred 10 m","RMSE"))
for (nm in names(FORMS)) {
  p <- FORMS[[nm]](oak_hi$h_m)*at2
  cat(sprintf("    %-10s %10.4f %10.4f %10.4f\n", nm, FORMS[[nm]](6)*at2,
              FORMS[[nm]](10)*at2, sqrt(mean((oak_hi$CH4_best_flux_nmol_m2_s-p)^2))))
}
cat(sprintf("    %-10s %10.4f %10.4f %10s\n","OBSERVED",
            mean(bo$CH4_best_flux_nmol_m2_s[abs(bo$h_m-6)<.01]),
            mean(bo$CH4_best_flux_nmol_m2_s[abs(bo$h_m-10)<.01]),"--"))
cat("
  Every decaying form under-predicts the oak above 2 m; constant is closest. Combined
  with (a) the impossibility of discriminating forms on 3 points and (b) the fact that
  constant flux makes the answer independent of the unmeasured vertical area
  distribution, this is the defensible basis for a single main-text value -- not that
  constant is 'true', but that it is the only choice that does not smuggle in an
  unconstrained assumption.
")

# ==============================================================================
rule("PART 4  DOES THE RF NEED EXPLICIT INTERACTIONS?")
# ==============================================================================
d <- tree_train_complete; y <- d$stem_flux_corrected; k <- is.finite(y); d <- d[k,]; y <- y[k]
base <- data.frame(species=factor(d$species_clean), dbh_m=d$dbh_m,
  soil_moisture_at_tree=d$soil_moisture_at_tree, soil_temp_C_mean=d$soil_temp_C_mean,
  air_temp_C_mean=d$air_temp_C_mean,
  height_cm=ifelse(is.na(d$measurement_height_cm),125,d$measurement_height_cm))
variants <- list(
  `base (current, 7 predictors)` = base,
  `+ moisture x height`     = transform(base, mh = soil_moisture_at_tree*height_cm),
  `+ moisture x height, deeper (node 2)` = transform(base, mh = soil_moisture_at_tree*height_cm),
  `+ mois x ht, + soilT x ht` = transform(base, mh = soil_moisture_at_tree*height_cm,
                                                th = soil_temp_C_mean*height_cm))
cat(sprintf("\n  %-38s %10s %10s %10s\n","variant","OOB R2","mtry","min.node"))
for (i in seq_along(variants)) {
  nm <- names(variants)[i]
  mn <- if (grepl("deeper", nm)) 2 else 5
  r2 <- mean(sapply(1:5, function(s) ranger(x=variants[[i]], y=y, num.trees=1000,
              min.node.size=mn, mtry=max(2,floor(sqrt(ncol(variants[[i]])))),
              num.threads=1, seed=s)$r.squared))
  cat(sprintf("  %-38s %10.4f %10d %10d\n", nm, r2,
              max(2,floor(sqrt(ncol(variants[[i]])))), mn))
}
cat("
  Random forests represent interactions implicitly through successive splits, so an
  explicit product term is usually redundant. The comparison above shows whether that
  holds here. If the interaction term does not raise OOB R2, the honest conclusion is
  that these data do not support a moisture-dependent height response strong enough to
  improve prediction -- not that the mechanism is absent, but that 434 height rows from
  a single campaign cannot resolve it. That is a limitation to state, not to engineer
  around: adding a term that does not improve OOB R2 would be fitting noise.
")
cat("\nDone.\n")
