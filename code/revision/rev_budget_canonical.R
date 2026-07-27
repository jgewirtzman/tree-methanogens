#!/usr/bin/env Rscript
# ==============================================================================
# rev_budget_canonical.R
# ------------------------------------------------------------------------------
# THE canonical budget. Every number downstream is computed here and written out;
# nothing is hardcoded in a figure script.
#
# TREE, MEASURED BAND (0-2 m)
#   The RF carries measurement height as a predictor, so the band is integrated
#   rather than represented by a single value: for each stem and each month the model
#   is evaluated at a series of heights spanning the measured range (50-200 cm), and
#   the per-stem profile is averaged over the band. Below 50 cm, outside the training
#   range, the 50 cm value is held. This replaces the earlier single-value-per-stem
#   scaling and is 1.4x larger, because flux declines with height and a value taken at
#   2 m understates the 0-2 m mean.
#
# GEOMETRY  matches compute_geometry() in code/04_scaling/02_rf_models.R:
#   lateral stem area = pi * DBH * 2 m, except Kalmia latifolia (a shrub) at 0.75 m.
#   soil area = plot area minus total basal area.
#
# SOIL  taken from the pipeline's monthly predictions (MONTHLY_FLUXES.csv), which are
#   current with the locked SoilRF. The two chambers removed by the C0 screen were
#   already absent from the RF training set, so that screen does not change the soil
#   budget; verified in rev_qc_c0_screen.R and the task #17 audit.
#
# Output: outputs/revision/canonical_budget.csv   (one row per quantity)
#         outputs/revision/canonical_monthly.csv  (monthly series, both components)
#         outputs/tables/tree_flux_predictions.csv  (per-stem, regenerated)
# ==============================================================================
suppressMessages({library(ranger); library(dplyr)})
set.seed(42)
outdir <- "outputs/revision"
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- rf_workflow_data$PLACEHOLDER_INVENTORY; DR <- rf_workflow_data$PLACEHOLDER_DRIVERS
d <- tree_train_complete; trained <- sort(unique(as.character(d$species_clean)))

PLOT_SIDE_M   <- 200
KALMIA_BAND_M <- 0.75      # shrub; matches compute_geometry()
STEM_BAND_M   <- 2.00
# The RF is a STEP function in height: it was trained at 50/125/200 cm and can only
# change value at split points, returning just three distinct levels across the band.
# The band integral is therefore a Riemann sum over a step function, and a coarse grid
# lands wherever the steps happen to fall (13 bins overestimated by 1.4%). A 1 cm grid
# resolves every step exactly and the integral is converged.
HEIGHT_BINS   <- seq(50, 200, by = 1)      # cm, 1 cm resolution over the measured range
A_PLOT <- PLOT_SIDE_M^2
SEC_PER_DAY <- 86400; DAYS_PER_MONTH <- 30.4
NMOL_TO_MG  <- 16 * 1e-6

fx <- function(dd,s){dd<-ifelse(!is.na(dd)&dd>3,dd/100,dd)
 dd<-ifelse(grepl("Betula",s)&!is.na(dd)&dd*100>200,dd/10,dd)
 dd<-ifelse(s=="Pinus strobus"&!is.na(dd)&dd*100>230,dd/10,dd)
 ifelse(s=="Kalmia latifolia"&!is.na(dd)&dd*100>100,dd/100,
 ifelse(s=="Kalmia latifolia"&!is.na(dd)&dd*100>10,dd/10,dd))}
INV$dbh <- fx(INV$dbh_m, INV$species); INV <- INV[is.finite(INV$dbh)&INV$dbh>0,]
INV <- INV %>% group_by(species) %>%
  mutate(dbh_within_z = if(n()>1 && sd(dbh,na.rm=TRUE)>0) as.numeric(scale(dbh)) else 0) %>%
  ungroup() %>% as.data.frame()
INV$sp     <- ifelse(INV$species %in% trained, INV$species, "SPECIES_OTHER")
INV$band_m <- ifelse(INV$species == "Kalmia latifolia", KALMIA_BAND_M, STEM_BAND_M)
INV$A_stem <- pi * INV$dbh * INV$band_m
INV$BA     <- pi * (INV$dbh/2)^2
A_SOIL <- A_PLOT - sum(INV$BA)

cat(sprintf("GEOMETRY\n  stems %d | stem area 0-2 m %.0f m2 (index %.4f)\n  basal area %.1f m2 | soil area %.0f m2 (%.3f of plot)\n",
  nrow(INV), sum(INV$A_stem), sum(INV$A_stem)/A_PLOT, sum(INV$BA), A_SOIL, A_SOIL/A_PLOT))

mm <- d %>% group_by(month) %>% summarise(m=mean(soil_moisture_at_tree,na.rm=TRUE),.groups="drop")
sm <- soil_train_complete %>% group_by(month) %>% summarise(ms=mean(soil_moisture_at_site,na.rm=TRUE),.groups="drop")
DR <- DR %>% left_join(mm,"month") %>% left_join(sm,"month") %>% mutate(m=ifelse(is.finite(m),m,ms))
fl <- function(v,mo) approx(mo[is.finite(v)],v[is.finite(v)],xout=mo,rule=2)$y
DR$m <- fl(DR$m,DR$month); DR$soil_temp_C_mean <- fl(DR$soil_temp_C_mean,DR$month)

cat(sprintf("\nPredicting %d stems x %d height bins x 12 months ...\n", nrow(INV), length(HEIGHT_BINS)))
P <- array(NA_real_, c(nrow(INV), length(HEIGHT_BINS), 12))
for (mo in 1:12) for (j in seq_along(HEIGHT_BINS))
  P[,j,mo] <- predict(TreeRF, data.frame(species=factor(INV$sp, levels=trained),
    dbh_m=INV$dbh, dbh_within_z=INV$dbh_within_z, soil_moisture_at_tree=DR$m[mo],
    soil_temp_C_mean=DR$soil_temp_C_mean[mo], air_temp_C_mean=DR$air_temp_C_mean[mo],
    height_cm=HEIGHT_BINS[j]), num.threads=1)$predictions

# Band integral per stem per month. Lateral area of a cylinder accumulates uniformly
# with height (dA/dz = pi*DBH is constant), so integrating flux x area over the band is
# the mean flux over height times the band area:
#     int_0^H f(z) * pi*D dz  =  pi*D*H * mean_z f(z)  =  A_stem * mean_z f(z)
# The mean is taken over a 1 cm grid from 0.5-2 m, with 0-0.5 m held at the 50 cm value
# since that is outside the RF's training range. For Kalmia the band is 0.75 m, so only
# the 0-0.75 m portion is used.
w_meas  <- (STEM_BAND_M-0.5)/STEM_BAND_M
w_below <- 0.5/STEM_BAND_M
band_month <- w_meas*apply(P, c(1,3), mean) + w_below*P[,1,]
# Kalmia: the 0.75 m band lies entirely below the 50 cm measurement, apart from 25 cm,
# so its band value is the 50 cm prediction throughout.
is_k <- INV$species == "Kalmia latifolia"
if (any(is_k)) band_month[is_k, ] <- P[is_k, 1, ]
tree_monthly <- colSums(band_month * INV$A_stem) / A_PLOT     # nmol m-2(ground) s-1

mf <- read.csv("outputs/tables/MONTHLY_FLUXES.csv")
soil_monthly <- mf$Phi_soil_umol_m2_s[match(1:12, mf$month)]  # already per m2 ground

M <- data.frame(month=1:12,
  tree_nmol_m2_s = tree_monthly, soil_nmol_m2_s = soil_monthly,
  tree_mg_m2_d = tree_monthly*SEC_PER_DAY*NMOL_TO_MG,
  soil_mg_m2_d = soil_monthly*SEC_PER_DAY*NMOL_TO_MG)
M$plot_mg_m2_d <- M$tree_mg_m2_d + M$soil_mg_m2_d
write.csv(M, file.path(outdir,"canonical_monthly.csv"), row.names=FALSE)

tree_ann <- sum(M$tree_mg_m2_d)*DAYS_PER_MONTH
soil_ann <- sum(M$soil_mg_m2_d)*DAYS_PER_MONTH
# annual mean of each stem's 2 m value, for comparison with the band-integrated basis
f2m_annual  <- rowMeans(P[, length(HEIGHT_BINS), ])
tree_ann_2m <- sum(f2m_annual*INV$A_stem)/A_PLOT*SEC_PER_DAY*365.25*NMOL_TO_MG

B <- data.frame(
  quantity = c("plot_area_m2","stem_area_0_2m_m2","stem_area_index","basal_area_m2",
               "soil_area_m2","soil_area_fraction","n_stems",
               "tree_measured_mg_m2_yr","tree_at_2m_basis_mg_m2_yr","soil_mg_m2_yr",
               "net_measured_mg_m2_yr","tree_pct_of_soil","tree_r2_oob","soil_r2_oob"),
  value = c(A_PLOT, sum(INV$A_stem), sum(INV$A_stem)/A_PLOT, sum(INV$BA),
            A_SOIL, A_SOIL/A_PLOT, nrow(INV),
            tree_ann, tree_ann_2m, soil_ann,
            tree_ann+soil_ann, 100*abs(tree_ann)/abs(soil_ann),
            TreeRF$r.squared, SoilRF$r.squared))
write.csv(B, file.path(outdir,"canonical_budget.csv"), row.names=FALSE)

cat("\nCANONICAL BUDGET\n")
for(i in seq_len(nrow(B))) cat(sprintf("  %-28s %14.4f\n", B$quantity[i], B$value[i]))

OUT <- data.frame(tree_id=INV$tree_id, species=INV$species, dbh_m=INV$dbh,
  x=INV$x, y=INV$y, band_m=INV$band_m, A_stem_m2=INV$A_stem,
  flux_band_nmol_m2_s=rowMeans(band_month), flux_2m_nmol_m2_s=f2m_annual)
OUT$flux_nmol_m2_s <- OUT$flux_band_nmol_m2_s
write.csv(OUT, "outputs/tables/tree_flux_predictions.csv", row.names=FALSE)
cat(sprintf("\n  Written: canonical_budget.csv, canonical_monthly.csv, tree_flux_predictions.csv (%d stems)\n", nrow(OUT)))
