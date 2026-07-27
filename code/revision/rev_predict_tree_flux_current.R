#!/usr/bin/env Rscript
# ==============================================================================
# rev_predict_tree_flux_current.R
# ------------------------------------------------------------------------------
# Per-stem CH4 flux predictions for the map panel of Figure 9.
#
# INVENTORY: the fg_final tree list (7,965 stems), which is what the map panels have
#   always used -- it carries the geodetic transform and local-outlier removal from
#   07_maps/01_forestgeo_alignment.R. The tree list and coordinates are taken from the
#   committed predictions file; only the FLUX values are recomputed.
#
#   An earlier version of this script substituted PLACEHOLDER_INVENTORY (7,010 stems,
#   the 200 m census core). That silently dropped 963 stems and left visible gaps in
#   the map. The two lists are not interchangeable: the core is the right basis for the
#   BUDGET, the full list is the right basis for the MAP. Both are produced here, with
#   `in_core` marking membership of the census core.
#
# MODEL: the locked 7-predictor TreeRF (species, dbh_m, dbh_within_z, soil moisture,
#   soil temperature, air temperature, height). The previous flux values came from the
#   superseded one-hot / SI_tree / chamber-flag feature matrix and spanned only
#   0.052-0.084 nmol m-2 s-1.
#
# HEIGHT: the RF is a step function in height with three levels across 50-200 cm, so
#   the band integral converges quickly; HEIGHT_BINS below is sufficient and the
#   residual discretisation error is under ~1%.
#
# Output: outputs/tables/tree_flux_predictions.csv
# ==============================================================================
suppressMessages({library(ranger); library(dplyr)})
set.seed(42)
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
CORE <- rf_workflow_data$PLACEHOLDER_INVENTORY     # 200 m census core, for the flag
DR   <- rf_workflow_data$PLACEHOLDER_DRIVERS
d <- tree_train_complete; trained <- sort(unique(as.character(d$species_clean)))

STEM_BAND_M <- 2.00; KALMIA_BAND_M <- 0.75
HEIGHT_BINS <- seq(50, 200, by = 12.5)

# --- fg_final tree list (coordinates + DBH), fluxes to be recomputed ----------
FG <- read.csv("outputs/tables/tree_flux_predictions.csv", stringsAsFactors=FALSE)
FG$code <- trimws(FG$species)
lk <- CORE %>% distinct(species_code, species) %>% filter(!is.na(species_code))
FG$species_name <- lk$species[match(FG$code, lk$species_code)]
FG$sp <- ifelse(!is.na(FG$species_name) & FG$species_name %in% trained,
                FG$species_name, "SPECIES_OTHER")
FG$dbh <- FG$dbh_m
FG <- FG[is.finite(FG$dbh) & FG$dbh > 0, ]
FG <- FG %>% group_by(sp) %>%
  mutate(dbh_within_z = if(n()>1 && sd(dbh,na.rm=TRUE)>0) as.numeric(scale(dbh)) else 0) %>%
  ungroup() %>% as.data.frame()
FG$band_m  <- ifelse(FG$species_name %in% "Kalmia latifolia", KALMIA_BAND_M, STEM_BAND_M)
FG$A_stem  <- pi * FG$dbh * FG$band_m
FG$in_core <- FG$tree_id %in% CORE$tree_id
cat(sprintf("fg_final stems %d | mapped to a trained species %d | SPECIES_OTHER %d | in 200 m core %d\n",
            nrow(FG), sum(FG$sp != "SPECIES_OTHER"), sum(FG$sp == "SPECIES_OTHER"), sum(FG$in_core)))

mm <- d %>% group_by(month) %>% summarise(m=mean(soil_moisture_at_tree,na.rm=TRUE),.groups="drop")
sm <- soil_train_complete %>% group_by(month) %>% summarise(ms=mean(soil_moisture_at_site,na.rm=TRUE),.groups="drop")
DR <- DR %>% left_join(mm,"month") %>% left_join(sm,"month") %>% mutate(m=ifelse(is.finite(m),m,ms))
fl <- function(v,mo) approx(mo[is.finite(v)],v[is.finite(v)],xout=mo,rule=2)$y
DR$m <- fl(DR$m,DR$month); DR$soil_temp_C_mean <- fl(DR$soil_temp_C_mean,DR$month)

cat(sprintf("predicting %d stems x %d height bins x 12 months ...\n", nrow(FG), length(HEIGHT_BINS)))
P <- array(NA_real_, c(nrow(FG), length(HEIGHT_BINS), 12))
for (mo in 1:12) for (j in seq_along(HEIGHT_BINS))
  P[,j,mo] <- predict(TreeRF, data.frame(species=factor(FG$sp, levels=trained),
    dbh_m=FG$dbh, dbh_within_z=FG$dbh_within_z, soil_moisture_at_tree=DR$m[mo],
    soil_temp_C_mean=DR$soil_temp_C_mean[mo], air_temp_C_mean=DR$air_temp_C_mean[mo],
    height_cm=HEIGHT_BINS[j]), num.threads=1)$predictions

w_meas <- (STEM_BAND_M-0.5)/STEM_BAND_M; w_below <- 0.5/STEM_BAND_M
band <- w_meas*apply(P,c(1,3),mean) + w_below*P[,1,]
is_k <- FG$species_name %in% "Kalmia latifolia"
if (any(is_k)) band[is_k,] <- P[is_k,1,]

OUT <- data.frame(tree_id=FG$tree_id, species=FG$code, species_name=FG$species_name,
  dbh_m=FG$dbh, x=FG$x, y=FG$y, in_core=FG$in_core, band_m=FG$band_m, A_stem_m2=FG$A_stem,
  flux_band_nmol_m2_s=rowMeans(band), flux_2m_nmol_m2_s=rowMeans(P[,length(HEIGHT_BINS),]))
OUT$flux_nmol_m2_s <- OUT$flux_band_nmol_m2_s
write.csv(OUT, "outputs/tables/tree_flux_predictions.csv", row.names=FALSE)

A_PLOT <- 200*200; CONV <- 86400*365.25*16e-6
cat(sprintf("
  stems written        %d  (map uses all; budget uses in_core = %d)
  band flux  median %.4f  range %.4f - %.4f nmol m-2 s-1
  band budget, core-only stems : %.3f mg CH4 m-2 yr-1
", nrow(OUT), sum(OUT$in_core), median(OUT$flux_band_nmol_m2_s),
   min(OUT$flux_band_nmol_m2_s), max(OUT$flux_band_nmol_m2_s),
   sum(OUT$flux_band_nmol_m2_s[OUT$in_core]*OUT$A_stem_m2[OUT$in_core])/A_PLOT*CONV))
