#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_height_extrapolation.R
# ------------------------------------------------------------------------------
# Can we extrapolate the RF's learned height effect above 2 m, as one more scenario?
#
# PART 1  Demonstrates that a random forest CANNOT extrapolate. Beyond the largest
#         split threshold every tree in the forest returns the same leaf value, so
#         predicting at 5 m, 10 m or 100 m returns EXACTLY the 2 m prediction. Naively
#         "extrapolating the RF" is therefore bit-for-bit identical to the constant
#         scenario -- it is not an independent scenario at all.
#
# PART 2  The coherent version. The RF's value is that it learns a height response
#         CONDITIONED on species, soil moisture, temperature and DBH. We can extract
#         each inventory tree's own RF-predicted profile across the measured band, fit
#         a functional form TO THAT PER-TREE PROFILE, and project it. The shape is then
#         tree-specific and RF-derived, rather than a single stand-mean curve applied
#         to everybody. This IS a distinct scenario and is what "extrapolating the
#         learned height effect" can legitimately mean.
#
#         Caveat that does not go away: the per-tree shape is still estimated only from
#         0.5-2.0 m, so projecting it to 25 m is extrapolation however it is dressed.
#         It belongs alongside the other functional forms as a scenario, not above them.
#
# Output: outputs/revision/rf_height_extrapolation.txt
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(tidyr)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- rf_workflow_data$PLACEHOLDER_INVENTORY; DR <- rf_workflow_data$PLACEHOLDER_DRIVERS
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

d <- tree_train_complete
trained <- sort(unique(as.character(d$species_clean)))

# ==============================================================================
rule("PART 1  THE RF IS FLAT ABOVE ITS TRAINING RANGE (demonstration)")
# ==============================================================================
probe <- data.frame(species=factor("Acer rubrum", levels=trained),
  dbh_m=0.3, dbh_within_z=0, soil_moisture_at_tree=0.25,
  soil_temp_C_mean=17, air_temp_C_mean=22,
  height_cm=c(50,125,200,250,500,1000,2500,10000))
probe$pred <- predict(TreeRF, probe, num.threads=1)$predictions
cat(sprintf("\n  training heights are 50 / 125 / 200 cm only\n\n  %10s %12s\n","height_cm","prediction"))
for (i in seq_len(nrow(probe)))
  cat(sprintf("  %10.0f %12.6f%s\n", probe$height_cm[i], probe$pred[i],
              if (probe$height_cm[i] > 200) "   <- identical to 200 cm" else ""))
cat(sprintf("\n  All predictions above 200 cm are identical: %s\n",
            all(abs(probe$pred[probe$height_cm>200] - probe$pred[probe$height_cm==200]) < 1e-12)))
cat("
  So 'extrapolate the RF upward' and 'hold the 2 m value constant' are THE SAME
  SCENARIO. Reporting them as two independent scenarios would be double-counting.
")

# ==============================================================================
rule("PART 2  PER-TREE RF-LEARNED SHAPE, PROJECTED")
# ==============================================================================
fix_dbh <- function(dd,s){ dd<-ifelse(!is.na(dd)&dd>3,dd/100,dd)
  dd<-ifelse(grepl("Betula",s)&!is.na(dd)&dd*100>200,dd/10,dd)
  dd<-ifelse(s=="Pinus strobus"&!is.na(dd)&dd*100>230,dd/10,dd)
  ifelse(s=="Kalmia latifolia"&!is.na(dd)&dd*100>100,dd/100,
  ifelse(s=="Kalmia latifolia"&!is.na(dd)&dd*100>10,dd/10,dd)) }
INV$dbh <- fix_dbh(INV$dbh_m, INV$species); INV <- INV[is.finite(INV$dbh)&INV$dbh>0,]
INV <- INV %>% group_by(species) %>%
  mutate(dbh_within_z = if (n()>1 && sd(dbh,na.rm=TRUE)>0) as.numeric(scale(dbh)) else 0) %>%
  ungroup() %>% as.data.frame()
INV$sp <- ifelse(INV$species %in% trained, INV$species, "SPECIES_OTHER")

mm <- d %>% group_by(month) %>% summarise(m=mean(soil_moisture_at_tree,na.rm=TRUE),.groups="drop")
sm <- soil_train_complete %>% group_by(month) %>% summarise(ms=mean(soil_moisture_at_site,na.rm=TRUE),.groups="drop")
DR <- DR %>% left_join(mm,"month") %>% left_join(sm,"month") %>% mutate(m=ifelse(is.finite(m),m,ms))
fl <- function(v,mo) approx(mo[is.finite(v)],v[is.finite(v)],xout=mo,rule=2)$y
DR$m <- fl(DR$m,DR$month); DR$soil_temp_C_mean <- fl(DR$soil_temp_C_mean,DR$month)

HG <- seq(50,200,by=12.5)
P <- array(NA_real_, c(nrow(INV), length(HG), 12))
for (mo in 1:12) for (j in seq_along(HG)) {
  nd <- data.frame(species=factor(INV$sp, levels=trained), dbh_m=INV$dbh,
    dbh_within_z=INV$dbh_within_z, soil_moisture_at_tree=DR$m[mo],
    soil_temp_C_mean=DR$soil_temp_C_mean[mo], air_temp_C_mean=DR$air_temp_C_mean[mo],
    height_cm=HG[j])
  P[,j,mo] <- predict(TreeRF, nd, num.threads=1)$predictions
}
PROF <- apply(P, c(1,2), mean)            # per-tree annual profile over 0.5-2.0 m
z <- HG/100

# per-tree log-linear slope of the RF's own learned profile
slope <- as.numeric(apply(PROF, 1, function(f) {
  if (any(!is.finite(f)) || all(f <= 0) || sd(f) == 0) return(0)
  as.numeric(coef(lm(log(pmax(f,1e-9)) ~ z))[2]) }))
f2 <- PROF[, length(HG)]

cat(sprintf("\n  Per-tree RF-learned log-linear height slope (per m), across %d stems:\n", nrow(INV)))
print(round(quantile(slope, c(0,.1,.25,.5,.75,.9,1)), 4))
cat(sprintf("\n  stems with a DECLINING learned profile: %d (%.0f%%)\n",
            sum(slope < 0), 100*mean(slope < 0)))
cat(sprintf("  implied ratio f(2m)/f(0.5m), median: %.3f\n", median(exp(slope*1.5))))

cat("\n  Slope by species (RF-learned, area-weighted by pi*DBH):\n\n")
bs <- data.frame(sp=INV$sp, sl=slope, wt=as.numeric(pi*INV$dbh)) %>% group_by(sp) %>%
  summarise(n=n(), slope=sum(sl*wt)/sum(wt), .groups="drop") %>%
  mutate(ratio_2_05=exp(slope*1.5)) %>% arrange(slope)
cat(sprintf("    %-24s %7s %10s %12s\n","species","n","slope","f2m/f0.5m"))
for (i in seq_len(nrow(bs)))
  cat(sprintf("    %-24s %7d %10.4f %12.3f\n", bs$sp[i], bs$n[i], bs$slope[i], bs$ratio_2_05[i]))

# --- project each tree's own learned shape and integrate to the canopy --------
A_PLOT <- 200*200; CONV <- 86400*365.25*16e-6
GY <- c("Pinus strobus","Tsuga canadensis")
D_ANCH <- quantile(INV$dbh[INV$dbh>0.10], 0.95, na.rm=TRUE)
INV$H <- 1.37 + ifelse(INV$species %in% GY, (25-1.37)/D_ANCH^0.60, (25-1.37)/D_ANCH^0.53) *
         INV$dbh^ifelse(INV$species %in% GY, 0.60, 0.53)
INV$H <- ifelse(INV$species=="Kalmia latifolia", pmin(INV$H, 3), INV$H)   # shrub cap
S_conic <- pi*INV$dbh*INV$H/2
for (BS in c(2.7, 4.0)) for (CC in c(1.20, 1.35)) {
  A_w <- S_conic*CC*(1+BS)
  whole_const <- sum(f2*A_w)/A_PLOT*CONV
  whole_rf <- sum(sapply(seq_len(nrow(INV)), function(i) {
    zz <- seq(0.05, INV$H[i], length.out=120)
    f2[i]*mean(exp(slope[i]*(zz-2))) * A_w[i] }))/A_PLOT*CONV
  if (CC==1.20 && BS==2.7) cat(sprintf("\n  %-46s %10s %10s %8s\n",
        "WAI configuration","constant","RF shape","ratio"))
  cat(sprintf("  branch:stem %.1f, taper %.2f  (WAI %.2f)        %10.2f %10.2f %7.2fx\n",
      BS, CC, sum(A_w)/A_PLOT, whole_const, whole_rf, whole_rf/whole_const))
}
cat("
  The RF-learned per-tree shape is a legitimate additional scenario -- the shape comes
  from species, moisture, temperature and DBH rather than being imposed -- but it is
  still an out-of-range projection of a slope estimated over 1.5 m of stem, and the
  cross-validation on the 143 three-height trees showed log-linear extrapolation to be
  the LEAST reliable form tested upward (RMSE 6.90 vs 0.43 for constant). Report it as
  a low-end bounding scenario, not as the best estimate.
")
cat("\nDone.\n")
