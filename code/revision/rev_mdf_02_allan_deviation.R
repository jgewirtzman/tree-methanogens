#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_02_allan_deviation.R      STEP 2 of the empirical-MDF rebuild
# ------------------------------------------------------------------------------
# Per-measurement instrument noise from the saved click.peak traces.
#
# METHOD (identical to ch4-data-filtering scripts/00_setup.R and 03_mdf_computation.R,
# so the two projects remain directly comparable)
#
#   Allan deviation at tau = 1 s:
#       sigma_allan = SD(diff(x)) / sqrt(2)
#   computed on the FIRST DIFFERENCE of the dry mole fraction within the measurement
#   window. First-differencing removes the smooth flux trend, so what remains is
#   instrument white noise rather than the signal we are trying to measure.
#   Restricted to flag == 1, i.e. the manually retained points from click.peak.
#
#   Christiansen et al. (2015) MDF:
#       MDF = (sigma_allan * 3 * t_crit) / t * flux.term
#   at 90 / 95 / 99% confidence, with t_crit from the t distribution on (t - 2) df.
#
# WHY THE SAVED TRACES AND NOT THE RAW ZIPS
#   The manID files written by click.peak contain the extracted 1 Hz trace for every
#   measurement (~500 points each) together with the MANUALLY CORRECTED window
#   (start.time_corr, obs.length_corr). Those manual adjustments represent a large
#   amount of irreplaceable human work; recomputing windows from the raw traces would
#   discard them. So noise is estimated on exactly the data goFlux fitted.
#
# COVERAGE (see rev_mdf_01_reconcile_windows.R)
#   height_2021     441/461  data/processed/flux/lgr_manual_identification_results_final.csv
#   semirigid_soil  288/296  data/processed/flux/lgr_manual_identification_results_soil.csv
#   semirigid_tree  369/396  deprecated/figures/flux_code_misc/lgr_manual_identification_results.csv
#   untagged_rescue  22/49   outputs/revision/untagged_manID.rds
#   static_2023       0/406  clicks lost to a filename collision between the 2021 and
#                            2023 processing scripts; handled separately.
#
# Output: outputs/revision/mdf_allan_per_measurement.csv
#         outputs/revision/fig_allan_deviation.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# --- method, verbatim from the companion project ------------------------------
allan_sd <- function(x) {
  x <- x[!is.na(x)]
  d <- diff(x)
  if (length(d) < 2) return(NA_real_)
  sd(d) / sqrt(2)
}

SRC <- list(
  height_2021    = list(path="data/processed/flux/lgr_manual_identification_results_final.csv", type="csv"),
  semirigid_soil = list(path="data/processed/flux/lgr_manual_identification_results_soil.csv",  type="csv"),
  semirigid_tree = list(path="deprecated/figures/flux_code_misc/lgr_manual_identification_results.csv", type="csv"),
  untagged       = list(path="outputs/revision/untagged_manID.rds", type="rds"))

rule("PER-MEASUREMENT ALLAN DEVIATION")
res <- list()
for (nm in names(SRC)) {
  s <- SRC[[nm]]
  if (!file.exists(s$path)) { cat(sprintf("  MISSING %s\n", s$path)); next }
  d <- if (s$type=="rds") readRDS(s$path) else read.csv(s$path, stringsAsFactors=FALSE)
  ch4 <- intersect(c("CH4dry_ppb","CH4d_ppb","CH4_ppb"), names(d))[1]
  co2 <- intersect(c("CO2dry_ppm","CO2d_ppm","CO2_ppm"), names(d))[1]
  if (is.na(ch4)) { cat(sprintf("  %-16s no CH4 column (%s)\n", nm, paste(head(names(d),6),collapse=","))); next }
  if ("flag" %in% names(d)) d <- d[!is.na(d$flag) & d$flag == 1, ]
  a <- d %>% group_by(UniqueID) %>%
    summarise(allan_sd_CH4 = allan_sd(.data[[ch4]]),
              allan_sd_CO2 = if (!is.na(co2)) allan_sd(.data[[co2]]) else NA_real_,
              n_pts = n(), .groups="drop") %>%
    mutate(campaign = nm)
  res[[nm]] <- a
  cat(sprintf("
  %-16s measurements %4d   points/meas median %4.0f
                   Allan CH4 (ppb)  median %6.3f   IQR %6.3f - %6.3f   range %6.3f - %7.3f
", nm, nrow(a), median(a$n_pts, na.rm=TRUE),
   median(a$allan_sd_CH4,na.rm=TRUE), quantile(a$allan_sd_CH4,.25,na.rm=TRUE),
   quantile(a$allan_sd_CH4,.75,na.rm=TRUE), min(a$allan_sd_CH4,na.rm=TRUE),
   max(a$allan_sd_CH4,na.rm=TRUE)))
}
A <- bind_rows(res)
# Save the FULL Allan table before any join. The joined table below keeps only the
# campaigns with a matching goFlux output, and an earlier version wrote only that,
# silently discarding the semirigid_tree and untagged Allan values.
write.csv(A, file.path(outdir,"mdf_allan_all_campaigns.csv"), row.names=FALSE)
cat(sprintf("\n  Full Allan table (all campaigns): %d measurements -> mdf_allan_all_campaigns.csv\n",
            nrow(A)))

rule("EMPIRICAL vs MANUFACTURER PRECISION")
SPEC <- 0.9
cat(sprintf("
  Manufacturer specification (ABB/LGR GLA131, 1-sigma at 1 s): %.2f ppb CH4
  Empirical, this dataset:
     median   %6.3f ppb   = %5.2fx the specification
     IQR      %6.3f - %6.3f ppb
     90th pct %6.3f ppb   = %5.2fx
     n        %d measurements
", SPEC, median(A$allan_sd_CH4,na.rm=TRUE), median(A$allan_sd_CH4,na.rm=TRUE)/SPEC,
   quantile(A$allan_sd_CH4,.25,na.rm=TRUE), quantile(A$allan_sd_CH4,.75,na.rm=TRUE),
   quantile(A$allan_sd_CH4,.90,na.rm=TRUE), quantile(A$allan_sd_CH4,.90,na.rm=TRUE)/SPEC,
   sum(is.finite(A$allan_sd_CH4))))
cat("\n  by campaign (median Allan SD, ppb, and multiple of spec):\n\n")
bc <- A %>% group_by(campaign) %>%
  summarise(n=n(), med=median(allan_sd_CH4,na.rm=TRUE),
            mult=median(allan_sd_CH4,na.rm=TRUE)/SPEC, .groups="drop") %>% arrange(med)
for (i in seq_len(nrow(bc))) with(bc[i,],
  cat(sprintf("    %-16s n=%4d   %6.3f ppb   %5.2fx spec\n", campaign, n, med, mult)))
cat("
  NOTE the analyser degrades over time, so precision must be estimated from THESE
  measurements in THIS period -- a value borrowed from another project or another
  year would not describe this instrument.
")

# --- rebuild MDF with empirical precision --------------------------------------
rule("CHRISTIANSEN MDF WITH EMPIRICAL PRECISION")
GF <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    transmute(UniqueID, best.flux, MDF_manuf=MDF, flux.term, nb.obs, LM.r2, quality.check),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    transmute(UniqueID, best.flux, MDF_manuf=MDF, flux.term, nb.obs, LM.r2, quality.check))
M <- A %>% inner_join(GF, by="UniqueID") %>%
  mutate(t_sec = nb.obs,
         df_m  = pmax(t_sec-2, 1),
         MDF_chr99 = (allan_sd_CH4*3*qt(.995, df_m))/t_sec*flux.term,
         MDF_chr95 = (allan_sd_CH4*3*qt(.975, df_m))/t_sec*flux.term,
         MDF_chr90 = (allan_sd_CH4*3*qt(.95,  df_m))/t_sec*flux.term,
         ratio = MDF_chr95/MDF_manuf)
cat(sprintf("
  matched to goFlux output: %d measurements

  MDF (nmol m-2 s-1)        median      IQR
    manufacturer          %8.5f   %8.5f - %8.5f
    Christiansen 95%%      %8.5f   %8.5f - %8.5f
    ratio (emp/manuf)     %8.2f   %8.2f - %8.2f
",
 nrow(M), median(M$MDF_manuf,na.rm=TRUE), quantile(M$MDF_manuf,.25,na.rm=TRUE),
 quantile(M$MDF_manuf,.75,na.rm=TRUE), median(M$MDF_chr95,na.rm=TRUE),
 quantile(M$MDF_chr95,.25,na.rm=TRUE), quantile(M$MDF_chr95,.75,na.rm=TRUE),
 median(M$ratio,na.rm=TRUE), quantile(M$ratio,.25,na.rm=TRUE), quantile(M$ratio,.75,na.rm=TRUE)))

cat("\n  DETECTION STATUS, manufacturer vs empirical:\n\n")
cat(sprintf("    %-18s %10s %10s\n","","manufacturer","Christiansen95"))
for (lab in c("below detection","negative","significantly negative")) {
  v1 <- switch(lab, `below detection`=mean(abs(M$best.flux)<M$MDF_manuf,na.rm=TRUE),
                    `negative`=mean(M$best.flux<0,na.rm=TRUE),
                    `significantly negative`=mean(M$best.flux < -M$MDF_manuf,na.rm=TRUE))
  v2 <- switch(lab, `below detection`=mean(abs(M$best.flux)<M$MDF_chr95,na.rm=TRUE),
                    `negative`=mean(M$best.flux<0,na.rm=TRUE),
                    `significantly negative`=mean(M$best.flux < -M$MDF_chr95,na.rm=TRUE))
  cat(sprintf("    %-18s %9.1f%% %9.1f%%\n", lab, 100*v1, 100*v2))
}
cat("\n  by campaign (Christiansen 95%):\n\n")
cat(sprintf("    %-16s %5s %10s %10s %12s\n","campaign","n","%below MDF","%negative","%sig.negative"))
for (cp in unique(M$campaign)) { q <- M[M$campaign==cp,]
  cat(sprintf("    %-16s %5d %9.1f%% %9.1f%% %11.1f%%\n", cp, nrow(q),
      100*mean(abs(q$best.flux)<q$MDF_chr95,na.rm=TRUE),
      100*mean(q$best.flux<0,na.rm=TRUE),
      100*mean(q$best.flux < -q$MDF_chr95,na.rm=TRUE))) }

write.csv(M, file.path(outdir,"mdf_allan_per_measurement.csv"), row.names=FALSE)
cat(sprintf("\n  Written: outputs/revision/mdf_allan_per_measurement.csv (%d rows)\n", nrow(M)))

# --- figure --------------------------------------------------------------------
p1 <- ggplot(A, aes(allan_sd_CH4, campaign, fill=campaign)) +
  geom_violin(scale="width", linewidth=.2, alpha=.75) +
  geom_vline(xintercept=SPEC, linetype="22", colour="#B2182B") +
  annotate("text", x=SPEC, y=0.6, label="manufacturer spec 0.9 ppb", hjust=-0.05,
           size=2.6, colour="#B2182B") +
  scale_x_log10() +
  labs(title="a  per-measurement Allan deviation by campaign",
       x=expression(sigma[Allan]~"at"~tau==1~s~(ppb~CH[4])), y=NULL) +
  theme_bw(base_size=8) + theme(legend.position="none")
p2 <- ggplot(M, aes(abs(best.flux))) +
  geom_point(aes(y=MDF_manuf, colour="manufacturer"), size=.7, alpha=.6) +
  geom_point(aes(y=MDF_chr95, colour="Christiansen 95%"), size=.7, alpha=.6) +
  geom_abline(slope=1, intercept=0, linetype="22") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c(manufacturer="grey55", `Christiansen 95%`="#2166AC"), name=NULL) +
  labs(title="b  detection limit vs measured flux",
       subtitle="points above the 1:1 line are indistinguishable from zero",
       x=expression("|"*CH[4]*" flux| (nmol m"^-2*" s"^-1*")"), y="MDF") +
  theme_bw(base_size=8)
suppressMessages(ggsave(file.path(outdir,"fig_allan_deviation.png"),
  gridExtra::grid.arrange(p1,p2,ncol=1), width=8, height=7, dpi=200, bg="white"))
cat("Written: outputs/revision/fig_allan_deviation.png\n")
