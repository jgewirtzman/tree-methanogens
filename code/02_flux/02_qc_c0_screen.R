#!/usr/bin/env Rscript
# ==============================================================================
# 02_qc_c0_screen.R -- procedural screen: chamber not at ambient at closure
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root.
#
# THE RULE
#   Reject a measurement if its chamber starting concentration C0 exceeds
#   C0_MULT x the median C0 of its campaign. A chamber that opens far above ambient
#   was not sampling ambient air at closure -- it was contaminated or carried over
#   from a previous high-concentration sample. This is a physical diagnostic of what
#   went wrong, not a statistical threshold on the flux value (the Tier 1 removal of
#   the framework in the companion ch4-data-filtering project). Plausible flux values
#   are not filtered statistically: that biases the mean upward by 17-135% (that
#   project, sec 5.3.4).
#
# WHY 1.5x
#   There is a 4.5-fold gap in the C0 distribution between the contaminated chambers
#   and everything else, so the threshold is insensitive across it:
#       49049 ppb (24.11x)  <- rejected
#       11687 ppb ( 5.74x)  <- rejected
#     ----------------------------- 4.5x gap, any threshold here gives the same answer
#        2570 ppb ( 1.26x)
#        2304 ppb ( 1.13x)
#        2274 ppb ( 1.12x)
#   Any multiplier from ~1.3x to ~5.7x rejects exactly the same two measurements;
#   1.5x sits in the middle of that plateau. The three next-highest C0 chambers
#   (1.12-1.26x) all lie in WS_1-1 / WS_1-2, the collars removed by
#   OUT_OF_STAND_COLLARS -- mildly elevated C0 and independently excluded.
#
# WHAT IT CATCHES
#   Soil (n=288): 2 measurements, both from the field session of 2020-06-17
#   (20200617_2-3-WS_19, C0 24.1x median; 20200617_2-1-WS_18, 5.7x). Each is far out of
#   line with its own collar. Height campaign (n=441): none (C0 1809-2156 ppb).
#
# COVERAGE
#   C0 is available only for the soil and 2021 height-campaign fits. The monthly stems
#   and the 2023 campaign are not screened (qc_c0_screened = FALSE), which is not the
#   same as passing.
#
# Inputs:  data/processed/flux/CH4_best_flux_lgr_results_soil.csv
#          data/processed/flux/CH4_best_flux_lgr_results_2021_multiheight.csv
#          data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv,
#          semirigid_tree_final_complete_dataset_soil.csv (summary only)
# Outputs: data/processed/flux/qc_excluded_measurements.csv  canonical exclusion list (the
#          qc_pass flag of code/tables/10_flux_stem.R and 11_flux_soil.R; rows are kept)
#          outputs/data/qc_c0_screened_ids.csv        deployments that were screened
# ==============================================================================

source("code/lib/outputs.R")
source("code/lib/monthly_cols.R")
suppressMessages({library(dplyr)})
C0_MULT <- 1.5

rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# --- campaigns that carry C0 --------------------------------------------------
SRC <- list(
  soil   = "data/processed/flux/CH4_best_flux_lgr_results_soil.csv",
  height = "data/processed/flux/CH4_best_flux_lgr_results_2021_multiheight.csv")

rule("C0 SCREEN")
EX <- list()
for (nm in names(SRC)) {
  d <- read.csv(SRC[[nm]], stringsAsFactors=FALSE)
  med <- median(d$C0, na.rm=TRUE); thr <- C0_MULT*med
  bad <- d[is.finite(d$C0) & d$C0 > thr, ]
  cat(sprintf("
  %-8s n = %3d   median C0 = %6.0f ppb   threshold = %6.0f ppb (%.1fx)
           C0 range %.0f - %.0f ppb   REJECTED: %d
", nm, nrow(d), med, thr, C0_MULT, min(d$C0,na.rm=TRUE), max(d$C0,na.rm=TRUE), nrow(bad)))
  if (nrow(bad)) {
    for (i in seq_len(nrow(bad))) cat(sprintf("             %-22s C0 %7.0f (%4.1fx)  flux %9.3f\n",
        bad$UniqueID[i], bad$C0[i], bad$C0[i]/med, bad$best.flux[i]))
    EX[[nm]] <- data.frame(campaign=nm, UniqueID=bad$UniqueID, C0_ppb=bad$C0,
      C0_multiple_of_median=round(bad$C0/med,2), flux_nmol_m2_s=bad$best.flux,
      goflux_quality=bad$quality.check,
      reason=sprintf("C0 %.0f ppb = %.1fx campaign median; chamber not at ambient at closure",
                     bad$C0, bad$C0/med), stringsAsFactors=FALSE)
  }
  # sensitivity of the threshold
  cat("             threshold sensitivity: ")
  cat(paste(sapply(c(1.2,1.5,2,5,10), function(m)
      sprintf("%.1fx->%d", m, sum(d$C0 > m*med, na.rm=TRUE))), collapse="  "), "\n")
}
EXCL <- bind_rows(EX)
write.csv(EXCL, "data/processed/flux/qc_excluded_measurements.csv", row.names=FALSE)
cat(sprintf("\n  Canonical exclusion list: %d measurements -> data/processed/flux/qc_excluded_measurements.csv\n",
            nrow(EXCL)))

# --- record what was screened ---------------------------------------------------
rule("SCREENED DEPLOYMENTS")
screened <- unique(unlist(lapply(SRC, function(p) read.csv(p, stringsAsFactors=FALSE)$UniqueID)))
write.csv(data.frame(UniqueID = screened), out_path("qc_c0_screened_ids.csv"), row.names = FALSE)
# annotate in memory for the summary below (the processed monthly-survey files;
# data/compiled/ is built later)
.zt <- monthly_plain(read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv", stringsAsFactors=FALSE, check.names=FALSE))
.zs <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv", stringsAsFactors=FALSE, check.names=FALSE)
Z <- rbind(
  data.frame(unique_id = .zt$UniqueID, measurement_type = "tree_stem", landscape_position = .zt$Plot.Letter,
             CH4_best_flux_nmol_m2_s = .zt$CH4_best.flux, stringsAsFactors = FALSE),
  data.frame(unique_id = .zs$UniqueID, measurement_type = "soil", landscape_position = .zs[["Plot letter"]],
             CH4_best_flux_nmol_m2_s = .zs$CH4_best.flux, stringsAsFactors = FALSE))
Z$qc_pass   <- !(Z$unique_id %in% EXCL$UniqueID)
Z$qc_reason <- EXCL$reason[match(Z$unique_id, EXCL$UniqueID)]
Z$qc_c0_screened <- Z$unique_id %in% screened
cat(sprintf("
  rows            %d
  qc_pass = TRUE  %d
  qc_pass = FALSE %d   %s
  C0-screened     %d (%.0f%%)   -- the remainder have no C0 and are UNSCREENED,
                                   not verified clean
",
 nrow(Z), sum(Z$qc_pass), sum(!Z$qc_pass),
 paste(Z$unique_id[!Z$qc_pass], collapse=", "),
 sum(Z$qc_c0_screened), 100*mean(Z$qc_c0_screened)))
# --- effect on the soil flux distribution -------------------------------------
rule("EFFECT")
s <- Z[Z$measurement_type=="soil" & is.finite(Z$CH4_best_flux_nmol_m2_s),]
s$plot   <- sub("^[0-9]+_([0-9]-[0-9])-.*$", "\\1", s$unique_id)
s$collar <- paste0(s$landscape_position, "_", s$plot)
CONV <- 86400*365.25*16e-6
show <- function(v, lab) cat(sprintf("  %-46s n=%3d  mean %8.4f  median %8.4f  -> %8.1f mg m-2 yr-1\n",
  lab, length(v), mean(v), median(v), mean(v)*CONV))
cat("\n  ALL soil collars:\n")
show(s$CH4_best_flux_nmol_m2_s, "before")
show(s$CH4_best_flux_nmol_m2_s[s$qc_pass], "after C0 screen")
k <- !(s$collar %in% c("WS_1-1","WS_1-2"))
cat("\n  SCALING SUBSET (OUT_OF_STAND_COLLARS excluded):\n")
show(s$CH4_best_flux_nmol_m2_s[k], "before")
show(s$CH4_best_flux_nmol_m2_s[k & s$qc_pass], "after C0 screen")
cat("
  The MEAN moves substantially while the MEDIAN barely shifts -- the signature of a
  small number of artefacts rather than a heavy-tailed real distribution.

  NOTE: these are raw measurement means, NOT the stand budget. The budget comes from
  SoilRF predictions over the surface and must be REFIT, not extrapolated from these
  numbers. Direction is certain (the sink becomes shallower); magnitude is not.
")

cat("\nDone.\n")
