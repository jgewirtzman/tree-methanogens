source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# rev_qc_c0_screen.R          TIER 1 HARD QC: chamber contamination screen
# ------------------------------------------------------------------------------
# REVISION override, following the rev_compile_zenodo_semirigid.R pattern: originals
# under code/ are left untouched; this script annotates the compiled products.
#
# THE RULE
#   Reject a measurement if its chamber starting concentration C0 exceeds
#   C0_MULT x the median C0 of its campaign. A chamber that opens far above ambient
#   was not sampling ambient air at closure -- it was contaminated or carried over
#   from a previous high-concentration sample. This is a PHYSICAL diagnostic of what
#   went wrong, not a statistical threshold on the flux value, and so is a legitimate
#   Tier 1 removal under the framework in the companion ch4-data-filtering project.
#   Statistical filtering of *plausible* flux values is NOT done here: it biases the
#   mean upward by 17-135% (that project, sec 5.3.4).
#
# WHY 1.5x
#   There is a 4.5-fold GAP in the C0 distribution between the contaminated chambers
#   and everything else, so the threshold is insensitive across it:
#       49049 ppb (24.11x)  <- rejected
#       11687 ppb ( 5.74x)  <- rejected
#     ----------------------------- 4.5x gap, any threshold here gives the same answer
#        2570 ppb ( 1.26x)
#        2304 ppb ( 1.13x)
#        2274 ppb ( 1.12x)
#   Any multiplier from ~1.3x to ~5.7x rejects exactly the same two measurements.
#   1.5x sits in the middle of that plateau. Documented as a choice nonetheless.
#   Consistency check: the three next-highest C0 chambers (1.12-1.26x, fluxes +25.7,
#   +11.6, +51.3) all lie in WS_1-1 / WS_1-2, the collars already removed by
#   OUT_OF_STAND_COLLARS -- mildly elevated C0 and independently excluded.
#
# WHAT IT CATCHES
#   Soil (n=288)            2 measurements, both from the field session of 2020-06-17:
#                             20200617_2-3-WS_19   C0 49049 ppb (24.1x median)  flux -176.11
#                             20200617_2-1-WS_18   C0 11687 ppb ( 5.7x median)  flux  +31.60
#   Height campaign (n=441) 0 measurements. C0 range 1809-2156 ppb; a clean campaign.
#
#   These are exactly the two soil values that lie outside the literature-plausible
#   range AND survive the OUT_OF_STAND_COLLARS exclusion, and each is far out of line
#   with its own collar (WS_2-3 median 0.074, WS_2-1 median -0.122). The screen was
#   derived from chamber physics, not tuned to remove them.
#
# COVERAGE GAP -- READ THIS
#   C0 is only available for the 288 soil and 441 height-campaign measurements. The
#   static-chamber campaign has no CH4 C0 in the processed outputs, so it CANNOT be
#   screened this way and its rows are marked qc_c0 = NA (not "pass"). Do not read a
#   missing value as a clean bill of health.
#
# OUTPUTS
#   outputs/revision/qc_excluded_measurements.csv   canonical exclusion list
#   data/compiled/semirigid_chamber_flux.csv        + qc_pass, qc_reason columns
#                                                   (rows RETAINED; filter on qc_pass)
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
C0_MULT <- 1.5

rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# --- campaigns that carry C0 --------------------------------------------------
SRC <- list(
  soil   = "data/processed/flux/CH4_best_flux_lgr_results_soil.csv",
  height = "data/processed/flux/CH4_best_flux_lgr_results.csv")

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
write.csv(EXCL, out_path("qc_excluded_measurements.csv"), row.names=FALSE)
cat(sprintf("\n  Canonical exclusion list: %d measurements -> outputs/revision/qc_excluded_measurements.csv\n",
            nrow(EXCL)))

# --- annotate the compiled (Zenodo) product -----------------------------------
rule("ANNOTATING data/compiled/semirigid_chamber_flux.csv")
f <- "data/compiled/semirigid_chamber_flux.csv"
Z <- read.csv(f, stringsAsFactors=FALSE)
Z$qc_pass   <- !(Z$unique_id %in% EXCL$UniqueID)
Z$qc_reason <- EXCL$reason[match(Z$unique_id, EXCL$UniqueID)]
# rows we could not screen: static campaign has no CH4 C0
screened <- unique(unlist(lapply(SRC, function(p) read.csv(p, stringsAsFactors=FALSE)$UniqueID)))
Z$qc_c0_screened <- Z$unique_id %in% screened
write.csv(Z, f, row.names=FALSE)
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
cat("
  Rows are RETAINED with a flag rather than deleted: this file is the public archive,
  and an excluded measurement that is visible and explained is more useful than one
  that silently vanished. Downstream analyses must filter on qc_pass.
")

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

# --- guard for downstream ------------------------------------------------------
cat(sprintf("
  Downstream revision scripts should assert the flag exists, e.g.
      stopifnot(\"qc_pass\" %%in%% names(d))
      d <- d[d$qc_pass, ]
  so that re-running the ORIGINAL compiler (which does not add qc_pass) fails loudly
  instead of silently reverting this screen.
"))
cat("\nDone.\n")
