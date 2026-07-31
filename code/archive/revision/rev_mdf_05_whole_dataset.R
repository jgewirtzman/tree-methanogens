#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_05_whole_dataset.R      STEP 5: detection status across EVERY campaign
# ------------------------------------------------------------------------------
# The earlier threshold comparison covered only the 727 measurements that had both a
# per-measurement Allan deviation AND a matched goFlux output. This assembles all
# five campaigns so the reported n and % describe the whole dataset.
#
# PRECISION SOURCE PER CAMPAIGN
#   height_2021     per-measurement Allan deviation   (clicks recovered)
#   semirigid_soil  per-measurement Allan deviation   (clicks recovered)
#   semirigid_tree  per-measurement Allan deviation   (clicks recovered from deprecated/)
#   untagged        per-measurement Allan deviation   (clicks in untagged_manID.rds)
#   static_2023     ROLLING-WINDOW per-instrument     (clicks destroyed by a filename
#                   1.298 ppb, scanned from the        collision; see task #21)
#                   2023 raw recordings themselves
#
#   The rolling-window scan reads the analyser's quietest periods, so it sits slightly
#   BELOW per-measurement Allan. Measured offset where both are available:
#       static_2021        Allan 1.038  rolling 0.935   ratio 0.90
#       semirigid_2020-21  Allan 2.000  rolling 1.634   ratio 0.82
#   The 2023 rolling value is therefore divided by 0.86 (the mean of those ratios) to
#   put it on the same footing as the Allan-based campaigns. This is an adjustment, and
#   it is applied only to 2023; it is reported rather than buried.
#
# THRESHOLD CHOICE
#   The manuscript makes claims in BOTH directions -- that stems emit, and that stems
#   do not take up. A liberal threshold flatters the emission claim; a conservative one
#   flatters the null. Neither extreme is neutral, which is the argument for the
#   conventional middle. All variants are reported so the choice is visible.
#
# Output: outputs/revision/mdf_whole_dataset.csv / .txt
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

ROLLING_2023 <- 1.298      # ppb, from rev_mdf_03
ROLL_TO_ALLAN <- 0.86      # measured rolling/Allan ratio where both exist
# FULL Allan table across all campaigns (1143), not the goFlux-joined subset.
ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
A   <- read.csv("outputs/revision/mdf_allan_per_measurement.csv", stringsAsFactors=FALSE)

# --- height_2021 + semirigid_soil: already joined to goFlux --------------------
P <- A %>% transmute(UniqueID, campaign, flux=best.flux, sigma=allan_sd_CH4,
                     t=nb.obs, flux.term, r2=LM.r2, prec_source="Allan (per-measurement)")

# --- semirigid_tree + untagged: goFlux columns live in the compiled tree dataset
# The ".y" columns are the ones the Zenodo compiler uses (rev_compile_zenodo_semirigid.R).
st <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE)
gcol <- function(base) {
  for (suf in c(".y", ".x", "")) {
    n <- paste0("CH4_", base, suf)
    if (n %in% names(st)) { v <- suppressWarnings(as.numeric(st[[n]]))
      if (any(is.finite(v))) return(v) }
  }
  rep(NA_real_, nrow(st))
}
allan_st <- ALL$allan_sd_CH4[match(st$UniqueID, ALL$UniqueID)]
camp_st  <- ALL$campaign[match(st$UniqueID, ALL$UniqueID)]
ST <- data.frame(UniqueID=st$UniqueID,
                 campaign=ifelse(is.na(camp_st), "semirigid_tree", camp_st),
                 flux=gcol("best.flux"), sigma=allan_st, t=gcol("nb.obs"),
                 flux.term=gcol("flux.term"), r2=gcol("LM.r2"),
                 prec_source="Allan (per-measurement)", stringsAsFactors=FALSE)
cat(sprintf("\n  semirigid tree file: %d rows | Allan matched %d | flux.term present %d | nb.obs present %d\n",
    nrow(ST), sum(is.finite(ST$sigma)), sum(is.finite(ST$flux.term)), sum(is.finite(ST$t))))

# --- static_2023: no Allan; use adjusted rolling-window precision --------------
s23 <- read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv",
                stringsAsFactors=FALSE)
S23 <- data.frame(UniqueID=s23$UniqueID, campaign="static_2023", flux=s23$best.flux,
                  sigma=ROLLING_2023/ROLL_TO_ALLAN, t=s23$nb.obs,
                  flux.term=s23$flux.term, r2=s23$LM.r2,
                  prec_source=sprintf("rolling-window %.2f ppb / %.2f", ROLLING_2023, ROLL_TO_ALLAN),
                  stringsAsFactors=FALSE)

D <- bind_rows(P, ST, S23) %>%
  filter(is.finite(flux), is.finite(sigma), is.finite(flux.term), is.finite(t), t > 2) %>%
  distinct(UniqueID, .keep_all=TRUE)
D$is_soil <- D$campaign == "semirigid_soil"

rule("COVERAGE")
cat("\n")
print(D %>% group_by(campaign, prec_source) %>%
  summarise(n=n(), median_sigma_ppb=round(median(sigma),3), .groups="drop") %>% as.data.frame(),
  row.names=FALSE)
cat(sprintf("\n  TOTAL %d measurements  (stems %d, soil %d)\n",
            nrow(D), sum(!D$is_soil), sum(D$is_soil)))

# --- thresholds ----------------------------------------------------------------
mk <- function(k, crit) D$sigma * k * crit / D$t * D$flux.term
V <- list(`manufacturer`      = 0.9/D$t*D$flux.term,
          `Wassmann 90%`      = mk(1, qnorm(.95)),
          `Wassmann 95%`      = mk(1, qnorm(.975)),
          `Wassmann 99%`      = mk(1, qnorm(.995)),
          `Christiansen 90%`  = mk(3, qt(.95,  pmax(D$t-2,1))),
          `Christiansen 95%`  = mk(3, qt(.975, pmax(D$t-2,1))),
          `Christiansen 99%`  = mk(3, qt(.995, pmax(D$t-2,1))))

rule("WHOLE DATASET: DETECTION STATUS")
cat(sprintf("\n  STEM MEASUREMENTS (n = %d)\n\n", sum(!D$is_soil)))
cat(sprintf("  %-18s %9s %9s %9s %9s %11s\n","threshold","DETECTED","%detected","sig.POS","sig.NEG","below MDF"))
R <- bind_rows(lapply(names(V), function(nm) {
  m <- V[[nm]]; tr <- !D$is_soil
  det <- abs(D$flux) >= m; sp <- D$flux >= m; sn <- D$flux <= -m
  cat(sprintf("  %-18s %9d %8.1f%% %9d %9d %11d\n", nm, sum(det[tr]), 100*mean(det[tr]),
      sum(sp[tr]), sum(sn[tr]), sum(!det[tr])))
  data.frame(threshold=nm,
    stem_n=sum(tr), stem_detected=sum(det[tr]), stem_pct_detected=100*mean(det[tr]),
    stem_sig_pos=sum(sp[tr]), stem_sig_neg=sum(sn[tr]),
    soil_n=sum(!tr), soil_detected=sum(det[!tr]), soil_pct_detected=100*mean(det[!tr]),
    soil_sig_pos=sum(sp[!tr]), soil_sig_neg=sum(sn[!tr]))
}))
cat(sprintf("\n  SOIL MEASUREMENTS (n = %d)\n\n", sum(D$is_soil)))
cat(sprintf("  %-18s %9s %9s %9s %9s %11s\n","threshold","DETECTED","%detected","sig.POS","sig.NEG","below MDF"))
for (i in seq_len(nrow(R))) with(R[i,],
  cat(sprintf("  %-18s %9d %8.1f%% %9d %9d %11d\n", threshold, soil_detected, soil_pct_detected,
      soil_sig_pos, soil_sig_neg, soil_n-soil_detected)))

cat("\n  BY CAMPAIGN, at Wassmann 95%:\n\n")
m95 <- V[["Wassmann 95%"]]
cat(sprintf("  %-16s %6s %9s %9s %9s\n","campaign","n","%detected","sig.POS","sig.NEG"))
for (cp in unique(D$campaign)) { k <- D$campaign==cp; m <- m95[k]; f <- D$flux[k]
  cat(sprintf("  %-16s %6d %8.1f%% %9d %9d\n", cp, sum(k), 100*mean(abs(f)>=m),
      sum(f>=m), sum(f<=-m))) }

write.csv(R, file.path(outdir,"mdf_whole_dataset.csv"), row.names=FALSE)
cat(sprintf("\n  Written: outputs/revision/mdf_whole_dataset.csv\n"))

rule("THE SYMMETRY ARGUMENT")
cat("
  The manuscript claims BOTH that stems emit AND that stems do not take up. Those two
  claims are flattered by OPPOSITE thresholds:

    a LIBERAL threshold  -> more measurements count as detected emissions (helps the
                            positive claim), but makes 'no significant uptake' easier
                            to achieve trivially
    a CONSERVATIVE one   -> the reverse

  So neither end is neutral for this manuscript, and the conventional middle is
  defensible precisely because it is not tuned to either claim. Whatever is chosen, the
  full table above should be reported so the choice is visible and the reader can see
  that the conclusions do not depend on it.
")
