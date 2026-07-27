#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_FINAL_precision_and_detection.R
# ------------------------------------------------------------------------------
# THE METHOD. Supersedes rev_mdf_02 through rev_mdf_13, which document the route.
#
# PRECISION
#   sigma = MAD(dx) / sqrt(2), where dx is the difference between consecutive 1 Hz
#   readings, taken over the ENTIRE record for each field period.
#
#   MAD (median absolute deviation) is a spread measure built on medians rather than
#   squares, so it is unaffected until more than half the data are contaminated.
#   Why this works without selecting quiet periods:
#     - At 1 s a real flux moves the headspace ~0.016 ppb (typical) to 0.13 ppb
#       (our largest observed), against 1-2 ppb of instrument noise. Every
#       second-to-second difference is therefore instrument-dominated regardless of
#       what the analyser is pointed at.
#     - A constant rate of change shifts every difference equally, and MAD measures
#       spread about the median, so a flux cancels exactly. Verified: sigma is
#       identical for drifts of 0, 0.13, 0.5 and 2.0 ppb/s.
#     - Moves, breaths and chamber changes create large isolated differences. MAD
#       tolerates these to ~40% contamination; SD does not (SD/MAD is 6-10x here).
#
#   VALIDATION: computed over the whole record it gives 1.200 ppb for the
#   Height+molecular period, against 1.038 ppb computed independently from
#   first differences inside manually-clicked closure windows only -- two disjoint
#   data subsets agreeing to 15%. Measured inside vs outside closures on the same
#   recordings, the ratio is 0.91, so the estimate does not depend on which part of
#   the record is used.
#
#   Approaches considered and rejected, documented in rev_mdf_02/03/13:
#     - quietest-5% of 30 s rolling windows: 30 s sits past the Allan minimum (~4 s)
#       so it captures real concentration change, and selecting an extreme quantile
#       biases the estimate low by ~23% even on pure white noise.
#     - a closed-loop or standard-gas run would be the textbook method, and is
#       recommended for future campaigns; none was recorded here.
#
# DETECTION
#   MDF = z * sigma / t * flux.term,  z at 90% confidence.
#   This is the conventional chamber form. It is CONSERVATIVE: it treats the
#   detectable change as a single-sample criterion and so discards the averaging
#   benefit of fitting a slope to t points. Regression theory gives a detection
#   limit about 4.8x smaller, which we verified reproduces the fitted slope standard
#   errors to within 2%. We retain the conventional form deliberately, and report the
#   regression test (p < 0.05) alongside it.
#
# Output: outputs/revision/FINAL_detection.txt, flux_FINAL.csv
# ==============================================================================
suppressMessages({library(dplyr)})
outdir <- "outputs/revision"
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

SIGMA <- c(`Height+molecular`=1.200, `Cross-species`=1.725, `Monthly survey`=2.181)
Z <- qnorm(0.95)

G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Height+molecular", type="stem"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Monthly survey", type="soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Cross-species", type="stem")) %>%
  filter(is.finite(best.flux),is.finite(flux.term),is.finite(nb.obs),nb.obs>2)
st <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE)
gc <- function(b){for(s in c(".y",".x","")){n<-paste0("CH4_",b,s)
  if(n %in% names(st)){v<-suppressWarnings(as.numeric(st[[n]])); if(any(is.finite(v))) return(v)}}
  rep(NA_real_,nrow(st))}
G <- bind_rows(G, data.frame(UniqueID=st$UniqueID, camp="Monthly survey", type="stem",
  best.flux=gc("best.flux"), flux.term=gc("flux.term"), nb.obs=gc("nb.obs"),
  LM.p.val=gc("LM.p.val"), LM.r2=gc("LM.r2"), stringsAsFactors=FALSE)) %>%
  filter(is.finite(best.flux),is.finite(flux.term),is.finite(nb.obs),nb.obs>2) %>%
  distinct(UniqueID,.keep_all=TRUE)
excl <- read.csv(file.path(outdir,"qc_excluded_measurements.csv"),stringsAsFactors=FALSE)$UniqueID
G <- G[!(G$UniqueID %in% excl),]
G$sigma <- SIGMA[G$camp]
G$MDF   <- Z*G$sigma/G$nb.obs*G$flux.term
G$detected <- abs(G$best.flux) > G$MDF
G$class <- ifelse(G$best.flux > G$MDF, "emission",
           ifelse(G$best.flux < -G$MDF, "uptake", "below detection"))
S <- G[G$type=="stem",]; SO <- G[G$type=="soil",]

rule("PRECISION AND DETECTION LIMITS")
cat("\n"); print(G %>% group_by(camp) %>% summarise(n=n(), sigma_ppb=sigma[1],
  x_spec=round(sigma[1]/0.9,2), median_MDF=round(median(MDF),4),
  closure_s=sprintf("%d-%d", min(nb.obs), max(nb.obs)), .groups="drop") %>%
  as.data.frame(), row.names=FALSE)

rule("DETECTION STATUS")
cat(sprintf("\n  STEMS n=%d\n    detected  %d (%.1f%%)\n    emissions %d\n    uptakes   %d\n    below MDF %d (%.1f%%)\n",
  nrow(S), sum(S$detected), 100*mean(S$detected), sum(S$class=="emission"),
  sum(S$class=="uptake"), sum(!S$detected), 100*mean(!S$detected)))
cat(sprintf("\n  SOIL  n=%d\n    detected  %d (%.1f%%)\n    uptakes   %d\n    emissions %d\n",
  nrow(SO), sum(SO$detected), 100*mean(SO$detected),
  sum(SO$class=="uptake"), sum(SO$class=="emission")))
pv <- S$LM.p.val<0.05; pv[is.na(pv)]<-FALSE
cat(sprintf("\n  regression p<0.05 (reported alongside): %d stems (%.1f%%)\n", sum(pv), 100*mean(pv)))

rule("FLUX DISTRIBUTIONS")
qs <- c(.05,.25,.5,.75,.95)
sh <- function(v,lab){ if(!length(v)){cat(sprintf("  %-24s none\n",lab));return(invisible())}
  cat(sprintf("\n  %-24s n=%4d  mean %8.4f  median %8.4f\n",lab,length(v),mean(v),median(v)))
  cat("    "); print(round(quantile(v,qs),4)) }
cat("\n  STEMS"); sh(S$best.flux,"all")
sh(S$best.flux[S$class=="emission"],"detected emissions")
sh(S$best.flux[S$class=="uptake"],"detected uptakes")
cat("\n  SOIL"); sh(SO$best.flux,"all")
sh(SO$best.flux[SO$class=="uptake"],"detected uptakes")

rule("BY CAMPAIGN")
cat("\n"); print(S %>% group_by(camp) %>% summarise(sigma=sigma[1], n=n(),
  pct_detected=round(100*mean(detected),1), n_uptake=sum(class=="uptake"),
  pct_uptake=round(100*mean(class=="uptake"),1), deepest=round(min(best.flux),4),
  med_r2_uptake=round(median(LM.r2[class=="uptake"],na.rm=TRUE),2), .groups="drop") %>%
  arrange(sigma) %>% as.data.frame(), row.names=FALSE)

rule("NEGATIVE BOUND FOR THE SCALING SCENARIO")
u <- S$best.flux[S$class=="uptake"]; su <- SO$best.flux[SO$class=="uptake"]
cat(sprintf("\n  detected stem uptakes n=%d\n    median %.4f | mean %.4f | 5th pct %.4f | deepest %.4f\n",
  length(u), median(u), mean(u), quantile(u,.05), min(u)))
hb <- S$best.flux[S$camp=="Height+molecular" & S$class=="uptake"]
cat(sprintf("    deepest in the lowest-sigma campaign: %s\n",
  if(length(hb)) sprintf("%.4f", min(hb)) else "none detected"))
cat(sprintf("    median detected soil uptake %.4f  -> stem uptake is %.1f%% of it\n",
  median(su), 100*abs(median(u))/abs(median(su))))
cat(sprintf("\n  >>> BOUND = %.4f nmol m-2 s-1 (median detected stem uptake)\n", median(u)))
write.csv(G, file.path(outdir,"flux_FINAL.csv"), row.names=FALSE)
cat("\n  Written: outputs/revision/flux_FINAL.csv\n")
