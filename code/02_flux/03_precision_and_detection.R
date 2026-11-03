#!/usr/bin/env Rscript
# ==============================================================================
# 03_precision_and_detection.R -- minimum detectable flux and detection status of
#                                 every CH4 flux
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root. Runs after 02_qc_c0_screen.R and
# 02b_instrument_precision.R.
#
# PRECISION (02b_instrument_precision.R)
#   sigma is the analyzer precision of the closure's field day: the median over that day's
#   closures of the second-difference precision inside each fitted window
#   (goFlux::empirical.prec, method "hadamard"). Every goFlux fit used it as prec.
#
# DETECTION
#   MDF = 1.96 * sigma / t * flux.term, t = closure length in SECONDS: goFlux's MDF column,
#   since the fits ran with conf = 0.95. t is recovered per closure from that column, so it
#   is exact for every logging interval. nb.obs is the number of logged samples and is NOT
#   t: the Height+molecular campaign logged every 5 s, so nb.obs there is 5x smaller than t
#   (nb.obs is used only for the n >= 3 filter).
#   This is the conventional chamber form. It is CONSERVATIVE: it treats the detectable
#   change as a single-sample criterion and so discards the averaging benefit of fitting a
#   slope to t points. We retain it deliberately, and report the regression test
#   (p < 0.05) alongside it.
#
# Measurements in qc_excluded_measurements.csv are left out.
#
# Inputs:  data/processed/flux/qc_excluded_measurements.csv (02)
#          data/processed/flux/CH4_best_flux_lgr_results_{2021_multiheight,soil,2023_cross_species}.csv
#          data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv
# Output:  data/processed/flux/flux_FINAL.csv  every CH4 flux with sigma, t_sec, MDF, detected,
#          class (read by the stage-B flux tables)
# ==============================================================================
source("code/lib/outputs.R")
source("code/lib/monthly_cols.R")
suppressMessages({library(dplyr)})
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# z for the reference MDF (two-sided 95%); every fit used its field day's sigma as prec with
# conf = 0.95 (code/lib/flux_windows.R), so goFlux's MDF column is this MDF already
Z <- qnorm(0.975)

G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_2021_multiheight.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Height+molecular", type="stem"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Monthly survey", type="soil"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_2023_cross_species.csv",stringsAsFactors=FALSE) %>%
    mutate(camp="Cross-species", type="stem")) %>%
  filter(is.finite(best.flux),is.finite(flux.term),is.finite(nb.obs),nb.obs>2)
st <- monthly_plain(read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE))
gc <- function(b){n<-paste0("CH4_",b); stopifnot(n %in% names(st)); suppressWarnings(as.numeric(st[[n]]))}
G <- bind_rows(G, data.frame(UniqueID=st$UniqueID, camp="Monthly survey", type="stem",
  best.flux=gc("best.flux"), flux.term=gc("flux.term"), nb.obs=gc("nb.obs"),
  MDF=gc("MDF"), prec=gc("prec"),
  LM.p.val=gc("LM.p.val"), LM.r2=gc("LM.r2"), stringsAsFactors=FALSE)) %>%
  filter(is.finite(best.flux),is.finite(flux.term),is.finite(nb.obs),nb.obs>2) %>%
  distinct(UniqueID,.keep_all=TRUE)
excl <- read.csv("data/processed/flux/qc_excluded_measurements.csv",stringsAsFactors=FALSE)$UniqueID
G <- G[!(G$UniqueID %in% excl),]
G$sigma <- G$prec                                 # the closure's field-day sigma (02b)
# closure length in seconds, from goFlux's MDF = 1.96 prec/t*flux.term (t in seconds).
# nb.obs is a sample count: equal to t only at 1 Hz, 5x too small at 5 s logging.
G$t_sec <- Z*G$prec*G$flux.term/G$MDF
stopifnot(all(is.finite(G$sigma)), all(G$sigma > 0), all(is.finite(G$t_sec)), all(G$t_sec > 0))
G$detected <- abs(G$best.flux) > G$MDF
G$class <- ifelse(G$best.flux > G$MDF, "emission",
           ifelse(G$best.flux < -G$MDF, "uptake", "below detection"))
S <- G[G$type=="stem",]; SO <- G[G$type=="soil",]

rule("PRECISION AND DETECTION LIMITS")
cat("\n"); print(G %>% group_by(camp) %>% summarise(n=n(), sigma_ppb=round(median(sigma),3),
  x_spec=round(median(sigma)/0.9,2), median_MDF=round(median(MDF),4),
  closure_s=sprintf("%.0f-%.0f", min(t_sec), max(t_sec)), .groups="drop") %>%
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
cat("\n"); print(S %>% group_by(camp) %>% summarise(sigma=round(median(sigma),3), n=n(),
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
write.csv(G, "data/processed/flux/flux_FINAL.csv", row.names=FALSE)
cat("\n  Written: data/processed/flux/flux_FINAL.csv\n")
