#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_12_final_filter.R      STEP 12: the filtering decision, and what it leaves
# ------------------------------------------------------------------------------
# TIER 1 -- hard QC. Removal requires evidence of PROCEDURAL failure, never a
# threshold on the flux value itself (filtering plausible values biases the mean;
# reproduced on our own data at +233%).
#   (a) C0 > 1.5x campaign median          chamber grossly contaminated
#   (b) C0 - local ambient > 50 ppb        chamber not at ambient at closure.
#       Verified on the raw trace for 20230623_974_091900 (CH4 spike to 3796 ppb
#       minutes before; chamber opened 152 ppb high; smooth decay fitted as -0.391,
#       r2 0.964). Dose-response: negative prevalence 16% -> 73% and significant
#       negatives 3.7% -> 40% from lowest to highest excess bin; |negative flux|
#       scales with excess (rho 0.324, p 0.0025). Threshold at 50 ppb, where
#       prevalence roughly triples over baseline.
#       MEASURABLE ONLY WHERE RAW TRACES WERE MATCHED (2023).
#
# TIER 2 -- detection: flag, do not remove. Wassmann 95% with campaign-level
# empirical precision. Regression p<0.05 reported alongside as the statistically
# correct test; MDF is uniformly stricter and is the conservative companion.
#
# Output: outputs/revision/final_filter.txt, flux_after_tier1.csv
# ==============================================================================
suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")
EXCESS_LIM <- 50

ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="height_2021"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="semirigid_soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="static_2023")) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs>2)
st <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE)
g <- function(b){for(s in c(".y",".x","")){n<-paste0("CH4_",b,s)
  if(n %in% names(st)){v<-suppressWarnings(as.numeric(st[[n]])); if(any(is.finite(v))) return(v)}}
  rep(NA_real_,nrow(st))}
G <- bind_rows(G, data.frame(UniqueID=st$UniqueID, campaign="semirigid_tree",
  best.flux=g("best.flux"), flux.term=g("flux.term"), nb.obs=g("nb.obs"),
  LM.p.val=g("LM.p.val"), LM.r2=g("LM.r2"), C0=g("C0"), stringsAsFactors=FALSE)) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs>2) %>%
  distinct(UniqueID, .keep_all=TRUE)
SIG <- ALL %>% mutate(campaign=ifelse(campaign=="untagged","semirigid_tree",campaign)) %>%
  group_by(campaign) %>% summarise(sigma=median(allan_sd_CH4,na.rm=TRUE), .groups="drop")
SIG <- bind_rows(SIG, data.frame(campaign="static_2023", sigma=1.298/0.86))
G <- G %>% left_join(SIG,"campaign") %>% filter(is.finite(sigma))
G$MDF95 <- G$sigma*qnorm(.975)/G$nb.obs*G$flux.term
G$is_soil <- G$campaign=="semirigid_soil"

excl_c0 <- read.csv("outputs/revision/qc_excluded_measurements.csv", stringsAsFactors=FALSE)$UniqueID
L <- read.csv("outputs/revision/local_ambient_check.csv", stringsAsFactors=FALSE)
G$excess <- L$excess[match(G$UniqueID, L$UniqueID)]
G$t1_c0 <- G$UniqueID %in% excl_c0
G$t1_excess <- is.finite(G$excess) & G$excess > EXCESS_LIM
G$tier1_fail <- G$t1_c0 | G$t1_excess

rule("TIER 1: GENUINELY EXCLUDABLE")
cat(sprintf("\n  (a) C0 > 1.5x campaign median            %2d\n  (b) excess above local ambient > %d ppb   %2d\n      total excluded  %2d of %d (%.1f%%)\n\n  carryover-screenable: %d of %d (%.0f%%), all 2023; the rest UNSCREENED\n",
 sum(G$t1_c0), EXCESS_LIM, sum(G$t1_excess), sum(G$tier1_fail), nrow(G),
 100*mean(G$tier1_fail), sum(is.finite(G$excess)), nrow(G), 100*mean(is.finite(G$excess))))
e <- G[G$t1_excess,] %>% arrange(desc(excess))
cat("\n  excluded by excess:\n\n")
cat(sprintf("    %-28s %9s %8s %8s\n","UniqueID","flux","excess","r2"))
for (i in seq_len(nrow(e))) with(e[i,], cat(sprintf("    %-28s %9.4f %8.1f %8.3f\n",
  substr(UniqueID,1,28), best.flux, excess, ifelse(is.na(LM.r2),NA,LM.r2))))
cat(sprintf("\n    %d of %d negative (%.0f%%) vs %.0f%% among retained 2023 measurements\n",
  sum(e$best.flux<0), nrow(e), 100*mean(e$best.flux<0),
  100*mean(G$best.flux[G$campaign=="static_2023" & !G$tier1_fail] < 0)))

K <- G[!G$tier1_fail,]; S <- K[!K$is_soil,]; SO <- K[K$is_soil,]
rule("FLUX DISTRIBUTIONS AFTER TIER 1")
cat(sprintf("\n  stems %d (from %d) | soil %d (from %d)\n", nrow(S), sum(!G$is_soil),
            nrow(SO), sum(G$is_soil)))
qs <- c(.01,.05,.10,.25,.50,.75,.90,.95,.99)
show <- function(v,lab){ if(!length(v)){cat(sprintf("  %-24s none\n",lab));return(invisible())}
  cat(sprintf("\n  %-24s n=%4d  mean %8.4f  median %8.4f\n", lab, length(v), mean(v), median(v)))
  cat("    "); print(round(quantile(v,qs),4)) }
cat("\n  STEMS"); show(S$best.flux,"all")
show(S$best.flux[S$best.flux>0],"positive only")
show(S$best.flux[S$best.flux<0],"negative only")
show(S$best.flux[S$best.flux>S$MDF95],"significant positive")
show(S$best.flux[S$best.flux< -S$MDF95],"significant negative")
cat("\n  SOIL"); show(SO$best.flux,"all")
show(SO$best.flux[SO$best.flux< -SO$MDF95],"significant negative")

rule("DETECTION SUMMARY AFTER TIER 1")
cat(sprintf("\n  %-26s %8s %8s %8s %8s\n","","stems","%","soil","%"))
for (nm in c("regression p<0.05","above MDF (Wassmann 95%)","sig. positive","sig. negative")) {
  fs <- switch(nm, `regression p<0.05`=S$LM.p.val<0.05,
    `above MDF (Wassmann 95%)`=abs(S$best.flux)>S$MDF95,
    `sig. positive`=S$best.flux>S$MDF95, `sig. negative`=S$best.flux< -S$MDF95)
  fo <- switch(nm, `regression p<0.05`=SO$LM.p.val<0.05,
    `above MDF (Wassmann 95%)`=abs(SO$best.flux)>SO$MDF95,
    `sig. positive`=SO$best.flux>SO$MDF95, `sig. negative`=SO$best.flux< -SO$MDF95)
  fs[is.na(fs)]<-FALSE; fo[is.na(fo)]<-FALSE
  cat(sprintf("  %-26s %8d %7.1f%% %8d %7.1f%%\n", nm, sum(fs), 100*mean(fs), sum(fo), 100*mean(fo)))
}

rule("A NEGATIVE BOUND FOR THE LINEAR SCENARIO?")
sn <- S[S$best.flux < -S$MDF95,]
cl <- sn[!is.finite(sn$excess) | sn$excess <= 10,]
cat(sprintf("\n  significant negatives surviving Tier 1 : %d\n    at excess <= 10 ppb (screened)      : %d\n    in UNSCREENED campaigns             : %d\n\n  candidate bounds (nmol m-2 s-1):\n    most negative surviving Tier 1      %s\n    median significant negative         %s\n    most negative at excess <= 10 ppb   %s\n  scale: soil median %.3f | stem positive median %.4f\n",
 nrow(sn), sum(is.finite(sn$excess) & sn$excess<=10), sum(!is.finite(sn$excess)),
 if(nrow(sn)) sprintf("%.4f",min(sn$best.flux)) else "-",
 if(nrow(sn)) sprintf("%.4f",median(sn$best.flux)) else "-",
 if(nrow(cl)) sprintf("%.4f",min(cl$best.flux)) else "-",
 median(SO$best.flux), median(S$best.flux[S$best.flux>0])))
cat("\n  surviving significant negatives:\n\n")
cat(sprintf("    %-28s %9s %8s %7s %-14s\n","UniqueID","flux","excess","r2","campaign"))
sn2 <- sn[order(sn$best.flux),]
for (i in seq_len(nrow(sn2))) with(sn2[i,], cat(sprintf("    %-28s %9.4f %8s %7.3f %-14s\n",
  substr(UniqueID,1,28), best.flux, ifelse(is.finite(excess),sprintf("%.1f",excess),"unscr"),
  ifelse(is.na(LM.r2),NA,LM.r2), campaign)))
write.csv(K, file.path(outdir,"flux_after_tier1.csv"), row.names=FALSE)
cat("\n  Written: outputs/revision/flux_after_tier1.csv\n")
