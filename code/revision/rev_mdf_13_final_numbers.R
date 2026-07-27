#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_13_final_numbers.R
# ------------------------------------------------------------------------------
# Every number reported in the detection/uptake SI, computed under the FINAL method:
#   precision   = residual SD about a local 30 s linear fit, median of the quietest
#                 5% of windows, computed separately for each instrument-period
#   MDF         = z * sigma / t * flux.term, z at 90% confidence
#   QC          = C0 screen only (2 soil chambers)
#   treatment   = retain all values, report detection status
# One method throughout. No per-measurement estimator, no correction factors.
# ==============================================================================
suppressMessages({library(dplyr)})
outdir <- "outputs/revision"
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

SIGMA <- c(`Height+molecular`=0.935, `Monthly survey`=1.634, `Cross-species`=1.298)
Z <- qnorm(0.95)

G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(camp="Height+molecular", type="stem"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    mutate(camp="Monthly survey", type="soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(camp="Cross-species", type="stem")) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs>2)
st <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE)
g <- function(b){for(s in c(".y",".x","")){n<-paste0("CH4_",b,s)
  if(n %in% names(st)){v<-suppressWarnings(as.numeric(st[[n]])); if(any(is.finite(v))) return(v)}}
  rep(NA_real_,nrow(st))}
G <- bind_rows(G, data.frame(UniqueID=st$UniqueID, camp="Monthly survey", type="stem",
  best.flux=g("best.flux"), flux.term=g("flux.term"), nb.obs=g("nb.obs"),
  LM.p.val=g("LM.p.val"), LM.r2=g("LM.r2"), stringsAsFactors=FALSE)) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs>2) %>%
  distinct(UniqueID,.keep_all=TRUE)
excl <- read.csv(file.path(outdir,"qc_excluded_measurements.csv"), stringsAsFactors=FALSE)$UniqueID
G <- G[!(G$UniqueID %in% excl),]
G$sigma <- SIGMA[G$camp]
G$MDF <- G$sigma*Z/G$nb.obs*G$flux.term
S <- G[G$type=="stem",]; SO <- G[G$type=="soil",]

rule("PRECISION AND DETECTION LIMITS")
cat("\n"); print(G %>% group_by(camp) %>% summarise(n=n(), sigma_ppb=sigma[1],
  x_spec=round(sigma[1]/0.9,2), median_MDF=round(median(MDF),4), .groups="drop") %>%
  as.data.frame(), row.names=FALSE)
cat(sprintf("\n  median detected stem emission / its own MDF: %.1f\n",
  median((S$best.flux/S$MDF)[S$best.flux > S$MDF])))

rule("DETECTION STATUS")
cat(sprintf("\n  stems n=%d | detected %d (%.1f%%) | emissions %d | uptakes %d\n",
  nrow(S), sum(abs(S$best.flux)>S$MDF), 100*mean(abs(S$best.flux)>S$MDF),
  sum(S$best.flux>S$MDF), sum(S$best.flux < -S$MDF)))
cat(sprintf("  soil  n=%d | detected %d (%.1f%%) | uptakes %d | emissions %d\n",
  nrow(SO), sum(abs(SO$best.flux)>SO$MDF), 100*mean(abs(SO$best.flux)>SO$MDF),
  sum(SO$best.flux < -SO$MDF), sum(SO$best.flux>SO$MDF)))
sp <- S$best.flux[S$best.flux>S$MDF]; sn <- S$best.flux[S$best.flux < -S$MDF]
son <- SO$best.flux[SO$best.flux < -SO$MDF]
cat(sprintf("\n  detected stem emission  median %.4f  IQR %.4f - %.4f\n",
  median(sp), quantile(sp,.25), quantile(sp,.75)))
cat(sprintf("  detected stem uptake    median %.4f  IQR %.4f - %.4f\n",
  median(sn), quantile(sn,.25), quantile(sn,.75)))
cat(sprintf("  detected soil uptake    median %.4f  IQR %.4f - %.4f  (%.0fx the stem uptake)\n",
  median(son), quantile(son,.25), quantile(son,.75), abs(median(son)/median(sn))))

rule("BY CAMPAIGN")
cat("\n"); print(S %>% group_by(camp) %>% summarise(sigma=sigma[1], n=n(),
  pct_detected=round(100*mean(abs(best.flux)>MDF),1),
  n_uptake=sum(best.flux < -MDF), pct_uptake=round(100*mean(best.flux < -MDF),1),
  deepest=round(min(best.flux),4),
  med_r2_uptake=round(median(LM.r2[best.flux < -MDF], na.rm=TRUE),2), .groups="drop") %>%
  arrange(sigma) %>% as.data.frame(), row.names=FALSE)

rule("FILTER SENSITIVITY")
cat("\n"); cat(sprintf("  %-24s %8s %9s %8s %9s %9s\n","criterion","kept","%kept","uptakes","%uptake","bias"))
CR <- list(`retain all`=rep(TRUE,nrow(S)),
           `p < 0.05`=S$LM.p.val<0.05,
           `above MDF (90%)`=abs(S$best.flux)>S$MDF,
           `above MDF (95%)`=abs(S$best.flux)>S$MDF*qnorm(.975)/Z,
           `r2 > 0.7`=S$LM.r2>0.7, `r2 > 0.9`=S$LM.r2>0.9)
for (nm in names(CR)) { m <- CR[[nm]]; m[is.na(m)] <- FALSE; f <- S$best.flux[m]
  cat(sprintf("  %-24s %8d %8.1f%% %8d %8.1f%% %8.0f%%\n", nm, sum(m), 100*mean(m),
      sum(f<0), 100*mean(f<0), 100*(mean(f)/mean(S$best.flux)-1))) }

rule("NEGATIVE BOUND")
cat(sprintf("\n  detected uptakes n=%d\n    median %.4f | mean %.4f | 5th pct %.4f | deepest %.4f\n",
  length(sn), median(sn), mean(sn), quantile(sn,.05), min(sn)))
hb <- S$best.flux[S$camp=="Height+molecular" & S$best.flux < -S$MDF]
cat(sprintf("    deepest in the lowest-detection-limit campaign: %s\n",
  if(length(hb)) sprintf("%.4f", min(hb)) else "none"))
cat(sprintf("    as %% of median soil uptake: %.1f%%\n", 100*abs(median(sn))/abs(median(son))))
write.csv(G, file.path(outdir,"flux_final_method.csv"), row.names=FALSE)
cat("\n  Written: outputs/revision/flux_final_method.csv\n")
