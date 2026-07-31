#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_07_criteria_comparison.R      STEP 7
# ------------------------------------------------------------------------------
# EVERY PLAUSIBLE DETECTABILITY CRITERION, SIDE BY SIDE.
#
# Four families of criterion are in common use, and they answer different questions:
#
#   REGRESSION   is the fitted slope distinguishable from zero?
#                The statistically correct question, using all n points. SE scales as
#                sigma/sqrt(n), so a long measurement gains real power.
#                Blind spot: a smooth artefact (drift, slow leak, baseline shift)
#                produces a highly significant slope that is not a flux.
#
#   MDF          could the instrument have resolved a change this size?
#                Formulated as sigma*z/t*flux.term, i.e. WITHOUT the sqrt(n) benefit,
#                so it is systematically stricter than the regression test -- by
#                roughly sqrt(n/12). It is a conservative screen, not a detection
#                limit in the statistical sense. Its virtue is catching the smooth
#                artefacts the p-value misses.
#
#   FIT QUALITY  does the model describe the data well? (r2)
#                Note r2 is NOT a detectability criterion: a tiny flux measured very
#                cleanly has high r2, and a large flux with one outlier has low r2.
#                Included because it is widely used, and to show what it does.
#
#   goFlux FLAGS the package's own diagnostics (SE, AICc, g-factor, intercept...).
#                Procedural rather than statistical; catches specific failure modes.
#
# For each criterion this reports how many measurements are retained, and of those,
# how many are positive vs negative -- so the effect on the SOURCE/SINK conclusion is
# visible, not just the sample size.
#
# Output: outputs/revision/mdf_criteria_comparison.csv / .txt
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="height_2021"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="semirigid_soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="static_2023")) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs > 2)

SIG <- ALL %>% mutate(campaign=ifelse(campaign=="untagged","semirigid_tree",campaign)) %>%
  group_by(campaign) %>% summarise(sigma=median(allan_sd_CH4,na.rm=TRUE), .groups="drop")
SIG <- bind_rows(SIG, data.frame(campaign="static_2023", sigma=1.298/0.86))
G <- G %>% left_join(SIG, by="campaign") %>% filter(is.finite(sigma))
G$is_soil <- G$campaign=="semirigid_soil"
b <- G$sigma/G$nb.obs*G$flux.term          # goFlux form with empirical sigma
G$r2 <- ifelse(G$model=="HM" & is.finite(G$HM.r2), G$HM.r2, G$LM.r2)
G$noflag <- !is.finite(match(TRUE, nchar(trimws(ifelse(is.na(G$quality.check),"",G$quality.check)))>0)) |
            nchar(trimws(ifelse(is.na(G$quality.check),"",G$quality.check)))==0

CRIT <- list(
  `-- none (retain all) --`          = rep(TRUE, nrow(G)),
  `REGRESSION p < 0.05`              = G$LM.p.val < 0.05,
  `REGRESSION p < 0.01`              = G$LM.p.val < 0.01,
  `REGRESSION p < 0.001`             = G$LM.p.val < 0.001,
  `MDF goFlux form, empirical`       = abs(G$best.flux) > b,
  `MDF Wassmann 90%`                 = abs(G$best.flux) > b*qnorm(.95),
  `MDF Wassmann 95%`                 = abs(G$best.flux) > b*qnorm(.975),
  `MDF Wassmann 99%`                 = abs(G$best.flux) > b*qnorm(.995),
  `MDF Christiansen 95% (x3)`        = abs(G$best.flux) > b*3*qt(.975, pmax(G$nb.obs-2,1)),
  `R2 > 0.5`                         = G$r2 > 0.5,
  `R2 > 0.7`                         = G$r2 > 0.7,
  `R2 > 0.9`                         = G$r2 > 0.9,
  `no goFlux quality flag`           = nchar(trimws(ifelse(is.na(G$quality.check),"",G$quality.check)))==0,
  `p<0.05 AND > MDF(W95)`            = (G$LM.p.val<0.05) & (abs(G$best.flux) > b*qnorm(.975)),
  `p<0.05 AND R2>0.7`                = (G$LM.p.val<0.05) & (G$r2>0.7),
  `p<0.05 AND >MDF(W95) AND no flag` = (G$LM.p.val<0.05) & (abs(G$best.flux) > b*qnorm(.975)) &
                                       (nchar(trimws(ifelse(is.na(G$quality.check),"",G$quality.check)))==0))

report <- function(mask, lab, sub) {
  m <- mask & sub; m[is.na(m)] <- FALSE
  f <- G$best.flux[m]
  data.frame(criterion=lab, n_total=sum(sub), n_retained=sum(m),
             pct_retained=100*sum(m)/sum(sub),
             n_pos=sum(f>0), n_neg=sum(f<0), pct_neg_of_retained=100*mean(f<0),
             mean_retained=mean(f), median_retained=median(f),
             mean_all=mean(G$best.flux[sub]),
             bias_pct=100*(mean(f)/mean(G$best.flux[sub])-1))
}

for (grp in c("STEMS","SOIL")) {
  sub <- if (grp=="STEMS") !G$is_soil else G$is_soil
  rule(sprintf("%s  (n = %d)", grp, sum(sub)))
  cat(sprintf("\n  %-34s %7s %7s %7s %7s %8s %10s %9s\n",
      "criterion","kept","%kept","n POS","n NEG","%neg","mean flux","bias"))
  R <- bind_rows(lapply(names(CRIT), function(nm) {
    r <- report(CRIT[[nm]], nm, sub)
    cat(sprintf("  %-34s %7d %6.1f%% %7d %7d %7.1f%% %10.4f %8.0f%%\n", nm,
        r$n_retained, r$pct_retained, r$n_pos, r$n_neg, r$pct_neg_of_retained,
        r$mean_retained, r$bias_pct))
    cbind(group=grp, r)
  }))
  assign(paste0("R_", grp), R)
}
write.csv(bind_rows(R_STEMS, R_SOIL), file.path(outdir,"mdf_criteria_comparison.csv"),
          row.names=FALSE)

rule("READING THIS TABLE")
cat("
  'bias' is the change in the MEAN flux caused by applying the criterion as a FILTER.
  It is shown to make the cost visible, not to endorse filtering: for budgets and
  scaling the correct treatment is to retain all values and report detection status
  separately. A criterion with a large positive bias is one that preferentially
  discards small and negative fluxes.

  Note what r2 does. It is not a detectability criterion and behaves unlike one: a
  small flux measured cleanly passes, a large flux with one bad point fails. Compare
  its retained-set bias with the MDF criteria.

  The p-value and the MDF criteria disagree by construction. The regression uses all
  n points, so its standard error carries a 1/sqrt(n) factor the MDF formula omits;
  the MDF is therefore stricter by roughly sqrt(n/12), which for these measurement
  lengths is a factor of about 5.
")
cat("\n  Written: outputs/revision/mdf_criteria_comparison.csv\n")
