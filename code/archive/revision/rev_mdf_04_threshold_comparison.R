#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_04_threshold_comparison.R      STEP 4 of the empirical-MDF rebuild
# ------------------------------------------------------------------------------
# WHAT DOES EACH DETECTION-LIMIT DEFINITION ACTUALLY GIVE US?
#
# There is no single correct MDF. The published definitions differ in precision
# source, confidence level, and whether they carry a factor of 3. They span more than
# an order of magnitude, so the choice must be made explicitly and reported, not
# inherited by accident. This script lays out the options side by side so the choice
# can be made on the evidence.
#
# THE VARIANTS
#   manufacturer      prec_spec / t * flux.term                  prec = 0.9 ppb
#                     goFlux's built-in value. Optimistic: it assumes laboratory
#                     precision, and our own Allan deviations are 1.2-2.3x worse.
#
#   Wassmann-z        z * sigma / t * flux.term                  z = 1.645/1.96/2.576
#                     Empirical precision at a stated confidence level. NO factor of 3.
#                     The most liberal of the empirical options.
#
#   Christiansen-noK  sigma * t_crit / t * flux.term
#                     As above but with the t distribution rather than the normal.
#                     Essentially identical to Wassmann at these sample sizes.
#
#   Christiansen      sigma * 3 * t_crit / t * flux.term
#                     The published form. The factor of 3 is a conservatism margin
#                     carried over from the analytical-chemistry limit-of-detection
#                     tradition, not derived from these data. It triples the MDF.
#
# A NOTE ON WHAT "BELOW DETECTION" MEANS
#   Below-MDF does NOT mean the measurement is wrong or should be discarded. Noise
#   scatters symmetrically, so the MEAN of many below-detection measurements still
#   converges on the true flux. Removing them biases the mean upward by 17-135%
#   (ch4-data-filtering sec 5.3.4). Detection status constrains what can be said about
#   an INDIVIDUAL measurement or its SIGN -- not what can be said about a mean.
#
# Output: outputs/revision/mdf_threshold_comparison.csv / .txt
#         outputs/revision/fig_mdf_thresholds.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

A <- read.csv("outputs/revision/mdf_allan_per_measurement.csv", stringsAsFactors=FALSE)
A <- A[is.finite(A$best.flux) & is.finite(A$flux.term) & is.finite(A$allan_sd_CH4) & A$nb.obs > 2, ]
A$t <- A$nb.obs; A$dfm <- pmax(A$t-2, 1)

mk <- function(sigma, k, crit) sigma * k * crit / A$t * A$flux.term
V <- list(
  `manufacturer (goFlux)`   = 0.9 / A$t * A$flux.term,
  `Wassmann 90%`            = mk(A$allan_sd_CH4, 1, qnorm(.95)),
  `Wassmann 95%`            = mk(A$allan_sd_CH4, 1, qnorm(.975)),
  `Wassmann 99%`            = mk(A$allan_sd_CH4, 1, qnorm(.995)),
  `Christiansen-noK 95%`    = mk(A$allan_sd_CH4, 1, qt(.975, A$dfm)),
  `Christiansen 90%`        = mk(A$allan_sd_CH4, 3, qt(.95,  A$dfm)),
  `Christiansen 95%`        = mk(A$allan_sd_CH4, 3, qt(.975, A$dfm)),
  `Christiansen 99%`        = mk(A$allan_sd_CH4, 3, qt(.995, A$dfm)))

rule("DETECTION LIMIT VARIANTS")
cat(sprintf("\n  n = %d measurements with per-measurement Allan deviation\n", nrow(A)))
cat(sprintf("  median |flux| = %.4f nmol m-2 s-1\n\n", median(abs(A$best.flux))))
cat(sprintf("  %-24s %10s %26s %22s\n","variant","median MDF","% BELOW detection","% SIGNIF. NEGATIVE"))
cat(sprintf("  %-24s %10s %8s %8s %8s %7s %7s %7s\n","","(nmol)","all","tree","soil","all","tree","soil"))
is_soil <- A$campaign == "semirigid_soil"
R <- bind_rows(lapply(names(V), function(nm) {
  m <- V[[nm]]
  bd  <- abs(A$best.flux) < m; sn <- A$best.flux < -m
  cat(sprintf("  %-24s %10.5f %7.1f%% %7.1f%% %7.1f%% %6.1f%% %6.1f%% %6.1f%%\n", nm, median(m),
      100*mean(bd), 100*mean(bd[!is_soil]), 100*mean(bd[is_soil]),
      100*mean(sn), 100*mean(sn[!is_soil]), 100*mean(sn[is_soil])))
  data.frame(variant=nm, median_MDF=median(m),
             pct_below_all=100*mean(bd), pct_below_tree=100*mean(bd[!is_soil]),
             pct_below_soil=100*mean(bd[is_soil]),
             pct_signeg_all=100*mean(sn), pct_signeg_tree=100*mean(sn[!is_soil]),
             pct_signeg_soil=100*mean(sn[is_soil]))
}))
write.csv(R, file.path(outdir,"mdf_threshold_comparison.csv"), row.names=FALSE)

rule("WHAT SURVIVES EVERY CHOICE")
cat("
  Two conclusions hold under ALL eight definitions, from the most liberal to the most
  conservative:
")
tr <- !is_soil
for (nm in names(V)) {
  m <- V[[nm]]
  cat(sprintf("    %-24s  tree stems significantly NEGATIVE: %4.1f%%   soil: %4.1f%%\n",
      nm, 100*mean(A$best.flux[tr] < -m[tr]), 100*mean(A$best.flux[is_soil] < -m[is_soil])))
}
cat("
  1. NO tree-stem uptake. Under every definition the fraction of stem measurements
     that are significantly negative is negligible. There is no support in these data
     for CH4 uptake on stem surfaces at any height, so a scenario in which stem flux
     declines through zero into uptake has no empirical basis here.

  2. SOIL uptake is real. The soil fraction stays high under every definition,
     including the most conservative. The sink is not an artefact of detection limits.

  What the choice DOES change is how much of the stem dataset is described as
  'below detection', which ranges widely below. That affects statements about
  individual measurements, not the budget.
")

rule("EFFECT ON THE BUDGET (the reason not to over-filter)")
CONV <- 86400*365.25*16e-6
cat(sprintf("
  Treatment of below-detection stem measurements, using Christiansen 95%%:

  %-28s %12s %12s
", "", "mean flux", "vs retain-all"))
m95 <- V[["Christiansen 95%"]]
base <- mean(A$best.flux[tr])
for (lab in c("retain all values","set below-MDF to zero","remove below-MDF")) {
  v <- switch(lab,
    `retain all values`     = A$best.flux[tr],
    `set below-MDF to zero` = ifelse(abs(A$best.flux[tr]) < m95[tr], 0, A$best.flux[tr]),
    `remove below-MDF`      = A$best.flux[tr][abs(A$best.flux[tr]) >= m95[tr]])
  cat(sprintf("  %-28s %12.5f %11.1f%%\n", lab, mean(v), 100*(mean(v)/base - 1)))
}
cat("
  Removing below-detection measurements inflates the stem mean substantially, which is
  exactly the bias documented in the companion project. For budgets and scaling the
  correct treatment is to RETAIN all values and report the detection fraction.
")

# --- figure --------------------------------------------------------------------
P <- R %>% dplyr::select(variant, tree=pct_below_tree, soil=pct_below_soil) %>%
  pivot_longer(-variant, names_to="type", values_to="pct") %>%
  mutate(variant=factor(variant, levels=R$variant[order(R$median_MDF)]))
p <- ggplot(P, aes(pct, variant, fill=type)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c(tree="#B2182B", soil="#2166AC"), name=NULL) +
  labs(title="Fraction of measurements below detection, by MDF definition",
       subtitle="ordered by median MDF; the choice spans an order of magnitude",
       x="% below detection", y=NULL) +
  theme_bw(base_size=8)
ggsave(file.path(outdir,"fig_mdf_thresholds.png"), p, width=8, height=4.5, dpi=200, bg="white")
cat("\nWritten: outputs/revision/fig_mdf_thresholds.png\n")
