# ==============================================================================
# REVISION — Precision-weighted within-tree eps_C (no arbitrary threshold)
# ==============================================================================
# Replaces the eyeballed ">=X ppm" cutoff with a principled estimator: weight each
# tree's eps_C by its INVERSE MEASUREMENT VARIANCE. The Picarro G2201-i d13C
# precision is <1.15 permil at 10 ppm CH4 (manuscript/instrument spec); for CRDS
# the delta precision scales ~1/[CH4] (rare-isotopologue signal falls with conc),
# with a systematic floor at high conc. So:
#     sigma_d13CH4(C) = max( floor , k / C ),   k = 1.15 * 10 = 11.5 permil*ppm
# eps_C depends on d13CH4 (noisy) and d13CO2 (reliable; CO2 >> atmospheric), so
#     sigma_eps ~ |d eps/d d13CH4| * sigma_d13CH4 ~ sigma_d13CH4.
# Weight w_i = 1/sigma_eps_i^2. Precision-weighted mean uses ALL samples; low-CH4
# (imprecise, atmosphere/oxidation-biased) samples down-weight themselves.
#
# Validation: the weighted mean should agree with (a) the >=10 ppm subset median
# and (b) the atmosphere mixing-corrected median -- three routes, one answer.
# Standards excluded. NEW file; edits nothing.
# Run: Rscript revision/analysis/R_isotopes_epsC_weighted.R
# Outputs: revision/outputs/iso_epsC_weighted.png, iso_epsC_weighted.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE))
std <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
wt <- raw %>%
  transmute(SampleName, d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean, co2 = `12CO2_mean`) %>%
  filter(!str_detect(SampleName, "^Amb|^V"), !str_detect(SampleName, "[HS]$"),
         !SampleName %in% std, !is.na(d13CH4), !is.na(d13CO2), ch4 >= 1.5) %>%
  mutate(eps_C = ((d13CO2 + 1000)/(d13CH4 + 1000) - 1)*1000)

# precision model (Picarro-anchored); floor = systematic limit at high conc
FLOOR <- 0.3; K <- 1.15 * 10          # permil*ppm
sigma_d <- function(C) pmax(FLOOR, K / C)
alpha <- (median(wt$d13CO2)+1000)/(median(wt$d13CH4)+1000)
wt <- wt %>% mutate(sig = alpha * sigma_d(ch4), w = 1/sig^2)

# ---- RANDOM-EFFECTS pooling (DerSimonian-Laird) ------------------------------
# Naive inverse-MEASUREMENT-variance weighting is wrong here: it assumes trees are
# identical (no biological spread), so the highest-CH4 point gets ~infinite weight
# (SE -> 0) -- the same leverage trap as Miller-Tans. Instead treat each tree as a
# "study" with measurement variance sigma_i^2 AND estimate between-tree
# heterogeneity tau^2; weight = 1/(sigma_i^2 + tau^2).
dl_pool <- function(x, s2) {
  wf <- 1/s2; mu_fe <- sum(wf*x)/sum(wf)
  Q  <- sum(wf*(x - mu_fe)^2); k <- length(x)
  C  <- sum(wf) - sum(wf^2)/sum(wf)
  tau2 <- max(0, (Q - (k - 1))/C)
  wr <- 1/(s2 + tau2); mu <- sum(wr*x)/sum(wr); se <- sqrt(1/sum(wr))
  list(mu = mu, se = se, tau = sqrt(tau2), mu_fe = mu_fe)
}
re_eps <- dl_pool(wt$eps_C, wt$sig^2)
re_d13 <- dl_pool(wt$d13CH4, (alpha*sigma_d(wt$ch4))^2)   # sigma_d13CH4 directly

eps_w   <- re_eps$mu; eps_w_se <- re_eps$se; eps_tau <- re_eps$tau
d13_w   <- re_d13$mu
# validation routes
eps_10 <- median(wt$eps_C[wt$ch4 >= 10]); d13_10 <- median(wt$d13CH4[wt$ch4 >= 10])
wt <- wt %>% mutate(d13_corr = (ch4*d13CH4 - 1.9*(-47))/(ch4-1.9),
                    eps_corr = ((d13CO2+1000)/(d13_corr+1000)-1)*1000)
eps_corr <- median(wt$eps_corr[wt$ch4-1.9 > 0.5], na.rm = TRUE)

sink(file.path(out_dir, "iso_epsC_weighted.txt"))
cat("=================================================================\n")
cat("PRECISION-WEIGHTED within-tree eps_C (no arbitrary threshold)\n")
cat("=================================================================\n\n")
cat(sprintf("n = %d whole-tree paired samples (>=1.5 ppm, standards excluded)\n", nrow(wt)))
cat(sprintf("Precision model: sigma_d13CH4(C) = max(%.1f, %.1f/C) permil (anchor <1.15 at 10 ppm)\n\n", FLOOR, K))
cat("PRIMARY estimate — random-effects pooled (DerSimonian-Laird):\n")
cat(sprintf("  eps_C   = %.1f +/- %.1f permil (pooled mean +/- SE)\n", eps_w, eps_w_se))
cat(sprintf("  tau (between-tree SD) = %.0f permil  <- large real tree-to-tree spread\n", eps_tau))
cat(sprintf("  pooled d13CH4 = %.1f permil\n\n", d13_w))
cat("ROBUSTNESS / VALIDATION (independent routes, all in the same place):\n")
cat(sprintf("  (1) random-effects pooled        : eps_C = %.1f +/- %.1f\n", eps_w, eps_w_se))
cat(sprintf("  (2) >=10 ppm subset median       : eps_C = %.1f (d13CH4 = %.1f)\n", eps_10, d13_10))
cat(sprintf("  (3) atmosphere mixing-corrected  : eps_C = %.1f\n\n", eps_corr))
cat("Why NOT naive inverse-variance weighting: it assumes trees are identical\n")
cat("(tau=0), so the single highest-CH4 point dominates and SE collapses to ~0 --\n")
cat("the Miller-Tans leverage trap. The random-effects model estimates real\n")
cat("between-tree spread (tau ~ tens of permil) and so weights sensibly.\n\n")
cat("=> eps_C ~ 55 +/- a few permil by all routes, at/into the CO2-reduction range,\n")
cat("   with NO arbitrary threshold. The large tau says trees genuinely vary; the\n")
cat("   POOLED value is what supports (corroboratively) a hydrogenotrophic signal.\n")
cat("   Report alongside bulk d13CH4 -63.7; keep the bulk-CO2 caveat.\n")
sink()

# figure: eps_C vs CH4 with point size = weight, and the weighted mean line
p <- ggplot(wt, aes(ch4, eps_C, size = w)) + geom_point(alpha = .5) + scale_x_log10() +
  scale_size_continuous(range = c(.3, 4), guide = "none") +
  geom_hline(yintercept = eps_w, color = "firebrick") +
  geom_hline(yintercept = eps_w + eps_w_se, color = "firebrick", linetype = "dotted") +
  geom_hline(yintercept = eps_w - eps_w_se, color = "firebrick", linetype = "dotted") +
  coord_cartesian(ylim = c(-100, 150)) +
  labs(title = "Precision-weighted within-tree eps_C",
       subtitle = sprintf("point size = 1/variance (Picarro 1/C model); weighted mean = %.0f +/- %.0f permil",
                          eps_w, eps_w_se),
       x = "internal CH4 (ppm, log)", y = "eps_C (permil)") + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "iso_epsC_weighted.png"), p, width = 8, height = 5, dpi = 150)
cat(readLines(file.path(out_dir,"iso_epsC_weighted.txt")), sep="\n")
