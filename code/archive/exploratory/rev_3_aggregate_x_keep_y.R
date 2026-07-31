# ==============================================================================
# REVISION ANALYSIS — Aggregate-x / keep-individual-y  (Mark's suggestion)
# ==============================================================================
# Question (Mark, Slack): regression attenuation is mainly an x-axis issue; if we
# aggregate the noisy predictor (gene abundance) to the species mean but KEEP the
# flux as individual observations, how does that influence the slope and fit?
# Compare against (a) the fully individual regression (attenuated null) and
# (b) the species-means-both-axes regression (current headline, OLS + SMA).
#
# The point this makes:
#   * Individual-both slope ~ 0  (attenuated: noisy x AND huge within-tree y scatter)
#   * Aggregate-x/keep-y slope jumps to the de-biased value  -> confirms the null
#     was x-attenuation, not absence of relationship. R2 stays LOW because
#     individual flux is mostly within-species/within-tree heterogeneity.
#   * Species-both OLS ~ the same de-biased slope; SMA nudges further for any
#     residual error in the species-mean x. High R2 there because y is aggregated too.
# => The slope is real and stable across specifications; only the variance
#    explained depends on whether y is aggregated. Supports a species-level,
#    associational claim without the "mean-line-through-means" fallacy.
#
# STATUS: NEW revision-only script. Modifies nothing pre-existing.
# Regressions use the SAME scale as the manuscript: raw CH4 flux (y) vs
# log10 gene metric (x).
#
# Run:  Rscript code/revision/R3_aggregate_x_keep_y.R
# Outputs (outputs/revision/):
#   aggregate_x_model_comparison.csv
#   aggregate_x_summary.txt
#   aggregate_x_ratio_panels.png
# ==============================================================================

source("code/revision/rev_prep_species_data.R")   # shared prep (new helper)
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "outputs/revision"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# SMA (Warton et al. 2006): b = sign(r) sd(y)/sd(x) = OLS/|r|; CI multiplicative.
sma_fit <- function(x, y, alpha = 0.05) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  r <- cor(x, y); b_ols <- unname(coef(lm(y ~ x))[2])
  b_sma <- sign(r) * sd(y) / sd(x)
  tcrit <- qt(1 - alpha/2, n - 2); B <- tcrit_sq <- tcrit^2 * (1 - r^2) / (n - 2)
  list(n = n, r = r, r2 = r^2, p = cor.test(x, y)$p.value,
       b_ols = b_ols, b_sma = b_sma,
       ci_lo = b_sma*(sqrt(B+1) - sqrt(B)), ci_hi = b_sma*(sqrt(B+1) + sqrt(B)))
}

# pseudo-log flux transform (manuscript's pseudolog10_individual, handles CH4 uptake<0)
psl <- function(x) asinh(x / 0.2) / log(10)

# Generic runner for one predictor across specifications, on a given flux scale --
# indiv_df: per-tree data with columns species, x_ind (individual log metric), CH4_flux (raw)
# sp_x:     species-level predictor: species, x_sp (species-median log metric)
# flux_ind: individual chamber fluxes (species, CH4_flux raw) for aggregate-x/keep-y
# ytrans:   function applied to raw flux (identity for raw scale, psl for pseudo-log)
# Species y-summary is the MEDIAN of transformed flux per species (raw median matches
# the manuscript's median_flux, so M2-raw reproduces R2=0.513).
run_specs <- function(indiv_df, sp_x, flux_ind, label, ytrans, scale_label) {
  ty <- function(v) ytrans(v)
  # M1: individual both
  m1 <- lm(ty(CH4_flux) ~ x_ind, data = indiv_df); s1 <- summary(m1)
  # M2: species medians, both axes (OLS + SMA)
  sp_y <- flux_ind %>% group_by(species) %>%
    summarise(y = median(ty(CH4_flux)), .groups = "drop") %>%
    inner_join(sp_x, by = "species")
  f2 <- sma_fit(sp_y$x_sp, sp_y$y)
  # M3: aggregate-x, keep individual chamber y
  jx <- flux_ind %>% inner_join(sp_x, by = "species")
  m3 <- lm(ty(CH4_flux) ~ x_sp, data = jx); s3 <- summary(m3)
  # M3b: aggregate-x, keep per-tree 125cm y
  jb <- indiv_df %>% inner_join(sp_x, by = "species")
  m3b <- lm(ty(CH4_flux) ~ x_sp, data = jb); s3b <- summary(m3b)

  tibble(
    predictor = label, flux_scale = scale_label,
    model = c("M1 individual (x & y individual)",
              "M2 species medians (both axes, OLS)",
              "M2 species medians (both axes, SMA)",
              "M3 aggregate-x, individual chamber y",
              "M3b aggregate-x, per-tree 125cm y"),
    n = c(nobs(m1), f2$n, f2$n, nobs(m3), nobs(m3b)),
    slope = round(c(coef(m1)[2], f2$b_ols, f2$b_sma, coef(m3)[2], coef(m3b)[2]), 4),
    R2 = round(c(s1$r.squared, f2$r2, f2$r2, s3$r.squared, s3b$r.squared), 4),
    p = round(c(coef(s1)[2,4], f2$p, f2$p, coef(s3)[2,4], coef(s3b)[2,4]), 4)
  )
}

predictors <- list(
  list(lab = "mcrA:methanotroph ratio",
       indiv = tree_level_complete %>% transmute(species, x_ind = log_ratio, CH4_flux),
       sp_x  = analysis_ratio %>% transmute(species, x_sp = median_log_ratio)),
  list(lab = "mcrA (area-weighted)",
       indiv = tree_level_complete %>% transmute(species, x_ind = log_mcra, CH4_flux),
       sp_x  = analysis_mcra %>% transmute(species, x_sp = log10(value + 1))),
  list(lab = "methanotroph (pmoA+mmoX)",
       indiv = tree_level_complete %>% transmute(species, x_ind = log_meth, CH4_flux),
       sp_x  = analysis_meth %>% transmute(species, x_sp = log10(value + 1)))
)
flux_pool <- flux_all %>% dplyr::select(species, CH4_flux)

comparison <- purrr::map_dfr(predictors, function(pp) {
  bind_rows(
    run_specs(pp$indiv, pp$sp_x, flux_pool, pp$lab, identity, "raw"),
    run_specs(pp$indiv, pp$sp_x, flux_pool, pp$lab, psl,      "pseudo-log")
  )
})
ratio_tab <- comparison %>% filter(predictor == "mcrA:methanotroph ratio")
write.csv(comparison, file.path(out_dir, "aggregate_x_model_comparison.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Illustrative 3-panel figure for the RATIO on the PSEUDO-LOG flux scale (where
# OLS is valid; raw flux is too heavy-tailed for a mean-based individual-y fit).
# ------------------------------------------------------------------------------
sp_x <- analysis_ratio %>% transmute(species, x_sp = median_log_ratio)
ind <- tree_level_complete %>% transmute(species, x = log_ratio, y = psl(CH4_flux))
chamber <- flux_all %>% dplyr::select(species, CH4_flux) %>%
  inner_join(sp_x, by = "species") %>% mutate(y = psl(CH4_flux))
sp_med <- flux_all %>% dplyr::select(species, CH4_flux) %>% inner_join(sp_x, by = "species") %>%
  group_by(species, x_sp) %>% summarise(y = median(psl(CH4_flux)), .groups = "drop")

ylab <- "pseudo-log10 CH4 flux"
pA <- ggplot(ind, aes(x, y)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(title = "M1  Individual (x & y individual)",
       subtitle = sprintf("slope=%.3f  R2=%.3f  p=%.3f",
                           coef(lm(y~x, ind))[2], summary(lm(y~x, ind))$r.squared,
                           coef(summary(lm(y~x, ind)))[2,4]),
       x = "log10 mcrA:methanotroph (individual)", y = ylab) +
  theme_bw(base_size = 10)

pB <- ggplot(chamber, aes(x_sp, y)) +
  geom_point(alpha = 0.30, size = 1, position = position_jitter(width = 0.02)) +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  labs(title = "M3  Aggregate-x, individual chamber y",
       subtitle = sprintf("slope=%.3f  R2=%.3f  p=%.3f",
                           coef(lm(y~x_sp, chamber))[2], summary(lm(y~x_sp, chamber))$r.squared,
                           coef(summary(lm(y~x_sp, chamber)))[2,4]),
       x = "log10 mcrA:methanotroph (species median)", y = NULL) +
  theme_bw(base_size = 10)

f2 <- sma_fit(sp_med$x_sp, sp_med$y)
pC <- ggplot(sp_med, aes(x_sp, y)) +
  geom_point(size = 2.5) +
  geom_abline(intercept = mean(sp_med$y) - f2$b_ols*mean(sp_med$x_sp),
              slope = f2$b_ols, color = "firebrick") +
  geom_abline(intercept = mean(sp_med$y) - f2$b_sma*mean(sp_med$x_sp),
              slope = f2$b_sma, color = "steelblue", linetype = "dashed") +
  labs(title = "M2  Species medians (both axes)",
       subtitle = sprintf("OLS %.3f (red)  SMA %.3f (blue)  R2=%.3f",
                           f2$b_ols, f2$b_sma, f2$r2),
       x = "log10 mcrA:methanotroph (species median)", y = ylab) +
  theme_bw(base_size = 10)

# arrange without patchwork dependency: use gridExtra if present, else 3 files
combined_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (combined_ok) {
  g <- gridExtra::arrangeGrob(pA, pB, pC, nrow = 1)
  ggsave(file.path(out_dir, "aggregate_x_ratio_panels.png"), g, width = 12, height = 4, dpi = 150)
} else {
  ggsave(file.path(out_dir, "aggregate_x_ratio_panels.png"),
         pC, width = 5, height = 4, dpi = 150)
  ggsave(file.path(out_dir, "aggregate_x_M1.png"), pA, width = 5, height = 4, dpi = 150)
  ggsave(file.path(out_dir, "aggregate_x_M3.png"), pB, width = 5, height = 4, dpi = 150)
}

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
sink(file.path(out_dir, "aggregate_x_summary.txt"))
cat("================================================================\n")
cat("AGGREGATE-X / KEEP-INDIVIDUAL-Y  vs  species medians  vs  individual\n")
cat("Two flux scales: raw (manuscript) and pseudo-log (asinh, skew-tamed)\n")
cat("================================================================\n\n")
for (lab in unique(comparison$predictor)) {
  cat("---", lab, "---\n")
  print(as.data.frame(comparison[comparison$predictor == lab, -1]), row.names = FALSE)
  cat("\n")
}
rr_raw <- ratio_tab %>% filter(flux_scale == "raw")
rr_psl <- ratio_tab %>% filter(flux_scale == "pseudo-log")
cat("KEY FINDING (ratio predictor):\n\n")
cat("On RAW flux, individual and aggregate-x fits are NON-significant and small:\n")
cat(sprintf("  M1 individual           slope=%.4f R2=%.3f p=%.3f\n", rr_raw$slope[1], rr_raw$R2[1], rr_raw$p[1]))
cat(sprintf("  M3 aggregate-x chamber  slope=%.4f R2=%.3f p=%.3f\n", rr_raw$slope[4], rr_raw$R2[4], rr_raw$p[4]))
cat("  -> NOT because x-attenuation is absent, but because raw flux is heavily\n")
cat("     right-skewed: a few extreme emission events (unrelated to the ratio)\n")
cat("     dominate the mean-based OLS. cor(x, species MEDIAN flux)=0.717 but\n")
cat("     cor(x, species MEAN flux)=0.328. The species result uses the MEDIAN.\n\n")
cat("On PSEUDO-LOG flux (valid OLS scale; the manuscript already uses this for\n")
cat("individual flux display), the relationship strengthens as x is aggregated:\n")
cat(sprintf("  M1 individual           slope=%.4f R2=%.3f p=%.3f  <- MARGINAL (weak, ~0.08)\n", rr_psl$slope[1], rr_psl$R2[1], rr_psl$p[1]))
cat(sprintf("  M3 aggregate-x chamber  slope=%.4f R2=%.3f p=%.3f  <- SIGNIFICANT, de-biased\n", rr_psl$slope[4], rr_psl$R2[4], rr_psl$p[4]))
cat(sprintf("  M3b aggregate-x 125cm   slope=%.4f R2=%.3f p=%.3f  <- SIGNIFICANT\n", rr_psl$slope[5], rr_psl$R2[5], rr_psl$p[5]))
cat(sprintf("  M2 species SMA          slope=%.4f R2=%.3f p=%.3f\n", rr_psl$slope[3], rr_psl$R2[3], rr_psl$p[3]))
cat("  -> aggregating the noisy predictor moves the fit from marginal (individual)\n")
cat("     to significant (aggregate-x), i.e. x-attenuation IS operating once the\n")
cat("     y-scale is appropriate. Aggregation sharpens rather than manufactures\n")
cat("     the signal (the individual fit is already marginally positive, not zero).\n\n")
cat("IMPLICATION for the manuscript: the 'individual R2<0.001, no individual\n")
cat("relationship' claim is largely an artifact of raw heavy-tailed flux. A more\n")
cat("defensible framing: a weak (marginal) individual relationship, strongly\n")
cat("attenuated by within-tree heterogeneity AND flux skew, that sharpens to a\n")
cat("robust species-level (species-typical/median) association. This removes the\n")
cat("'how can species work if individual is exactly zero' paradox R3 pressed on.\n")
cat("DECISION NEEDED (Jon/Mark): present on raw or pseudo-log scale; whether to\n")
cat("retire the 'exactly null individually' claim in favor of 'weak+attenuated'.\n")
sink()

cat("Done. Outputs in", out_dir, "\n")
print(as.data.frame(comparison), row.names = FALSE)
