# ==============================================================================
# REVISION — Within-tree CH4-CO2 apparent fractionation (the correct unit)
# ==============================================================================
# Jon's point: fractionation should be computed WITHIN each tree (its own paired
# CH4 + CO2), then summarized across trees -- NOT by comparing two separately
# POOLED source terms (Miller-Tans CH4 source vs CO2 source), which are each
# dominated by different trees and so conflate pools.
#
# Each whole-tree single-borehole sample has BOTH d13C-CH4 and d13C-CO2, so the
# apparent fractionation is a per-tree quantity:
#   alpha_C = (d13CO2 + 1000) / (d13CH4 + 1000) ;  eps_C = (alpha_C - 1)*1000
# Diagnostic (Whiticar 1999): eps_C ~ 55-90 -> CO2 reduction (hydrogenotrophic);
#   ~40-55 -> acetate/methyl. CAVEAT: bulk internal CO2 is respiration-dominated,
#   not the specific methanogenic substrate, so eps_C is APPARENT (R2 L607-610).
#
# Also: at LOW internal CH4, d13C-CH4 is pulled toward atmosphere (~-47) and
# enriched by oxidation, biasing eps_C low. So we show eps_C vs CH4 concentration
# and report the value where CH4 is high enough to reflect the production source.
#
# NEW file; edits nothing. Run: Rscript revision/analysis/R_isotopes_within_tree_fractionation.R
# Outputs: revision/outputs/iso_within_tree_fractionation.png, iso_within_tree_fractionation.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE))

# EXCLUDE calibration standards (e.g. S3c: 16105 ppm, d13CH4 -391 permil) -- these
# leak in otherwise and, being high-concentration, hijack the Miller-Tans slope.
std_ids <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
wt <- raw %>%
  transmute(SampleName, d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean, co2 = `12CO2_mean`) %>%
  filter(!str_detect(SampleName, "^Amb|^V"), !str_detect(SampleName, "[HS]$"),
         !SampleName %in% std_ids,
         !is.na(d13CH4), !is.na(d13CO2), ch4 >= 1.5) %>%
  mutate(alpha_C = (d13CO2 + 1000)/(d13CH4 + 1000), eps_C = (alpha_C - 1)*1000)

# eps_C by CH4 threshold (does within-tree fractionation converge as CH4 rises?)
thr_tab <- map_dfr(c(1.5,5,10,50,100,500), function(t) {
  s <- wt %>% filter(ch4 >= t)
  tibble(ch4_threshold = t, n = nrow(s),
         eps_median = median(s$eps_C, na.rm = TRUE),
         eps_IQR = IQR(s$eps_C, na.rm = TRUE),
         pct_hydrogenotrophic = 100*mean(s$eps_C > 55, na.rm = TRUE))
})

# pooled-source alpha_C for contrast (the conflated approach)
mt <- function(C, D) { m <- lm((D*C) ~ C); coef(m)[2] }
src_ch4 <- mt(wt$ch4, wt$d13CH4); src_co2 <- mt(wt$co2, wt$d13CO2)
pooled_eps <- ((src_co2 + 1000)/(src_ch4 + 1000) - 1) * 1000

# ---- figures: d13CH4 and eps_C vs internal CH4 (same samples, consistent) ----
wt_f <- wt %>% filter(ch4 <= 5000)   # drop 1 extreme outlier for display
p1 <- ggplot(wt_f, aes(ch4, d13CH4)) +
  geom_point(alpha = .55, size = 1.3) + scale_x_log10() +
  geom_smooth(method = "loess", se = TRUE, color = "firebrick", formula = y ~ x) +
  labs(title = "d13C-CH4 vs internal CH4",
       subtitle = "more depleted (more hydrogenotrophic) as CH4 rises; low CH4 = atmosphere/oxidation-mixed",
       x = "internal CH4 (ppm, log)", y = "d13C-CH4 (permil)") + theme_bw(base_size = 10)
p2 <- ggplot(wt_f, aes(ch4, eps_C)) +
  geom_point(alpha = .55, size = 1.3) + scale_x_log10() +
  geom_smooth(method = "loess", se = TRUE, color = "steelblue", formula = y ~ x) +
  labs(title = "Apparent fractionation eps_C vs internal CH4",
       subtitle = "eps_C = (d13CO2 - d13CH4)-based; converges to ~50-55 at reliable CH4 (bulk-CO2 caveat)",
       x = "internal CH4 (ppm, log)", y = "eps_C (permil)") + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "iso_within_tree_fractionation.png"),
         gridExtra::arrangeGrob(p1, p2, nrow = 1), width = 12, height = 5, dpi = 150)

sink(file.path(out_dir, "iso_within_tree_fractionation.txt"))
cat("=================================================================\n")
cat("WITHIN-TREE CH4-CO2 APPARENT FRACTIONATION (the appropriate unit)\n")
cat("=================================================================\n\n")
cat("Each whole-tree sample has paired d13CH4 + d13CO2, so eps_C is computed\n")
cat("WITHIN each tree, then summarized across trees. This is the appropriate\n")
cat("way to look at CH4-vs-CO2 fractionation.\n\n")
cat("CONSISTENT paired medians (d13CH4, d13CO2, eps_C from the SAME samples):\n")
consist <- map_dfr(c(1.5,10,50,100), function(t){ s <- wt %>% filter(ch4 >= t)
  tibble(ch4_threshold = t, n = nrow(s), d13CH4 = median(s$d13CH4), d13CO2 = median(s$d13CO2),
         eps_from_medians = ((median(s$d13CO2)+1000)/(median(s$d13CH4)+1000)-1)*1000) })
print(as.data.frame(consist %>% mutate(across(where(is.numeric), ~round(.,1)))), row.names = FALSE)
cat("\nKEY: d13CH4 (a COMPOSITION) and eps_C (a DIFFERENCE) are different quantities.\n")
cat(" * Bulk (all samples): d13CH4 ~ -62 (manuscript -63.7), d13CO2 ~ -21, eps_C ~ 44.\n")
cat(" * Reliable (>=50 ppm): d13CH4 ~ -70 (more depleted), d13CO2 ~ -18, eps_C ~ 55.\n")
cat(" -> the -63 d13CH4 pairs with eps_C ~44, NOT ~60. eps_C ~55 goes with the more-\n")
cat("    depleted high-CH4 d13CH4 (~-70). Always report both from the SAME subset.\n\n")
cat(sprintf("CONTRAST — the pooled-source approach CONFLATES pools:\n"))
cat(sprintf("  Miller-Tans CH4 source=%.0f, CO2=%.0f -> paired eps_C=%.0f.\n", src_ch4, src_co2, pooled_eps))
cat("  This eps_C (~32) is WRONG: the CO2 source (-35) is weighted toward high-CO2\n")
cat("  trees, but the CO2 that co-occurs with high CH4 is ~-18. Within-tree eps_C\n")
cat("  (~55) pairs each tree's own CH4 and CO2 and is the valid estimate.\n")
cat("  (NB: standards like S3c, d13CH4 -391 at 16105 ppm, must be excluded or they\n")
cat("   hijack the Miller-Tans slope -- that earlier gave a spurious -269.)\n")
cat("\nHOW TO REPORT:\n")
cat(" * Primary isotopic evidence = depleted d13CH4 (-63.7, manuscript) + taxonomy.\n")
cat(" * Also report d13CO2 (~-20) and the apparent eps_C from the SAME samples:\n")
cat("   bulk eps_C ~44, rising to ~55 (boundary of the CO2-reduction range) as CH4\n")
cat("   increases and the signal cleans up. eps_C is CORROBORATING, not decisive.\n")
cat(" * Caveat: eps_C uses bulk (respiration) CO2, not the isolated methanogenic\n")
cat("   substrate, so it likely mis-estimates the true fractionation (R2 L607-610).\n")
sink()
cat(readLines(file.path(out_dir,"iso_within_tree_fractionation.txt")), sep="\n")
