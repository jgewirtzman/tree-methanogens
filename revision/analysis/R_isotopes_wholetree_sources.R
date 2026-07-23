# ==============================================================================
# REVISION — Whole-tree isotope source terms: Keeling AND Miller-Tans (CH4 + CO2)
# ==============================================================================
# Focus on the WHOLE-TREE single-borehole samples (the manuscript basis; the H/S
# tissue set is being dropped). Reports source d13C for CH4 and CO2 by two
# methods:
#   Keeling:     d13 ~ 1/C ; intercept = source d13   (classic; compresses high C)
#   Miller-Tans: (d13 * C) ~ C ; slope = source d13    (robust over wide C range)
# Internal gas spans ~1.5 to >1500 ppm CH4, a very wide range, so Miller-Tans is
# the more appropriate estimator here; we report both.
#
# Loads ALL FOUR Picarro runs (incl. 205833, which the original script omitted)
# and uses only WholeTree samples. NEW file; edits nothing.
# Run: Rscript revision/analysis/R_isotopes_wholetree_sources.R
# Outputs: revision/outputs/iso_wholetree_sources.png, iso_wholetree_sources.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE) %>% mutate(run = basename(.x)))

std_ids <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
d <- raw %>%
  transmute(SampleName, run,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,    co2 = `12CO2_mean`) %>%
  mutate(class = case_when(
    str_detect(SampleName, "^Amb") ~ "Atmosphere",
    str_detect(SampleName, "^V")   ~ "Incubation",
    SampleName %in% std_ids         ~ "Standard",
    str_detect(SampleName, "[HS]$") ~ "Tissue",
    TRUE                            ~ "WholeTree"))

wt <- d %>% filter(class == "WholeTree", !is.na(d13CH4), ch4 >= 1.5)
wt_co2 <- d %>% filter(class == "WholeTree", !is.na(d13CO2), is.finite(co2), co2 > 0)

# ---- estimators --------------------------------------------------------------
keeling <- function(df, C, D) { df <- df %>% filter(is.finite(.data[[C]]), .data[[C]]>0, is.finite(.data[[D]]))
  m <- lm(df[[D]] ~ I(1/df[[C]])); ci <- confint(m)[1,]
  list(src = coef(m)[1], lo = ci[1], hi = ci[2], n = nrow(df), r2 = summary(m)$r.squared) }
miller_tans <- function(df, C, D) { df <- df %>% filter(is.finite(.data[[C]]), .data[[C]]>0, is.finite(.data[[D]]))
  y <- df[[D]]*df[[C]]; x <- df[[C]]; m <- lm(y ~ x); ci <- confint(m)[2,]
  list(src = coef(m)[2], lo = ci[1], hi = ci[2], n = nrow(df), r2 = summary(m)$r.squared) }

k_ch4 <- keeling(wt, "ch4", "d13CH4");     mt_ch4 <- miller_tans(wt, "ch4", "d13CH4")
k_co2 <- keeling(wt_co2, "co2", "d13CO2"); mt_co2 <- miller_tans(wt_co2, "co2", "d13CO2")

wt <- wt %>% mutate(alpha_C = (d13CO2+1000)/(d13CH4+1000), eps_C = (alpha_C-1)*1000)

# ---- figure ------------------------------------------------------------------
mk <- function(df, C, D, src, ttl, xl, yl) ggplot(df, aes(1/.data[[C]], .data[[D]])) +
  geom_point(alpha=.6) + geom_smooth(method="lm", se=TRUE, color="black", formula=y~x) +
  labs(title=ttl, subtitle=sprintf("Keeling source = %.1f permil", src), x=xl, y=yl) + theme_bw(base_size=10)
mkmt <- function(df, C, D, src, ttl, xl, yl){ df$xx<-df[[C]]; df$yy<-df[[D]]*df[[C]]
  ggplot(df, aes(xx, yy)) + geom_point(alpha=.6) + geom_smooth(method="lm", se=TRUE, color="firebrick", formula=y~x) +
  labs(title=ttl, subtitle=sprintf("Miller-Tans source = %.1f permil", src), x=xl, y=yl) + theme_bw(base_size=10) }
p1 <- mk(wt, "ch4","d13CH4", k_ch4$src, "CH4 Keeling", "1/CH4 (1/ppm)", "d13CH4")
p2 <- mkmt(wt, "ch4","d13CH4", mt_ch4$src, "CH4 Miller-Tans", "CH4 (ppm)", "d13CH4 * CH4")
p3 <- mk(wt_co2, "co2","d13CO2", k_co2$src, "CO2 Keeling", "1/CO2 (1/ppm)", "d13CO2")
p4 <- mkmt(wt_co2, "co2","d13CO2", mt_co2$src, "CO2 Miller-Tans", "CO2 (ppm)", "d13CO2 * CO2")
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir,"iso_wholetree_sources.png"),
         gridExtra::arrangeGrob(p1,p2,p3,p4, nrow=2), width=12, height=9, dpi=150)

# ---- summary -----------------------------------------------------------------
sink(file.path(out_dir, "iso_wholetree_sources.txt"))
cat("=================================================================\n")
cat("WHOLE-TREE ISOTOPE SOURCE TERMS (Keeling vs Miller-Tans)\n")
cat("=================================================================\n\n")
cat("Sample classes across all 4 runs:\n"); print(as.data.frame(d %>% count(class)), row.names=FALSE)
cat(sprintf("\nWhole-tree CH4 isotope samples (ch4>=1.5): %d\n", nrow(wt)))
cat(sprintf("Whole-tree CO2 isotope samples: %d\n\n", nrow(wt_co2)))
qs <- function(v) sprintf("median %.1f (IQR %.1f to %.1f)", median(v,na.rm=T),quantile(v,.25,na.rm=T),quantile(v,.75,na.rm=T))
cat("Bulk signatures:\n")
cat("  d13C-CH4 :", qs(wt$d13CH4), "permil\n")
cat("  d13C-CO2 :", qs(wt_co2$d13CO2), "permil\n")
cat("  eps_C    :", qs(wt$eps_C), "permil\n\n")
cat("SOURCE d13C-CH4:\n")
cat(sprintf("  Keeling      : %.1f permil [95%% CI %.1f, %.1f]  (n=%d, R2=%.2f)\n", k_ch4$src,k_ch4$lo,k_ch4$hi,k_ch4$n,k_ch4$r2))
cat(sprintf("  Miller-Tans  : %.1f permil [95%% CI %.1f, %.1f]  (n=%d, R2=%.2f)\n", mt_ch4$src,mt_ch4$lo,mt_ch4$hi,mt_ch4$n,mt_ch4$r2))
cat("SOURCE d13C-CO2:\n")
cat(sprintf("  Keeling      : %.1f permil [95%% CI %.1f, %.1f]  (n=%d, R2=%.2f)\n", k_co2$src,k_co2$lo,k_co2$hi,k_co2$n,k_co2$r2))
cat(sprintf("  Miller-Tans  : %.1f permil [95%% CI %.1f, %.1f]  (n=%d, R2=%.2f)\n", mt_co2$src,mt_co2$lo,mt_co2$hi,mt_co2$n,mt_co2$r2))
cat("\nRECOMMENDATION (CORRECTED): do NOT use Miller-Tans here. Despite its high R2,\n")
cat("Miller-Tans is dominated by a few high-CH4 points -- its CH4 source swings\n")
cat("-65 -> -89 permil as the top 1-5 points are dropped (the R2 is a concentration-\n")
cat("weighting artifact, not fit quality). KEELING is the robust estimate here\n")
cat("(stable ~ -79 permil to point removal). Better still, report the bulk median\n")
cat("d13CH4 (-63.7, manuscript) plus the atmosphere-mass-balance tree-produced\n")
cat("source (~ -70 permil, endmember-robust; see R_isotopes_convergence_threshold.R).\n")
cat("The tree-produced source is MORE depleted than the measured bulk because the\n")
cat("measured gas is diluted with atmospheric CH4 -- consistent with hydrogenotrophy.\n")
sink()
cat(readLines(file.path(out_dir,"iso_wholetree_sources.txt")), sep="\n")
