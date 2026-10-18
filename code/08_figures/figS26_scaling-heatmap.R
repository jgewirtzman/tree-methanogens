#!/usr/bin/env Rscript
# ==============================================================================
# figS26_scaling-heatmap.R -- Figure S29: whole-woody-surface stem emission across the
# scaling grid
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) median total for each flux form above 2 m x woody area index, over the area
#       shapes
#   (b) median total for each flux form x branch-area placement, over WAI and bole
#       shape (branch placement is collapsed into the cells of (a))
#   (c) leverage: for each assumption, the range of the total across its own levels
#       with every other assumption held fixed, over all combinations of the others
#       (box). A variance share is not used: it would rank WAI, which doubles the
#       answer, at a few percent.
# Diverging colour centred on zero: with the uptake form the stem term can be a net
# sink. Counts are read from the grid, never typed.
# Input: outputs/data/scaling_full_grid.csv (06_upscale/07_scaling_full_grid.R)
# Outputs: outputs/figures/generated/fig_scaling_heatmap.png,
#          outputs/data/scaling_leverage.csv (both leverage statistics, for the text)
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
source("code/lib/outputs.R")
R <- read.csv(out_path("scaling_full_grid.csv"), stringsAsFactors = FALSE)

FO <- c("constant","exp_band_slope","power","exponential","linear_floored",
        "linear_bounded_median")
stopifnot(setequal(FO, unique(R$flux)))
R$flux <- factor(R$flux, levels = rev(FO))
# linear_bounded_median (linear into uptake) is contradicted by the only direct evidence
# above 2 m: every climbed-tree measurement there is positive and larger than the sink it
# assumes. Kept, because one tree cannot exclude uptake elsewhere, but marked "*".
R$WAI    <- factor(R$WAI,    levels = sort(unique(R$WAI)))
R$branch <- factor(R$branch, levels = sort(unique(R$branch)))
R$bole   <- factor(R$bole,   levels = sort(unique(R$bole)))

N_COMB <- nrow(R)
LIM  <- max(abs(R$total_mg)) * c(-1, 1)

DIV <- scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = LIM,
                            name = expression("mg CH"[4]*" m"^-2*" yr"^-1))
# readable names for the grid's codes; * marks the form the climbed tree contradicts
FLUX_LAB <- c(constant = "Constant", exp_band_slope = "Per-stem exponential", power = "Power",
              exponential = "Exponential", linear_floored = "Linear, floored at zero",
              linear_bounded_median = "Linear into uptake *")
BRANCH_LAB <- c(uniform_all = "Uniform, whole stem", uniform_top50 = "Uniform, upper half",
                gaussian_50 = "Gaussian at 50% height", gaussian_75 = "Gaussian at 75% height")
stopifnot(setequal(names(FLUX_LAB), FO), setequal(names(BRANCH_LAB), levels(R$branch)))
WAI_LAB <- function(v) gsub(" \\((.*)\\)", ", \\1", sub("^([0-9.]+) (.*)$", "\\1\n\\2", v))
th <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9.5),
        axis.text.y = element_text(size = 10))

# ---- a: flux x WAI -----------------------------------------------------------
A <- R %>% group_by(flux, WAI) %>% summarise(m = median(total_mg), .groups = "drop")
pa <- ggplot(A, aes(WAI, flux, fill = m)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", m)), size = 3.6) + DIV +
  scale_y_discrete(labels = FLUX_LAB) + scale_x_discrete(labels = WAI_LAB) +
  labs(x = "Woody area index", y = "Flux form above 2 m") + th

# ---- b: flux x branch placement ---------------------------------------------
B <- R %>% group_by(flux, branch) %>% summarise(m = median(total_mg), .groups = "drop")
pb <- ggplot(B, aes(branch, flux, fill = m)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", m)), size = 3.6) + DIV +
  scale_y_discrete(labels = NULL) + scale_x_discrete(labels = BRANCH_LAB) +
  labs(x = "Branch-area placement", y = NULL) + th

# ---- c: leverage, the range each assumption spans with the others held fixed -----
FACT <- c("flux","WAI","branch","bole")
LEV <- bind_rows(lapply(FACT, function(v) {
  others <- setdiff(FACT, v)
  R %>% group_by(across(all_of(others))) %>%
    summarise(range_mg = max(total_mg) - min(total_mg),
              lo = min(total_mg), hi = max(total_mg), .groups = "drop") %>%
    mutate(assumption = v, n_levels = nlevels(R[[v]]),
           fold = ifelse(lo > 0, hi/lo, NA_real_))
}))
LEVS <- LEV %>% group_by(assumption, n_levels) %>%
  summarise(med = median(range_mg), q1 = quantile(range_mg, .25),
            q3 = quantile(range_mg, .75), mn = min(range_mg), mx = max(range_mg),
            med_fold = suppressWarnings(median(fold, na.rm = TRUE)), .groups = "drop") %>%
  arrange(med)
ASSUMP_LAB <- c(flux = "Flux form above 2 m", WAI = "Woody area index", branch = "Branch placement", bole = "Bole shape")
LEVS$assumption <- factor(LEVS$assumption, levels = LEVS$assumption)
LEV$assumption  <- factor(LEV$assumption,  levels = levels(LEVS$assumption))
pc <- ggplot(LEV, aes(assumption, range_mg)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, fill = "#d1e5f0", colour = "grey35", linewidth = 0.4) +
  coord_flip() + scale_x_discrete(labels = ASSUMP_LAB) +
  labs(x = NULL, y = expression("Range in whole-surface total (mg CH"[4]*" m"^-2*" yr"^-1*")")) +
  theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10.5))

fig <- (pa | pb) / pc + plot_layout(heights = c(1, 0.42), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 14, face = "bold"), legend.position = "right")
ggsave(out_path("fig_scaling_heatmap.png"), fig, width = 13, height = 9, dpi = 300, bg = "white")

# Both leverage statistics, for the text: the range of the median total across an
# assumption's levels (main text, Methods S7), and the median range with all other
# assumptions held fixed (panel c).
RM <- sapply(FACT, function(v) { m <- tapply(R$total_mg, R[[v]], median); max(m) - min(m) })
write.csv(data.frame(assumption = FACT, range_of_medians_mg = round(RM[FACT], 2),
                     median_range_others_fixed_mg = round(LEVS$med[match(FACT, LEVS$assumption)], 2),
                     iqr_lo = round(LEVS$q1[match(FACT, LEVS$assumption)], 2),
                     iqr_hi = round(LEVS$q3[match(FACT, LEVS$assumption)], 2)),
          out_path("scaling_leverage.csv"), row.names = FALSE)

cat(sprintf("=== LEVERAGE: range spanned by each assumption, others held fixed (%d combinations) ===\n", N_COMB))
print(as.data.frame(LEVS %>% arrange(-med) %>%
      transmute(assumption, levels = n_levels,
                median_range_mg = round(med,1), IQR = sprintf("%.1f-%.1f", q1, q3),
                full_range = sprintf("%.1f-%.1f", mn, mx),
                median_fold = ifelse(is.na(med_fold), "-", sprintf("%.2fx", med_fold)))),
      row.names = FALSE)
cat("\n=== branch placement, the axis (a) hides ===\n")
print(as.data.frame(R %>% group_by(branch) %>%
  summarise(median_total = round(median(total_mg), 1),
            lo = round(min(total_mg), 1), hi = round(max(total_mg), 1), .groups = "drop")),
  row.names = FALSE)
cat("\nwritten: outputs/figures/generated/fig_scaling_heatmap.png\n")
