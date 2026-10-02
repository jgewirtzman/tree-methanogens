source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# figS26_scaling-heatmap.R
# ------------------------------------------------------------------------------
# Summary view of the scaling sensitivity grid.
#
#   a  flux form x WAI, median over the area shapes -- the two axes carrying the
#      headline range
#   b  flux form x BRANCH placement, median over WAI and bole. Branch placement is
#      the least intuitive of the geometry assumptions and its effect is invisible
#      in (a), where it is collapsed into the cell. Bole shape and WAI get no panel
#      of their own: bole shape is negligible (see d) and WAI is a simple
#      multiplier that (a) already shows.
#   c  LEVERAGE as first-order variance indices. The grid is a complete factorial,
#      so this decomposition is exact rather than approximate: for each assumption
#      the variance of its group means over the total variance is the share of
#      spread attributable to it alone, and the remainder is interaction. This
#      replaces a max-minus-min "fold change", which is undefined here because the
#      grid spans zero.
#
# Diverging colour centred on zero throughout: with the bounded negative form the
# tree term can be a net sink, so sign matters and a sequential scale would hide it.
# (2026-10-02: the all-combinations heatmap and the titles were removed; the full
# grid is outputs/data/scaling_full_grid.csv and the archive.) Counts are read from the grid, never typed -- they had
# previously drifted to "280 combinations" and "7 flux forms" against a grid with
# 240 and 6, the same failure as the hardcoded "216" in the Figure 9 caption.
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs"
R <- read.csv(out_path("scaling_full_grid.csv"), stringsAsFactors = FALSE)

FO <- c("constant","exp_band_slope","power","exponential","linear_floored",
        "linear_bounded_median")
R$flux <- factor(R$flux, levels = rev(FO))
# linear_bounded_median is CONTRADICTED by the only direct evidence above 2 m: all
# four climbed-tree measurements there are positive (0.031-0.354 nmol m-2 s-1),
# each larger in magnitude than the (median detected stem uptake) sink it assumes, and the profile's only
# negative value is at 0.5 m where the flux routine flags it below detection.
# Retained because n = 1 tree cannot exclude uptake in other species or conditions,
# but marked so it is not read as an equal member of the set.
CONTRADICTED <- "linear_bounded_median"
R$WAI    <- factor(R$WAI,    levels = sort(unique(R$WAI)))
R$branch <- factor(R$branch, levels = sort(unique(R$branch)))
R$bole   <- factor(R$bole,   levels = sort(unique(R$bole)))

N_COMB <- nrow(R); N_WAI <- nlevels(R$WAI); N_BOLE <- nlevels(R$bole)
N_BRANCH <- nlevels(R$branch); N_FLUX <- nlevels(R$flux)
N_SHAPE <- N_BOLE * N_BRANCH
MEAS <- unique(round(R$measured_mg, 2))[1]
LIM  <- max(abs(R$total_mg)) * c(-1, 1)
star <- function(v) ifelse(v == CONTRADICTED, paste0(v, " *"), as.character(v))

DIV <- scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = LIM,
                            name = expression("mg CH"[4]*" m"^-2*" yr"^-1))
# readable names for the grid's codes (2026-10-02); * marks the form the climbed tree contradicts
FLUX_LAB <- c(constant = "Constant", exp_band_slope = "Per-stem exponential (named)", power = "Power",
              exponential = "Exponential", linear_floored = "Linear, floored at zero",
              linear_bounded_median = "Linear into uptake *")
BRANCH_LAB <- c(uniform_all = "Uniform, whole stem", uniform_top50 = "Uniform, upper half",
                gaussian_50 = "Gaussian at 50% height", gaussian_75 = "Gaussian at 75% height")
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
# (a variance share was rejected: the flux axis takes 89% and WAI, which doubles the
#  answer, reads as 3%). For each assumption, hold every other fixed and take the
# range across its own levels; over every combination of the others this gives a
# distribution of ranges, shown as a box.
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
