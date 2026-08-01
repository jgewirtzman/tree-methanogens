#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_heatmap.R
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
#   c  ALL combinations, nothing collapsed: flux form against WAI nested within
#      branch placement, split by bole shape. Every cell of the grid appears once,
#      so the reader sees the thing itself rather than a summary of it.
#   d  LEVERAGE as first-order variance indices. The grid is a complete factorial,
#      so this decomposition is exact rather than approximate: for each assumption
#      the variance of its group means over the total variance is the share of
#      spread attributable to it alone, and the remainder is interaction. This
#      replaces a max-minus-min "fold change", which is undefined here because the
#      grid spans zero.
#
# Diverging colour centred on zero throughout: with the bounded negative form the
# tree term can be a net sink, so sign matters and a sequential scale would hide it.
# Every count and value in the titles is read from the grid, never typed -- they had
# previously drifted to "280 combinations" and "7 flux forms" against a grid with
# 240 and 6, the same failure as the hardcoded "216" in the Figure 9 caption.
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"
R <- read.csv(file.path(outdir, "scaling_full_grid.csv"), stringsAsFactors = FALSE)

FO <- c("constant","exp_band_slope","power","exponential","linear_floored",
        "linear_bounded_median")
R$flux <- factor(R$flux, levels = rev(FO))
# linear_bounded_median is CONTRADICTED by the only direct evidence above 2 m: all
# four climbed-tree measurements there are positive (0.031-0.354 nmol m-2 s-1),
# each larger in magnitude than the -0.0262 sink it assumes, and the profile's only
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
th <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.4, colour = "grey30"),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 6.2),
        axis.text.y = element_text(size = 7))

# ---- a: flux x WAI -----------------------------------------------------------
A <- R %>% group_by(flux, WAI) %>% summarise(m = median(total_mg), .groups = "drop")
pa <- ggplot(A, aes(WAI, flux, fill = m)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", m)), size = 2.7) + DIV +
  scale_y_discrete(labels = star) +
  labs(title = "a   flux form x woody area index",
       subtitle = sprintf("median over the %d area shapes; blue = net sink; * contradicted by the climbed tree", N_SHAPE),
       x = NULL, y = NULL) + th

# ---- b: flux x branch placement ---------------------------------------------
B <- R %>% group_by(flux, branch) %>% summarise(m = median(total_mg), .groups = "drop")
pb <- ggplot(B, aes(branch, flux, fill = m)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f", m)), size = 2.7) + DIV +
  scale_y_discrete(labels = star) +
  labs(title = "b   flux form x branch placement",
       subtitle = "median over WAI and bole shape; the axis collapsed inside (a)",
       x = NULL, y = NULL) + th

# ---- c: ONE heatmap, every combination, one scale ---------------------------
# Previously this was facet_grid(bole ~ branch), which ggplot renders as eight
# panels. Faceting was the wrong structure: the eye cannot compare a cell in one
# panel with a cell in another, so the "all combinations" panel showed less than
# the two collapsed panels above it. Now a single heatmap on one scale, with the
# two structural axes CROSSED into the rows and columns rather than split across
# panels: rows are flux form x bole shape, columns are woody area index x branch
# placement. Every one of the combinations is a cell, and any cell is directly
# comparable to any other.
CC <- R %>% mutate(
  row = interaction(flux, bole, sep = "  |  ", lex.order = TRUE),
  col = interaction(WAI, branch, sep = "  |  ", lex.order = TRUE))
# order rows by flux (already levelled) then bole; columns by WAI then branch
CC$row <- factor(CC$row, levels = rev(levels(CC$row)))
pc <- ggplot(CC, aes(col, row, fill = total_mg)) +
  geom_tile(colour = "white", linewidth = 0.35) + DIV +
  labs(title = sprintf("c   all %d combinations on one scale", N_COMB),
       subtitle = "rows = flux form x bole shape; columns = woody area index x branch placement",
       x = NULL, y = NULL) +
  th + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 4.6),
             axis.text.y = element_text(size = 5.2),
             legend.position = "right")

# ---- d: leverage as the range each assumption actually spans ------------------
# A variance share was tried and rejected. On this grid the flux axis takes 89% of
# the variance and everything else reads as 0-3%, which is arithmetically true and
# substantively misleading: WAI is a multiplier that DOUBLES the answer, and "3%"
# reads as "irrelevant". The problem is that a share of total variance is dominated
# by whichever factor happens to have the widest range, so it cannot express what a
# reader wants to know -- how much does THIS assumption move the answer.
#
# Instead: for each assumption, hold every other assumption fixed and take the
# range it spans across its own levels. Doing that for every combination of the
# others gives a DISTRIBUTION of ranges rather than one number, which also shows
# how much each factor's leverage depends on the others (an interaction, made
# visible instead of relegated to a residual).
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
LEVS$assumption <- factor(LEVS$assumption, levels = LEVS$assumption)
LEV$assumption  <- factor(LEV$assumption,  levels = levels(LEVS$assumption))

pd <- ggplot(LEV, aes(assumption, range_mg)) +
  geom_boxplot(width = 0.5, outlier.size = 0.7, fill = "#d1e5f0",
               colour = "grey35", linewidth = 0.4) +
  geom_text(data = LEVS, aes(assumption, mx,
              label = sprintf("  median %.0f mg%s", med,
                ifelse(is.na(med_fold), "", sprintf("   (%.1fx)", med_fold)))),
            hjust = 0, size = 2.9) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, max(LEV$range_mg)*1.55), expand = c(0, 0)) +
  labs(title = "d   leverage: how far each assumption moves the answer",
       subtitle = "range spanned by that assumption with all others held fixed, over every combination of the others;\nbox = IQR across those combinations. Fold change shown where the range stays positive.",
       x = NULL, y = expression("range in whole-surface total (mg CH"[4]*" m"^-2*" yr"^-1*")")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.4, colour = "grey30"),
        plot.margin = margin(6, 92, 6, 6))

fig <- (pa | pb) / pc / pd + plot_layout(heights = c(1, 1.85, 0.62)) +
  plot_annotation(
    title = sprintf("Stem CH4 upscaling: %d combinations", N_COMB),
    subtitle = sprintf("%d woody area index x %d bole shapes x %d branch placements x %d vertical flux forms   |   measured 0-2 m band = %.2f mg CH4 m-2 yr-1 in every combination",
                       N_WAI, N_BOLE, N_BRANCH, N_FLUX, MEAS),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.4, colour = "grey25")))
ggsave(file.path(outdir, "fig_scaling_heatmap.png"), fig,
       width = 13, height = 15.5, dpi = 200, bg = "white")

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
cat("\nwritten: outputs/revision/fig_scaling_heatmap.png\n")
