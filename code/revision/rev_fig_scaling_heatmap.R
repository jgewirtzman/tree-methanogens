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

FO <- c("constant","rf_learned_capped","power","exponential","linear_floored",
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

# ---- c: every combination ----------------------------------------------------
pc <- ggplot(R, aes(WAI, flux, fill = total_mg)) +
  geom_tile(colour = "white", linewidth = 0.25) + DIV +
  scale_y_discrete(labels = star) +
  facet_grid(bole ~ branch) +
  labs(title = sprintf("c   all %d combinations, nothing collapsed", N_COMB),
       subtitle = "columns = woody area index within branch placement; rows = bole shape",
       x = NULL, y = NULL) +
  th + theme(strip.text = element_text(size = 7),
             axis.text.x = element_text(angle = 55, hjust = 1, size = 5),
             axis.text.y = element_text(size = 6))

# ---- d: leverage, first-order variance indices -------------------------------
vt <- var(R$total_mg); gm <- mean(R$total_mg)
S <- bind_rows(lapply(c("flux","WAI","branch","bole"), function(v) {
  g <- R %>% group_by(.data[[v]]) %>% summarise(m = mean(total_mg), n = dplyr::n(), .groups = "drop")
  data.frame(assumption = v,
             S1 = sum(g$n * (g$m - gm)^2)/(nrow(R) - 1)/vt,
             lo = min(g$m), hi = max(g$m), spread = max(g$m) - min(g$m))
}))
S$assumption <- factor(S$assumption, levels = S$assumption[order(S$S1)])
inter <- 1 - sum(S$S1)
# A variance share is a share of the WHOLE grid's spread, which the flux axis
# dominates, so WAI's index understates what it does in practice: it is a
# multiplier, and within any one flux form it moves the total by a factor of two.
# Both readings are reported, because quoting only the 3% would be misleading.
wai_within <- R %>% group_by(flux) %>%
  summarise(lo = min(total_mg), hi = max(total_mg), .groups = "drop") %>%
  mutate(fold = ifelse(lo > 0, hi/lo, NA_real_))
wf <- wai_within$fold[wai_within$flux == "constant"]
pd <- ggplot(S, aes(assumption, S1)) +
  geom_col(fill = "#4393c3", width = 0.62) +
  geom_text(aes(label = sprintf("%.0f%%    spread %.0f mg", 100*S1, spread)),
            hjust = -0.05, size = 2.9) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15), expand = c(0, 0)) +
  labs(title = "d   leverage: share of grid variance attributable to each assumption",
       subtitle = sprintf("first-order indices on a complete factorial, so exact; interactions %.0f%%.\nA share of the WHOLE grid understates WAI: it is a multiplier, spanning %.1fx within the constant-flux row alone.",
                          100*inter, wf),
       x = NULL, y = "share of variance in the whole-surface total") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.4, colour = "grey30"),
        plot.margin = margin(6, 74, 6, 6))

fig <- (pa | pb) / pc / pd + plot_layout(heights = c(1, 1.3, 0.72)) +
  plot_annotation(
    title = sprintf("Stem CH4 upscaling: %d combinations", N_COMB),
    subtitle = sprintf("%d woody area index x %d bole shapes x %d branch placements x %d vertical flux forms   |   measured 0-2 m band = %.2f mg CH4 m-2 yr-1 in every combination",
                       N_WAI, N_BOLE, N_BRANCH, N_FLUX, MEAS),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.4, colour = "grey25")))
ggsave(file.path(outdir, "fig_scaling_heatmap.png"), fig,
       width = 11.5, height = 13.5, dpi = 200, bg = "white")

cat(sprintf("=== LEVERAGE (first-order variance indices, %d combinations) ===\n", N_COMB))
print(as.data.frame(S %>% arrange(-S1) %>%
      transmute(assumption, share = sprintf("%.1f%%", 100*S1),
                mean_lo = round(lo, 1), mean_hi = round(hi, 1),
                spread_mg = round(spread, 1))), row.names = FALSE)
cat(sprintf("  interactions: %.1f%%\n", 100*inter))
cat("\n  NOTE: the variance share is a share of the whole grid, which the flux axis\n")
cat("  dominates. WAI is a multiplier; within a single flux form it spans:\n")
print(as.data.frame(wai_within %>%
  transmute(flux, lo = round(lo,1), hi = round(hi,1),
            fold = ifelse(is.na(fold), "-", sprintf("%.2fx", fold)))), row.names = FALSE)

cat("\n=== branch placement, the axis (a) hides ===\n")
print(as.data.frame(R %>% group_by(branch) %>%
  summarise(median_total = round(median(total_mg), 1),
            lo = round(min(total_mg), 1), hi = round(max(total_mg), 1), .groups = "drop")),
  row.names = FALSE)
cat("\nwritten: outputs/revision/fig_scaling_heatmap.png\n")
