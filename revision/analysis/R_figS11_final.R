# ==============================================================================
# REVISION — S11 rebuilt: 3-column scale-dependent gene-flux figure.
# Columns = aggregation ladder:
#   M1 Individual        : individual x (gene), individual y (flux)
#   M3 Aggregate predictor: species-median x (gene), INDIVIDUAL y (flux)   [NEW]
#   M2 Species-level     : species-median x & y
# Rows = 5 predictors (mcrA, pmoA, mmoX, pmoA+mmoX, ratio). Bottom row = R2 bars.
#
# Every panel reports slope + R2 + p from a GENE-ONLY simple linear model fit to the
# plotted response (so the reported stats match the drawn line, and the R2 bars are
# the gene's R2 — NOT the old full species+gene model). M1/M3 responses are on the
# pseudolog flux axis (individual flux); M2 on median flux. This makes the ladder
# honest: R2 is ~0 at the individual and aggregate-predictor levels and only emerges
# once the RESPONSE is also aggregated (M2) — i.e. much of the species-level signal
# is response-aggregation, not the predictor.
#
# Sources the original generator for data/theme/palette/pseudolog only; builds ALL
# panels fresh. ORIGINAL SCRIPT UNTOUCHED; the two outputs it re-renders (figS11,
# figS14) are git-restored by the wrapper after. Tags a..o (row-major) + p,q,r bars.
# Output: revision/outputs/figS11_final.png
# ==============================================================================
suppressWarnings(suppressMessages(
  source("code/05_gene_flux_analysis/02_scale_dependent_gene_patterns.R")
))
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(grid) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

PRED <- tibble(
  key   = c("mcrA", "pmoA", "mmoX", "methanotroph", "ratio"),
  label = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  ix    = c("log_tree_mcra", "log_pmoa", "log_tree_mmox", "log_methanotroph", "log_ratio"),
  col   = c("#E74C3C", "#3498DB", "#5DADE2", "#1F618D", "#9B59B6"))
xlab_i <- list(mcrA = expression(log[10]~mcrA), pmoA = expression(log[10]~pmoA), mmoX = expression(log[10]~mmoX),
               methanotroph = expression(log[10]~(pmoA+mmoX)), ratio = expression(log[10]~ratio))
xlab_s <- list(mcrA = expression(log[10]~median~mcrA), pmoA = expression(log[10]~median~pmoA), mmoX = expression(log[10]~median~mmoX),
               methanotroph = expression(log[10]~median~(pmoA+mmoX)), ratio = expression(log[10]~median~ratio))
# species-median predictor + flux (matches the original M2 x definitions)
A <- list(
  mcrA         = analysis_mcra         %>% mutate(sx = log10(median_mcra + 1)),
  pmoA         = analysis_pmoa         %>% mutate(sx = log10(median_pmoa + 1)),
  mmoX         = analysis_mmox         %>% mutate(sx = log10(median_mmox + 1)),
  methanotroph = analysis_methanotroph %>% mutate(sx = log10(median_methanotroph + 1)),
  ratio        = analysis_ratio        %>% mutate(sx = median_log_ratio))
spx <- lapply(A, function(x) x %>% transmute(species, sx))

mk <- function(k, level, tag, ylab = "") {
  col <- PRED$col[match(k, PRED$key)]
  if (level == "ind") {
    d <- tree_level_complete %>% mutate(xv = .data[[PRED$ix[match(k, PRED$key)]]], yv = pseudolog10_individual(CH4_flux)) %>% filter(is.finite(xv), is.finite(yv))
    xlab <- xlab_i[[k]]
  } else if (level == "agg") {
    d <- tree_level_complete %>% inner_join(spx[[k]], by = "species") %>% mutate(xv = sx, yv = pseudolog10_individual(CH4_flux)) %>% filter(is.finite(xv), is.finite(yv))
    xlab <- xlab_s[[k]]
  } else {
    d <- A[[k]] %>% mutate(xv = sx, yv = median_flux) %>% filter(is.finite(xv), is.finite(yv))
    xlab <- xlab_s[[k]]
  }
  f <- lm(yv ~ xv, d); s <- coef(f)[2]; r2 <- summary(f)$r.squared; p <- summary(f)$coef[2, 4]
  g <- ggplot(d, aes(xv, yv, color = species)) +
    geom_point(size = if (level == "sp") 3.1 else 1.9, alpha = if (level == "sp") 0.85 else 0.6) +
    geom_smooth(aes(group = 1), method = "lm", formula = y~x, se = TRUE, color = col, fill = col, alpha = 0.15, linewidth = 1) +
    annotate("label", x = Inf, y = Inf, label = sprintf("slope=%.3f\nR2=%.3f  p=%.3f", s, r2, p),
             hjust = 1.04, vjust = 1.04, size = 2.55, fill = "white", alpha = 0.9, label.size = 0.15, lineheight = 0.92) +
    scale_color_manual(values = species_palette, name = "Species") +
    labs(tag = tag, x = xlab, y = ylab) + theme_pub_gene + theme(legend.position = "none")
  if (level == "sp") {
    g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  } else {
    g <- g + geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
      scale_y_continuous(breaks = pseudolog10_individual(c(0, 0.1, 1)), labels = c("0", "0.1", "1"))
  }
  list(panel = g, r2 = r2, p = p)
}

tags <- matrix(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)"), ncol = 3, byrow = TRUE)
levels3 <- c("ind", "agg", "sp"); ylabs3 <- list(ind = y_lab_individual, agg = "", sp = y_lab_species)
res <- list()
for (r in 1:5) for (cc in 1:3) {
  k <- PRED$key[r]; lv <- levels3[cc]
  res[[paste(k, lv)]] <- mk(k, lv, tags[r, cc], ylab = ylabs3[[lv]])
}

# ---- R2 bars: gene-only R2 per column, shared y-limit ------------------------
bar_df <- function(lv) data.frame(Model = PRED$label,
  R2 = sapply(PRED$key, function(k) res[[paste(k, lv)]]$r2),
  P  = sapply(PRED$key, function(k) res[[paste(k, lv)]]$p)) %>%
  mutate(Significant = P < 0.05, Model = factor(Model, levels = Model[order(R2)]))
d1 <- bar_df("ind"); d3 <- bar_df("agg"); d2 <- bar_df("sp")
y_lim <- max(c(d1$R2, d3$R2, d2$R2)) * 1.2
bar_panel <- function(dat, tag, ylab = "") ggplot(dat, aes(Model, R2, fill = Significant)) +
  geom_col(alpha = 0.85) + geom_text(aes(label = sprintf("%.3f", R2)), vjust = -0.3, size = 2.9) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#27AE60"), labels = c("NS", "p < 0.05"), name = "") +
  ylim(0, y_lim) + labs(tag = tag, x = "", y = ylab) + theme_pub_gene +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
p_b1 <- bar_panel(d1, "(p)", expression(italic(R)^2~"(gene only)")); p_b3 <- bar_panel(d3, "(q)"); p_b2 <- bar_panel(d2, "(r)")

# ---- assemble ----------------------------------------------------------------
hh <- function(txt) wrap_elements(full = textGrob(txt, gp = gpar(fontface = "bold", fontsize = 12)))
h1 <- hh("Individual\n(individual x, individual y)")
h2 <- hh("Aggregate predictor\n(species-median x, individual y)")
h3 <- hh("Species-level\n(species-median x & y)")
hts <- c(0.13, 1, 1, 1, 1, 1, 0.92)
P <- function(k, lv) res[[paste(k, lv)]]$panel
col1 <- h1 / P("mcrA","ind") / P("pmoA","ind") / P("mmoX","ind") / P("methanotroph","ind") / P("ratio","ind") / p_b1 + plot_layout(ncol = 1, heights = hts)
col2 <- h2 / P("mcrA","agg") / P("pmoA","agg") / P("mmoX","agg") / P("methanotroph","agg") / P("ratio","agg") / p_b3 + plot_layout(ncol = 1, heights = hts)
col3 <- h3 / P("mcrA","sp")  / P("pmoA","sp")  / P("mmoX","sp")  / P("methanotroph","sp")  / P("ratio","sp")  / p_b2 + plot_layout(ncol = 1, heights = hts)
combined <- (col1 | col2 | col3) / legend_panel + plot_layout(heights = c(1, 0.05))
ggsave(file.path(out, "figS11_final.png"), combined, width = 17, height = 16, dpi = 300, bg = "white")
cat("M1 gene R2:", paste(sprintf("%s=%.3f", d1$Model, d1$R2), collapse="  "), "\n")
cat("M3 gene R2:", paste(sprintf("%s=%.3f", d3$Model, d3$R2), collapse="  "), "\n")
cat("M2 gene R2:", paste(sprintf("%s=%.3f", d2$Model, d2$R2), collapse="  "), "\n")
cat("Wrote figS11_final.png\n")
