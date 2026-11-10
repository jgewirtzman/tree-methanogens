# ==============================================================================
# fig_si_heatmap.R -- feature x sample heatmap for the SI. SOURCED after fig_style.R.
# ------------------------------------------------------------------------------
# Used by figS14_taxonomy-mcra.R and figS09_taxonomy-pmoa.R (16S families, Figs S9-S10)
# and figS11-13_picrust-heatmaps.R (MetaCyc pathways, Figs S12-S14). It draws what
# pheatmap(scale = "row") drew before, in ggplot2 so the SI house style applies:
#   cells   each feature (row) scaled across samples to mean 0, SD 1 (pheatmap's
#           scale = "row"); diverging colours symmetric about 0, limits +/- max |z|
#   top     log10 of the gene and the compartment of each sample
#   left    an optional class per feature (association FDR, or > 50% from gene ASVs)
# Columns and rows are drawn in the order given; nothing is clustered.
# ==============================================================================
suppressPackageStartupMessages(library(ggnewscale))

HM_COLS   <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea", "#f7f7f7",
                                "#f5c5a3", "#e8845a", "#b40426"))(100)   # row z-scores
GENE_RAMP <- colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100)   # log10 gene

# mat         features x samples (raw relative abundance), rows and columns in display order
# row_labels  character, one per row of mat (drawn top to bottom)
# gene        "mcrA" / "pmoA" (italic in labels)
# log_gene    log10(1 + gene), one per column;  compartment  "Heartwood"/"Sapwood", per column
# row_class   optional character, one per row; row_class_cols named colours (legend order);
#             row_class_title  legend title (expression allowed); row_class_labels  legend labels
# row_in      inches per feature row at the intended figure height (sizes the top bars)
# row_text    font size (pt) of the row labels
si_feature_heatmap <- function(mat, row_labels, gene, log_gene, compartment,
                               row_class = NULL, row_class_cols = NULL,
                               row_class_title = NULL, row_class_labels = waiver(),
                               row_in = 0.2, row_text = 8.5) {
  n <- nrow(mat); m <- ncol(mat)
  cbar <- function(o) guide_colourbar(order = o, theme = theme(legend.key.width = unit(7, "lines"),
                                                                legend.key.height = unit(0.7, "lines")))
  z <- t(scale(t(mat)))                       # (x - row mean) / row SD, as pheatmap scale = "row"
  zmax <- max(abs(range(z, na.rm = TRUE)))

  cells <- data.frame(x = rep(seq_len(m), each = n), y = rep(n:1, times = m), z = as.vector(z))

  # top bars, in row units: compartment directly above the cells, log10 gene above it
  bar_h <- 0.15 / row_in; gap <- 0.05 / row_in
  y_comp <- n + 0.5 + gap + bar_h / 2
  y_gene <- y_comp + bar_h + gap
  top <- data.frame(x = seq_len(m), comp = factor(compartment, levels = c("Heartwood", "Sapwood")),
                    lg = log_gene)

  p <- ggplot() +
    geom_raster(data = cells, aes(x, y, fill = z)) +
    scale_fill_gradientn(colours = HM_COLS, limits = c(-zmax, zmax), na.value = "grey90",
                         name = "Row z-score", guide = cbar(1)) +
    new_scale_fill() +
    geom_tile(data = top, aes(x, y_comp, fill = comp), height = bar_h, width = 1) +
    scale_fill_manual(values = COMP_COLS[c("Heartwood", "Sapwood")], name = "Compartment",
                      guide = guide_legend(order = 2, ncol = 1)) +
    new_scale_fill() +
    geom_tile(data = top, aes(x, y_gene, fill = lg), height = bar_h, width = 1) +
    scale_fill_gradientn(colours = GENE_RAMP, name = lab_log10(gene),
                         guide = cbar(3))

  x_left <- 0.5
  if (!is.null(row_class)) {
    bar_w <- 0.022 * m; bar_gap <- 0.008 * m
    x_left <- 0.5 - bar_gap - bar_w
    left <- data.frame(y = n:1, cls = row_class, w = bar_w)
    absent <- setdiff(names(row_class_cols), row_class)   # zero-width tiles: legend keys for
    if (length(absent)) left <- rbind(left, data.frame(y = 1, cls = absent, w = 0))   # absent classes
    left$cls <- factor(left$cls, levels = names(row_class_cols))
    p <- p + new_scale_fill() +
      geom_tile(data = left, aes(x = 0.5 - bar_gap - bar_w / 2, y = y, fill = cls, width = w),
                height = 1) +
      scale_fill_manual(values = row_class_cols, labels = row_class_labels, drop = FALSE,
                        name = row_class_title, guide = guide_legend(order = 4, ncol = 1))
  }

  labs_y <- do.call(expression, c(as.list(row_labels), list("Compartment", lab_log10(gene))))
  p +
    scale_x_continuous(expand = c(0, 0), limits = c(x_left, m + 0.5)) +
    scale_y_continuous(breaks = c(n:1, y_comp, y_gene), labels = labs_y, position = "right",
                       expand = c(0, 0), limits = c(0.5, y_gene + bar_h / 2)) +
    labs(x = NULL, y = NULL) +
    theme_si() +
    theme(panel.border = element_blank(), panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y.right = element_text(size = row_text, colour = "black",
                                           margin = margin(l = 2)),
          legend.position = "bottom", legend.box = "horizontal",
          legend.justification = "left", legend.title.position = "top",
          legend.title = element_text(size = SI_BASE * 0.9), legend.text = element_text(size = SI_BASE * 0.8),
          legend.spacing.x = unit(16, "pt"),
          legend.key.size = unit(0.8, "lines"),
          plot.margin = margin(4, 4, 2, 4))
}
