# ==============================================================================
# pmoA/mmoX Ratio Analysis (Figure S12)
# ==============================================================================
# Purpose: Analyzes pMMO:sMMO balance (pmoA:mmoX ratio) and produces a
#   2-panel figure: (A) pmoA vs mmoX abundance, (B) ratio vs total abundance.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - Metadata: data/raw/picrust/16S_tree_sample_table_with_meta.csv
#
# Outputs:
#   - figS12_methanotroph_abundance_patterns.pdf
#
# Required packages: ggplot2, cowplot, gridExtra
# ==============================================================================

library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)

# ==============================================================================
# STEP 0: Load and filter data
# ==============================================================================

meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(mcra_probe_loose)) %>%
  filter(X16S_per_ul >= 100)

cat("Wood samples with mcrA + 16S >= 100:", nrow(wood), "\n")

## For pmoA heatmap (legacy, kept for reference)
annotation_colors_pmoa = list(
  core_type = c("heartwood" = "#8B4513", "sapwood" = "#DEB887"),
  log16S = colorRampPalette(c("white", "forestgreen"))(100),
  log10_pmoa_qpcr = colorRampPalette(c("white", "navy"))(100),  # white to dark blue
  `-logFDR` = colorRampPalette(c("gray90", "black"))(100),
  mean_counts = colorRampPalette(c("white", "orange"))(100)
)

# pdf("outputs/figures/pheatmap_picrust_pmoa.pdf", height = 8, width = 16)
# pheatmap(t(sig_pathways_percent_pmoa),
#          color = colorRampPalette(c("dodgerblue", "white", "red"))(100),
#          scale = "row",
#          annotation_col = annotation_col_pmoa,
#          annotation_row = annotation_row_pmoa,
#          annotation_colors = annotation_colors_pmoa,
#          show_colnames = F,
#          treeheight_row = 0,
#          cluster_cols = F,
#          cluster_rows = F,
#          labels_row = sig_pathways_pmoa$description,
#          annotation_names_row = F)
# dev.off()




library(ggplot2)

## Prepare data
wood$core_type_label = factor(wood$core_type, 
                              levels = c("Inner", "Outer"),
                              labels = c("Heartwood", "Sapwood"))

## Simple scatter plot
ggplot(wood, aes(x = log10(1 + pmoa_loose), 
                 y = log10(1 + mmox_loose),
                 color = core_type_label)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "log10(pmoA + 1)", 
       y = "log10(mmoX + 1)",
       color = "Core Type") +
  theme_bw() +
  theme(legend.position = "top")

# ggsave("outputs/figures/pmoa_vs_mmox.pdf", width = 6, height = 5)


## Only samples where both are detected
wood_both = wood[wood$pmoa_loose > 0 & wood$mmox_loose > 0,]

## Scatter plot with abundance, colored by environmental factors
ggplot(wood_both, aes(x = log10(pmoa_loose), 
                      y = log10(mmox_loose),
                      color = core_type_label)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  labs(x = "log10(pmoA)", 
       y = "log10(mmoX)",
       color = "Core Type",
       title = "pmoA vs mmoX abundance (both detected)") +
  theme_bw()
# ggsave("outputs/figures/pmoa_vs_mmox_abundance.pdf", width = 6, height = 5)

## Check correlation
cat("\n=== Correlation between pmoA and mmoX abundance ===\n")
cat("Pearson correlation:", cor(log10(wood_both$pmoa_loose), log10(wood_both$mmox_loose)), "\n")

## Look at ratio across environmental gradients
wood_both$ratio = log10(wood_both$pmoa_loose / wood_both$mmox_loose)

## Plot ratio vs environmental factors
if("mean_orp" %in% names(wood)) {
  ggplot(wood_both, aes(x = mean_orp, y = ratio)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    labs(x = "Redox potential (mV)", 
         y = "log10(pmoA / mmoX)",
         title = "Does redox affect pmoA:mmoX ratio?") +
    theme_bw()
  # ggsave("outputs/figures/ratio_vs_orp.pdf", width = 6, height = 4)
}

if("CH4_int" %in% names(wood)) {
  ggplot(wood_both, aes(x = log10(1 + CH4_int), y = ratio)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    labs(x = "log10(CH4 + 1)", 
         y = "log10(pmoA / mmoX)",
         title = "Does CH4 affect pmoA:mmoX ratio?") +
    theme_bw()
  # ggsave("outputs/figures/ratio_vs_ch4.pdf", width = 6, height = 4)
}

## Summary stats on the ratio
cat("\n=== pmoA:mmoX ratio distribution ===\n")
cat("Mean log10(pmoA/mmoX):", mean(wood_both$ratio, na.rm = TRUE), "\n")
cat("Median log10(pmoA/mmoX):", median(wood_both$ratio, na.rm = TRUE), "\n")
cat("SD log10(pmoA/mmoX):", sd(wood_both$ratio, na.rm = TRUE), "\n")

## Test if ratio differs by core type
if(sum(!is.na(wood_both$ratio)) > 0) {
  t_test = t.test(ratio ~ core_type, data = wood_both)
  cat("\n=== pmoA:mmoX ratio by core type ===\n")
  print(t_test)
}


## Does the ratio correlate with total methanotroph abundance?
ggplot(wood_both, aes(x = log10(pmoa_loose + mmox_loose), y = ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  labs(x = "log10(total pmoA + mmoX)", 
       y = "log10(pmoA / mmoX)",
       title = "Does total methanotroph abundance affect enzyme ratio?") +
  theme_bw()
# ggsave("outputs/figures/ratio_vs_total_methanotroph.pdf", width = 6, height = 4)

## Correlation with methanogen abundance
cor.test(wood_both$ratio, log10(1 + wood_both$mcra_probe_loose))

## Correlation between total methanotroph abundance and the pmoA:mmoX ratio
cat("\n=== Total methanotroph abundance vs pmoA:mmoX ratio ===\n")
cor.test(log10(wood_both$pmoa_loose + wood_both$mmox_loose), wood_both$ratio)










library(ggplot2)
library(cowplot)

## Prepare data
wood$core_type_label = factor(wood$core_type, 
                              levels = c("Inner", "Outer"),
                              labels = c("Heartwood", "Sapwood"))

wood_both = wood[wood$pmoa_loose > 0 & wood$mmox_loose > 0,]
wood_both$ratio = log10(wood_both$pmoa_loose / wood_both$mmox_loose)
wood_both$total_methanotroph = log10(wood_both$pmoa_loose + wood_both$mmox_loose)

## Calculate stats for Panel A
pmoa_log = log10(1 + wood$pmoa_loose)
mmox_log = log10(1 + wood$mmox_loose)
fit1 = lm(mmox_log ~ pmoa_log)
r2_1 = summary(fit1)$r.squared
slope_1 = coef(fit1)[2]
p_1 = summary(fit1)$coefficients[2,4]

## Calculate stats for Panel B
fit2 = lm(ratio ~ total_methanotroph, data = wood_both)
r2_2 = summary(fit2)$r.squared
slope_2 = coef(fit2)[2]
p_2 = summary(fit2)$coefficients[2,4]

## Format p-values
p1_text = ifelse(p_1 < 0.001, "p < 0.001", sprintf("p = %.3f", p_1))
p2_text = ifelse(p_2 < 0.001, "p < 0.001", sprintf("p = %.3f", p_2))

## Determine axis limits for equal scaling in Panel A
axis_max = max(c(pmoa_log, mmox_log), na.rm = TRUE)
axis_min = min(c(pmoa_log, mmox_log), na.rm = TRUE)

## Panel A: pmoA vs mmoX
p1 = ggplot(wood, aes(x = log10(1 + pmoa_loose), 
                      y = log10(1 + mmox_loose),
                      color = core_type_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  coord_fixed(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
  annotate("text", x = axis_min + 0.2, y = axis_max - 0.5,
           label = sprintf("atop(atop(slope == %.2f, italic(R)^2 == %.3f), '%s')",
                           slope_1, r2_1, p1_text),
           parse = TRUE, hjust = 0, size = 3.5) +
  labs(x = "log10(pmoA + 1)",
       y = "log10(mmoX + 1)",
       color = "Core Type") +
  theme_bw() +
  theme(legend.position = "none")

## Panel B: Ratio vs total abundance
## Get the actual range for proper placement
x_range = range(wood_both$total_methanotroph, na.rm = TRUE)
y_range = range(wood_both$ratio, na.rm = TRUE)

p2 = ggplot(wood_both, aes(x = total_methanotroph, 
                           y = ratio,
                           color = core_type_label)) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  annotate("text", x = x_range[1] + 0.1,
           y = y_range[2] - 0.3,
           label = sprintf("atop(atop(slope == %.2f, italic(R)^2 == %.3f), '%s')",
                           slope_2, r2_2, p2_text),
           parse = TRUE, hjust = 0, size = 3.5) +
  labs(x = "log10(total pmoA + mmoX)",
       y = "log10(pmoA / mmoX)",
       color = "Core Type") +
  theme_bw() +
  theme(legend.position = "none")

## Extract legend
p_legend = ggplot(wood, aes(x = log10(1 + pmoa_loose), 
                            y = log10(1 + mmox_loose),
                            color = core_type_label)) +
  geom_point() +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

legend = get_legend(p_legend)

## Combine with shared legend
combined = plot_grid(p1, p2, nrow = 1, labels = c("A", "B"), align = "h")

## Extract the legend
legend = get_legend(
  ggplot(wood, aes(x = log10(1 + pmoa_loose), 
                   y = log10(1 + mmox_loose),
                   color = core_type_label)) +
    geom_point() +
    scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
    labs(color = "Core Type") +
    theme(legend.position = "bottom")
)

## Combine: plots on top, legend on bottom
plots = plot_grid(p1, p2, nrow = 1, labels = c("A", "B"))
combined = plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.1))
combined






library(ggplot2)
library(gridExtra)

## Filter out NA core types FIRST
wood = wood[!is.na(wood$core_type),]

## Then prepare data
wood$core_type_label = factor(wood$core_type, 
                              levels = c("Inner", "Outer"),
                              labels = c("Heartwood", "Sapwood"))

## Now create wood_both from the already-filtered wood
wood_both = wood[wood$pmoa_loose > 0 & wood$mmox_loose > 0,]
wood_both$ratio = log10(wood_both$pmoa_loose / wood_both$mmox_loose)
wood_both$total_methanotroph = log10(wood_both$pmoa_loose + wood_both$mmox_loose)

## Double check - remove any remaining NAs in wood_both just to be safe
wood_both = wood_both[!is.na(wood_both$core_type_label),]

## Calculate stats for Panel A
pmoa_log = log10(1 + wood$pmoa_loose)
mmox_log = log10(1 + wood$mmox_loose)
fit1 = lm(mmox_log ~ pmoa_log)
r2_1 = summary(fit1)$r.squared
slope_1 = coef(fit1)[2]
p_1 = summary(fit1)$coefficients[2,4]

## Calculate stats for Panel B
fit2 = lm(ratio ~ total_methanotroph, data = wood_both)
r2_2 = summary(fit2)$r.squared
slope_2 = coef(fit2)[2]
p_2 = summary(fit2)$coefficients[2,4]

## Format p-values
p1_text = ifelse(p_1 < 0.001, "p < 0.001", sprintf("p = %.3f", p_1))
p2_text = ifelse(p_2 < 0.001, "p < 0.001", sprintf("p = %.3f", p_2))

## Axis limits for Panel A
axis_max = max(c(pmoa_log, mmox_log), na.rm = TRUE)
axis_min = min(c(pmoa_log, mmox_log), na.rm = TRUE)

## Panel A
p1 = ggplot(wood, aes(x = log10(1 + pmoa_loose), 
                      y = log10(1 + mmox_loose),
                      color = core_type_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  coord_fixed(xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
  annotate("text", x = axis_min + 0.2, y = axis_max - 0.5,
           label = sprintf("atop(atop(slope == %.2f, italic(R)^2 == %.3f), '%s')",
                           slope_1, r2_1, p1_text),
           parse = TRUE, hjust = 0, size = 3.5) +
  labs(x = "log10(pmoA + 1)", y = "log10(mmoX + 1)", color = "Core Type") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

## Panel B
x_range = range(wood_both$total_methanotroph, na.rm = TRUE)
y_range = range(wood_both$ratio, na.rm = TRUE)

p2 = ggplot(wood_both, aes(x = total_methanotroph, y = ratio, color = core_type_label)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Heartwood" = "#8B4513", "Sapwood" = "#DEB887")) +
  annotate("text", x = x_range[1] + 0.1, y = y_range[2] - 0.3,
           label = sprintf("atop(atop(slope == %.2f, italic(R)^2 == %.3f), '%s')",
                           slope_2, r2_2, p2_text),
           parse = TRUE, hjust = 0, size = 3.5) +
  labs(x = "log10(total pmoA + mmoX)", y = "log10(pmoA / mmoX)", color = "Core Type") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

## Combine: plots on top, legend on bottom
plots = plot_grid(p1, p2, nrow = 1)
combined = plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.1))
combined


## Add panel labels
library(patchwork)
p1_labeled <- p1 + labs(tag = "(a)") + theme(plot.tag = element_text(size = 11, face = "bold"))
p2_labeled <- p2 + labs(tag = "(b)") + theme(plot.tag = element_text(size = 11, face = "bold"))

## Save
ggsave("outputs/figures/supplementary/figS12_methanotroph_abundance_patterns.pdf",
       arrangeGrob(p1_labeled, p2_labeled, ncol = 2),
       width = 11, height = 5)
