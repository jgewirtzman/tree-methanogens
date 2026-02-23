# ==============================================================================
# PICRUSt2 MetaCyc Pathway–mcrA Association Heatmap
# ==============================================================================
# Purpose: Visualizes MetaCyc pathways significantly associated with mcrA
#   (methanogen) abundance. Uses results from the no-mcrA-OTU analysis:
#   pathway abundances reconstructed after removing methanogen contributions,
#   then tested for association with mcrA via LMER. This reveals community-
#   level functional shifts associated with methanogens beyond the methanogens
#   themselves.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - Pre-computed results: data/processed/molecular/metacyc_pathways_association_w_mcra.csv
#   - Pathway abundances:   data/raw/picrust/path_abun_unstrat_descrip.tsv
#   - Metadata:             data/raw/picrust/16S_tree_sample_table_with_meta.csv
#
# Outputs:
#   - fig_picrust_mcra_pathways.png
#
# Required packages: tidyverse, pheatmap, grDevices
# ==============================================================================

library(tidyverse)
library(pheatmap)

# ==============================================================================
# STEP 1: Load metadata and apply consistent filtering
# ==============================================================================

meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(mcra_probe_loose)) %>%
  filter(X16S_per_ul >= 100)

cat("Wood samples with mcrA + 16S >= 100:", nrow(wood), "\n")

# ==============================================================================
# STEP 2: Load pathway abundances and filter to matching samples
# ==============================================================================

pathways_raw <- read.table("data/raw/picrust/path_abun_unstrat_descrip.tsv",
                            sep = "\t", header = TRUE, quote = "")
descriptions <- pathways_raw[, 1:2]
rownames(descriptions) <- pathways_raw$pathway

abundance <- pathways_raw[, 3:ncol(pathways_raw)]
rownames(abundance) <- pathways_raw$pathway
abundance <- t(abundance)

## Filter low-depth samples
abundance <- abundance[rowSums(abundance) > 10000, ]

## Match samples
samples <- intersect(rownames(abundance), rownames(wood))
wood <- wood[samples, ]
abundance <- abundance[samples, ]

cat("Samples after matching:", nrow(wood), "\n")

## Relative abundance
percent <- abundance / rowSums(abundance)
percent <- percent[, colSums(percent > 0) > 20]

cat("Pathways passing prevalence filter (> 20 samples):", ncol(percent), "\n")

# ==============================================================================
# STEP 3: Load pre-computed association results
# ==============================================================================

pvals_combined <- read.csv("data/processed/molecular/metacyc_pathways_association_w_mcra.csv",
                            stringsAsFactors = FALSE)

cat("\nTotal pathways tested:", nrow(pvals_combined), "\n")
cat("Significant at FDR.no_mcra < 0.01:", sum(pvals_combined$FDR.no_mcra < 0.01, na.rm = TRUE), "\n")
cat("Significant at FDR.all < 0.01:", sum(pvals_combined$FDR.all < 0.01, na.rm = TRUE), "\n")

# ==============================================================================
# STEP 4: Build heatmap from no-mcrA-OTU results (FDR < 0.01)
# ==============================================================================

sig_pathways <- pvals_combined %>%
  filter(FDR.no_mcra < 0.01) %>%
  arrange(t.no_mcra)

rownames(sig_pathways) <- sig_pathways$pathway
sig_pathways$logFDR <- -log10(sig_pathways$FDR.no_mcra)

cat("\nPathways in heatmap:", nrow(sig_pathways), "\n")

## Match to available pathway abundances
available <- intersect(sig_pathways$pathway, colnames(percent))
sig_pathways <- sig_pathways[sig_pathways$pathway %in% available, ]

cat("Pathways with abundance data:", nrow(sig_pathways), "\n")

## Order samples by mcrA abundance
ordered_samples <- rownames(wood[order(wood$mcra_probe_loose), ])
sig_percent <- percent[ordered_samples, sig_pathways$pathway]

## Column annotations (samples)
wood$log10_mcra <- log10(1 + wood$mcra_probe_loose)
ann_col <- wood[ordered_samples, c("core_type", "log16S", "log10_mcra")]
colnames(ann_col) <- c("Type", "log10(16S)", "log10(mcrA)")

## Row annotations (pathways)
ann_row <- sig_pathways[, c("logFDR", "mean_counts_per_million"), drop = FALSE]
colnames(ann_row) <- c("-log10(FDR)", "Mean CPM")

## Color schemes
ann_colors <- list(
  Type = c(Inner = "#a6611a", Outer = "#dfc27d")
)

## Figure dimensions based on number of pathways
fig_height <- max(8, nrow(sig_pathways) * 0.3 + 3)

# ==============================================================================
# STEP 5: Generate and save heatmap
# ==============================================================================

png("outputs/figures/fig_picrust_mcra_pathways.png",
    width = 16, height = fig_height, units = "in", res = 300)
pheatmap(t(sig_percent),
         color = colorRampPalette(c("dodgerblue", "white", "red"))(100),
         scale = "row",
         annotation_col = ann_col,
         annotation_row = ann_row,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         treeheight_row = 0,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         labels_row = sig_pathways$description,
         annotation_names_row = FALSE,
         fontsize_row = 8,
         fontsize = 10,
         main = "MetaCyc pathways associated with mcrA (non-methanogen OTUs, FDR < 0.01)")
dev.off()

cat("\nHeatmap saved: fig_picrust_mcra_pathways.png\n")
cat("Done.\n")
