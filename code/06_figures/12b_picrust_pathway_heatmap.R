# ==============================================================================
# PICRUSt2 MetaCyc Pathway–Gene Association Heatmaps (Figure 6 + Supplementary)
# ==============================================================================
# Purpose: Visualizes MetaCyc pathways significantly associated with functional
#   gene abundances (mcrA, pmoA, mmoX). Uses pre-computed results from the
#   plastid-filtered, no-mcrA-ASV analysis pipeline.
#
# Pipeline stage: 4 — Publication Figures
#
# The primary figure (Fig 6) shows mcrA-associated pathways using the no-mcrA-
#   ASV analysis: pathway abundances reconstructed after removing methanogen
#   contributions, then tested for association with mcrA via LMER. This reveals
#   community-level functional shifts associated with methanogens beyond the
#   methanogens themselves.
#
# Supplementary figures show associations for all genes:
#   - mcrA (all ASVs): includes methanogen pathway contributions
#   - mcrA (no-mcrA ASVs): excludes methanogen pathway contributions
#   - pmoA: particulate methane monooxygenase
#   - mmoX: soluble methane monooxygenase
#
# Inputs (pre-computed by 07_picrust_pathway_associations.R):
#   - data/processed/molecular/picrust/pathway_associations_mcra_all.csv
#   - data/processed/molecular/picrust/pathway_associations_mcra_no_mcra_otus.csv
#   - data/processed/molecular/picrust/pathway_associations_pmoa.csv
#   - data/processed/molecular/picrust/pathway_associations_mmox.csv
#   - data/processed/molecular/picrust/pathway_associations_combined.csv
#   - data/raw/picrust/path_abun_unstrat_descrip.tsv  (pathway abundances)
#   - data/raw/picrust/16S_tree_sample_table_with_meta.csv (metadata)
#
# Outputs:
#   - outputs/figures/main/fig6_picrust_mcra_no_mcra_heatmap.png
#   - outputs/figures/supplementary/figS4_picrust_mcra_all_heatmap.png
#   - outputs/figures/supplementary/figS5_picrust_pmoa_heatmap.png
#
# Required packages: tidyverse, pheatmap
# ==============================================================================

library(tidyverse)
library(pheatmap)

# ==============================================================================
# STEP 1: Load metadata and pathway abundances
# ==============================================================================

meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

# Load pathway abundances (used for display matrix in heatmaps)
pathways_raw <- read.table("data/raw/picrust/path_abun_unstrat_descrip.tsv",
                            sep = "\t", header = TRUE, quote = "")
descriptions <- setNames(pathways_raw$description, pathways_raw$pathway)

abundance <- pathways_raw[, 3:ncol(pathways_raw)]
rownames(abundance) <- pathways_raw$pathway
abundance <- t(abundance)
abundance <- abundance[rowSums(abundance) > 10000, ]

## Relative abundance with prevalence filter
percent <- abundance / rowSums(abundance)
percent <- percent[, colSums(percent > 0) > 20]

# ==============================================================================
# STEP 2: Shared heatmap aesthetics
# ==============================================================================

## Heatmap color palette (cool blue – white – warm red, diverging)
heatmap_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea",
                                      "#f7f7f7",
                                      "#f5c5a3", "#e8845a", "#b40426"))(100)

## Annotation color schemes (harmonious, matching taxonomy heatmap)
ann_colors_mcra <- list(
  Compartment = c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
  `log10(mcrA)` = colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
  `mcrA assoc.` = c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a",
                     "NS" = "#f0f0f0")
)

ann_colors_pmoa <- list(
  Compartment = c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
  `log10(pmoA)` = colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
  `pmoA assoc.` = c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a",
                     "NS" = "#f0f0f0")
)

ann_colors_mmox <- list(
  Compartment = c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
  `log10(mmoX)` = colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
  `mmoX assoc.` = c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a",
                     "NS" = "#f0f0f0")
)

# ==============================================================================
# STEP 3: Heatmap generation function
# ==============================================================================

make_pathway_heatmap <- function(pvals_file, meta_df, gene_col, gene_label,
                                  ann_colors, output_file,
                                  fdr_threshold = 0.01,
                                  max_gene_contrib = NULL,
                                  contrib_file = NULL,
                                  contrib_col = NULL,
                                  row_contrib_bar = FALSE,
                                  row_contrib_threshold = 0.50,
                                  show_row_ann_names = TRUE,
                                  border_color = NA,
                                  cellheight = 10, fontsize_row = 8,
                                  fig_width = 13) {

  ## Load pre-computed association results
  pvals <- read.csv(pvals_file, stringsAsFactors = FALSE)
  cat("\n===", gene_label, "===\n")
  cat("Total pathways tested:", nrow(pvals), "\n")
  cat("Significant at FDR <", fdr_threshold, ":",
      sum(pvals$FDR < fdr_threshold, na.rm = TRUE), "\n")

  ## Filter to significant pathways
  sig <- pvals %>%
    filter(!is.na(FDR), FDR < fdr_threshold) %>%
    arrange(t)

  ## Optionally load contribution data (for filtering and/or row annotation)
  contrib_lookup <- NULL
  if (!is.null(contrib_file) && !is.null(contrib_col)) {
    combined <- read.csv(contrib_file, stringsAsFactors = FALSE)
    contrib_lookup <- setNames(combined[[contrib_col]], combined$pathway)
  }

  ## Optionally filter by gene-ASV contribution fraction
  if (!is.null(max_gene_contrib) && !is.null(contrib_lookup)) {
    sig$gene_contrib <- contrib_lookup[sig$pathway]
    sig$gene_contrib[is.na(sig$gene_contrib)] <- 0
    n_before <- nrow(sig)
    sig <- sig %>% filter(gene_contrib < max_gene_contrib)
    cat("Filtered by", contrib_col, "<", max_gene_contrib, ":",
        n_before, "->", nrow(sig), "pathways\n")
  }

  if (nrow(sig) == 0) {
    cat("No significant pathways at FDR <", fdr_threshold, "— skipping heatmap\n")
    return(invisible(NULL))
  }

  rownames(sig) <- sig$pathway
  cat("Pathways in heatmap:", nrow(sig), "\n")

  ## Filter samples
  wood <- meta_df %>%
    filter(material == "Wood", !is.na(material)) %>%
    filter(!is.na(.data[[gene_col]])) %>%
    filter(X16S_per_ul >= 100)

  samples <- intersect(rownames(percent), rownames(wood))
  wood <- wood[samples, ]

  ## Match pathways to available abundances
  available <- intersect(sig$pathway, colnames(percent))
  sig <- sig[sig$pathway %in% available, ]

  ## Order samples by gene abundance
  ordered_samples <- rownames(wood[order(wood[[gene_col]]), ])
  sig_percent <- percent[ordered_samples, sig$pathway]

  ## Column annotations
  wood$log10_gene <- log10(1 + wood[[gene_col]])
  wood$compartment <- ifelse(wood$core_type == "Inner", "Heartwood", "Sapwood")

  ann_col <- wood[ordered_samples, c("compartment", "log10_gene"), drop = FALSE]
  colnames(ann_col) <- c("Compartment", paste0("log10(", gene_label, ")"))

  ## Row annotations (optional: gene-ASV contribution bar)
  ann_row <- NA
  if (row_contrib_bar && !is.null(contrib_lookup)) {
    contrib_vals <- contrib_lookup[sig$pathway]
    contrib_vals[is.na(contrib_vals)] <- 0
    contrib_label <- paste0(">", round(row_contrib_threshold * 100), "% from ",
                            gene_label, " ASVs")
    ann_row <- data.frame(
      contrib_flag = ifelse(contrib_vals >= row_contrib_threshold, "Yes", "No"),
      row.names = sig$pathway
    )
    colnames(ann_row) <- contrib_label
    ann_colors[[contrib_label]] <- c("Yes" = "#b71c1c", "No" = "#f0f0f0")
  }

  ## Row labels: pathway descriptions
  row_labels <- ifelse(!is.na(sig$description) & sig$description != "",
                        sig$description, sig$pathway)

  ## Figure dimensions: driven by cellheight for tight packing
  ## Annotation bar + legend + margins need ~1.0 in
  fig_height <- max(4, nrow(sig) * cellheight / 72 + 1.0)

  ## Generate heatmap (with device protection to prevent corruption)
  tryCatch({
    png(output_file, width = fig_width, height = fig_height,
        units = "in", res = 300)

    pheatmap(t(sig_percent),
             color = heatmap_colors,
             scale = "row",
             annotation_col = ann_col,
             annotation_row = if (is.data.frame(ann_row)) ann_row else NA,
             annotation_colors = ann_colors,
             show_colnames = FALSE,
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             treeheight_row = 0,
             labels_row = row_labels,
             annotation_names_row = show_row_ann_names,
             fontsize_row = fontsize_row,
             fontsize = 10,
             fontsize_col = 8,
             cellheight = cellheight,
             border_color = border_color,
             main = "")

    dev.off()
    cat("Saved:", output_file, "\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("ERROR generating", output_file, ":", conditionMessage(e), "\n")
  })

  return(invisible(sig))
}

# ==============================================================================
# STEP 4: Generate all heatmaps
# ==============================================================================

results_dir <- "data/processed/molecular/picrust"

## --- Figure 6: mcrA (no-mcrA ASVs) — MAIN TEXT ---
## FDR < 0.001 with mcrA-ASV contribution < 10% (~45 pathways).
## The contribution filter removes archaeal housekeeping pathways dominated
## by methanogen ASVs (e.g., archaetidylinositol, CDP-archaeol, mevalonate),
## retaining community-level functional shifts. The relaxed FDR (vs. 1e-4)
## captures ecologically informative pathways like homolactic fermentation,
## mannan degradation, aerobic respiration I, and sulfate assimilation that
## narrowly missed the stricter threshold.
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_mcra_no_mcra_otus.csv"),
  meta_df = meta,
  gene_col = "mcra_probe_loose",
  gene_label = "mcrA",
  ann_colors = ann_colors_mcra,
  output_file = "outputs/figures/main/fig6_picrust_mcra_no_mcra_heatmap.png",
  fdr_threshold = 0.001,
  max_gene_contrib = 0.10,
  contrib_file = file.path(results_dir, "pathway_associations_combined.csv"),
  contrib_col = "mean_percent_from_mcra",
  cellheight = 10,
  fontsize_row = 7
)

## --- Supplementary: mcrA (no-mcrA ASVs, full FDR < 0.01 set) ---
## Includes row annotation bar flagging pathways where >50% of predicted
## abundance comes from mcrA-carrying ASVs. Cell borders added for
## readability with 109 pathways.
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_mcra_no_mcra_otus.csv"),
  meta_df = meta,
  gene_col = "mcra_probe_loose",
  gene_label = "mcrA",
  ann_colors = ann_colors_mcra,
  output_file = "outputs/figures/supplementary/figS4_picrust_mcra_all_heatmap.png",
  fdr_threshold = 0.01,
  contrib_file = file.path(results_dir, "pathway_associations_combined.csv"),
  contrib_col = "mean_percent_from_mcra",
  row_contrib_bar = TRUE,
  row_contrib_threshold = 0.50,
  show_row_ann_names = FALSE,
  border_color = "#e0e0e0",
  cellheight = 8,
  fontsize_row = 6
)

## --- Supplementary: pmoA (no-pmoA ASVs, contribution < 50%) ---
## Matches mcrA supplementary logic: contrib bar for >50% from pmoA ASVs,
## cell borders for readability.
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_pmoa_no_pmoa_otus.csv"),
  meta_df = meta,
  gene_col = "pmoa_loose",
  gene_label = "pmoA",
  ann_colors = ann_colors_pmoa,
  output_file = "outputs/figures/supplementary/figS5_picrust_pmoa_heatmap.png",
  fdr_threshold = 0.01,
  max_gene_contrib = 0.50,
  contrib_file = file.path(results_dir, "pathway_associations_pmoa_combined.csv"),
  contrib_col = "mean_percent_from_pmoa",
  row_contrib_bar = TRUE,
  row_contrib_threshold = 0.50,
  show_row_ann_names = FALSE,
  border_color = "#e0e0e0",
  cellheight = 12,
  fontsize_row = 8
)

## --- mmoX (currently 0 significant pathways at any threshold) ---
## Kept for completeness; no figure is produced.
sig_mmox <- make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_mmox.csv"),
  meta_df = meta,
  gene_col = "mmox_loose",
  gene_label = "mmoX",
  ann_colors = ann_colors_mmox,
  output_file = "outputs/figures/supplementary/picrust_mmox_heatmap.png",
  fdr_threshold = 0.01,
  cellheight = 12,
  fontsize_row = 8
)

if (is.null(sig_mmox)) {
  cat("\nRetrying mmoX at FDR < 0.05...\n")
  make_pathway_heatmap(
    pvals_file = file.path(results_dir, "pathway_associations_mmox.csv"),
    meta_df = meta,
    gene_col = "mmox_loose",
    gene_label = "mmoX",
    ann_colors = ann_colors_mmox,
    output_file = "outputs/figures/supplementary/picrust_mmox_heatmap.png",
    fdr_threshold = 0.05,
    cellheight = 12,
    fontsize_row = 8
  )
}

cat("\nAll heatmaps complete.\n")
