# ==============================================================================
# figS11-13_picrust-heatmaps.R -- Figures S13, S14, S15: MetaCyc pathways associated
# with mcrA and pmoA
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Row-scaled relative abundance of PICRUSt2 MetaCyc pathways associated with the gene
# (LMER, code/07_molecular/02_picrust_pathway_associations.R) across 2021 wood samples
# ordered by the gene; pathways ordered by t. Pathway abundances are reconstructed
# without the gene-carrying ASVs ("no_mcra_otus", "no_pmoa_otus").
#   Figure S13  mcrA, FDR < 0.001, pathways with < 10% of their abundance from
#               mcrA-carrying ASVs (community-level shifts, not methanogen housekeeping)
#   Figure S14  mcrA, FDR < 0.01, all; row bar: > 50% from mcrA-carrying ASVs
#   Figure S15  pmoA, FDR < 0.01, all; row bar: > 50% from pmoA-carrying ASVs
# Inputs
#   outputs/data/picrust/pathway_associations_mcra_no_mcra_otus.csv,
#     pathway_associations_combined.csv (mean_percent_from_mcra),
#     pathway_associations_pmoa_no_pmoa_otus.csv,
#     pathway_associations_pmoa_combined.csv (mean_percent_from_pmoa)
#   data/tables/picrust_pathways.csv  pathway abundances of the non-organelle ASVs
#   survey_16s_meta() (code/lib/tables_read.R): ddPCR and 16S qPCR per sample
# Outputs
#   outputs/figures/original/main/fig6_picrust_mcra_no_mcra_heatmap.png          (S13)
#   outputs/figures/original/supplementary/figS4_picrust_mcra_all_heatmap.png    (S14)
#   outputs/figures/original/supplementary/figS5_picrust_pmoa_heatmap.png        (S15)
# ==============================================================================

library(tidyverse)
library(pheatmap)

# ==============================================================================
# STEP 1: Load metadata and pathway abundances
# ==============================================================================

source("code/lib/tables_read.R")
meta <- survey_16s_meta()   # per sequenced 2021 sample (data/tables/: samples, samples_16S, ddpcr, qpcr, trees)
meta$log16S <- log10(1 + meta$X16S_per_ul)

# Pathway abundances for the display matrix: data/tables/picrust_pathways.csv, the
# summed PICRUSt2 contributions of the non-organelle ASVs per library (the abundances
# the associations are tested on). Libraries with > 10,000; pathways in > 20 libraries.
HTML_ENT <- c("&alpha;" = "\u03b1", "&beta;" = "\u03b2", "&gamma;" = "\u03b3", "&delta;" = "\u03b4", "&omega;" = "\u03c9", "&amp;" = "&")
abundance <- read_table("picrust_pathways") %>%
  dplyr::select(library_id, pathway, abundance) %>%
  pivot_wider(names_from = pathway, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("library_id") %>% as.matrix()
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

## Annotation colours (as the family heatmaps, Figures S10-S11)
ann_colors_for <- function(gene_label) setNames(list(
  c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
  colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
  c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a", "NS" = "#f0f0f0")),
  c("Compartment", paste0("log10(", gene_label, ")"), paste0(gene_label, " assoc.")))

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
    stopifnot(!anyNA(sig$gene_contrib))   # every tested pathway has a contribution
    n_before <- nrow(sig)
    sig <- sig %>% filter(gene_contrib < max_gene_contrib)
    cat("Filtered by", contrib_col, "<", max_gene_contrib, ":",
        n_before, "->", nrow(sig), "pathways\n")
  }

  if (nrow(sig) == 0) stop("no pathways at FDR < ", fdr_threshold, " for ", gene_label, call. = FALSE)

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
  # MetaCyc descriptions carry HTML entities (e.g. "&beta;-alanine"); render the characters
  for (ent in names(HTML_ENT)) row_labels <- gsub(ent, HTML_ENT[[ent]], row_labels, fixed = TRUE)

  ## Figure dimensions: driven by cellheight for tight packing
  ## Annotation bar + legend + margins need ~1.0 in
  fig_height <- max(4, nrow(sig) * cellheight / 72 + 1.0)

  ## Generate the heatmap. On error, close every device this call may have left open and
  ## delete the half-written file, then stop: a later script must not draw into an open
  ## device, and the assembler must not copy a stale PNG.
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
    while (dev.cur() != 1L) try(dev.off(), silent = TRUE)
    if (file.exists(output_file)) unlink(output_file)
    stop("failed generating ", output_file, ": ", conditionMessage(e), call. = FALSE)
  })

  return(invisible(sig))
}

# ==============================================================================
# STEP 4: Generate all heatmaps
# ==============================================================================

results_dir <- "outputs/data/picrust"

## --- Figure S13: mcrA (no-mcrA ASVs), FDR < 0.001, mcrA-ASV contribution < 10% ---
## The contribution filter removes archaeal housekeeping pathways dominated by
## methanogen ASVs (e.g. archaetidylinositol, CDP-archaeol, mevalonate).
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_mcra_no_mcra_otus.csv"),
  meta_df = meta,
  gene_col = "mcra_probe_loose",
  gene_label = "mcrA",
  ann_colors = ann_colors_for("mcrA"),
  output_file = "outputs/figures/original/main/fig6_picrust_mcra_no_mcra_heatmap.png",
  fdr_threshold = 0.001,
  max_gene_contrib = 0.10,
  contrib_file = file.path(results_dir, "pathway_associations_combined.csv"),
  contrib_col = "mean_percent_from_mcra",
  cellheight = 10,
  fontsize_row = 7
)

## --- Figure S14: mcrA (no-mcrA ASVs), FDR < 0.01 ---
## Row bar: pathways with > 50% of predicted abundance from mcrA-carrying ASVs.
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_mcra_no_mcra_otus.csv"),
  meta_df = meta,
  gene_col = "mcra_probe_loose",
  gene_label = "mcrA",
  ann_colors = ann_colors_for("mcrA"),
  output_file = "outputs/figures/original/supplementary/figS4_picrust_mcra_all_heatmap.png",
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

## --- Figure S15: pmoA (no-pmoA ASVs), FDR < 0.01 ---
## Row bar: pathways with > 50% of predicted abundance from pmoA-carrying ASVs.
make_pathway_heatmap(
  pvals_file = file.path(results_dir, "pathway_associations_pmoa_no_pmoa_otus.csv"),
  meta_df = meta,
  gene_col = "pmoa_loose",
  gene_label = "pmoA",
  ann_colors = ann_colors_for("pmoA"),
  output_file = "outputs/figures/original/supplementary/figS5_picrust_pmoa_heatmap.png",
  fdr_threshold = 0.01,
  contrib_file = file.path(results_dir, "pathway_associations_pmoa_combined.csv"),
  contrib_col = "mean_percent_from_pmoa",
  row_contrib_bar = TRUE,
  row_contrib_threshold = 0.50,
  show_row_ann_names = FALSE,
  border_color = "#e0e0e0",
  cellheight = 12,
  fontsize_row = 8
)

cat("\nAll heatmaps complete.\n")
