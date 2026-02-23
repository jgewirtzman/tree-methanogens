# ==============================================================================
# Taxonomy–mcrA Association Heatmaps (Family-level 16S)
# ==============================================================================
# Purpose: Tests which bacterial/archaeal families (16S) are associated with
#   mcrA (methanogen) abundance using LMER mixed-effects models, then
#   visualizes significant families as publication-quality heatmaps.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - OTU table:     data/raw/16s/OTU_table.txt
#   - Taxonomy:      data/raw/16s/taxonomy_table.txt
#   - Metadata:      data/raw/picrust/16S_tree_sample_table_with_meta.csv
#
# Outputs:
#   - fig_taxonomy_mcra_strict.png  (prevalence > 20, p < 0.01)
#   - fig_taxonomy_mcra_loose.png   (prevalence > 2,  p < 0.05)
#
# Required packages: tidyverse, lme4, pheatmap, grDevices
# ==============================================================================

library(tidyverse)
library(lme4)
library(pheatmap)

# ==============================================================================
# STEP 1: Load and prepare data
# ==============================================================================

## Taxonomy
taxonomy <- read.table("data/raw/16s/taxonomy_table.txt",
                        sep = "\t", header = TRUE, row.names = 1, quote = "",
                        stringsAsFactors = FALSE)

## OTU abundance table → numeric matrix, transposed to samples × OTUs
otu_table <- read.table("data/raw/16s/OTU_table.txt",
                         sep = "\t", header = TRUE, row.names = 1, quote = "",
                         check.names = FALSE)
otu_matrix <- apply(as.matrix(otu_table), 2, as.numeric)
rownames(otu_matrix) <- rownames(otu_table)
abundance <- t(otu_matrix)

## Metadata — shared with PICRUSt analysis
meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

# ==============================================================================
# STEP 2: Filter samples (consistent with all heatmap scripts)
# ==============================================================================

wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(mcra_probe_loose)) %>%
  filter(X16S_per_ul >= 100)

## Match samples present in both OTU table and metadata
samples <- intersect(rownames(abundance), rownames(wood))
wood <- wood[samples, ]
abundance <- abundance[samples, ]

cat("Wood samples with mcrA + 16S >= 100:", nrow(wood), "\n")

# ==============================================================================
# STEP 3: Parse taxonomy → aggregate to Family level
# ==============================================================================

tax_split <- str_split(taxonomy$Taxon, "; ", simplify = TRUE)
families <- tax_split[, 5]
families <- gsub("^\\s+|\\s+$", "", families)
families[families == "none" | families == ""] <- NA
names(families) <- rownames(taxonomy)

## Keep OTUs present in both abundance and taxonomy
common_otus <- intersect(colnames(abundance), names(families))
abundance <- abundance[, common_otus]
families_filtered <- families[common_otus]

## Sum OTU abundances by family
unique_families <- unique(families_filtered)
unique_families <- unique_families[!is.na(unique_families) &
                                     unique_families != "" &
                                     unique_families != "none"]

family_abundance <- matrix(0, nrow = nrow(abundance),
                            ncol = length(unique_families))
rownames(family_abundance) <- rownames(abundance)
colnames(family_abundance) <- unique_families

for (i in seq_len(ncol(abundance))) {
  fam <- families_filtered[colnames(abundance)[i]]
  if (!is.na(fam) && fam != "" && fam != "none" && fam %in% colnames(family_abundance)) {
    family_abundance[, fam] <- family_abundance[, fam] + as.numeric(abundance[, i])
  }
}

cat("Unique families:", length(unique_families), "\n")

# ==============================================================================
# STEP 4: Run LMER association tests (two prevalence thresholds)
# ==============================================================================

run_lmer_associations <- function(percent_mat, wood_df) {
  pvals_list <- lapply(colnames(percent_mat), function(family) {
    meta_copy <- wood_df
    meta_copy$family <- scale(percent_mat[, family])
    meta_copy$mcra <- scale(meta_copy$mcra_probe_loose)

    fit0 <- lmer(mcra ~ core_type + log16S + (1 | seq_id), data = meta_copy)
    fit1 <- lmer(mcra ~ family + core_type + log16S + (1 | seq_id), data = meta_copy)
    p <- anova(fit0, fit1)[2, "Pr(>Chisq)"]
    t_val <- coef(summary(fit1))["family", "t value"]
    data.frame(family = family, p = p, t = t_val, stringsAsFactors = FALSE)
  })
  pvals <- do.call("rbind", pvals_list)
  pvals$FDR <- p.adjust(pvals$p, "BH")
  return(pvals)
}

## Relative abundance
percent_all <- family_abundance / rowSums(family_abundance)

## --- Strict version: prevalence > 20 ---
percent_strict <- percent_all[, colSums(percent_all > 0) > 20]
cat("\nStrict: testing", ncol(percent_strict), "families (prevalence > 20)\n")
pvals_strict <- run_lmer_associations(percent_strict, wood)

## --- Loose version: prevalence > 2 ---
percent_loose <- percent_all[, colSums(percent_all > 0) > 2]
cat("Loose: testing", ncol(percent_loose), "families (prevalence > 2)\n")
pvals_loose <- run_lmer_associations(percent_loose, wood)

# ==============================================================================
# STEP 5: Build pheatmap function
# ==============================================================================

make_taxonomy_heatmap <- function(pvals, percent_mat, wood_df, p_cutoff,
                                   filename, fig_title) {

  sig <- pvals %>%
    filter(!is.na(p), p < p_cutoff) %>%
    arrange(t)

  cat(fig_title, ": found", nrow(sig), "significant families at p <", p_cutoff, "\n")
  if (nrow(sig) == 0) {
    cat("  No significant families — skipping heatmap.\n")
    return(invisible(NULL))
  }

  rownames(sig) <- sig$family

  ## Order samples by mcrA abundance
  ordered_samples <- rownames(wood_df[order(wood_df$mcra_probe_loose), ])
  sig_percent <- percent_mat[ordered_samples, sig$family]

  ## Column annotations (samples)
  wood_df$log10_mcra <- log10(1 + wood_df$mcra_probe_loose)
  ann_col <- wood_df[ordered_samples, c("core_type", "log16S", "log10_mcra")]
  colnames(ann_col) <- c("Type", "log10(16S)", "log10(mcrA)")

  ## Row annotations (families)
  ann_row <- data.frame(
    Significance = ifelse(sig$p < 0.01, "p < 0.01",
                    ifelse(sig$p < 0.05, "p < 0.05", "N.S."))
  )
  rownames(ann_row) <- sig$family

  ## Color schemes
  ann_colors <- list(
    Type = c(Inner = "#a6611a", Outer = "#dfc27d"),
    Significance = c("p < 0.01" = "#d7191c", "p < 0.05" = "#fdae61", "N.S." = "#abd9e9")
  )

  ## Determine height based on number of families
  fig_height <- max(6, nrow(sig) * 0.3 + 3)

  ## Save as PNG
  png(paste0("outputs/figures/", filename), width = 12, height = fig_height,
      units = "in", res = 300)
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
           labels_row = sig$family,
           annotation_names_row = FALSE,
           fontsize_row = 9,
           fontsize = 10,
           main = fig_title)
  dev.off()
  cat("  Saved:", filename, "\n")
}

# ==============================================================================
# STEP 6: Generate heatmaps
# ==============================================================================

make_taxonomy_heatmap(pvals_strict, percent_strict, wood, p_cutoff = 0.01,
                       filename = "main/fig6_taxonomy_mcra_heatmap.png",
                       fig_title = "Family-level 16S associations with mcrA (prevalence > 20, p < 0.01)")

# make_taxonomy_heatmap(pvals_loose, percent_loose, wood, p_cutoff = 0.05,
#                        filename = "fig_taxonomy_mcra_loose.png",
#                        fig_title = "Family-level 16S associations with mcrA (prevalence > 2, p < 0.05)")

## Save association tables
write.csv(arrange(pvals_strict, p),
          "outputs/tables/family_mcra_associations_strict.csv", row.names = FALSE)
write.csv(arrange(pvals_loose, p),
          "outputs/tables/family_mcra_associations_loose.csv", row.names = FALSE)

cat("\nDone. Taxonomy–mcrA heatmap figures saved.\n")
