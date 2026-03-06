# ==============================================================================
# Methanotroph ddPCR Gene–Function (PICRUSt2) Relationships
# ==============================================================================
# Purpose: Visualize how ddPCR pmoA and mmoX relate to PICRUSt2-predicted
#   methanotrophic pathway abundances. Includes overall volcano plots from
#   pre-computed results plus new compartment-specific LMER analysis.
#
# Analyses:
#   A. Volcano plots: pmoA and mmoX pathway associations (existing results)
#   B. pmoA vs mmoX t-statistic comparison across pathways
#   C. Compartment-specific LMER re-analysis (Heartwood vs Sapwood)
#   D. Compartment-specific heatmap
#
# Inputs:
#   - data/processed/molecular/picrust/pathway_associations_pmoa.csv (existing)
#   - data/processed/molecular/picrust/pathway_associations_mmox.csv (existing)
#   - data/raw/picrust/path_abun_contrib.tsv (6 GB, for compartment-specific)
#   - data/raw/picrust/path_abun_unstrat_descrip.tsv (pathway descriptions)
#   - data/raw/picrust/KO_predicted.tsv (KEGG orthologs)
#   - data/raw/picrust/16S_tree_sample_table_with_meta.csv (metadata)
#   - data/raw/16s/taxonomy_table.txt (for plastid filtering)
#
# Outputs:
#   - fig_gene_vs_function.png
#   - compartment_specific_pathway_results.csv
#
# Required packages: tidyverse, patchwork, lme4, data.table, pheatmap
# ==============================================================================

library(tidyverse)
library(patchwork)
library(lme4)

# ==============================================================================
# STEP 1: Load pre-computed pathway association results
# ==============================================================================

pvals_pmoa <- read.csv("data/processed/molecular/picrust/pathway_associations_pmoa.csv",
                        stringsAsFactors = FALSE)
pvals_mmox <- read.csv("data/processed/molecular/picrust/pathway_associations_mmox.csv",
                        stringsAsFactors = FALSE)

cat("pmoA pathway results:", nrow(pvals_pmoa), "pathways\n")
cat("  Significant (FDR < 0.05):", sum(pvals_pmoa$FDR < 0.05, na.rm = TRUE), "\n")
cat("mmoX pathway results:", nrow(pvals_mmox), "pathways\n")
cat("  Significant (FDR < 0.05):", sum(pvals_mmox$FDR < 0.05, na.rm = TRUE), "\n")

# ==============================================================================
# STEP 2: Volcano plots (A & B)
# ==============================================================================

# Highlight methane/C1-related pathways
c1_pathways <- c("RUMP-PWY", "METHANOGENESIS-PWY", "METH-ACETATE-PWY",
                  "P241-PWY", "P261-PWY", "PWY-6148", "CODH-PWY",
                  "PWY-6167", "PWY-6141", "PWY-6174")

make_volcano <- function(pvals_df, gene_label) {
  pvals_df$neg_log_fdr <- -log10(pvals_df$FDR)
  pvals_df$sig <- case_when(
    pvals_df$FDR < 0.01 ~ "FDR < 0.01",
    pvals_df$FDR < 0.05 ~ "FDR < 0.05",
    TRUE ~ "NS"
  )
  pvals_df$label <- ifelse(pvals_df$FDR < 0.01, pvals_df$description, "")

  # Truncate long labels
  pvals_df$label <- ifelse(
    nchar(pvals_df$label) > 40,
    paste0(substr(pvals_df$label, 1, 37), "..."),
    pvals_df$label
  )

  ggplot(pvals_df, aes(x = t, y = neg_log_fdr, color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text(aes(label = label), size = 2.2, hjust = -0.05, vjust = 0.5,
              check_overlap = TRUE, color = "black") +
    scale_color_manual(
      values = c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a", "NS" = "grey70"),
      name = ""
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    labs(
      x = "t-statistic",
      y = expression(-log[10]*"(FDR)"),
      title = paste(gene_label, "pathway associations")
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom")
}

p_volcano_pmoa <- make_volcano(pvals_pmoa, "pmoA")
p_volcano_mmox <- make_volcano(pvals_mmox, "mmoX")

# ==============================================================================
# STEP 3: pmoA vs mmoX t-statistic comparison (panel C)
# ==============================================================================

# Merge on pathway
comparison <- merge(
  pvals_pmoa %>% select(pathway, t, FDR, description),
  pvals_mmox %>% select(pathway, t, FDR),
  by = "pathway", suffixes = c(".pmoa", ".mmox")
)

comparison$sig_class <- case_when(
  comparison$FDR.pmoa < 0.05 & comparison$FDR.mmox < 0.05 ~ "Both",
  comparison$FDR.pmoa < 0.05 ~ "pmoA only",
  comparison$FDR.mmox < 0.05 ~ "mmoX only",
  TRUE ~ "Neither"
)

comparison$label <- ifelse(
  comparison$sig_class != "Neither",
  ifelse(nchar(comparison$description) > 35,
         paste0(substr(comparison$description, 1, 32), "..."),
         comparison$description),
  ""
)

p_comparison <- ggplot(comparison, aes(x = t.pmoa, y = t.mmox, color = sig_class)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), size = 2, hjust = -0.05, vjust = 0.5,
            check_overlap = TRUE, color = "black") +
  scale_color_manual(
    values = c("Both" = "#7b2d8e", "pmoA only" = "#e8845a",
               "mmoX only" = "#3b4cc0", "Neither" = "grey80"),
    name = "Significant for"
  ) +
  labs(
    x = "t-statistic (pmoA)",
    y = "t-statistic (mmoX)",
    title = "Pathway associations: pmoA vs mmoX"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ==============================================================================
# STEP 4: Compartment-specific LMER analysis
# ==============================================================================

cat("\n=== Loading data for compartment-specific analysis ===\n")

# Load metadata
meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv", row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

# Load pathway descriptions
pathways_raw <- read.table("data/raw/picrust/path_abun_unstrat_descrip.tsv",
                            sep = "\t", header = TRUE, quote = "")
descriptions <- setNames(pathways_raw$description, pathways_raw$pathway)

# Load taxonomy for plastid filtering
taxonomy <- read.table("data/raw/16s/taxonomy_table.txt",
                        sep = "\t", header = TRUE, row.names = 1,
                        quote = "", stringsAsFactors = FALSE)
tax_split <- str_split(taxonomy$Taxon, "; ", simplify = TRUE)
rownames(tax_split) <- rownames(taxonomy)
mito_otus <- rownames(tax_split)[which(trimws(tax_split[, 5]) == "Mitochondria")]
chloro_otus <- rownames(tax_split)[which(trimws(tax_split[, 4]) == "Chloroplast")]
plastid_otus <- unique(c(mito_otus, chloro_otus))

# Base filter: wood samples, 16S >= 100
wood_base <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(X16S_per_ul >= 100)

# Split by compartment
heartwood <- wood_base %>%
  filter(core_type == "Inner", !is.na(pmoa_loose), !is.na(mmox_loose))
sapwood <- wood_base %>%
  filter(core_type == "Outer", !is.na(pmoa_loose), !is.na(mmox_loose))

cat("Heartwood samples:", nrow(heartwood), "\n")
cat("Sapwood samples:", nrow(sapwood), "\n")

# Load pathway contribution table
cat("Loading pathway contributions (this may take a few minutes)...\n")
library(data.table)
contrib <- fread("data/raw/picrust/path_abun_contrib.tsv", header = TRUE)
cat("Loaded", nrow(contrib), "rows\n")

# Remove plastid contributions
contrib <- contrib[!(contrib$taxon %in% plastid_otus), ]

# Helper: rebuild and normalize abundance matrix
rebuild_abundance <- function(contrib_dt, sample_names) {
  sub <- contrib_dt[contrib_dt$sample %in% sample_names, ]
  agg <- sub %>%
    group_by(sample, `function`) %>%
    summarize(abundance = sum(taxon_function_abun), .groups = "drop") %>%
    pivot_wider(names_from = `function`, values_from = abundance, values_fill = 0) %>%
    as.data.frame()
  rownames(agg) <- agg$sample
  agg <- agg[, -1, drop = FALSE]
  return(agg)
}

normalize_and_filter <- function(abun_mat, min_depth = 10000, min_prevalence = 20) {
  abun_mat <- abun_mat[rowSums(abun_mat) > min_depth, , drop = FALSE]
  pct <- abun_mat / rowSums(abun_mat)
  pct <- pct[, colSums(pct > 0) > min(min_prevalence, nrow(pct) * 0.1), drop = FALSE]
  return(pct)
}

# Build per-compartment matrices
all_comp_samples <- c(rownames(heartwood), rownames(sapwood))
abun_all <- rebuild_abundance(contrib, all_comp_samples)
cat("Abundance matrix:", nrow(abun_all), "samples x", ncol(abun_all), "pathways\n")

rm(contrib)
gc()

# Per-compartment LMER (no core_type covariate since within-compartment)
run_lmer_compartment <- function(pct_mat, meta_df, gene_col, gene_label) {
  shared <- intersect(rownames(pct_mat), rownames(meta_df))
  meta_df <- meta_df[shared, ]

  cat("Running LMER for", gene_label, "(", nrow(meta_df), "samples,",
      ncol(pct_mat), "pathways)...\n")

  pvals_list <- lapply(colnames(pct_mat), function(pathway) {
    meta_copy <- meta_df
    meta_copy$pathway <- scale(pct_mat[rownames(meta_copy), pathway])
    meta_copy$gene <- scale(meta_copy[[gene_col]])

    tryCatch({
      fit0 <- lmer(gene ~ log16S + (1 | seq_id), data = meta_copy)
      fit1 <- lmer(gene ~ pathway + log16S + (1 | seq_id), data = meta_copy)
      p <- anova(fit0, fit1)[2, "Pr(>Chisq)"]
      t_val <- coef(summary(fit1))["pathway", "t value"]
      data.frame(pathway = pathway, p = p, t = t_val, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(pathway = pathway, p = NA, t = NA, stringsAsFactors = FALSE)
    })
  })

  pvals <- do.call("rbind", pvals_list)
  pvals$FDR <- p.adjust(pvals$p, "BH")
  pvals$description <- descriptions[pvals$pathway]
  return(pvals)
}

# Heartwood analysis
hw_samples <- intersect(rownames(abun_all), rownames(heartwood))
pct_hw <- normalize_and_filter(abun_all[hw_samples, ])
pvals_hw_pmoa <- run_lmer_compartment(pct_hw, heartwood, "pmoa_loose", "pmoA-Heartwood")
pvals_hw_mmox <- run_lmer_compartment(pct_hw, heartwood, "mmox_loose", "mmoX-Heartwood")

# Sapwood analysis
sw_samples <- intersect(rownames(abun_all), rownames(sapwood))
pct_sw <- normalize_and_filter(abun_all[sw_samples, ])
pvals_sw_pmoa <- run_lmer_compartment(pct_sw, sapwood, "pmoa_loose", "pmoA-Sapwood")
pvals_sw_mmox <- run_lmer_compartment(pct_sw, sapwood, "mmox_loose", "mmoX-Sapwood")

rm(abun_all)
gc()

# Save compartment-specific results
comp_results <- bind_rows(
  pvals_hw_pmoa %>% mutate(gene = "pmoA", compartment = "Heartwood"),
  pvals_hw_mmox %>% mutate(gene = "mmoX", compartment = "Heartwood"),
  pvals_sw_pmoa %>% mutate(gene = "pmoA", compartment = "Sapwood"),
  pvals_sw_mmox %>% mutate(gene = "mmoX", compartment = "Sapwood")
)
write.csv(comp_results,
          "outputs/figures/new_paper/compartment_specific_pathway_results.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 5: Compartment-specific heatmap (panel D)
# ==============================================================================

# Select pathways significant in at least one gene × compartment combination
sig_pathways <- comp_results %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  pull(pathway) %>%
  unique()

cat("\nPathways significant in at least one combination:", length(sig_pathways), "\n")

if (length(sig_pathways) > 0) {
  # Build t-statistic matrix: pathways × (gene-compartment)
  heatmap_data <- comp_results %>%
    filter(pathway %in% sig_pathways) %>%
    mutate(combo = paste(gene, compartment, sep = "-")) %>%
    select(pathway, description, combo, t) %>%
    pivot_wider(names_from = combo, values_from = t)

  # Order by mean absolute t
  heatmap_data$mean_abs_t <- rowMeans(abs(heatmap_data[, 3:ncol(heatmap_data)]), na.rm = TRUE)
  heatmap_data <- heatmap_data %>% arrange(desc(mean_abs_t))

  # Truncate to top 30 if too many
  if (nrow(heatmap_data) > 30) {
    heatmap_data <- heatmap_data[1:30, ]
  }

  # Prepare for pheatmap
  hm_mat <- as.matrix(heatmap_data[, grep("^(pmoA|mmoX)", colnames(heatmap_data))])
  rownames(hm_mat) <- ifelse(
    nchar(heatmap_data$description) > 50,
    paste0(substr(heatmap_data$description, 1, 47), "..."),
    heatmap_data$description
  )
  hm_mat[is.na(hm_mat)] <- 0

  library(pheatmap)

  hm_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea",
                                    "#f7f7f7",
                                    "#f5c5a3", "#e8845a", "#b40426"))(100)

  # Column annotations
  ann_col <- data.frame(
    Gene = gsub("-.*", "", colnames(hm_mat)),
    Compartment = gsub(".*-", "", colnames(hm_mat))
  )
  rownames(ann_col) <- colnames(hm_mat)

  ann_colors <- list(
    Gene = c("pmoA" = "#e8845a", "mmoX" = "#3b4cc0"),
    Compartment = c("Heartwood" = "#a6611a", "Sapwood" = "#dfc27d")
  )

  png("outputs/figures/new_paper/fig_compartment_pathway_heatmap.png",
      width = 10, height = max(7, nrow(hm_mat) * 0.35 + 3),
      units = "in", res = 300)

  pheatmap(hm_mat,
           color = hm_colors,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           fontsize_row = 8,
           fontsize = 10,
           border_color = NA,
           main = "Compartment-specific pathway associations (t-statistic)")

  dev.off()
  cat("Saved: fig_compartment_pathway_heatmap.png\n")
} else {
  cat("No significant pathways found at FDR < 0.05 in compartment-specific analysis\n")
}

# ==============================================================================
# STEP 6: Combine volcano + comparison figure
# ==============================================================================

fig_function <- (p_volcano_pmoa + p_volcano_mmox) / p_comparison +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_gene_vs_function.png",
       fig_function, width = 14, height = 12, dpi = 300)
cat("Saved: fig_gene_vs_function.png\n")
cat("\nDone.\n")
