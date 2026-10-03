# ==============================================================================
# PICRUSt2 MetaCyc Pathway–Gene Abundance Associations
# ==============================================================================
# Purpose: Tests associations between MetaCyc pathway abundances and functional
#   gene abundances (mcrA, pmoA, mmoX) using LMER mixed-effects models.
#   Filters plastid OTUs (mitochondria, chloroplasts) from pathway predictions
#   and optionally removes mcrA- or pmoA-carrying OTU contributions.
#
# Pipeline stage: 3 — Analysis
#
# Approach:
#   1. Load pathway contribution table (per-OTU contributions to each pathway)
#   2. Identify and remove plastid OTUs (mitochondria + chloroplasts)
#   3. Identify mcrA-carrying OTUs via KO_predicted.tsv (K00399) and
#      pmoA-carrying OTUs (K10944)
#   4. Reconstruct pathway abundances: (a) all non-plastid OTUs, (b) excluding
#      mcrA OTUs, (c) excluding pmoA OTUs
#   5. Run LMER for each pathway against mcrA (both versions), pmoA, and mmoX
#   6. Save association results and pre-computed abundance matrices
#
# Filtering (consistent with all 16S/PICRUSt scripts):
#   - Remove mitochondria (Family == "Mitochondria") and chloroplasts
#     (Order == "Chloroplast") from pathway contributions
#   - Wood samples only, non-NA gene abundance, 16S >= 100
#   - Pathway depth > 10,000
#   - Prevalence filter: pathways present in > 20 samples
#
# Inputs:
#   - data/raw/picrust/path_abun_contrib.tsv       (per-OTU pathway contributions)
#   - data/raw/picrust/path_abun_unstrat_descrip.tsv (pathway descriptions)
#   - data/raw/picrust/KO_predicted.tsv              (KEGG orthologs per OTU)
#   - data/raw/picrust/16S_tree_sample_table_with_meta.csv (metadata)
#   - data/raw/16s/taxonomy_table.txt                (OTU taxonomy for plastid ID)
#
# Outputs (data/processed/molecular/picrust/):
#   - pathway_associations_mcra_all.csv
#   - pathway_associations_mcra_no_mcra_otus.csv
#   - pathway_associations_pmoa.csv
#   - pathway_associations_pmoa_no_pmoa_otus.csv
#   - pathway_associations_mmox.csv
#   - pathway_associations_combined.csv       (merged mcrA all + no-mcrA results)
#   - pathway_associations_pmoa_combined.csv  (merged pmoA all + no-pmoA results)
#
# Required packages: tidyverse, lme4, data.table
# ==============================================================================

library(tidyverse)
library(lme4)
library(data.table)   # for fast reading of the 6 GB contrib file

# ==============================================================================
# STEP 1: Load metadata and identify sample sets for each gene
# ==============================================================================

meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

# Base wood filter (material + 16S threshold)
wood_base <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(X16S_per_ul >= 100)

# Gene-specific sample sets
wood_mcra <- wood_base %>% filter(!is.na(mcra_probe_loose))
wood_pmoa <- wood_base %>% filter(!is.na(pmoa_loose))
wood_mmox <- wood_base %>% filter(!is.na(mmox_loose))

cat("=== Sample sizes ===\n")
cat("Wood samples (16S >= 100):", nrow(wood_base), "\n")
cat("  with mcrA:", nrow(wood_mcra), "\n")
cat("  with pmoA:", nrow(wood_pmoa), "\n")
cat("  with mmoX:", nrow(wood_mmox), "\n")

# ==============================================================================
# STEP 2: Identify plastid OTUs from taxonomy
# ==============================================================================

taxonomy <- read.table("data/raw/16s/taxonomy_table.txt",
                        sep = "\t", header = TRUE, row.names = 1,
                        quote = "", stringsAsFactors = FALSE)

tax_split <- str_split(taxonomy$Taxon, "; ", simplify = TRUE)
rownames(tax_split) <- rownames(taxonomy)

# Mitochondria at Family level (column 5), Chloroplast at Order level (column 4)
mito_otus <- rownames(tax_split)[which(trimws(tax_split[, 5]) == "Mitochondria")]
chloro_otus <- rownames(tax_split)[which(trimws(tax_split[, 4]) == "Chloroplast")]
plastid_otus <- unique(c(mito_otus, chloro_otus))

cat("\n=== Plastid OTUs ===\n")
cat("Mitochondria:", length(mito_otus), "\n")
cat("Chloroplast:", length(chloro_otus), "\n")
cat("Total plastid OTUs:", length(plastid_otus), "\n")

# ==============================================================================
# STEP 3: Identify mcrA- and pmoA-carrying OTUs from KO predictions
# ==============================================================================

KO <- fread("data/raw/picrust/KO_predicted.tsv", header = TRUE)
KO_df <- as.data.frame(KO)
rownames(KO_df) <- KO_df[[1]]
KO_df <- KO_df[, -1]

mcra_ko <- "K00399"
if (mcra_ko %in% colnames(KO_df)) {
  mcra_otus <- rownames(KO_df)[KO_df[, mcra_ko] > 0]
  cat("mcrA-carrying OTUs (K00399):", length(mcra_otus), "\n")
} else {
  stop("K00399 (mcrA) not found in KO_predicted.tsv")
}

pmoa_ko <- "K10944"
if (pmoa_ko %in% colnames(KO_df)) {
  pmoa_otus <- rownames(KO_df)[KO_df[, pmoa_ko] > 0]
  cat("pmoA-carrying OTUs (K10944):", length(pmoa_otus), "\n")
} else {
  stop("K10944 (pmoA) not found in KO_predicted.tsv")
}

# ==============================================================================
# STEP 4: Load pathway contributions and filter plastids
# ==============================================================================

cat("\nLoading pathway contributions (this may take a few minutes)...\n")
contrib <- fread("data/raw/picrust/path_abun_contrib.tsv", header = TRUE)
cat("Loaded", nrow(contrib), "rows\n")

# Load pathway descriptions for later use
pathways_raw <- read.table("data/raw/picrust/path_abun_unstrat_descrip.tsv",
                            sep = "\t", header = TRUE, quote = "")
descriptions <- setNames(pathways_raw$description, pathways_raw$pathway)

# --- Remove plastid OTU contributions ---
n_before <- nrow(contrib)
contrib <- contrib[!(contrib$taxon %in% plastid_otus), ]
n_after <- nrow(contrib)
cat("Removed", n_before - n_after, "plastid contribution rows (",
    round(100 * (n_before - n_after) / n_before, 1), "%)\n")

# ==============================================================================
# STEP 5: Reconstruct pathway abundances (plastid-filtered)
# ==============================================================================

## Helper function: rebuild abundance matrix from contrib table
rebuild_abundance <- function(contrib_dt, sample_names) {
  # Filter to relevant samples first to reduce memory
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

## Helper: normalize to relative abundance with prevalence filter
normalize_and_filter <- function(abun_mat, min_depth = 10000, min_prevalence = 20) {
  abun_mat <- abun_mat[rowSums(abun_mat) > min_depth, , drop = FALSE]
  pct <- abun_mat / rowSums(abun_mat)
  pct <- pct[, colSums(pct > 0) > min_prevalence, drop = FALSE]
  return(pct)
}

# --- Build abundance matrices for each gene's sample set ---
cat("\n=== Rebuilding pathway abundances (plastid-filtered) ===\n")

# All non-plastid OTUs (for mcrA-all, pmoA, mmoX)
all_gene_samples <- unique(c(rownames(wood_mcra), rownames(wood_pmoa), rownames(wood_mmox)))
abun_all <- rebuild_abundance(contrib, all_gene_samples)
cat("All-OTU matrix:", nrow(abun_all), "samples x", ncol(abun_all), "pathways\n")

# No-mcrA OTUs (for mcrA no-mcrA analysis)
contrib_no_mcra <- contrib[!(contrib$taxon %in% mcra_otus), ]
cat("Contribution rows after removing mcrA OTUs:", nrow(contrib_no_mcra), "\n")
abun_no_mcra <- rebuild_abundance(contrib_no_mcra, rownames(wood_mcra))
cat("No-mcrA-OTU matrix:", nrow(abun_no_mcra), "samples x", ncol(abun_no_mcra), "pathways\n")
rm(contrib_no_mcra)

# No-pmoA OTUs (for pmoA no-pmoA analysis)
contrib_no_pmoa <- contrib[!(contrib$taxon %in% pmoa_otus), ]
cat("Contribution rows after removing pmoA OTUs:", nrow(contrib_no_pmoa), "\n")
abun_no_pmoa <- rebuild_abundance(contrib_no_pmoa, rownames(wood_pmoa))
cat("No-pmoA-OTU matrix:", nrow(abun_no_pmoa), "samples x", ncol(abun_no_pmoa), "pathways\n")

# Free memory
rm(contrib, contrib_no_pmoa)
gc()

# ==============================================================================
# STEP 6: Match samples and normalize
# ==============================================================================

## mcrA (all OTUs, plastid-filtered)
samples_mcra <- intersect(rownames(abun_all), rownames(wood_mcra))
wood_mcra <- wood_mcra[samples_mcra, ]
pct_mcra_all <- normalize_and_filter(abun_all[samples_mcra, ])
cat("\nmcrA (all OTUs): ", nrow(pct_mcra_all), "samples,", ncol(pct_mcra_all), "pathways\n")

## mcrA (no-mcrA OTUs)
samples_mcra_nm <- intersect(rownames(abun_no_mcra), rownames(wood_mcra))
pct_mcra_no <- normalize_and_filter(abun_no_mcra[samples_mcra_nm, ])
cat("mcrA (no-mcrA OTUs):", nrow(pct_mcra_no), "samples,", ncol(pct_mcra_no), "pathways\n")

## pmoA (all OTUs, plastid-filtered)
samples_pmoa <- intersect(rownames(abun_all), rownames(wood_pmoa))
wood_pmoa <- wood_pmoa[samples_pmoa, ]
pct_pmoa <- normalize_and_filter(abun_all[samples_pmoa, ])
cat("pmoA (all OTUs):", nrow(pct_pmoa), "samples,", ncol(pct_pmoa), "pathways\n")

## pmoA (no-pmoA OTUs)
samples_pmoa_np <- intersect(rownames(abun_no_pmoa), rownames(wood_pmoa))
pct_pmoa_no <- normalize_and_filter(abun_no_pmoa[samples_pmoa_np, ])
cat("pmoA (no-pmoA OTUs):", nrow(pct_pmoa_no), "samples,", ncol(pct_pmoa_no), "pathways\n")

## mmoX
samples_mmox <- intersect(rownames(abun_all), rownames(wood_mmox))
wood_mmox <- wood_mmox[samples_mmox, ]
pct_mmox <- normalize_and_filter(abun_all[samples_mmox, ])
cat("mmoX:", nrow(pct_mmox), "samples,", ncol(pct_mmox), "pathways\n")

# ==============================================================================
# STEP 7: Compute gene-OTU contribution fractions per pathway
# ==============================================================================

## Helper: compute fraction of each pathway's abundance from a set of OTUs
compute_contrib_frac <- function(abun_all_mat, abun_no_gene_mat, sample_set) {
  shared_samp <- intersect(rownames(abun_all_mat), rownames(abun_no_gene_mat))
  shared_samp <- intersect(shared_samp, sample_set)
  shared_pw <- intersect(colnames(abun_all_mat), colnames(abun_no_gene_mat))

  sapply(shared_pw, function(pw) {
    total <- sum(abun_all_mat[shared_samp, pw])
    no_gene <- sum(abun_no_gene_mat[shared_samp, pw])
    if (total == 0) return(0)
    return((total - no_gene) / total)
  })
}

cat("\n=== Computing OTU contribution fractions ===\n")

# mcrA-OTU contributions
mcra_contrib_frac <- compute_contrib_frac(
  abun_all[rownames(wood_mcra), ], abun_no_mcra, rownames(wood_mcra))
cat("mcrA-OTU contributions computed for", length(mcra_contrib_frac), "pathways\n")

# pmoA-OTU contributions
pmoa_contrib_frac <- compute_contrib_frac(
  abun_all[rownames(wood_pmoa), ], abun_no_pmoa, rownames(wood_pmoa))
cat("pmoA-OTU contributions computed for", length(pmoa_contrib_frac), "pathways\n")

rm(abun_all, abun_no_mcra, abun_no_pmoa)
gc()

# ==============================================================================
# STEP 8: Run LMER association tests
# ==============================================================================

run_lmer <- function(pct_mat, meta_df, gene_col, gene_label) {
  # Subset metadata to only samples present in the abundance matrix
  shared <- intersect(rownames(pct_mat), rownames(meta_df))
  meta_df <- meta_df[shared, ]

  cat("\nRunning LMER for", gene_label, "(",
      nrow(pct_mat), "samples,", ncol(pct_mat), "pathways)...\n")

  pvals_list <- lapply(colnames(pct_mat), function(pathway) {
    meta_copy <- meta_df
    meta_copy$pathway <- scale(pct_mat[rownames(meta_copy), pathway])
    meta_copy$gene <- scale(meta_copy[[gene_col]])

    tryCatch({
      fit0 <- lmer(gene ~ core_type + log16S + (1 | seq_id), data = meta_copy)
      fit1 <- lmer(gene ~ pathway + core_type + log16S + (1 | seq_id), data = meta_copy)
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

  # Add mean CPM
  pathway_cpm <- 1e6 * colMeans(pct_mat)
  pvals$mean_counts_per_million <- pathway_cpm[pvals$pathway]

  cat("  Significant at FDR < 0.01:", sum(pvals$FDR < 0.01, na.rm = TRUE), "\n")
  cat("  Significant at FDR < 0.05:", sum(pvals$FDR < 0.05, na.rm = TRUE), "\n")

  return(pvals)
}

## mcrA — all OTUs (plastid-filtered)
pvals_mcra_all <- run_lmer(pct_mcra_all, wood_mcra, "mcra_probe_loose", "mcrA (all OTUs)")

## mcrA — no mcrA OTUs
pvals_mcra_no <- run_lmer(pct_mcra_no, wood_mcra[rownames(pct_mcra_no), ],
                           "mcra_probe_loose", "mcrA (no-mcrA OTUs)")

## pmoA — all OTUs (plastid-filtered)
pvals_pmoa <- run_lmer(pct_pmoa, wood_pmoa, "pmoa_loose", "pmoA (all OTUs)")

## pmoA — no pmoA OTUs
pvals_pmoa_no <- run_lmer(pct_pmoa_no, wood_pmoa[rownames(pct_pmoa_no), ],
                           "pmoa_loose", "pmoA (no-pmoA OTUs)")

## mmoX
pvals_mmox <- run_lmer(pct_mmox, wood_mmox, "mmox_loose", "mmoX")

# ==============================================================================
# STEP 9: Build combined tables (all + no-gene results)
# ==============================================================================

## --- mcrA combined ---
pvals_combined <- merge(
  pvals_mcra_all %>% select(pathway, p, t, FDR, description, mean_counts_per_million),
  pvals_mcra_no %>% select(pathway, p, t, FDR),
  by = "pathway", suffixes = c(".all", ".no_mcra")
)
pvals_combined$mean_percent_from_mcra <- mcra_contrib_frac[pvals_combined$pathway]
pvals_combined$mean_percent_from_mcra[is.na(pvals_combined$mean_percent_from_mcra)] <- 0

## --- pmoA combined ---
pvals_pmoa_combined <- merge(
  pvals_pmoa %>% select(pathway, p, t, FDR, description, mean_counts_per_million),
  pvals_pmoa_no %>% select(pathway, p, t, FDR),
  by = "pathway", suffixes = c(".all", ".no_pmoa")
)
pvals_pmoa_combined$mean_percent_from_pmoa <- pmoa_contrib_frac[pvals_pmoa_combined$pathway]
pvals_pmoa_combined$mean_percent_from_pmoa[is.na(pvals_pmoa_combined$mean_percent_from_pmoa)] <- 0

# ==============================================================================
# STEP 10: Save all results
# ==============================================================================

out_dir <- "data/processed/molecular/picrust"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(arrange(pvals_mcra_all, p),
          file.path(out_dir, "pathway_associations_mcra_all.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_mcra_no, p),
          file.path(out_dir, "pathway_associations_mcra_no_mcra_otus.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_pmoa, p),
          file.path(out_dir, "pathway_associations_pmoa.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_pmoa_no, p),
          file.path(out_dir, "pathway_associations_pmoa_no_pmoa_otus.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_mmox, p),
          file.path(out_dir, "pathway_associations_mmox.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_combined, p.all),
          file.path(out_dir, "pathway_associations_combined.csv"),
          row.names = FALSE)

write.csv(arrange(pvals_pmoa_combined, p.all),
          file.path(out_dir, "pathway_associations_pmoa_combined.csv"),
          row.names = FALSE)

cat("\n=== All results saved to", out_dir, "===\n")
cat("  pathway_associations_mcra_all.csv\n")
cat("  pathway_associations_mcra_no_mcra_otus.csv\n")
cat("  pathway_associations_pmoa.csv\n")
cat("  pathway_associations_pmoa_no_pmoa_otus.csv\n")
cat("  pathway_associations_mmox.csv\n")
cat("  pathway_associations_combined.csv\n")
cat("  pathway_associations_pmoa_combined.csv\n")
cat("\nDone.\n")
