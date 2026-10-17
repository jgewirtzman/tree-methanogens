# ==============================================================================
# PICRUSt2 MetaCyc Pathway–Gene Abundance Associations
# ==============================================================================
# Purpose: Tests associations between MetaCyc pathway abundances and functional
#   gene abundances (mcrA, pmoA, mmoX) using LMER mixed-effects models, with the
#   pathway abundances from all non-organelle ASVs and, for mcrA and pmoA, without
#   the ASVs predicted to carry the gene.
#
# Pipeline stage: D (runs every time, from data/tables/)
#
# Inputs (data/tables/, via code/lib/tables_read.R):
#   - picrust_pathways.csv      per library x pathway: summed PICRUSt2 contributions of
#                               the non-organelle ASVs (abundance), and without the
#                               mcrA- (K00399) or pmoA- (K10944) carrying ASVs; built in
#                               stage A by code/01_import/06_picrust_pathways.R from the
#                               PICRUSt2 predictions (data/raw/picrust/)
#   - picrust_pathway_info.csv  pathway descriptions
#   - sample metadata: survey_16s_meta() (samples, samples_16S, ddpcr, qpcr, trees)
#
# Approach:
#   1. Sample sets: wood samples, 16S >= 100 copies per uL, non-NA gene abundance
#   2. Abundance matrices: (a) all non-organelle ASVs, (b) without mcrA ASVs,
#      (c) without pmoA ASVs
#   3. Relative abundance; depth > 10,000; pathways present in > 20 samples
#   4. LMER for each pathway against mcrA (both versions), pmoA (both), mmoX
#   5. Share of each pathway's abundance from the gene-carrying ASVs
#
# Outputs (outputs/data/picrust/):
#   - pathway_associations_mcra_all.csv
#   - pathway_associations_mcra_no_mcra_otus.csv
#   - pathway_associations_pmoa.csv
#   - pathway_associations_pmoa_no_pmoa_otus.csv
#   - pathway_associations_mmox.csv
#   - pathway_associations_combined.csv       (merged mcrA all + no-mcrA results)
#   - pathway_associations_pmoa_combined.csv  (merged pmoA all + no-pmoA results)
#
# Required packages: tidyverse, lme4
# ==============================================================================

library(tidyverse)
library(lme4)

# ==============================================================================
# STEP 1: Load metadata and identify sample sets for each gene
# ==============================================================================

source("code/lib/tables_read.R")
meta <- survey_16s_meta()   # per sequenced 2021 sample (data/tables/: samples, samples_16S, ddpcr, qpcr, trees)
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
# STEP 2: Pathway abundances per library (data/tables/picrust_pathways.csv)
# ==============================================================================
# Organelle ASVs are already out (code/01_import/06_picrust_pathways.R).
PW <- read_table("picrust_pathways")
descriptions <- with(read_table("picrust_pathway_info"), setNames(description, pathway))

## Helper: library x pathway matrix of one abundance version (pathways with a
## non-zero abundance in the libraries; 0 where a library has none)
rebuild_abundance <- function(col, sample_names) {
  sub <- PW[PW$library_id %in% sample_names & PW[[col]] > 0, c("library_id", "pathway", col)]
  agg <- sub %>%
    pivot_wider(names_from = pathway, values_from = all_of(col), values_fill = 0) %>%
    as.data.frame()
  rownames(agg) <- agg$library_id
  agg[, -1, drop = FALSE]
}

## Helper: normalize to relative abundance with prevalence filter
normalize_and_filter <- function(abun_mat, min_depth = 10000, min_prevalence = 20) {
  abun_mat <- abun_mat[rowSums(abun_mat) > min_depth, , drop = FALSE]
  pct <- abun_mat / rowSums(abun_mat)
  pct <- pct[, colSums(pct > 0) > min_prevalence, drop = FALSE]
  return(pct)
}

# --- Build abundance matrices for each gene's sample set ---
cat("\n=== Pathway abundance matrices (non-organelle ASVs) ===\n")

# All non-organelle ASVs (for mcrA-all, pmoA, mmoX)
all_gene_samples <- unique(c(rownames(wood_mcra), rownames(wood_pmoa), rownames(wood_mmox)))
abun_all <- rebuild_abundance("abundance", all_gene_samples)
cat("All-ASV matrix:", nrow(abun_all), "samples x", ncol(abun_all), "pathways\n")

# Without the mcrA-carrying ASVs (for mcrA no-mcrA analysis)
abun_no_mcra <- rebuild_abundance("abundance_no_mcra_asvs", rownames(wood_mcra))
cat("No-mcrA-ASV matrix:", nrow(abun_no_mcra), "samples x", ncol(abun_no_mcra), "pathways\n")

# Without the pmoA-carrying ASVs (for pmoA no-pmoA analysis)
abun_no_pmoa <- rebuild_abundance("abundance_no_pmoa_asvs", rownames(wood_pmoa))
cat("No-pmoA-ASV matrix:", nrow(abun_no_pmoa), "samples x", ncol(abun_no_pmoa), "pathways\n")
rm(PW)

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

out_dir <- "outputs/data/picrust"
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
