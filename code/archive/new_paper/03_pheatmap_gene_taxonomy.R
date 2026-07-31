# ==============================================================================
# Taxonomy–Gene Association Heatmaps (Family-level 16S) — pmoA and mmoX
# ==============================================================================
# Purpose: Pheatmap-style heatmaps showing which bacterial/archaeal families
#   (all families, not just methanotrophs) are significantly associated with
#   pmoA and mmoX ddPCR gene abundances via LMER mixed-effects models.
#
# Adapted from: code/06_figures/12c_taxonomy_pmoa_heatmap.R
#
# Approach:
#   - Wood samples only, 16S qPCR >= 100 copies/µL
#   - Family-level aggregation from rarefied 16S data
#   - LMER: gene ~ family_abundance + core_type + log16S + (1 | seq_id)
#   - Top 20 families by mean abundance + any FDR < 0.05 families
#   - Families ordered by t-statistic; samples ordered by gene abundance
#
# Outputs:
#   - fig_pheatmap_pmoa.png         (pmoA associations)
#   - fig_pheatmap_mmox.png         (mmoX associations)
#   - fig_pheatmap_comparison.png   (pmoA vs mmoX t-stat comparison)
#   - family_pmoa_associations_new.csv
#   - family_mmox_associations_new.csv
#
# Required packages: phyloseq, tidyverse, lme4, pheatmap, RColorBrewer
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(lme4)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# STEP 1: Build phyloseq object (consistent with 08b_ and 12c_ scripts)
# ==============================================================================

ddpcr <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
water <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv")
qpcr_16s <- subset(qpcr_16s, qpcr_16s$Sample.ID != "None")
qpcr_16s <- qpcr_16s[, c(3, 4, 6)]

water$seq_id <- toupper(water$seq_id)
ddpcr <- merge(ddpcr, water, by = c("core_type", "seq_id"), all.x = TRUE)
ddpcr <- merge(ddpcr, qpcr_16s, by = c("core_type", "seq_id"), all.x = TRUE)

otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)
bastard_tax <- otu_tab[, 590:596]
bastard_tax[bastard_tax == ""] <- NA
tax_tab_pre <- tax_table(bastard_tax)
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))

otu_tab_corr <- otu_tab[, 1:589]
otu_table_pre <- otu_table(otu_tab_corr, taxa_are_rows = TRUE)

phylo_tree <- read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
samp_data$RowName <- row.names(samp_data)

samp_data$seq_id <- sub("prime", "'", samp_data$seq_id)
samp_data$seq_id <- sub("star", "*", samp_data$seq_id)
samp_data$seq_id <- sub("HM", "H", samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
dups <- which(duplicated(samp_data_merged$RowName) == TRUE)
samp_data_merged <- samp_data_merged[-c(dups), ]
row.names(samp_data_merged) <- samp_data_merged$RowName

raw_ps <- phyloseq(tax_tab_pre, otu_table_pre, phylo_tree, sample_data(samp_data_merged))

# ==============================================================================
# STEP 2: Remove mitochondria and chloroplasts
# ==============================================================================

pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
badTaxa <- c(mitochondria, chloroplast)
no_mito <- pop_taxa(raw_ps, badTaxa)

cat("Removed", length(mitochondria), "mitochondrial and",
    length(chloroplast), "chloroplast taxa\n")

# ==============================================================================
# STEP 3: Rarefy and transform to relative abundance
# ==============================================================================

taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)
ps.ra <- transform_sample_counts(ps.rare, function(x) 100 * x / sum(x))

colnames(tax_table(ps.ra)) <- c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species")

cat("After rarefaction:", nsamples(ps.ra), "samples,", ntaxa(ps.ra), "taxa\n")

# ==============================================================================
# STEP 4: Filter to wood samples with gene data and 16S >= 100
# ==============================================================================

ps.wood <- subset_samples(ps.ra, material == "Wood")

meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)
meta$compartment <- ifelse(meta$core_type == "Inner", "Heartwood", "Sapwood")

# Well-sampled species only
well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                   "FAGR", "FRAM", "PIST", "QURU", "TSCA")

# Filter: wood, has pmoA, has mmoX, 16S >= 100, well-sampled species
wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(pmoa_loose), !is.na(mmox_loose)) %>%
  filter(X16S_per_ul >= 100) %>%
  filter(species.x %in% well_sampled)

ps_samples <- sample_names(ps.wood)
meta_samples <- rownames(wood)
shared_samples <- intersect(ps_samples, meta_samples)

wood <- wood[shared_samples, ]
ps.wood <- prune_samples(shared_samples, ps.wood)

cat("Wood samples with pmoA + mmoX + 16S >= 100:", nrow(wood), "\n")
cat("By compartment:\n")
print(table(wood$compartment))

# ==============================================================================
# STEP 5: Aggregate to Family level
# ==============================================================================

otu_ra <- as.data.frame(t(otu_table(ps.wood)))  # samples x taxa
tax_df <- as.data.frame(tax_table(ps.wood))

families <- tax_df$Family
families <- gsub("^\\s+|\\s+$", "", families)
families[is.na(families) | families == "" | families == "none"] <- NA
names(families) <- rownames(tax_df)

valid_families <- unique(families[!is.na(families)])
family_abundance <- matrix(0, nrow = nrow(otu_ra), ncol = length(valid_families))
rownames(family_abundance) <- rownames(otu_ra)
colnames(family_abundance) <- valid_families

for (i in seq_len(ncol(otu_ra))) {
  taxon <- colnames(otu_ra)[i]
  fam <- families[taxon]
  if (!is.na(fam) && fam %in% colnames(family_abundance)) {
    family_abundance[, fam] <- family_abundance[, fam] + otu_ra[, i]
  }
}

cat("Unique families:", length(valid_families), "\n")

# ==============================================================================
# STEP 6: Run LMER association tests for pmoA AND mmoX
# ==============================================================================

## Relative abundance at family level (re-normalize within families)
percent_all <- family_abundance / rowSums(family_abundance)

## Prevalence filter: families present in > 20 samples
percent_filt <- percent_all[, colSums(percent_all > 0) > 20]
cat("\nFamilies passing prevalence filter (> 20 samples):", ncol(percent_filt), "\n")

run_lmer_associations <- function(percent_mat, wood_df, gene_col) {
  pvals_list <- lapply(colnames(percent_mat), function(family) {
    meta_copy <- wood_df
    meta_copy$family <- scale(percent_mat[, family])
    meta_copy$gene <- scale(meta_copy[[gene_col]])

    tryCatch({
      fit0 <- lmer(gene ~ core_type + log16S + (1 | seq_id), data = meta_copy)
      fit1 <- lmer(gene ~ family + core_type + log16S + (1 | seq_id), data = meta_copy)
      p <- anova(fit0, fit1)[2, "Pr(>Chisq)"]
      t_val <- coef(summary(fit1))["family", "t value"]
      data.frame(family = family, p = p, t = t_val, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(family = family, p = NA, t = NA, stringsAsFactors = FALSE)
    })
  })
  pvals <- do.call("rbind", pvals_list)
  pvals$FDR <- p.adjust(pvals$p, "BH")
  return(pvals)
}

cat("Running pmoA LMER associations...\n")
pvals_pmoa <- run_lmer_associations(percent_filt, wood, "pmoa_loose")
cat("  Done:", nrow(pvals_pmoa), "families tested\n")

cat("Running mmoX LMER associations...\n")
pvals_mmox <- run_lmer_associations(percent_filt, wood, "mmox_loose")
cat("  Done:", nrow(pvals_mmox), "families tested\n")

# Add mean relative abundance
mean_abundance <- colMeans(percent_filt)
pvals_pmoa$mean_pct <- mean_abundance[pvals_pmoa$family]
pvals_mmox$mean_pct <- mean_abundance[pvals_mmox$family]

# Save association tables
write.csv(arrange(pvals_pmoa, p),
          "outputs/figures/new_paper/family_pmoa_associations_new.csv", row.names = FALSE)
write.csv(arrange(pvals_mmox, p),
          "outputs/figures/new_paper/family_mmox_associations_new.csv", row.names = FALSE)

# Report results
cat("\n=== pmoA associations ===\n")
cat("  FDR < 0.01:", sum(pvals_pmoa$FDR < 0.01, na.rm = TRUE), "\n")
cat("  FDR < 0.05:", sum(pvals_pmoa$FDR < 0.05, na.rm = TRUE), "\n")
cat("\n=== mmoX associations ===\n")
cat("  FDR < 0.01:", sum(pvals_mmox$FDR < 0.01, na.rm = TRUE), "\n")
cat("  FDR < 0.05:", sum(pvals_mmox$FDR < 0.05, na.rm = TRUE), "\n")

# ==============================================================================
# STEP 7: Select families for each heatmap
# ==============================================================================

select_heatmap_families <- function(pvals_df, mean_abund, n_top = 20) {
  top_n <- names(sort(mean_abund, decreasing = TRUE))[1:n_top]
  sig_005 <- pvals_df %>% filter(!is.na(FDR), FDR < 0.05) %>% pull(family)
  display <- union(top_n, sig_005)
  cat("  Top", n_top, "by abundance:", length(top_n), "\n")
  cat("  Significant (FDR < 0.05):", length(sig_005), "\n")
  cat("  Additional significant not in top", n_top, ":",
      length(setdiff(sig_005, top_n)), "\n")
  cat("  Total families to display:", length(display), "\n")
  return(display)
}

cat("\npmoA heatmap families:\n")
pmoa_families <- select_heatmap_families(pvals_pmoa, mean_abundance)

cat("\nmmoX heatmap families:\n")
mmox_families <- select_heatmap_families(pvals_mmox, mean_abundance)

# ==============================================================================
# STEP 8: Build and save pheatmaps
# ==============================================================================

build_pheatmap <- function(percent_mat, wood_df, pvals_df, display_families,
                           gene_col, gene_label, pvals_other = NULL,
                           other_label = NULL, out_file) {

  ## Order samples by gene abundance
  ordered_samples <- rownames(wood_df[order(wood_df[[gene_col]]), ])
  heatmap_data <- percent_mat[ordered_samples, display_families, drop = FALSE]

  ## Column annotations: compartment, species, gene level
  wood_df$log10_gene <- log10(1 + wood_df[[gene_col]])

  species_mapping <- c(
    "ACRU" = "A. rubrum", "ACSA" = "A. saccharum",
    "BEAL" = "B. alleghaniensis", "BELE" = "B. lenta",
    "BEPA" = "B. papyrifera", "FAGR" = "F. grandifolia",
    "FRAM" = "F. americana", "PIST" = "P. strobus",
    "QURU" = "Q. rubra", "TSCA" = "T. canadensis"
  )
  wood_df$species_short <- species_mapping[wood_df$species.x]

  ann_col <- wood_df[ordered_samples, c("compartment", "species_short", "log10_gene"),
                     drop = FALSE]
  colnames(ann_col) <- c("Compartment", "Species", paste0("log10(", gene_label, ")"))

  ## Row annotations: significance for this gene (and other gene if provided)
  pvals_lookup <- pvals_df %>% filter(family %in% display_families)
  fdr_vec <- setNames(pvals_lookup$FDR, pvals_lookup$family)

  sig_status <- sapply(display_families, function(f) {
    fdr <- fdr_vec[f]
    if (is.na(fdr)) return("NS")
    if (fdr < 0.01) return("FDR < 0.01")
    if (fdr < 0.05) return("FDR < 0.05")
    return("NS")
  })

  ann_row <- data.frame(
    sig = sig_status,
    check.names = FALSE
  )
  colnames(ann_row) <- paste0(gene_label, " assoc.")
  rownames(ann_row) <- display_families

  ## Add other gene's significance if provided
  if (!is.null(pvals_other) && !is.null(other_label)) {
    pvals_lookup2 <- pvals_other %>% filter(family %in% display_families)
    fdr_vec2 <- setNames(pvals_lookup2$FDR, pvals_lookup2$family)
    sig_status2 <- sapply(display_families, function(f) {
      fdr <- fdr_vec2[f]
      if (is.na(fdr)) return("NS")
      if (fdr < 0.01) return("FDR < 0.01")
      if (fdr < 0.05) return("FDR < 0.05")
      return("NS")
    })
    ann_row[[paste0(other_label, " assoc.")]] <- sig_status2
  }

  ## Color schemes
  # Species colors (10 species)
  species_cols <- setNames(
    c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
      "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#469990"),
    sort(unique(na.omit(wood_df$species_short)))
  )

  ann_colors <- list(
    Compartment = c(Heartwood = "#a6611a", Sapwood = "#dfc27d"),
    Species = species_cols
  )
  ann_colors[[paste0("log10(", gene_label, ")")]] <-
    colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100)
  ann_colors[[paste0(gene_label, " assoc.")]] <-
    c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a", "NS" = "#f0f0f0")

  if (!is.null(other_label)) {
    ann_colors[[paste0(other_label, " assoc.")]] <-
      c("FDR < 0.01" = "#1a237e", "FDR < 0.05" = "#90caf9", "NS" = "#f0f0f0")
  }

  ## Heatmap body colors
  heatmap_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea",
                                        "#f7f7f7",
                                        "#f5c5a3", "#e8845a", "#b40426"))(100)

  ## Order families by t-statistic
  pvals_ordered <- pvals_df %>%
    filter(family %in% display_families) %>%
    arrange(t)
  ordered_families <- pvals_ordered$family

  heatmap_data <- heatmap_data[, ordered_families]
  ann_row <- ann_row[ordered_families, , drop = FALSE]

  ## Dimensions
  fig_height <- max(7, length(display_families) * 0.35 + 3)
  fig_width <- 14

  ## Draw heatmap
  png(out_file, width = fig_width, height = fig_height, units = "in", res = 300)

  pheatmap(t(heatmap_data),
           color = heatmap_colors,
           scale = "row",
           annotation_col = ann_col,
           annotation_row = ann_row,
           annotation_colors = ann_colors,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           treeheight_row = 0,
           labels_row = ordered_families,
           annotation_names_row = TRUE,
           fontsize_row = 9,
           fontsize = 10,
           fontsize_col = 8,
           border_color = NA,
           main = "")

  dev.off()
  cat("  Saved:", out_file, "\n")
}

## Build pmoA heatmap (with mmoX cross-annotations)
cat("\nBuilding pmoA pheatmap...\n")
build_pheatmap(percent_filt, wood, pvals_pmoa, pmoa_families,
               gene_col = "pmoa_loose", gene_label = "pmoA",
               pvals_other = pvals_mmox, other_label = "mmoX",
               out_file = "outputs/figures/new_paper/fig_pheatmap_pmoa.png")

## Build mmoX heatmap (with pmoA cross-annotations)
cat("Building mmoX pheatmap...\n")
build_pheatmap(percent_filt, wood, pvals_mmox, mmox_families,
               gene_col = "mmox_loose", gene_label = "mmoX",
               pvals_other = pvals_pmoa, other_label = "pmoA",
               out_file = "outputs/figures/new_paper/fig_pheatmap_mmox.png")

# ==============================================================================
# STEP 9: pmoA vs mmoX comparison heatmap (t-statistics side by side)
# ==============================================================================

# Merge t-statistics for all families tested in both
comparison <- merge(
  pvals_pmoa %>% select(family, t_pmoa = t, fdr_pmoa = FDR, mean_pct),
  pvals_mmox %>% select(family, t_mmox = t, fdr_mmox = FDR),
  by = "family"
)

# Significance categories
comparison$sig_cat <- with(comparison, case_when(
  fdr_pmoa < 0.05 & fdr_mmox < 0.05 ~ "Both",
  fdr_pmoa < 0.05 ~ "pmoA only",
  fdr_mmox < 0.05 ~ "mmoX only",
  TRUE ~ "Neither"
))

cat("\n=== Comparison of pmoA vs mmoX family associations ===\n")
print(table(comparison$sig_cat))

# Select families: significant for either gene, or top 20 by abundance
top20_comp <- comparison %>% arrange(desc(mean_pct)) %>% slice(1:20) %>% pull(family)
sig_comp <- comparison %>%
  filter(fdr_pmoa < 0.05 | fdr_mmox < 0.05) %>%
  pull(family)
display_comp <- union(top20_comp, sig_comp)

comp_disp <- comparison %>%
  filter(family %in% display_comp) %>%
  arrange(t_pmoa)

# Build comparison heatmap data: families × 2 columns (pmoA t, mmoX t)
comp_mat <- as.matrix(comp_disp[, c("t_pmoa", "t_mmox")])
rownames(comp_mat) <- comp_disp$family
colnames(comp_mat) <- c("pmoA t-stat", "mmoX t-stat")

# Row annotations: significance category
ann_row_comp <- data.frame(
  `Sig. gene` = comp_disp$sig_cat,
  check.names = FALSE
)
rownames(ann_row_comp) <- comp_disp$family

ann_colors_comp <- list(
  `Sig. gene` = c("Both" = "#7b1fa2", "pmoA only" = "#b71c1c",
                   "mmoX only" = "#1a237e", "Neither" = "#f0f0f0")
)

comp_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea",
                                   "#f7f7f7",
                                   "#f5c5a3", "#e8845a", "#b40426"))(100)

fig_height_comp <- max(7, length(display_comp) * 0.35 + 3)

png("outputs/figures/new_paper/fig_pheatmap_comparison.png",
    width = 7, height = fig_height_comp, units = "in", res = 300)

pheatmap(comp_mat,
         color = comp_colors,
         scale = "none",
         annotation_row = ann_row_comp,
         annotation_colors = ann_colors_comp,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         treeheight_row = 0,
         fontsize_row = 9,
         fontsize = 11,
         fontsize_col = 11,
         border_color = NA,
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "black",
         main = "Family–gene LMER associations (t-statistics)")

dev.off()
cat("Saved: fig_pheatmap_comparison.png\n")

# ==============================================================================
# STEP 10: Print top associations for interpretation
# ==============================================================================

cat("\n=== Top 10 pmoA-associated families (by |t|) ===\n")
top_pmoa <- pvals_pmoa %>%
  filter(!is.na(FDR)) %>%
  arrange(desc(abs(t))) %>%
  slice(1:10)
print(top_pmoa[, c("family", "t", "FDR", "mean_pct")], row.names = FALSE)

cat("\n=== Top 10 mmoX-associated families (by |t|) ===\n")
top_mmox <- pvals_mmox %>%
  filter(!is.na(FDR)) %>%
  arrange(desc(abs(t))) %>%
  slice(1:10)
print(top_mmox[, c("family", "t", "FDR", "mean_pct")], row.names = FALSE)

cat("\n=== Families significant for BOTH genes ===\n")
both_sig <- comparison %>% filter(sig_cat == "Both")
if (nrow(both_sig) > 0) {
  print(both_sig[, c("family", "t_pmoa", "fdr_pmoa", "t_mmox", "fdr_mmox", "mean_pct")],
        row.names = FALSE)
} else {
  cat("  None\n")
}

cat("\nDone.\n")
