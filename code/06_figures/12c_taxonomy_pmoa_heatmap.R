# ==============================================================================
# Taxonomy–pmoA Association Heatmaps (Family-level 16S) — Figure S15
# ==============================================================================
# Purpose: Visualizes the top bacterial/archaeal families (16S) across wood
#   samples, highlighting those significantly associated with pmoA (methanotroph)
#   abundance via LMER mixed-effects models.
#
# Pipeline stage: 4 — Publication Figures
#
# Approach: Shows the top 20 most abundant families (by mean relative
#   abundance) + any additional families significantly associated with pmoA
#   (FDR < 0.05) that are not already in the top 20. Families are ordered by
#   t-statistic; samples are ordered by pmoA abundance.
#
# Filtering (consistent with all 16S sequencing scripts):
#   1. Remove mitochondria and chloroplasts (plastid filtering)
#   2. Rarefy to 3,500 reads per sample (seed = 46814)
#   3. Transform to relative abundance
#   4. Filter to wood samples with pmoA data and 16S >= 100
#
# Inputs:
#   - OTU table:     data/raw/16s/OTU_table.txt
#   - Phylo tree:    data/raw/16s/unrooted_tree.nwk
#   - Sample map:    data/raw/16s/tree_16s_mapping_dada2_corrected.txt
#   - ddPCR meta:    data/raw/ddpcr/ddPCR_meta_all_data.csv
#   - Core data:     data/raw/tree_cores/Tree_Core_Sectioning_Data.csv
#   - 16S qPCR:      data/raw/16s/16s_w_metadata.csv
#   - PICRUSt meta:  data/raw/picrust/16S_tree_sample_table_with_meta.csv
#
# Outputs:
#   - figS15_taxonomy_pmoa_heatmap.png
#   - family_pmoa_associations.csv
#
# Required packages: phyloseq, tidyverse, lme4, pheatmap, RColorBrewer
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(lme4)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# STEP 1: Build phyloseq object (consistent with 08_/08b_/08c_/08d_ scripts)
# ==============================================================================

# Load and merge metadata
ddpcr <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
water <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv")
qpcr_16s <- subset(qpcr_16s, qpcr_16s$Sample.ID != "None")
qpcr_16s <- qpcr_16s[, c(3, 4, 6)]

water$seq_id <- toupper(water$seq_id)
ddpcr <- merge(ddpcr, water, by = c("core_type", "seq_id"), all.x = TRUE)
ddpcr <- merge(ddpcr, qpcr_16s, by = c("core_type", "seq_id"), all.x = TRUE)

# Parse OTU table (taxonomy embedded in columns 590-596)
otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)
bastard_tax <- otu_tab[, 590:596]
bastard_tax[bastard_tax == ""] <- NA
tax_tab_pre <- tax_table(bastard_tax)
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))

otu_tab_corr <- otu_tab[, 1:589]
otu_table_pre <- otu_table(otu_tab_corr, taxa_are_rows = TRUE)

# Load tree and sample metadata
phylo_tree <- read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
samp_data$RowName <- row.names(samp_data)

# Fix sample names
samp_data$seq_id <- sub("prime", "'", samp_data$seq_id)
samp_data$seq_id <- sub("star", "*", samp_data$seq_id)
samp_data$seq_id <- sub("HM", "H", samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

# Merge and create phyloseq object
samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
dups <- which(duplicated(samp_data_merged$RowName) == TRUE)
samp_data_merged <- samp_data_merged[-c(dups), ]
row.names(samp_data_merged) <- samp_data_merged$RowName

raw_ps <- phyloseq(tax_tab_pre, otu_table_pre, phylo_tree, sample_data(samp_data_merged))

# ==============================================================================
# STEP 2: Remove mitochondria and chloroplasts (plastid filtering)
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

# Transform to relative abundance (0-100% scale)
ps.ra <- transform_sample_counts(ps.rare, function(x) 100 * x / sum(x))

# Set proper taxonomy column names
colnames(tax_table(ps.ra)) <- c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species")

cat("After rarefaction:", nsamples(ps.ra), "samples,", ntaxa(ps.ra), "taxa\n")

# ==============================================================================
# STEP 4: Filter to wood samples and prepare metadata for LMER
# ==============================================================================

ps.wood <- subset_samples(ps.ra, material == "Wood")

# Load PICRUSt metadata for pmoA + 16S qPCR values used in LMER
meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$log16S <- log10(1 + meta$X16S_per_ul)

wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(pmoa_loose)) %>%
  filter(X16S_per_ul >= 100)

# Match samples present in both phyloseq and metadata
ps_samples <- sample_names(ps.wood)
meta_samples <- rownames(wood)
shared_samples <- intersect(ps_samples, meta_samples)

wood <- wood[shared_samples, ]
ps.wood <- prune_samples(shared_samples, ps.wood)

cat("Wood samples with pmoA + 16S >= 100:", nrow(wood), "\n")

# ==============================================================================
# STEP 5: Aggregate to Family level from rarefied/filtered phyloseq
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
# STEP 6: Run LMER association tests (pmoA)
# ==============================================================================

## Relative abundance at family level
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

pvals <- run_lmer_associations(percent_filt, wood, "pmoa_loose")
cat("LMER tests complete for", nrow(pvals), "families\n")

# Compute mean relative abundance per family
mean_abundance <- colMeans(percent_filt)
pvals$mean_pct <- mean_abundance[pvals$family]

## Save full association table
write.csv(arrange(pvals, p),
          "outputs/tables/family_pmoa_associations.csv", row.names = FALSE)

# ==============================================================================
# STEP 7: Select families for heatmap (top 20 + significant)
# ==============================================================================

# Top 20 families by mean relative abundance
top20 <- names(sort(mean_abundance, decreasing = TRUE))[1:20]

# Families significantly associated with pmoA (FDR < 0.05)
sig_001 <- pvals %>% filter(!is.na(FDR), FDR < 0.01) %>% pull(family)
sig_005 <- pvals %>% filter(!is.na(FDR), FDR < 0.05) %>% pull(family)

# Union: top 20 + any significant families (FDR < 0.05) not already in top 20
display_families <- union(top20, sig_005)

cat("\nFamilies in heatmap:", length(display_families), "\n")
cat("  Top 20 by abundance:", length(top20), "\n")
cat("  Significant (FDR < 0.01):", length(sig_001), "\n")
cat("  Significant (FDR < 0.05):", length(sig_005), "\n")
cat("  Additional significant not in top 20:",
    length(setdiff(sig_005, top20)), "\n")

# ==============================================================================
# STEP 8: Build heatmap data and annotations
# ==============================================================================

## Order samples by pmoA abundance
ordered_samples <- rownames(wood[order(wood$pmoa_loose), ])

## Subset abundance matrix to display families
heatmap_data <- percent_filt[ordered_samples, display_families]

## --- Column annotations (samples) ---
wood$log10_pmoa <- log10(1 + wood$pmoa_loose)
wood$compartment <- ifelse(wood$core_type == "Inner", "Heartwood", "Sapwood")

ann_col <- wood[ordered_samples, c("compartment", "log10_pmoa"), drop = FALSE]
colnames(ann_col) <- c("Compartment", "log10(pmoA)")

## --- Row annotations (families): just significance ---
pvals_lookup <- pvals %>% filter(family %in% display_families)
fdr_vec <- setNames(pvals_lookup$FDR, pvals_lookup$family)

sig_status <- sapply(display_families, function(f) {
  fdr <- fdr_vec[f]
  if (is.na(fdr)) return("NS")
  if (fdr < 0.01) return("FDR < 0.01")
  if (fdr < 0.05) return("FDR < 0.05")
  return("NS")
})

ann_row <- data.frame(
  `pmoA assoc.` = sig_status,
  check.names = FALSE
)
rownames(ann_row) <- display_families

## --- Color schemes (clean, harmonious) ---
ann_colors <- list(
  Compartment = c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
  `log10(pmoA)` = colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
  `pmoA assoc.` = c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a", "NS" = "#f0f0f0")
)

## --- Heatmap color palette (cool blue – white – warm red, diverging) ---
heatmap_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea",
                                      "#f7f7f7",
                                      "#f5c5a3", "#e8845a", "#b40426"))(100)

## Figure dimensions
fig_height <- max(7, length(display_families) * 0.35 + 3)
fig_width <- 13

# ==============================================================================
# STEP 9: Generate and save heatmap
# ==============================================================================

## Order families by mean t-statistic (negative to positive association)
pvals_ordered <- pvals %>%
  filter(family %in% display_families) %>%
  arrange(t)
ordered_families <- pvals_ordered$family

## Reorder heatmap data and annotations to match
heatmap_data <- heatmap_data[, ordered_families]
ann_row <- ann_row[ordered_families, , drop = FALSE]

png("outputs/figures/supplementary/figS15_taxonomy_pmoa_heatmap.png",
    width = fig_width, height = fig_height, units = "in", res = 300)

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

cat("\nHeatmap saved: figS15_taxonomy_pmoa_heatmap.png\n")
cat("Done.\n")
