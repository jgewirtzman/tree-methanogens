# ==============================================================================
# figS09_taxonomy-pmoa.R -- Figure S11: 16S families associated with pmoA in wood
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Row-scaled relative abundance of the 20 most abundant families, plus any other family
# associated with pmoA at FDR < 0.05, across wood samples ordered by pmoA; families
# ordered by the association t. Samples: 2021 wood samples with pmoA and 16S qPCR
# >= 100 copies per uL. 16S: organelles removed, rarefied (rarefy_16s()), percent.
# pmoA (pmoa_loose): ddPCR copies per 20 uL well, survey_16s_meta().
# The same code as figS14_taxonomy-mcra.R; only the parameters below differ.
# Inputs: the 16S phyloseq object (code/lib/build_phyloseq.R: data/tables/ and the
#   ASV tree data/raw/16s/unrooted_tree.nwk); data/tables/ via code/lib/tables_read.R
# Outputs: outputs/tables/family_pmoa_associations.csv
#          outputs/figures/original/supplementary/figS2_taxonomy_pmoa_heatmap.png
# ==============================================================================
GENE_COL   <- "pmoa_loose"
GENE_LABEL <- "pmoA"
OUT_TABLE  <- "outputs/tables/family_pmoa_associations.csv"
OUT_PNG    <- "outputs/figures/original/supplementary/figS2_taxonomy_pmoa_heatmap.png"

library(phyloseq)
library(tidyverse)
library(lme4)
library(pheatmap)
source("code/lib/build_phyloseq.R")
source("code/lib/tables_read.R")

# ==============================================================================
# STEP 1: 16S relative abundance (rarefied), wood samples with the gene and 16S >= 100
# ==============================================================================

.ps16   <- build_16s_phyloseq()
no_mito <- .ps16$ps
cat("Removed", .ps16$n_mitochondria, "mitochondrial and",
    .ps16$n_chloroplast, "chloroplast taxa\n")
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
ps.rare <- rarefy_16s(no_mito)
ps.ra <- transform_sample_counts(ps.rare, function(x) 100 * x / sum(x))   # percent
cat("After rarefaction:", nsamples(ps.ra), "samples,", ntaxa(ps.ra), "taxa\n")

ps.wood <- subset_samples(ps.ra, material == "Wood")

# per sequenced 2021 sample: ddPCR (copies per 20 uL well) and 16S qPCR used in the LMER
meta <- survey_16s_meta()
meta$log16S <- log10(1 + meta$X16S_per_ul)
wood <- meta %>%
  filter(material == "Wood", !is.na(material)) %>%
  filter(!is.na(.data[[GENE_COL]])) %>%
  filter(X16S_per_ul >= 100)

shared_samples <- intersect(sample_names(ps.wood), rownames(wood))
wood <- wood[shared_samples, ]
ps.wood <- prune_samples(shared_samples, ps.wood)
cat("Wood samples with", GENE_LABEL, "+ 16S >= 100:", nrow(wood), "\n")

# ==============================================================================
# STEP 2: Aggregate to family level
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
  fam <- families[colnames(otu_ra)[i]]
  if (!is.na(fam)) family_abundance[, fam] <- family_abundance[, fam] + otu_ra[, i]
}
cat("Unique families:", length(valid_families), "\n")

# ==============================================================================
# STEP 3: LMER association of each family with the gene
# ==============================================================================
# gene ~ family + core_type + log16S + (1 | tree), against the model without the
# family (likelihood-ratio test); both variables scaled. BH across families.

percent_all <- family_abundance / rowSums(family_abundance)
percent_filt <- percent_all[, colSums(percent_all > 0) > 20]   # present in > 20 samples
cat("\nFamilies passing prevalence filter (> 20 samples):", ncol(percent_filt), "\n")

pvals <- do.call(rbind, lapply(colnames(percent_filt), function(family) {
  d <- wood
  d$family <- scale(percent_filt[, family])
  d$gene <- scale(d[[GENE_COL]])
  fit0 <- lmer(gene ~ core_type + log16S + (1 | seq_id), data = d)
  fit1 <- lmer(gene ~ family + core_type + log16S + (1 | seq_id), data = d)
  data.frame(family = family, p = anova(fit0, fit1)[2, "Pr(>Chisq)"],
             t = coef(summary(fit1))["family", "t value"], stringsAsFactors = FALSE)
}))
pvals$FDR <- p.adjust(pvals$p, "BH")
cat("LMER tests complete for", nrow(pvals), "families\n")

mean_abundance <- colMeans(percent_filt)
pvals$mean_pct <- mean_abundance[pvals$family]
write.csv(arrange(pvals, p), OUT_TABLE, row.names = FALSE)

# ==============================================================================
# STEP 4: Families shown: the 20 most abundant + any other with FDR < 0.05
# ==============================================================================

top20 <- names(sort(mean_abundance, decreasing = TRUE))[1:20]
sig_001 <- pvals %>% filter(!is.na(FDR), FDR < 0.01) %>% pull(family)
sig_005 <- pvals %>% filter(!is.na(FDR), FDR < 0.05) %>% pull(family)
display_families <- union(top20, sig_005)

cat("\nFamilies in heatmap:", length(display_families), "\n")
cat("  Top 20 by abundance:", length(top20), "\n")
cat("  Significant (FDR < 0.01):", length(sig_001), "\n")
cat("  Significant (FDR < 0.05):", length(sig_005), "\n")
cat("  Additional significant not in top 20:",
    length(setdiff(sig_005, top20)), "\n")

# ==============================================================================
# STEP 5: Heatmap: samples ordered by the gene, families by t
# ==============================================================================

ordered_samples <- rownames(wood[order(wood[[GENE_COL]]), ])
heatmap_data <- percent_filt[ordered_samples, display_families]

LOG_LABEL <- paste0("log10(", GENE_LABEL, ")"); ASSOC_LABEL <- paste0(GENE_LABEL, " assoc.")
wood$log10_gene <- log10(1 + wood[[GENE_COL]])
wood$compartment <- ifelse(wood$core_type == "Inner", "Heartwood", "Sapwood")
ann_col <- wood[ordered_samples, c("compartment", "log10_gene"), drop = FALSE]
colnames(ann_col) <- c("Compartment", LOG_LABEL)

fdr_vec <- setNames(pvals$FDR, pvals$family)[display_families]
ann_row <- data.frame(ifelse(is.na(fdr_vec), "NS", ifelse(fdr_vec < 0.01, "FDR < 0.01",
                             ifelse(fdr_vec < 0.05, "FDR < 0.05", "NS"))),
                      row.names = display_families)
colnames(ann_row) <- ASSOC_LABEL

ann_colors <- setNames(list(c(Heartwood = "#5D4037", Sapwood = "#BCAAA4"),
                            colorRampPalette(c("#f0f0f0", "#7b9fc2", "#1a3a5c"))(100),
                            c("FDR < 0.01" = "#b71c1c", "FDR < 0.05" = "#ef9a9a", "NS" = "#f0f0f0")),
                       c("Compartment", LOG_LABEL, ASSOC_LABEL))
heatmap_colors <- colorRampPalette(c("#3b4cc0", "#7092d5", "#c5d5ea", "#f7f7f7",
                                     "#f5c5a3", "#e8845a", "#b40426"))(100)   # diverging, row z-scores

ordered_families <- pvals %>% filter(family %in% display_families) %>% arrange(t) %>% pull(family)
heatmap_data <- heatmap_data[, ordered_families]
ann_row <- ann_row[ordered_families, , drop = FALSE]

png(OUT_PNG, width = 13, height = max(7, length(display_families) * 0.35 + 3), units = "in", res = 300)
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
cat("\nHeatmap saved:", OUT_PNG, "\n")
