# ==============================================================================
# Methanotroph Community Ordination & Clustering
# ==============================================================================
# Purpose: Dimensional reduction and clustering of methanotroph communities
#   overlaid with ddPCR pmoA/mmoX gene abundances. Tests whether community
#   structure maps onto gene abundance patterns.
#
# Analyses:
#   A. RDA biplot (constrained ordination with pmoA/mmoX as explanatory)
#   B. NMDS (unconstrained, colored by gene abundance)
#   C. PERMANOVA (community ~ pmoA + mmoX + compartment)
#   D. Hierarchical clustering with gene abundance comparison
#
# Inputs: Same as 01_methanotroph_gene_taxonomy.R (runs same phyloseq build)
#
# Outputs:
#   - fig_methanotroph_ordination.png (multi-panel)
#   - permanova_results.txt
#
# Required packages: phyloseq, tidyverse, vegan, patchwork, ggrepel
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)
library(ggrepel)

# ==============================================================================
# STEP 1: Build phyloseq object (identical to Script 01)
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

pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
no_mito <- pop_taxa(raw_ps, c(mitochondria, chloroplast))

taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)
colnames(tax_table(ps.rare)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# ==============================================================================
# STEP 2: Filter and assign compartments
# ==============================================================================

samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)

well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                   "FAGR", "FRAM", "PIST", "QURU", "TSCA")
samp_df <- samp_df %>%
  filter(!is.na(compartment), !is.na(species.x), species.x %in% well_sampled,
         !is.na(pmoa_loose), !is.na(mmox_loose))

samp_df$log_pmoa <- log10(samp_df$pmoa_loose + 1)
samp_df$log_mmox <- log10(samp_df$mmox_loose + 1)
samp_df$pct_pmoa <- ifelse(
  samp_df$pmoa_loose + samp_df$mmox_loose > 0,
  samp_df$pmoa_loose / (samp_df$pmoa_loose + samp_df$mmox_loose) * 100,
  NA
)

sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(sample_names(ps.ra) %in% rownames(samp_df), ps.ra)

# ==============================================================================
# STEP 3: Subset to methanotroph ASVs only
# ==============================================================================

source("code/00_harmonization/load_methanotroph_definitions.R")
mt_defs <- load_methanotroph_defs()
tax_df <- as.data.frame(tax_table(ps.filt))
tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs, include_conditional = TRUE)

mt_asvs <- rownames(tax_df)[!is.na(tax_df$mt_status)]
cat("Methanotroph ASVs in filtered dataset:", length(mt_asvs), "\n")

# Subset phyloseq to methanotroph ASVs
ps.mt <- prune_taxa(mt_asvs, ps.filt)

# Remove samples with zero methanotroph reads
ps.mt <- prune_samples(sample_sums(ps.mt) > 0, ps.mt)
cat("Samples with methanotroph reads:", nsamples(ps.mt), "\n")

compartment_colors <- c(
  "Heartwood" = "#a6611a", "Sapwood" = "#dfc27d",
  "Mineral Soil" = "#80cdc1", "Organic Soil" = "#018571"
)

# ==============================================================================
# STEP 4A: RDA biplot (constrained ordination)
# ==============================================================================

# Build family-level abundance matrix for RDA
mt_tax <- as.data.frame(tax_table(ps.mt))
mt_otu <- as.data.frame(t(otu_table(ps.mt)))  # samples x ASVs
mt_meta <- data.frame(sample_data(ps.mt))

# Aggregate to family level
families <- mt_tax$Family
families[is.na(families) | families == ""] <- "Unclassified"
names(families) <- rownames(mt_tax)

unique_fams <- unique(families)
fam_mat <- matrix(0, nrow = nrow(mt_otu), ncol = length(unique_fams))
rownames(fam_mat) <- rownames(mt_otu)
colnames(fam_mat) <- unique_fams

for (i in seq_len(ncol(mt_otu))) {
  fam <- families[colnames(mt_otu)[i]]
  fam_mat[, fam] <- fam_mat[, fam] + mt_otu[, i]
}

# Remove families with near-zero variance
fam_mat <- fam_mat[, colSums(fam_mat) > 0, drop = FALSE]

# Match samples
shared <- intersect(rownames(fam_mat), rownames(mt_meta))
fam_mat <- fam_mat[shared, ]
mt_meta <- mt_meta[shared, ]

# Hellinger transform for RDA
fam_hell <- decostand(fam_mat, method = "hellinger")

# Run RDA with pmoA and mmoX as constraints
rda_result <- rda(fam_hell ~ log_pmoa + log_mmox, data = mt_meta)

# Significance test
rda_anova <- anova.cca(rda_result, permutations = 999)
cat("\n=== RDA ANOVA ===\n")
print(rda_anova)

rda_axis_anova <- anova.cca(rda_result, by = "axis", permutations = 999)
cat("\n=== RDA Axis ANOVA ===\n")
print(rda_axis_anova)

# Extract scores for plotting
site_scores <- as.data.frame(scores(rda_result, display = "sites", choices = 1:2))
site_scores$compartment <- mt_meta$compartment[match(rownames(site_scores), rownames(mt_meta))]

species_scores <- as.data.frame(scores(rda_result, display = "species", choices = 1:2))
species_scores$family <- rownames(species_scores)

biplot_scores <- as.data.frame(scores(rda_result, display = "bp", choices = 1:2))
biplot_scores$variable <- rownames(biplot_scores)

# Variance explained
eig <- rda_result$CCA$eig
prop_explained <- eig / sum(rda_result$CCA$eig + rda_result$CA$eig) * 100

p_rda <- ggplot() +
  geom_point(data = site_scores,
             aes(x = RDA1, y = RDA2, color = compartment),
             alpha = 0.5, size = 1.5) +
  geom_text_repel(data = species_scores,
                   aes(x = RDA1, y = RDA2, label = family),
                   size = 3, fontface = "italic", max.overlaps = 15) +
  geom_segment(data = biplot_scores,
               aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text(data = biplot_scores,
            aes(x = RDA1 * 2.2, y = RDA2 * 2.2, label = variable),
            size = 3.5, fontface = "bold") +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(
    x = sprintf("RDA1 (%.1f%%)", prop_explained[1]),
    y = sprintf("RDA2 (%.1f%%)", prop_explained[2]),
    title = "RDA: Methanotroph families constrained by pmoA & mmoX"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ==============================================================================
# STEP 4B: NMDS (unconstrained)
# ==============================================================================

# Bray-Curtis on methanotroph ASV-level data
set.seed(42)
nmds_result <- ordinate(ps.mt, method = "NMDS", distance = "bray", trymax = 100)
cat("\nNMDS Stress:", nmds_result$stress, "\n")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$compartment <- mt_meta$compartment[match(rownames(nmds_scores), rownames(mt_meta))]
nmds_scores$log_pmoa <- mt_meta$log_pmoa[match(rownames(nmds_scores), rownames(mt_meta))]
nmds_scores$log_mmox <- mt_meta$log_mmox[match(rownames(nmds_scores), rownames(mt_meta))]
nmds_scores$pct_pmoa <- mt_meta$pct_pmoa[match(rownames(nmds_scores), rownames(mt_meta))]

# NMDS colored by compartment
p_nmds_comp <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = compartment)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(title = sprintf("NMDS (stress = %.3f)", nmds_result$stress)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# NMDS colored by pmoA
p_nmds_pmoa <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = log_pmoa)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis_c(name = expression(log[10]*"(pmoA+1)"), option = "C") +
  labs(title = "NMDS colored by pmoA abundance") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# NMDS colored by mmoX
p_nmds_mmox <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = log_mmox)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis_c(name = expression(log[10]*"(mmoX+1)"), option = "D") +
  labs(title = "NMDS colored by mmoX abundance") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# NMDS colored by %pmoA
p_nmds_pct <- ggplot(nmds_scores %>% filter(!is.na(pct_pmoa)),
                      aes(x = NMDS1, y = NMDS2, color = pct_pmoa)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_viridis_c(name = "% pmoA", option = "B") +
  labs(title = "NMDS colored by pmoA/(pmoA+mmoX)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ==============================================================================
# STEP 4C: PERMANOVA
# ==============================================================================

bc_dist <- vegdist(fam_mat[shared, ], method = "bray")

perm_result <- adonis2(
  bc_dist ~ compartment + log_pmoa + log_mmox,
  data = mt_meta, permutations = 999, method = "bray"
)

cat("\n=== PERMANOVA ===\n")
print(perm_result)

# With interaction
perm_interact <- adonis2(
  bc_dist ~ compartment * log_pmoa + compartment * log_mmox,
  data = mt_meta, permutations = 999, method = "bray"
)

cat("\n=== PERMANOVA (with interactions) ===\n")
print(perm_interact)

# Save PERMANOVA results
sink("outputs/figures/new_paper/permanova_results.txt")
cat("=== PERMANOVA: Methanotroph community ~ compartment + pmoA + mmoX ===\n\n")
print(perm_result)
cat("\n=== PERMANOVA with interactions ===\n\n")
print(perm_interact)
cat("\n=== RDA ANOVA ===\n\n")
print(rda_anova)
cat("\n=== RDA Axis ANOVA ===\n\n")
print(rda_axis_anova)
sink()

# ==============================================================================
# STEP 4D: Hierarchical clustering
# ==============================================================================

hc <- hclust(bc_dist, method = "ward.D2")

# Determine k by silhouette
library(cluster)
sil_widths <- sapply(2:8, function(k) {
  cl <- cutree(hc, k = k)
  mean(silhouette(cl, bc_dist)[, 3])
})
best_k <- which.max(sil_widths) + 1
cat("\nBest k by silhouette:", best_k, "\n")

clusters <- cutree(hc, k = best_k)
mt_meta$cluster <- factor(clusters[rownames(mt_meta)])

p_cluster_pmoa <- ggplot(mt_meta, aes(x = cluster, y = log_pmoa, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(x = "Community cluster", y = expression(log[10]*"(pmoA + 1)"),
       title = "pmoA by methanotroph community cluster") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

p_cluster_mmox <- ggplot(mt_meta, aes(x = cluster, y = log_mmox, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  labs(x = "Community cluster", y = expression(log[10]*"(mmoX + 1)"),
       title = "mmoX by methanotroph community cluster") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# Cluster composition (which compartments in each cluster)
p_cluster_comp <- ggplot(mt_meta, aes(x = cluster, fill = compartment)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = compartment_colors, name = "Compartment") +
  labs(x = "Community cluster", y = "Proportion",
       title = "Compartment composition by cluster") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# ==============================================================================
# STEP 5: Combine and save
# ==============================================================================

# Top row: RDA + NMDS by compartment
# Middle row: NMDS by pmoA, mmoX, %pmoA
# Bottom row: Clustering

fig_top <- p_rda + p_nmds_comp + plot_layout(widths = c(1, 1))
fig_mid <- p_nmds_pmoa + p_nmds_mmox + p_nmds_pct + plot_layout(widths = c(1, 1, 1))
fig_bot <- p_cluster_pmoa + p_cluster_mmox + p_cluster_comp + plot_layout(widths = c(1, 1, 1))

fig_ordination <- fig_top / fig_mid / fig_bot +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_methanotroph_ordination.png",
       fig_ordination, width = 16, height = 18, dpi = 300)
cat("\nSaved: fig_methanotroph_ordination.png\n")
cat("Saved: permanova_results.txt\n")
cat("Done.\n")
