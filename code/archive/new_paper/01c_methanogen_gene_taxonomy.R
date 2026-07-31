# ==============================================================================
# Methanogen ddPCR Gene–Taxonomy Correlations (mcrA)
# ==============================================================================
# Purpose: Systematically test how ddPCR mcrA gene abundance relates to
#   16S-inferred methanogen taxonomy across compartments and species.
#   Parallel analysis to 01_methanotroph_gene_taxonomy.R but for methanogens.
#
# Key comparison: Relative abundance (% of 16S community) vs Estimated absolute
#   abundance (relative abundance × total 16S qPCR copies/µL).
#
# Methanogen families: Methanobacteriaceae, Methanomassiliicoccaceae,
#   Methanoregulaceae, Methanocellaceae, Methanosaetaceae, Methanomicrobiaceae,
#   Methanosarcinaceae, Methanomethyliaceae, Methanocorpusculaceae
#
# Metabolic groupings:
#   Hydrogenotrophic: Methanobacteriaceae, Methanoregulaceae,
#                     Methanomicrobiaceae, Methanocorpusculaceae
#   Acetoclastic:     Methanosaetaceae, Methanosarcinaceae (also versatile)
#   Methylotrophic:   Methanomethyliaceae, Methanomassiliicoccaceae (methyl-reducing)
#
# Inputs:
#   - OTU_table.txt, unrooted_tree.nwk (16S data)
#   - tree_16s_mapping_dada2_corrected.txt (16S sample mapping)
#   - ddPCR_meta_all_data.csv (ddPCR with mcrA)
#   - 16S_tree_sample_table_with_meta.csv (16S qPCR absolute abundance)
#   - Tree_Core_Sectioning_Data.csv, 16s_w_metadata.csv (metadata)
#
# Outputs:
#   - fig_mcra_correlation_heatmap_paired.png (rel + abs heatmaps)
#   - fig_mcra_vs_taxonomy_paired.png (scatter plots, rel vs abs)
#   - fig_mcra_vs_taxonomy_combined.png (combined scatter, absolute)
#   - fig_mcra_vs_taxonomy_species.png (species-level)
#
# Required packages: phyloseq, tidyverse, patchwork, ggpubr
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)

# ==============================================================================
# STEP 1: Build phyloseq object (same as methanotroph script)
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

# Remove mitochondria and chloroplasts
pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
badTaxa <- c(mitochondria, chloroplast)
no_mito <- pop_taxa(raw_ps, badTaxa)

taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)
colnames(tax_table(ps.rare)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# ==============================================================================
# STEP 2: Assign compartments, filter samples
# ==============================================================================

samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)

samp_df <- samp_df %>%
  filter(!is.na(compartment), !is.na(species.x), species.x != "")

well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                   "FAGR", "FRAM", "PIST", "QURU", "TSCA")
samp_df <- samp_df %>% filter(species.x %in% well_sampled)

sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(sample_names(ps.ra) %in% rownames(samp_df), ps.ra)

cat("Filtered samples:", nsamples(ps.filt), "\n")

# ==============================================================================
# STEP 3: Define methanogen families and taxonomic subsets
# ==============================================================================

otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))

# Core methanogen families
methanogen_families <- c(
  "Methanobacteriaceae",
  "Methanomassiliicoccaceae",
  "Methanoregulaceae",
  "Methanocellaceae",
  "Methanosaetaceae",
  "Methanomicrobiaceae",
  "Methanosarcinaceae",
  "Methanomethyliaceae",
  "Methanocorpusculaceae"
)

# Identify methanogen ASVs
all_mg_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]
cat("Methanogen ASVs:", length(all_mg_asvs), "\n")

# Helper functions
sum_asvs <- function(asv_set, otu_mat) {
  if (length(asv_set) == 0) return(rep(0, ncol(otu_mat)))
  colSums(otu_mat[asv_set, , drop = FALSE], na.rm = TRUE)
}

family_asvs <- function(fam_name, tax_df) {
  rownames(tax_df)[which(tax_df$Family == fam_name)]
}

# --- Build taxonomy subsets ---
subsets <- list()

# All methanogens
subsets[["All methanogens"]] <- all_mg_asvs

# By family
for (fam in methanogen_families) {
  fam_asv <- family_asvs(fam, tax_df)
  if (length(fam_asv) > 0) subsets[[fam]] <- fam_asv
}

# By metabolic type
hydrogenotrophic_fams <- c("Methanobacteriaceae", "Methanoregulaceae",
                            "Methanomicrobiaceae", "Methanocorpusculaceae")
acetoclastic_fams <- c("Methanosaetaceae", "Methanosarcinaceae")
methylotrophic_fams <- c("Methanomethyliaceae", "Methanomassiliicoccaceae")

subsets[["Hydrogenotrophic"]] <- unlist(lapply(hydrogenotrophic_fams, family_asvs, tax_df = tax_df))
subsets[["Acetoclastic"]] <- unlist(lapply(acetoclastic_fams, family_asvs, tax_df = tax_df))
subsets[["Methylotrophic/methyl-reducing"]] <- unlist(lapply(methylotrophic_fams, family_asvs, tax_df = tax_df))

# Report subset sizes
cat("\n=== Methanogen taxonomic subset ASV counts ===\n")
for (nm in names(subsets)) {
  cat(sprintf("  %-40s %d ASVs\n", nm, length(subsets[[nm]])))
}

# ==============================================================================
# STEP 4: Compute per-sample abundances (relative and absolute)
# ==============================================================================

samp_meta <- data.frame(sample_data(ps.filt), stringsAsFactors = FALSE)

# Relative abundances
for (nm in names(subsets)) {
  col_name <- paste0("tx_", gsub("[^A-Za-z0-9]", "_", nm))
  samp_meta[[col_name]] <- sum_asvs(subsets[[nm]], otu_df)[rownames(samp_meta)]
}

# Load 16S qPCR for absolute abundance
picrust_meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                          row.names = 1)
samp_meta$X16S_per_ul <- picrust_meta[rownames(samp_meta), "X16S_per_ul"]

cat("\n=== 16S qPCR coverage ===\n")
cat("Samples with X16S_per_ul:", sum(!is.na(samp_meta$X16S_per_ul)), "/", nrow(samp_meta), "\n")

# Compute absolute abundance
tx_cols <- grep("^tx_", colnames(samp_meta), value = TRUE)
for (tc in tx_cols) {
  abs_col <- sub("^tx_", "abs_", tc)
  samp_meta[[abs_col]] <- samp_meta[[tc]] * samp_meta$X16S_per_ul / 100
}

# Log-transform absolute abundances
abs_cols <- grep("^abs_", colnames(samp_meta), value = TRUE)
for (ac in abs_cols) {
  log_col <- paste0("log_", ac)
  samp_meta[[log_col]] <- log10(samp_meta[[ac]] + 1)
}

# Add mcrA gene metrics
samp_meta$mcra <- samp_meta$mcra_probe_loose
samp_meta$log_mcra <- log10(samp_meta$mcra + 1)

# Filter to samples with ddPCR mcrA data
plot_data <- samp_meta %>% filter(!is.na(mcra))
cat("\nSamples with mcrA ddPCR and 16S:", nrow(plot_data), "\n")
cat("  Of those, with 16S qPCR:", sum(!is.na(plot_data$X16S_per_ul)), "\n")
cat("By compartment:\n")
print(table(plot_data$compartment))

# ==============================================================================
# STEP 5: Species mapping & shared aesthetics
# ==============================================================================

species_mapping <- c(
  "ACRU" = "A. rubrum", "ACSA" = "A. saccharum",
  "BEAL" = "B. alleghaniensis", "BELE" = "B. lenta",
  "BEPA" = "B. papyrifera", "FAGR" = "F. grandifolia",
  "FRAM" = "F. americana", "PIST" = "P. strobus",
  "QURU" = "Q. rubra", "TSCA" = "T. canadensis"
)
plot_data$species_label <- species_mapping[plot_data$species.x]
plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

compartment_colors <- c(
  "Heartwood" = "#a6611a", "Sapwood" = "#dfc27d",
  "Mineral Soil" = "#80cdc1", "Organic Soil" = "#018571"
)

# ==============================================================================
# STEP 6: Compute correlation matrices — BOTH relative and absolute
# ==============================================================================

tx_cols <- grep("^tx_", colnames(plot_data), value = TRUE)
tx_labels <- gsub("^tx_", "", tx_cols)
tx_labels <- gsub("_", " ", tx_labels)

log_abs_cols <- paste0("log_abs_", gsub("^tx_", "", tx_cols))
stopifnot(all(log_abs_cols %in% colnames(plot_data)))

# For methanogens, single gene metric: mcrA
gene_metrics <- c("log_mcra")
gene_labels <- c("mcrA")

compartments <- levels(plot_data$compartment)

# Generic correlation runner
run_correlations <- function(data, y_cols, gene_metrics, compartments, labels) {
  results <- list()
  for (comp in compartments) {
    comp_data <- data %>% filter(compartment == comp)
    rho_mat <- matrix(NA, nrow = length(y_cols), ncol = length(gene_metrics))
    p_mat   <- matrix(NA, nrow = length(y_cols), ncol = length(gene_metrics))
    for (i in seq_along(y_cols)) {
      for (j in seq_along(gene_metrics)) {
        x <- comp_data[[y_cols[i]]]
        y <- comp_data[[gene_metrics[j]]]
        valid <- complete.cases(x, y)
        if (sum(valid) >= 5) {
          ct <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
          rho_mat[i, j] <- ct$estimate
          p_mat[i, j] <- ct$p.value
        }
      }
    }
    rownames(rho_mat) <- labels
    colnames(rho_mat) <- gene_labels
    rownames(p_mat) <- labels
    colnames(p_mat) <- gene_labels
    results[[comp]] <- list(rho = rho_mat, p = p_mat, n = nrow(comp_data))
  }
  return(results)
}

cor_rel <- run_correlations(plot_data, tx_cols, gene_metrics, compartments, tx_labels)
cor_abs <- run_correlations(plot_data, log_abs_cols, gene_metrics, compartments, tx_labels)

# ==============================================================================
# STEP 7: Correlation heatmaps — paired (relative vs absolute)
# ==============================================================================

cor_to_long <- function(cor_results, compartments, tx_labels, gene_labels, type_label) {
  cor_long <- data.frame()
  for (comp in compartments) {
    rho <- cor_results[[comp]]$rho
    pval <- cor_results[[comp]]$p
    for (i in seq_len(nrow(rho))) {
      for (j in seq_len(ncol(rho))) {
        sig <- ""
        if (!is.na(pval[i, j])) {
          if (pval[i, j] < 0.001) sig <- "***"
          else if (pval[i, j] < 0.01) sig <- "**"
          else if (pval[i, j] < 0.05) sig <- "*"
        }
        cor_long <- rbind(cor_long, data.frame(
          compartment = comp,
          taxon_subset = rownames(rho)[i],
          gene_metric = colnames(rho)[j],
          rho = rho[i, j],
          p = pval[i, j],
          sig = sig,
          abundance_type = type_label,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  cor_long$compartment <- factor(cor_long$compartment, levels = compartments)
  cor_long$taxon_subset <- factor(cor_long$taxon_subset, levels = rev(tx_labels))
  cor_long$gene_metric <- factor(cor_long$gene_metric, levels = gene_labels)
  return(cor_long)
}

cor_long_rel <- cor_to_long(cor_rel, compartments, tx_labels, gene_labels, "Relative abundance (%)")
cor_long_abs <- cor_to_long(cor_abs, compartments, tx_labels, gene_labels, "Estimated absolute abundance")

make_heatmap <- function(cor_df, title = "") {
  ggplot(cor_df, aes(x = gene_metric, y = taxon_subset, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig), size = 3.5, vjust = 0.5) +
    facet_wrap(~ compartment, nrow = 1) +
    scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426",
                          midpoint = 0, limits = c(-1, 1),
                          name = expression(rho)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold")
    )
}

p_heatmap_rel <- make_heatmap(cor_long_rel, "mcrA vs methanogen taxonomy — Relative abundance")
p_heatmap_abs <- make_heatmap(cor_long_abs, "mcrA vs methanogen taxonomy — Estimated absolute abundance")

ggsave("outputs/figures/new_paper/fig_mcra_correlation_heatmap_relative.png",
       p_heatmap_rel, width = 10, height = 8, dpi = 300)
ggsave("outputs/figures/new_paper/fig_mcra_correlation_heatmap_absolute.png",
       p_heatmap_abs, width = 10, height = 8, dpi = 300)

fig_paired <- p_heatmap_rel / p_heatmap_abs +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_mcra_correlation_heatmap_paired.png",
       fig_paired, width = 10, height = 15, dpi = 300)
cat("Saved: fig_mcra_correlation_heatmap_paired.png\n")

# ==============================================================================
# STEP 8: Key scatter plots — relative and absolute side by side
# ==============================================================================

make_scatter <- function(data, x_col, y_col, x_lab, y_lab, title = "") {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = compartment)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
    stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
    facet_wrap(~ compartment, nrow = 1, scales = "free") +
    scale_color_manual(values = compartment_colors) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 9, face = "bold"),
      strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.title = element_text(size = 10, face = "bold")
    )
}

# Paired scatter: mcrA vs All methanogens
s1_rel <- make_scatter(plot_data, "log_mcra", "tx_All_methanogens",
                        expression(log[10]*"(mcrA + 1)"),
                        "Rel. abundance (%)",
                        "mcrA vs All methanogens — Relative")
s1_abs <- make_scatter(plot_data, "log_mcra", "log_abs_All_methanogens",
                        expression(log[10]*"(mcrA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mcrA vs All methanogens — Absolute")

# Paired scatter: mcrA vs Hydrogenotrophic
s2_rel <- make_scatter(plot_data, "log_mcra", "tx_Hydrogenotrophic",
                        expression(log[10]*"(mcrA + 1)"),
                        "Rel. abundance (%)",
                        "mcrA vs Hydrogenotrophic — Relative")
s2_abs <- make_scatter(plot_data, "log_mcra", "log_abs_Hydrogenotrophic",
                        expression(log[10]*"(mcrA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mcrA vs Hydrogenotrophic — Absolute")

# Paired scatter: mcrA vs Acetoclastic
s3_rel <- make_scatter(plot_data, "log_mcra", "tx_Acetoclastic",
                        expression(log[10]*"(mcrA + 1)"),
                        "Rel. abundance (%)",
                        "mcrA vs Acetoclastic — Relative")
s3_abs <- make_scatter(plot_data, "log_mcra", "log_abs_Acetoclastic",
                        expression(log[10]*"(mcrA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mcrA vs Acetoclastic — Absolute")

# Paired scatter: mcrA vs Methylotrophic
s4_rel <- make_scatter(plot_data, "log_mcra", "tx_Methylotrophic_methyl_reducing",
                        expression(log[10]*"(mcrA + 1)"),
                        "Rel. abundance (%)",
                        "mcrA vs Methylotrophic — Relative")
s4_abs <- make_scatter(plot_data, "log_mcra", "log_abs_Methylotrophic_methyl_reducing",
                        expression(log[10]*"(mcrA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mcrA vs Methylotrophic — Absolute")

fig_scatter_paired <- (s1_rel / s1_abs / s2_rel / s2_abs /
                        s3_rel / s3_abs / s4_rel / s4_abs) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_mcra_vs_taxonomy_paired.png",
       fig_scatter_paired, width = 14, height = 32, dpi = 300)
cat("Saved: fig_mcra_vs_taxonomy_paired.png\n")

# ==============================================================================
# STEP 9: Combined scatter (absolute abundance, single panels)
# ==============================================================================

make_scatter_combined <- function(data, x_col, y_col, x_lab, y_lab, title = "") {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = compartment)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
    scale_color_manual(values = compartment_colors, name = "Compartment") +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 11, face = "bold")
    )
}

ac1 <- make_scatter_combined(plot_data, "log_mcra", "log_abs_All_methanogens",
                              expression(log[10]*"(mcrA + 1)"),
                              expression(log[10]*"(est. absolute methanogen abund. + 1)"),
                              "mcrA vs All methanogens (absolute)")
ac2 <- make_scatter_combined(plot_data, "log_mcra", "log_abs_Hydrogenotrophic",
                              expression(log[10]*"(mcrA + 1)"),
                              expression(log[10]*"(est. absolute hydrogenotrophic + 1)"),
                              "mcrA vs Hydrogenotrophic (absolute)")
ac3 <- make_scatter_combined(plot_data, "log_mcra", "log_abs_Acetoclastic",
                              expression(log[10]*"(mcrA + 1)"),
                              expression(log[10]*"(est. absolute acetoclastic + 1)"),
                              "mcrA vs Acetoclastic (absolute)")
ac4 <- make_scatter_combined(plot_data, "log_mcra", "log_abs_Methylotrophic_methyl_reducing",
                              expression(log[10]*"(mcrA + 1)"),
                              expression(log[10]*"(est. absolute methylotrophic + 1)"),
                              "mcrA vs Methylotrophic (absolute)")

fig_abs_combined <- (ac1 + ac2) / (ac3 + ac4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_mcra_vs_taxonomy_combined.png",
       fig_abs_combined, width = 12, height = 10, dpi = 300)
cat("Saved: fig_mcra_vs_taxonomy_combined.png\n")

# ==============================================================================
# STEP 10: Species-level aggregation
# ==============================================================================

species_summary <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(
    across(starts_with("tx_"), ~ mean(.x, na.rm = TRUE)),
    mean_log_abs_all = mean(log_abs_All_methanogens, na.rm = TRUE),
    mean_mcra = mean(log_mcra, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

ps1_rel <- ggplot(species_summary, aes(x = mean_mcra, y = tx_All_methanogens,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(mcrA + 1)"),
       y = "Mean methanogen rel. abund. (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ps1_abs <- ggplot(species_summary, aes(x = mean_mcra, y = mean_log_abs_all,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(mcrA + 1)"),
       y = expression("Mean "*log[10]*"(est. abs. abund. + 1)")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

fig_species <- ps1_rel / ps1_abs +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_mcra_vs_taxonomy_species.png",
       fig_species, width = 8, height = 10, dpi = 300)
cat("Saved: fig_mcra_vs_taxonomy_species.png\n")

# ==============================================================================
# STEP 11: Save correlation tables
# ==============================================================================

write.csv(cor_long_rel, "outputs/figures/new_paper/mcra_correlation_results_relative.csv", row.names = FALSE)
write.csv(cor_long_abs, "outputs/figures/new_paper/mcra_correlation_results_absolute.csv", row.names = FALSE)
cat("\nAll mcrA-methanogen gene-taxonomy figures complete.\n")
