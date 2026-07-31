# ==============================================================================
# Methanotroph ddPCR Gene–Taxonomy Correlations
# ==============================================================================
# Purpose: Systematically test how ddPCR pmoA and mmoX gene abundances relate
#   to 16S-inferred methanotroph taxonomy across compartments and species.
#   Tests correlations against multiple taxonomic subsets (by definition level,
#   family, ecological grouping, MMO type).
#
# Key comparison: Relative abundance (% of 16S community) vs Estimated absolute
#   abundance (relative abundance × total 16S qPCR copies/µL). Soil has ~100–
#   800× more total bacteria than wood, so relative abundance can produce
#   misleading negative correlations when comparing to absolute ddPCR targets.
#
# Key biological question: Do pmoA (pMMO, obligate methanotrophs) and mmoX
#   (sMMO, facultative methanotrophs) track with the expected taxonomic groups?
#
# Note on taxonomic resolution: Without pmoA amplicon sequencing, we cannot
#   resolve Upland Soil Clusters (USCα, USCγ) directly. Putative methanotroph
#   assignments are inferred from 16S taxonomy using Knief (2015) definitions.
#
# Inputs:
#   - OTU_table.txt, taxonomy_table.txt, unrooted_tree.nwk (16S data)
#   - tree_16s_mapping_dada2_corrected.txt (16S sample mapping)
#   - ddPCR_meta_all_data.csv (ddPCR with pmoA/mmoX)
#   - 16S_tree_sample_table_with_meta.csv (16S qPCR absolute abundance)
#   - Tree_Core_Sectioning_Data.csv, 16s_w_metadata.csv (metadata)
#   - methanotroph_definitions.csv (Knief 2015)
#
# Outputs:
#   - fig_correlation_heatmap_relative.png (relative abundance correlations)
#   - fig_correlation_heatmap_absolute.png (absolute abundance correlations)
#   - fig_correlation_heatmap_paired.png (both side by side)
#   - fig_gene_vs_taxonomy_faceted.png (scatter plots, faceted by compartment)
#   - fig_gene_vs_taxonomy_combined.png (scatter plots, colored by compartment)
#   - fig_biomass_context.png (16S qPCR by compartment — context panel)
#
# Required packages: phyloseq, tidyverse, patchwork, vegan, ggpubr
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)

# ==============================================================================
# STEP 1: Build phyloseq object (from 08b_methanotroph_16s_composition.R)
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

# Filter to samples with compartment, species, and ddPCR data
samp_df <- samp_df %>%
  filter(!is.na(compartment), !is.na(species.x), species.x != "")

# Exclude under-sampled species
well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                   "FAGR", "FRAM", "PIST", "QURU", "TSCA")
samp_df <- samp_df %>% filter(species.x %in% well_sampled)

sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(sample_names(ps.ra) %in% rownames(samp_df), ps.ra)

cat("Filtered samples:", nsamples(ps.filt), "\n")

# ==============================================================================
# STEP 3: Classify methanotroph ASVs and compute per-sample abundances
# ==============================================================================

source("code/00_harmonization/load_methanotroph_definitions.R")
mt_defs <- load_methanotroph_defs()

otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))

tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs, include_conditional = TRUE)

known_asvs <- rownames(tax_df)[tax_df$mt_status == "Known" & !is.na(tax_df$mt_status)]
putative_asvs <- rownames(tax_df)[tax_df$mt_status == "Putative" & !is.na(tax_df$mt_status)]
all_mt_asvs <- c(known_asvs, putative_asvs)

cat("Known methanotroph ASVs:", length(known_asvs), "\n")
cat("Putative methanotroph ASVs:", length(putative_asvs), "\n")

# ==============================================================================
# STEP 4: Define taxonomic subsets and compute per-sample abundances
# ==============================================================================

# Helper: sum relative abundance for a set of ASVs per sample
sum_asvs <- function(asv_set, otu_mat) {
  if (length(asv_set) == 0) return(rep(0, ncol(otu_mat)))
  colSums(otu_mat[asv_set, , drop = FALSE], na.rm = TRUE)
}

# Helper: get ASVs matching a family
family_asvs <- function(fam_name, tax_df, mt_asvs = NULL) {
  idx <- which(tax_df$Family == fam_name)
  asvs <- rownames(tax_df)[idx]
  if (!is.null(mt_asvs)) asvs <- intersect(asvs, mt_asvs)
  return(asvs)
}

# Helper: get ASVs matching genus names
genus_asvs <- function(genus_names, tax_df) {
  idx <- which(tax_df$Genus %in% genus_names)
  rownames(tax_df)[idx]
}

# --- Build taxonomy subsets ---
subsets <- list()

# By definition level
subsets[["All methanotrophs"]] <- all_mt_asvs
subsets[["Known only"]] <- known_asvs
subsets[["Putative only"]] <- putative_asvs

# By family (individual) — restrict to methanotroph ASVs within each family
subsets[["Methylococcaceae"]] <- family_asvs("Methylococcaceae", tax_df, all_mt_asvs)
subsets[["Methylomonadaceae"]] <- family_asvs("Methylomonadaceae", tax_df, all_mt_asvs)
subsets[["Methylothermaceae"]] <- family_asvs("Methylothermaceae", tax_df, all_mt_asvs)
subsets[["Methylacidiphilaceae"]] <- family_asvs("Methylacidiphilaceae", tax_df, all_mt_asvs)
subsets[["Methylocystaceae"]] <- family_asvs("Methylocystaceae", tax_df, all_mt_asvs)
subsets[["Beijerinckiaceae"]] <- family_asvs("Beijerinckiaceae", tax_df, all_mt_asvs)

# By ecological grouping
type_i_fams <- c("Methylococcaceae", "Methylothermaceae", "Methylomonadaceae")
subsets[["Type I Gamma"]] <- unlist(lapply(type_i_fams, family_asvs, tax_df = tax_df))

type_ii_fams <- c("Methylocystaceae", "Beijerinckiaceae")
type_ii_genera <- c("Methylocystis", "Methylosinus", "Methylocapsa", "Methylocella", "Methyloferula")
subsets[["Type II Alpha"]] <- union(
  unlist(lapply(type_ii_fams, family_asvs, tax_df = tax_df, mt_asvs = all_mt_asvs)),
  genus_asvs(type_ii_genera, tax_df)
)

pmmo_genera <- c("Methylocystis", "Methylosinus")
subsets[["pMMO-expected"]] <- union(
  subsets[["Type I Gamma"]],
  union(family_asvs("Methylacidiphilaceae", tax_df),
        genus_asvs(pmmo_genera, tax_df))
)

smmo_only_genera <- c("Methylocella", "Methyloferula")
subsets[["sMMO-only (Methylocella/Methyloferula)"]] <- genus_asvs(smmo_only_genera, tax_df)

# Report subset sizes
cat("\n=== Taxonomic subset ASV counts ===\n")
for (nm in names(subsets)) {
  cat(sprintf("  %-40s %d ASVs\n", nm, length(subsets[[nm]])))
}

# --- Compute per-sample RELATIVE abundances for each subset ---
samp_meta <- data.frame(sample_data(ps.filt), stringsAsFactors = FALSE)

for (nm in names(subsets)) {
  col_name <- paste0("tx_", gsub("[^A-Za-z0-9]", "_", nm))
  samp_meta[[col_name]] <- sum_asvs(subsets[[nm]], otu_df)[rownames(samp_meta)]
}

# ==============================================================================
# STEP 4b: Load 16S qPCR and compute ABSOLUTE abundances
# ==============================================================================

picrust_meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                          row.names = 1)
# Merge X16S_per_ul into samp_meta by matching row names
samp_meta$X16S_per_ul <- picrust_meta[rownames(samp_meta), "X16S_per_ul"]

cat("\n=== 16S qPCR coverage ===\n")
cat("Samples with X16S_per_ul:", sum(!is.na(samp_meta$X16S_per_ul)), "/", nrow(samp_meta), "\n")
cat("Median 16S copies/µL by compartment:\n")
print(tapply(samp_meta$X16S_per_ul, samp_meta$compartment, median, na.rm = TRUE))

# Compute absolute abundance: relative (%) × 16S copies/µL / 100
# Units: estimated 16S copies/µL attributable to each taxon subset
tx_cols <- grep("^tx_", colnames(samp_meta), value = TRUE)
for (tc in tx_cols) {
  abs_col <- sub("^tx_", "abs_", tc)
  samp_meta[[abs_col]] <- samp_meta[[tc]] * samp_meta$X16S_per_ul / 100
}

# Log-transform absolute abundances for plotting
abs_cols <- grep("^abs_", colnames(samp_meta), value = TRUE)
for (ac in abs_cols) {
  log_col <- paste0("log_", ac)
  samp_meta[[log_col]] <- log10(samp_meta[[ac]] + 1)
}

# Add gene metrics
samp_meta$pmoa <- samp_meta$pmoa_loose
samp_meta$mmox <- samp_meta$mmox_loose
samp_meta$pmoa_plus_mmox <- samp_meta$pmoa_loose + samp_meta$mmox_loose
samp_meta$pct_pmoa <- ifelse(
  samp_meta$pmoa_plus_mmox > 0,
  samp_meta$pmoa_loose / samp_meta$pmoa_plus_mmox * 100,
  NA
)

# Log-transform gene abundances
samp_meta$log_pmoa <- log10(samp_meta$pmoa + 1)
samp_meta$log_mmox <- log10(samp_meta$mmox + 1)
samp_meta$log_combined <- log10(samp_meta$pmoa_plus_mmox + 1)

# Also log-transform 16S total for context
samp_meta$log_16S_qpcr <- log10(samp_meta$X16S_per_ul + 1)

# Filter to samples with ddPCR data
plot_data <- samp_meta %>% filter(!is.na(pmoa) & !is.na(mmox))
cat("\nSamples with both ddPCR and 16S:", nrow(plot_data), "\n")
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
# STEP 5b: Biomass context figure — 16S qPCR by compartment
# ==============================================================================

biomass_df <- plot_data %>% filter(!is.na(X16S_per_ul))

p_biomass <- ggplot(biomass_df, aes(x = compartment, y = X16S_per_ul, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = compartment_colors) +
  labs(x = NULL, y = "Total bacterial 16S copies / µL (qPCR)",
       title = "Total bacterial biomass by compartment") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

ggsave("outputs/figures/new_paper/fig_biomass_context.png",
       p_biomass, width = 5, height = 4, dpi = 300)
cat("Saved: fig_biomass_context.png\n")

# ==============================================================================
# STEP 6: Compute correlation matrices — BOTH relative and absolute
# ==============================================================================

# Identify the tx_ (relative) columns and their matching abs_ (absolute) and
# log_abs_ (log-absolute) columns
tx_cols <- grep("^tx_", colnames(plot_data), value = TRUE)
tx_labels <- gsub("^tx_", "", tx_cols)
tx_labels <- gsub("_", " ", tx_labels)

# Absolute versions: log10(abs + 1)
log_abs_cols <- paste0("log_abs_", gsub("^tx_", "", tx_cols))
# Verify they exist
stopifnot(all(log_abs_cols %in% colnames(plot_data)))

gene_metrics <- c("log_pmoa", "log_mmox", "log_combined", "pct_pmoa")
gene_labels <- c("pmoA", "mmoX", "pmoA + mmoX", "% pmoA")

compartments <- levels(plot_data$compartment)

# Generic correlation runner
run_correlations <- function(data, y_cols, gene_metrics, compartments) {
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
    rownames(rho_mat) <- tx_labels
    colnames(rho_mat) <- gene_labels
    rownames(p_mat) <- tx_labels
    colnames(p_mat) <- gene_labels
    results[[comp]] <- list(rho = rho_mat, p = p_mat, n = nrow(comp_data))
  }
  return(results)
}

# Run for both relative and absolute
cor_rel <- run_correlations(plot_data, tx_cols, gene_metrics, compartments)
cor_abs <- run_correlations(plot_data, log_abs_cols, gene_metrics, compartments)

# ==============================================================================
# STEP 7: Correlation heatmaps — relative, absolute, and paired
# ==============================================================================

# Converter: correlation results → long format
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
cor_long_all <- rbind(cor_long_rel, cor_long_abs)
cor_long_all$abundance_type <- factor(cor_long_all$abundance_type,
  levels = c("Relative abundance (%)", "Estimated absolute abundance"))

# Heatmap helper
make_heatmap <- function(cor_df, title = "") {
  ggplot(cor_df, aes(x = gene_metric, y = taxon_subset, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig), size = 3, vjust = 0.5) +
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

# Individual heatmaps
p_heatmap_rel <- make_heatmap(cor_long_rel, "Relative abundance (% of 16S community)")
p_heatmap_abs <- make_heatmap(cor_long_abs, "Estimated absolute abundance (rel. abund. × 16S qPCR)")

ggsave("outputs/figures/new_paper/fig_correlation_heatmap_relative.png",
       p_heatmap_rel, width = 14, height = 8, dpi = 300)
cat("Saved: fig_correlation_heatmap_relative.png\n")

ggsave("outputs/figures/new_paper/fig_correlation_heatmap_absolute.png",
       p_heatmap_abs, width = 14, height = 8, dpi = 300)
cat("Saved: fig_correlation_heatmap_absolute.png\n")

# Paired heatmap (stacked)
fig_paired <- p_heatmap_rel / p_heatmap_abs +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_correlation_heatmap_paired.png",
       fig_paired, width = 14, height = 15, dpi = 300)
cat("Saved: fig_correlation_heatmap_paired.png\n")

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

# Paired scatter: pmoA vs All methanotrophs (relative then absolute)
s1_rel <- make_scatter(plot_data, "log_pmoa", "tx_All_methanotrophs",
                        expression(log[10]*"(pmoA + 1)"),
                        "Rel. abundance (%)",
                        "pmoA vs All methanotrophs — Relative")
s1_abs <- make_scatter(plot_data, "log_pmoa", "log_abs_All_methanotrophs",
                        expression(log[10]*"(pmoA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "pmoA vs All methanotrophs — Absolute")

# Paired scatter: mmoX vs All methanotrophs
s2_rel <- make_scatter(plot_data, "log_mmox", "tx_All_methanotrophs",
                        expression(log[10]*"(mmoX + 1)"),
                        "Rel. abundance (%)",
                        "mmoX vs All methanotrophs — Relative")
s2_abs <- make_scatter(plot_data, "log_mmox", "log_abs_All_methanotrophs",
                        expression(log[10]*"(mmoX + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mmoX vs All methanotrophs — Absolute")

# Paired scatter: pmoA vs Methylacidiphilaceae
s3_rel <- make_scatter(plot_data, "log_pmoa", "tx_Methylacidiphilaceae",
                        expression(log[10]*"(pmoA + 1)"),
                        "Rel. abundance (%)",
                        "pmoA vs Methylacidiphilaceae — Relative")
s3_abs <- make_scatter(plot_data, "log_pmoa", "log_abs_Methylacidiphilaceae",
                        expression(log[10]*"(pmoA + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "pmoA vs Methylacidiphilaceae — Absolute")

# Paired scatter: mmoX vs Beijerinckiaceae
s4_rel <- make_scatter(plot_data, "log_mmox", "tx_Beijerinckiaceae",
                        expression(log[10]*"(mmoX + 1)"),
                        "Rel. abundance (%)",
                        "mmoX vs Beijerinckiaceae — Relative")
s4_abs <- make_scatter(plot_data, "log_mmox", "log_abs_Beijerinckiaceae",
                        expression(log[10]*"(mmoX + 1)"),
                        expression(log[10]*"(est. abs. abund. + 1)"),
                        "mmoX vs Beijerinckiaceae — Absolute")

fig_scatter_paired <- (s1_rel / s1_abs / s2_rel / s2_abs /
                        s3_rel / s3_abs / s4_rel / s4_abs) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_gene_vs_taxonomy_paired.png",
       fig_scatter_paired, width = 14, height = 32, dpi = 300)
cat("Saved: fig_gene_vs_taxonomy_paired.png\n")

# ==============================================================================
# STEP 9: Combined scatter (single panels, colored by compartment)
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

# Absolute abundance combined scatter — the key summary figure
ac1 <- make_scatter_combined(plot_data, "log_pmoa", "log_abs_All_methanotrophs",
                              expression(log[10]*"(pmoA + 1)"),
                              expression(log[10]*"(est. absolute methanotroph abund. + 1)"),
                              "pmoA vs All methanotrophs (absolute)")
ac2 <- make_scatter_combined(plot_data, "log_mmox", "log_abs_All_methanotrophs",
                              expression(log[10]*"(mmoX + 1)"),
                              expression(log[10]*"(est. absolute methanotroph abund. + 1)"),
                              "mmoX vs All methanotrophs (absolute)")
ac3 <- make_scatter_combined(plot_data, "log_pmoa", "log_abs_Methylacidiphilaceae",
                              expression(log[10]*"(pmoA + 1)"),
                              expression(log[10]*"(est. absolute Methylacidiphilaceae + 1)"),
                              "pmoA vs Methylacidiphilaceae (absolute)")
ac4 <- make_scatter_combined(plot_data, "log_mmox", "log_abs_Beijerinckiaceae",
                              expression(log[10]*"(mmoX + 1)"),
                              expression(log[10]*"(est. absolute Beijerinckiaceae + 1)"),
                              "mmoX vs Beijerinckiaceae (absolute)")

fig_abs_combined <- (ac1 + ac2) / (ac3 + ac4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_gene_vs_taxonomy_absolute_combined.png",
       fig_abs_combined, width = 12, height = 10, dpi = 300)
cat("Saved: fig_gene_vs_taxonomy_absolute_combined.png\n")

# ==============================================================================
# STEP 10: Species-level aggregation (both relative and absolute)
# ==============================================================================

species_summary <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(
    across(starts_with("tx_"), ~ mean(.x, na.rm = TRUE)),
    mean_log_abs_all = mean(log_abs_All_methanotrophs, na.rm = TRUE),
    mean_pmoa = mean(log_pmoa, na.rm = TRUE),
    mean_mmox = mean(log_mmox, na.rm = TRUE),
    mean_combined = mean(log_combined, na.rm = TRUE),
    mean_pct_pmoa = mean(pct_pmoa, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Relative abundance species plots
ps1_rel <- ggplot(species_summary, aes(x = mean_pmoa, y = tx_All_methanotrophs,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(pmoA + 1)"),
       y = "Mean methanotroph rel. abund. (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ps2_rel <- ggplot(species_summary, aes(x = mean_mmox, y = tx_All_methanotrophs,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(mmoX + 1)"),
       y = "Mean methanotroph rel. abund. (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

# Absolute abundance species plots
ps1_abs <- ggplot(species_summary, aes(x = mean_pmoa, y = mean_log_abs_all,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(pmoA + 1)"),
       y = expression("Mean "*log[10]*"(est. abs. abund. + 1)")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ps2_abs <- ggplot(species_summary, aes(x = mean_mmox, y = mean_log_abs_all,
                                     color = compartment, label = species_label)) +
  geom_point(size = 3) +
  geom_text(size = 2.5, hjust = -0.1, vjust = -0.3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_color_manual(values = compartment_colors, name = "Compartment") +
  labs(x = expression("Mean "*log[10]*"(mmoX + 1)"),
       y = expression("Mean "*log[10]*"(est. abs. abund. + 1)")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

fig_species <- (ps1_rel + ps2_rel) / (ps1_abs + ps2_abs) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_gene_vs_taxonomy_species.png",
       fig_species, width = 12, height = 10, dpi = 300)
cat("Saved: fig_gene_vs_taxonomy_species.png\n")

# ==============================================================================
# STEP 11: Save correlation tables for reference
# ==============================================================================

write.csv(cor_long_rel, "outputs/figures/new_paper/correlation_results_relative.csv", row.names = FALSE)
write.csv(cor_long_abs, "outputs/figures/new_paper/correlation_results_absolute.csv", row.names = FALSE)
cat("\nAll gene-taxonomy figures complete.\n")
