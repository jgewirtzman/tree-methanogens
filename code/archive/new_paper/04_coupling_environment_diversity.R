# ==============================================================================
# Methanotroph–Methanogen Coupling, Environmental Drivers, and Diversity
# ==============================================================================
# Purpose: Expanded analyses for the new paper, focusing on:
#   (A) Methanotroph–methanogen spatial coupling (pmoA/mmoX vs mcrA)
#   (B) Environmental drivers of gene abundance (VWC, ORP, soil temp)
#   (C) Gene abundance vs internal CH4 concentrations
#   (D) Methanotroph community diversity vs dominance by compartment
#
# Inputs:
#   - 16S_tree_sample_table_with_meta.csv (main metadata: ddPCR + env + gas)
#   - OTU_table.txt + phyloseq pipeline (for diversity calculations)
#   - methanotroph_definitions.csv (Knief 2015)
#
# Outputs:
#   - fig_coupling_pmoa_mcra.png
#   - fig_environmental_drivers.png
#   - fig_gene_vs_ch4_internal.png
#   - fig_methanotroph_diversity.png
#   - coupling_stats.csv
#   - environmental_driver_stats.csv
#
# Required packages: phyloseq, tidyverse, patchwork, ggpubr, lme4, vegan
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(lme4)
library(vegan)

# ==============================================================================
# STEP 1: Build phyloseq + load metadata (same as Script 01)
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

# Remove plastids
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
colnames(tax_table(ps.rare)) <- c("Kingdom", "Phylum", "Class", "Order",
                                    "Family", "Genus", "Species")
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# Load metadata with environmental variables
meta <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv",
                  row.names = 1)
meta$compartment <- case_when(
  meta$core_type == "Inner" ~ "Heartwood",
  meta$core_type == "Outer" ~ "Sapwood",
  meta$core_type == "Mineral" ~ "Mineral Soil",
  meta$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)

well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                   "FAGR", "FRAM", "PIST", "QURU", "TSCA")

meta_filt <- meta %>%
  filter(!is.na(compartment), species.x %in% well_sampled)

# Log-transform gene abundances
meta_filt$log_pmoa <- log10(meta_filt$pmoa_loose + 1)
meta_filt$log_mmox <- log10(meta_filt$mmox_loose + 1)
meta_filt$log_mcra <- log10(meta_filt$mcra_probe_loose + 1)
meta_filt$log_combined <- log10(meta_filt$pmoa_loose + meta_filt$mmox_loose + 1)
meta_filt$pct_pmoa <- ifelse(
  (meta_filt$pmoa_loose + meta_filt$mmox_loose) > 0,
  meta_filt$pmoa_loose / (meta_filt$pmoa_loose + meta_filt$mmox_loose) * 100,
  NA
)

meta_filt$compartment <- factor(meta_filt$compartment,
  levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

compartment_colors <- c(
  "Heartwood" = "#a6611a", "Sapwood" = "#dfc27d",
  "Mineral Soil" = "#80cdc1", "Organic Soil" = "#018571"
)

cat("Metadata samples:", nrow(meta_filt), "\n")
cat("By compartment:\n")
print(table(meta_filt$compartment))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ ANALYSIS A: Methanotroph–Methanogen Coupling (pmoA/mmoX vs mcrA)          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n=== ANALYSIS A: Methanotroph–Methanogen Coupling ===\n")

# Filter to samples with both methanotroph and methanogen gene data
coupling_data <- meta_filt %>%
  filter(!is.na(pmoa_loose), !is.na(mmox_loose), !is.na(mcra_probe_loose))

cat("Samples with pmoA + mmoX + mcrA:", nrow(coupling_data), "\n")

# Scatter: pmoA vs mcrA, faceted by compartment
p_coupling_1 <- ggplot(coupling_data, aes(x = log_mcra, y = log_pmoa, color = compartment)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
  facet_wrap(~ compartment, nrow = 1, scales = "free") +
  scale_color_manual(values = compartment_colors) +
  labs(x = expression(log[10]*"(mcrA + 1)"),
       y = expression(log[10]*"(pmoA + 1)"),
       title = "pmoA vs mcrA (methanotrophy vs methanogenesis)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 10, face = "bold"))

# Scatter: mmoX vs mcrA
p_coupling_2 <- ggplot(coupling_data, aes(x = log_mcra, y = log_mmox, color = compartment)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
  facet_wrap(~ compartment, nrow = 1, scales = "free") +
  scale_color_manual(values = compartment_colors) +
  labs(x = expression(log[10]*"(mcrA + 1)"),
       y = expression(log[10]*"(mmoX + 1)"),
       title = "mmoX vs mcrA") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 10, face = "bold"))

# Combined: pmoA+mmoX vs mcrA
p_coupling_3 <- ggplot(coupling_data, aes(x = log_mcra, y = log_combined, color = compartment)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
  facet_wrap(~ compartment, nrow = 1, scales = "free") +
  scale_color_manual(values = compartment_colors) +
  labs(x = expression(log[10]*"(mcrA + 1)"),
       y = expression(log[10]*"(pmoA + mmoX + 1)"),
       title = "Total methanotroph genes vs mcrA") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 10, face = "bold"))

# %pmoA vs mcrA (does methanogen abundance predict which MMO type dominates?)
p_coupling_4 <- ggplot(coupling_data %>% filter(!is.na(pct_pmoa)),
                        aes(x = log_mcra, y = pct_pmoa, color = compartment)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
  facet_wrap(~ compartment, nrow = 1, scales = "free") +
  scale_color_manual(values = compartment_colors) +
  labs(x = expression(log[10]*"(mcrA + 1)"),
       y = "% pmoA / (pmoA + mmoX)",
       title = "% pmoA vs mcrA (gene ratio vs methanogen abundance)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 10, face = "bold"))

fig_coupling <- p_coupling_1 / p_coupling_2 / p_coupling_3 / p_coupling_4 +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_coupling_pmoa_mcra.png",
       fig_coupling, width = 14, height = 16, dpi = 300)
cat("Saved: fig_coupling_pmoa_mcra.png\n")

# Coupling statistics table
coupling_stats <- data.frame()
for (comp in levels(coupling_data$compartment)) {
  cd <- coupling_data %>% filter(compartment == comp)
  for (pair in list(
    c("log_pmoa", "log_mcra", "pmoA vs mcrA"),
    c("log_mmox", "log_mcra", "mmoX vs mcrA"),
    c("log_combined", "log_mcra", "pmoA+mmoX vs mcrA")
  )) {
    valid <- complete.cases(cd[[pair[1]]], cd[[pair[2]]])
    if (sum(valid) >= 5) {
      ct <- cor.test(cd[[pair[1]]][valid], cd[[pair[2]]][valid],
                     method = "spearman", exact = FALSE)
      coupling_stats <- rbind(coupling_stats, data.frame(
        compartment = comp, comparison = pair[3],
        n = sum(valid), rho = ct$estimate, p = ct$p.value
      ))
    }
  }
}
coupling_stats$FDR <- p.adjust(coupling_stats$p, "BH")
write.csv(coupling_stats, "outputs/figures/new_paper/coupling_stats.csv", row.names = FALSE)
cat("\nCoupling statistics:\n")
print(coupling_stats, row.names = FALSE)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ ANALYSIS B: Environmental Drivers of Gene Abundance                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n=== ANALYSIS B: Environmental Drivers ===\n")

# Available: mean_vwc, mean_orp, mean_st (soil temp)
env_data <- meta_filt %>%
  filter(!is.na(pmoa_loose), !is.na(mmox_loose))

cat("Samples with environmental data:\n")
cat("  VWC:", sum(!is.na(env_data$mean_vwc)), "\n")
cat("  ORP:", sum(!is.na(env_data$mean_orp)), "\n")
cat("  Soil temp:", sum(!is.na(env_data$mean_st)), "\n")

# Helper for environmental scatter plots
make_env_scatter <- function(data, x_col, y_col, x_lab, y_lab, title = "") {
  ggplot(data %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])),
         aes(x = .data[[x_col]], y = .data[[y_col]], color = compartment)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
    stat_cor(method = "spearman", size = 2.5, label.x.npc = 0.02, label.y.npc = 0.95) +
    facet_wrap(~ compartment, nrow = 1, scales = "free") +
    scale_color_manual(values = compartment_colors) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          strip.text = element_text(size = 9, face = "bold"),
          strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          plot.title = element_text(size = 10, face = "bold"))
}

# VWC vs genes
e1 <- make_env_scatter(env_data, "mean_vwc", "log_pmoa",
                        "Mean VWC (%)", expression(log[10]*"(pmoA + 1)"),
                        "VWC vs pmoA")
e2 <- make_env_scatter(env_data, "mean_vwc", "log_mmox",
                        "Mean VWC (%)", expression(log[10]*"(mmoX + 1)"),
                        "VWC vs mmoX")

# ORP vs genes
e3 <- make_env_scatter(env_data, "mean_orp", "log_pmoa",
                        "Mean ORP (mV)", expression(log[10]*"(pmoA + 1)"),
                        "ORP vs pmoA")
e4 <- make_env_scatter(env_data, "mean_orp", "log_mmox",
                        "Mean ORP (mV)", expression(log[10]*"(mmoX + 1)"),
                        "ORP vs mmoX")

# Soil temp vs genes
e5 <- make_env_scatter(env_data, "mean_st", "log_pmoa",
                        "Mean soil temp (C)", expression(log[10]*"(pmoA + 1)"),
                        "Soil temp vs pmoA")
e6 <- make_env_scatter(env_data, "mean_st", "log_mmox",
                        "Mean soil temp (C)", expression(log[10]*"(mmoX + 1)"),
                        "Soil temp vs mmoX")

# VWC vs %pmoA (does moisture predict which MMO dominates?)
e7 <- make_env_scatter(env_data %>% filter(!is.na(pct_pmoa)),
                        "mean_vwc", "pct_pmoa",
                        "Mean VWC (%)", "% pmoA / (pmoA + mmoX)",
                        "VWC vs % pmoA")

# ORP vs %pmoA
e8 <- make_env_scatter(env_data %>% filter(!is.na(pct_pmoa)),
                        "mean_orp", "pct_pmoa",
                        "Mean ORP (mV)", "% pmoA / (pmoA + mmoX)",
                        "ORP vs % pmoA")

fig_env <- (e1 + e2) / (e3 + e4) / (e5 + e6) / (e7 + e8) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_environmental_drivers.png",
       fig_env, width = 14, height = 20, dpi = 300)
cat("Saved: fig_environmental_drivers.png\n")

# Environmental driver statistics
env_stats <- data.frame()
env_vars <- c("mean_vwc", "mean_orp", "mean_st")
env_labels <- c("VWC", "ORP", "Soil temp")
gene_vars <- c("log_pmoa", "log_mmox", "log_mcra", "pct_pmoa")
gene_labels_v <- c("pmoA", "mmoX", "mcrA", "%pmoA")

for (comp in levels(env_data$compartment)) {
  cd <- env_data %>% filter(compartment == comp)
  for (i in seq_along(env_vars)) {
    for (j in seq_along(gene_vars)) {
      x <- cd[[env_vars[i]]]
      y <- cd[[gene_vars[j]]]
      valid <- complete.cases(x, y)
      if (sum(valid) >= 5) {
        ct <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
        env_stats <- rbind(env_stats, data.frame(
          compartment = comp, env_var = env_labels[i], gene = gene_labels_v[j],
          n = sum(valid), rho = ct$estimate, p = ct$p.value
        ))
      }
    }
  }
}
env_stats$FDR <- p.adjust(env_stats$p, "BH")
write.csv(env_stats, "outputs/figures/new_paper/environmental_driver_stats.csv", row.names = FALSE)

cat("\nSignificant environmental associations (FDR < 0.05):\n")
sig_env <- env_stats %>% filter(FDR < 0.05) %>% arrange(FDR)
if (nrow(sig_env) > 0) {
  print(sig_env, row.names = FALSE)
} else {
  cat("  None significant after FDR correction\n")
}

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ ANALYSIS C: Gene Abundance vs Internal CH4 Concentrations                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n=== ANALYSIS C: Gene vs Internal CH4 ===\n")

# CH4_int = internal stem CH4 concentration
# ch4_50, ch4_125, ch4_200 = CH4 at 50, 125, 200 cm heights
ch4_data <- meta_filt %>%
  filter(!is.na(pmoa_loose), !is.na(mmox_loose))

cat("Samples with CH4 data:\n")
cat("  CH4_int:", sum(!is.na(ch4_data$CH4_int)), "\n")
cat("  ch4_50:", sum(!is.na(ch4_data$ch4_50)), "\n")
cat("  ch4_125:", sum(!is.na(ch4_data$ch4_125)), "\n")
cat("  ch4_200:", sum(!is.na(ch4_data$ch4_200)), "\n")

# Log-transform CH4 (add 1 to handle zeros, use absolute value for negative)
ch4_data$log_ch4_int <- log10(abs(ch4_data$CH4_int) + 1) * sign(ch4_data$CH4_int)
ch4_data$log_ch4_50 <- log10(abs(ch4_data$ch4_50) + 1) * sign(ch4_data$ch4_50)

# Scatter: pmoA vs internal CH4
c1 <- make_env_scatter(ch4_data, "log_ch4_int", "log_pmoa",
                        expression(log[10]*"(|CH"[4]*" internal| + 1)"),
                        expression(log[10]*"(pmoA + 1)"),
                        "Internal CH4 vs pmoA")

c2 <- make_env_scatter(ch4_data, "log_ch4_int", "log_mmox",
                        expression(log[10]*"(|CH"[4]*" internal| + 1)"),
                        expression(log[10]*"(mmoX + 1)"),
                        "Internal CH4 vs mmoX")

c3 <- make_env_scatter(ch4_data, "log_ch4_int", "log_mcra",
                        expression(log[10]*"(|CH"[4]*" internal| + 1)"),
                        expression(log[10]*"(mcrA + 1)"),
                        "Internal CH4 vs mcrA")

# CH4 at 50cm vs genes (for soil-proximal data)
c4 <- make_env_scatter(ch4_data, "log_ch4_50", "log_pmoa",
                        expression(log[10]*"(|CH"[4]*" 50cm| + 1)"),
                        expression(log[10]*"(pmoA + 1)"),
                        "CH4 at 50cm vs pmoA")

c5 <- make_env_scatter(ch4_data, "log_ch4_50", "log_mmox",
                        expression(log[10]*"(|CH"[4]*" 50cm| + 1)"),
                        expression(log[10]*"(mmoX + 1)"),
                        "CH4 at 50cm vs mmoX")

# %pmoA vs CH4_int (does CH4 availability predict which MMO dominates?)
c6 <- make_env_scatter(ch4_data %>% filter(!is.na(pct_pmoa)),
                        "log_ch4_int", "pct_pmoa",
                        expression(log[10]*"(|CH"[4]*" internal| + 1)"),
                        "% pmoA / (pmoA + mmoX)",
                        "Internal CH4 vs % pmoA")

fig_ch4 <- (c1 + c2) / (c3 + c4) / (c5 + c6) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_gene_vs_ch4_internal.png",
       fig_ch4, width = 14, height = 15, dpi = 300)
cat("Saved: fig_gene_vs_ch4_internal.png\n")

# CH4 statistics
ch4_stats <- data.frame()
ch4_vars <- c("log_ch4_int", "log_ch4_50")
ch4_labels <- c("CH4_int", "CH4_50cm")
for (comp in levels(ch4_data$compartment)) {
  cd <- ch4_data %>% filter(compartment == comp)
  for (i in seq_along(ch4_vars)) {
    for (j in seq_along(gene_vars)) {
      x <- cd[[ch4_vars[i]]]
      y <- cd[[gene_vars[j]]]
      valid <- complete.cases(x, y)
      if (sum(valid) >= 5) {
        ct <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
        ch4_stats <- rbind(ch4_stats, data.frame(
          compartment = comp, ch4_var = ch4_labels[i], gene = gene_labels_v[j],
          n = sum(valid), rho = ct$estimate, p = ct$p.value
        ))
      }
    }
  }
}
ch4_stats$FDR <- p.adjust(ch4_stats$p, "BH")
write.csv(ch4_stats, "outputs/figures/new_paper/ch4_gene_stats.csv", row.names = FALSE)

cat("\nSignificant CH4–gene associations (FDR < 0.05):\n")
sig_ch4 <- ch4_stats %>% filter(FDR < 0.05) %>% arrange(FDR)
if (nrow(sig_ch4) > 0) {
  print(sig_ch4, row.names = FALSE)
} else {
  cat("  None significant after FDR correction\n")
}

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ ANALYSIS D: Methanotroph Community Diversity vs Dominance                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

cat("\n=== ANALYSIS D: Methanotroph Diversity ===\n")

# Classify methanotroph ASVs
source("code/00_harmonization/load_methanotroph_definitions.R")
mt_defs <- load_methanotroph_defs()

tax_df <- as.data.frame(tax_table(ps.ra))
tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs, include_conditional = TRUE)

mt_asvs <- rownames(tax_df)[tax_df$mt_status %in% c("Known", "Putative") & !is.na(tax_df$mt_status)]
cat("Methanotroph ASVs:", length(mt_asvs), "\n")

# Subset phyloseq to methanotroph ASVs only (use rarefied counts, not RA)
ps.mt <- prune_taxa(mt_asvs, ps.rare)

# Compute alpha diversity within the methanotroph community
# Use the rarefied count data (ps.rare subset to mt ASVs)
mt_otu <- as.data.frame(t(otu_table(ps.mt)))  # samples x mt ASVs

# Shannon diversity
mt_shannon <- diversity(mt_otu, index = "shannon")
# Richness (observed ASVs)
mt_richness <- specnumber(mt_otu)
# Simpson (1-D, higher = more even)
mt_simpson <- diversity(mt_otu, index = "simpson")
# Dominance: max proportion (how much the top taxon contributes)
mt_ra <- mt_otu / rowSums(mt_otu)
mt_ra[is.nan(as.matrix(mt_ra))] <- 0
mt_dominance <- apply(mt_ra, 1, max)

# Build diversity data frame
div_df <- data.frame(
  sample = names(mt_shannon),
  mt_shannon = mt_shannon,
  mt_richness = mt_richness,
  mt_simpson = mt_simpson,
  mt_dominance = mt_dominance,
  stringsAsFactors = FALSE
)

# Merge with metadata
div_meta <- merge(div_df, meta_filt, by.x = "sample", by.y = "row.names", all.x = TRUE)
div_meta <- div_meta %>% filter(!is.na(compartment), !is.na(pmoa_loose))

div_meta$compartment <- factor(div_meta$compartment,
  levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

cat("Diversity samples with gene data:", nrow(div_meta), "\n")

# Panel 1: Shannon diversity by compartment
d1 <- ggplot(div_meta, aes(x = compartment, y = mt_shannon, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = compartment_colors) +
  labs(x = NULL, y = "Methanotroph Shannon diversity",
       title = "Methanotroph alpha diversity by compartment") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# Panel 2: Richness by compartment
d2 <- ggplot(div_meta, aes(x = compartment, y = mt_richness, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = compartment_colors) +
  labs(x = NULL, y = "Methanotroph ASV richness",
       title = "Methanotroph richness by compartment") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# Panel 3: Shannon vs pmoA (does diversity predict gene abundance?)
d3 <- make_env_scatter(div_meta, "mt_shannon", "log_pmoa",
                        "Methanotroph Shannon diversity",
                        expression(log[10]*"(pmoA + 1)"),
                        "Methanotroph diversity vs pmoA")

# Panel 4: Shannon vs mmoX
d4 <- make_env_scatter(div_meta, "mt_shannon", "log_mmox",
                        "Methanotroph Shannon diversity",
                        expression(log[10]*"(mmoX + 1)"),
                        "Methanotroph diversity vs mmoX")

# Panel 5: Dominance vs pmoA (does a single dominant taxon predict high pmoA?)
d5 <- make_env_scatter(div_meta, "mt_dominance", "log_pmoa",
                        "Methanotroph dominance (max prop.)",
                        expression(log[10]*"(pmoA + 1)"),
                        "Methanotroph dominance vs pmoA")

# Panel 6: Shannon vs %pmoA (does diversity predict gene ratio?)
d6 <- make_env_scatter(div_meta %>% filter(!is.na(pct_pmoa)),
                        "mt_shannon", "pct_pmoa",
                        "Methanotroph Shannon diversity",
                        "% pmoA / (pmoA + mmoX)",
                        "Methanotroph diversity vs % pmoA")

# Wood-only focus: Shannon vs pmoA in wood compartments
wood_div <- div_meta %>% filter(compartment %in% c("Heartwood", "Sapwood"))

d7 <- ggplot(wood_div, aes(x = mt_shannon, y = log_pmoa, color = compartment)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.02, label.y.npc = 0.95) +
  facet_wrap(~ compartment, scales = "free") +
  scale_color_manual(values = compartment_colors) +
  labs(x = "Methanotroph Shannon diversity",
       y = expression(log[10]*"(pmoA + 1)"),
       title = "Wood focus: methanotroph diversity vs pmoA") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 11, face = "bold"))

# Identify dominant methanotroph family per sample
mt_fam_ra <- mt_ra
mt_tax <- tax_df[mt_asvs, ]
mt_fam_labels <- mt_tax$Family
names(mt_fam_labels) <- rownames(mt_tax)

# Aggregate to family level within methanotroph RA
mt_fams <- unique(na.omit(mt_fam_labels))
mt_fam_mat <- matrix(0, nrow = nrow(mt_fam_ra), ncol = length(mt_fams))
rownames(mt_fam_mat) <- rownames(mt_fam_ra)
colnames(mt_fam_mat) <- mt_fams

for (asv in colnames(mt_fam_ra)) {
  fam <- mt_fam_labels[asv]
  if (!is.na(fam) && fam %in% mt_fams) {
    mt_fam_mat[, fam] <- mt_fam_mat[, fam] + mt_fam_ra[, asv]
  }
}

# Find dominant family per sample
dom_family <- apply(mt_fam_mat, 1, function(x) {
  if (all(x == 0)) return(NA)
  colnames(mt_fam_mat)[which.max(x)]
})

div_meta$dominant_mt_family <- dom_family[div_meta$sample]

d8 <- ggplot(div_meta %>% filter(!is.na(dominant_mt_family),
                                  compartment %in% c("Heartwood", "Sapwood")),
             aes(x = dominant_mt_family, y = log_pmoa, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ compartment) +
  scale_fill_manual(values = compartment_colors) +
  labs(x = "Dominant methanotroph family",
       y = expression(log[10]*"(pmoA + 1)"),
       title = "Wood focus: pmoA by dominant methanotroph family") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.title = element_text(size = 10, face = "bold"))

fig_div <- (d1 + d2) / (d3 + d4) / (d5 + d6) / (d7 + d8) +
  plot_annotation(tag_levels = "a")

ggsave("outputs/figures/new_paper/fig_methanotroph_diversity.png",
       fig_div, width = 14, height = 20, dpi = 300)
cat("Saved: fig_methanotroph_diversity.png\n")

# Diversity stats
cat("\nMethanotroph diversity summary by compartment:\n")
div_summary <- div_meta %>%
  group_by(compartment) %>%
  summarize(
    n = n(),
    mean_shannon = mean(mt_shannon, na.rm = TRUE),
    sd_shannon = sd(mt_shannon, na.rm = TRUE),
    mean_richness = mean(mt_richness, na.rm = TRUE),
    mean_dominance = mean(mt_dominance, na.rm = TRUE),
    .groups = "drop"
  )
print(div_summary, n = Inf)

# Kruskal-Wallis test for Shannon by compartment
kw <- kruskal.test(mt_shannon ~ compartment, data = div_meta)
cat("\nKruskal-Wallis test (Shannon ~ compartment):\n")
cat("  chi-squared =", round(kw$statistic, 2), ", p =", format(kw$p.value, digits = 3), "\n")

# Dominant family frequency table (wood only)
cat("\nDominant methanotroph family frequency (wood):\n")
wood_dom <- div_meta %>%
  filter(compartment %in% c("Heartwood", "Sapwood"), !is.na(dominant_mt_family)) %>%
  count(compartment, dominant_mt_family) %>%
  arrange(compartment, desc(n)) %>%
  as.data.frame()
print(wood_dom)

cat("\nAll analyses complete.\n")
