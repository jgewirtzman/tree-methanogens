# ==============================================================================
# Methanogenic Taxa Relative Abundance (Figure 5)
# ==============================================================================
# Purpose: Two-panel figure showing (a) total methanogen relative abundance
#   (% of prokaryotic community) and (b) methanogen taxonomic composition
#   (family-level) across tree species and forest compartments (heartwood,
#   sapwood, mineral soil, organic soil) from 16S rRNA gene sequencing.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - OTU_table.txt (from data/raw/16s/)
#   - taxonomy_table.txt (from data/raw/16s/)
#   - unrooted_tree.nwk (from data/raw/16s/)
#   - tree_16s_mapping_dada2_corrected.txt (from data/raw/16s/)
#   - ddPCR_meta_all_data.csv (from data/raw/ddpcr/)
#   - Tree_Core_Sectioning_Data.csv (from data/raw/tree_cores/)
#   - 16s_w_metadata.csv (from data/raw/16s/)
#
# Outputs:
#   - fig5_methanogen_abundance_composition.png
#
# Required packages: phyloseq, tidyverse, patchwork, RColorBrewer
# ==============================================================================

# Load required libraries
library(phyloseq)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

# ==============================================================================
# STEP 1: Build phyloseq object (same pipeline as other 16S scripts)
# ==============================================================================

#### Import Metadata ####
ddpcr <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
water <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv")
qpcr_16s <- subset(qpcr_16s, qpcr_16s$Sample.ID != "None")
qpcr_16s <- qpcr_16s[, c(3, 4, 6)]

# Process metadata
water$seq_id <- toupper(water$seq_id)
ddpcr <- merge(ddpcr, water, by = c("core_type", "seq_id"), all.x = TRUE)
ddpcr <- merge(ddpcr, qpcr_16s, by = c("core_type", "seq_id"), all.x = TRUE)

#### Process 16S Data ####
otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)

# Process taxonomy table (embedded in last 7 columns)
bastard_tax <- otu_tab[, 590:596]
bastard_tax[bastard_tax == ""] <- NA
tax_tab_pre <- tax_table(bastard_tax)
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))

# Process OTU table
otu_tab_corr <- otu_tab[, 1:589]
otu_table_pre <- otu_table(otu_tab_corr, taxa_are_rows = TRUE)

# Import sample data and phylogenetic tree
phylo_tree <- read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
samp_data$RowName <- row.names(samp_data)

# Fix sample names
samp_data$seq_id <- sub("prime", "'", samp_data$seq_id)
samp_data$seq_id <- sub("star", "*", samp_data$seq_id)
samp_data$seq_id <- sub("HM", "H", samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

# Merge metadata
samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)

# Remove duplicates
dups <- which(duplicated(samp_data_merged$RowName) == TRUE)
samp_data_merged <- samp_data_merged[-c(dups), ]
row.names(samp_data_merged) <- samp_data_merged$RowName

# Create phyloseq object
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

# Rename ASVs and rarefy
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)

# Assign clean taxonomy column names
colnames(tax_table(ps.rare)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Transform to relative abundance (proportion of total community)
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# ==============================================================================
# STEP 2: Define compartments and filter samples
# ==============================================================================

# Create compartment labels from core_type and material
samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)

# Keep only samples with defined compartments and known species
# Note: after merge, species column becomes species.x
samp_df <- samp_df %>% filter(!is.na(compartment), !is.na(species.x), species.x != "")
sample_data(ps.ra) <- sample_data(samp_df)

# Subset phyloseq to matching samples
ps.filt <- prune_samples(!is.na(sample_data(ps.ra)$compartment), ps.ra)

# ==============================================================================
# STEP 3: Define methanogen families
# ==============================================================================

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

# ==============================================================================
# STEP 4: Calculate per-sample methanogen abundance
# ==============================================================================

# Extract OTU and taxonomy tables
otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))
samp_meta <- data.frame(sample_data(ps.filt))
sample_names <- colnames(otu_df)  # These are the phyloseq sample names

# Identify methanogen ASVs
methanogen_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]

# Per-sample total methanogen relative abundance
methanogen_abund <- colSums(otu_df[methanogen_asvs, , drop = FALSE], na.rm = TRUE)
samp_meta$methanogen_pct <- methanogen_abund[rownames(samp_meta)]

# Per-sample methanogen abundance by family — build directly into samp_meta
for (fam in methanogen_families) {
  fam_asvs <- rownames(tax_df)[which(tax_df$Family == fam)]
  if (length(fam_asvs) > 0) {
    fam_sums <- colSums(otu_df[fam_asvs, , drop = FALSE], na.rm = TRUE)
    samp_meta[[fam]] <- fam_sums[rownames(samp_meta)]
  } else {
    samp_meta[[fam]] <- 0
  }
}

plot_data <- samp_meta

# ==============================================================================
# STEP 5: Species name mapping and ordering
# ==============================================================================

species_mapping <- c(
  "ACRU" = "A. rubrum",
  "ACSA" = "A. saccharum",
  "BEAL" = "B. alleghaniensis",
  "BELE" = "B. lenta",
  "BEPA" = "B. papyrifera",
  "FAGR" = "F. grandifolia",
  "FRAM" = "F. americana",
  "PIST" = "P. strobus",
  "QURU" = "Q. rubra",
  "TSCA" = "T. canadensis",
  "CAOV" = "C. ovata",
  "KALA" = "K. latifolia",
  "PRSE" = "P. serotina",
  "QUAL" = "Q. alba",
  "QUVE" = "Q. velutina",
  "SAAL" = "S. albidum"
)

plot_data$species_label <- species_mapping[plot_data$species.x]
# Drop samples with unmapped species
plot_data <- plot_data %>% filter(!is.na(species_label))

# Order compartments
plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

# Order species alphabetically by full name
sp_order <- sort(unique(plot_data$species_label))
plot_data$species_label <- factor(plot_data$species_label, levels = sp_order)

# ==============================================================================
# STEP 6: Panel (a) — Total methanogen relative abundance
# ==============================================================================

# Summarize by species × compartment
summary_a <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(
    mean_pct = mean(methanogen_pct, na.rm = TRUE),
    se_pct = sd(methanogen_pct, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_a <- ggplot(summary_a, aes(x = species_label, y = mean_pct, fill = compartment)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_pct - se_pct, 0), ymax = mean_pct + se_pct),
                width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(
    values = c("Heartwood" = "#a6611a", "Sapwood" = "#dfc27d",
               "Mineral Soil" = "#80cdc1", "Organic Soil" = "#018571"),
    name = "Compartment"
  ) +
  labs(
    x = NULL,
    y = "Methanogen rel. abundance (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.margin = margin(5, 10, 5, 10)
  )

# ==============================================================================
# STEP 7: Panel (b) — Methanogen family composition (stacked bars)
# ==============================================================================

# Reshape family abundance data to long format
family_cols <- intersect(methanogen_families, colnames(plot_data))

composition_data <- plot_data %>%
  select(species_label, compartment, all_of(family_cols)) %>%
  pivot_longer(cols = all_of(family_cols), names_to = "Family", values_to = "Abundance")

# Summarize mean abundance by species × compartment × family
summary_b <- composition_data %>%
  group_by(species_label, compartment, Family) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Calculate total methanogen abundance per group for normalization
summary_b <- summary_b %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# For species×compartment combos with zero methanogens, keep them so bars appear
# (they'll be empty/grey indicating no detection)

# Viridis discrete color palette for methanogen families
library(viridis)
n_families <- length(methanogen_families)
family_colors <- setNames(viridis(n_families, option = "D"), methanogen_families)

p_b <- ggplot(summary_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = family_colors, name = "Methanogen Family") +
  labs(
    x = NULL,
    y = "Proportion of methanogens (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 8),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.margin = margin(5, 10, 5, 10)
  )

# ==============================================================================
# STEP 8: Combine and save
# ==============================================================================

fig5 <- p_a / p_b +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a")

print(fig5)

# ggsave("outputs/figures/fig5_methanogen_abundance_composition.png",
#        fig5, width = 14, height = 10, dpi = 300)

cat("Figure 5 saved to outputs/figures/fig5_methanogen_abundance_composition.png\n")
