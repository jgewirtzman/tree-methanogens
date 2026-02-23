# ==============================================================================
# Combined Methanogen + Methanotroph Composition (4-panel layout)
# ==============================================================================
# Purpose: Combines the methanogen and methanotroph 16S composition figures
#   into a single 4-panel layout:
#     (a) Methanogen relative abundance
#     (b) Methanogen taxonomic composition
#     (c) Methanotroph relative abundance
#     (d) Methanotroph taxonomic composition
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   Same as 08_methanogen_16s_composition.R and 08b_methanotroph_16s_composition.R
#
# Outputs:
#   - fig5_combined_methane_cycling_composition.png
#
# Required packages: phyloseq, tidyverse, patchwork, viridis
# ==============================================================================

# Load required libraries
library(phyloseq)
library(tidyverse)
library(patchwork)
library(viridis)

# ==============================================================================
# STEP 1: Build phyloseq object (shared for both methanogen and methanotroph)
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

# Filter samples
samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)
samp_df <- samp_df %>% filter(!is.na(compartment), !is.na(species.x), species.x != "")
sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(!is.na(sample_data(ps.ra)$compartment), ps.ra)

# Extract tables
otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))
samp_meta <- data.frame(sample_data(ps.filt))

# ==============================================================================
# STEP 2: Species mapping (shared)
# ==============================================================================

species_mapping <- c(
  "ACRU" = "A. rubrum", "ACSA" = "A. saccharum",
  "BEAL" = "B. alleghaniensis", "BELE" = "B. lenta", "BEPA" = "B. papyrifera",
  "FAGR" = "F. grandifolia", "FRAM" = "F. americana",
  "PIST" = "P. strobus", "QURU" = "Q. rubra", "TSCA" = "T. canadensis",
  "CAOV" = "C. ovata", "KALA" = "K. latifolia", "PRSE" = "P. serotina",
  "QUAL" = "Q. alba", "QUVE" = "Q. velutina", "SAAL" = "S. albidum"
)

compartment_colors <- c("Heartwood" = "#a6611a", "Sapwood" = "#dfc27d",
                          "Mineral Soil" = "#80cdc1", "Organic Soil" = "#018571")

# ==============================================================================
# STEP 3: Methanogen calculations
# ==============================================================================

methanogen_families <- c(
  "Methanobacteriaceae", "Methanomassiliicoccaceae", "Methanoregulaceae",
  "Methanocellaceae", "Methanosaetaceae", "Methanomicrobiaceae",
  "Methanosarcinaceae", "Methanomethyliaceae", "Methanocorpusculaceae"
)

mg_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]
samp_meta$methanogen_pct <- colSums(otu_df[mg_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]

for (fam in methanogen_families) {
  fam_asvs <- rownames(tax_df)[which(tax_df$Family == fam)]
  if (length(fam_asvs) > 0) {
    samp_meta[[paste0("mg_", fam)]] <- colSums(otu_df[fam_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
  } else {
    samp_meta[[paste0("mg_", fam)]] <- 0
  }
}

# ==============================================================================
# STEP 4: Methanotroph calculations
# ==============================================================================

mt_exclusive_families <- c("Methylococcaceae", "Methylacidiphilaceae",
                            "Methylomonadaceae", "Methylomirabilaceae")
beij_mt_genera <- c("Methylocapsa", "Methylocella", "Methylorosula", "Methylovirgula",
                     "Methyloferula", "Methylobacterium-Methylorubrum", "1174-901-12", "Roseiarcus")
mcyst_mt_genera <- c("Methylosinus", "Methylocystis")
other_mt_families <- c("Methylophilaceae", "Methylopilaceae", "Methyloligellaceae")

mt_asvs_excl <- rownames(tax_df)[tax_df$Family %in% mt_exclusive_families]
mt_asvs_beij <- rownames(tax_df)[tax_df$Family == "Beijerinckiaceae" & tax_df$Genus %in% beij_mt_genera]
mt_asvs_mcyst <- rownames(tax_df)[tax_df$Family == "Methylocystaceae" & tax_df$Genus %in% mcyst_mt_genera]
mt_asvs_other <- rownames(tax_df)[tax_df$Family %in% other_mt_families]
mt_asvs_creno <- rownames(tax_df)[tax_df$Genus == "Crenothrix"]
all_mt_asvs <- unique(c(mt_asvs_excl, mt_asvs_beij, mt_asvs_mcyst, mt_asvs_other, mt_asvs_creno))

samp_meta$methanotroph_pct <- colSums(otu_df[all_mt_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]

# Classify ASVs into display families
tax_df$mt_family <- NA_character_
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_excl] <- tax_df$Family[rownames(tax_df) %in% mt_asvs_excl]
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_beij] <- "Beijerinckiaceae"
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_mcyst] <- "Methylocystaceae"
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_other] <- tax_df$Family[rownames(tax_df) %in% mt_asvs_other]
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_creno] <- "Crenotrichaceae"

mt_display_families <- sort(unique(tax_df$mt_family[!is.na(tax_df$mt_family)]))

for (dfam in mt_display_families) {
  dfam_asvs <- rownames(tax_df)[which(tax_df$mt_family == dfam)]
  if (length(dfam_asvs) > 0) {
    samp_meta[[paste0("mt_", dfam)]] <- colSums(otu_df[dfam_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
  } else {
    samp_meta[[paste0("mt_", dfam)]] <- 0
  }
}

# ==============================================================================
# STEP 5: Prepare plot data
# ==============================================================================

plot_data <- samp_meta
plot_data$species_label <- species_mapping[plot_data$species.x]
plot_data <- plot_data %>% filter(!is.na(species_label))
plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
sp_order <- sort(unique(plot_data$species_label))
plot_data$species_label <- factor(plot_data$species_label, levels = sp_order)

# ==============================================================================
# STEP 6: Panel (a) — Methanogen relative abundance
# ==============================================================================

summary_mg_a <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(methanogen_pct, na.rm = TRUE),
            se_pct = sd(methanogen_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# ==============================================================================
# Panel (a) — Methanogen relative abundance
# ==============================================================================

p_mg_a <- ggplot(summary_mg_a, aes(x = species_label, y = mean_pct, fill = compartment)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_pct - se_pct, 0), ymax = mean_pct + se_pct),
                width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = compartment_colors, name = "Compartment") +
  labs(x = NULL, y = "Methanogen rel.\nabundance (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(5, 10, 2, 10))

# ==============================================================================
# Panel (b) — Methanogen family composition
# ==============================================================================

mg_family_cols <- paste0("mg_", methanogen_families)
mg_family_cols <- intersect(mg_family_cols, colnames(plot_data))

mg_comp <- plot_data %>%
  select(species_label, compartment, all_of(mg_family_cols)) %>%
  pivot_longer(cols = all_of(mg_family_cols), names_to = "Family", values_to = "Abundance") %>%
  mutate(Family = sub("^mg_", "", Family))

summary_mg_b <- mg_comp %>%
  group_by(species_label, compartment, Family) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Sort families alphabetically so color ramp matches legend order
mg_families_sorted <- sort(methanogen_families)
mg_colors <- setNames(viridis(length(mg_families_sorted), option = "D"), mg_families_sorted)
summary_mg_b$Family <- factor(summary_mg_b$Family, levels = mg_families_sorted)

p_mg_b <- ggplot(summary_mg_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mg_colors, name = "Methanogen\nFamily") +
  labs(x = NULL, y = "Proportion of\nmethanogens (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "right", legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(2, 10, 2, 10))

# ==============================================================================
# Panel (c) — Methanotroph relative abundance
# ==============================================================================

summary_mt_a <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(methanotroph_pct, na.rm = TRUE),
            se_pct = sd(methanotroph_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

p_mt_a <- ggplot(summary_mt_a, aes(x = species_label, y = mean_pct, fill = compartment)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_pct - se_pct, 0), ymax = mean_pct + se_pct),
                width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = compartment_colors, name = "Compartment") +
  labs(x = NULL, y = "Methanotroph rel.\nabundance (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(2, 10, 2, 10))

# ==============================================================================
# Panel (d) — Methanotroph family composition
# ==============================================================================

mt_family_cols <- paste0("mt_", mt_display_families)
mt_family_cols <- intersect(mt_family_cols, colnames(plot_data))

mt_comp <- plot_data %>%
  select(species_label, compartment, all_of(mt_family_cols)) %>%
  pivot_longer(cols = all_of(mt_family_cols), names_to = "Family", values_to = "Abundance") %>%
  mutate(Family = sub("^mt_", "", Family))

summary_mt_b <- mt_comp %>%
  group_by(species_label, compartment, Family) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Sort families alphabetically so color ramp matches legend order
mt_families_sorted <- sort(mt_display_families)
mt_colors <- setNames(viridis(length(mt_families_sorted), option = "D"), mt_families_sorted)
summary_mt_b$Family <- factor(summary_mt_b$Family, levels = mt_families_sorted)

p_mt_b <- ggplot(summary_mt_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mt_colors, name = "Methanotroph\nFamily") +
  labs(x = NULL, y = "Proportion of\nmethanotrophs (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "right", legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(2, 10, 5, 10))

# ==============================================================================
# STEP 10: Combine all panels vertically with shared x-axis
# ==============================================================================

# Vertical stack: a (methanogen abundance) / b (methanogen composition) /
#                 c (methanotroph abundance) / d (methanotroph composition)
# Compartment facet headers only on panel (a), species labels only on panel (d)
# Add bold tags to each panel individually (ggplot2 4.x compatible)
tag_theme <- theme(plot.tag = element_text(size = 16, face = "bold"))
p_mg_a <- p_mg_a + tag_theme
p_mg_b <- p_mg_b + tag_theme
p_mt_a <- p_mt_a + tag_theme
p_mt_b <- p_mt_b + tag_theme

fig_combined <- p_mg_a / p_mg_b / p_mt_a / p_mt_b +
  plot_layout(heights = c(1, 1.2, 1, 1.2)) +
  plot_annotation(tag_levels = "a")

print(fig_combined)

ggsave("../../outputs/figures/main/fig5_combined_methane_cycling_composition.png",
       fig_combined, width = 12, height = 12, dpi = 300)

cat("Combined figure saved to outputs/figures/fig5_combined_methane_cycling_composition.png\n")
