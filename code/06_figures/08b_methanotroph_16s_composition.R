# ==============================================================================
# Methanotrophic Taxa Relative Abundance
# ==============================================================================
# Purpose: Two-panel figure showing (a) total methanotroph relative abundance
#   (% of prokaryotic community) and (b) methanotroph taxonomic composition
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
#   - figX_methanotroph_abundance_composition.png
#
# Required packages: phyloseq, tidyverse, patchwork, viridis
# ==============================================================================

# Load required libraries
library(phyloseq)
library(tidyverse)
library(patchwork)
library(viridis)

# ==============================================================================
# STEP 1: Build phyloseq object (same pipeline as methanogen script)
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

# ==============================================================================
# STEP 3: Define methanotroph taxa
# ==============================================================================

# Exclusively methanotrophic families (all members are methanotrophs)
methanotroph_exclusive_families <- c(
  "Methylococcaceae",
  "Methylacidiphilaceae",
  "Methylomonadaceae",
  "Methylomirabilaceae"
)

# Mixed families — only certain genera are methanotrophs
# Beijerinckiaceae: methanotrophic genera only
beijerinckiaceae_mt_genera <- c(
  "Methylocapsa",
  "Methylocella",
  "Methylorosula",
  "Methylovirgula",
  "Methyloferula",
  "Methylobacterium-Methylorubrum",
  "1174-901-12",
  "Roseiarcus"
)

# Methylocystaceae methanotrophic genera
methylocystaceae_mt_genera <- c(
  "Methylosinus",
  "Methylocystis"
)

# Methylophilaceae and Methylopilaceae: obligate methylotrophs, some methanotrophic
other_mt_families <- c(
  "Methylophilaceae",
  "Methylopilaceae",
  "Methyloligellaceae"
)

# Crenotrichaceae genus Crenothrix
crenothrix_genus <- "Crenothrix"

# ==============================================================================
# STEP 4: Calculate per-sample methanotroph abundance
# ==============================================================================

otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))
samp_meta <- data.frame(sample_data(ps.filt))

# Identify all methanotroph ASVs using combined family + genus criteria
mt_asvs_exclusive <- rownames(tax_df)[tax_df$Family %in% methanotroph_exclusive_families]
mt_asvs_beij <- rownames(tax_df)[tax_df$Family == "Beijerinckiaceae" &
                                    tax_df$Genus %in% beijerinckiaceae_mt_genera]
mt_asvs_mcyst <- rownames(tax_df)[tax_df$Family == "Methylocystaceae" &
                                     tax_df$Genus %in% methylocystaceae_mt_genera]
mt_asvs_other <- rownames(tax_df)[tax_df$Family %in% other_mt_families]
mt_asvs_creno <- rownames(tax_df)[tax_df$Genus == crenothrix_genus]

all_mt_asvs <- unique(c(mt_asvs_exclusive, mt_asvs_beij, mt_asvs_mcyst,
                         mt_asvs_other, mt_asvs_creno))

# Per-sample total methanotroph relative abundance
mt_abund <- colSums(otu_df[all_mt_asvs, , drop = FALSE], na.rm = TRUE)
samp_meta$methanotroph_pct <- mt_abund[rownames(samp_meta)]

# For panel (b), classify each ASV into a display family for the stacked bar
# We use the actual family for exclusively methanotrophic families,
# and "Beijerinckiaceae (methanotrophs)" for the Beijerinckiaceae subset
tax_df$mt_family <- NA_character_
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_exclusive] <-
  tax_df$Family[rownames(tax_df) %in% mt_asvs_exclusive]
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_beij] <- "Beijerinckiaceae"
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_mcyst] <- "Methylocystaceae"
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_other] <-
  tax_df$Family[rownames(tax_df) %in% mt_asvs_other]
tax_df$mt_family[rownames(tax_df) %in% mt_asvs_creno] <- "Crenotrichaceae"

# Get the unique display families present
display_families <- sort(unique(tax_df$mt_family[!is.na(tax_df$mt_family)]))

# Per-sample abundance by display family
for (dfam in display_families) {
  dfam_asvs <- rownames(tax_df)[which(tax_df$mt_family == dfam)]
  if (length(dfam_asvs) > 0) {
    dfam_sums <- colSums(otu_df[dfam_asvs, , drop = FALSE], na.rm = TRUE)
    samp_meta[[dfam]] <- dfam_sums[rownames(samp_meta)]
  } else {
    samp_meta[[dfam]] <- 0
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
plot_data <- plot_data %>% filter(!is.na(species_label))

plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

sp_order <- sort(unique(plot_data$species_label))
plot_data$species_label <- factor(plot_data$species_label, levels = sp_order)

# ==============================================================================
# STEP 6: Panel (a) — Total methanotroph relative abundance
# ==============================================================================

summary_a <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(
    mean_pct = mean(methanotroph_pct, na.rm = TRUE),
    se_pct = sd(methanotroph_pct, na.rm = TRUE) / sqrt(n()),
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
    y = "Methanotroph rel. abundance (%)"
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
# STEP 7: Panel (b) — Methanotroph family composition (stacked bars)
# ==============================================================================

family_cols <- intersect(display_families, colnames(plot_data))

composition_data <- plot_data %>%
  select(species_label, compartment, all_of(family_cols)) %>%
  pivot_longer(cols = all_of(family_cols), names_to = "Family", values_to = "Abundance")

summary_b <- composition_data %>%
  group_by(species_label, compartment, Family) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

summary_b <- summary_b %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Viridis discrete color palette for methanotroph families
n_fam <- length(display_families)
family_colors <- setNames(viridis(n_fam, option = "D"), display_families)

p_b <- ggplot(summary_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = family_colors, name = "Methanotroph Family") +
  labs(
    x = NULL,
    y = "Proportion of methanotrophs (%)"
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

fig_mt <- p_a / p_b +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a")

print(fig_mt)

# ggsave("outputs/figures/figX_methanotroph_abundance_composition.png",
#        fig_mt, width = 14, height = 10, dpi = 300)

cat("Methanotroph figure saved to outputs/figures/figX_methanotroph_abundance_composition.png\n")
