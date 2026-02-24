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
# STEP 4: Methanotroph calculations (using canonical definitions)
# ==============================================================================

source("code/00_harmonization/load_methanotroph_definitions.R")
mt_defs <- load_methanotroph_defs()

# Classify each ASV as Known, Putative, or NA
tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs, include_conditional = TRUE)

known_asvs   <- rownames(tax_df)[tax_df$mt_status == "Known" & !is.na(tax_df$mt_status)]
putative_asvs <- rownames(tax_df)[tax_df$mt_status == "Putative" & !is.na(tax_df$mt_status)]
all_mt_asvs  <- c(known_asvs, putative_asvs)

cat("Known methanotroph ASVs:", length(known_asvs), "\n")
cat("Putative methanotroph ASVs:", length(putative_asvs), "\n")

# Per-sample abundance: total, known, putative
samp_meta$methanotroph_pct <- colSums(otu_df[all_mt_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
samp_meta$mt_known_pct     <- colSums(otu_df[known_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
samp_meta$mt_putative_pct  <- colSums(otu_df[putative_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]

# Assign display family (using the ASV's Family) + status for grouping
tax_df$mt_family <- NA_character_
tax_df$mt_family[!is.na(tax_df$mt_status)] <- tax_df$Family[!is.na(tax_df$mt_status)]

# Build per-family per-status abundance columns
mt_asvs_all <- rownames(tax_df)[!is.na(tax_df$mt_status)]
mt_fam_status <- unique(tax_df[mt_asvs_all, c("mt_family", "mt_status")])
mt_fam_status <- mt_fam_status[order(mt_fam_status$mt_status, mt_fam_status$mt_family), ]

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
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(5, 10, 0, 10))

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
mg_colors <- setNames(viridis(length(mg_families_sorted), option = "D", direction = -1, alpha = 0.7), mg_families_sorted)
summary_mg_b$Family <- factor(summary_mg_b$Family, levels = mg_families_sorted)

p_mg_b <- ggplot(summary_mg_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mg_colors, name = "Methanogen\nFamily") +
  labs(x = NULL, y = "Proportion of\nmethanogens (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 0, 10))

# ==============================================================================
# Panel (c) — Methanotroph relative abundance (Known + Putative stacked)
# ==============================================================================

summary_mt_known <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(mt_known_pct, na.rm = TRUE),
            se_pct = sd(mt_known_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(Status = "Known")

summary_mt_putative <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(mt_putative_pct, na.rm = TRUE),
            se_pct = sd(mt_putative_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(Status = "Putative")

summary_mt_a <- bind_rows(summary_mt_known, summary_mt_putative) %>%
  mutate(Status = factor(Status, levels = c("Known", "Putative")))

# Error bars on total (Known + Putative combined)
summary_mt_total <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_total = mean(methanotroph_pct, na.rm = TRUE),
            se_total = sd(methanotroph_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Known bars match panel (a) compartment colors; Putative in grey
# Create a composite fill variable: Known bars get compartment name, Putative bars get "Putative"
summary_mt_a$fill_group <- ifelse(summary_mt_a$Status == "Known",
                                   as.character(summary_mt_a$compartment),
                                   "Putative")

# Build fill palette: compartment colors for Known + grey for Putative
panel_c_fills <- c(compartment_colors, "Putative" = "grey70")

# Factor ordering: compartment levels first (so Known stacks on bottom), then Putative on top
summary_mt_a$fill_group <- factor(summary_mt_a$fill_group,
                                   levels = c(names(compartment_colors), "Putative"))

p_mt_a <- ggplot(summary_mt_a, aes(x = species_label, y = mean_pct, fill = fill_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(data = summary_mt_total,
                aes(x = species_label, y = mean_total,
                    ymin = pmax(mean_total - se_total, 0),
                    ymax = mean_total + se_total),
                inherit.aes = FALSE, width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = panel_c_fills,
                    breaks = "Putative",
                    labels = "Putative",
                    name = "Classification") +
  labs(x = NULL, y = "Methanotroph rel.\nabundance (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.justification = c(0, 0.5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 0, 10))

# ==============================================================================
# Panel (d) — Methanotroph family composition (Known families first, then Putative)
# ==============================================================================

# Build per-sample abundance by family + status
mt_fam_status_unique <- unique(tax_df[!is.na(tax_df$mt_status),
                                       c("mt_family", "mt_status")])

for (i in seq_len(nrow(mt_fam_status_unique))) {
  fam <- mt_fam_status_unique$mt_family[i]
  stat <- mt_fam_status_unique$mt_status[i]
  col_name <- paste0("mt_", stat, "_", fam)
  fam_stat_asvs <- rownames(tax_df)[which(tax_df$mt_family == fam &
                                            tax_df$mt_status == stat)]
  if (length(fam_stat_asvs) > 0) {
    samp_meta[[col_name]] <- colSums(otu_df[fam_stat_asvs, , drop = FALSE],
                                     na.rm = TRUE)[rownames(samp_meta)]
  } else {
    samp_meta[[col_name]] <- 0
  }
}

# Rebuild plot_data with new columns
plot_data <- samp_meta
plot_data$species_label <- species_mapping[plot_data$species.x]
plot_data <- plot_data %>% filter(!is.na(species_label))
plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
plot_data$species_label <- factor(plot_data$species_label, levels = sp_order)

mt_fam_stat_cols <- grep("^mt_(Known|Putative)_", colnames(plot_data), value = TRUE)

mt_comp <- plot_data %>%
  select(species_label, compartment, all_of(mt_fam_stat_cols)) %>%
  pivot_longer(cols = all_of(mt_fam_stat_cols),
               names_to = "Family_Status", values_to = "Abundance") %>%
  mutate(
    Status = sub("^mt_(Known|Putative)_.*", "\\1", Family_Status),
    Family = sub("^mt_(Known|Putative)_", "", Family_Status)
  )

summary_mt_b <- mt_comp %>%
  group_by(species_label, compartment, Family, Status) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Drop families with zero total abundance (not present in plotted samples)
fam_totals <- summary_mt_b %>% group_by(Family, Status) %>%
  summarize(total_abund = sum(mean_abund), .groups = "drop") %>%
  filter(total_abund > 0)
summary_mt_b <- summary_mt_b %>%
  semi_join(fam_totals, by = c("Family", "Status"))

# Distinct color for each Family x Status combination
all_mt_fams <- sort(unique(summary_mt_b$Family))
summary_mt_b$Family_status <- paste0(summary_mt_b$Family, " (", summary_mt_b$Status, ")")

# Order levels: for each family, Known first then Putative
status_levels <- c()
for (fam in all_mt_fams) {
  if (any(summary_mt_b$Family_status == paste0(fam, " (Known)")))
    status_levels <- c(status_levels, paste0(fam, " (Known)"))
  if (any(summary_mt_b$Family_status == paste0(fam, " (Putative)")))
    status_levels <- c(status_levels, paste0(fam, " (Putative)"))
}
summary_mt_b$Family_status <- factor(summary_mt_b$Family_status, levels = status_levels)

mt_colors <- setNames(viridis(length(status_levels), option = "D", direction = -1, alpha = 0.7), status_levels)

p_mt_b <- ggplot(summary_mt_b, aes(x = species_label, y = proportion, fill = Family_status)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mt_colors, name = "Methanotroph\nFamily") +
  labs(x = NULL, y = "Proportion of\nmethanotrophs (%)") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 9),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 5, 10))

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

ggsave("outputs/figures/main/fig5_combined_methane_cycling_composition.png",
       fig_combined, width = 14, height = 12, dpi = 300)

cat("Combined figure saved to outputs/figures/fig5_combined_methane_cycling_composition.png\n")
