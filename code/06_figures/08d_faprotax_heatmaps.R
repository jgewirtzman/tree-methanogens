# ==============================================================================
# FAPROTAX Functional Prediction Heatmaps (Figure S3)
# ==============================================================================
# Purpose: Creates a two-panel figure showing FAPROTAX functional predictions
#   from 16S rRNA amplicon data:
#     Panel (a): Average relative abundance of metabolisms per core type
#                (Inner vs Outer), grouped by functional category
#     Panel (b): Average relative abundance of selected metabolisms per species
#                (heartwood only)
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - OTU_table.txt (from data/raw/16s/)
#   - unrooted_tree.nwk (from data/raw/16s/)
#   - tree_16s_mapping_dada2_corrected.txt (from data/raw/16s/)
#   - ddPCR_meta_all_data.csv (from data/raw/ddpcr/)
#   - Tree_Core_Sectioning_Data.csv (from data/raw/tree_cores/)
#   - 16s_w_metadata.csv (from data/raw/16s/)
#
# Outputs:
#   - figS3_faprotax_heatmaps.png (in outputs/figures/supplementary/)
#
# Required packages: phyloseq, microeco, tidyverse, patchwork
# ==============================================================================

library(phyloseq)
library(microeco)
library(tidyverse)
library(patchwork)

# ==============================================================================
# STEP 1: Build phyloseq object (same pipeline as 08c_combined_...)
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

# Remove mitochondria and chloroplasts
pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
no_mito <- pop_taxa(raw_ps, c(mitochondria, chloroplast))

# Rename ASVs and rarefy
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)

# Transform to relative abundance (0-100% scale)
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# Set proper taxonomy column names
colnames(tax_table(ps.ra)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ==============================================================================
# STEP 2: Filter for wood samples (keep ALL species)
# ==============================================================================

ps.wood <- subset_samples(ps.ra, material == "Wood")

# ==============================================================================
# STEP 3: Prepare taxonomy for microeco (add QIIME-style prefixes)
# ==============================================================================

tax_table(ps.wood)[, 1] <- paste0("k__", tax_table(ps.wood)[, 1])
tax_table(ps.wood)[, 2] <- paste0("p__", tax_table(ps.wood)[, 2])
tax_table(ps.wood)[, 3] <- paste0("c__", tax_table(ps.wood)[, 3])
tax_table(ps.wood)[, 4] <- paste0("o__", tax_table(ps.wood)[, 4])
tax_table(ps.wood)[, 5] <- paste0("f__", tax_table(ps.wood)[, 5])
tax_table(ps.wood)[, 6] <- paste0("g__", tax_table(ps.wood)[, 6])
tax_table(ps.wood)[, 7] <- paste0("s__", tax_table(ps.wood)[, 7])

# ==============================================================================
# STEP 4: Convert to microeco and run FAPROTAX
# ==============================================================================

# --- Panel (a): Merge by core_type (Inner / Outer) ---
meco_all <- microtable$new(
  otu_table = as.data.frame(otu_table(ps.wood)),
  tax_table = noquote(as.data.frame(tax_table(ps.wood))),
  sample_table = as.data.frame(as.matrix(sample_data(ps.wood))),
  phylo_tree = phy_tree(ps.wood)
)

meco_core <- clone(meco_all)$merge_samples("core_type")

cat("Running FAPROTAX for core type grouping...\n")
t_func_core <- trans_func$new(meco_core)
t_func_core$cal_spe_func(prok_database = "FAPROTAX")
t_func_core$cal_spe_func_perc(abundance_weighted = TRUE, dec = 2)
perc_core <- t_func_core$res_spe_func_perc
# perc_core: rows = "Inner","Outer"; cols = function names

# --- Panel (b): Heartwood only, merge by species ---
ps.hw <- subset_samples(ps.wood, core_type == "Inner")

meco_hw <- microtable$new(
  otu_table = as.data.frame(otu_table(ps.hw)),
  tax_table = noquote(as.data.frame(tax_table(ps.hw))),
  sample_table = as.data.frame(as.matrix(sample_data(ps.hw))),
  phylo_tree = phy_tree(ps.hw)
)

meco_species <- clone(meco_hw)$merge_samples("species.x")

cat("Running FAPROTAX for species grouping (heartwood)...\n")
t_func_sp <- trans_func$new(meco_species)
t_func_sp$cal_spe_func(prok_database = "FAPROTAX")
t_func_sp$cal_spe_func_perc(abundance_weighted = TRUE, dec = 2)
perc_species <- t_func_sp$res_spe_func_perc
# perc_species: rows = species codes; cols = function names

# ==============================================================================
# STEP 5: Curate metabolisms and categories
# ==============================================================================

# Species code → full species name mapping
species_names <- c(
  "ACRU" = "A. rubrum",       "ACSA" = "A. saccharum",
  "BEAL" = "B. alleghaniensis", "BELE" = "B. lenta",
  "BEPA" = "B. papyrifera",  "CAOV" = "C. ovata",
  "FAGR" = "F. grandifolia",  "FRAM" = "F. americana",
  "KALA" = "K. latifolia",    "PIST" = "P. strobus",
  "PRSE" = "P. serotina",     "QUAL" = "Q. alba",
  "QURU" = "Q. rubra",        "QUVE" = "Q. velutina",
  "SAAL" = "S. albidum",      "TSCA" = "T. canadensis"
)

# Define the curated set of metabolisms and their categories for panel (a)
# Includes all methane-cycling pathways
metabolism_info <- tribble(
  ~metabolism,                                           ~category,
  "photoheterotrophy",                                   "Energy source",
  "aerobic_chemoheterotrophy",                           "Energy source",
  "methanotrophy",                                       "CH4 cycling",
  "methanogenesis",                                      "CH4 cycling",
  "methanogenesis_by_CO2_reduction_with_H2",             "CH4 cycling",
  "methanogenesis_by_reduction_of_methyl_compounds_with_H2", "CH4 cycling",
  "methanol_oxidation",                                  "CH4 cycling",
  "methylotrophy",                                       "CH4 cycling",
  "fermentation",                                        "C-cycle",
  "cellulolysis",                                        "C-cycle",
  "nitrate_respiration",                                 "N-cycle",
  "nitrate_reduction",                                   "N-cycle",
  "aerobic_ammonia_oxidation",                           "N-cycle",
  "nitrogen_fixation",                                   "N-cycle",
  "sulfur_respiration",                                  "S-cycle",
  "sulfate_respiration",                                 "S-cycle",
  "manganese_oxidation",                                 "Others",
  "dark_hydrogen_oxidation",                             "Others"
)

# Clean display names (remove underscores)
metabolism_info <- metabolism_info %>%
  mutate(display_name = str_replace_all(metabolism, "_", " "))

# Category and metabolism order (top to bottom)
cat_order <- c("Energy source", "CH4 cycling", "C-cycle", "N-cycle", "S-cycle", "Others")
metabolism_info$category <- factor(metabolism_info$category, levels = cat_order)

# Create ordered factor for metabolism display names
metabolism_info$display_name <- factor(
  metabolism_info$display_name,
  levels = rev(metabolism_info$display_name)
)

# ==============================================================================
# STEP 6: Build panel (a) data — Heartwood vs Sapwood
# ==============================================================================

# Extract values: perc_core rows = groups, cols = functions
# Display as Heartwood / Sapwood
panel_a_data <- metabolism_info %>%
  rowwise() %>%
  mutate(
    Heartwood = if (metabolism %in% colnames(perc_core)) perc_core["Inner", metabolism] else 0,
    Sapwood   = if (metabolism %in% colnames(perc_core)) perc_core["Outer", metabolism] else 0
  ) %>%
  ungroup()

# Pivot to long format for ggplot
panel_a_long <- panel_a_data %>%
  pivot_longer(cols = c(Heartwood, Sapwood), names_to = "Core_Type", values_to = "Percentage") %>%
  mutate(Core_Type = factor(Core_Type, levels = c("Heartwood", "Sapwood")))

# ==============================================================================
# STEP 7: Build panel (b) data — species x methane metabolisms (heartwood)
# ==============================================================================

# Panel (b): all methane-related pathways + fermentation + dark H2 oxidation
panel_b_metabolisms <- c(
  "methanotrophy", "methanogenesis",
  "methanogenesis_by_CO2_reduction_with_H2",
  "methanogenesis_by_reduction_of_methyl_compounds_with_H2",
  "methanol_oxidation", "methylotrophy",
  "fermentation", "dark_hydrogen_oxidation"
)

# Build info for panel (b)
panel_b_info <- tibble(
  metabolism = panel_b_metabolisms,
  category = c(rep("CH4 cycling", 6), "C-cycle", "Others"),
  display_name = str_replace_all(panel_b_metabolisms, "_", " ")
)
panel_b_info$category <- factor(panel_b_info$category, levels = cat_order)
panel_b_info$display_name <- factor(
  panel_b_info$display_name,
  levels = rev(panel_b_info$display_name)
)

# Get all species present, sorted alphabetically by code
all_species <- sort(rownames(perc_species))

# Map species codes to full names for display
species_display <- species_names[all_species]

# Build panel (b) data
panel_b_data <- panel_b_info %>%
  crossing(species = all_species) %>%
  rowwise() %>%
  mutate(
    Percentage = if (metabolism %in% colnames(perc_species)) perc_species[species, metabolism] else 0
  ) %>%
  ungroup() %>%
  mutate(
    species_label = species_names[species],
    species_label = factor(species_label, levels = species_display)
  )

# ==============================================================================
# STEP 8: Create panel (a) heatmap
# ==============================================================================

# Log-transform percentages for fill so low values are visible
# Raw percentages kept as text labels
panel_a_long <- panel_a_long %>%
  mutate(log_pct = log10(Percentage + 1))
panel_b_data <- panel_b_data %>%
  mutate(log_pct = log10(Percentage + 1))

log_max <- max(c(panel_a_long$log_pct, panel_b_data$log_pct), na.rm = TRUE)

p_a <- ggplot(panel_a_long, aes(x = Core_Type, y = display_name, fill = log_pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Percentage, 2)), size = 3.5, color = "white") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colours = c("#1a1a2e", "#16213e", "#533483", "#e94560", "#f5a623"),
    limits = c(0, log_max),
    name = expression(log[10]*"(% + 1)")
  ) +
  scale_x_discrete(position = "bottom") +
  labs(x = "Core Type", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# ==============================================================================
# STEP 9: Create panel (b) heatmap
# ==============================================================================

p_b <- ggplot(panel_b_data, aes(x = species_label, y = display_name, fill = log_pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Percentage, 2)), size = 2.2, color = "white") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colours = c("#1a1a2e", "#16213e", "#533483", "#e94560", "#f5a623"),
    limits = c(0, log_max),
    name = expression(log[10]*"(% + 1)")
  ) +
  labs(x = "Species (heartwood)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# ==============================================================================
# STEP 10: Combine panels and save
# ==============================================================================

tag_theme <- theme(plot.tag = element_text(size = 12, face = "bold"))
p_a <- p_a + tag_theme
p_b <- p_b + tag_theme

combined <- p_a / p_b +
  plot_layout(heights = c(2, 1.2), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

print(combined)

ggsave("outputs/figures/supplementary/figS3_faprotax_heatmaps.png",
       combined, width = 14, height = 12, dpi = 300)

cat("\nSaved figS3_faprotax_heatmaps.png to outputs/figures/supplementary/\n")
