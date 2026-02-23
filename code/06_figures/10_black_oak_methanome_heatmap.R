# ==============================================================================
# Black Oak Methanogen/Methanotroph Heatmap (Figure S12)
# ==============================================================================
# Purpose: Creates a heatmap of methanogenic and methanotrophic taxa across
#   tissue types in a felled black oak (QUVE), based on 16S rRNA amplicon
#   data. Also computes alpha diversity metrics for methane-cycling taxa.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - OTU_table.txt (from data/raw/16s/black_oak/)
#
# Outputs:
#   - figS12_black_oak_methanome.png (absolute abundance heatmap)
#   - figS12_black_oak_methanome_relative.png (relative abundance heatmap)
#
# Required packages: tidyverse, viridis, cluster
# ==============================================================================

# Load libraries
library(tidyverse)
library(viridis)
library(cluster)

# ==============================================================================
# STEP 1: Load and Prepare Data
# ==============================================================================
otu_table <- read.table("data/raw/16s/black_oak/OTU_table.txt",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define methanogen families
methanogen_families <- c(
  "Bathyarchaeia Class",
  "Methanobacteriaceae",
  "Methanocellaceae",
  "Methanocorpusculaceae",
  "Methanomassiliicoccaceae",
  "Methanomethylophilaceae",
  "Methanomicrobiaceae",
  "Methanoregulaceae",
  "Methanosaetaceae",
  "Methanosarcinaceae",
  "Rice Cluster II",
  "Methanopyraceae",
  "Methanothermaceae",
  "Methanocaldococcaceae",
  "Methanoplanaceae",
  "Methanospirillaceae"
)

# Define methanotroph families
methanotroph_families <- c(
  "Methylococcaceae",
  "Methylacidiphilaceae",
  "Beijerinckiaceae",
  "Methylocystaceae",
  "Methylomonadaceae",
  "Methylophilaceae",
  "Methylopilaceae",
  "Methyloligellaceae",
  "Methylomirabilaceae"
)

# Define methanotroph genera
methanotroph_genera <- c(
  "Methylobacterium-Methylorubrum",
  "Methylocella",
  "Methylorosula",
  "Methylocapsa",
  "1174-901-12",
  "Roseiarcus"
)

# ==============================================================================
# STEP 2: Filter and classify taxa
# ==============================================================================
otu_table <- otu_table %>%
  filter(!is.na(Family)) %>%
  mutate(Group = case_when(
    Family %in% methanogen_families ~ "Methanogen",
    (Family %in% methanotroph_families &
       !(Family %in% c("Methylopilaceae", "Beijerinckiaceae") &
           Genus %in% c("Bosea", "Microvirga", "Psychroglaciecola"))) |
      Genus %in% methanotroph_genera ~ "Methanotroph",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# Pivot to long format and assign tissue types
data_cols <- otu_table %>%
  select(where(is.numeric)) %>%
  colnames()

otu_table_long <- otu_table %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(Media = case_when(
    str_detect(Sample, "HEART") ~ "Heartwood",
    str_detect(Sample, "SAP") ~ "Sapwood",
    str_detect(Sample, "LITTER") ~ "Leaf Litter",
    str_detect(Sample, "BARK") ~ "Bark",
    str_detect(Sample, "FOLIAGE") ~ "Foliage",
    str_detect(Sample, "MINERAL") ~ "Mineral Soil",
    str_detect(Sample, "ORGANIC") ~ "Organic Soil",
    str_detect(Sample, "BRANCH") ~ "Branch",
    str_detect(Sample, "COARSE") ~ "Coarse Root",
    str_detect(Sample, "FINE") ~ "Fine Root",
    str_detect(Sample, "ROT") ~ "Rot",
    TRUE ~ "Other"
  )) %>%
  filter(Media != "Other")

# ==============================================================================
# STEP 3: Absolute abundance heatmap
# ==============================================================================
aggregated_data <- otu_table_long %>%
  group_by(Group, Family, Genus, Media) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Media, values_from = Total_Abundance, values_fill = 0)

# Log transform
log_transformed_data <- aggregated_data %>%
  mutate(across(where(is.numeric), ~ log10(. + 1)))

# Prepare for ggplot
heatmap_data <- log_transformed_data %>%
  pivot_longer(cols = -c(Group, Family, Genus), names_to = "Media", values_to = "Abundance") %>%
  mutate(Row_Label = paste(Family, Genus, sep = " | "))

# Ecological order: tree interior -> exterior -> ground
media_order <- c("Heartwood", "Sapwood", "Rot", "Bark", "Branch",
                 "Foliage", "Leaf Litter", "Coarse Root", "Fine Root",
                 "Organic Soil", "Mineral Soil")

heatmap_data <- heatmap_data %>%
  mutate(Media = factor(Media, levels = media_order))

# Split into methanogen / methanotroph for separate panels with square cells
library(patchwork)

n_cols <- length(media_order)  # 11
heatmap_mg <- heatmap_data %>% filter(Group == "Methanogen")
heatmap_mt <- heatmap_data %>% filter(Group == "Methanotroph")
n_rows_mg <- n_distinct(heatmap_mg$Row_Label)  # 1
n_rows_mt <- n_distinct(heatmap_mt$Row_Label)  # 9

shared_heatmap_theme <- theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )

# Methanogen panel (top)
p_mg <- ggplot(heatmap_mg, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile(color = "white") +
  coord_fixed(ratio = 1) +
  shared_heatmap_theme +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5, 5, 0, 5)) +
  labs(subtitle = "Methanogen") +
  scale_fill_viridis_c(option = "C", name = expression(log[10]~"(Abundance + 1)"))

# Methanotroph panel (bottom)
p_mt <- ggplot(heatmap_mt, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile(color = "white") +
  coord_fixed(ratio = 1) +
  shared_heatmap_theme +
  theme(plot.margin = margin(0, 5, 5, 5)) +
  labs(y = "Taxa (Family | Genus)", subtitle = "Methanotroph") +
  scale_fill_viridis_c(option = "C", name = expression(log[10]~"(Abundance + 1)"))

# Stack with heights proportional to row counts
p_mg <- p_mg + theme(legend.position = "none")
p_abs <- p_mg / p_mt +
  plot_layout(heights = c(n_rows_mg, n_rows_mt))

print(p_abs)

ggsave("outputs/figures/figS12_black_oak_methanome.png",
       p_abs, width = 12, height = 8, dpi = 300)

# ==============================================================================
# STEP 4: Relative abundance heatmap
# ==============================================================================
otu_table_long_relative <- otu_table_long %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

aggregated_relative_data <- otu_table_long_relative %>%
  group_by(Group, Family, Genus, Media) %>%
  summarize(Total_Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Media, values_from = Total_Relative_Abundance, values_fill = 0)

log_transformed_relative_data <- aggregated_relative_data %>%
  mutate(across(where(is.numeric), ~ log10(. + 1)))

heatmap_relative_data <- log_transformed_relative_data %>%
  pivot_longer(cols = -c(Group, Family, Genus), names_to = "Media", values_to = "Abundance") %>%
  mutate(Row_Label = paste(Family, Genus, sep = " | "))

# Same ecological order for relative heatmap
heatmap_relative_data <- heatmap_relative_data %>%
  mutate(Media = factor(Media, levels = media_order))

# Split relative heatmap similarly
heatmap_rel_mg <- heatmap_relative_data %>% filter(Group == "Methanogen")
heatmap_rel_mt <- heatmap_relative_data %>% filter(Group == "Methanotroph")

p_rel_mg <- ggplot(heatmap_rel_mg, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile(color = "white") +
  coord_fixed(ratio = 1) +
  shared_heatmap_theme +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(5, 5, 0, 5)) +
  labs(subtitle = "Methanogen") +
  scale_fill_viridis_c(option = "C", name = expression(log[10]~"(Rel. Abundance + 1)"))

p_rel_mt <- ggplot(heatmap_rel_mt, aes(x = Media, y = Row_Label, fill = Abundance)) +
  geom_tile(color = "white") +
  coord_fixed(ratio = 1) +
  shared_heatmap_theme +
  theme(plot.margin = margin(0, 5, 5, 5)) +
  labs(y = "Taxa (Family | Genus)", subtitle = "Methanotroph") +
  scale_fill_viridis_c(option = "C", name = expression(log[10]~"(Rel. Abundance + 1)"))

p_rel_mg <- p_rel_mg + theme(legend.position = "none")
p_rel <- p_rel_mg / p_rel_mt +
  plot_layout(heights = c(n_rows_mg, n_rows_mt))

print(p_rel)

ggsave("outputs/figures/figS12_black_oak_methanome_relative.png",
       p_rel, width = 12, height = 8, dpi = 300)

# ==============================================================================
# STEP 5: Alpha Diversity Metrics
# ==============================================================================
library(vegan)

otu_table_long$OTU_ID <- otu_table_long$X

# Calculate alpha diversity for each sample and media type
diversity_data_otu <- otu_table_long %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Media) %>%
  summarize(
    Total_OTUs = n_distinct(OTU_ID),
    Chao1 = estimateR(Abundance)[1],
    Shannon = diversity(Abundance, index = "shannon"),
    Simpson = diversity(Abundance, index = "simpson"),
    .groups = "drop"
  )

# Total OTUs across media
total_otus_by_media <- otu_table_long %>%
  filter(Abundance > 0) %>%
  group_by(Media) %>%
  summarize(
    Total_OTUs = n_distinct(OTU_ID),
    .groups = "drop"
  )

# Plot alpha diversity
diversity_long_otu <- diversity_data_otu %>%
  pivot_longer(cols = c(Total_OTUs, Chao1, Shannon, Simpson),
               names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(diversity_long_otu, aes(x = Media, y = Value, fill = Media)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Media",
    y = "Diversity Metric Value",
    title = "Alpha Diversity Metrics for Methane-Cycling Taxa (OTU Level)"
  )

print(p_alpha)

# Plot total OTUs
p_total <- ggplot(total_otus_by_media, aes(x = Media, y = Total_OTUs, fill = Media)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Media",
    y = "Total OTUs",
    title = "Total Methane-Cycling OTUs Detected Across Tissues"
  )

print(p_total)
