# ==============================================================================
# δ13CH4 Isotope Figures — Single Per Tree (Dataset 1)
# ==============================================================================
# Purpose: Creates three figures from Picarro δ13CH4 isotope data for the main
#   study trees (one sample per tree). Annotated with methanogenesis pathway
#   isotope ranges.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - Picarro results: data/raw/internal_gas/picarro/20251128_211030_results.csv
#   - Picarro results: data/raw/internal_gas/picarro/20251128_213226_results.csv
#   - Picarro results: data/raw/internal_gas/picarro/20251128_215521_results.csv
#   - ddPCR metadata:  data/raw/ddpcr/ddPCR_meta_all_data.csv (for species mapping)
#
# Outputs:
#   - d13ch4_single_overall.png
#   - d13ch4_single_by_species.png
#   - d13ch4_single_vs_concentration.png
#   - d13ch4_single_rainfall.png
#   - d13ch4_single_keeling.png
#
# Required packages: tidyverse, ggridges
# ==============================================================================

library(tidyverse)
library(ggridges)

# y-axis label for δ13CH4
# Note: ‰ (permil) not renderable in R plotmath; VPDB implies per mille
d13_label <- expression(delta^{13}*CH[4]~"(VPDB)")

# ==============================================================================
# STEP 1: Load and combine Picarro data from Runs 2, 3, 4
# ==============================================================================

pic2 <- read_csv("data/raw/internal_gas/picarro/20251128_211030_results.csv",
                 show_col_types = FALSE)
pic3 <- read_csv("data/raw/internal_gas/picarro/20251128_213226_results.csv",
                 show_col_types = FALSE)
pic4 <- read_csv("data/raw/internal_gas/picarro/20251128_215521_results.csv",
                 show_col_types = FALSE)

# Run 2: keep only samples WITHOUT H/S suffix (single per tree) + Amb atmosphere
pic2_single <- pic2 %>%
  filter(!str_detect(SampleName, "[HS]$")) %>%
  filter(!SampleName %in% c("SB1", "SB3a", "SB3b", "SB3", "SB4a", "SB4",
                              "SB5a", "SB5", "S3a", "S3b", "S3c"))

# Run 3: exclude standards, keep all tree + atmosphere samples
pic3_single <- pic3 %>%
  filter(!SampleName %in% c("SB1", "SB3a", "SB3b", "SB3", "SB4a", "SB4",
                              "SB5a", "SB5", "S3a", "S3b", "S3c"))

# Run 4: keep only RM trees (exclude V_ vial experiments and control samples)
pic4_single <- pic4 %>%
  filter(!str_detect(SampleName, "^V|^FLUSH|^RINS|^ROOM"))

# Combine all
all_picarro <- bind_rows(pic2_single, pic3_single, pic4_single)

# ==============================================================================
# STEP 2: Map sample names to species via ddPCR metadata
# ==============================================================================

ddpcr <- read_csv("data/raw/ddpcr/ddPCR_meta_all_data.csv", show_col_types = FALSE)
species_map <- ddpcr %>%
  distinct(seq_id, species) %>%
  rename(SampleName = seq_id)

# Join species
iso_data <- all_picarro %>%
  left_join(species_map, by = "SampleName")

# Classify atmosphere samples
iso_data <- iso_data %>%
  mutate(species = case_when(
    str_detect(SampleName, "^Amb") ~ "Atmosphere",
    TRUE ~ species
  ))

# Drop samples with no species match (standards, unknowns)
iso_data <- iso_data %>% filter(!is.na(species))

# ==============================================================================
# STEP 3: Extract key columns and clean
# ==============================================================================

iso_clean <- iso_data %>%
  select(SampleName, species,
         d13CH4 = HR_Delta_iCH4_Raw_mean,
         ch4_ppm = HR_12CH4_dry_mean) %>%
  filter(!is.na(d13CH4)) %>%
  # Remove any sample below 1.5 ppm (unreliable isotope at very low conc.)
  filter(ch4_ppm >= 1.5) %>%
  # Remove atmosphere samples above 5 ppm (contaminated with tree gas)
  filter(!(species == "Atmosphere" & ch4_ppm > 5))

# Species label mapping (genus abbreviations)
species_labels <- c(
  "ACRU" = "A. rubrum", "ACSA" = "A. saccharum",
  "BEAL" = "B. alleghaniensis", "BELE" = "B. lenta",
  "BEPA" = "B. papyrifera", "CAOV" = "C. ovata",
  "FAGR" = "F. grandifolia", "FRAM" = "F. americana",
  "KALA" = "K. latifolia", "PIST" = "P. strobus",
  "PRSE" = "P. serotina", "QUAL" = "Q. alba",
  "QURU" = "Q. rubra", "QUVE" = "Q. velutina",
  "SAAL" = "S. albidum", "TSCA" = "T. canadensis",
  "Atmosphere" = "Atmosphere"
)

iso_clean <- iso_clean %>%
  mutate(species_label = species_labels[species]) %>%
  filter(!is.na(species_label))

# Classify as internal vs atmosphere
iso_clean <- iso_clean %>%
  mutate(sample_type = if_else(species == "Atmosphere", "Atmosphere", "Internal"))

cat("Dataset 1 summary:\n")
cat("Total samples:", nrow(iso_clean), "\n")
print(table(iso_clean$species))

# ==============================================================================
# STEP 4: Methanogenesis pathway annotation layers
# ==============================================================================

# Standard literature ranges (δ13CH4, ‰ VPDB)
# Hydrogenotrophic: −110 to −60‰
# Acetoclastic:     −65 to −50‰
# Methylotrophic:   −70 to −50‰

pathway_rects <- function(xmin_val, xmax_val) {
  list(
    annotate("rect", xmin = xmin_val, xmax = xmax_val,
             ymin = -110, ymax = -60,
             fill = "#4393c3", alpha = 0.12),
    annotate("rect", xmin = xmin_val, xmax = xmax_val,
             ymin = -65, ymax = -50,
             fill = "#d6604d", alpha = 0.12),
    annotate("rect", xmin = xmin_val, xmax = xmax_val,
             ymin = -70, ymax = -50,
             fill = "#b2abd2", alpha = 0.12)
  )
}

pathway_labels <- function(x_pos) {
  list(
    annotate("text", x = x_pos, y = -85, label = "Hydrogenotrophic",
             size = 3, fontface = "italic", color = "#2166ac", hjust = 1),
    annotate("text", x = x_pos, y = -56, label = "Acetoclastic",
             size = 3, fontface = "italic", color = "#b2182b", hjust = 1),
    annotate("text", x = x_pos, y = -72, label = "Methylotrophic",
             size = 3, fontface = "italic", color = "#7a5195", hjust = 1)
  )
}

# Shared theme
iso_theme <- theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

# ==============================================================================
# FIGURE 1: Overall δ13CH4 distribution (Internal vs Atmosphere)
# ==============================================================================

p1 <- ggplot(iso_clean, aes(x = sample_type, y = d13CH4)) +
  pathway_rects(0.4, 2.6) +
  geom_violin(aes(fill = sample_type), alpha = 0.4, color = NA) +
  geom_jitter(aes(color = sample_type), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("Atmosphere" = "#7fbf7b", "Internal" = "#af8dc3")) +
  scale_color_manual(values = c("Atmosphere" = "#4d9221", "Internal" = "#762a83")) +
  iso_theme +
  theme(legend.position = "none") +
  labs(x = NULL,
       y = d13_label) +
  pathway_labels(2.55)

print(p1)
# ggsave("outputs/figures/d13ch4_single_overall.png",
#        p1, width = 6, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 2: δ13CH4 by species
# ==============================================================================

# Order species alphabetically, atmosphere last
sp_order <- c(sort(unique(iso_clean$species_label[iso_clean$species != "Atmosphere"])),
              "Atmosphere")
iso_clean$species_label <- factor(iso_clean$species_label, levels = sp_order)

n_species <- length(sp_order)

p2 <- ggplot(iso_clean, aes(x = species_label, y = d13CH4)) +
  pathway_rects(0.4, n_species + 0.6) +
  geom_boxplot(aes(fill = sample_type), alpha = 0.4, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = sample_type), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("Atmosphere" = "#7fbf7b", "Internal" = "#af8dc3")) +
  scale_color_manual(values = c("Atmosphere" = "#4d9221", "Internal" = "#762a83")) +
  iso_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 10),
        legend.position = "none") +
  labs(x = NULL,
       y = d13_label) +
  pathway_labels(n_species + 0.5)

print(p2)
# ggsave("outputs/figures/d13ch4_single_by_species.png",
#        p2, width = 14, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 3: CH4 concentration vs δ13CH4
# ==============================================================================

# Filter to samples with reliable isotope data (CH4 > 5 ppm)
iso_conc <- iso_clean %>%
  filter(species != "Atmosphere", ch4_ppm > 5)

# Get x-axis range for pathway rectangles
xrange <- range(log10(iso_conc$ch4_ppm), na.rm = TRUE)

p3 <- ggplot(iso_conc, aes(x = ch4_ppm, y = d13CH4)) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -110, ymax = -60, fill = "#4393c3", alpha = 0.12) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -65, ymax = -50, fill = "#d6604d", alpha = 0.12) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -70, ymax = -50, fill = "#b2abd2", alpha = 0.12) +
  geom_point(aes(color = species_label), size = 2.5, alpha = 0.7) +
  scale_x_log10() +
  scale_color_viridis_d(option = "D", name = "Species") +
  iso_theme +
  theme(legend.text = element_text(face = "italic", size = 9),
        legend.title = element_text(size = 10)) +
  labs(x = expression(CH[4]~"(ppm)"),
       y = d13_label) +
  annotate("text", x = 10^(xrange[2] + 0.15), y = -85,
           label = "Hydrogenotrophic", size = 3, fontface = "italic",
           color = "#2166ac", hjust = 1) +
  annotate("text", x = 10^(xrange[2] + 0.15), y = -56,
           label = "Acetoclastic", size = 3, fontface = "italic",
           color = "#b2182b", hjust = 1) +
  annotate("text", x = 10^(xrange[2] + 0.15), y = -72,
           label = "Methylotrophic", size = 3, fontface = "italic",
           color = "#7a5195", hjust = 1)

print(p3)
# ggsave("outputs/figures/d13ch4_single_vs_concentration.png",
#        p3, width = 10, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 4: Rainfall plot — δ13CH4 distribution with pathway brackets
# ==============================================================================

# Internal samples only (drop atmosphere)
iso_rain <- iso_clean %>%
  filter(sample_type == "Internal")

# Compute kernel density manually so we control the height
dens <- density(iso_rain$d13CH4, na.rm = TRUE)
dens_df <- data.frame(x = dens$x, y = dens$y)
# Scale density so max height fills a nice region (peak at y ~ 0.55)
dens_df$y_scaled <- dens_df$y / max(dens_df$y) * 0.55

# Jittered points: random y between -0.05 and -0.28
set.seed(42)
iso_rain$jitter_y <- runif(nrow(iso_rain), -0.28, -0.05)

# Pathway bracket positions (below jittered points, in clip-off margin)
by <- -0.45  # base y for top bracket
bs <- 0.16   # spacing between bracket rows
bc <- 0.04   # end-cap height

p4 <- ggplot() +
  # Atmospheric δ13CH4 reference line (behind all data layers)
  geom_vline(xintercept = -47, linetype = "dashed", color = "gray30", linewidth = 0.5) +
  annotate("text", x = -47, y = 0.53, label = "Atmosphere", color = "gray30",
           size = 3, fontface = "italic", hjust = -0.1) +
  # Density ridge (filled polygon above y=0)
  geom_ribbon(data = dens_df, aes(x = x, ymin = 0, ymax = y_scaled),
              fill = "#af8dc3", color = "#762a83", alpha = 0.5, linewidth = 0.6) +
  # Jittered points below
  geom_point(data = iso_rain, aes(x = d13CH4, y = jitter_y),
             color = "#762a83", fill = "#af8dc3", shape = 21,
             size = 2, alpha = 0.6, stroke = 0.3) +
  # Hydrogenotrophic bracket: -110 to -60
  annotate("segment", x = -110, xend = -60, y = by, yend = by,
           color = "gray30", linewidth = 1) +
  annotate("segment", x = -110, xend = -110, y = by, yend = by + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("segment", x = -60, xend = -60, y = by, yend = by + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("text", x = -85, y = by - 0.04, label = "Hydrogenotrophic",
           color = "gray30", size = 3, fontface = "italic") +
  # Acetoclastic bracket: -65 to -50
  annotate("segment", x = -65, xend = -50, y = by - bs, yend = by - bs,
           color = "gray30", linewidth = 1) +
  annotate("segment", x = -65, xend = -65, y = by - bs, yend = by - bs + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("segment", x = -50, xend = -50, y = by - bs, yend = by - bs + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("text", x = -57.5, y = by - bs - 0.04, label = "Acetoclastic",
           color = "gray30", size = 3, fontface = "italic") +
  # Methylotrophic bracket: -70 to -50
  annotate("segment", x = -70, xend = -50, y = by - 2*bs, yend = by - 2*bs,
           color = "gray30", linewidth = 1) +
  annotate("segment", x = -70, xend = -70, y = by - 2*bs, yend = by - 2*bs + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("segment", x = -50, xend = -50, y = by - 2*bs, yend = by - 2*bs + bc,
           color = "gray30", linewidth = 0.6) +
  annotate("text", x = -60, y = by - 2*bs - 0.04, label = "Methylotrophic",
           color = "gray30", size = 3, fontface = "italic") +
  # Zoom x-axis to show core data + pathway ranges; outliers beyond are clipped
  coord_cartesian(xlim = c(-125, 40), clip = "off") +
  iso_theme +
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(t = 5, r = 10, b = 55, l = 5, unit = "pt")) +
  labs(x = d13_label, y = NULL)

print(p4)
ggsave("outputs/figures/supplementary/figS11_d13ch4_rainfall.png",
       p4, width = 10, height = 5, dpi = 300)

# ==============================================================================
# FIGURE 5: Keeling plot — 1/[CH4] vs δ13CH4
# ==============================================================================

# Prepare Keeling data from internal samples (already filtered to ch4 > 5 ppm)
keeling_data <- iso_conc %>%
  mutate(inv_ch4 = 1 / ch4_ppm)

# Atmosphere reference points
atm_keeling <- iso_clean %>%
  filter(species == "Atmosphere", ch4_ppm > 0) %>%
  mutate(inv_ch4 = 1 / ch4_ppm)

# Fit overall linear model
keeling_lm <- lm(d13CH4 ~ inv_ch4, data = keeling_data)
keeling_intercept <- coef(keeling_lm)[1]
keeling_r2 <- summary(keeling_lm)$r.squared

# x-range for pathway background rectangles
xrange_k <- range(c(keeling_data$inv_ch4, atm_keeling$inv_ch4), na.rm = TRUE)
xpad <- diff(xrange_k) * 0.05

# Intercept annotation text
intercept_text <- paste0("y-int = ", round(keeling_intercept, 1),
                          "  (R\u00B2 = ", round(keeling_r2, 3), ")")

p5 <- ggplot() +
  # Pathway shaded bands
  annotate("rect", xmin = xrange_k[1] - xpad, xmax = xrange_k[2] + xpad,
           ymin = -110, ymax = -60, fill = "#4393c3", alpha = 0.10) +
  annotate("rect", xmin = xrange_k[1] - xpad, xmax = xrange_k[2] + xpad,
           ymin = -65, ymax = -50, fill = "#d6604d", alpha = 0.10) +
  annotate("rect", xmin = xrange_k[1] - xpad, xmax = xrange_k[2] + xpad,
           ymin = -70, ymax = -50, fill = "#b2abd2", alpha = 0.10) +
  # Internal samples colored by species
  geom_point(data = keeling_data,
             aes(x = inv_ch4, y = d13CH4, color = species_label),
             size = 2.5, alpha = 0.7) +
  # Atmosphere reference points (open diamonds)
  geom_point(data = atm_keeling,
             aes(x = inv_ch4, y = d13CH4),
             shape = 5, size = 3, color = "#4d9221", stroke = 1.2) +
  # Overall regression line with CI
  geom_smooth(data = keeling_data,
              aes(x = inv_ch4, y = d13CH4),
              method = "lm", se = TRUE,
              color = "black", fill = "gray80", alpha = 0.3,
              linewidth = 0.8) +
  # Species color scale
  scale_color_viridis_d(option = "D", name = "Species") +
  # Y-intercept marker on y-axis
  annotate("point", x = 0, y = keeling_intercept,
           shape = 18, size = 4, color = "black") +
  # Y-intercept annotation
  annotate("text", x = xrange_k[2] * 0.5, y = keeling_intercept + 8,
           label = intercept_text,
           size = 3.5, fontface = "bold") +
  # Pathway labels
  annotate("text", x = xrange_k[2] + xpad, y = -85,
           label = "Hydrogenotrophic", size = 2.8, fontface = "italic",
           color = "#2166ac", hjust = 1) +
  annotate("text", x = xrange_k[2] + xpad, y = -56,
           label = "Acetoclastic", size = 2.8, fontface = "italic",
           color = "#b2182b", hjust = 1) +
  annotate("text", x = xrange_k[2] + xpad, y = -72,
           label = "Methylotrophic", size = 2.8, fontface = "italic",
           color = "#7a5195", hjust = 1) +
  iso_theme +
  theme(legend.text = element_text(face = "italic", size = 9),
        legend.title = element_text(size = 10)) +
  labs(x = expression(1/CH[4]~"(ppm"^{-1}*")"),
       y = d13_label)

print(p5)
# ggsave("outputs/figures/d13ch4_single_keeling.png",
#        p5, width = 10, height = 7, dpi = 300)

cat("Dataset 1 figures saved.\n")
