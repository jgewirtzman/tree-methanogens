# ==============================================================================
# δ13CH4 Isotope Figures — Paired Heartwood/Sapwood (Dataset 2)
# ==============================================================================
# Purpose: Creates three figures from Picarro δ13CH4 isotope data for the paired
#   heartwood/sapwood sampling (H/S suffix samples). Annotated with methanogenesis
#   pathway isotope ranges.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - GC metadata: data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv
#   - Picarro Run 1: data/raw/internal_gas/picarro/20251128_205833_results.csv
#   - Picarro Run 2: data/raw/internal_gas/picarro/20251128_211030_results.csv
#
# Outputs:
#   - d13ch4_paired_overall.png
#   - d13ch4_paired_by_species.png
#   - d13ch4_paired_vs_concentration.png
#   - d13ch4_paired_rainfall.png
#   - d13ch4_paired_keeling.png
#
# Required packages: tidyverse, ggridges
# ==============================================================================

library(tidyverse)
library(ggridges)

# y-axis label for δ13CH4
# Note: ‰ (permil) not renderable in R plotmath; VPDB implies per mille
d13_label <- expression(delta^{13}*CH[4]~"(VPDB)")

# ==============================================================================
# STEP 1: Load GC metadata (master sample table)
# ==============================================================================

gc_meta <- read.csv("data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv",
                     fileEncoding = "UTF-8-BOM")

# Keep only biological samples and atmosphere
gc_samples <- gc_meta %>%
  filter(Sample.Type %in% c("Sample", "Atmosphere")) %>%
  select(Sample.ID, Species, Tree.No, Tissue, CH4_concentration)

# ==============================================================================
# STEP 2: Load and combine Picarro data from Runs 1 and 2
# ==============================================================================

pic1 <- read_csv("data/raw/internal_gas/picarro/20251128_205833_results.csv",
                 show_col_types = FALSE)
pic2 <- read_csv("data/raw/internal_gas/picarro/20251128_211030_results.csv",
                 show_col_types = FALSE)

# Run 1: all samples are paired H/S or atmosphere (no standards in this run)
# Run 2: keep only H/S-suffixed samples + Atm atmosphere (exclude standards and
#         single-per-tree samples that belong to Dataset 1)
pic2_paired <- pic2 %>%
  filter(str_detect(SampleName, "[HS]$") | str_detect(SampleName, "^Atm"))

# Combine
all_picarro <- bind_rows(pic1, pic2_paired)

# Exclude standards that may have H/S-like patterns
all_picarro <- all_picarro %>%
  filter(!SampleName %in% c("SB1", "SB3a", "SB3b", "SB3", "SB4a", "SB4",
                              "SB5a", "SB5", "S3a", "S3b", "S3c"))

# ==============================================================================
# STEP 3: Join Picarro isotope data to GC metadata
# ==============================================================================

iso_data <- all_picarro %>%
  left_join(gc_samples, by = c("SampleName" = "Sample.ID"))

# Check for unmatched samples
unmatched <- iso_data %>% filter(is.na(Species) & !str_detect(SampleName, "^Atm"))
if (nrow(unmatched) > 0) {
  cat("Warning: unmatched Picarro samples:\n")
  print(unmatched$SampleName)
}

# Classify atmosphere samples
iso_data <- iso_data %>%
  mutate(
    Species = case_when(
      str_detect(SampleName, "^Atm") ~ "Atmosphere",
      TRUE ~ Species
    ),
    Tissue = case_when(
      str_detect(SampleName, "^Atm") ~ "Atmosphere",
      TRUE ~ Tissue
    )
  )

# Drop unmatched (standards, unknowns)
iso_data <- iso_data %>% filter(!is.na(Species))

# ==============================================================================
# STEP 4: Extract key columns and clean
# ==============================================================================

iso_clean <- iso_data %>%
  select(SampleName, Species, Tree.No, Tissue,
         d13CH4 = HR_Delta_iCH4_Raw_mean,
         ch4_ppm_picarro = HR_12CH4_dry_mean,
         ch4_ppm_gc = CH4_concentration) %>%
  filter(!is.na(d13CH4)) %>%
  # Use coalesced concentration for filtering (GC preferred, Picarro fallback)
  mutate(ch4_ppm_filter = coalesce(ch4_ppm_gc, ch4_ppm_picarro)) %>%
  # Remove any sample below 1.5 ppm (unreliable isotope at very low conc.)
  filter(ch4_ppm_filter >= 1.5) %>%
  # Remove atmosphere samples above 5 ppm (contaminated with tree gas)
  filter(!(Species == "Atmosphere" & ch4_ppm_filter > 5)) %>%
  select(-ch4_ppm_filter)

# Species label mapping (genus abbreviations, matching Fig 5 convention)
species_labels <- c(
  "SM" = "A. saccharum", "BB" = "B. lenta",
  "WP" = "P. strobus", "H" = "T. canadensis",
  "RO" = "Q. rubra", "RM" = "A. rubrum",
  "Atmosphere" = "Atmosphere"
)

iso_clean <- iso_clean %>%
  mutate(species_label = species_labels[Species]) %>%
  filter(!is.na(species_label))

# Tissue labels for plotting
tissue_labels <- c("H" = "Heartwood", "S" = "Sapwood", "Atmosphere" = "Atmosphere")
iso_clean <- iso_clean %>%
  mutate(tissue_label = tissue_labels[Tissue])

cat("Dataset 2 summary:\n")
cat("Total samples:", nrow(iso_clean), "\n")
print(table(iso_clean$Species, iso_clean$Tissue))

# ==============================================================================
# STEP 5: Methanogenesis pathway annotation layers
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

# Color palette: Heartwood / Sapwood / Atmosphere
tissue_fill <- c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4", "Atmosphere" = "#7fbf7b")
tissue_color <- c("Heartwood" = "#8c510a", "Sapwood" = "#08519c", "Atmosphere" = "#4d9221")

# ==============================================================================
# FIGURE 1: Overall δ13CH4 distribution (Heartwood vs Sapwood vs Atmosphere)
# ==============================================================================

# Order tissue groups
iso_clean$tissue_label <- factor(iso_clean$tissue_label,
                                  levels = c("Heartwood", "Sapwood", "Atmosphere"))

n_groups <- length(levels(iso_clean$tissue_label))

p1 <- ggplot(iso_clean, aes(x = tissue_label, y = d13CH4)) +
  pathway_rects(0.4, n_groups + 0.6) +
  geom_violin(aes(fill = tissue_label), alpha = 0.4, color = NA) +
  geom_jitter(aes(color = tissue_label), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = tissue_fill) +
  scale_color_manual(values = tissue_color) +
  iso_theme +
  theme(legend.position = "none") +
  labs(x = NULL,
       y = d13_label) +
  pathway_labels(n_groups + 0.55)

print(p1)
# ggsave("outputs/figures/d13ch4_paired_overall.png",
#        p1, width = 6, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 2: δ13CH4 by species, colored by tissue
# ==============================================================================

# Order species alphabetically, atmosphere last
sp_order <- c(sort(unique(iso_clean$species_label[iso_clean$Species != "Atmosphere"])),
              "Atmosphere")
iso_clean$species_label <- factor(iso_clean$species_label, levels = sp_order)

n_species <- length(sp_order)

p2 <- ggplot(iso_clean, aes(x = species_label, y = d13CH4)) +
  pathway_rects(0.4, n_species + 0.6) +
  geom_boxplot(aes(fill = tissue_label), alpha = 0.4, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = tissue_label), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = tissue_fill, name = "Tissue") +
  scale_color_manual(values = tissue_color, name = "Tissue") +
  iso_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 10),
        legend.position = "bottom") +
  labs(x = NULL,
       y = d13_label) +
  pathway_labels(n_species + 0.5)

print(p2)
# ggsave("outputs/figures/d13ch4_paired_by_species.png",
#        p2, width = 10, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 3: CH4 concentration vs δ13CH4
# ==============================================================================

# Use GC CH4 concentration for paired samples; filter > 5 ppm for reliable isotopes
iso_conc <- iso_clean %>%
  filter(Species != "Atmosphere") %>%
  mutate(ch4_ppm = coalesce(ch4_ppm_gc, ch4_ppm_picarro)) %>%
  filter(ch4_ppm > 5)

# Get x-axis range for pathway rectangles
xrange <- range(log10(iso_conc$ch4_ppm), na.rm = TRUE)

p3 <- ggplot(iso_conc, aes(x = ch4_ppm, y = d13CH4)) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -110, ymax = -60, fill = "#4393c3", alpha = 0.12) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -65, ymax = -50, fill = "#d6604d", alpha = 0.12) +
  annotate("rect", xmin = 10^(xrange[1] - 0.2), xmax = 10^(xrange[2] + 0.2),
           ymin = -70, ymax = -50, fill = "#b2abd2", alpha = 0.12) +
  geom_point(aes(color = species_label, shape = tissue_label), size = 2.5, alpha = 0.7) +
  scale_x_log10() +
  scale_color_viridis_d(option = "D", name = "Species") +
  scale_shape_manual(values = c("Heartwood" = 16, "Sapwood" = 17), name = "Tissue") +
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
# ggsave("outputs/figures/d13ch4_paired_vs_concentration.png",
#        p3, width = 10, height = 7, dpi = 300)

# ==============================================================================
# FIGURE 4: Rainfall plot — δ13CH4 distribution with pathway brackets
# ==============================================================================

# Prepare factor ordering for y-axis
iso_rain <- iso_clean %>%
  mutate(tissue_label = factor(tissue_label,
                                levels = c("Heartwood", "Sapwood", "Atmosphere")))

# Pathway bracket positions below the data
# Groups at y=1 (Heartwood), y=2 (Sapwood), y=3 (Atmosphere); brackets below y=1
by <- 0.2    # base y for top bracket
bs <- 0.35   # spacing between bracket rows
bc <- 0.08   # end-cap height

p4 <- ggplot(iso_rain, aes(x = d13CH4, y = tissue_label, fill = tissue_label,
                            color = tissue_label)) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_raincloud(ygap = 0.05, adjust_vlines = TRUE),
    point_size = 2, point_alpha = 0.6,
    alpha = 0.4, scale = 0.7
  ) +
  # Hydrogenotrophic bracket: -110 to -60
  annotate("segment", x = -110, xend = -60, y = by, yend = by,
           color = "#2166ac", linewidth = 1.2) +
  annotate("segment", x = -110, xend = -110, y = by, yend = by + bc,
           color = "#2166ac", linewidth = 0.7) +
  annotate("segment", x = -60, xend = -60, y = by, yend = by + bc,
           color = "#2166ac", linewidth = 0.7) +
  annotate("text", x = -85, y = by - 0.08, label = "Hydrogenotrophic",
           color = "#2166ac", size = 3, fontface = "italic") +
  # Acetoclastic bracket: -65 to -50
  annotate("segment", x = -65, xend = -50, y = by - bs, yend = by - bs,
           color = "#b2182b", linewidth = 1.2) +
  annotate("segment", x = -65, xend = -65, y = by - bs, yend = by - bs + bc,
           color = "#b2182b", linewidth = 0.7) +
  annotate("segment", x = -50, xend = -50, y = by - bs, yend = by - bs + bc,
           color = "#b2182b", linewidth = 0.7) +
  annotate("text", x = -57.5, y = by - bs - 0.08, label = "Acetoclastic",
           color = "#b2182b", size = 3, fontface = "italic") +
  # Methylotrophic bracket: -70 to -50
  annotate("segment", x = -70, xend = -50, y = by - 2*bs, yend = by - 2*bs,
           color = "#7a5195", linewidth = 1.2) +
  annotate("segment", x = -70, xend = -70, y = by - 2*bs, yend = by - 2*bs + bc,
           color = "#7a5195", linewidth = 0.7) +
  annotate("segment", x = -50, xend = -50, y = by - 2*bs, yend = by - 2*bs + bc,
           color = "#7a5195", linewidth = 0.7) +
  annotate("text", x = -60, y = by - 2*bs - 0.08, label = "Methylotrophic",
           color = "#7a5195", size = 3, fontface = "italic") +
  scale_fill_manual(values = tissue_fill) +
  scale_color_manual(values = tissue_color) +
  coord_cartesian(clip = "off") +
  iso_theme +
  theme(legend.position = "none",
        plot.margin = margin(t = 5, r = 10, b = 60, l = 5, unit = "pt")) +
  labs(x = d13_label,
       y = NULL)

print(p4)
# ggsave("outputs/figures/d13ch4_paired_rainfall.png",
#        p4, width = 10, height = 6, dpi = 300)

# ==============================================================================
# FIGURE 5: Keeling plot — 1/[CH4] vs δ13CH4
# ==============================================================================

# Prepare Keeling data from internal samples (already filtered to ch4 > 5 ppm)
keeling_data <- iso_conc %>%
  mutate(inv_ch4 = 1 / ch4_ppm)

# Atmosphere reference points
atm_keeling <- iso_clean %>%
  filter(Species == "Atmosphere") %>%
  mutate(ch4_ppm = coalesce(ch4_ppm_gc, ch4_ppm_picarro)) %>%
  filter(ch4_ppm > 0) %>%
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
  # Internal samples: color = species, shape = tissue
  geom_point(data = keeling_data,
             aes(x = inv_ch4, y = d13CH4,
                 color = species_label, shape = tissue_label),
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
  # Tissue shapes
  scale_shape_manual(values = c("Heartwood" = 16, "Sapwood" = 17), name = "Tissue") +
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
# ggsave("outputs/figures/d13ch4_paired_keeling.png",
#        p5, width = 10, height = 7, dpi = 300)

cat("Dataset 2 figures saved.\n")
