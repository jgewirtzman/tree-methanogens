# ==============================================================================
# Radial Wood Cross-Section Plots (Figure 8, S9)
# ==============================================================================
# Purpose: Publication-quality radial wood cross-section plots showing spatial
#   gene distribution.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - flux data (from data/processed/flux/)
#
# Outputs:
#   - tree_mcra_cross_sections_compressed.pdf (to outputs/figures/)
# ==============================================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggforce)
  library(patchwork)
})

# ---- Load data ----
# Update these paths as needed
ymf2023 <- read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")

merged_final <- ymf2021  # Adjust if ddPCR data is elsewhere

# ---- Species mapping for labels ----
species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "SAAL" = "Sassafras albidum"
)

# ============================================================
# 1) Parse ddPCR data into long format
# ============================================================
prepare_long_mcra <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    pivot_longer(
      cols = starts_with("ddpcr_"),
      names_to = "measurement_type",
      values_to = "gene_copies"
    ) %>%
    filter(!is.na(gene_copies)) %>%
    separate(
      measurement_type,
      into = c("method", "gene", "part1", "part2", "part3"),
      sep = "_", extra = "merge", fill = "right"
    ) %>%
    mutate(
      is_probe = (part1 == "probe"),
      location = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(
        location %in% c("Inner","inner") ~ "Inner",
        location %in% c("Outer","outer") ~ "Outer",
        location %in% c("Mineral","mineral") ~ "Mineral",
        location %in% c("Organic","organic") ~ "Organic",
        TRUE ~ location
      ),
      gene = case_when(
        gene == "mcra" ~ "mcrA",
        gene == "mmox" ~ "mmoX",
        gene == "pmoa" ~ "pmoA",
        TRUE ~ gene
      ),
      sample_type = case_when(
        location == "Inner" ~ "Heartwood",
        location == "Outer" ~ "Sapwood",
        location == "Mineral" ~ "Mineral",
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene == "mcrA", is_probe) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

long_data <- prepare_long_mcra(merged_final)

# ============================================================
# 2) Area-weighted function (closed-form, no density)
# ============================================================
area_weighted_mcra <- function(mcra_inner, mcra_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  mcra_outer + (mcra_inner - mcra_outer) * (S / (3 * R^2))
}

# Vectorized version for use in mutate
area_weighted_mcra_vec <- Vectorize(area_weighted_mcra)

# ============================================================
# 3) Calculate per-tree metrics
# ============================================================
mcra_with_species <- long_data %>%
  left_join(
    merged_final %>% dplyr::select(tree_id, species_id),
    by = "tree_id"
  ) %>%
  mutate(species = unname(species_mapping[species_id])) %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, dbh),
    by = "tree_id"
  ) %>%
  filter(is.finite(dbh))

# Calculate inner/outer means and area-weighted value per tree
tree_weighted <- mcra_with_species %>%
  group_by(tree_id, species_id, species, dbh) %>%
  summarise(
    mcra_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    mcra_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mcra_inner), is.finite(mcra_outer), is.finite(dbh)) %>%
  mutate(
    mcra_area_weighted = area_weighted_mcra_vec(mcra_inner, mcra_outer, dbh),
    species_label = ifelse(is.na(species) | species == "", species_id, species)
  )

# ============================================================
# 4) Summary statistics
# ============================================================
cat("\n=============================\n")
cat("AREA-WEIGHTED mcrA SUMMARY\n")
cat("=============================\n")
cat(sprintf("Total trees with complete data: %d\n", nrow(tree_weighted)))
cat(sprintf("Unique species: %d\n", n_distinct(tree_weighted$species_id)))
cat(sprintf("Median log10(mcrA+1): %.3f\n",
            median(log10(tree_weighted$mcra_area_weighted + 1), na.rm = TRUE)))

# Species counts
species_summary <- tree_weighted %>%
  group_by(species_label) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcra_area_weighted, na.rm = TRUE)
  ) %>%
  arrange(desc(n_trees))

cat("\nSpecies representation:\n")
print(species_summary)

# ============================================================
# 5) Final visualization: All trees by species
# ============================================================

# Filter for species with enough trees (adjust threshold as needed)
min_trees_per_species <- 3  # Changed to 3 as requested

species_to_include <- tree_weighted %>%
  group_by(species_label) %>%
  summarise(n_trees = n()) %>%
  filter(n_trees >= min_trees_per_species) %>%
  pull(species_label)

# Prepare trees for plotting - order species by mean area-weighted mcrA
trees_to_plot <- tree_weighted %>%
  filter(species_label %in% species_to_include) %>%
  group_by(species_label) %>%
  mutate(
    species_mean_mcra = mean(mcra_area_weighted, na.rm = TRUE)
  ) %>%
  arrange(desc(species_mean_mcra), mcra_area_weighted) %>%
  mutate(tree_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    # Order species from highest to lowest mean mcrA
    species_label = reorder(species_label, -species_mean_mcra)
  )

cat(sprintf("\nTrees included in visualization: %d\n", nrow(trees_to_plot)))
cat(sprintf("Species included: %d\n", n_distinct(trees_to_plot$species_label)))

# Generate ring cloud data with consistent grid spacing
generate_tree_grid <- function(tree_row, grid_resolution = 0.5) {
  R <- tree_row$dbh / 2
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  # Create rectangular grid
  x_seq <- seq(-R, R, by = grid_resolution)
  y_seq <- seq(-R, R, by = grid_resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  # Calculate radius for each point
  grid$r <- sqrt(grid$x^2 + grid$y^2)
  
  # Keep only points inside circle
  grid <- grid[grid$r <= R, ]
  
  # First convert to log scale for both inner and outer values
  log_inner <- log10(tree_row$mcra_inner + 1)
  log_outer <- log10(tree_row$mcra_outer + 1)
  
  # Calculate log-scale interpolation
  grid$Clog <- ifelse(grid$r <= r1, log_inner,
                      ifelse(grid$r >= r2, log_outer,
                             log_inner + (log_outer - log_inner) *
                               (grid$r - r1) / max(r2 - r1, 1e-9)))
  
  # Back-calculate C for reference (if needed)
  grid$C <- 10^grid$Clog - 1
  
  grid$species_label <- tree_row$species_label
  grid$tree_rank <- tree_row$tree_rank
  grid$tree_id <- tree_row$tree_id
  grid$dbh <- tree_row$dbh
  grid$r1 <- r1
  grid$r2 <- r2
  grid$R <- R
  
  grid
}

# Generate data for all trees
all_cloud <- bind_rows(lapply(1:nrow(trees_to_plot), function(i) {
  generate_tree_grid(trees_to_plot[i,])
}))

# Main visualization
p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression("log"[10]*"(mcrA)"),
    option = "plasma"
  ) +
  geom_circle(
    data = all_cloud %>% distinct(species_label, tree_rank, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.3
  ) +
  facet_grid(species_label ~ tree_rank, switch = "both") +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "italic", size = 8),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.05, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  ) #+
  #labs(
   # title = "mcrA distribution across tree cross-sections",
    #subtitle = sprintf("Rows = species (n≥%d trees), columns = individual trees (ordered by area-weighted mcrA)", 
     #                  min_trees_per_species)
  #)

print(p_main)
# 
# # ============================================================
# # 6) Optional: Summary plot showing species-level patterns
# # ============================================================
# species_summary_plot <- tree_weighted %>%
#   filter(species_label %in% species_to_include) %>%
#   mutate(
#     # Use same ordering as main plot
#     species_label = factor(species_label, 
#                            levels = rev(levels(trees_to_plot$species_label)))
#   ) %>%
#   ggplot(aes(x = species_label, 
#              y = log10(mcra_area_weighted + 1))) +
#   geom_boxplot(aes(fill = species_label), alpha = 0.7, show.legend = FALSE) +
#   geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
#   scale_fill_viridis_d(option = "plasma") +
#   coord_flip() +
#   theme_minimal(base_size = 10) +
#   labs(
#     title = "Area-weighted mcrA by species",
#     subtitle = "Species ordered by mean mcrA (highest to lowest)",
#     x = "",
#     y = expression("log"[10]*"(mcrA + 1) copies/g")
#   )
# 
# print(species_summary_plot)

# ============================================================
# 7) Export results
# ============================================================
# write.csv(tree_weighted, "tree_mcra_area_weighted_clean.csv", row.names = FALSE)
# ggsave("tree_mcra_cross_sections.pdf", p_main, width = 12, height = 8)




# Check for trees with higher sapwood than heartwood
higher_sapwood <- tree_weighted %>%
  filter(mcra_outer > mcra_inner) %>%
  dplyr::select(tree_id, species_label, mcra_inner, mcra_outer, dbh) %>%
  mutate(
    ratio = mcra_outer / mcra_inner,
    diff = mcra_outer - mcra_inner
  )

# Summary
cat(sprintf("Trees with higher sapwood mcrA: %d out of %d (%.1f%%)\n", 
            nrow(higher_sapwood), 
            nrow(tree_weighted),
            100 * nrow(higher_sapwood) / nrow(tree_weighted)))

# Show some examples
if(nrow(higher_sapwood) > 0) {
  print(higher_sapwood %>% arrange(desc(ratio)) %>% head(10))
  
  # By species
  cat("\nBy species:\n")
  higher_sapwood %>%
    group_by(species_label) %>%
    summarise(
      n = n(),
      median_ratio = median(ratio),
      max_ratio = max(ratio)
    ) %>%
    arrange(desc(n)) %>%
    print()
}


# Main visualization with reduced vertical spacing
p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression("log"[10]*"(mcrA)"),
    option = "plasma"
  ) +
  geom_circle(
    data = all_cloud %>% distinct(species_label, tree_rank, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.2  # Thinner borders
  ) +
  facet_grid(species_label ~ tree_rank, switch = "both") +
  coord_equal() +
  theme_minimal(base_size = 8) +  # Smaller base size
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    # Reduce species label size and margins
    strip.text.y.left = element_text(
      angle = 0, 
      face = "italic", 
      size = 6.5,  # Smaller font
      margin = margin(0, 0.5, 0, 0.5)  # Very minimal margins
    ),
    strip.text.x = element_text(
      size = 5,  # Smaller column headers
      margin = margin(0.5, 0, 0.5, 0)  # Minimal margins
    ),
    strip.background = element_rect(fill = NA, color = NA),  # Remove strip backgrounds
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),  # Remove legend margins
    legend.box.margin = margin(-3, 0, 0, 0),  # Pull legend closer
    # Key adjustment: small positive spacing to avoid cutoff
    panel.spacing.y = unit(0.02, "lines"),  # Very small spacing
    panel.spacing.x = unit(0.01, "lines"),  # Minimal horizontal spacing
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1)  # Minimal plot margins
  )

print(p_main)

# Optional: If you want even more compression, save with adjusted dimensions
# This creates a wider, shorter plot that emphasizes the grid layout
ggsave("outputs/figures/supplementary/figS11_tree_radial_sections.pdf",
       p_main,
       width = 12,
       height = 5.1,
       dpi = 300)



# Main visualization with reduced vertical spacing
p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression("log"[10]*"(mcrA)"),
    option = "plasma"
  ) +
  geom_circle(
    data = all_cloud %>% distinct(species_label, tree_rank, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.2  # Thinner borders
  ) +
  facet_grid(species_label ~ tree_rank, switch = "both") +
  coord_equal() +
  theme_minimal(base_size = 8) +  # Smaller base size
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    # Reduce species label size and margins
    strip.text.y.left = element_text(
      angle = 0, 
      face = "italic", 
      size = 6.5,  # Smaller font
      margin = margin(0, 0.5, 0, 0.5)  # Very minimal margins
    ),
    strip.text.x = element_text(
      size = 5,  # Smaller column headers
      margin = margin(0.5, 0, 0.5, 0)  # Minimal margins
    ),
    strip.background = element_rect(fill = NA, color = NA),  # Remove strip backgrounds
    legend.position = c(0.98, 0.98),  # Top right corner
    legend.justification = c("right", "top"),  # Anchor to top-right
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),  # White background with border
    legend.margin = margin(2, 2, 2, 2),  # Small margin inside legend box
    legend.box.margin = margin(0, 0, 0, 0),  # No extra margin around legend
    # Key adjustment: small positive spacing to avoid cutoff
    panel.spacing.y = unit(0.02, "lines"),  # Very small spacing
    panel.spacing.x = unit(0.01, "lines"),  # Minimal horizontal spacing
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1)  # Minimal plot margins
  )
print(p_main)


# Main visualization with reduced vertical spacing
p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression("log"[10]*"(mcrA)"),
    option = "plasma"
  ) +
  geom_circle(
    data = all_cloud %>% distinct(species_label, tree_rank, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.2  # Thinner borders
  ) +
  facet_grid(species_label ~ tree_rank, switch = "both") +
  coord_equal() +
  theme_minimal(base_size = 8) +  # Smaller base size
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    # Increase species label size
    strip.text.y.left = element_text(
      angle = 0, 
      face = "italic", 
      size = 8,  # Increased from 6.5 to 8
      margin = margin(0, 0.5, 0, 0.5)  # Very minimal margins
    ),
    strip.text.x = element_blank(),  # Remove the top labels (numbers)
    strip.background = element_rect(fill = NA, color = NA),  # Remove strip backgrounds
    legend.position = c(0.98, 0.98),  # Top right corner
    legend.justification = c("right", "top"),  # Anchor to top-right
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),  # White background with border
    legend.margin = margin(2, 2, 2, 2),  # Small margin inside legend box
    legend.box.margin = margin(0, 0, 0, 0),  # No extra margin around legend
    # Key adjustment: small positive spacing to avoid cutoff
    panel.spacing.y = unit(0.02, "lines"),  # Very small spacing
    panel.spacing.x = unit(0.01, "lines"),  # Minimal horizontal spacing
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1)  # Minimal plot margins
  )
print(p_main)

# 
# # Main visualization with true-size DBH legend
# p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
#   geom_raster(interpolate = TRUE) +
#   scale_fill_viridis_c(
#     name = expression("log"[10]*"(mcrA)"),
#     option = "plasma"
#   ) +
#   # Keep original circles WITHOUT size aesthetic
#   geom_circle(
#     data = all_cloud %>% distinct(species_label, tree_rank, R),
#     aes(x0 = 0, y0 = 0, r = R),
#     inherit.aes = FALSE, color = "black", linewidth = 0.2
#   ) +
#   # Add invisible points for legend with area scaling
#   geom_point(
#     data = data.frame(x = NA, y = NA, dbh = c(20, 30, 40, 50, 60)),
#     aes(x = x, y = y, size = dbh),
#     na.rm = TRUE, shape = 21, fill = NA, color = "black", stroke = 0.5
#   ) +
#   # Use scale_size_area for true proportional circles
#   scale_size_area(
#     name = "DBH (cm)",
#     breaks = c(20, 30, 40, 50, 60),
#     labels = c("20", "30", "40", "50", "60"),
#     max_size = 12  # Adjust this to match your actual circle sizes
#   ) +
#   facet_grid(species_label ~ tree_rank, switch = "both") +
#   coord_equal() +
#   theme_minimal(base_size = 8) +
#   theme(
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank(),
#     strip.text.y.left = element_text(
#       angle = 0, 
#       face = "italic", 
#       size = 8,
#       margin = margin(0, 0.5, 0, 0.5)
#     ),
#     strip.text.x = element_blank(),
#     strip.background = element_rect(fill = NA, color = NA),
#     legend.position = c(0.98, 0.98),
#     legend.justification = c("right", "top"),
#     legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
#     legend.margin = margin(2, 2, 2, 2),
#     legend.box.margin = margin(0, 0, 0, 0),
#     legend.box = "horizontal",
#     legend.spacing.x = unit(0.3, "cm"),
#     legend.key.size = unit(1.2, "cm"),  # Larger to accommodate biggest circle
#     legend.title = element_text(size = 7),
#     legend.text = element_text(size = 6),
#     panel.spacing.y = unit(0.02, "lines"),
#     panel.spacing.x = unit(0.01, "lines"),
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.margin = margin(1, 1, 1, 1)
#   )
# 
# print(p_main)
# 
# 
# # --- 1) Put the mcrA colorbar at the far top-right
# p_main <- p_main +
#   theme(
#     legend.position      = c(0.985, 0.985),  # top-right corner
#     legend.justification = c(1, 1),
#     legend.background    = element_rect(fill = "white", colour = NA),
#     legend.box.margin    = margin(0, 0, 0, 0),
#     legend.margin        = margin(0, 0, 0, 0)
#   )
# 
# # --- 2) Build a vertical (true-size) DBH legend
# dbh_breaks <- c(20, 30, 40, 50, 60)   # adjust as you like
# gap <- 3                               # spacing between circles (cm)
# 
# legend_df <- tibble::tibble(dbh = rev(dbh_breaks), R = dbh/2) |>
#   dplyr::mutate(
#     y0 = cumsum(dplyr::lag(2*R + gap, default = 0)) + R,
#     x0 = 0
#   )
# 
# pad <- 6
# legend_xlim <- c(-max(legend_df$R) - pad, max(legend_df$R) + pad)
# legend_ylim <- c(-2, max(legend_df$y0 + legend_df$R) + 6)
# 
# p_dbh_legend <- ggplot(legend_df) +
#   ggforce::geom_circle(aes(x0 = x0, y0 = y0, r = R),
#                        fill = NA, colour = "black", linewidth = 0.3) +
#   geom_text(aes(x = max(R) + 5, y = y0, label = dbh), hjust = 0, size = 2.6) +
#   annotate("text",
#            x = legend_xlim[1], y = max(legend_df$y0 + legend_df$R) + 4,
#            label = "DBH (cm)", hjust = 0, vjust = 0.5, size = 2.8) +
#   coord_equal(xlim = legend_xlim, ylim = legend_ylim, expand = FALSE) +
#   theme_void() +
#   theme(
#     plot.background  = element_rect(fill = NA, colour = NA),
#     panel.background = element_rect(fill = NA, colour = NA),
#     plot.margin      = margin(0, 0, 0, 0)
#   )
# 
# # --- 3) Inset the DBH legend to the LEFT of the colorbar (side-by-side)
# # Tweak the box if your figure size/layout changes.
# p_final <- p_main +
#   patchwork::inset_element(
#     p_dbh_legend,
#     left = 0.86,   # left edge of the DBH column
#     right = 0.955, # just left of the colorbar
#     bottom = 0.62, # vertical bounds of the open area
#     top = 0.985,
#     align_to = "panel",
#     clip = FALSE
#   )
# 
# print(p_final)
