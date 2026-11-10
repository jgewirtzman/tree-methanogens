# ==============================================================================
# figS22_radial-sections.R -- Figure S23: radial mcrA distribution in individual trees
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Each circle is one tree, drawn to its DBH: log10(mcrA + 1) in heartwood within 5 cm
# of the centre, sapwood within 5 cm of the bark, interpolated linearly on the log
# scale between. Rows: species with at least three trees, ordered by mean area-weighted
# mcrA (area_weighted_gene(), code/lib/species_gene_data.R); columns: trees.
# Inputs: data/tables/ via code/lib/species_gene_data.R: ymf2021 (2021 survey per tree;
#   species and DBH from trees.csv) and long_all (wood genes, mcrA probe assay, loose
#   threshold; ddpcr.csv + samples.csv)
# Output: outputs/figures/original/supplementary/figS13_tree_radial_sections.png;
#   per-tree summaries to the console
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
})
source("code/lib/fig_style.R")   # SI house style: theme_si(), si_save()

# ---- Load data (canonical tables) ----
# 2021 survey per tree (species and DBH from trees.csv), its wood genes long (loose
# threshold; mcrA probe assay only) and the area weighting: code/lib/species_gene_data.R
source("code/lib/species_gene_data.R")
merged_final <- ymf2021

# ============================================================
# 1) mcrA (probe assay), heartwood and sapwood
# ============================================================
long_data <- long_all %>% filter(gene == "mcrA")

# ============================================================
# 2) Area weighting: area_weighted_gene() (closed-form, no density)
# ============================================================
area_weighted_mcra_vec <- Vectorize(area_weighted_gene)

# ============================================================
# 3) Calculate per-tree metrics
# ============================================================
mcra_with_species <- long_data %>%
  left_join(
    merged_final %>% dplyr::select(tree_id, species_id),
    by = "tree_id"
  ) %>%
  mutate(species = species_name(species_id)) %>%
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

# Species with at least MIN_TREES_PER_SPECIES trees
MIN_TREES_PER_SPECIES <- 3

species_to_include <- tree_weighted %>%
  group_by(species_label) %>%
  summarise(n_trees = n()) %>%
  filter(n_trees >= MIN_TREES_PER_SPECIES) %>%
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

# Cross-sections: rows = species, columns = trees ordered by area-weighted mcrA
# (SI house style: code/lib/fig_style.R)
p_main <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression(log[10] ~ italic(mcrA) ~ "(copies" ~ g^{-1} * ")"),
    option = "plasma",
    guide = guide_colourbar(theme = theme(legend.key.width = unit(10, "lines"),
                                          legend.key.height = unit(0.7, "lines")))
  ) +
  geom_circle(
    data = all_cloud %>% distinct(species_label, tree_rank, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.2
  ) +
  facet_grid(species_label ~ tree_rank, switch = "both") +
  coord_equal() +
  theme_si() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "italic",
                                     size = SI_BASE * 0.8, margin = margin(0, 3, 0, 0)),
    strip.text.x.bottom = element_text(size = SI_BASE * 0.65, colour = "grey30",
                                       margin = margin(1, 0, 0, 0)),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    panel.spacing.y = unit(1, "pt"),
    panel.spacing.x = unit(1, "pt"),
    plot.margin = margin(3, 3, 3, 3)
  )

si_save("outputs/figures/original/supplementary/figS13_tree_radial_sections.png",
        p_main, height = 4.9)
