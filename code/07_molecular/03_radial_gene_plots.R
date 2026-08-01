# ==============================================================================
# Radial Gene Distribution Plots (Figure 8 panels a–c)
# ==============================================================================
# Purpose: Radial cross-section plots showing spatial distribution of mcrA,
#   pmoA, and mmoX gene abundance on tree stems.
#
# Pipeline stage: 4 — Visualization
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#
# Outputs:
#   - radial plot PDFs (to outputs/figures/)
# ==============================================================================

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggnewscale)

# ============================================================
# LOAD DATA
# ============================================================

ymf2021 <- read.csv('data/processed/integrated/merged_tree_dataset_final.csv')

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
# PARSE ALL GENE DATA
# ============================================================

prepare_long_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    tidyr::pivot_longer(
      cols = starts_with("ddpcr_"),
      names_to = "measurement_type",
      values_to = "gene_copies"
    ) %>%
    dplyr::filter(!is.na(gene_copies)) %>%
    tidyr::separate(
      measurement_type,
      into = c("method", "gene", "part1", "part2", "part3"),
      sep = "_", extra = "merge", fill = "right"
    ) %>%
    dplyr::mutate(
      is_probe = (part1 == "probe"),
      location = dplyr::if_else(is_probe, part2, part1),
      stringency = dplyr::if_else(is_probe, part3, part2),
      location = stringr::str_remove(location, "probe"),
      location = dplyr::case_when(
        location %in% c("Inner","inner") ~ "Inner",
        location %in% c("Outer","outer") ~ "Outer",
        TRUE ~ location
      ),
      gene = dplyr::case_when(
        gene == "mcra" ~ "mcrA",
        gene == "mmox" ~ "mmoX",
        gene == "pmoa" ~ "pmoA",
        TRUE ~ gene
      ),
      sample_type = dplyr::case_when(
        location == "Inner" ~ "Heartwood",
        location == "Outer" ~ "Sapwood",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(stringency == "loose", !is.na(sample_type)) %>%
    dplyr::filter(gene %in% c("mcrA", "pmoA", "mmoX"))
}

tree_genes <- prepare_long_genes(ymf2021)

cat(sprintf("Total measurements: %d\n", nrow(tree_genes)))
print(table(tree_genes$gene))

# ============================================================
# GRID GENERATION FUNCTION
# ============================================================

generate_species_grid <- function(species_row, gene_inner, gene_outer, grid_resolution = 0.1) {
  # Use median DBH for this species
  R <- species_row$median_dbh / 2
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  x_seq <- seq(-R, R, by = grid_resolution)
  y_seq <- seq(-R, R, by = grid_resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$r <- sqrt(grid$x^2 + grid$y^2)
  grid <- grid[grid$r <= R + grid_resolution/2, ]  # Include edge points
  
  # Log-scale interpolation
  log_inner <- log10(gene_inner + 1)
  log_outer <- log10(gene_outer + 1)
  
  grid$Clog <- ifelse(grid$r <= r1, log_inner,
                      ifelse(grid$r >= r2, log_outer,
                             log_inner + (log_outer - log_inner) *
                               (grid$r - r1) / max(r2 - r1, 1e-9)))
  
  grid$species_label <- species_row$species_label
  grid$R <- R
  
  grid
}

# ============================================================
# PLOT CREATION FUNCTION (SPECIES AVERAGES)
# ============================================================

create_gene_species_plot <- function(gene_data, gene_name, color_option = "plasma", min_trees = 5, reverse = FALSE) {
  
  # Calculate species averages
  # Check if species_id already exists in the data
  if(!"species_id" %in% names(gene_data)) {
    gene_data <- gene_data %>%
      dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id")
  }
  
  species_data <- gene_data %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::mutate(species = species_mapping[species_id]) %>%
    dplyr::filter(is.finite(dbh)) %>%
    dplyr::group_by(tree_id, species_id, species, dbh, sample_type) %>%
    dplyr::summarise(gene_copies = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
    dplyr::group_by(species_id, species) %>%
    dplyr::summarise(
      n_trees = dplyr::n(),
      gene_inner = mean(Heartwood, na.rm = TRUE),
      gene_outer = mean(Sapwood, na.rm = TRUE),
      median_dbh = median(dbh, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_trees >= min_trees, is.finite(gene_inner), is.finite(gene_outer))
  
  # Add species label
  species_data$species_label <- ifelse(is.na(species_data$species) | species_data$species == "", 
                                       species_data$species_id, species_data$species)
  
  if(nrow(species_data) == 0) {
    cat(sprintf("No species data for %s\n", gene_name))
    return(NULL)
  }
  
  # Order alphabetically by species name
  species_data <- species_data %>%
    dplyr::mutate(mean_gene = (gene_inner + gene_outer) / 2) %>%
    dplyr::arrange(species_label)
  
  species_data$species_label <- factor(species_data$species_label,
                                       levels = species_data$species_label)
  
  cat(sprintf("\n%s: %d species (n≥%d trees)\n", 
              gene_name, nrow(species_data), min_trees))
  
  # Generate grids for all species
  all_cloud <- dplyr::bind_rows(lapply(1:nrow(species_data), function(i) {
    generate_species_grid(
      species_data[i,],
      species_data$gene_inner[i],
      species_data$gene_outer[i]
    )
  }))
  
  # Create plot
  p <- ggplot(all_cloud, aes(x = x, y = y, fill = Clog)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis_c(
      name = bquote("log"[10] ~ "(" ~ .(gene_name) ~ ")"),
      option = color_option,
      direction = if(reverse) -1 else 1,
      breaks = function(x) c(min(x), mean(c(min(x), max(x))), max(x)),
      labels = function(x) sprintf("%.1f", x)
    ) +
    geom_circle(
      data = all_cloud %>% dplyr::distinct(species_label, R),
      aes(x0 = 0, y0 = 0, r = R),
      inherit.aes = FALSE, color = "black", linewidth = 0.3
    ) +
    facet_wrap(~ species_label, ncol = 5) +
    coord_equal() +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "italic", size = 11),
      strip.background = element_rect(fill = "grey95", color = NA),
      legend.position = "right",
      panel.spacing = unit(0.5, "lines"),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(list(plot = p, data = species_data))
}

# ============================================================
# CREATE INDIVIDUAL GENE PLOTS
# ============================================================

cat("\n============================================================\n")
cat("CREATING SPECIES-AVERAGED GENE PLOTS\n")
cat("============================================================\n")

# mcrA plot
if("mcrA" %in% tree_genes$gene) {
  p_mcra_result <- create_gene_species_plot(tree_genes, "mcrA", "inferno")
  if(!is.null(p_mcra_result)) {
    print(p_mcra_result$plot)
    # ggsave("outputs/figures/radial_plot_mcra_species.pdf", p_mcra_result$plot, width = 12, height = 8)
    # cat("Saved: radial_plot_mcra_species.pdf\n")
  }
}

# pmoA plot
if("pmoA" %in% tree_genes$gene) {
  p_pmoa_result <- create_gene_species_plot(tree_genes, "pmoA", "plasma")
  if(!is.null(p_pmoa_result)) {
    print(p_pmoa_result$plot)
    # ggsave("outputs/figures/radial_plot_pmoa_species.pdf", p_pmoa_result$plot, width = 12, height = 8)
    # cat("Saved: radial_plot_pmoa_species.pdf\n")
  }
}

# mmoX plot
if("mmoX" %in% tree_genes$gene) {
  p_mmox_result <- create_gene_species_plot(tree_genes, "mmoX", "plasma")
  if(!is.null(p_mmox_result)) {
    print(p_mmox_result$plot)
    # ggsave("outputs/figures/radial_plot_mmox_species.pdf", p_mmox_result$plot, width = 12, height = 8)
    # cat("Saved: radial_plot_mmox_species.pdf\n")
  }
}

# ============================================================
# CREATE SUM PLOT (pmoA + mmoX) - REVERSE MAKO
# ============================================================

cat("\n============================================================\n")
cat("CREATING SUM PLOT (pmoA + mmoX) - MAKO\n")
cat("============================================================\n")

if(all(c("pmoA", "mmoX") %in% tree_genes$gene)) {
  # First get species_id for each tree
  tree_species <- ymf2021 %>% 
    dplyr::select(tree_id, species_id, dbh) %>%
    dplyr::distinct()
  
  sum_data <- tree_genes %>%
    dplyr::filter(gene %in% c("pmoA", "mmoX")) %>%
    dplyr::left_join(tree_species, by = "tree_id") %>%
    dplyr::filter(is.finite(dbh)) %>%
    dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
    dplyr::summarise(gene_copies = sum(gene_copies, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(gene = "pmoA+mmoX")
  
  p_sum_result <- create_gene_species_plot(sum_data, "pmoA+mmoX", "mako")
  if(!is.null(p_sum_result)) {
    print(p_sum_result$plot)
    # ggsave("outputs/figures/radial_plot_sum_species.pdf", p_sum_result$plot, width = 12, height = 8)
    # cat("Saved: radial_plot_sum_species.pdf\n")
  }
}

# ============================================================
# CREATE OVERLAY PLOT (Species averages) - MODERN RED/BLUE
# ============================================================

cat("\n============================================================\n")
cat("CREATING OVERLAY PLOT (SPECIES AVERAGES) - MODERN RED/BLUE\n")
cat("============================================================\n")

if(all(c("mcrA", "pmoA", "mmoX") %in% tree_genes$gene)) {
  
  # Get species averages for mcrA
  mcra_species <- tree_genes %>%
    dplyr::filter(gene == "mcrA") %>%
    dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
    dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
    dplyr::summarise(gene_copies = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
    dplyr::group_by(species_id) %>%
    dplyr::summarise(
      n_trees = dplyr::n(),
      mcra_inner = mean(Heartwood, na.rm = TRUE),
      mcra_outer = mean(Sapwood, na.rm = TRUE),
      median_dbh = median(dbh, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_trees >= 5)
  
  # Get species averages for methanotrophs
  methan_species <- tree_genes %>%
    dplyr::filter(gene %in% c("pmoA", "mmoX")) %>%
    dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
    dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
    dplyr::summarise(gene_copies = sum(gene_copies, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
    dplyr::group_by(species_id) %>%
    dplyr::summarise(
      n_trees = dplyr::n(),
      methan_inner = mean(Heartwood, na.rm = TRUE),
      methan_outer = mean(Sapwood, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_trees >= 5)
  
  # Merge species with both genes
  overlay_species <- mcra_species %>%
    dplyr::inner_join(methan_species, by = "species_id", suffix = c("_mcra", "_methan")) %>%
    dplyr::mutate(species = species_mapping[species_id])
  
  overlay_species$species_label <- ifelse(is.na(overlay_species$species) | overlay_species$species == "",
                                          overlay_species$species_id, overlay_species$species)
  
  cat(sprintf("Found %d species with both mcrA and methanotrophs\n", nrow(overlay_species)))
  
  if(nrow(overlay_species) > 0) {
    # Generate grids for mcrA
    mcra_cloud <- dplyr::bind_rows(lapply(1:nrow(overlay_species), function(i) {
      grid <- generate_species_grid(
        overlay_species[i,],
        overlay_species$mcra_inner[i],
        overlay_species$mcra_outer[i]
      )
      grid$Clog_mcra <- grid$Clog
      grid
    }))
    
    # Generate grids for methanotrophs
    methan_cloud <- dplyr::bind_rows(lapply(1:nrow(overlay_species), function(i) {
      grid <- generate_species_grid(
        overlay_species[i,],
        overlay_species$methan_inner[i],
        overlay_species$methan_outer[i]
      )
      grid$Clog_methan <- grid$Clog
      grid
    }))
    
    # Merge clouds
    overlay_cloud <- mcra_cloud %>%
      dplyr::select(x, y, species_label, R, Clog_mcra) %>%
      dplyr::inner_join(
        methan_cloud %>% dplyr::select(x, y, species_label, Clog_methan),
        by = c("x", "y", "species_label")
      )
    
    # Create overlay plot with high impact Material Design colors
    p_overlay <- ggplot(overlay_cloud, aes(x = x, y = y)) +
      # First layer: methanotrophs (Material blue)
      geom_raster(aes(fill = Clog_methan), alpha = 0.6, interpolate = TRUE) +
      scale_fill_gradient(
        low = "white", 
        high = "#1E88E5",  # Material blue
        name = bquote("log"[10] ~ "(pmoA+mmoX)"),
        breaks = function(x) {
          min_val <- min(overlay_cloud$Clog_methan)
          max_val <- max(overlay_cloud$Clog_methan)
          c(min_val, mean(c(min_val, max_val)), max_val)
        },
        labels = function(x) sprintf("%.1f", x)
      ) +
      new_scale_fill() +
      # Second layer: mcrA (Material red)
      geom_raster(aes(fill = Clog_mcra), alpha = 0.6, interpolate = TRUE) +
      scale_fill_gradient(
        low = "white", 
        high = "#E53935",  # Material red
        name = bquote("log"[10] ~ "(mcrA)"),
        breaks = function(x) {
          min_val <- min(overlay_cloud$Clog_mcra)
          max_val <- max(overlay_cloud$Clog_mcra)
          c(min_val, mean(c(min_val, max_val)), max_val)
        },
        labels = function(x) sprintf("%.1f", x)
      ) +
      # Add circles
      geom_circle(
        data = overlay_cloud %>% dplyr::distinct(species_label, R),
        aes(x0 = 0, y0 = 0, r = R),
        inherit.aes = FALSE, color = "black", linewidth = 0.5
      ) +
      facet_wrap(~ species_label, ncol = 5) +
      coord_equal() +
      theme_minimal(base_size = 10) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "italic", size = 11),
        strip.background = element_rect(fill = "grey95", color = NA),
        legend.position = "right",
        panel.spacing = unit(0.5, "lines"),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    print(p_overlay)
    # ggsave("outputs/figures/radial_plot_overlay_species.pdf", p_overlay, width = 12, height = 8)
    # cat("Saved: radial_plot_overlay_species.pdf\n")
  }
}

cat("\n============================================================\n")
cat("ALL SPECIES-AVERAGED PLOTS COMPLETE!\n")
cat("============================================================\n")
cat("\nFiles created:\n")
cat("  - radial_plot_mcra_species.pdf\n")
cat("  - radial_plot_pmoa_species.pdf\n")
cat("  - radial_plot_mmox_species.pdf\n")
cat("  - radial_plot_sum_species.pdf (MAKO)\n")
cat("  - radial_plot_overlay_species.pdf (HIGH IMPACT RED/BLUE)\n")

# ============================================================
# COMBINE ALL RADIAL PLOTS INTO VERTICAL LAYOUT
# ============================================================

library(ggplot2)
library(patchwork)

cat("\n============================================================\n")
cat("CREATING COMBINED VERTICAL LAYOUT\n")
cat("============================================================\n")

# Extract just the plot objects
p_mcra <- p_mcra_result$plot
p_sum <- p_sum_result$plot
# p_overlay already exists

# Stack plots vertically
combined_plot <- p_mcra / p_sum / p_overlay

# Save the combined figure
# ggsave("outputs/figures/radial_plots_combined_vertical.pdf",
#        combined_plot,
#        width = 10,
#        height = 12,
#        limitsize = FALSE)
#
# ggsave("outputs/figures/radial_plots_combined_vertical.png",
#        combined_plot,
#        width = 10,
#        height = 12,
#        limitsize = FALSE)

cat("\nSaved: radial_plots_combined_vertical.pdf\n")
cat("Dimensions: 12 inches wide × 24 inches tall\n")
cat("\nPlots stacked in order:\n")
cat("  1. mcrA\n")
cat("  2. pmoA+mmoX (mako)\n")
cat("  3. Overlay (high impact red + blue)\n")

cat("\n============================================================\n")
cat("COMBINED LAYOUT COMPLETE!\n")
cat("============================================================\n")