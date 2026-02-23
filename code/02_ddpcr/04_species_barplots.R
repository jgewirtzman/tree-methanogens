# ==============================================================================
# Species Gene Abundance Barplots (Figure 4)
# ==============================================================================
# Purpose: Creates barplots of gene abundance by species with phylogenetic
#   coloring.
#
# Pipeline stage: 03 Analysis
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#
# Outputs:
#   - mcra_barplot.png
# ==============================================================================

library(tidyverse)
library(scales)
library(RColorBrewer)
library(ape)
library(phytools)
library(viridis)

# Species name mapping
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

# Data Quality Check Function
check_data_quality <- function(merged_final) {
  # Reshape data (same as in your plotting function)
  long_data <- merged_final %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select_if(~ !all(is.na(.))) %>%
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
      is_probe   = (part1 == "probe"),
      location   = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(
        location %in% c("Inner","inner")     ~ "Inner",
        location %in% c("Outer","outer")     ~ "Outer",
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
        location == "Inner"   ~ "Heartwood",
        location == "Outer"   ~ "Sapwood",
        location == "Mineral" ~ "Mineral",
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      ),
      sample_id = paste(tree_id, sample_type, sep = "_")
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX"))
  
  # Check mcrA probe data
  mcra_probe_data <- long_data %>%
    filter(gene == "mcrA", is_probe)
  
  # Summary statistics for mcrA probe
  cat("=== mcrA PROBE DATA QUALITY CHECK ===\n")
  cat("Total mcrA probe measurements:", nrow(mcra_probe_data), "\n")
  cat("Unique samples with mcrA probe data:", length(unique(mcra_probe_data$sample_id)), "\n")
  
  mcra_summary <- mcra_probe_data %>%
    summarise(
      min_val = min(gene_copies, na.rm = TRUE),
      max_val = max(gene_copies, na.rm = TRUE),
      mean_val = mean(gene_copies, na.rm = TRUE),
      median_val = median(gene_copies, na.rm = TRUE),
      n_zeros = sum(gene_copies == 0, na.rm = TRUE),
      n_negatives = sum(gene_copies < 0, na.rm = TRUE),
      n_na = sum(is.na(gene_copies))
    )
  
  print(mcra_summary)
  
  if(mcra_summary$n_zeros > 0) {
    cat("\nSamples with zero mcrA values:", mcra_summary$n_zeros, "\n")
  }
  
  if(mcra_summary$n_negatives > 0) {
    cat("Samples with negative mcrA values:", mcra_summary$n_negatives, "\n")
  }
  
  return(list(
    mcra_summary = mcra_summary
  ))
}

# Main Bar Plot Function with Fixed Ordering
create_mcra_barplot_by_species <- function(merged_final, species_mapping, 
                                           tree_file = "data/processed/metadata/PhytoPhylo",
                                           error_bar_width = 0.2) {
  
  # Process data
  long_data <- merged_final %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select_if(~ !all(is.na(.))) %>%
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
      is_probe   = (part1 == "probe"),
      location   = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(
        location %in% c("Inner","inner")     ~ "Inner",
        location %in% c("Outer","outer")     ~ "Outer", 
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
        location == "Inner"   ~ "Heartwood",
        location == "Outer"   ~ "Sapwood",
        location == "Mineral" ~ "Mineral", 
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      ),
      sample_id = paste(tree_id, sample_type, sep = "_")
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene == "mcrA", is_probe)  # Only mcrA probe data
  
  # Map species_id to full species name (keep full names for phylogenetic analysis)
  mcra_data <- long_data %>%
    left_join(
      merged_final %>% select(tree_id, species_id), 
      by = "tree_id"
    ) %>%
    mutate(
      # Map species_id to full species name 
      species = case_when(
        species_id %in% names(species_mapping) ~ species_mapping[species_id],
        !is.na(species_id) ~ paste("Unknown species:", species_id),
        TRUE ~ "No species ID found"
      ),
      material = case_when(
        sample_type == "Heartwood" ~ "Heartwood",
        sample_type == "Sapwood" ~ "Sapwood", 
        sample_type == "Mineral" ~ "Mineral Soil",
        sample_type == "Organic" ~ "Organic Soil",
        TRUE ~ sample_type
      ),
      # Set factor levels to control facet order
      material = factor(material, levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
    ) %>%
    filter(!is.na(species)) %>%
    # CHANGED: compute per-row log10(x+1) instead of (x+1)
    mutate(mcra_log10p1 = log10(pmax(gene_copies, 0) + 1)) %>%
    select(species_id, species, material, mcra_copies = gene_copies, mcra_log10p1)
  
  # Calculate summary statistics
  summary_stats <- mcra_data %>%
    group_by(species, material) %>%
    summarise(
      n = n(),
      mean_mcra = mean(mcra_log10p1, na.rm = TRUE),
      se_mcra = sd(mcra_log10p1, na.rm = TRUE) / sqrt(n()),
      median_mcra = median(mcra_log10p1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Handle cases where SE is NA (only 1 observation)
    mutate(se_mcra = ifelse(is.na(se_mcra), 0, se_mcra))
  
  # Calculate overall mean for each material type (for horizontal line)
  overall_means <- mcra_data %>%
    group_by(material) %>%
    summarise(overall_mean = mean(mcra_log10p1, na.rm = TRUE), .groups = "drop")
  
  # Add overall means and line colors
  summary_stats <- summary_stats %>%
    left_join(overall_means, by = "material") %>%
    mutate(
      line_color = case_when(
        material == "Heartwood" ~ "#a6611a",
        material == "Sapwood" ~ "#dfc27d", 
        material == "Mineral Soil" ~ "#80cdc1",
        material == "Organic Soil" ~ "#018571",
        TRUE ~ "red"
      )
    )
  
  # KEY FIX: Create position variable for ordering within each facet
  summary_stats <- summary_stats %>%
    group_by(material) %>%
    arrange(mean_mcra) %>%
    mutate(
      position = row_number(),
      # Create a unique identifier that combines position and material
      species_position = paste0(material, "_", sprintf("%02d", position), "_", species)
    ) %>%
    ungroup() %>%
    # Make species_position a factor with correct ordering
    mutate(species_position = factor(species_position, levels = species_position))
  
  # --- PHYLOGENETIC DISTANCE CALCULATION ---
  phylo_colors <- NULL
  
  if (file.exists(tree_file)) {
    tryCatch({
      tree_scenario1 <- read.tree(tree_file)
      present_species <- unique(summary_stats$species)
      present_species_clean <- gsub(" ", "_", present_species)
      tree_species <- intersect(present_species_clean, tree_scenario1$tip.label)
      
      if (length(tree_species) > 1) {
        pruned_tree <- keep.tip(tree_scenario1, tree_species)
        dist_matrix <- cophenetic(pruned_tree)
        average_distances <- rowMeans(dist_matrix)
        
        distance_df <- data.frame(
          species_clean = names(average_distances),
          distance = as.numeric(average_distances)
        ) %>%
          mutate(species = gsub("_", " ", species_clean)) %>%
          arrange(distance)
        
        phylo_colors <- setNames(
          viridis(nrow(distance_df)), 
          distance_df$species
        )
        
        cat("=== PHYLOGENETIC DISTANCE INFO ===\n")
        cat("Species with phylogenetic data:", length(tree_species), "\n")
        cat("Using phylogenetic distance-based colors\n\n")
        
      } else {
        cat("Warning: Not enough species overlap with phylogenetic tree\n")
      }
    }, error = function(e) {
      cat("Warning: Could not calculate phylogenetic distances:", e$message, "\n")
    })
  } else {
    cat("Warning: Tree file not found at:", tree_file, "\n")
  }
  
  # Set up colors and phylogenetic distance info for continuous scale
  phylo_distance_data <- NULL
  
  if (is.null(phylo_colors)) {
    present_species <- unique(summary_stats$species)
    n_species <- length(present_species)
    
    if(n_species <= 12) {
      plot_colors <- RColorBrewer::brewer.pal(min(n_species, 12), "Set3")[1:n_species]
    } else {
      plot_colors <- rainbow(n_species)
    }
    names(plot_colors) <- present_species
    
    # No phylogenetic data available - use discrete legend
    use_continuous_legend <- FALSE
    legend_order <- sort(present_species)
  } else {
    plot_colors <- phylo_colors
    
    # Extract phylogenetic distance data for continuous scale
    phylo_distance_data <- data.frame(
      species = names(phylo_colors),
      color = phylo_colors,
      stringsAsFactors = FALSE
    )
    
    # Get the original distance values from the tree calculation
    # We need to recalculate to get actual distance values
    if (file.exists(tree_file)) {
      tree_scenario1 <- read.tree(tree_file)
      present_species <- unique(summary_stats$species)
      present_species_clean <- gsub(" ", "_", present_species)
      tree_species <- intersect(present_species_clean, tree_scenario1$tip.label)
      
      if (length(tree_species) > 1) {
        pruned_tree <- keep.tip(tree_scenario1, tree_species)
        dist_matrix <- cophenetic(pruned_tree)
        average_distances <- rowMeans(dist_matrix)
        
        phylo_distance_data <- data.frame(
          species = gsub("_", " ", names(average_distances)),
          phylo_distance = as.numeric(average_distances),
          stringsAsFactors = FALSE
        )
      }
    }
    
    use_continuous_legend <- TRUE
  }
  
  # Add phylogenetic distance to summary_stats if available (using full names)
  if (!is.null(phylo_distance_data) && "phylo_distance" %in% names(phylo_distance_data)) {
    summary_stats <- summary_stats %>%
      left_join(phylo_distance_data, by = "species")
  }
  
  # Add abbreviated species names AFTER all phylogenetic calculations
  summary_stats <- summary_stats %>%
    mutate(
      species_abbrev = case_when(
        species == "Acer rubrum" ~ "A. rubrum",
        species == "Acer saccharum" ~ "A. saccharum", 
        species == "Betula alleghaniensis" ~ "B. alleghaniensis",
        species == "Betula lenta" ~ "B. lenta",
        species == "Betula papyrifera" ~ "B. papyrifera",
        species == "Fagus grandifolia" ~ "F. grandifolia",
        species == "Fraxinus americana" ~ "Fr. americana",
        species == "Pinus strobus" ~ "P. strobus",
        species == "Quercus rubra" ~ "Q. rubra",
        species == "Tsuga canadensis" ~ "T. canadensis",
        species == "Carya ovata" ~ "C. ovata",
        species == "Kalmia latifolia" ~ "K. latifolia",
        species == "Prunus serotina" ~ "Pr. serotina",
        species == "Quercus alba" ~ "Q. alba",
        species == "Quercus velutina" ~ "Q. velutina",
        species == "Sassafras albidum" ~ "S. albidum",
        TRUE ~ species
      )
    )
  
  # Print diagnostics
  cat("=== DIAGNOSTIC INFO ===\n")
  cat("Original mcrA range:", 
      min(mcra_data$mcra_copies, na.rm = TRUE), "to", 
      max(mcra_data$mcra_copies, na.rm = TRUE), "\n")
  cat("log10p1 range (per-row):", 
      min(mcra_data$mcra_log10p1, na.rm = TRUE), "to", 
      max(mcra_data$mcra_log10p1, na.rm = TRUE), "\n")
  cat("Mean(log10p1) range:", 
      min(summary_stats$mean_mcra, na.rm = TRUE), "to", 
      max(summary_stats$mean_mcra, na.rm = TRUE), "\n")
  
  n_zeros <- sum(mcra_data$mcra_copies == 0, na.rm = TRUE)
  cat("Number of zero values:", n_zeros, "out of", nrow(mcra_data), "total observations\n")
  
  if(n_zeros / nrow(mcra_data) > 0.8) {
    cat("WARNING: Over 80% of values are zero. Consider using a different transformation.\n")
  }
  
  # Prepare individual data points for jittering
  individual_points <- mcra_data %>%
    left_join(phylo_distance_data, by = "species") %>%
    left_join(summary_stats %>% select(species, material, species_position), 
              by = c("species", "material"))
  
  # Create annotation data for material types
  annotation_data <- summary_stats %>%
    group_by(material) %>%
    summarise(
      max_y = max(mean_mcra + se_mcra, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Position annotations at top-left corner of each facet
      y_pos = Inf,  # Use Inf to position at top
      x_pos = -Inf,  # Use -Inf to position at left
      # Clean up material names for display
      material_label = case_when(
        material == "Heartwood" ~ "Heartwood",
        material == "Sapwood" ~ "Sapwood",
        material == "Mineral Soil" ~ "Mineral Soil", 
        material == "Organic Soil" ~ "Organic Soil",
        TRUE ~ as.character(material)
      )
    )
  
  # Get min and max species for legend labels (convert to shortened names)
  min_distance_species_full <- phylo_distance_data$species[which.min(phylo_distance_data$phylo_distance)]
  max_distance_species_full <- phylo_distance_data$species[which.max(phylo_distance_data$phylo_distance)]
  
  # Convert to shortened names for legend
  min_distance_species <- case_when(
    min_distance_species_full == "Tsuga canadensis" ~ "Tsuga",
    min_distance_species_full == "Quercus velutina" ~ "Quercus", 
    TRUE ~ format(min(phylo_distance_data$phylo_distance, na.rm = TRUE), digits = 2)
  )
  
  max_distance_species <- case_when(
    max_distance_species_full == "Tsuga canadensis" ~ "Tsuga",
    max_distance_species_full == "Quercus velutina" ~ "Quercus",
    TRUE ~ format(max(phylo_distance_data$phylo_distance, na.rm = TRUE), digits = 2)
  )
  
  p <- ggplot(summary_stats, aes(x = species_position, y = mean_mcra, fill = phylo_distance)) +
    geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_jitter(
      data = individual_points,
      aes(x = species_position, y = mcra_log10p1, fill = phylo_distance),
      width = 0.1, height = 0, alpha = 0.6, size = 1, shape = 21, color = "black"
    ) +
    geom_errorbar(
      aes(ymin = pmax(mean_mcra - se_mcra, 0.1), ymax = mean_mcra + se_mcra),
      width = error_bar_width, 
      linewidth = 0.5,
      color = "black"
    ) +
    geom_hline(
      aes(yintercept = overall_mean), 
      color = "grey50", alpha = 0.6, linewidth = 0.6
    ) +
    geom_point(
      data = summary_stats,
      aes(x = 1, y = overall_mean),
      size = 3, shape = 23, color = "grey50", alpha=0.9, fill = summary_stats$line_color
    ) +
    # Add material type annotations
    geom_text(
      data = annotation_data,
      aes(x = x_pos, y = y_pos, label = material_label),
      inherit.aes = FALSE,  # Don't inherit fill aesthetic
      #fontface = "bold",
      size = 4,
      hjust = 0,  # Left-align text
      vjust = 1   # Top-align text
    ) +
    facet_wrap(~ material, scales = "free_x", ncol = 2) +
    scale_fill_viridis_c(
      name = "Phylogenetic Distance",
      labels = c(min_distance_species, "", "", "", max_distance_species),
      breaks = seq(min(summary_stats$phylo_distance, na.rm = TRUE), 
                   max(summary_stats$phylo_distance, na.rm = TRUE), 
                   length.out = 5),
      trans = "log10"
    ) +
    scale_y_continuous(
      breaks = c(0, 1, 2, 3, 4, 5, 6),
      labels = c(expression(10^0), expression(10^1), expression(10^2), 
                 expression(10^3), expression(10^4), expression(10^5), expression(10^6)),
      expand = expansion(mult = c(0.02, 0.2))  # Increased top expansion to accommodate annotations
    ) +
    scale_x_discrete(labels = function(x) {
      # Extract species name from species_position
      species_names <- sapply(strsplit(as.character(x), "_", fixed = TRUE), function(parts) {
        paste(parts[3:length(parts)], collapse = "_")
      })
      
      # Return just abbreviated species names (no sample sizes)
      sapply(species_names, function(sp_name) {
        matching_row <- summary_stats[summary_stats$species == sp_name, ]
        if(nrow(matching_row) > 0) {
          matching_row$species_abbrev[1]  # Just the abbreviated name
        } else {
          sp_name
        }
      })
    }) +
    labs(
      x = NULL,
      y = expression("mcrA gene abundance (copies g"^-1*")"),
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(margin = margin(r = 2)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),  # Hide facet titles since we're using annotations
      legend.position = "bottom",
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      
      # Minimize spacing between facets
      panel.spacing = unit(0.15, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"),  # Small top margin for annotations
      
      # Reduce legend spacing
      legend.margin = margin(t = 1, r = 0, b = 1, l = -50),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = -50),
      legend.spacing = unit(0.5, "lines"),
      legend.box.spacing = unit(0, "lines")
    )
  
  return(list(
    plot = p,
    data = mcra_data,
    summary = summary_stats,
    phylo_colors = plot_colors
  ))
}

# Run the analysis
cat("Running data quality check...\n")
quality_check <- check_data_quality(merged_final)

cat("Creating mcrA bar plot...\n")
result <- create_mcra_barplot_by_species(merged_final, species_mapping)
print(result$plot)

# ggsave("outputs/figures/mcra_barplot.png", result$plot, width = 12, height = 8, dpi = 300)