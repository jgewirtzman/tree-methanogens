invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
rm(list = ls())

# Complete mcrA Analysis Script with Fixed Ordering
# Load all required libraries
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ape)
library(phytools)
library(viridis)

merged_final<-read.csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/merged_tree_dataset_final.csv")

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
                                           tree_file = "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/metadata/PhytoPhylo",
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
    # Add 1 to handle negatives/zeros for log transformation
    mutate(mcra_copies_log = gene_copies + 1) %>%
    select(species_id, species, material, mcra_copies = gene_copies, mcra_copies_log)
  
  # Calculate summary statistics
  summary_stats <- mcra_data %>%
    group_by(species, material) %>%
    summarise(
      n = n(),
      mean_mcra = mean(mcra_copies_log, na.rm = TRUE),
      se_mcra = sd(mcra_copies_log, na.rm = TRUE) / sqrt(n()),
      median_mcra = median(mcra_copies_log, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Handle cases where SE is NA (only 1 observation)
    mutate(se_mcra = ifelse(is.na(se_mcra), 0, se_mcra))
  
  # Calculate overall mean for each material type (for horizontal line)
  overall_means <- mcra_data %>%
    group_by(material) %>%
    summarise(overall_mean = mean(mcra_copies_log, na.rm = TRUE), .groups = "drop")
  
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
  cat("Original mcrA range (before +1):", 
      min(mcra_data$mcra_copies, na.rm = TRUE), "to", 
      max(mcra_data$mcra_copies, na.rm = TRUE), "\n")
  cat("Transformed range (after +1):", 
      min(mcra_data$mcra_copies_log, na.rm = TRUE), "to", 
      max(mcra_data$mcra_copies_log, na.rm = TRUE), "\n")
  cat("Mean values range:", 
      min(summary_stats$mean_mcra, na.rm = TRUE), "to", 
      max(summary_stats$mean_mcra, na.rm = TRUE), "\n")
  
  n_zeros <- sum(mcra_data$mcra_copies == 0, na.rm = TRUE)
  cat("Number of zero values:", n_zeros, "out of", nrow(mcra_data), "total observations\n")
  
  if(n_zeros / nrow(mcra_data) > 0.8) {
    cat("WARNING: Over 80% of values are zero. Consider using a different transformation.\n")
  }
  
  # Create the plot with conditional legend type
  if (use_continuous_legend && !is.null(phylo_distance_data) && "phylo_distance" %in% names(phylo_distance_data)) {
    # Get min and max species for legend labels (convert to abbreviated names for display)
    min_distance_species_full <- phylo_distance_data$species[which.min(phylo_distance_data$phylo_distance)]
    max_distance_species_full <- phylo_distance_data$species[which.max(phylo_distance_data$phylo_distance)]
    
    # Convert to abbreviated names for legend (with special cases for legend-only)
    min_distance_species <- case_when(
      min_distance_species_full == "Acer rubrum" ~ "A. rubrum",
      min_distance_species_full == "Acer saccharum" ~ "A. saccharum", 
      min_distance_species_full == "Betula alleghaniensis" ~ "B. alleghaniensis",
      min_distance_species_full == "Betula lenta" ~ "B. lenta",
      min_distance_species_full == "Betula papyrifera" ~ "B. papyrifera",
      min_distance_species_full == "Fagus grandifolia" ~ "F. grandifolia",
      min_distance_species_full == "Fraxinus americana" ~ "Fr. americana",
      min_distance_species_full == "Pinus strobus" ~ "P. strobus",
      min_distance_species_full == "Quercus rubra" ~ "Q. rubra",
      min_distance_species_full == "Tsuga canadensis" ~ "Tsuga",
      min_distance_species_full == "Carya ovata" ~ "C. ovata",
      min_distance_species_full == "Kalmia latifolia" ~ "K. latifolia",
      min_distance_species_full == "Prunus serotina" ~ "Pr. serotina",
      min_distance_species_full == "Quercus alba" ~ "Q. alba",
      min_distance_species_full == "Quercus velutina" ~ "Quercus",
      min_distance_species_full == "Sassafras albidum" ~ "S. albidum",
      TRUE ~ min_distance_species_full
    )
    
    max_distance_species <- case_when(
      max_distance_species_full == "Acer rubrum" ~ "A. rubrum",
      max_distance_species_full == "Acer saccharum" ~ "A. saccharum", 
      max_distance_species_full == "Betula alleghaniensis" ~ "B. alleghaniensis",
      max_distance_species_full == "Betula lenta" ~ "B. lenta",
      max_distance_species_full == "Betula papyrifera" ~ "B. papyrifera",
      max_distance_species_full == "Fagus grandifolia" ~ "F. grandifolia",
      max_distance_species_full == "Fraxinus americana" ~ "Fr. americana",
      max_distance_species_full == "Pinus strobus" ~ "P. strobus",
      max_distance_species_full == "Quercus rubra" ~ "Q. rubra",
      max_distance_species_full == "Tsuga canadensis" ~ "Tsuga",
      max_distance_species_full == "Carya ovata" ~ "C. ovata",
      max_distance_species_full == "Kalmia latifolia" ~ "K. latifolia",
      max_distance_species_full == "Prunus serotina" ~ "Pr. serotina",
      max_distance_species_full == "Quercus alba" ~ "Q. alba",
      max_distance_species_full == "Quercus velutina" ~ "Quercus",
      max_distance_species_full == "Sassafras albidum" ~ "S. albidum",
      TRUE ~ max_distance_species_full
    )
    
    p <- ggplot(summary_stats, aes(x = species_position, y = mean_mcra, fill = phylo_distance)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      geom_errorbar(
        aes(ymin = pmax(mean_mcra - se_mcra, 0.1), ymax = mean_mcra + se_mcra),
        width = error_bar_width, 
        size = 0.5,
        color = "black"
      ) +
      geom_hline(
        aes(yintercept = overall_mean, color = I(line_color)), 
        linetype = "solid", 
        size = 1
      ) +
      # geom_text(
      #   aes(label = paste0("n=", n)), 
      #   vjust = -0.3,
      #   hjust = 0.5,
      #   size = 2.5,
      #   y = summary_stats$mean_mcra + summary_stats$se_mcra + 0.1
      # ) +
      facet_wrap(~ material, scales = "free_x", ncol = 2) +
      scale_fill_viridis_c(
        name = "Phylogenetic Distance",
        labels = c(min_distance_species, "", "", "", max_distance_species),
        breaks = seq(min(summary_stats$phylo_distance, na.rm = TRUE), 
                     max(summary_stats$phylo_distance, na.rm = TRUE), 
                     length.out = 5),
        trans = "log10"
      )
  } else {
    # Fallback to discrete colors if no phylogenetic data
    p <- ggplot(summary_stats, aes(x = species_position, y = mean_mcra, fill = species)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      geom_errorbar(
        aes(ymin = pmax(mean_mcra - se_mcra, 0.1), ymax = mean_mcra + se_mcra),
        width = error_bar_width, 
        size = 0.5,
        color = "black"
      ) +
      geom_hline(
        aes(yintercept = overall_mean, color = I(line_color)), 
        linetype = "dashed", 
        size = 1
      ) +
      # geom_text(
      #   aes(label = paste0("n=", n)), 
      #   vjust = -0.3,
      #   hjust = 0.5,
      #   size = 2.5,
      #   y = summary_stats$mean_mcra + summary_stats$se_mcra + 0.1
      # ) +
      facet_wrap(~ material, scales = "free_x", ncol = 2) +
      scale_fill_manual(values = plot_colors, breaks = legend_order, name = "Species")
  }
  
  # Add common plot elements
  p <- p +
    scale_y_log10(
      labels = trans_format("log10", math_format(10^.x)),
      expand = expansion(mult = c(0.02, 0.15))
    ) +
    # Use abbreviated species names with sample sizes as x-axis labels
    scale_x_discrete(labels = function(x) {
      # Extract species name from species_position
      species_names <- sapply(strsplit(as.character(x), "_", fixed = TRUE), function(parts) {
        paste(parts[3:length(parts)], collapse = "_")
      })
      
      # Add sample sizes to abbreviated species names
      sapply(species_names, function(sp_name) {
        # Find the sample size for this species in summary_stats
        matching_row <- summary_stats[summary_stats$species == sp_name, ]
        if(nrow(matching_row) > 0) {
          paste0(matching_row$species_abbrev[1], " (n=", matching_row$n[1], ")")
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
      axis.title.y = element_text(margin = margin(r = 20)),  # Add more space to the right of y-axis title
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(0.8, "lines"),  # Make legend smaller
      legend.text = element_text(size = 10),  # Smaller legend text
      legend.title = element_text(size = 10)  # Smaller legend title
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

ggsave("mcra_barplot.png", result$plot, width = 12, height = 8, dpi = 300)

