# Load required libraries
library(ggplot2)
library(dplyr)
library(ape)
library(phytools)
library(viridis)
library(RColorBrewer)

fg<-read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW.csv')
head(fg)
str(fg)

# Convert coordinates to numeric
fg$Latitude <- as.numeric(fg$Latitude)
fg$Longitude <- as.numeric(fg$Longitude)
fg$PX <- as.numeric(fg$PX)
fg$PY <- as.numeric(fg$PY)

# Remove rows with missing coordinates
fg_clean <- fg %>%
  filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(PX) & !is.na(PY))

# Calculate basal area for each stem (π × (DBH/2)²) in cm², then convert to m²
fg_clean$BasalArea_cm2 <- pi * (fg_clean$DBH / 2)^2
fg_clean$BasalArea_m2 <- fg_clean$BasalArea_cm2 / 10000

# Species name mapping for phylogenetic analysis
species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "ACPE" = "Acer pensylvanicum",
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "SAAL" = "Sassafras albidum",
  "HAVI" = "Hamamelis virginiana",
  "LITU" = "Liriodendron tulipifera",
  "VACO" = "Vaccinium corymbosum"
)

# Map species codes to full names
fg_clean$Species_Name <- case_when(
  fg_clean$Species_Code %in% names(species_mapping) ~ species_mapping[fg_clean$Species_Code],
  fg_clean$Species_Code == "" ~ "Unknown",
  !is.na(fg_clean$Species_Code) ~ paste("Unknown (", fg_clean$Species_Code, ")", sep = ""),
  TRUE ~ "Unknown"
)

cat("=== ROTATED COORDINATE TRANSFORMATION WITH OUTLIER HANDLING ===\n")
cat("Starting with", nrow(fg_clean), "stems with complete coordinates\n\n")

# OUTLIER DETECTION FUNCTION for local coordinates only
detect_outliers <- function(x, multiplier = 0.75) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - multiplier * IQR
  upper_bound <- Q3 + multiplier * IQR
  return(x < lower_bound | x > upper_bound)
}

# DETECT OUTLIERS IN LOCAL COORDINATES
fg_clean$outlier_PX <- detect_outliers(fg_clean$PX)
fg_clean$outlier_PY <- detect_outliers(fg_clean$PY)
fg_clean$local_outlier <- fg_clean$outlier_PX | fg_clean$outlier_PY

# Count outliers
n_outliers <- sum(fg_clean$local_outlier)
n_good <- sum(!fg_clean$local_outlier)

cat("=== LOCAL OUTLIER DETECTION (0.75 × IQR threshold) ===\n")
cat("Good local coordinates:", n_good, "(", round(n_good/nrow(fg_clean)*100, 1), "%)\n")
cat("Local outliers:", n_outliers, "(", round(n_outliers/nrow(fg_clean)*100, 1), "%)\n\n")

if(n_outliers > 0) {
  outlier_details <- fg_clean %>%
    filter(local_outlier) %>%
    dplyr::select(Tag, PX, PY, Latitude, Longitude, outlier_PX, outlier_PY) %>%
    arrange(Tag)
  
  cat("=== LOCAL OUTLIER DETAILS ===\n")
  print(outlier_details)
  cat("\n")
} else {
  cat("No local outliers detected with current threshold\n\n")
}

# Show coordinate ranges for diagnostics
cat("=== COORDINATE RANGES ===\n")
cat("PX range:", min(fg_clean$PX), "to", max(fg_clean$PX), "\n")
cat("PY range:", min(fg_clean$PY), "to", max(fg_clean$PY), "\n")
cat("Latitude range:", min(fg_clean$Latitude), "to", max(fg_clean$Latitude), "\n")
cat("Longitude range:", min(fg_clean$Longitude), "to", max(fg_clean$Longitude), "\n\n")

# AFFINE TRANSFORMATION FUNCTION
affine_transform <- function(px, py, params) {
  # params: c(scale_x, scale_y, rotation_angle, offset_x, offset_y, translate_lon, translate_lat)
  scale_x <- params[1]
  scale_y <- params[2] 
  theta <- params[3]  # rotation angle in radians
  offset_x <- params[4]  # offset in local coordinates
  offset_y <- params[5]
  translate_lon <- params[6]  # translation to GPS
  translate_lat <- params[7]
  
  # Apply offset
  px_adj <- px - offset_x
  py_adj <- py - offset_y
  
  # Apply scaling
  px_scaled <- px_adj * scale_x
  py_scaled <- py_adj * scale_y
  
  # Apply rotation
  lon_local <- px_scaled * cos(theta) - py_scaled * sin(theta)
  lat_local <- px_scaled * sin(theta) + py_scaled * cos(theta)
  
  # Translate to GPS coordinates
  longitude <- lon_local + translate_lon
  latitude <- lat_local + translate_lat
  
  return(list(longitude = longitude, latitude = latitude))
}

# OBJECTIVE FUNCTION for optimization
objective_function <- function(params, px, py, true_lon, true_lat) {
  predicted <- affine_transform(px, py, params)
  
  # Calculate RMSE
  lon_error <- sum((predicted$longitude - true_lon)^2)
  lat_error <- sum((predicted$latitude - true_lat)^2)
  
  return(sqrt((lon_error + lat_error) / length(px)))
}

# INITIAL PARAMETER ESTIMATES using ALL points (including outliers for rotation detection)
px_range <- max(fg_clean$PX) - min(fg_clean$PX)
py_range <- max(fg_clean$PY) - min(fg_clean$PY)
lon_range <- max(fg_clean$Longitude) - min(fg_clean$Longitude)
lat_range <- max(fg_clean$Latitude) - min(fg_clean$Latitude)

# Estimate scales (convert plot units to GPS degrees)
initial_scale_x <- lon_range / px_range
initial_scale_y <- lat_range / py_range

# Initial parameters: scale_x, scale_y, rotation, offset_x, offset_y, translate_lon, translate_lat
initial_params <- c(
  initial_scale_x,  # scale_x
  initial_scale_y,  # scale_y  
  0,                # rotation (start with 0)
  mean(fg_clean$PX), # offset_x (center of PX)
  mean(fg_clean$PY), # offset_y (center of PY)
  mean(fg_clean$Longitude), # translate_lon
  mean(fg_clean$Latitude)   # translate_lat
)

cat("=== INITIAL PARAMETER ESTIMATES ===\n")
cat("Scale X:", initial_params[1], "\n")
cat("Scale Y:", initial_params[2], "\n")
cat("Rotation:", initial_params[3], "radians\n")
cat("Offset X:", initial_params[4], "\n")
cat("Offset Y:", initial_params[5], "\n")
cat("Translate Lon:", initial_params[6], "\n")
cat("Translate Lat:", initial_params[7], "\n\n")

# OPTIMIZE TRANSFORMATION PARAMETERS using ALL points (including outliers)
cat("Optimizing transformation parameters using ALL points (including outliers)...\n")
result <- optim(
  par = initial_params,
  fn = objective_function,
  px = fg_clean$PX,
  py = fg_clean$PY,
  true_lon = fg_clean$Longitude,
  true_lat = fg_clean$Latitude,
  method = "BFGS"
)

optimal_params <- result$par
cat("Optimization converged:", result$convergence == 0, "\n")
cat("Final RMSE on all points:", result$value, "\n\n")

# DISPLAY OPTIMAL PARAMETERS
cat("=== OPTIMAL TRANSFORMATION PARAMETERS ===\n")
cat("Scale X:", optimal_params[1], "\n")
cat("Scale Y:", optimal_params[2], "\n")
cat("Rotation:", optimal_params[3], "radians (", optimal_params[3] * 180 / pi, "degrees)\n")
cat("Offset X:", optimal_params[4], "\n")
cat("Offset Y:", optimal_params[5], "\n")
cat("Translate Longitude:", optimal_params[6], "\n")
cat("Translate Latitude:", optimal_params[7], "\n\n")

# APPLY TRANSFORMATION TO ALL POINTS
transformed <- affine_transform(fg_clean$PX, fg_clean$PY, optimal_params)
fg_clean$Longitude_converted <- transformed$longitude
fg_clean$Latitude_converted <- transformed$latitude

# CREATE FINAL DATASET - REMOVE LOCAL OUTLIERS
fg_final <- fg_clean %>% 
  filter(!local_outlier) %>%
  mutate(
    Longitude_final = Longitude_converted,
    Latitude_final = Latitude_converted
  )

# CALCULATE PHYLOGENETIC COLORS (DISCRETE)
phylo_colors <- NULL
phylo_distance_data <- NULL
tree_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/metadata/PhytoPhylo"

# Get all species present in final dataset
all_species <- unique(fg_final$Species_Name)
known_species <- all_species[!grepl("Unknown", all_species)]
unknown_species <- all_species[grepl("Unknown", all_species)]

if (file.exists(tree_file) && length(known_species) > 1) {
  tryCatch({
    tree_scenario1 <- read.tree(tree_file)
    present_species_clean <- gsub(" ", "_", known_species)
    tree_species <- intersect(present_species_clean, tree_scenario1$tip.label)
    
    if (length(tree_species) > 1) {
      pruned_tree <- keep.tip(tree_scenario1, tree_species)
      dist_matrix <- cophenetic(pruned_tree)
      average_distances <- rowMeans(dist_matrix)
      
      # Create discrete color assignments based on phylogenetic distance
      distance_df <- data.frame(
        species_clean = names(average_distances),
        distance = as.numeric(average_distances)
      ) %>%
        mutate(species = gsub("_", " ", species_clean)) %>%
        arrange(distance) %>%
        mutate(
          # Create discrete groups based on distance quantiles
          distance_group = cut(distance, 
                               breaks = quantile(distance, probs = seq(0, 1, length.out = min(8, length(distance)))), 
                               include.lowest = TRUE, 
                               labels = FALSE),
          # Assign viridis colors to each group
          color = viridis(max(distance_group, na.rm = TRUE))[distance_group]
        )
      
      # Create color mapping for known species
      phylo_colors <- setNames(distance_df$color, distance_df$species)
      
      # Store phylogenetic distance data for legend ordering
      phylo_distance_data <- distance_df %>%
        dplyr::select(species, distance, distance_group) %>%
        arrange(distance)
      
      cat("=== PHYLOGENETIC COLOR MAPPING ===\n")
      cat("Species with phylogenetic data:", length(tree_species), "\n")
      cat("Using discrete phylogenetic distance-based colors\n")
      print(phylo_distance_data)
      cat("\n")
      
    } else {
      cat("Warning: Not enough species overlap with phylogenetic tree\n")
    }
  }, error = function(e) {
    cat("Warning: Could not calculate phylogenetic distances:", e$message, "\n")
  })
} else {
  cat("Warning: Tree file not found or insufficient known species, using discrete colors\n")
}

# Create complete color mapping including unknowns
if (is.null(phylo_colors)) {
  # Fallback: use discrete colors for known species
  if(length(known_species) <= 12) {
    known_colors <- RColorBrewer::brewer.pal(min(max(3, length(known_species)), 12), "Set3")[1:length(known_species)]
  } else {
    known_colors <- rainbow(length(known_species))
  }
  names(known_colors) <- known_species
  phylo_colors <- known_colors
}

# Add grey colors for unknown species
unknown_colors <- rep("grey80", length(unknown_species))
names(unknown_colors) <- unknown_species

# Combine all colors
final_colors <- c(phylo_colors, unknown_colors)

# Create legend order: phylogenetic species first (ordered by distance), then unknowns
if (!is.null(phylo_distance_data)) {
  legend_order <- c(phylo_distance_data$species, unknown_species)
} else {
  legend_order <- c(sort(known_species), unknown_species)
}

cat("=== FINAL DATASET ===\n")
cat("Original points:", nrow(fg_clean), "\n")
cat("After removing local outliers:", nrow(fg_final), "\n")
cat("Removed:", nrow(fg_clean) - nrow(fg_final), "local outliers\n\n")

# VALIDATION STATISTICS on good points only (non-outliers)
good_points <- fg_clean %>% filter(!local_outlier)
rmse_lon <- sqrt(mean((good_points$Longitude - good_points$Longitude_converted)^2))
rmse_lat <- sqrt(mean((good_points$Latitude - good_points$Latitude_converted)^2))
mae_lon <- mean(abs(good_points$Longitude - good_points$Longitude_converted))
mae_lat <- mean(abs(good_points$Latitude - good_points$Latitude_converted))

cat("=== VALIDATION RESULTS (Non-outlier points) ===\n")
cat("Longitude RMSE:", rmse_lon, "\n")
cat("Latitude RMSE:", rmse_lat, "\n")
cat("Longitude MAE:", mae_lon, "\n")
cat("Latitude MAE:", mae_lat, "\n")
cat("Approximate accuracy: ±", round(max(rmse_lon, rmse_lat) * 111320, 2), "meters\n\n")

# CONVERSION FUNCTION FOR NEW DATA
cat("=== CONVERSION FUNCTION FOR NEW DATA ===\n")
cat("# Use this function to convert new PX/PY coordinates:\n")
cat("convert_coordinates <- function(px, py) {\n")
cat("  params <- c(", paste(sprintf("%.10f", optimal_params), collapse = ", "), ")\n")
cat("  \n")
cat("  scale_x <- params[1]\n")
cat("  scale_y <- params[2]\n") 
cat("  theta <- params[3]\n")
cat("  offset_x <- params[4]\n")
cat("  offset_y <- params[5]\n")
cat("  translate_lon <- params[6]\n")
cat("  translate_lat <- params[7]\n")
cat("  \n")
cat("  px_adj <- px - offset_x\n")
cat("  py_adj <- py - offset_y\n")
cat("  px_scaled <- px_adj * scale_x\n")
cat("  py_scaled <- py_adj * scale_y\n")
cat("  \n")
cat("  longitude <- px_scaled * cos(theta) - py_scaled * sin(theta) + translate_lon\n")
cat("  latitude <- px_scaled * sin(theta) + py_scaled * cos(theta) + translate_lat\n")
cat("  \n")
cat("  return(data.frame(longitude = longitude, latitude = latitude))\n")
cat("}\n\n")

# VISUALIZATIONS
# 1. Local outlier identification
outlier_plot <- ggplot(fg_clean, aes(x = PX, y = PY, color = local_outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Keep in final dataset", "TRUE" = "Remove from final dataset"),
    name = "Status"
  ) +
  labs(title = "Local Coordinate Outlier Detection (0.75 × IQR)",
       subtitle = "Red points used for rotation calculation but excluded from final dataset",
       x = "PX (local)", y = "PY (local)") +
  theme_minimal() +
  coord_fixed()

print(outlier_plot)

# 2. Transformation validation using all points
validation_plot <- ggplot(fg_clean) +
  geom_point(aes(x = Longitude, y = Latitude, color = "Original GPS"), 
             alpha = 0.6, size = 1.5) +
  geom_point(aes(x = Longitude_converted, y = Latitude_converted, color = "Converted from PX/PY"), 
             alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Original GPS" = "blue", "Converted from PX/PY" = "red")) +
  labs(title = "Transformation Results (All Points Used for Rotation)",
       subtitle = paste("Rotation angle:", round(optimal_params[3] * 180 / pi, 1), "degrees"),
       x = "Longitude", y = "Latitude", color = "Source") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(validation_plot)

# 3. Final dataset map with basal area sizing and discrete phylogenetic coloring
final_plot <- ggplot(fg_final, aes(x = Longitude_final, y = Latitude_final)) +
  geom_point(aes(size = BasalArea_m2, color = Species_Name), alpha = 0.8) +
  scale_color_manual(
    values = final_colors,
    breaks = legend_order,
    name = "Species\n(phylogenetic order)"
  ) +
  scale_size_continuous(
    name = "Basal Area\n(m²)",
    range = c(0.5, 4),
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.2"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    title = "Final Forest Map: GPS Coordinates from Local Transformation",
    subtitle = paste("N =", nrow(fg_final), "stems, colored by phylogenetic distance (discrete), sized by basal area"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

print(final_plot)

# 4. Error visualization for non-outliers
good_points$lon_error <- good_points$Longitude - good_points$Longitude_converted
good_points$lat_error <- good_points$Latitude - good_points$Latitude_converted
good_points$total_error <- sqrt(good_points$lon_error^2 + good_points$lat_error^2)

error_plot <- ggplot(good_points, aes(x = PX, y = PY, color = total_error)) +
  geom_point(size = 2) +
  scale_color_viridis_c(name = "Error\n(degrees)") +
  labs(title = "Transformation Error by Plot Position (Non-outliers)",
       x = "PX (local)", y = "PY (local)") +
  theme_minimal() +
  coord_fixed()

print(error_plot)

# SUMMARY STATISTICS
cat("=== TRANSFORMATION SUMMARY ===\n")
cat("Plot is rotated", round(optimal_params[3] * 180 / pi, 1), "degrees from North-South orientation\n")
cat("Used ALL", nrow(fg_clean), "points for rotation calculation (including", n_outliers, "local outliers)\n")
cat("Final clean dataset:", nrow(fg_final), "points (local outliers removed)\n")
cat("Transformation accuracy on non-outliers: ±", round(max(rmse_lon, rmse_lat) * 111320, 2), "meters (approx.)\n")
cat("Scale factors suggest plot units are approximately", 
    round(1/mean(c(abs(optimal_params[1]), abs(optimal_params[2]))), 0), 
    "units per degree\n\n")

# BASAL AREA SUMMARY
cat("=== BASAL AREA SUMMARY ===\n")
total_basal_area <- sum(fg_final$BasalArea_m2, na.rm = TRUE)
cat("Total basal area (final dataset):", round(total_basal_area, 3), "m²\n")

species_basal_area <- fg_final %>%
  group_by(Species_Code, Species_Name) %>%
  summarise(
    n_stems = n(),
    total_basal_area_m2 = sum(BasalArea_m2, na.rm = TRUE),
    mean_dbh_cm = mean(DBH, na.rm = TRUE),
    mean_basal_area_m2 = mean(BasalArea_m2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_basal_area_m2))

print(species_basal_area)