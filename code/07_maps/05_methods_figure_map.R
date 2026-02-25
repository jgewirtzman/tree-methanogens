# ==============================================================================
# Methods Figure Map
# ==============================================================================
# Purpose: Study site methods figure map for publication.
#
# Pipeline stage: 4 — Visualization
#
# Outputs:
#   - methods map PNGs (to outputs/figures/)
# ==============================================================================

# Source upstream dependencies (must run from project root)
source("code/07_maps/01_forestgeo_alignment.R")
source("code/07_maps/02_spatial_interpolation.R")

# Overlay stem map on interpolated moisture map

# Create overlay plot with stem map on moisture/hillshade background
stem_moisture_overlay <- ggplot() +
  # Base hillshade layer for terrain effect
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation background
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.1) +
  # Moisture overlay (includes river influence)
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.6) +
  # ForestGEO trees - sized by basal area, colored by species
  geom_point(data = fg_final,
             aes(x = Longitude_final, y = Latitude_final,
                 size = BasalArea_m2, color = Species_Name),
             alpha = 0.8, stroke = 0.3) +
  # Research plots as reference points
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude),
             color = "white", fill = "black", size = 2, stroke = 1, shape = 21, alpha = 0.9) +
  # Plot ellipses for context
  geom_polygon(data = plot_tree_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "white", size = 0.8, alpha = 0.7, linetype = "dashed") +
  
  # Color and size scales
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "mako", direction = -1) +
  scale_color_manual(
    values = final_colors,
    breaks = legend_order,
    name = "Tree Species"
  ) +
  scale_size_continuous(
    name = "Basal Area\n(m²)",
    range = c(0.5, 4),
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.2"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_alpha_identity() +
  
  # Formatting
  coord_equal() +
  labs(
    title = "ForestGEO Stem Map Overlaid on Soil Moisture + Terrain",
    subtitle = paste("N =", nrow(fg_final), "trees on moisture gradient with research plots"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
    panel.background = element_rect(fill = "grey90")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1),
    size = guide_legend(ncol = 1)
  )

print(stem_moisture_overlay)

# Create a cleaner version focusing on tree distribution patterns
tree_distribution_plot <- ggplot() +
  # Simplified moisture background
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  # Trees with larger points for better visibility
  geom_point(data = fg_final, 
             aes(x = Longitude_final, y = Latitude_final, 
                 size = BasalArea_m2, color = Species_Name), 
             alpha = 0.9) +
  # Research plots as reference
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "black", fill = "white", size = 3, stroke = 1.5, shape = 21) +
  
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "viridis", direction = -1) +
  scale_color_manual(
    values = final_colors,
    breaks = legend_order,
    name = "Tree Species"
  ) +
  scale_size_continuous(
    name = "Basal Area\n(m²)",
    range = c(1, 5),
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.2")
  ) +
  
  coord_equal() +
  labs(
    title = "Tree Species Distribution Across Moisture Gradient",
    subtitle = "Tree size = basal area, colors = phylogenetic grouping",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    axis.text = element_text(size = 10)
  )

print(tree_distribution_plot)

# Create transformation method comparison on moisture background
transformation_comparison <- ggplot() +
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.5) +
  geom_point(data = fg_final, 
             aes(x = Longitude_final, y = Latitude_final, 
                 color = transformation_method, shape = Coordinates_Estimated), 
             size = 2, alpha = 0.8) +
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 2, stroke = 1, shape = 21) +
  
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "plasma", direction = -1) +
  scale_color_manual(
    values = c("Geodetic Transform" = "red", "Original GPS" = "blue"),
    name = "Coordinate Method"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Original GPS", "TRUE" = "Estimated GPS"),
    name = "GPS Source"
  ) +
  
  coord_equal() +
  labs(
    title = "Tree Coordinate Sources on Moisture Gradient",
    subtitle = "Red = geodetic transform from PX/PY, Blue = original GPS coordinates",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(transformation_comparison)

# Note: figS1 is saved at the end of this script after computing the
# clipped/aligned final_extended_plot (the publication version)

# Print summary of overlay
cat("=== STEM MAP - MOISTURE OVERLAY SUMMARY ===\n")
cat("Trees in ForestGEO dataset:", nrow(fg_final), "\n")
cat("Trees with geodetic transformation:", sum(fg_final$transformation_method == "Geodetic Transform"), "\n")
cat("Trees with original GPS:", sum(fg_final$transformation_method == "Original GPS"), "\n")
cat("Research plots for reference:", nrow(plots_data), "\n")
cat("Moisture interpolation grid points:", nrow(moisture_df), "\n\n")

# Check coordinate alignment
lon_overlap <- range(fg_final$Longitude_final)[1] >= range(moisture_df$Longitude)[1] && 
  range(fg_final$Longitude_final)[2] <= range(moisture_df$Longitude)[2]
lat_overlap <- range(fg_final$Latitude_final)[1] >= range(moisture_df$Latitude)[1] && 
  range(fg_final$Latitude_final)[2] <= range(moisture_df$Latitude)[2]

cat("Coordinate alignment check:\n")
cat("Trees longitude range:", round(range(fg_final$Longitude_final), 6), "\n")
cat("Moisture longitude range:", round(range(moisture_df$Longitude), 6), "\n")
cat("Trees latitude range:", round(range(fg_final$Latitude_final), 6), "\n")
cat("Moisture latitude range:", round(range(moisture_df$Latitude), 6), "\n")
cat("Longitude overlap:", lon_overlap, "\n")
cat("Latitude overlap:", lat_overlap, "\n\n")

if(!lon_overlap || !lat_overlap) {
  cat("WARNING: Tree coordinates may not fully overlap with moisture interpolation area\n")
  cat("Consider checking coordinate reference systems or expanding interpolation bounds\n")
}












# Alternative methods to extend moisture interpolation beyond convex hull
# Load additional required libraries
library(fields)
library(MBA)
library(mgcv)

# Method 1: Add boundary points to extend the interpolation domain
extend_moisture_with_boundary_points <- function(combined_moisture, stem_extent, buffer = 0.1) {
  
  # Get current moisture data extent
  moisture_extent <- list(
    lon = range(combined_moisture$Longitude),
    lat = range(combined_moisture$Latitude)
  )
  
  # Calculate expanded extent to cover stem map
  extended_extent <- list(
    lon = c(min(stem_extent$lon[1], moisture_extent$lon[1]), 
            max(stem_extent$lon[2], moisture_extent$lon[2])),
    lat = c(min(stem_extent$lat[1], moisture_extent$lat[1]), 
            max(stem_extent$lat[2], moisture_extent$lat[2]))
  )
  
  # Add buffer
  lon_buffer <- diff(extended_extent$lon) * buffer
  lat_buffer <- diff(extended_extent$lat) * buffer
  extended_extent$lon <- extended_extent$lon + c(-lon_buffer, lon_buffer)
  extended_extent$lat <- extended_extent$lat + c(-lat_buffer, lat_buffer)
  
  # Create boundary points with extrapolated moisture values
  # Use edge points to estimate boundary values
  boundary_points <- expand.grid(
    Longitude = c(extended_extent$lon[1], extended_extent$lon[2]),
    Latitude = seq(extended_extent$lat[1], extended_extent$lat[2], length.out = 10)
  )
  boundary_points <- rbind(boundary_points,
                           expand.grid(
                             Longitude = seq(extended_extent$lon[1], extended_extent$lon[2], length.out = 10),
                             Latitude = c(extended_extent$lat[1], extended_extent$lat[2])
                           )
  )
  
  # Estimate VWC for boundary points using inverse distance weighting
  boundary_points$VWC <- NA
  for(i in 1:nrow(boundary_points)) {
    distances <- sqrt((combined_moisture$Longitude - boundary_points$Longitude[i])^2 + 
                        (combined_moisture$Latitude - boundary_points$Latitude[i])^2)
    weights <- 1 / (distances^2 + 1e-6)  # Add small constant to avoid division by zero
    boundary_points$VWC[i] <- sum(combined_moisture$VWC * weights) / sum(weights)
  }
  
  # Combine original data with boundary points
  extended_data <- rbind(combined_moisture, boundary_points)
  
  return(list(data = extended_data, extent = extended_extent))
}

# Method 2: Use Thin Plate Splines (fields package)
create_tps_interpolation <- function(data, extent, resolution = 150) {
  
  # Create interpolation grid
  grid_x <- seq(extent$lon[1], extent$lon[2], length.out = resolution)
  grid_y <- seq(extent$lat[1], extent$lat[2], length.out = resolution)
  grid_locations <- expand.grid(Longitude = grid_x, Latitude = grid_y)
  
  # Fit thin plate spline
  tps_fit <- Tps(x = cbind(data$Longitude, data$Latitude), Y = data$VWC)
  
  # Predict on grid
  grid_locations$VWC <- predict(tps_fit, x = cbind(grid_locations$Longitude, grid_locations$Latitude))
  
  return(grid_locations)
}

# Method 3: Use GAM (mgcv package) 
create_gam_interpolation <- function(data, extent, resolution = 150) {
  
  # Fit GAM model
  gam_fit <- gam(VWC ~ s(Longitude, Latitude, k = 20), data = data)
  
  # Create prediction grid
  grid_x <- seq(extent$lon[1], extent$lon[2], length.out = resolution)
  grid_y <- seq(extent$lat[1], extent$lat[2], length.out = resolution)
  grid_locations <- expand.grid(Longitude = grid_x, Latitude = grid_y)
  
  # Predict on grid
  grid_locations$VWC <- predict(gam_fit, newdata = grid_locations)
  
  return(grid_locations)
}

# Get stem map extent
stem_extent <- list(
  lon = range(fg_final$Longitude_final, na.rm = TRUE),
  lat = range(fg_final$Latitude_final, na.rm = TRUE)
)

cat("=== EXTENDED INTERPOLATION METHODS ===\n")
cat("Stem map extent:\n")
cat("  Longitude:", round(stem_extent$lon, 6), "\n")
cat("  Latitude:", round(stem_extent$lat, 6), "\n\n")

# Method 1: Boundary point extension
cat("Method 1: Adding boundary points...\n")
extended_result <- extend_moisture_with_boundary_points(combined_moisture, stem_extent)
extended_moisture_akima <- interp(x = extended_result$data$Longitude,
                                  y = extended_result$data$Latitude,
                                  z = extended_result$data$VWC,
                                  xo = seq(extended_result$extent$lon[1], extended_result$extent$lon[2], length = 150),
                                  yo = seq(extended_result$extent$lat[1], extended_result$extent$lat[2], length = 150),
                                  duplicate = "mean")

extended_akima_df <- expand.grid(Longitude = extended_moisture_akima$x, 
                                 Latitude = extended_moisture_akima$y)
extended_akima_df$VWC <- as.vector(extended_moisture_akima$z)
extended_akima_df <- extended_akima_df[!is.na(extended_akima_df$VWC), ]
cat("  Grid points:", nrow(extended_akima_df), "\n")

# Method 2: Thin Plate Splines
cat("Method 2: Thin Plate Splines...\n")
tps_df <- create_tps_interpolation(combined_moisture, extended_result$extent)
tps_df <- tps_df[!is.na(tps_df$VWC), ]
cat("  Grid points:", nrow(tps_df), "\n")

# Method 3: GAM interpolation
cat("Method 3: GAM interpolation...\n")
gam_df <- create_gam_interpolation(combined_moisture, extended_result$extent)
gam_df <- gam_df[!is.na(gam_df$VWC), ]
cat("  Grid points:", nrow(gam_df), "\n\n")

# Check tree coverage for each method
check_coverage <- function(trees, interp_df, method_name) {
  in_bounds <- trees$Longitude_final >= min(interp_df$Longitude) & 
    trees$Longitude_final <= max(interp_df$Longitude) &
    trees$Latitude_final >= min(interp_df$Latitude) & 
    trees$Latitude_final <= max(interp_df$Latitude)
  
  cat(method_name, "coverage:", sum(in_bounds), "/", nrow(trees), 
      "(", round(sum(in_bounds)/nrow(trees)*100, 1), "%)\n")
  return(in_bounds)
}

cat("=== COVERAGE ANALYSIS ===\n")
akima_coverage <- check_coverage(fg_final, extended_akima_df, "Extended Akima")
tps_coverage <- check_coverage(fg_final, tps_df, "Thin Plate Splines")
gam_coverage <- check_coverage(fg_final, gam_df, "GAM")
cat("\n")

# Create comparison plots
create_method_plot <- function(interp_df, trees, title, method_name) {
  ggplot() +
    geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
    geom_point(data = trees, aes(x = Longitude_final, y = Latitude_final, size = BasalArea_m2), 
               alpha = 0.6, color = "black", stroke = 0.2) +
    geom_point(data = combined_moisture, aes(x = Longitude, y = Latitude), 
               color = "red", size = 1, alpha = 0.8) +
    scale_fill_viridis_c(name = "VWC %", option = "viridis", direction = -1) +
    scale_size_continuous(name = "Basal Area\n(m²)", range = c(0.3, 2)) +
    coord_equal() +
    labs(title = title,
         subtitle = paste("Red dots = moisture data points, Black dots = trees"),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "right")
}

# Create plots for each method
p1 <- create_method_plot(extended_akima_df, fg_final, 
                         "Extended Akima (Boundary Points)", "akima")
p2 <- create_method_plot(tps_df, fg_final, 
                         "Thin Plate Splines", "tps")
p3 <- create_method_plot(gam_df, fg_final, 
                         "GAM Interpolation", "gam")

print(p1)
print(p2)
print(p3)

# Choose the best method (highest coverage) for final overlay
best_coverage <- which.max(c(sum(akima_coverage), sum(tps_coverage), sum(gam_coverage)))
best_method <- c("Extended Akima", "Thin Plate Splines", "GAM")[best_coverage]
best_df <- list(extended_akima_df, tps_df, gam_df)[[best_coverage]]

cat("=== RECOMMENDED METHOD ===\n")
cat("Best coverage achieved by:", best_method, "\n")
cat("Using this method for final overlay plot...\n\n")

# Function to calculate minimum area bounding box with rotation
calculate_minimum_bounding_box <- function(points, buffer_pct = 0.01) {
  
  # Extract coordinates
  x <- points$x
  y <- points$y
  
  # Test rotation angles from 0 to 180 degrees (every 1 degree)
  angles <- seq(0, 179, by = 1) * pi / 180
  min_area <- Inf
  best_angle <- 0
  best_box <- NULL
  
  for(angle in angles) {
    # Rotate points
    cos_a <- cos(angle)
    sin_a <- sin(angle)
    
    x_rot <- x * cos_a - y * sin_a
    y_rot <- x * sin_a + y * cos_a
    
    # Calculate bounding box in rotated space
    x_range <- range(x_rot)
    y_range <- range(y_rot)
    
    # Calculate area
    area <- diff(x_range) * diff(y_range)
    
    # Keep track of minimum area
    if(area < min_area) {
      min_area <- area
      best_angle <- angle
      
      # Add buffer
      x_buffer <- diff(x_range) * buffer_pct
      y_buffer <- diff(y_range) * buffer_pct
      x_range_buffered <- x_range + c(-x_buffer, x_buffer)
      y_range_buffered <- y_range + c(-y_buffer, y_buffer)
      
      # Create corner points in rotated space
      corners_rot <- expand.grid(x = x_range_buffered, y = y_range_buffered)
      
      # Rotate back to original space
      corners_orig <- data.frame(
        x = corners_rot$x * cos(-angle) - corners_rot$y * sin(-angle),
        y = corners_rot$x * sin(-angle) + corners_rot$y * cos(-angle)
      )
      
      best_box <- list(
        corners = corners_orig,
        angle = best_angle * 180 / pi,
        area = min_area,
        rotated_ranges = list(x = x_range_buffered, y = y_range_buffered)
      )
    }
  }
  
  return(best_box)
}

# Function to check if points are inside rotated bounding box
point_in_rotated_box <- function(test_points, box_info) {
  
  angle <- box_info$angle * pi / 180
  cos_a <- cos(angle)
  sin_a <- sin(angle)
  
  # Rotate test points
  x_rot <- test_points$x * cos_a - test_points$y * sin_a
  y_rot <- test_points$x * sin_a + test_points$y * cos_a
  
  # Check if in rotated bounding box
  x_in <- x_rot >= box_info$rotated_ranges$x[1] & x_rot <= box_info$rotated_ranges$x[2]
  y_in <- y_rot >= box_info$rotated_ranges$y[1] & y_rot <= box_info$rotated_ranges$y[2]
  
  return(x_in & y_in)
}

# Get all feature coordinates
all_features_coords <- data.frame(
  x = c(fg_final$Longitude_final, plot_tree_ellipses$Longitude),
  y = c(fg_final$Latitude_final, plot_tree_ellipses$Latitude)
)

cat("=== CALCULATING MINIMUM AREA BOUNDING BOX ===\n")
cat("Total feature points:", nrow(all_features_coords), "\n")

# Calculate minimum bounding box
min_bbox <- calculate_minimum_bounding_box(all_features_coords, buffer_pct = 0.01)

cat("Optimal rotation angle:", round(min_bbox$angle, 2), "degrees\n")
cat("Minimum bounding box area:", round(min_bbox$area, 8), "\n")

# Compare with axis-aligned bounding box
axis_aligned_area <- diff(range(all_features_coords$x)) * diff(range(all_features_coords$y))
area_reduction <- (1 - min_bbox$area / axis_aligned_area) * 100

cat("Axis-aligned area:", round(axis_aligned_area, 8), "\n")
cat("Area reduction:", round(area_reduction, 1), "%\n\n")

# Clip moisture raster using rotated bounding box
moisture_coords <- data.frame(x = best_df$Longitude, y = best_df$Latitude)
inside_box <- point_in_rotated_box(moisture_coords, min_bbox)
clipped_moisture_df <- best_df[inside_box, ]

cat("=== CLIPPING MOISTURE RASTER ===\n")
cat("Original moisture grid points:", nrow(best_df), "\n")
cat("Clipped moisture grid points:", nrow(clipped_moisture_df), "\n")
cat("Points removed:", nrow(best_df) - nrow(clipped_moisture_df), "\n")
cat("Reduction:", round((1 - nrow(clipped_moisture_df)/nrow(best_df)) * 100, 1), "%\n\n")

# Create bounding box polygon for visualization
bbox_polygon <- data.frame(
  Longitude = c(min_bbox$corners$x, min_bbox$corners$x[1]),  # Close the polygon
  Latitude = c(min_bbox$corners$y, min_bbox$corners$y[1])
)

# Create clean publication figure (figS1) with clipped moisture overlay
# Three point types: All Trees (filled circle), Measured Trees (asterisk),
#   Research Plots (open circle) — matching the "Point Type" legend

# Add point_type column to each dataset for unified legend
fg_plot_data <- fg_final
fg_plot_data$point_type <- "All Trees"

measured_plot_data <- trees_with_plots
measured_plot_data$point_type <- "Measured Trees"

plots_plot_data <- plots_data
plots_plot_data$point_type <- "Research Plots"

final_extended_plot <- ggplot() +
  # Clipped moisture interpolation background
  geom_raster(data = clipped_moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  # Plot ellipses
  geom_polygon(data = plot_tree_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot),
               fill = NA, color = "black", linewidth = 0.8, alpha = 0.8) +
  # ForestGEO trees (all stems) — filled circles colored by species
  geom_point(data = fg_plot_data,
             aes(x = Longitude_final, y = Latitude_final,
                 size = BasalArea_m2, color = Species_Name, shape = point_type),
             alpha = 0.8, stroke = 0.3) +
  # Measured/sampled trees — asterisk symbol
  geom_point(data = measured_plot_data,
             aes(x = Longitude, y = Latitude, shape = point_type),
             color = "black", size = 3, stroke = 0.8, alpha = 0.9) +
  # Soil collars — white open circles
  geom_point(data = plots_plot_data,
             aes(x = Longitude, y = Latitude, shape = point_type),
             color = "black", fill = NA, size = 2.5, stroke = 1, alpha = 0.9) +

  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "mako", direction = -1) +
  scale_color_manual(values = final_colors, breaks = legend_order, name = "Tree Species") +
  scale_size_continuous(name = "Basal Area\n(m\u00b2)", range = c(0.5, 4),
                        breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
                        guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_shape_manual(name = "Point Type",
                     values = c("All Trees" = 16, "Measured Trees" = 8, "Research Plots" = 1),
                     guide = guide_legend(override.aes = list(size = c(3, 3, 3),
                                                               color = c("black", "black", "black"),
                                                               alpha = 1))) +

  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1),
    size = guide_legend(ncol = 1)
  )

print(final_extended_plot)

# Save figS1 publication figure
ggsave("outputs/figures/supplementary/figS1_moisture_overlay.png", final_extended_plot, width = 15, height = 9, dpi = 300)

# Store the best interpolation for future use
best_moisture_df <- best_df