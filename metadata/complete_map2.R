# First, identify matching trees between YM_trees_measured and fg_final datasets
# Assuming YM_trees_measured uses 'Label' and fg_final uses 'Tag'
matching_trees <- fg_final[fg_final$Tag %in% trees_data$Label, ]

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

# Get all feature coordinates (trees + ellipses for clipping)
all_features_coords <- data.frame(
  x = c(fg_final$Longitude_final, plot_tree_ellipses$Longitude),
  y = c(fg_final$Latitude_final, plot_tree_ellipses$Latitude)
)

# Calculate minimum bounding box for clipping
min_bbox <- calculate_minimum_bounding_box(all_features_coords, buffer_pct = 0.02)

# Get JUST tree coordinates for tree bounding box
tree_coords <- data.frame(
  x = fg_final$Longitude_final,
  y = fg_final$Latitude_final
)

# Calculate minimum bounding box for trees only (no buffer)
tree_bbox <- calculate_minimum_bounding_box(tree_coords, buffer_pct = 0.00)

# Create tree bounding box as proper rectangle corners only
tree_bbox_rect <- data.frame(
  x = tree_bbox$corners$x[c(1,2,4,3,1)],  # Connect corners in rectangle order
  y = tree_bbox$corners$y[c(1,2,4,3,1)]   # Connect corners in rectangle order
)

# Clip moisture raster using rotated bounding box
moisture_coords <- data.frame(x = best_df$Longitude, y = best_df$Latitude)
inside_box <- point_in_rotated_box(moisture_coords, min_bbox)
clipped_moisture_df <- best_df[inside_box, ]

# Create a dummy dataset for the shape legend using clipped data range
shape_legend_data <- data.frame(
  x = rep(mean(range(clipped_moisture_df$Longitude)), 3),
  y = rep(mean(range(clipped_moisture_df$Latitude)), 3),
  shape_type = factor(c("All Trees", "Measured Trees", "Research Plots"), 
                      levels = c("All Trees", "Measured Trees", "Research Plots"))
)

# Create the modified final plot
final_extended_plot_modified <- ggplot() +
  # Extended moisture interpolation background (now clipped)
  geom_raster(data = clipped_moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  
  # Plot ellipses from previous analysis
  geom_polygon(data = plot_tree_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "black", size = 1, alpha = 0.8) +
  
  # Thin grey rotated bounding box around tree points only
  geom_path(data = tree_bbox_rect, aes(x = x, y = y),
            color = "grey50", size = 0.5, alpha = 1) +
  
  # ForestGEO trees with thin dark grey outlines
  geom_point(data = fg_final, 
             aes(x = Longitude_final, y = Latitude_final, 
                 size = BasalArea_m2, color = Species_Name), 
             alpha = 0.8, stroke = 0.2, colour = "darkgrey") +
  geom_point(data = fg_final, 
             aes(x = Longitude_final, y = Latitude_final, 
                 size = BasalArea_m2, color = Species_Name), 
             alpha = 0.8, stroke = 0) +
  
  # Highlighted trees that match between datasets (stars) - same species colors
  geom_point(data = matching_trees, 
             aes(x = Longitude_final, y = Latitude_final, 
                 size = BasalArea_m2, color = Species_Name), 
             shape = 8, stroke = 1, alpha = 1) +
  
  # Research plots as open white circles
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = NA, size = 2, stroke = 1, shape = 1, alpha = 1) +
  
  # Blue Xs where soil moisture was collected
  geom_point(data = combined_moisture[combined_moisture$Type == "Soil", ], 
             aes(x = Longitude, y = Latitude), 
             color = "blue", size = 2, shape = 4, stroke = 1, alpha = 0.8) +
  
  # Add the shape legend layer (invisible points just for legend)
  geom_point(data = shape_legend_data, aes(x = x, y = y, shape = shape_type),
             alpha = 0, size = 3) +
  
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "mako", direction = -1) +
  scale_color_manual(values = final_colors, breaks = legend_order, name = "Tree Species",
                     guide = guide_legend(override.aes = list(shape = 16, size = 3, stroke = 0))) +
  scale_size_continuous(name = "Basal Area\n(m²)", range = c(0.5, 4), 
                        breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
                        guide = guide_legend(override.aes = list(shape = 16))) +
  scale_shape_manual(name = "Point Type",
                     values = c("All Trees" = 16, "Measured Trees" = 8, "Research Plots" = 1),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 3, 
                                                              color = c("black", "black", "black"),
                                                              fill = c("black", "black", NA)))) +
  
  coord_equal() +
  labs(
    #title = paste("Forest Plot with Soil Moisture Collection Points and Highlighted Trees"),
    #subtitle = paste("Blue X = soil moisture points, Yellow stars = YM measured trees, Grey outlines on all trees"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm")
  )

print(final_extended_plot_modified)

# Print summary of matching trees
cat("Summary of matching trees between datasets:\n")
cat("Total trees in fg_final:", nrow(fg_final), "\n")
cat("Total trees in YM_trees_measured:", nrow(trees_data), "\n")
cat("Matching trees (highlighted with stars):", nrow(matching_trees), "\n")

if(nrow(matching_trees) > 0) {
  cat("Matching tree species:\n")
  print(table(matching_trees$Species_Name))
} else {
  cat("No matching trees found. Check if column names are 'Tag' in fg_final and 'Label' in trees_data\n")
}

ggsave("final_extended_moisture_overlay_modified.png", final_extended_plot_modified, width = 10, height = 8, dpi = 300)
