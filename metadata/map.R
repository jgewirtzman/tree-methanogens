# Load required libraries
library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(akima)
library(readxl)
library(sf)
library(geosphere)

# Read the soil moisture data
data <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/soil_moisture_20201216.csv')

# Read the river data
river_data <- read_excel('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/River.xlsx')

# Read the research plots data
plots_data <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/plots.csv')

# Remove rows with missing VWC values
data_clean <- data[!is.na(data$VWC), ]

# Prepare river data - add VWC = 100% for river points
river_points <- data.frame(
  Longitude = river_data$Longitude,
  Latitude = river_data$Latitude,
  VWC = 100  # Assume 100% moisture at river
)

# Combine soil moisture data with river data for interpolation
combined_data <- rbind(
  data.frame(Longitude = data_clean$Longitude, 
             Latitude = data_clean$Latitude, 
             VWC = data_clean$VWC,
             Type = "Soil"),
  data.frame(Longitude = river_points$Longitude,
             Latitude = river_points$Latitude,
             VWC = river_points$VWC,
             Type = "River")
)

# Check data
print(paste("Number of soil data points:", nrow(data_clean)))
print(paste("Number of river data points:", nrow(river_points)))
print(paste("Number of research plots:", nrow(plots_data)))
print(paste("Combined data points:", nrow(combined_data)))
print(paste("VWC range:", round(min(combined_data$VWC), 2), "to", round(max(combined_data$VWC), 2)))

# Get extent of combined data and add small buffer
lon_range <- range(combined_data$Longitude)
lat_range <- range(combined_data$Latitude)
lon_buffer <- diff(lon_range) * 0.1  # 10% buffer
lat_buffer <- diff(lat_range) * 0.1

# Create interpolation using akima package with combined data (soil + river)
# Increase resolution for smoother map
interp_result <- interp(x = combined_data$Longitude, 
                        y = combined_data$Latitude, 
                        z = combined_data$VWC,
                        xo = seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, length = 100),
                        yo = seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, length = 100),
                        duplicate = "mean")

# Convert to dataframe for ggplot
interp_df <- expand.grid(Longitude = interp_result$x, 
                         Latitude = interp_result$y)
interp_df$VWC_predicted <- as.vector(interp_result$z)

# Remove NA values
interp_df <- interp_df[!is.na(interp_df$VWC_predicted), ]

# Function to create bounding ellipse for each plot
library(cluster)
create_plot_ellipses <- function(plots_data) {
  ellipses_list <- list()
  
  # Create unique plot identifier combining Site and Plot
  plots_data$Site_Plot <- paste(plots_data$Site, plots_data$Plot, sep = "_")
  
  for(site_plot in unique(plots_data$Site_Plot)) {
    plot_subset <- plots_data[plots_data$Site_Plot == site_plot, ]
    
    if(nrow(plot_subset) >= 3) {  # Need at least 3 points for ellipse
      # Calculate ellipse parameters
      center_lon <- mean(plot_subset$Longitude)
      center_lat <- mean(plot_subset$Latitude)
      
      # Calculate covariance matrix
      coords <- cbind(plot_subset$Longitude, plot_subset$Latitude)
      cov_matrix <- cov(coords)
      
      # Create ellipse points (99% confidence ellipse)
      angles <- seq(0, 2*pi, length.out = 100)
      eigenvals <- eigen(cov_matrix)$values
      eigenvecs <- eigen(cov_matrix)$vectors
      
      # Scale factor for 99% confidence ellipse
      scale_factor <- sqrt(qchisq(0.99, df = 2))  # Approximately 3.03
      
      ellipse_x <- center_lon + scale_factor * sqrt(eigenvals[1]) * cos(angles) * eigenvecs[1,1] + 
        scale_factor * sqrt(eigenvals[2]) * sin(angles) * eigenvecs[1,2]
      ellipse_y <- center_lat + scale_factor * sqrt(eigenvals[1]) * cos(angles) * eigenvecs[2,1] + 
        scale_factor * sqrt(eigenvals[2]) * sin(angles) * eigenvecs[2,2]
      
      ellipses_list[[site_plot]] <- data.frame(
        Longitude = ellipse_x,
        Latitude = ellipse_y,
        Site_Plot = site_plot,
        Site = plot_subset$Site[1],
        Plot = plot_subset$Plot[1]
      )
    }
  }
  
  do.call(rbind, ellipses_list)
}

# Create ellipses for plots
plot_ellipses <- create_plot_ellipses(plots_data)

# Function to calculate plot-level VWC three ways
calculate_plot_vwc <- function(plots_data, interp_df, interp_result, plot_ellipses) {
  results <- data.frame()
  
  # Create unique plot identifier
  plots_data$Site_Plot <- paste(plots_data$Site, plots_data$Plot, sep = "_")
  
  for(site_plot in unique(plots_data$Site_Plot)) {
    plot_subset <- plots_data[plots_data$Site_Plot == site_plot, ]
    
    # Method 1: Mean of all interpolated values within bounding ellipse
    if(site_plot %in% plot_ellipses$Site_Plot) {
      ellipse_coords <- plot_ellipses[plot_ellipses$Site_Plot == site_plot, ]
      
      # Create polygon from ellipse
      ellipse_poly <- st_polygon(list(as.matrix(ellipse_coords[, c("Longitude", "Latitude")])))
      ellipse_sf <- st_sfc(ellipse_poly, crs = 4326)
      
      # Convert interpolated data to sf points
      interp_sf <- st_as_sf(interp_df, coords = c("Longitude", "Latitude"), crs = 4326)
      
      # Find points within ellipse
      within_ellipse <- st_within(interp_sf, ellipse_sf, sparse = FALSE)
      method1_vwc <- mean(interp_df$VWC_predicted[within_ellipse], na.rm = TRUE)
    } else {
      method1_vwc <- NA
    }
    
    # Method 2: Mean of interpolated values at exact subplot locations
    method2_values <- c()
    for(i in 1:nrow(plot_subset)) {
      # Find closest interpolated point to subplot location
      distances <- sqrt((interp_df$Longitude - plot_subset$Longitude[i])^2 + 
                          (interp_df$Latitude - plot_subset$Latitude[i])^2)
      closest_idx <- which.min(distances)
      method2_values <- c(method2_values, interp_df$VWC_predicted[closest_idx])
    }
    method2_vwc <- mean(method2_values, na.rm = TRUE)
    
    # Method 3: Mean of interpolated values within 5m buffer of each subplot
    method3_values <- c()
    for(i in 1:nrow(plot_subset)) {
      # Convert 5m to degrees (approximate: 1 degree ≈ 111,000m at equator)
      buffer_deg <- 5 / 111000
      
      # Find points within 5m buffer
      distances <- sqrt((interp_df$Longitude - plot_subset$Longitude[i])^2 + 
                          (interp_df$Latitude - plot_subset$Latitude[i])^2)
      within_buffer <- distances <= buffer_deg
      
      if(sum(within_buffer) > 0) {
        method3_values <- c(method3_values, mean(interp_df$VWC_predicted[within_buffer], na.rm = TRUE))
      }
    }
    method3_vwc <- mean(method3_values, na.rm = TRUE)
    
    # Store results
    results <- rbind(results, data.frame(
      Site_Plot = site_plot,
      Site = plot_subset$Site[1],
      Plot = plot_subset$Plot[1],
      n_subplots = nrow(plot_subset),
      Method1_Ellipse_Mean = round(method1_vwc, 2),
      Method2_Subplot_Points = round(method2_vwc, 2),
      Method3_5m_Buffer = round(method3_vwc, 2)
    ))
  }
  
  results
}

# Calculate plot-level VWC using all three methods
plot_vwc_results <- calculate_plot_vwc(plots_data, interp_df, interp_result, plot_ellipses)

# Create the final map
p1 <- ggplot() +
  geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = VWC_predicted)) +
  geom_polygon(data = plot_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "black", size = 1, alpha = 0.8) +
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "black", fill = "white", size = 2, stroke = 1, shape = 21) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "viridis", direction = -1) +
  coord_equal() +
  labs(title = "Soil Moisture Map with Research Plots and Bounding Ellipses",
       subtitle = "Black ellipses: plot boundaries, White dots: subplot locations",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2))

print(p1)

# Create contour version
p2 <- ggplot() +
  geom_raster(data = interp_df, aes(x = Longitude, y = Latitude, fill = VWC_predicted), alpha = 0.7) +
  geom_contour(data = interp_df, aes(x = Longitude, y = Latitude, z = VWC_predicted),
               color = "white", size = 0.5, bins = 10) +
  geom_polygon(data = plot_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "black", size = 1, alpha = 0.8) +
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "black", fill = "white", size = 2, stroke = 1, shape = 21) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "turbo", direction = -1) +
  coord_equal() +
  labs(title = "Soil Moisture Contour Map with Research Plots",
       subtitle = "Black ellipses: plot boundaries, White contours: moisture zones",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2))

print(p2)

# Save plots
ggsave("soil_moisture_plots_with_ellipses.png", p1, width = 12, height = 8, dpi = 300)
ggsave("soil_moisture_contour_with_ellipses.png", p2, width = 12, height = 8, dpi = 300)

# Print results
cat("Plot-Level VWC Analysis Results:\n")
cat("Method 1: Mean of all interpolated values within bounding ellipse\n")
cat("Method 2: Mean of interpolated values at exact subplot locations\n")
cat("Method 3: Mean of interpolated values within 5m buffer of each subplot\n\n")

print(plot_vwc_results)

# Save results to CSV
write.csv(plot_vwc_results, "plot_level_vwc_analysis.csv", row.names = FALSE)

# Summary statistics by method
cat("\nSummary Statistics by Method:\n")
summary_stats <- plot_vwc_results %>%
  summarise(
    Method1_Mean = round(mean(Method1_Ellipse_Mean, na.rm = TRUE), 2),
    Method1_SD = round(sd(Method1_Ellipse_Mean, na.rm = TRUE), 2),
    Method2_Mean = round(mean(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method2_SD = round(sd(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method3_Mean = round(mean(Method3_5m_Buffer, na.rm = TRUE), 2),
    Method3_SD = round(sd(Method3_5m_Buffer, na.rm = TRUE), 2)
  )

print(summary_stats)

# Correlation between methods
cat("\nCorrelations between methods:\n")
cor_matrix <- cor(plot_vwc_results[, c("Method1_Ellipse_Mean", "Method2_Subplot_Points", "Method3_5m_Buffer")], 
                  use = "complete.obs")
print(round(cor_matrix, 3))


# Load required libraries
library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(akima)
library(readxl)

# Set file path base
base_path <- '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data'

# Read all data files
soil_moisture <- read.csv(file.path(base_path, 'soil_moisture_20201216.csv'))
plots_data <- read.csv(file.path(base_path, 'plots.csv'))
wells_data <- read.csv(file.path(base_path, 'wells.csv'))
trees_data <- read.csv(file.path(base_path, 'YM_trees_measured.csv.csv'))
river_data <- read_excel(file.path(base_path, 'River.xlsx'))

# Extract elevation data from each file
# Soil moisture points
soil_elev <- data.frame(
  Longitude = soil_moisture$Longitude,
  Latitude = soil_moisture$Latitude,
  Elevation = soil_moisture$Elevation,
  Source = "Soil_Moisture",
  stringsAsFactors = FALSE
)

# Research plots
plots_elev <- data.frame(
  Longitude = plots_data$Longitude,
  Latitude = plots_data$Latitude,
  Elevation = plots_data$Elevation,
  Source = "Research_Plots",
  stringsAsFactors = FALSE
)

# Wells
wells_elev <- data.frame(
  Longitude = wells_data$Longitude,
  Latitude = wells_data$Latitude,
  Elevation = wells_data$Elevation,
  Source = "Wells",
  stringsAsFactors = FALSE
)

# Trees
trees_elev <- data.frame(
  Longitude = trees_data$Longitude,
  Latitude = trees_data$Latitude,
  Elevation = trees_data$Elevation,
  Source = "Trees",
  stringsAsFactors = FALSE
)

# Check for outliers in tree data
cat("Tree coordinate ranges:\n")
cat("Longitude range:", range(trees_elev$Longitude, na.rm = TRUE), "\n")
cat("Latitude range:", range(trees_elev$Latitude, na.rm = TRUE), "\n")

# Calculate median coordinates for trees
median_lon <- median(trees_elev$Longitude, na.rm = TRUE)
median_lat <- median(trees_elev$Latitude, na.rm = TRUE)

# Calculate distances from median center
trees_elev$dist_from_center <- sqrt((trees_elev$Longitude - median_lon)^2 + 
                                      (trees_elev$Latitude - median_lat)^2)

# Find potential outliers (points > 2 standard deviations from center)
dist_threshold <- mean(trees_elev$dist_from_center, na.rm = TRUE) + 2 * sd(trees_elev$dist_from_center, na.rm = TRUE)
outliers <- trees_elev[trees_elev$dist_from_center > dist_threshold, ]

if(nrow(outliers) > 0) {
  cat("Potential tree outliers found:\n")
  print(outliers)
  
  # Remove outliers
  trees_elev <- trees_elev[trees_elev$dist_from_center <= dist_threshold, ]
  cat("Removed", nrow(outliers), "outlier tree point(s)\n")
  cat("Remaining tree points:", nrow(trees_elev), "\n")
}

# Remove the distance column for final dataset
trees_elev <- trees_elev[, c("Longitude", "Latitude", "Elevation", "Source")]

# River points
river_elev <- data.frame(
  Longitude = river_data$Longitude,
  Latitude = river_data$Latitude,
  Elevation = river_data$Elevation,
  Source = "River",
  stringsAsFactors = FALSE
)

# Combine all elevation data
all_elevation <- rbind(soil_elev, plots_elev, wells_elev, trees_elev, river_elev)

# Remove any rows with missing elevation data
all_elevation <- all_elevation[!is.na(all_elevation$Elevation), ]

# Print summary
cat("Elevation Data Summary:\n")
print(table(all_elevation$Source))
cat("Total elevation points:", nrow(all_elevation), "\n")
cat("Elevation range:", round(min(all_elevation$Elevation), 2), "to", round(max(all_elevation$Elevation), 2), "m\n")

# Get extent of data and add buffer
lon_range <- range(all_elevation$Longitude)
lat_range <- range(all_elevation$Latitude)
lon_buffer <- diff(lon_range) * 0.1
lat_buffer <- diff(lat_range) * 0.1

# Create interpolation using akima
elev_interp <- interp(x = all_elevation$Longitude, 
                      y = all_elevation$Latitude, 
                      z = all_elevation$Elevation,
                      xo = seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, length = 150),
                      yo = seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, length = 150),
                      duplicate = "mean")

# Convert to dataframe for ggplot
elev_df <- expand.grid(Longitude = elev_interp$x, 
                       Latitude = elev_interp$y)
elev_df$Elevation <- as.vector(elev_interp$z)

# Remove NA values
elev_df <- elev_df[!is.na(elev_df$Elevation), ]

# Create raster for hillshade calculation
elev_raster <- raster(elev_interp)

# Calculate hillshade
# Sun angle: azimuth = 315° (northwest), altitude = 45°
hillshade <- hillShade(terrain(elev_raster, opt = "slope"),
                       terrain(elev_raster, opt = "aspect"),
                       angle = 45, direction = 315)

# Convert hillshade to dataframe
hillshade_df <- as.data.frame(rasterToPoints(hillshade))
names(hillshade_df) <- c("Longitude", "Latitude", "Hillshade")

# Create hillshade elevation map
p1 <- ggplot() +
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, fill = Hillshade), alpha = 0.6) +
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude, fill = Elevation), alpha = 0.7) +
  geom_point(data = all_elevation, aes(x = Longitude, y = Latitude, color = Source), 
             size = 1.5, alpha = 0.9) +
  scale_fill_gradient(name = "Elevation\n(meters)", low = "black", high = "white") +
  scale_color_manual(name = "Data Source",
                     values = c("Soil_Moisture" = "red",
                                "Research_Plots" = "blue", 
                                "Wells" = "green",
                                "Trees" = "orange",
                                "River" = "cyan")) +
  coord_equal() +
  labs(title = "Hillshade Elevation Map",
       subtitle = "3D terrain effect: Black = Low, White = High elevation",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p1)

# Create contour version
p2 <- ggplot() +
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude, fill = Elevation), alpha = 0.8) +
  geom_contour(data = elev_df, aes(x = Longitude, y = Latitude, z = Elevation),
               color = "white", size = 0.5, bins = 15) +
  geom_point(data = all_elevation, aes(x = Longitude, y = Latitude, color = Source), 
             size = 1.5, alpha = 0.8) +
  scale_fill_gradient(name = "Elevation\n(meters)", low = "black", high = "white") +
  scale_color_manual(name = "Data Source",
                     values = c("Soil_Moisture" = "red",
                                "Research_Plots" = "blue", 
                                "Wells" = "green",
                                "Trees" = "orange",
                                "River" = "cyan")) +
  coord_equal() +
  labs(title = "Elevation Map with Contour Lines",
       subtitle = "Black = Low, White = High, White lines = Elevation contours",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p2)

# Create clean hillshade version without points
p3 <- ggplot() +
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude, fill = Elevation), alpha = 0.8) +
  scale_fill_gradient(name = "Elevation\n(meters)", low = "black", high = "white") +
  scale_alpha_identity() +
  coord_equal() +
  labs(title = "Clean Hillshade Elevation Map",
       subtitle = "3D terrain effect: Black = Low, White = High elevation",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p3)

# Save plots
ggsave("elevation_hillshade_with_points.png", p1, width = 12, height = 8, dpi = 300)
ggsave("elevation_map_with_contours.png", p2, width = 12, height = 8, dpi = 300)
ggsave("elevation_hillshade_clean.png", p3, width = 12, height = 8, dpi = 300)

# Print detailed statistics
cat("\nDetailed Elevation Statistics by Source:\n")
elev_stats <- all_elevation %>%
  group_by(Source) %>%
  summarise(
    n_points = n(),
    min_elev = round(min(Elevation), 2),
    max_elev = round(max(Elevation), 2),
    mean_elev = round(mean(Elevation), 2),
    sd_elev = round(sd(Elevation), 2),
    .groups = 'drop'
  )
print(elev_stats)

# Overall statistics
cat("\nOverall Elevation Statistics:\n")
cat("Total data points:", nrow(all_elevation), "\n")
cat("Elevation range:", round(min(all_elevation$Elevation), 2), "to", round(max(all_elevation$Elevation), 2), "meters\n")
cat("Mean elevation:", round(mean(all_elevation$Elevation), 2), "meters\n")
cat("Standard deviation:", round(sd(all_elevation$Elevation), 2), "meters\n")
cat("Study area extent:\n")
cat("  Longitude:", round(lon_range, 5), "\n")
cat("  Latitude:", round(lat_range, 5), "\n")

# Save elevation data
write.csv(all_elevation, "combined_elevation_data.csv", row.names = FALSE)

# Load the moisture interpolation we created earlier
# This assumes you've already run the soil moisture interpolation script
# If not, we'll need to recreate it here with river data included

# Read moisture data and prepare it (same as before)
moisture_clean <- soil_moisture[!is.na(soil_moisture$VWC), ]

# Prepare river data for moisture interpolation  
river_moisture <- data.frame(
  Longitude = river_data$Longitude,
  Latitude = river_data$Latitude,
  VWC = 100  # River at 100% moisture
)

# Combine soil + river moisture data (same as the earlier analysis)
combined_moisture <- rbind(
  data.frame(Longitude = moisture_clean$Longitude, 
             Latitude = moisture_clean$Latitude, 
             VWC = moisture_clean$VWC),
  river_moisture
)

# Create moisture interpolation on same grid as elevation for perfect overlay
moisture_interp <- interp(x = combined_moisture$Longitude, 
                          y = combined_moisture$Latitude, 
                          z = combined_moisture$VWC,
                          xo = elev_interp$x,  # Same grid as elevation
                          yo = elev_interp$y,  # Same grid as elevation
                          duplicate = "mean")

# Convert to dataframe
moisture_df <- expand.grid(Longitude = moisture_interp$x, 
                           Latitude = moisture_interp$y)
moisture_df$VWC <- as.vector(moisture_interp$z)
moisture_df <- moisture_df[!is.na(moisture_df$VWC), ]

# Create combined hillshade + moisture map (same moisture as earlier analysis)
p4 <- ggplot() +
  # Base hillshade layer for 3D terrain effect
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation for context
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.2) +
  # Moisture overlay (includes river influence)
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "viridis", direction = -1) +
  scale_alpha_identity() +
  coord_equal() +
  labs(title = "Soil Moisture (with River) on Hillshade Terrain",
       subtitle = "3D terrain with transparent moisture overlay including river influence",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p4)

# Create version with research plots overlaid
p5 <- ggplot() +
  # Base hillshade layer
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation background
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.15) +
  # Moisture overlay (same as our earlier analysis with river)
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.8) +
  # Research plots as black dots with white outline
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 2.5, stroke = 1, shape = 21, alpha = 0.9) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "viridis", direction = -1) +
  scale_alpha_identity() +
  coord_equal() +
  labs(title = "Soil Moisture + Hillshade + Research Plots",
       subtitle = "Complete analysis: terrain relief + moisture patterns + study locations",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p5)

# Save the combined maps
ggsave("moisture_hillshade_overlay.png", p4, width = 12, height = 8, dpi = 300)
ggsave("moisture_hillshade_plots.png", p5, width = 12, height = 8, dpi = 300)

# Assign each tree to nearest research plot
assign_trees_to_plots <- function(trees_data, plots_data) {
  # First remove outliers from trees_data (same method as elevation)
  trees_clean <- trees_data
  
  # Calculate median coordinates for trees
  median_lon <- median(trees_clean$Longitude, na.rm = TRUE)
  median_lat <- median(trees_clean$Latitude, na.rm = TRUE)
  
  # Calculate distances from median center
  trees_clean$dist_from_center <- sqrt((trees_clean$Longitude - median_lon)^2 + 
                                         (trees_clean$Latitude - median_lat)^2)
  
  # Find potential outliers (points > 2 standard deviations from center)
  dist_threshold <- mean(trees_clean$dist_from_center, na.rm = TRUE) + 2 * sd(trees_clean$dist_from_center, na.rm = TRUE)
  outliers <- trees_clean[trees_clean$dist_from_center > dist_threshold, ]
  
  if(nrow(outliers) > 0) {
    cat("Removing", nrow(outliers), "outlier tree(s):\n")
    print(outliers[, c("Label", "Longitude", "Latitude")])
    
    # Remove outliers
    trees_clean <- trees_clean[trees_clean$dist_from_center <= dist_threshold, ]
    cat("Remaining trees after outlier removal:", nrow(trees_clean), "\n")
  }
  
  # Remove the distance column
  trees_clean <- trees_clean[, !names(trees_clean) %in% "dist_from_center"]
  
  # Now assign remaining trees to plots
  trees_assigned <- trees_clean
  trees_assigned$Nearest_Plot <- NA
  trees_assigned$Nearest_Site <- NA
  trees_assigned$Distance_to_Plot <- NA
  
  for(i in 1:nrow(trees_clean)) {
    # Calculate distances to all plot points
    distances <- sqrt((trees_clean$Longitude[i] - plots_data$Longitude)^2 + 
                        (trees_clean$Latitude[i] - plots_data$Latitude)^2)
    
    # Find nearest plot
    nearest_idx <- which.min(distances)
    
    # Assign tree to nearest plot
    trees_assigned$Nearest_Plot[i] <- paste(plots_data$Site[nearest_idx], 
                                            plots_data$Plot[nearest_idx], 
                                            sep = "_")
    trees_assigned$Nearest_Site[i] <- plots_data$Site[nearest_idx]
    trees_assigned$Distance_to_Plot[i] <- round(distances[nearest_idx] * 111000, 2) # Convert to meters
  }
  
  return(trees_assigned)
}

# Assign trees to plots (after removing outliers)
trees_with_plots <- assign_trees_to_plots(trees_data, plots_data)

# Print assignment summary
cat("Tree-Plot Assignment Summary:\n")
assignment_summary <- trees_with_plots %>%
  group_by(Nearest_Site, Nearest_Plot) %>%
  summarise(
    n_trees = n(),
    mean_distance_m = round(mean(Distance_to_Plot), 2),
    max_distance_m = round(max(Distance_to_Plot), 2),
    .groups = 'drop'
  )
print(assignment_summary)

# Create final map with trees as small black dots
p6 <- ggplot() +
  # Base hillshade layer
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation background
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.15) +
  # Moisture overlay
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  # Trees as small black dots
  geom_point(data = trees_with_plots, aes(x = Longitude, y = Latitude), 
             color = "black", size = 1, alpha = 0.8) +
  # Research plots as black dots with white outline
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 3, stroke = 1, shape = 21, alpha = 0.9) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "viridis", direction = -1) +
  scale_alpha_identity() +
  coord_equal() +
  labs(title = "Moisture + Terrain + Research Plots + Trees",
       subtitle = "Small black dots: trees, Large dots: research plots",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p6)

# Function to create ellipses including both plot points and assigned trees
create_plot_tree_ellipses <- function(plots_data, trees_with_plots) {
  ellipses_list <- list()
  
  # Create unique plot identifier
  plots_data$Site_Plot <- paste(plots_data$Site, plots_data$Plot, sep = "_")
  
  for(site_plot in unique(plots_data$Site_Plot)) {
    # Get plot points for this site_plot
    plot_subset <- plots_data[plots_data$Site_Plot == site_plot, ]
    
    # Get trees assigned to this site_plot
    tree_subset <- trees_with_plots[trees_with_plots$Nearest_Plot == site_plot, ]
    
    # Combine plot points and trees
    if(nrow(tree_subset) > 0) {
      combined_points <- rbind(
        data.frame(Longitude = plot_subset$Longitude, 
                   Latitude = plot_subset$Latitude,
                   Type = "Plot"),
        data.frame(Longitude = tree_subset$Longitude,
                   Latitude = tree_subset$Latitude, 
                   Type = "Tree")
      )
    } else {
      combined_points <- data.frame(Longitude = plot_subset$Longitude, 
                                    Latitude = plot_subset$Latitude,
                                    Type = "Plot")
    }
    
    if(nrow(combined_points) >= 3) {  # Need at least 3 points for ellipse
      # Calculate ellipse parameters
      center_lon <- mean(combined_points$Longitude)
      center_lat <- mean(combined_points$Latitude)
      
      # Calculate covariance matrix
      coords <- cbind(combined_points$Longitude, combined_points$Latitude)
      cov_matrix <- cov(coords)
      
      # Create ellipse points (95% confidence ellipse)
      angles <- seq(0, 2*pi, length.out = 100)
      eigenvals <- eigen(cov_matrix)$values
      eigenvecs <- eigen(cov_matrix)$vectors
      
      # Scale factor for 95% confidence ellipse
      scale_factor <- sqrt(qchisq(0.95, df = 2))  # Approximately 2.45
      
      ellipse_x <- center_lon + scale_factor * sqrt(eigenvals[1]) * cos(angles) * eigenvecs[1,1] + 
        scale_factor * sqrt(eigenvals[2]) * sin(angles) * eigenvecs[1,2]
      ellipse_y <- center_lat + scale_factor * sqrt(eigenvals[1]) * cos(angles) * eigenvecs[2,1] + 
        scale_factor * sqrt(eigenvals[2]) * sin(angles) * eigenvecs[2,2]
      
      ellipses_list[[site_plot]] <- data.frame(
        Longitude = ellipse_x,
        Latitude = ellipse_y,
        Site_Plot = site_plot,
        Site = plot_subset$Site[1],
        Plot = plot_subset$Plot[1],
        n_points = nrow(plot_subset),
        n_trees = nrow(tree_subset)
      )
    }
  }
  
  do.call(rbind, ellipses_list)
}

# Create new ellipses including trees
plot_tree_ellipses <- create_plot_tree_ellipses(plots_data, trees_with_plots)

# Print ellipse summary
cat("Plot-Tree Ellipse Summary:\n")
if(nrow(plot_tree_ellipses) > 0) {
  ellipse_summary <- plot_tree_ellipses %>%
    group_by(Site_Plot, Site, Plot, n_points, n_trees) %>%
    summarise(.groups = 'drop') %>%
    distinct()
  print(ellipse_summary)
} else {
  cat("No ellipses could be created (need at least 3 points per plot)\n")
}

# Create map with new ellipses including trees
p7 <- ggplot() +
  # Base hillshade layer
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation background
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.15) +
  # Moisture overlay
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  # Plot ellipses (now including trees)
  geom_polygon(data = plot_tree_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "black", size = 1, alpha = 0.8) +
  # Trees as small black dots
  geom_point(data = trees_with_plots, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 1, stroke = 1, shape = 21, alpha = 0.9) +
  # Research plots as larger black dots with white outline
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 3, stroke = 1, shape = 21, alpha = 0.7) +
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "mako", direction = -1) +
  scale_alpha_identity() +
  coord_equal() +
  labs(
    #title = "Moisture + Terrain + Plot Ellipses Including Trees",
    #subtitle = "Black ellipses: 95% confidence bounds around plots + assigned trees",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p7)

# Save the tree assignment maps
ggsave("moisture_hillshade_plots_trees.png", p6, width = 14, height = 8, dpi = 300)
ggsave("moisture_hillshade_ellipses_with_trees.png", p7, width = 14, height = 8, dpi = 300)

# Calculate plot-level VWC using three methods with new ellipses including trees
calculate_plot_vwc_with_trees <- function(plots_data, trees_with_plots, moisture_df, moisture_interp, plot_tree_ellipses) {
  results <- data.frame()
  
  # Create unique plot identifier
  plots_data$Site_Plot <- paste(plots_data$Site, plots_data$Plot, sep = "_")
  
  for(site_plot in unique(plots_data$Site_Plot)) {
    plot_subset <- plots_data[plots_data$Site_Plot == site_plot, ]
    tree_subset <- trees_with_plots[trees_with_plots$Nearest_Plot == site_plot, ]
    
    # Method 1: Mean of all interpolated values within new bounding ellipse (includes trees)
    if(site_plot %in% plot_tree_ellipses$Site_Plot) {
      ellipse_coords <- plot_tree_ellipses[plot_tree_ellipses$Site_Plot == site_plot, ]
      
      # Create polygon from ellipse
      ellipse_poly <- st_polygon(list(as.matrix(ellipse_coords[, c("Longitude", "Latitude")])))
      ellipse_sf <- st_sfc(ellipse_poly, crs = 4326)
      
      # Convert interpolated data to sf points
      moisture_sf <- st_as_sf(moisture_df, coords = c("Longitude", "Latitude"), crs = 4326)
      
      # Find points within ellipse
      within_ellipse <- st_within(moisture_sf, ellipse_sf, sparse = FALSE)
      method1_vwc <- mean(moisture_df$VWC[within_ellipse], na.rm = TRUE)
    } else {
      method1_vwc <- NA
    }
    
    # Method 2: Mean of interpolated values at exact subplot locations (plot points only)
    method2_values <- c()
    for(i in 1:nrow(plot_subset)) {
      # Find closest interpolated point to subplot location
      distances <- sqrt((moisture_df$Longitude - plot_subset$Longitude[i])^2 + 
                          (moisture_df$Latitude - plot_subset$Latitude[i])^2)
      closest_idx <- which.min(distances)
      method2_values <- c(method2_values, moisture_df$VWC[closest_idx])
    }
    method2_vwc <- mean(method2_values, na.rm = TRUE)
    
    # Method 3: Mean of interpolated values within 5m buffer of each subplot (plot points only)
    method3_values <- c()
    for(i in 1:nrow(plot_subset)) {
      # Convert 5m to degrees (approximate: 1 degree ≈ 111,000m at equator)
      buffer_deg <- 5 / 111000
      
      # Find points within 5m buffer
      distances <- sqrt((moisture_df$Longitude - plot_subset$Longitude[i])^2 + 
                          (moisture_df$Latitude - plot_subset$Latitude[i])^2)
      within_buffer <- distances <= buffer_deg
      
      if(sum(within_buffer) > 0) {
        method3_values <- c(method3_values, mean(moisture_df$VWC[within_buffer], na.rm = TRUE))
      }
    }
    method3_vwc <- mean(method3_values, na.rm = TRUE)
    
    # Store results
    results <- rbind(results, data.frame(
      Site_Plot = site_plot,
      Site = plot_subset$Site[1],
      Plot = plot_subset$Plot[1],
      n_subplots = nrow(plot_subset),
      n_trees = nrow(tree_subset),
      Method1_Ellipse_with_Trees = round(method1_vwc, 2),
      Method2_Subplot_Points = round(method2_vwc, 2),
      Method3_5m_Buffer = round(method3_vwc, 2)
    ))
  }
  
  results
}

# Calculate plot-level VWC using all three methods with new ellipses
plot_vwc_final <- calculate_plot_vwc_with_trees(plots_data, trees_with_plots, moisture_df, moisture_interp, plot_tree_ellipses)

# Print results
cat("\nFinal Plot-Level VWC Analysis Results (with tree-inclusive ellipses):\n")
cat("Method 1: Mean of all interpolated values within ellipse (including trees)\n")
cat("Method 2: Mean of interpolated values at exact subplot locations\n")
cat("Method 3: Mean of interpolated values within 5m buffer of each subplot\n\n")

print(plot_vwc_final)

# Save results to CSV
write.csv(plot_vwc_final, "final_plot_level_vwc_analysis.csv", row.names = FALSE)

# Summary statistics by method
cat("\nSummary Statistics by Method:\n")
summary_stats_final <- plot_vwc_final %>%
  summarise(
    Method1_Mean = round(mean(Method1_Ellipse_with_Trees, na.rm = TRUE), 2),
    Method1_SD = round(sd(Method1_Ellipse_with_Trees, na.rm = TRUE), 2),
    Method2_Mean = round(mean(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method2_SD = round(sd(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method3_Mean = round(mean(Method3_5m_Buffer, na.rm = TRUE), 2),
    Method3_SD = round(sd(Method3_5m_Buffer, na.rm = TRUE), 2)
  )

print(summary_stats_final)

# Correlation between methods
cat("\nCorrelations between methods:\n")
cor_matrix_final <- cor(plot_vwc_final[, c("Method1_Ellipse_with_Trees", "Method2_Subplot_Points", "Method3_5m_Buffer")], 
                        use = "complete.obs")
print(round(cor_matrix_final, 3))

# Summary by site type
cat("\nVWC Summary by Site Type:\n")
site_summary <- plot_vwc_final %>%
  group_by(Site) %>%
  summarise(
    n_plots = n(),
    total_trees = sum(n_trees),
    Method1_Mean = round(mean(Method1_Ellipse_with_Trees, na.rm = TRUE), 2),
    Method1_SD = round(sd(Method1_Ellipse_with_Trees, na.rm = TRUE), 2),
    Method2_Mean = round(mean(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method2_SD = round(sd(Method2_Subplot_Points, na.rm = TRUE), 2),
    Method3_Mean = round(mean(Method3_5m_Buffer, na.rm = TRUE), 2),
    Method3_SD = round(sd(Method3_5m_Buffer, na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(site_summary)

# Save tree assignment data
write.csv(trees_with_plots, "trees_assigned_to_plots.csv", row.names = FALSE)



# Calculate ellipse centers for all plots
ellipse_centers <- plot_tree_ellipses %>%
  group_by(Site_Plot, Site) %>%
  summarise(
    center_lon = mean(Longitude),
    center_lat = mean(Latitude),
    .groups = 'drop'
  )

# Get all ellipses for each habitat type
wetland_ellipses <- ellipse_centers[ellipse_centers$Site == "Wetland", ]
intermediate_ellipses <- ellipse_centers[ellipse_centers$Site == "Intermediate", ]
upland_ellipses <- ellipse_centers[ellipse_centers$Site == "Upland", ]

# Calculate smart label positions
lon_range <- range(plot_tree_ellipses$Longitude)
lat_range <- range(plot_tree_ellipses$Latitude)
offset_x <- diff(lon_range) * 0.2
offset_y <- diff(lat_range) * 0.15

# Position labels as requested
wetland_label_pos <- data.frame(
  x = mean(wetland_ellipses$center_lon),  # Halfway between the two wetland plots
  y = mean(wetland_ellipses$center_lat) - offset_y/2  # Nudged downward
)

intermediate_label_pos <- data.frame(
  x = mean(range(plot_tree_ellipses$Longitude)),  # Keep as is
  y = max(lat_range) + offset_y
)

upland_label_pos <- data.frame(
  x = mean(upland_ellipses$center_lon),  # Halfway between the two upland plots
  y = mean(upland_ellipses$center_lat) - offset_y*2  # Nudged downward
)

# Create the annotated map
p7 <- ggplot() +
  # Base hillshade layer
  geom_raster(data = hillshade_df, aes(x = Longitude, y = Latitude, alpha = Hillshade)) +
  # Subtle elevation background
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude), fill = "grey50", alpha = 0.15) +
  # Moisture overlay
  geom_raster(data = moisture_df, aes(x = Longitude, y = Latitude, fill = VWC), alpha = 0.7) +
  # Plot ellipses
  geom_polygon(data = plot_tree_ellipses, aes(x = Longitude, y = Latitude, group = Site_Plot), 
               fill = NA, color = "black", size = 1, alpha = 0.8) +
  # Trees as small black dots
  geom_point(data = trees_with_plots, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 1, stroke = 1, shape = 21, alpha = 0.9) +
  # Research plots as larger black dots
  geom_point(data = plots_data, aes(x = Longitude, y = Latitude), 
             color = "white", fill = "black", size = 3, stroke = 1, shape = 21, alpha = 0.7) +
  
  # Wetland label with lines to both wetland ellipses
  annotate("segment", 
           x = wetland_label_pos$x, y = wetland_label_pos$y,
           xend = wetland_ellipses$center_lon[1] - offset_x/6, 
           yend = wetland_ellipses$center_lat[1],
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("segment", 
           x = wetland_label_pos$x, y = wetland_label_pos$y,
           xend = wetland_ellipses$center_lon[2] - offset_x/6, 
           yend = wetland_ellipses$center_lat[2],
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", 
           x = wetland_label_pos$x, y = wetland_label_pos$y,
           label = "Wetland", color = "black", size = 4, fontface = "bold",
           hjust = 1, vjust = 0.5) +
  
  # Intermediate label with lines to both intermediate ellipses
  annotate("segment", 
           x = intermediate_label_pos$x, y = intermediate_label_pos$y,
           xend = intermediate_ellipses$center_lon[1], 
           yend = intermediate_ellipses$center_lat[1] + offset_y/6,
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("segment", 
           x = intermediate_label_pos$x, y = intermediate_label_pos$y,
           xend = intermediate_ellipses$center_lon[2], 
           yend = intermediate_ellipses$center_lat[2] + offset_y/6,
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", 
           x = intermediate_label_pos$x, y = intermediate_label_pos$y,
           label = "Intermediate", color = "black", size = 4, fontface = "bold",
           hjust = 0.5, vjust = 0) +
  
  # Upland label with lines to both upland ellipses
  annotate("segment", 
           x = upland_label_pos$x, y = upland_label_pos$y,
           xend = upland_ellipses$center_lon[1] + offset_x/6, 
           yend = upland_ellipses$center_lat[1],
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("segment", 
           x = upland_label_pos$x, y = upland_label_pos$y,
           xend = upland_ellipses$center_lon[2] + offset_x/6, 
           yend = upland_ellipses$center_lat[2],
           color = "black", size = 1.2, arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", 
           x = upland_label_pos$x, y = upland_label_pos$y,
           label = "Upland", color = "black", size = 4, fontface = "bold",
           hjust = 0, vjust = 0.5) +
  
  scale_fill_viridis_c(name = "Soil Moisture\n(VWC %)", option = "mako", direction = -1) +
  scale_alpha_identity() +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey80", linewidth = 0.2),
        panel.background = element_rect(fill = "grey90"))

print(p7)

# Save the annotated plot
ggsave("annotated_moisture_hillshade_ellipses.png", p7, width = 14, height = 8, dpi = 300)

