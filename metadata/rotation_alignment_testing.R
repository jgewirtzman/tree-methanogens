# LEAFLET MAP COMPARING DIFFERENT ROTATION METHODS
# This creates an interactive map where you can toggle between different coordinate transformations

library(leaflet)
library(dplyr)

cat("=== CREATING MULTI-ROTATION COMPARISON MAP ===\n")

# Use trees with PX/PY coordinates
trees_with_local <- fg_clean %>% filter(!is.na(PX) & !is.na(PY) & PX != 0 & PY != 0)

# Reference point and Earth radius
lat1_deg <- 41.989211
lon1_deg <- -72.13092
lat1_rad <- lat1_deg * pi/180
lon1_rad <- lon1_deg * pi/180
earth_radius <- 6371000

# Function for geodetic transformation (reusing from earlier)
geodetic_transform <- function(px, py, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  distance <- sqrt(px^2 + py^2)
  original_bearing <- atan2(px, py)
  rotated_bearing <- original_bearing + rotation_angle
  
  px_rotated <- distance * sin(rotated_bearing)
  py_rotated <- distance * cos(rotated_bearing)
  
  bearing <- atan2(px_rotated, py_rotated)
  distance_m <- sqrt(px_rotated^2 + py_rotated^2)
  
  lat2_rad <- asin(sin(lat1_rad) * cos(distance_m/earth_radius) + 
                     cos(lat1_rad) * sin(distance_m/earth_radius) * cos(bearing))
  
  lon2_rad <- lon1_rad + atan2(sin(bearing) * sin(distance_m/earth_radius) * cos(lat1_rad),
                               cos(distance_m/earth_radius) - sin(lat1_rad) * sin(lat2_rad))
  
  return(list(
    longitude = lon2_rad * 180/pi,
    latitude = lat2_rad * 180/pi
  ))
}

# Create datasets with different rotation angles
rotation_methods <- list(
  "Excel Method (10°)" = 10 * pi/180,
  "R Correlation (-8°)" = -8 * pi/180,
  "Corner Method (-10°)" = -10 * pi/180,
  "No Rotation (0°)" = 0,
  "PCA Method" = -135.37 * pi/180  # From your earlier diagnostics
)

# Apply each transformation
transformed_datasets <- list()

for(method_name in names(rotation_methods)) {
  rotation_angle <- rotation_methods[[method_name]]
  
  cat("Applying", method_name, "with rotation", rotation_angle * 180/pi, "degrees\n")
  
  result <- geodetic_transform(
    trees_with_local$PX, 
    trees_with_local$PY,
    lat1_rad, lon1_rad, rotation_angle, earth_radius
  )
  
  transformed_datasets[[method_name]] <- trees_with_local %>%
    mutate(
      longitude = result$longitude,
      latitude = result$latitude,
      method = method_name
    )
}

# Also add original GPS coordinates for comparison
if(sum(!is.na(trees_with_local$Longitude)) > 0) {
  transformed_datasets[["Original GPS"]] <- trees_with_local %>%
    filter(!is.na(Longitude) & !is.na(Latitude)) %>%
    mutate(
      longitude = Longitude,
      latitude = Latitude,
      method = "Original GPS"
    )
}

cat("Created", length(transformed_datasets), "coordinate datasets\n\n")

# Create color palette for different methods
method_colors <- c(
  "Excel Method (10°)" = "red",
  "R Correlation (-8°)" = "blue", 
  "Corner Method (-10°)" = "green",
  "No Rotation (0°)" = "purple",
  "PCA Method" = "orange",
  "Original GPS" = "black"
)

# Calculate overall bounding box for all methods
all_coords <- do.call(rbind, transformed_datasets)
center_lon <- mean(all_coords$longitude, na.rm = TRUE)
center_lat <- mean(all_coords$latitude, na.rm = TRUE)

cat("Map center:", round(center_lat, 6), "°N,", round(center_lon, 6), "°W\n")

# Create the interactive comparison map
comparison_map <- leaflet() %>%
  # Add satellite imagery
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addProviderTiles(providers$OpenStreetMap, group = "Street Map") %>%
  addProviderTiles(providers$Esri.WorldTopoMap, group = "Topographic")

# Add each transformation method as a separate layer
for(method_name in names(transformed_datasets)) {
  data <- transformed_datasets[[method_name]]
  color <- method_colors[[method_name]]
  
  comparison_map <- comparison_map %>%
    addCircleMarkers(
      data = data,
      lng = ~longitude,
      lat = ~latitude,
      radius = ~ifelse(!is.na(BasalArea_m2), sqrt(BasalArea_m2) * 8, 3),
      color = color,
      fillColor = color,
      fillOpacity = 0.6,
      stroke = TRUE,
      weight = 1,
      group = method_name,  # This creates the layer group
      popup = ~paste(
        "Method:", method, "<br/>",
        "Tag:", Tag, "<br/>",
        if("DBH" %in% names(.)) paste("DBH:", round(DBH, 1), "cm<br/>") else "",
        "Coordinates:", round(longitude, 6), ",", round(latitude, 6)
      )
    )
}

# Add layer controls to toggle between methods
comparison_map <- comparison_map %>%
  addLayersControl(
    baseGroups = c("Satellite", "Street Map", "Topographic"),
    overlayGroups = names(transformed_datasets),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  
  # Set the view
  setView(
    lng = center_lon,
    lat = center_lat,
    zoom = 17
  ) %>%
  
  # Add a legend
  addLegend(
    colors = method_colors[names(transformed_datasets)],
    labels = names(transformed_datasets),
    title = "Rotation Methods",
    position = "bottomright",
    opacity = 0.8
  )

# Display the map
print(comparison_map)

cat("\n=== INTERACTIVE COMPARISON MAP CREATED ===\n")
cat("Use the layer controls on the top-right to toggle different rotation methods on/off\n")
cat("Switch between satellite/street/topographic basemaps\n")
cat("Click on individual trees to see coordinates and method used\n")
cat("Zoom in to see how different methods align with forest features\n\n")

# Create a summary table of coordinate ranges for each method
cat("=== COORDINATE RANGES BY METHOD ===\n")
for(method_name in names(transformed_datasets)) {
  data <- transformed_datasets[[method_name]]
  lon_range <- range(data$longitude, na.rm = TRUE)
  lat_range <- range(data$latitude, na.rm = TRUE)
  
  cat(sprintf("%-20s: Lon %.6f to %.6f, Lat %.6f to %.6f\n", 
              method_name, lon_range[1], lon_range[2], lat_range[1], lat_range[2]))
}

cat("\nLook for the method that best aligns trees with visible forest boundaries in satellite imagery\n")