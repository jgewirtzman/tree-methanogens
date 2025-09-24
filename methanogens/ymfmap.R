# R script to create a map of US forest types with Yale-Myers Forest location
# Using the USFS Forest Type Groups data available from Google Earth Engine

# Load required packages
library(rgee)       # For Google Earth Engine integration
library(sf)         # For spatial data handling
library(raster)     # For raster data processing
library(ggplot2)    # For mapping
library(dplyr)      # For data manipulation
library(ggspatial)  # For north arrow and scale bar
library(RColorBrewer) # For color palettes
library(maps)       # For base maps

# Initialize Earth Engine (only needs to be done once per session)
ee_Initialize()

# Define Yale-Myers Forest location (northeastern Connecticut)
yale_myers <- data.frame(
  name = "Yale-Myers Forest",
  long = -72.1245,  # Longitude
  lat = 41.9359     # Latitude
)

# Load the forest type group data from GEE
forest_group <- ee$ImageCollection("projects/sat-io/open-datasets/USFS/national-forest-group")
forest_group_image <- forest_group$mosaic()

# Define the area of interest - focus on Connecticut and surrounding region
states <- ee$FeatureCollection('TIGER/2018/States')
CT <- states$filter(ee$Filter$eq('NAME', 'Connecticut'))
export_region <- CT$geometry()$buffer(50000)

# Define the bounding box for our map in lon/lat (WGS84)
bbox <- c(-73.7, 41.0, -71.8, 42.1)  # xmin, ymin, xmax, ymax

# Export the forest type data for our region of interest
forest_data <- ee_as_raster(
  image = forest_group_image, 
  region = export_region$bounds(),
  scale = 250,  # 250m resolution is sufficient for a regional overview
  via = "drive"  # Export via Google Drive
)

# Convert the raster to a data frame for ggplot
forest_df <- as.data.frame(forest_data, xy = TRUE)
names(forest_df) <- c("x", "y", "forest_group")

# Remove NA values (areas without forest data)
forest_df <- forest_df %>% filter(!is.na(forest_group))

# Define the forest type groups and codes
forest_codes <- data.frame(
  code = c(0, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 
           360, 370, 400, 500, 600, 700, 800, 900, 920, 940, 950, 999),
  group = c("Non-forest", "White/Red/Jack Pine", "Spruce/Fir", "Longleaf/Slash Pine", 
            "Loblolly/Shortleaf Pine", "Other Eastern Softwoods", "Pinyon/Juniper", 
            "Douglas-fir", "Ponderosa Pine", "Western White Pine", 
            "Fir/Spruce/Mountain Hemlock", "Lodgepole Pine", "Hemlock/Sitka Spruce", 
            "Other Western Softwoods", "California Mixed Conifer", "Exotic Softwoods", 
            "Oak/Pine", "Oak/Hickory", "Oak/Gum/Cypress", "Elm/Ash/Cottonwood", 
            "Maple/Beech/Birch", "Aspen/Birch", "Other Hardwoods", 
            "Other Western Hardwoods", "Exotic Hardwoods", "Nonstocked")
)

# Join the forest code descriptions to our data
forest_df$code <- forest_df$forest_group
forest_df <- left_join(forest_df, forest_codes, by = c("code" = "code"))

# Define colors for forest types found in the northeastern US
# Using a subset of colors that correspond to the forest types in our region
ne_forest_colors <- c(
  "0" = "#FFFFFF",     # Non-forest - White
  "100" = "#005E00",   # White/Red/Jack Pine - Dark Green
  "120" = "#008600",   # Spruce/Fir - Bright Green
  "400" = "#CF5E00",   # Oak/Pine - Orange-Brown
  "500" = "#A86400",   # Oak/Hickory - Medium Brown
  "700" = "#9AC0D8",   # Elm/Ash/Cottonwood - Light Blue
  "800" = "#FF8CA6",   # Maple/Beech/Birch - Pink
  "900" = "#FFBCD8"    # Aspen/Birch - Light Pink
)

# Get state boundaries for New England
states_map <- map_data("state")
ne_states <- subset(states_map, region %in% c("connecticut", "massachusetts", 
                                              "rhode island", "new hampshire", 
                                              "vermont", "maine", "new york"))

# Get county boundaries for CT
counties_map <- map_data("county")
ct_counties <- subset(counties_map, region == "connecticut")

# Create the main map
main_map <- ggplot() +
  # Add the forest type raster
  geom_raster(data = forest_df, 
              aes(x = x, y = y, fill = as.factor(code)), alpha = 0.8) +
  # Custom color scale for forest types
  scale_fill_manual(
    values = ne_forest_colors,
    name = "Forest Type Group",
    labels = function(x) {
      sapply(x, function(code) {
        group_name <- forest_codes$group[forest_codes$code == as.numeric(code)]
        return(paste0(code, ": ", group_name))
      })
    },
    drop = FALSE
  ) +
  # Add state boundaries
  geom_polygon(data = ne_states, 
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray40", size = 0.2) +
  # Add county boundaries for CT
  geom_polygon(data = ct_counties, 
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray60", size = 0.1) +
  # Highlight Connecticut with slightly darker outline
  geom_polygon(data = subset(states_map, region == "connecticut"),
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.5) +
  # Add Yale-Myers Forest point
  geom_point(data = yale_myers, 
             aes(x = long, y = lat),
             color = "red", size = 3) +
  # Add label for Yale-Myers Forest
  geom_text(data = yale_myers, 
            aes(x = long, y = lat + 0.05, label = name),
            size = 4, fontface = "bold") +
  # Set the map extent to our bbox
  coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
  # Add title and caption
  labs(title = "Location of Yale-Myers Forest",
       subtitle = "Forest Type Groups in Northeastern Connecticut",
       caption = "Data Source: USFS National Forest Type Groups via Google Earth Engine") +
  # Add north arrow and scale bar
  annotation_north_arrow(
    location = "bl", 
    which_north = "true",
    pad_x = unit(0.1, "in"), 
    pad_y = unit(0.1, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    bar_cols = c("black", "white"),
    text_family = "sans",
    text_cex = 0.7,
    width_hint = 0.3
  ) +
  # Clean theme
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE))

# Create inset map showing Yale-Myers Forest in the context of New England
inset_map <- ggplot() +
  # Add New England states
  geom_polygon(data = ne_states, 
               aes(x = long, y = lat, group = group),
               fill = "white", color = "gray40", size = 0.2) +
  # Highlight Connecticut
  geom_polygon(data = subset(states_map, region == "connecticut"),
               aes(x = long, y = lat, group = group),
               fill = "lightgreen", color = "black", size = 0.5, alpha = 0.3) +
  # Add Yale-Myers Forest location
  geom_point(data = yale_myers, 
             aes(x = long, y = lat),
             color = "red", size = 2) +
  # Set the map extent to show all of New England
  coord_fixed(xlim = c(-73.7, -66.9), ylim = c(41, 47.5), ratio = 1.3) +
  # Clean theme for inset
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA, alpha = 0.7)
  )

# Combine the main map and inset map
# Calculate the dimensions for the inset placement
inset_width <- 0.25  # Width of inset as proportion of main plot
inset_height <- 0.25 # Height of inset as proportion of main plot
inset_x <- 0.02      # X position of inset from left edge
inset_y <- 0.02      # Y position of inset from bottom edge

# Plot with inset
final_map <- main_map +
  annotation_custom(
    grob = ggplotGrob(inset_map),
    xmin = bbox[1] + inset_x * (bbox[3] - bbox[1]),
    xmax = bbox[1] + (inset_x + inset_width) * (bbox[3] - bbox[1]),
    ymin = bbox[2] + inset_y * (bbox[4] - bbox[2]),
    ymax = bbox[2] + (inset_y + inset_height) * (bbox[4] - bbox[2])
  )

# Display the final map
print(final_map)

# Save the map
ggsave("yale_myers_forest_map.png", final_map, width = 10, height = 7, dpi = 300)

# Alternative approach without rgee (if you've already downloaded the data)
# If you've already downloaded the forest type raster from GEE, you can use this code:

# forest_data <- raster("path/to/downloaded/forest_type.tif")
# forest_df <- as.data.frame(forest_data, xy = TRUE)
# names(forest_df) <- c("x", "y", "forest_group")
# forest_df <- forest_df %>% filter(!is.na(forest_group))
# Then continue with the rest of the code above starting from the forest codes definition