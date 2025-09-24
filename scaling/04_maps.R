# =============================================================================
# TREE CH4 FLUX MAP VISUALIZATION
# Creates maps of trees colored by predicted flux rates
# Run AFTER completing the RF workflow analysis
# =============================================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(scales)
library(sf)
library(tidyr)

cat("=== TREE CH4 FLUX MAP VISUALIZATION ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD DATA AND MODELS
# =============================================================================

cat("Loading RF models and data...\n")

# Load the saved models and data
load("RF_MODELS.RData")  # Contains TreeRF and SoilRF
load("rf_workflow_input_data_clean.RData")  # Contains all input data
load("TAXONOMY_PRIORS.RData")  # Contains taxonomy priors

# Extract needed components
INVENTORY <- rf_workflow_data$PLACEHOLDER_INVENTORY
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
TAXONOMY <- rf_workflow_data$PLACEHOLDER_TAXONOMY
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY

# Load SI tables and moisture calibration
SI_TABLES <- read.csv("SI_TABLES.csv")
MOISTURE_AFFINE_TABLE <- read.csv("MOISTURE_AFFINE_TABLE.csv")

cat("✓ Data loaded\n")
cat("  Trees in inventory:", nrow(INVENTORY), "\n\n")

# =============================================================================
# PREDICT FLUX FOR A SPECIFIC MONTH (or annual average)
# =============================================================================

# Choose which month to visualize (1-12, or 0 for annual average)
MONTH_TO_MAP <- 7  # July by default - change as needed

cat("Predicting fluxes for visualization...\n")

if(MONTH_TO_MAP == 0) {
  cat("  Calculating annual average fluxes...\n")
} else {
  cat("  Predicting for month", MONTH_TO_MAP, "\n")
}

# Function to predict tree flux for a given month
predict_tree_flux_month <- function(inventory_df, month_val, tree_rf, drivers, si_tables, 
                                    moisture_affine, moisture_fn, taxonomy_priors) {
  
  # Get monthly drivers
  air_temp_t <- drivers$air_temp_C_mean[month_val]
  soil_temp_t <- drivers$soil_temp_C_mean[month_val]
  if(is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  # Get seasonal index
  si_tree_t <- si_tables %>% 
    filter(group == "tree", month == month_val) %>% 
    pull(SI)
  if(length(si_tree_t) == 0) si_tree_t <- 0
  
  # Prepare predictions dataframe
  predictions <- inventory_df %>%
    mutate(
      # Moisture at tree location
      moisture_raw = moisture_fn(x, y),
      soil_moisture_abs = moisture_affine$alpha_t[month_val] + 
        moisture_affine$beta_t[month_val] * moisture_raw,
      soil_moisture_abs = pmax(0, pmin(0.75, soil_moisture_abs)),
      
      # Environmental variables
      air_temp_C = air_temp_t,
      soil_temp_C = soil_temp_t,
      SI_tree = si_tree_t,
      
      # Month as cyclic features
      month_sin = sin(2 * pi * month_val / 12),
      month_cos = cos(2 * pi * month_val / 12),
      
      # Chamber type (default to rigid for inventory)
      chamber_rigid = 1,
      chamber_semirigid = 0,
      
      # Interactions
      moisture_x_airT = soil_moisture_abs * air_temp_C,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      
      # Taxonomy
      genus = sub(" .*", "", species),
      family = case_when(
        genus == "Acer" ~ "Sapindaceae",
        genus %in% c("Betula") ~ "Betulaceae",
        genus %in% c("Fagus", "Quercus") ~ "Fagaceae",
        genus == "Fraxinus" ~ "Oleaceae",
        genus %in% c("Pinus", "Tsuga") ~ "Pinaceae",
        genus == "Carya" ~ "Juglandaceae",
        TRUE ~ "OTHER"
      )
    )
  
  # Add taxonomy priors
  predictions$taxon_prior_asinh <- 0
  
  # Extract priors as named vectors
  species_priors <- setNames(taxonomy_priors$species$med_resid, 
                             taxonomy_priors$species$species)
  genus_priors <- setNames(taxonomy_priors$genus$med_resid, 
                           taxonomy_priors$genus$genus)
  family_priors <- setNames(taxonomy_priors$family$med_resid, 
                            taxonomy_priors$family$family)
  
  for(i in 1:nrow(predictions)) {
    sp <- predictions$species[i]
    gn <- predictions$genus[i]
    fm <- predictions$family[i]
    
    if(sp %in% names(species_priors)) {
      predictions$taxon_prior_asinh[i] <- species_priors[sp]
    } else if(gn %in% names(genus_priors)) {
      predictions$taxon_prior_asinh[i] <- genus_priors[gn]
    } else if(fm %in% names(family_priors)) {
      predictions$taxon_prior_asinh[i] <- family_priors[fm]
    }
  }
  
  # Build feature matrix
  numeric_features <- predictions %>%
    select(dbh_m, air_temp_C, soil_temp_C, soil_moisture_abs, SI_tree, 
           moisture_x_airT, moisture_x_soilT, taxon_prior_asinh,
           chamber_rigid, chamber_semirigid, month_sin, month_cos) %>%
    as.matrix()
  
  # Add categorical dummies
  species_dummies <- model.matrix(~ species - 1, data = predictions)
  genus_dummies <- model.matrix(~ genus - 1, data = predictions)
  family_dummies <- model.matrix(~ family - 1, data = predictions)
  
  X_pred <- cbind(numeric_features, species_dummies, genus_dummies, family_dummies)
  
  # Align with training features
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred), 
                           ncol = length(tree_rf$forest$independent.variable.names))
  colnames(X_pred_aligned) <- tree_rf$forest$independent.variable.names
  
  for(col_name in colnames(X_pred_aligned)) {
    if(col_name %in% colnames(X_pred)) {
      X_pred_aligned[, col_name] <- X_pred[, col_name]
    }
  }
  
  # Predict
  pred_asinh <- predict(tree_rf, X_pred_aligned)$predictions
  pred_flux <- sinh(pred_asinh)
  
  return(pred_flux)
}

# Predict fluxes
if(MONTH_TO_MAP == 0) {
  # Calculate annual average
  all_months_flux <- matrix(NA, nrow = nrow(INVENTORY), ncol = 12)
  
  for(m in 1:12) {
    cat("    Month", m, "...")
    all_months_flux[, m] <- predict_tree_flux_month(
      INVENTORY, m, TreeRF, DRIVERS, SI_TABLES, 
      MOISTURE_AFFINE_TABLE, moisture_lookup_xy, TAXONOMY_PRIORS
    )
    cat(" done\n")
  }
  
  INVENTORY$flux_umol_m2_s <- rowMeans(all_months_flux, na.rm = TRUE)
  map_title <- "Annual Average Tree CH4 Flux"
  
} else {
  # Single month
  INVENTORY$flux_umol_m2_s <- predict_tree_flux_month(
    INVENTORY, MONTH_TO_MAP, TreeRF, DRIVERS, SI_TABLES, 
    MOISTURE_AFFINE_TABLE, moisture_lookup_xy, TAXONOMY_PRIORS
  )
  
  month_names <- c("January", "February", "March", "April", "May", "June",
                   "July", "August", "September", "October", "November", "December")
  map_title <- paste(month_names[MONTH_TO_MAP], "Tree CH4 Flux")
}

# Add flux categories
INVENTORY <- INVENTORY %>%
  mutate(
    flux_category = case_when(
      flux_umol_m2_s < 0 ~ "Sink (<0)",
      flux_umol_m2_s < 0.001 ~ "Very Low\n(0-0.001)",
      flux_umol_m2_s < 0.005 ~ "Low\n(0.001-0.005)",
      flux_umol_m2_s < 0.01 ~ "Moderate\n(0.005-0.01)",
      flux_umol_m2_s < 0.02 ~ "High\n(0.01-0.02)",
      TRUE ~ "Very High\n(>0.02)"
    ),
    flux_category = factor(flux_category, 
                           levels = c("Sink (<0)", "Very Low\n(0-0.001)", 
                                      "Low\n(0.001-0.005)", "Moderate\n(0.005-0.01)",
                                      "High\n(0.01-0.02)", "Very High\n(>0.02)"))
  )

cat("✓ Flux predictions complete\n")
cat("  Flux range:", round(min(INVENTORY$flux_umol_m2_s), 6), "to", 
    round(max(INVENTORY$flux_umol_m2_s), 6), "μmol m-2 s-1\n\n")

# =============================================================================
# CREATE MAPS
# =============================================================================

cat("Creating maps...\n")

# Map 1: Continuous color scale
p1 <- ggplot(INVENTORY, aes(x = x, y = y)) +  # x = longitude, y = latitude
  geom_point(aes(color = flux_umol_m2_s, size = dbh_m * 100), alpha = 0.7) +
  scale_color_viridis(
    name = "CH4 Flux\n(μmol m⁻² s⁻¹)",
    option = "turbo",
    trans = "pseudo_log",
    breaks = c(-0.001, 0, 0.001, 0.005, 0.01, 0.02),
    labels = c("-0.001", "0", "0.001", "0.005", "0.01", "0.02")
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.5, 4),
    breaks = c(5, 10, 20, 40, 60),
    limits = c(0, 80)
  ) +
  coord_equal() +
  labs(
    title = paste(map_title, "- Continuous Scale"),
    subtitle = paste("N =", nrow(INVENTORY), "trees | GPS coordinates"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    panel.grid = element_line(color = "grey90", size = 0.5),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("tree_flux_map_continuous.png", p1, width = 12, height = 9, dpi = 300)

# Map 2: Categorical color scale
p2 <- ggplot(INVENTORY, aes(x = x, y = y)) +
  geom_point(aes(color = flux_category, size = dbh_m * 100), alpha = 0.8) +
  scale_color_manual(
    name = "CH4 Flux Category\n(μmol m⁻² s⁻¹)",
    values = c("Sink (<0)" = "blue",
               "Very Low\n(0-0.001)" = "lightgreen",
               "Low\n(0.001-0.005)" = "yellow",
               "Moderate\n(0.005-0.01)" = "orange",
               "High\n(0.01-0.02)" = "darkorange",
               "Very High\n(>0.02)" = "red")
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.5, 4),
    breaks = c(5, 10, 20, 40, 60),
    limits = c(0, 80)
  ) +
  coord_equal() +
  labs(
    title = paste(map_title, "- Categorical Scale"),
    subtitle = paste("N =", nrow(INVENTORY), "trees | GPS coordinates"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    panel.grid = element_line(color = "grey90", size = 0.5),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("tree_flux_map_categorical.png", p2, width = 12, height = 9, dpi = 300)

# Map 3: Species-specific flux
# Calculate mean flux by species
species_flux <- INVENTORY %>%
  group_by(species) %>%
  summarise(
    mean_flux = mean(flux_umol_m2_s),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 5) %>%  # Only species with 5+ individuals
  arrange(desc(mean_flux))

# Keep only top emitting species
top_species <- head(species_flux$species, 10)

p3 <- INVENTORY %>%
  filter(species %in% top_species) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(data = INVENTORY, color = "grey80", size = 0.5, alpha = 0.3) +  # Background
  geom_point(aes(color = flux_umol_m2_s, size = dbh_m * 100), alpha = 0.8) +
  facet_wrap(~ species, ncol = 3) +
  scale_color_viridis(
    name = "CH4 Flux\n(μmol m⁻² s⁻¹)",
    option = "plasma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.5, 3),
    breaks = c(10, 20, 40)
  ) +
  coord_equal() +
  labs(
    title = "Tree CH4 Flux by Species",
    subtitle = "Top 10 species by abundance (≥5 individuals)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 10)
  )

ggsave("tree_flux_map_by_species.png", p3, width = 14, height = 12, dpi = 300)

# Map 4: Spatial hotspot analysis (hexbin)
library(hexbin)

# Create hexbin plot
p4 <- ggplot(INVENTORY, aes(x = x, y = y)) +
  stat_summary_hex(aes(z = flux_umol_m2_s), fun = mean, bins = 30, alpha = 0.8) +
  geom_point(data = INVENTORY %>% filter(flux_umol_m2_s > quantile(flux_umol_m2_s, 0.95)),
             aes(size = dbh_m * 100), color = "red", alpha = 0.7) +
  scale_fill_viridis(
    name = "Mean CH4 Flux\n(μmol m⁻² s⁻¹)",
    option = "magma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(
    name = "DBH (cm)\nTop 5% emitters",
    range = c(1, 4),
    breaks = c(10, 20, 40)
  ) +
  coord_equal() +
  labs(
    title = "CH4 Flux Spatial Hotspots",
    subtitle = "Hexagonal bins show mean flux; red points = top 5% emitters",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid = element_line(color = "grey90", size = 0.5)
  )

ggsave("tree_flux_map_hotspots.png", p4, width = 12, height = 9, dpi = 300)

cat("✓ Maps saved\n\n")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("=== FLUX STATISTICS ===\n")

# Overall statistics
cat("\nOverall:\n")
cat("  Mean flux:", round(mean(INVENTORY$flux_umol_m2_s), 6), "μmol m⁻² s⁻¹\n")
cat("  Median flux:", round(median(INVENTORY$flux_umol_m2_s), 6), "μmol m⁻² s⁻¹\n")
cat("  SD:", round(sd(INVENTORY$flux_umol_m2_s), 6), "μmol m⁻² s⁻¹\n")
cat("  Range:", round(min(INVENTORY$flux_umol_m2_s), 6), "to", 
    round(max(INVENTORY$flux_umol_m2_s), 6), "μmol m⁻² s⁻¹\n")
cat("  Trees as CH4 sinks:", sum(INVENTORY$flux_umol_m2_s < 0), 
    "(", round(100*mean(INVENTORY$flux_umol_m2_s < 0), 1), "%)\n")

# By species
cat("\nTop 10 emitting species (mean flux):\n")
species_summary <- INVENTORY %>%
  group_by(species) %>%
  summarise(
    n = n(),
    mean_flux = mean(flux_umol_m2_s),
    median_flux = median(flux_umol_m2_s),
    sd_flux = sd(flux_umol_m2_s),
    .groups = "drop"
  ) %>%
  filter(n >= 3) %>%
  arrange(desc(mean_flux))

print(head(species_summary, 10))

# By size class
cat("\nFlux by size class:\n")
size_summary <- INVENTORY %>%
  mutate(
    size_class = cut(dbh_m * 100, 
                     breaks = c(0, 10, 20, 40, 100),
                     labels = c("Small (<10cm)", "Medium (10-20cm)", 
                                "Large (20-40cm)", "Very Large (>40cm)"))
  ) %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    mean_flux = mean(flux_umol_m2_s),
    median_flux = median(flux_umol_m2_s),
    .groups = "drop"
  )

print(size_summary)

# Save results
write.csv(INVENTORY %>% select(tree_id, species, dbh_m, x, y, flux_umol_m2_s, flux_category),
          paste0("tree_flux_predictions_month", MONTH_TO_MAP, ".csv"), 
          row.names = FALSE)

cat("\n=== VISUALIZATION COMPLETE ===\n")
cat("Generated files:\n")
cat("  - tree_flux_map_continuous.png\n")
cat("  - tree_flux_map_categorical.png\n")
cat("  - tree_flux_map_by_species.png\n")
cat("  - tree_flux_map_hotspots.png\n")
cat("  - tree_flux_predictions_month", MONTH_TO_MAP, ".csv\n")




















# =============================================================================
# SPATIAL MAPS OF PREDICTED FLUXES WITH GEODETIC TRANSFORMATION
# Run this after completing the main RF workflow
# Properly transforms PX/PY coordinates to GPS using geodetic transformation
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)
library(scales)

cat("\n=== CREATING SPATIAL FLUX MAPS WITH GEODETIC TRANSFORMATION ===\n")

# =============================================================================
# GEODETIC TRANSFORMATION SETUP
# =============================================================================

# Parameters from your workflow
ref_lat <- 41.989211
ref_lon <- -72.13092
lat1_rad <- ref_lat * pi/180
lon1_rad <- ref_lon * pi/180
rotation_angle <- 0.174533  # 10 degrees in radians
earth_radius <- 6371000

# Geodetic transformation function (exact from your workflow)
geodetic_transform <- function(px, py, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  
  # Step 1: Apply rotation to PX/PY coordinates
  distance <- sqrt(px^2 + py^2)
  original_bearing <- atan2(px, py)
  rotated_bearing <- original_bearing + rotation_angle
  
  px_rotated <- distance * sin(rotated_bearing)
  py_rotated <- distance * cos(rotated_bearing)
  
  # Step 2: Convert to GPS using geodetic formulas
  bearing <- atan2(px_rotated, py_rotated)
  distance_m <- sqrt(px_rotated^2 + py_rotated^2)
  
  # Geodetic forward calculation
  lat2_rad <- asin(sin(lat1_rad) * cos(distance_m/earth_radius) + 
                     cos(lat1_rad) * sin(distance_m/earth_radius) * cos(bearing))
  
  lon2_rad <- lon1_rad + atan2(sin(bearing) * sin(distance_m/earth_radius) * cos(lat1_rad),
                               cos(distance_m/earth_radius) - sin(lat1_rad) * sin(lat2_rad))
  
  return(list(
    longitude = lon2_rad * 180/pi,
    latitude = lat2_rad * 180/pi
  ))
}

# Inverse geodetic transformation
geodetic_inverse <- function(lat, lon, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  lat_rad <- lat * pi/180
  lon_rad <- lon * pi/180
  
  # Calculate bearing and distance from reference point
  delta_lon <- lon_rad - lon1_rad
  
  y_comp <- sin(delta_lon) * cos(lat_rad)
  x_comp <- cos(lat1_rad) * sin(lat_rad) - sin(lat1_rad) * cos(lat_rad) * cos(delta_lon)
  bearing <- atan2(y_comp, x_comp)
  
  # Calculate distance
  a <- sin((lat_rad - lat1_rad)/2)^2 + cos(lat1_rad) * cos(lat_rad) * sin(delta_lon/2)^2
  c <- 2 * asin(sqrt(a))
  distance_m <- earth_radius * c
  
  # Convert to local coordinates with rotation
  bearing_adjusted <- bearing - rotation_angle
  x <- distance_m * sin(bearing_adjusted)
  y <- distance_m * cos(bearing_adjusted)
  
  return(list(x = x, y = y))
}

# =============================================================================
# TRANSFORM INVENTORY COORDINATES
# =============================================================================

cat("Transforming inventory coordinates...\n")

# Check what columns we have
cat("Available columns in INVENTORY:", names(INVENTORY), "\n")

# Check if we need to transform from x,y to GPS or if GPS is already correct
if("x" %in% names(INVENTORY) && "y" %in% names(INVENTORY)) {
  cat("  Found local x,y coordinates\n")
  
  # Check if GPS coordinates already exist and might be correct
  if("longitude" %in% names(INVENTORY) && "latitude" %in% names(INVENTORY)) {
    # Check if they look like they've already been transformed
    lon_range <- range(INVENTORY$longitude, na.rm = TRUE)
    lat_range <- range(INVENTORY$latitude, na.rm = TRUE)
    
    if(lon_range[1] > -180 && lon_range[2] < 180 && 
       lat_range[1] > -90 && lat_range[2] < 90 &&
       lon_range[1] > -73 && lon_range[2] < -72 &&  # Near expected area
       lat_range[1] > 41.5 && lat_range[2] < 42.5) {
      cat("  GPS coordinates appear to be already corrected\n")
      cat("  Longitude range:", lon_range, "\n")
      cat("  Latitude range:", lat_range, "\n")
      # Keep existing GPS coordinates
    } else {
      cat("  GPS coordinates don't look right - applying geodetic transformation\n")
      # Transform x,y to GPS
      trees_with_xy <- !is.na(INVENTORY$x) & !is.na(INVENTORY$y)
      
      if(sum(trees_with_xy) > 0) {
        cat("  Transforming", sum(trees_with_xy), "trees from x,y to GPS\n")
        
        # Transform x,y to GPS using geodetic transform
        result <- geodetic_transform(
          INVENTORY$x[trees_with_xy], 
          INVENTORY$y[trees_with_xy],
          lat1_rad, lon1_rad, rotation_angle, earth_radius
        )
        
        # Store transformed coordinates
        INVENTORY$longitude[trees_with_xy] <- result$longitude
        INVENTORY$latitude[trees_with_xy] <- result$latitude
      }
    }
  } else {
    # No GPS coordinates exist, create them from x,y
    cat("  No GPS coordinates found - transforming from x,y\n")
    trees_with_xy <- !is.na(INVENTORY$x) & !is.na(INVENTORY$y)
    
    if(sum(trees_with_xy) > 0) {
      cat("  Transforming", sum(trees_with_xy), "trees from x,y to GPS\n")
      
      # Transform x,y to GPS
      result <- geodetic_transform(
        INVENTORY$x[trees_with_xy], 
        INVENTORY$y[trees_with_xy],
        lat1_rad, lon1_rad, rotation_angle, earth_radius
      )
      
      # Store transformed coordinates
      INVENTORY$longitude <- NA
      INVENTORY$latitude <- NA
      INVENTORY$longitude[trees_with_xy] <- result$longitude
      INVENTORY$latitude[trees_with_xy] <- result$latitude
    }
  }
} else if("longitude" %in% names(INVENTORY) && "latitude" %in% names(INVENTORY)) {
  # Only GPS coordinates exist, use them
  cat("  Using existing longitude and latitude columns\n")
  # GPS coordinates already in place, nothing to do
  
  # Create local x,y for moisture lookup
  cat("  Creating local x,y coordinates for moisture lookup\n")
  xy_result <- geodetic_inverse(
    INVENTORY$latitude, INVENTORY$longitude,
    lat1_rad, lon1_rad, rotation_angle, earth_radius
  )
  INVENTORY$x <- xy_result$x
  INVENTORY$y <- xy_result$y
} else {
  stop("No usable coordinate columns found in INVENTORY!")
}

# Report coordinate ranges
cat("\nFinal GPS coordinate ranges:\n")
cat("  Longitude:", range(INVENTORY$longitude, na.rm=TRUE), "\n")
cat("  Latitude:", range(INVENTORY$latitude, na.rm=TRUE), "\n")
cat("  Local x:", range(INVENTORY$x, na.rm=TRUE), "\n")
cat("  Local y:", range(INVENTORY$y, na.rm=TRUE), "\n\n")

# =============================================================================
# FUNCTION TO CREATE FLUX MAP FOR A SPECIFIC MONTH
# =============================================================================

create_monthly_flux_maps <- function(month_num, 
                                     inventory_df = INVENTORY,
                                     soil_grid_res = 50,
                                     plot_area = PLOT_AREA) {
  
  cat("\nGenerating flux maps for month", month_num, "(", month.abb[month_num], ")...\n")
  
  # Get monthly drivers for this month
  air_temp_t <- DRIVERS$air_temp_C_mean[month_num]
  soil_temp_t <- DRIVERS$soil_temp_C_mean[month_num]
  if (is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  # =============================================================================
  # 1. TREE FLUX PREDICTIONS
  # =============================================================================
  
  cat("  Predicting tree fluxes...\n")
  
  # Get seasonal index for trees
  SI_tree_t <- SI_TABLES %>% 
    filter(group == "tree", month == month_num) %>% 
    pull(SI) %>% 
    {if(length(.) > 0) .[1] else 0}
  
  # Predict for each tree in inventory (using LOCAL x,y for moisture lookup)
  tree_predictions <- inventory_df %>%
    mutate(
      # Environmental conditions - using LOCAL x,y for moisture lookup
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = clip_to_range(
        MOISTURE_AFFINE_TABLE$alpha_t[month_num] + 
          MOISTURE_AFFINE_TABLE$beta_t[month_num] * moisture_raw,
        training_ranges$moisture
      ),
      air_temp_C = clip_to_range(air_temp_t, training_ranges$air_temp),
      soil_temp_C = clip_to_range(soil_temp_t, training_ranges$soil_temp),
      dbh_m_clipped = clip_to_range(dbh_m, training_ranges$dbh),
      
      # Features
      SI_tree = SI_tree_t,
      chamber_rigid = 1,
      chamber_semirigid = 0,
      moisture_x_airT = soil_moisture_abs * air_temp_C,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      
      # Taxonomy
      genus = sub(" .*", "", species),
      family = case_when(
        genus == "Acer" ~ "Sapindaceae",
        genus %in% c("Betula") ~ "Betulaceae", 
        genus %in% c("Fagus", "Quercus") ~ "Fagaceae",
        genus == "Fraxinus" ~ "Oleaceae",
        genus %in% c("Pinus", "Tsuga") ~ "Pinaceae",
        genus == "Carya" ~ "Juglandaceae",
        genus == "Kalmia" ~ "Ericaceae",
        TRUE ~ "OTHER"
      )
    )
  
  # Assign taxon priors
  tree_predictions$taxon_prior_asinh <- 0
  for (i in 1:nrow(tree_predictions)) {
    sp <- tree_predictions$species[i]
    if (exists("species_priors") && sp %in% names(species_priors)) {
      tree_predictions$taxon_prior_asinh[i] <- species_priors[sp]
    } else if (exists("genus_priors")) {
      gn <- tree_predictions$genus[i]
      if (gn %in% names(genus_priors)) {
        tree_predictions$taxon_prior_asinh[i] <- genus_priors[gn]
      }
    }
  }
  
  # Build feature matrix for prediction
  numeric_features <- tree_predictions %>%
    select(dbh_m_clipped, air_temp_C, soil_temp_C, 
           soil_moisture_abs, SI_tree, moisture_x_airT, 
           moisture_x_soilT, taxon_prior_asinh,
           chamber_rigid, chamber_semirigid) %>%
    as.matrix()
  
  # Add categorical features
  species_dummies <- model.matrix(~ species - 1, data = tree_predictions)
  genus_dummies <- model.matrix(~ genus - 1, data = tree_predictions)
  family_dummies <- model.matrix(~ family - 1, data = tree_predictions)
  
  X_pred_tree <- cbind(numeric_features, species_dummies, genus_dummies, family_dummies)
  
  # Align with training features
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred_tree), ncol = ncol(X_tree))
  colnames(X_pred_aligned) <- colnames(X_tree)
  
  for (col_name in colnames(X_tree)) {
    if (col_name %in% colnames(X_pred_tree)) {
      X_pred_aligned[, col_name] <- X_pred_tree[, col_name]
    }
  }
  
  # Predict and back-transform
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  tree_predictions$flux_umol_m2_s <- sinh(pred_asinh)
  tree_predictions$flux_nmol_m2_s <- tree_predictions$flux_umol_m2_s * 1000
  
  # Calculate size for plotting (proportional to DBH)
  tree_predictions$tree_size <- (tree_predictions$dbh_m * 100)^2 / 50
  
  # =============================================================================
  # 2. SOIL FLUX PREDICTIONS ON GPS GRID
  # =============================================================================
  
  cat("  Predicting soil fluxes on GPS grid...\n")
  
  # Create GPS grid
  lon_range <- range(inventory_df$longitude, na.rm = TRUE)
  lat_range <- range(inventory_df$latitude, na.rm = TRUE)
  
  # Add small buffer
  lon_buffer <- diff(lon_range) * 0.05
  lat_buffer <- diff(lat_range) * 0.05
  
  lon_grid <- seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, length.out = soil_grid_res)
  lat_grid <- seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, length.out = soil_grid_res)
  
  soil_grid_gps <- expand.grid(longitude = lon_grid, latitude = lat_grid)
  
  # Convert GPS grid points to local x,y for moisture lookup
  xy_result <- geodetic_inverse(
    soil_grid_gps$latitude, soil_grid_gps$longitude,
    lat1_rad, lon1_rad, rotation_angle, earth_radius
  )
  soil_grid_gps$x <- xy_result$x
  soil_grid_gps$y <- xy_result$y
  
  # Get seasonal index for soil
  SI_soil_t <- SI_TABLES %>% 
    filter(group == "soil", month == month_num) %>% 
    pull(SI) %>% 
    {if(length(.) > 0) .[1] else 0}
  
  # Predict for each grid point
  soil_predictions <- soil_grid_gps %>%
    mutate(
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = clip_to_range(
        MOISTURE_AFFINE_TABLE$alpha_t[month_num] + 
          MOISTURE_AFFINE_TABLE$beta_t[month_num] * moisture_raw,
        training_ranges$moisture
      ),
      air_temp_C = clip_to_range(air_temp_t, training_ranges$air_temp),
      soil_temp_C = clip_to_range(soil_temp_t, training_ranges$soil_temp),
      SI_soil = SI_soil_t,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      moisture_x_airT = soil_moisture_abs * air_temp_C
    )
  
  # Build feature matrix
  X_pred_soil <- soil_predictions %>%
    select(soil_temp_C, soil_moisture_abs, SI_soil, 
           moisture_x_soilT, air_temp_C, moisture_x_airT) %>%
    as.matrix()
  
  # Align with training features
  X_pred_soil_aligned <- matrix(0, nrow = nrow(X_pred_soil), ncol = ncol(X_soil))
  colnames(X_pred_soil_aligned) <- colnames(X_soil)
  
  for (col_name in colnames(X_soil)) {
    if (col_name %in% colnames(X_pred_soil)) {
      X_pred_soil_aligned[, col_name] <- X_pred_soil[, col_name]
    }
  }
  
  # Predict and back-transform
  pred_asinh_soil <- predict(SoilRF, X_pred_soil_aligned)$predictions
  soil_predictions$flux_umol_m2_s <- sinh(pred_asinh_soil)
  soil_predictions$flux_nmol_m2_s <- soil_predictions$flux_umol_m2_s * 1000
  
  # =============================================================================
  # 3. CREATE PLOTS IN GPS COORDINATES
  # =============================================================================
  
  cat("  Creating visualizations...\n")
  
  # Common theme elements
  theme_map <- theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  # SOIL FLUX MAP
  p_soil <- ggplot() +
    # Soil flux as raster
    geom_raster(data = soil_predictions,
                aes(x = longitude, y = latitude, fill = flux_nmol_m2_s), 
                alpha = 0.8) +
    # Add tree locations as points for context
    geom_point(data = inventory_df,
               aes(x = longitude, y = latitude),
               size = 0.3, alpha = 0.3, color = "black") +
    # Color scale - diverging for soil (can be negative)
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = expression(paste("Flux\n(nmol m"^-2, " s"^-1, ")")),
      breaks = pretty_breaks(n = 5)
    ) +
    coord_equal() +
    labs(
      title = paste("Soil CH₄ Flux -", month.abb[month_num]),
      subtitle = "Negative values indicate CH₄ uptake",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_map
  
  # TREE FLUX MAP
  flux_range_tree <- range(tree_predictions$flux_nmol_m2_s, na.rm = TRUE)
  
  if (flux_range_tree[1] < 0 && flux_range_tree[2] > 0) {
    color_scale_tree <- scale_color_gradient2(
      low = "blue", mid = "white", high = "darkgreen",
      midpoint = 0,
      name = expression(paste("Flux\n(nmol m"^-2, " s"^-1, ")")),
      breaks = pretty_breaks(n = 5)
    )
  } else {
    color_scale_tree <- scale_color_viridis_c(
      name = expression(paste("Flux\n(nmol m"^-2, " s"^-1, ")")),
      breaks = pretty_breaks(n = 5),
      option = "plasma"
    )
  }
  
  p_tree <- ggplot() +
    # Tree points colored by flux, sized by DBH
    geom_point(data = tree_predictions,
               aes(x = longitude, y = latitude, 
                   color = flux_nmol_m2_s,
                   size = tree_size),
               alpha = 0.8) +
    color_scale_tree +
    scale_size_continuous(
      name = "DBH (cm)",
      range = c(1, 8),
      breaks = c(10, 50, 100, 200)^2 / 50,
      labels = c("10", "50", "100", "200")
    ) +
    coord_equal() +
    labs(
      title = paste("Tree Stem CH₄ Flux -", month.abb[month_num]),
      subtitle = paste("n =", nrow(tree_predictions), "trees (geodetically transformed)"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_map
  
  # COMBINED VIEW - Trees overlaid on soil
  p_combined <- ggplot() +
    # Soil flux as base layer
    geom_raster(data = soil_predictions,
                aes(x = longitude, y = latitude, fill = flux_nmol_m2_s),
                alpha = 0.7) +
    # Trees as points on top
    geom_point(data = tree_predictions,
               aes(x = longitude, y = latitude, size = tree_size),
               color = "black", fill = "white", shape = 21,
               alpha = 0.8) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = expression(paste("Soil Flux\n(nmol m"^-2, " s"^-1, ")")),
      breaks = pretty_breaks(n = 5)
    ) +
    scale_size_continuous(
      name = "Tree DBH (cm)",
      range = c(0.5, 6),
      breaks = c(10, 50, 100, 200)^2 / 50,
      labels = c("10", "50", "100", "200")
    ) +
    coord_equal() +
    labs(
      title = paste("Combined Flux Map -", month.abb[month_num]),
      subtitle = "Soil flux (colors) with tree locations (circles)",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_map
  
  # Return all three plots and the data
  return(list(
    plots = list(
      soil = p_soil,
      tree = p_tree,
      combined = p_combined
    ),
    data = list(
      soil = soil_predictions,
      tree = tree_predictions
    ),
    month = month_num
  ))
}

# =============================================================================
# GENERATE MAPS FOR SPECIFIC MONTHS
# =============================================================================

# Select months to visualize
months_to_plot <- c(1, 4, 7, 10)  # Winter, Spring, Summer, Fall
cat("\nGenerating maps for months:", paste(month.abb[months_to_plot], collapse = ", "), "\n")

# Generate maps for each selected month
monthly_maps <- lapply(months_to_plot, create_monthly_flux_maps)

# Save individual month plots
for (i in 1:length(monthly_maps)) {
  month_num <- monthly_maps[[i]]$month
  month_name <- month.abb[month_num]
  
  # Save soil map
  ggsave(
    filename = paste0("MAP_SOIL_FLUX_", month_name, ".png"),
    plot = monthly_maps[[i]]$plots$soil,
    width = 10, height = 8
  )
  
  # Save tree map
  ggsave(
    filename = paste0("MAP_TREE_FLUX_", month_name, ".png"),
    plot = monthly_maps[[i]]$plots$tree,
    width = 10, height = 8
  )
  
  # Save combined map
  ggsave(
    filename = paste0("MAP_COMBINED_FLUX_", month_name, ".png"),
    plot = monthly_maps[[i]]$plots$combined,
    width = 10, height = 8
  )
  
  cat("  Saved maps for", month_name, "\n")
}

# =============================================================================
# CREATE MULTI-PANEL SEASONAL COMPARISON
# =============================================================================

cat("\nCreating seasonal comparison plots...\n")

# Soil flux seasonal comparison
soil_plots <- lapply(monthly_maps, function(x) x$plots$soil)
p_soil_seasonal <- grid.arrange(
  grobs = soil_plots,
  ncol = 2,
  top = "Soil CH₄ Flux - Seasonal Variation"
)
ggsave("MAP_SOIL_SEASONAL.png", p_soil_seasonal, width = 16, height = 14)

# Tree flux seasonal comparison
tree_plots <- lapply(monthly_maps, function(x) x$plots$tree)
p_tree_seasonal <- grid.arrange(
  grobs = tree_plots,
  ncol = 2,
  top = "Tree Stem CH₄ Flux - Seasonal Variation"
)
ggsave("MAP_TREE_SEASONAL.png", p_tree_seasonal, width = 16, height = 14)

# =============================================================================
# ANNUAL AVERAGE MAPS
# =============================================================================

cat("\nGenerating annual average maps...\n")

# Calculate annual averages by running all 12 months
all_months_data <- lapply(1:12, function(m) {
  result <- create_monthly_flux_maps(m)
  return(result$data)
})

# Average soil flux across all months
soil_annual <- all_months_data[[1]]$soil %>%
  select(longitude, latitude) %>%
  mutate(flux_sum = 0, flux_count = 0)

for (m in 1:12) {
  month_data <- all_months_data[[m]]$soil
  soil_annual$flux_sum <- soil_annual$flux_sum + month_data$flux_nmol_m2_s
  soil_annual$flux_count <- soil_annual$flux_count + 1
}

soil_annual$flux_nmol_m2_s <- soil_annual$flux_sum / soil_annual$flux_count

# Average tree flux across all months
tree_annual <- all_months_data[[1]]$tree %>%
  select(tree_id, species, longitude, latitude, dbh_m) %>%
  mutate(flux_sum = 0, flux_count = 0)

for (m in 1:12) {
  month_data <- all_months_data[[m]]$tree
  tree_annual$flux_sum <- tree_annual$flux_sum + month_data$flux_nmol_m2_s
  tree_annual$flux_count <- tree_annual$flux_count + 1
}

tree_annual$flux_nmol_m2_s <- tree_annual$flux_sum / tree_annual$flux_count
tree_annual$tree_size <- (tree_annual$dbh_m * 100)^2 / 50

# Create annual average plots
p_soil_annual <- ggplot() +
  geom_raster(data = soil_annual,
              aes(x = longitude, y = latitude, fill = flux_nmol_m2_s),
              alpha = 0.8) +
  geom_point(data = INVENTORY,
             aes(x = longitude, y = latitude),
             size = 0.3, alpha = 0.3, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = expression(paste("Flux\n(nmol m"^-2, " s"^-1, ")")),
    breaks = pretty_breaks(n = 5)
  ) +
  coord_equal() +
  labs(
    title = "Annual Average Soil CH₄ Flux",
    subtitle = "Mean of monthly predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  )

p_tree_annual <- ggplot() +
  geom_point(data = tree_annual,
             aes(x = longitude, y = latitude, 
                 color = flux_nmol_m2_s,
                 size = tree_size),
             alpha = 0.8) +
  scale_color_viridis_c(
    name = expression(paste("Flux\n(nmol m"^-2, " s"^-1, ")")),
    breaks = pretty_breaks(n = 5),
    option = "plasma"
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(1, 8),
    breaks = c(10, 50, 100, 200)^2 / 50,
    labels = c("10", "50", "100", "200")
  ) +
  coord_equal() +
  labs(
    title = "Annual Average Tree Stem CH₄ Flux",
    subtitle = "Mean of monthly predictions (geodetically transformed)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  )

ggsave("MAP_SOIL_ANNUAL_AVERAGE.png", p_soil_annual, width = 10, height = 8)
ggsave("MAP_TREE_ANNUAL_AVERAGE.png", p_tree_annual, width = 10, height = 8)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== SPATIAL FLUX SUMMARY ===\n")
cat("Maps created with geodetic transformation:\n")
cat("  - PX/PY rotated by", rotation_angle * 180/pi, "degrees\n")
cat("  - Transformed to GPS using reference point:", ref_lat, "°N,", ref_lon, "°W\n")
cat("  - Ready for satellite imagery overlay\n\n")

# Check transformation success
n_transformed <- sum(!is.na(INVENTORY$longitude))
n_total <- nrow(INVENTORY)
cat("Coordinate transformation summary:\n")
cat("  Trees successfully mapped:", n_transformed, "out of", n_total, "\n")
if("PX_original" %in% names(INVENTORY)) {
  n_from_px_py <- sum(!is.na(INVENTORY$PX_original) & !is.na(INVENTORY$longitude))
  cat("  Trees transformed from PX/PY:", n_from_px_py, "\n")
}

# Soil statistics
soil_stats <- soil_annual %>%
  summarise(
    mean_flux = mean(flux_nmol_m2_s, na.rm = TRUE),
    median_flux = median(flux_nmol_m2_s, na.rm = TRUE),
    min_flux = min(flux_nmol_m2_s, na.rm = TRUE),
    max_flux = max(flux_nmol_m2_s, na.rm = TRUE),
    pct_negative = 100 * mean(flux_nmol_m2_s < 0, na.rm = TRUE)
  )

cat("\nSoil Flux Statistics (annual average, nmol m-2 s-1):\n")
cat("  Mean:", round(soil_stats$mean_flux, 3), "\n")
cat("  Median:", round(soil_stats$median_flux, 3), "\n")
cat("  Range: [", round(soil_stats$min_flux, 3), ",", round(soil_stats$max_flux, 3), "]\n")
cat("  % area with uptake (negative flux):", round(soil_stats$pct_negative, 1), "%\n")

# Tree statistics by species
tree_stats_species <- tree_annual %>%
  group_by(species) %>%
  summarise(
    n_trees = n(),
    mean_flux = mean(flux_nmol_m2_s, na.rm = TRUE),
    median_flux = median(flux_nmol_m2_s, na.rm = TRUE),
    total_dbh_m2 = sum(pi * (dbh_m/2)^2),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_flux))

cat("\nTree Flux by Species (annual average, nmol m-2 s-1):\n")
print(tree_stats_species %>% head(10))

# Overall tree statistics
tree_stats <- tree_annual %>%
  summarise(
    n_trees = n(),
    mean_flux = mean(flux_nmol_m2_s, na.rm = TRUE),
    median_flux = median(flux_nmol_m2_s, na.rm = TRUE),
    min_flux = min(flux_nmol_m2_s, na.rm = TRUE),
    max_flux = max(flux_nmol_m2_s, na.rm = TRUE),
    pct_negative = 100 * mean(flux_nmol_m2_s < 0, na.rm = TRUE)
  )

cat("\nOverall Tree Statistics (annual average, nmol m-2 s-1):\n")
cat("  Number of trees:", tree_stats$n_trees, "\n")
cat("  Mean flux:", round(tree_stats$mean_flux, 3), "\n")
cat("  Median flux:", round(tree_stats$median_flux, 3), "\n")
cat("  Range: [", round(tree_stats$min_flux, 3), ",", round(tree_stats$max_flux, 3), "]\n")
cat("  % trees with uptake:", round(tree_stats$pct_negative, 1), "%\n")

cat("\n✓ Spatial flux maps complete with geodetic transformation\n")
cat("Tree coordinates properly transformed from PX/PY to GPS\n")