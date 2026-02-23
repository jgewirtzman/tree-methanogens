# ==============================================================================
# Predict Tree CH4 Flux
# ==============================================================================
# Purpose: Applies trained RF models to predict tree CH4 flux across the
#   ForestGEO inventory for each month and species.
#
# Pipeline stage: 3 — Upscaling
# Run after: 02_rf_models.R
#
# Inputs:
#   - RF_MODELS.RData (from 02_rf_models)
#   - rf_workflow_input_data_with_2023.RData (from 01_load_and_prep_data)
#
# Outputs:
#   - tree_monthly_predictions.RData (to outputs/)
#   - tree_flux_predictions.csv (to outputs/tables/)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ranger)
library(patchwork)

cat("=== TREE CH4 FLUX MAPPING ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD MODELS AND DATA
# =============================================================================

cat("Loading RF models and data...\n")

# Load the trained models
load("../../outputs/models/RF_MODELS.RData")  # Contains TreeRF
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

# Extract what we need
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
TAXONOMY <- rf_workflow_data$PLACEHOLDER_TAXONOMY
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY

# Load tables
MOISTURE_AFFINE_TABLE <- read.csv("../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")
SI_TABLES <- read.csv("../../outputs/tables/SI_TABLES.csv")

cat("✓ Models and data loaded\n\n")

# =============================================================================
# SOURCE MAP SCRIPTS FOR CONSISTENT SPATIAL SETUP
# =============================================================================

cat("Sourcing map scripts for spatial setup and tree data...\n")

# # Source to get fg_final with all coordinate corrections and outlier removal
# invisible(capture.output(source("fg_aligned.R")))
# invisible(capture.output(source("map.R")))
# invisible(capture.output(source("complete_map2.R")))

# Use fg_final as our inventory - it has:
# - Geodetic transformation applied
# - Outlier removal
# - Proper coordinates (Longitude_final, Latitude_final)
# - Species names and codes
# - BasalArea_m2 calculated
INVENTORY <- fg_final %>%
  rename(
    tree_id = Tag,
    species = Species_Code,
    x = Longitude_final,
    y = Latitude_final,
    dbh_m = DBH  # DBH is in cm in fg_final
  ) %>%
  mutate(
    dbh_m = dbh_m / 100  # Convert cm to m for model
  )

# Add BasalArea if not already present (fg_final should have it)
if(!"BasalArea_m2" %in% names(INVENTORY)) {
  INVENTORY$BasalArea_m2 <- pi * (INVENTORY$dbh_m / 2)^2
} else {
  # Rename if it exists
  INVENTORY <- INVENTORY %>%
    rename(BasalArea_m2 = BasalArea_m2)
}

cat("  Using fg_final as inventory (", nrow(INVENTORY), "trees)\n")
cat("  DBH range (cm):", round(range(INVENTORY$dbh_m * 100), 1), "\n")
cat("  Using clipped moisture surface from complete_map2.R\n")
cat("  Using min_bbox for consistent boundaries\n\n")

# =============================================================================
# PREDICT FLUX FOR EACH TREE FOR A SPECIFIC MONTH
# =============================================================================

# Choose month to display (July = 7 for peak growing season)
DISPLAY_MONTH <- 7

cat("Predicting tree flux for month", month.name[DISPLAY_MONTH], "\n")

# Get drivers for this month
month_drivers <- DRIVERS %>% filter(month == DISPLAY_MONTH)

if (nrow(month_drivers) == 0 || is.na(month_drivers$soil_temp_C_mean)) {
  stop("No temperature data available for month", DISPLAY_MONTH)
}

# Get moisture calibration for this month
month_affine <- MOISTURE_AFFINE_TABLE %>% filter(month == DISPLAY_MONTH)

if (nrow(month_affine) == 0 || is.na(month_affine$alpha_t) || is.na(month_affine$beta_t)) {
  stop("No moisture calibration available for month", DISPLAY_MONTH)
}

# Get SI value for trees this month
si_month <- SI_TABLES %>% filter(group == "tree", month == DISPLAY_MONTH)
si_value <- ifelse(nrow(si_month) > 0, si_month$SI[1], 0)

cat("  Temperature data available\n")
cat("  Moisture calibration available\n")
cat("  SI value:", si_value, "\n\n")

# =============================================================================
# BUILD FEATURES FOR EACH INVENTORY TREE
# =============================================================================

cat("Building features for", nrow(INVENTORY), "inventory trees...\n")

# Get December moisture at each tree location
INVENTORY$moisture_dec <- moisture_lookup_xy(INVENTORY$x, INVENTORY$y)

# Apply affine transformation to get monthly moisture
INVENTORY$soil_moisture_at_tree <- month_affine$alpha_t[1] + 
  month_affine$beta_t[1] * INVENTORY$moisture_dec

# Clip moisture to reasonable bounds
INVENTORY$soil_moisture_at_tree <- pmax(0, pmin(0.6, INVENTORY$soil_moisture_at_tree))

# Add taxonomy info
INVENTORY <- INVENTORY %>%
  left_join(TAXONOMY %>% dplyr::select(species, genus, family), by = "species")

# Handle missing taxonomy
INVENTORY <- INVENTORY %>%
  mutate(
    genus = ifelse(is.na(genus), "Unknown", genus),
    family = ifelse(is.na(family), "Unknown", family)
  )

# Create one-hot encoding for species (matching training)
species_levels <- unique(c(INVENTORY$species, "Unknown", "OTHER"))
INVENTORY$species_factor <- factor(INVENTORY$species, levels = species_levels)
species_onehot <- model.matrix(~ species_factor - 1, data = INVENTORY)
colnames(species_onehot) <- paste0("species_", colnames(species_onehot))

# Create one-hot encoding for genus
genus_levels <- unique(c(INVENTORY$genus, "Unknown", "GENUS_OTHER"))
INVENTORY$genus_factor <- factor(INVENTORY$genus, levels = genus_levels)
genus_onehot <- model.matrix(~ genus_factor - 1, data = INVENTORY)
colnames(genus_onehot) <- paste0("genus_", colnames(genus_onehot))

# Create one-hot encoding for family
family_levels <- unique(c(INVENTORY$family, "Unknown", "FAMILY_OTHER"))
INVENTORY$family_factor <- factor(INVENTORY$family, levels = family_levels)
family_onehot <- model.matrix(~ family_factor - 1, data = INVENTORY)
colnames(family_onehot) <- paste0("family_", colnames(family_onehot))

# Build feature matrix
features <- data.frame(
  dbh_m = INVENTORY$dbh_m,
  air_temp_C_mean = month_drivers$air_temp_C_mean[1],
  soil_temp_C_mean = month_drivers$soil_temp_C_mean[1],
  soil_moisture_at_tree = INVENTORY$soil_moisture_at_tree,
  SI_tree = si_value,
  moisture_x_airT = INVENTORY$soil_moisture_at_tree * month_drivers$air_temp_C_mean[1],
  moisture_x_soilT = INVENTORY$soil_moisture_at_tree * month_drivers$soil_temp_C_mean[1],
  month_sin = sin(2 * pi * DISPLAY_MONTH / 12),
  month_cos = cos(2 * pi * DISPLAY_MONTH / 12),
  chamber_rigid = 0,  # Assume rigid chamber for predictions
  chamber_semirigid = 1  # Assume semirigid (most common)
)

# Add taxonomic one-hot features
features <- cbind(features, species_onehot, genus_onehot, family_onehot)

# Add moisture × species interactions
for(col in colnames(species_onehot)) {
  features[[paste0("moisture_x_", col)]] <- INVENTORY$soil_moisture_at_tree * species_onehot[, col]
}

cat("  Feature matrix created with", ncol(features), "features\n")

# =============================================================================
# PREDICT FLUX
# =============================================================================

cat("\nPredicting flux for all trees...\n")

# Get the feature names the model expects
model_features <- TreeRF$forest$independent.variable.names

# Create empty dataframe with all expected features
features_aligned <- data.frame(matrix(0, nrow = nrow(features), ncol = length(model_features)))
colnames(features_aligned) <- model_features

# Fill in the features we have
for(col in names(features)) {
  if(col %in% model_features) {
    features_aligned[[col]] <- features[[col]]
  }
}

# Predict (on asinh scale)
pred_asinh <- predict(TreeRF, features_aligned)$predictions

# Back-transform to get flux in μmol/m²/s
pred_flux_umol <- sinh(pred_asinh)

# Convert to nmol for display
INVENTORY$flux_nmol_m2_s <- pred_flux_umol * 1000

cat("  Predictions complete\n")
cat("  Flux range (nmol/m²/s):", round(range(INVENTORY$flux_nmol_m2_s), 3), "\n")
cat("  Mean flux:", round(mean(INVENTORY$flux_nmol_m2_s), 3), "\n")
cat("  Median flux:", round(median(INVENTORY$flux_nmol_m2_s), 3), "\n\n")

# =============================================================================
# CREATE FLUX MAP (SIMPLE VERSION)
# =============================================================================

cat("Creating tree flux map...\n")

# Create simple plot with all trees colored by flux
p_tree_flux <- ggplot() +
  # All trees colored by flux
  geom_point(data = INVENTORY, 
             aes(x = x, y = y, 
                 size = BasalArea_m2,
                 color = flux_nmol_m2_s), 
             alpha = 0.8, 
             stroke = 0.2) +
  
  scale_color_viridis_c(
    option = "plasma",
    name = paste("Tree CH₄ Flux\n(nmol/m²/s)\n", month.name[DISPLAY_MONTH])
  ) +
  
  scale_size_continuous(
    name = "Basal Area\n(m²)", 
    range = c(0.5, 4),
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
    guide = guide_legend(override.aes = list(shape = 16))
  ) +
  
  coord_equal() +
  labs(
    x = "Longitude", 
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm")
  )

print(p_tree_flux)
ggsave("../../outputs/figures/tree_flux_map.png", p_tree_flux, width = 12, height = 8, dpi = 300)
cat("  Saved: tree_flux_map.png\n")

# =============================================================================
# CREATE SUMMARY BY SPECIES
# =============================================================================

cat("\nCreating species summary...\n")

species_summary <- INVENTORY %>%
  group_by(species) %>%
  summarise(
    n_trees = n(),
    mean_flux = mean(flux_nmol_m2_s),
    median_flux = median(flux_nmol_m2_s),
    sd_flux = sd(flux_nmol_m2_s),
    mean_dbh_cm = mean(dbh_m * 100),
    total_basal_area_m2 = sum(pi * (dbh_m/2)^2),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_flux))

print(head(species_summary, 10))

# Save results
write.csv(INVENTORY %>% dplyr::select(tree_id, species, dbh_m, x, y, flux_nmol_m2_s),
          "../../outputs/flux_predictions/tree_flux_predictions.csv", row.names = FALSE)
write.csv(species_summary, "../../outputs/tables/tree_flux_by_species.csv", row.names = FALSE)

cat("\n=== TREE FLUX MAPPING COMPLETE ===\n")
cat("Map shows", month.name[DISPLAY_MONTH], "flux predictions\n")
cat("Blue = low flux, Yellow = medium, Red = high flux\n")
cat("Point size represents DBH\n")











# =============================================================================
# TREE CH4 FLUX MAPPING - EXTENDED WITH MONTHLY PREDICTIONS
# Matches soil code functionality with species-specific analysis
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ranger)
library(patchwork)

cat("=== TREE CH4 FLUX MAPPING - EXTENDED ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD MODELS AND DATA
# =============================================================================

cat("Loading RF models and data...\n")

# Load the trained models
load("../../outputs/models/RF_MODELS.RData")  # Contains TreeRF
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

# Extract what we need
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
TAXONOMY <- rf_workflow_data$PLACEHOLDER_TAXONOMY
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY

# Load tables
MOISTURE_AFFINE_TABLE <- read.csv("../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")
SI_TABLES <- read.csv("../../outputs/tables/SI_TABLES.csv")

cat("✓ Models and data loaded\n\n")

# =============================================================================
# SOURCE MAP SCRIPTS FOR CONSISTENT SPATIAL SETUP
# =============================================================================

cat("Sourcing map scripts for spatial setup and tree data...\n")

# Source to get fg_final with all coordinate corrections
invisible(capture.output(source('../../code/07_maps/01_forestgeo_alignment.R')))
invisible(capture.output(source('../../code/07_maps/02_spatial_interpolation.R')))
invisible(capture.output(source('../../code/07_maps/03_interpolation_methods.R')))
invisible(capture.output(source('../../code/07_maps/04_seasonal_flux_maps.R')))

# Use fg_final as our inventory
INVENTORY <- fg_final %>%
  rename(
    tree_id = Tag,
    species = Species_Code,
    x = Longitude_final,
    y = Latitude_final,
    dbh_m = DBH
  ) %>%
  mutate(
    dbh_m = dbh_m / 100  # Convert cm to m for model
  )

# Add BasalArea
if(!"BasalArea_m2" %in% names(INVENTORY)) {
  INVENTORY$BasalArea_m2 <- pi * (INVENTORY$dbh_m / 2)^2
}

# Add taxonomy info
INVENTORY <- INVENTORY %>%
  left_join(TAXONOMY %>% dplyr::select(species, genus, family), by = "species") %>%
  mutate(
    genus = ifelse(is.na(genus), "Unknown", genus),
    family = ifelse(is.na(family), "Unknown", family)
  )

# Get December moisture at each tree location once (for base moisture pattern)
INVENTORY$moisture_dec <- moisture_lookup_xy(INVENTORY$x, INVENTORY$y)

cat("  Using fg_final as inventory (", nrow(INVENTORY), "trees)\n")
cat("  DBH range (cm):", round(range(INVENTORY$dbh_m * 100), 1), "\n")
cat("  Species count:", length(unique(INVENTORY$species)), "\n\n")

# =============================================================================
# PREDICTION FUNCTION FOR TREES BY MONTH
# =============================================================================

predict_tree_month <- function(month_val, inventory_df, moisture_affine, drivers, model, si_table) {
  
  cat("  Month", month_val, "...")
  
  # Get climate drivers
  month_drivers <- drivers %>% filter(month == month_val)
  
  if (nrow(month_drivers) == 0 || is.na(month_drivers$soil_temp_C_mean)) {
    cat(" skipped (no temperature data)\n")
    return(NULL)
  }
  
  # Get moisture calibration
  month_affine <- moisture_affine %>% filter(month == month_val)
  
  if (nrow(month_affine) == 0 || is.na(month_affine$alpha_t) || is.na(month_affine$beta_t)) {
    cat(" skipped (no moisture calibration)\n")
    return(NULL)
  }
  
  # Get SI value for trees
  si_month <- si_table %>% filter(group == "tree", month == month_val)
  si_value <- ifelse(nrow(si_month) > 0, si_month$SI[1], 0)
  
  # Apply affine transformation to get monthly moisture
  inventory_df$soil_moisture_at_tree <- month_affine$alpha_t[1] + 
    month_affine$beta_t[1] * inventory_df$moisture_dec
  
  # Clip moisture to reasonable bounds
  inventory_df$soil_moisture_at_tree <- pmax(0, pmin(0.6, inventory_df$soil_moisture_at_tree))
  
  # Create one-hot encoding for species
  species_levels <- unique(c(inventory_df$species, "Unknown", "OTHER"))
  inventory_df$species_factor <- factor(inventory_df$species, levels = species_levels)
  species_onehot <- model.matrix(~ species_factor - 1, data = inventory_df)
  colnames(species_onehot) <- paste0("species_", gsub("species_factor", "", colnames(species_onehot)))
  
  # Create one-hot encoding for genus
  genus_levels <- unique(c(inventory_df$genus, "Unknown", "GENUS_OTHER"))
  inventory_df$genus_factor <- factor(inventory_df$genus, levels = genus_levels)
  genus_onehot <- model.matrix(~ genus_factor - 1, data = inventory_df)
  colnames(genus_onehot) <- paste0("genus_", gsub("genus_factor", "", colnames(genus_onehot)))
  
  # Create one-hot encoding for family
  family_levels <- unique(c(inventory_df$family, "Unknown", "FAMILY_OTHER"))
  inventory_df$family_factor <- factor(inventory_df$family, levels = family_levels)
  family_onehot <- model.matrix(~ family_factor - 1, data = inventory_df)
  colnames(family_onehot) <- paste0("family_", gsub("family_factor", "", colnames(family_onehot)))
  
  # Build feature matrix
  features <- data.frame(
    dbh_m = inventory_df$dbh_m,
    air_temp_C_mean = month_drivers$air_temp_C_mean[1],
    soil_temp_C_mean = month_drivers$soil_temp_C_mean[1],
    soil_moisture_at_tree = inventory_df$soil_moisture_at_tree,
    SI_tree = si_value,
    moisture_x_airT = inventory_df$soil_moisture_at_tree * month_drivers$air_temp_C_mean[1],
    moisture_x_soilT = inventory_df$soil_moisture_at_tree * month_drivers$soil_temp_C_mean[1],
    month_sin = sin(2 * pi * month_val / 12),
    month_cos = cos(2 * pi * month_val / 12),
    chamber_rigid = 0,
    chamber_semirigid = 1
  )
  
  # Add taxonomic one-hot features
  features <- cbind(features, species_onehot, genus_onehot, family_onehot)
  
  # Add moisture × species interactions
  for(col in colnames(species_onehot)) {
    features[[paste0("moisture_x_", col)]] <- inventory_df$soil_moisture_at_tree * species_onehot[, col]
  }
  
  # Align features with model expectations
  model_features <- model$forest$independent.variable.names
  features_aligned <- data.frame(matrix(0, nrow = nrow(features), ncol = length(model_features)))
  colnames(features_aligned) <- model_features
  
  for(col in names(features)) {
    if(col %in% model_features) {
      features_aligned[[col]] <- features[[col]]
    }
  }
  
  # Predict
  tryCatch({
    pred_asinh <- predict(model, features_aligned)$predictions
    pred_flux <- sinh(pred_asinh)
    
    result <- inventory_df %>%
      dplyr::select(tree_id, species, x, y, dbh_m, BasalArea_m2) %>%
      mutate(
        month = month_val,
        flux_umol_m2_s = pred_flux,
        flux_nmol_m2_s = pred_flux * 1000,
        flux_mg_m2_d = pred_flux * 86400 * 16 * 1e-3,
        moisture_pct = inventory_df$soil_moisture_at_tree * 100
      )
    
    cat(" done (mean:", round(mean(pred_flux * 1000), 3), "nmol)\n")
    return(result)
    
  }, error = function(e) {
    cat(" failed:", e$message, "\n")
    return(NULL)
  })
}

# =============================================================================
# GENERATE PREDICTIONS FOR ALL MONTHS
# =============================================================================

cat("\nGenerating monthly tree flux predictions...\n")

# Get usable months (same as soil code)
complete_drivers <- DRIVERS %>%
  filter(!is.na(soil_temp_C_mean) & !is.na(air_temp_C_mean))
valid_affine <- MOISTURE_AFFINE_TABLE %>%
  filter(!is.na(alpha_t) & !is.na(beta_t))
usable_months <- intersect(complete_drivers$month, valid_affine$month)

cat("  Usable months:", paste(usable_months, collapse = ", "), "\n\n")

# Generate predictions for each month
prediction_list <- list()
for (m in usable_months) {
  result <- predict_tree_month(m, INVENTORY, MOISTURE_AFFINE_TABLE, 
                               DRIVERS, TreeRF, SI_TABLES)
  if (!is.null(result)) {
    prediction_list[[length(prediction_list) + 1]] <- result
  }
}

if (length(prediction_list) == 0) {
  stop("No successful predictions")
}

monthly_predictions <- bind_rows(prediction_list)
cat("\n✓ Predictions complete for", length(unique(monthly_predictions$month)), "months\n")

# =============================================================================
# CREATE SEASONAL MAPS (ALL TREES)
# =============================================================================

cat("\nCreating seasonal tree flux maps...\n")

# Function to create month map for all trees
create_tree_month_map <- function(month_val, pred_df) {
  month_data <- pred_df %>% filter(month == month_val)
  
  # Stats
  mean_flux <- mean(month_data$flux_nmol_m2_s)
  median_flux <- median(month_data$flux_nmol_m2_s)
  
  # Color limits based on all months
  flux_limits <- quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98))
  
  ggplot(month_data, aes(x = x, y = y)) +
    geom_point(aes(size = BasalArea_m2, color = flux_nmol_m2_s), 
               alpha = 0.7, stroke = 0.1) +
    scale_color_viridis_c(
      option = "plasma",
      limits = flux_limits,
      oob = scales::squish,
      name = "nmol/m²/s"
    ) +
    scale_size_continuous(
      name = "Basal Area (m²)",
      range = c(0.3, 3),
      guide = "none"  # Hide size legend for cleaner look
    ) +
    coord_equal() +
    labs(
      title = paste(month.abb[month_val], "- Tree CH₄ Flux"),
      subtitle = paste("Mean:", round(mean_flux, 2), "| Median:", 
                       round(median_flux, 2), "nmol/m²/s"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}

# dplyr::select seasonal months (same as soil code)
available_months <- sort(unique(monthly_predictions$month))
season_months <- c(
  winter = intersect(c(12, 1, 2), available_months)[1],
  spring = intersect(c(3, 4, 5), available_months)[1], 
  summer = intersect(c(6, 7, 8), available_months)[1],
  fall = intersect(c(9, 10, 11), available_months)[1]
)
season_months <- na.omit(season_months)

# Create seasonal plots
if (length(season_months) >= 2) {
  plot_list <- lapply(season_months, create_tree_month_map, monthly_predictions)
  
  if (length(plot_list) == 4) {
    combined <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]])
  } else if (length(plot_list) == 3) {
    combined <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[3]])
  } else if (length(plot_list) == 2) {
    combined <- plot_list[[1]] | plot_list[[2]]
  }
  
  ggsave("../../outputs/figures/tree_flux_seasonal.png", combined, width = 14, height = 10, dpi = 150)
  cat("  Saved: tree_flux_seasonal.png\n")
}

# =============================================================================
# ANNUAL AVERAGE MAP (ALL TREES)
# =============================================================================

cat("\nCreating annual average map...\n")

annual_avg <- monthly_predictions %>%
  group_by(tree_id, species, x, y, dbh_m, BasalArea_m2) %>%
  summarise(
    mean_flux_nmol = mean(flux_nmol_m2_s),
    mean_flux_mg = mean(flux_mg_m2_d),
    mean_moisture = mean(moisture_pct),
    .groups = "drop"
  )

p_annual <- ggplot(annual_avg, aes(x = x, y = y)) +
  geom_point(aes(size = BasalArea_m2, color = mean_flux_nmol), 
             alpha = 0.7, stroke = 0.1) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Annual Mean\nFlux (nmol/m²/s)"
  ) +
  scale_size_continuous(
    name = "Basal Area (m²)",
    range = c(0.3, 3),
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Annual Average Tree CH₄ Flux",
    subtitle = paste("Mean:", round(mean(annual_avg$mean_flux_nmol), 3), 
                     "nmol/m²/s | N =", nrow(annual_avg), "trees"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/tree_flux_annual.png", p_annual, width = 10, height = 8, dpi = 150)
cat("  Saved: tree_flux_annual.png\n")

# =============================================================================
# SPECIES-SPECIFIC FACETED MAPS (SEASONAL)
# =============================================================================

cat("\nCreating species-specific seasonal maps...\n")

# Get top species by abundance
species_counts <- annual_avg %>%
  count(species) %>%
  arrange(desc(n))

# dplyr::select top 12 species for faceting
top_species <- head(species_counts$species, 12)

# Filter data for top species
species_seasonal_data <- monthly_predictions %>%
  filter(species %in% top_species,
         month %in% season_months)

if (nrow(species_seasonal_data) > 0) {
  # Create faceted plot by species and season
  species_seasonal_data$month_label <- factor(
    month.abb[species_seasonal_data$month],
    levels = month.abb[season_months]
  )
  
  p_species_seasonal <- ggplot(species_seasonal_data, 
                               aes(x = x, y = y, color = flux_nmol_m2_s)) +
    geom_point(aes(size = BasalArea_m2), alpha = 0.6, stroke = 0) +
    facet_grid(species ~ month_label, scales = "free") +
    scale_color_viridis_c(
      option = "plasma",
      name = "Flux\n(nmol/m²/s)"
    ) +
    scale_size_continuous(
      range = c(0.2, 2),
      guide = "none"
    ) +
    coord_equal() +
    labs(
      title = "Seasonal Tree CH₄ Flux by Species",
      subtitle = "Top 12 species by abundance",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 7),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave("../../outputs/figures/tree_flux_species_seasonal.png", p_species_seasonal, 
         width = 12, height = 14, dpi = 150)
  cat("  Saved: tree_flux_species_seasonal.png\n")
}

# =============================================================================
# SPECIES-SPECIFIC FACETED MAPS (ANNUAL AVERAGE)
# =============================================================================

cat("\nCreating species-specific annual average maps...\n")

# Filter annual data for top species
species_annual_data <- annual_avg %>%
  filter(species %in% top_species)

if (nrow(species_annual_data) > 0) {
  p_species_annual <- ggplot(species_annual_data, 
                             aes(x = x, y = y, color = mean_flux_nmol)) +
    geom_point(aes(size = BasalArea_m2), alpha = 0.7, stroke = 0) +
    facet_wrap(~ species, ncol = 4) +
    scale_color_viridis_c(
      option = "plasma",
      name = "Annual Mean\nFlux (nmol/m²/s)"
    ) +
    scale_size_continuous(
      range = c(0.2, 2),
      guide = "none"
    ) +
    coord_equal() +
    labs(
      title = "Annual Average Tree CH₄ Flux by Species",
      subtitle = paste("Top 12 species | Overall mean:", 
                       round(mean(species_annual_data$mean_flux_nmol), 3), "nmol/m²/s"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave("../../outputs/figures/tree_flux_species_annual.png", p_species_annual, 
         width = 12, height = 9, dpi = 150)
  cat("  Saved: tree_flux_species_annual.png\n")
}

# =============================================================================
# DETAILED SPECIES COMPARISON (HIGHLIGHTED MONTHS)
# =============================================================================

cat("\nCreating detailed species comparison for highlighted months...\n")

# dplyr::select 6 most abundant species for detailed view
top_6_species <- head(species_counts$species, 6)

species_highlight_data <- monthly_predictions %>%
  filter(species %in% top_6_species,
         month %in% season_months)

if (nrow(species_highlight_data) > 0) {
  # Calculate species means by month
  species_month_means <- species_highlight_data %>%
    group_by(species, month) %>%
    summarise(
      mean_flux = mean(flux_nmol_m2_s),
      .groups = "drop"
    )
  
  # Create individual maps for each species-month combination
  for (sp in top_6_species) {
    sp_data <- species_highlight_data %>% filter(species == sp)
    
    if (nrow(sp_data) > 0) {
      plot_list <- list()
      
      for (i in 1:length(season_months)) {
        m <- season_months[i]
        month_sp_data <- sp_data %>% filter(month == m)
        
        if (nrow(month_sp_data) > 0) {
          mean_flux <- mean(month_sp_data$flux_nmol_m2_s)
          
          plot_list[[i]] <- ggplot(month_sp_data, aes(x = x, y = y)) +
            geom_point(aes(size = BasalArea_m2, color = flux_nmol_m2_s), 
                       alpha = 0.8) +
            scale_color_viridis_c(
              option = "plasma",
              limits = quantile(sp_data$flux_nmol_m2_s, c(0.02, 0.98)),
              oob = scales::squish,
              name = "nmol/m²/s"
            ) +
            scale_size_continuous(range = c(0.5, 3), guide = "none") +
            coord_equal() +
            labs(
              title = month.abb[m],
              subtitle = paste("Mean:", round(mean_flux, 2))
            ) +
            theme_minimal() +
            theme(legend.position = "right")
        }
      }
      
      if (length(plot_list) >= 2) {
        if (length(plot_list) == 4) {
          combined_sp <- (plot_list[[1]] | plot_list[[2]]) / 
            (plot_list[[3]] | plot_list[[4]]) +
            plot_annotation(
              title = paste("Species:", sp, "- Seasonal CH₄ Flux"),
              subtitle = paste("N =", nrow(sp_data %>% distinct(tree_id)), "trees")
            )
        } else if (length(plot_list) == 2) {
          combined_sp <- plot_list[[1]] | plot_list[[2]] +
            plot_annotation(
              title = paste("Species:", sp, "- Seasonal CH₄ Flux")
            )
        }
        
        ggsave(paste0("../../outputs/figures/tree_flux_", sp, "_seasonal.png"), combined_sp, 
               width = 12, height = 8, dpi = 150)
        cat("  Saved: tree_flux_", sp, "_seasonal.png\n")
      }
    }
  }
}

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

# Monthly summary
month_summary <- monthly_predictions %>%
  group_by(month) %>%
  summarise(
    mean_flux_nmol = mean(flux_nmol_m2_s),
    median_flux_nmol = median(flux_nmol_m2_s),
    sd_flux_nmol = sd(flux_nmol_m2_s),
    min_flux_nmol = min(flux_nmol_m2_s),
    max_flux_nmol = max(flux_nmol_m2_s),
    .groups = "drop"
  ) %>%
  mutate(month_name = month.abb[month])

cat("\nMonthly flux summary:\n")
print(month_summary)

# Species summary
species_summary <- annual_avg %>%
  group_by(species) %>%
  summarise(
    n_trees = n(),
    mean_flux = mean(mean_flux_nmol),
    median_flux = median(mean_flux_nmol),
    sd_flux = sd(mean_flux_nmol),
    mean_dbh_cm = mean(dbh_m * 100),
    total_basal_area_m2 = sum(BasalArea_m2),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_flux))

cat("\nTop 10 species by mean annual flux:\n")
print(head(species_summary, 10))

# Species by month summary
species_month_summary <- monthly_predictions %>%
  group_by(species, month) %>%
  summarise(
    n_trees = n_distinct(tree_id),
    mean_flux = mean(flux_nmol_m2_s),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = month, values_from = mean_flux, names_prefix = "month_")

# Save all summaries
write.csv(month_summary, "../../outputs/tables/tree_flux_monthly_summary.csv", row.names = FALSE)
write.csv(species_summary, "../../outputs/tables/tree_flux_species_summary.csv", row.names = FALSE)
write.csv(species_month_summary, "../../outputs/tables/tree_flux_species_by_month.csv", row.names = FALSE)
write.csv(annual_avg, "../../outputs/flux_predictions/tree_flux_annual_all_trees.csv", row.names = FALSE)

cat("\n=== TREE FLUX MAPPING COMPLETE ===\n")
cat("Generated predictions for", length(unique(monthly_predictions$month)), "months\n")
cat("Created seasonal and annual maps for all trees\n")
cat("Created species-specific faceted maps\n")
cat("All output files saved\n")

# Save the tree monthly predictions object for later use  
save(monthly_predictions, file = "../../outputs/models/tree_monthly_predictions.RData")











# =============================================================================
# FACETED ANNUAL FLUX DISTRIBUTION BY SPECIES
# =============================================================================

cat("\nCreating faceted annual flux distribution by species...\n")

# Filter for species with n > 50 and exclude KALA
species_for_dist <- species_counts %>%
  filter(n > 50, species != "KALA") %>%
  pull(species)

cat("  Species with n>50 (excluding KALA):", length(species_for_dist), "species\n")
cat("  Species included:", paste(species_for_dist, collapse = ", "), "\n")

# Filter annual data for selected species
species_dist_data <- annual_avg %>%
  filter(species %in% species_for_dist) %>%
  # Add species count for ordering
  group_by(species) %>%
  mutate(n_trees = n()) %>%
  ungroup() %>%
  # Create a label with species code and count
  mutate(species_label = paste0(species, " (n=", n_trees, ")"))

# Create faceted histogram/density plot
p_species_dist <- ggplot(species_dist_data, 
                         aes(x = mean_flux_nmol)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 15, 
                 fill = "steelblue", 
                 alpha = 0.7,
                 color = "white") +
  geom_density(color = "darkred", 
               linewidth = 0.8,
               alpha = 0.5) +
  facet_wrap(~ species_label, 
             scales = "free_y", 
             ncol = 4) +  # Adjusted ncol based on expected number of species
  labs(
    title = "Distribution of Mean Annual CH₄ Flux by Species",
    subtitle = paste(length(species_for_dist), 
                     "species with n>50 trees (KALA excluded)"),
    x = "Mean Annual Flux (nmol/m²/s)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank()
  )

ggsave("../../outputs/figures/tree_flux_distribution_by_species.png", p_species_dist, 
       width = 14, height = 10, dpi = 150)
cat("  Saved: tree_flux_distribution_by_species.png\n")

# Alternative: Box plot version for cleaner comparison
p_species_box <- ggplot(species_dist_data, 
                        aes(x = mean_flux_nmol, y = reorder(species, mean_flux_nmol))) +
  geom_boxplot(aes(fill = species), 
               alpha = 0.7, 
               show.legend = FALSE) +
  geom_point(aes(color = species), 
             position = position_jitter(height = 0.2), 
             size = 1, 
             alpha = 0.5,
             show.legend = FALSE) +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Mean Annual CH₄ Flux Distribution by Species",
    subtitle = paste(length(species_for_dist), 
                     "species with n>50 trees (KALA excluded) | Points show individual trees"),
    x = "Mean Annual Flux (nmol/m²/s)",
    y = "Species Code"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, face = "bold")
  )

ggsave("../../outputs/figures/tree_flux_boxplot_by_species.png", p_species_box, 
       width = 10, height = 8, dpi = 150)  # Adjusted height based on expected species count
cat("  Saved: tree_flux_boxplot_by_species.png\n")

# Create a faceted scatter plot showing spatial distribution with flux
p_species_spatial <- ggplot(species_dist_data, 
                            aes(x = x, y = y)) +
  geom_point(aes(color = mean_flux_nmol, size = BasalArea_m2), 
             alpha = 0.7) +
  facet_wrap(~ species, ncol = 4) +  # Adjusted ncol
  scale_color_viridis_c(
    option = "plasma",
    name = "Annual Mean\nFlux (nmol/m²/s)"
  ) +
  scale_size_continuous(
    range = c(0.5, 2),
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Spatial Distribution of Mean Annual CH₄ Flux by Species",
    subtitle = paste(length(species_for_dist), 
                     "species with n>50 trees (KALA excluded)"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave("../../outputs/figures/tree_flux_spatial_faceted_species.png", p_species_spatial, 
       width = 14, height = 10, dpi = 150)
cat("  Saved: tree_flux_spatial_faceted_species.png\n")

# Summary statistics for the faceted species
faceted_species_stats <- species_dist_data %>%
  group_by(species) %>%
  summarise(
    n_trees = n(),
    mean_flux = mean(mean_flux_nmol),
    median_flux = median(mean_flux_nmol),
    sd_flux = sd(mean_flux_nmol),
    cv_flux = sd_flux / mean_flux * 100,  # Coefficient of variation
    min_flux = min(mean_flux_nmol),
    max_flux = max(mean_flux_nmol),
    q25_flux = quantile(mean_flux_nmol, 0.25),
    q75_flux = quantile(mean_flux_nmol, 0.75),
    iqr_flux = q75_flux - q25_flux,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_flux))

cat("\nStatistics for species with n>50 (excluding KALA):\n")
print(faceted_species_stats)
write.csv(faceted_species_stats, "../../outputs/tables/tree_flux_faceted_species_stats.csv", row.names = FALSE)

# Also create a simple bar chart showing mean flux by species
p_species_bar <- ggplot(faceted_species_stats, 
                        aes(x = reorder(species, mean_flux), y = mean_flux)) +
  geom_bar(stat = "identity", aes(fill = mean_flux), show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_flux - sd_flux, ymax = mean_flux + sd_flux),
                width = 0.3, alpha = 0.5) +
  scale_fill_viridis_c(option = "plasma") +
  coord_flip() +
  labs(
    title = "Mean Annual CH₄ Flux by Species (±1 SD)",
    subtitle = paste(length(species_for_dist), 
                     "species with n>50 trees (KALA excluded)"),
    x = "Species Code",
    y = "Mean Annual Flux (nmol/m²/s)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, face = "bold")
  )

ggsave("../../outputs/figures/tree_flux_barplot_by_species.png", p_species_bar, 
       width = 10, height = 8, dpi = 150)
cat("  Saved: tree_flux_barplot_by_species.png\n")

cat("\n✓ Faceted annual flux plots by species complete\n")
cat("  Total trees in analysis:", nrow(species_dist_data), "\n")
cat("  Number of species analyzed:", length(species_for_dist), "\n")

