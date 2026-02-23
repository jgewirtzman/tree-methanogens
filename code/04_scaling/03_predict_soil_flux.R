# ==============================================================================
# Predict Soil CH4 Flux
# ==============================================================================
# Purpose: Applies trained RF models to predict soil CH4 flux with spatial
#   interpolation across the study area.
#
# Pipeline stage: 3 — Upscaling
# Run after: 02_rf_models.R
#
# Inputs:
#   - RF_MODELS.RData (from 02_rf_models)
#   - rf_workflow_input_data_with_2023.RData (from 01_load_and_prep_data)
#
# Outputs:
#   - soil_monthly_predictions.RData (to outputs/)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ranger)
library(akima)
library(patchwork)

cat("=== SOIL CH4 FLUX MAPPING - EXTENDED SURFACE ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD REQUIRED DATA
# =============================================================================

cat("Loading models and data...\n")

load("../../outputs/models/RF_MODELS.RData")
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

# Extract needed components
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
INVENTORY <- rf_workflow_data$PLACEHOLDER_INVENTORY
MOISTURE_DEC_RASTER <- rf_workflow_data$PLACEHOLDER_MOISTURE_DEC_RASTER

MOISTURE_AFFINE_TABLE <- read.csv("../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")
SI_TABLES <- read.csv("../../outputs/tables/SI_TABLES.csv")

cat("✓ Models and data loaded\n\n")

# =============================================================================
# SOURCE EXISTING SCRIPTS TO GET ALL OBJECTS
# =============================================================================

cat("Sourcing existing scripts...\n")

# Source your scripts that create all the objects we need
# This gives us best_df (extended moisture) and min_bbox (exact bounding box)

cat("  Sourcing fg_aligned.R...\n")
invisible(capture.output(source('../../code/07_maps/01_forestgeo_alignment.R')))

cat("  Sourcing map.R...\n")
invisible(capture.output(source('../../code/07_maps/02_spatial_interpolation.R')))

cat("  Sourcing complete_map2.R...\n") 
invisible(capture.output(source('../../code/07_maps/03_interpolation_methods.R'
)))

cat("  Sourcing complete_map2.R...\n") 
invisible(capture.output(source('../../code/07_maps/04_seasonal_flux_maps.R')))


cat("  Objects loaded successfully\n")
cat("  Using best_df as extended moisture surface\n")
cat("  Using min_bbox for clipping\n\n")

# =============================================================================
# APPLY EXISTING CLIPPING TO MOISTURE SURFACE
# =============================================================================

cat("Applying clipping boundaries from complete_map2.R...\n")

# We already have clipped_moisture_df from complete_map2.R which is exactly what we need!
# It's best_df (extended Akima with river) clipped to min_bbox boundaries

cat("  Using pre-clipped moisture surface from complete_map2.R\n")
cat("  Original grid points (best_df):", nrow(best_df), "\n")
cat("  Clipped grid points:", nrow(clipped_moisture_df), "\n")
cat("  Retained:", round(100 * nrow(clipped_moisture_df) / nrow(best_df), 1), "%\n")
cat("  Bounding box angle:", round(min_bbox$angle, 1), "degrees\n\n")

# =============================================================================
# CREATE PREDICTION GRID FOR SOIL FLUX
# =============================================================================

cat("Creating soil flux prediction grid...\n")

# Use the clipped moisture grid as base for predictions
pred_grid <- clipped_moisture_df %>%
  dplyr::select(Longitude, Latitude, VWC) %>%
  rename(x = Longitude, y = Latitude, moisture_dec = VWC)

cat("  Prediction points:", nrow(pred_grid), "\n\n")

# =============================================================================
# PREDICTION FUNCTION FOR SOIL FLUX
# =============================================================================

predict_soil_month_extended <- function(month_val, grid_df, moisture_affine, drivers, model, si_table) {
  
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
  
  # Get SI value
  si_month <- si_table %>% filter(group == "soil", month == month_val)
  si_value <- ifelse(nrow(si_month) > 0, si_month$SI[1], 0)
  
  # Apply affine transformation to December moisture
  # Note: grid_df$moisture_dec is the extended Akima interpolation (not actual December)
  # We'll treat it as the base moisture pattern to be transformed
  grid_df$soil_moisture_at_site <- month_affine$alpha_t[1] + 
    month_affine$beta_t[1] * (grid_df$moisture_dec / 100)  # Normalize to 0-1
  
  # Clip moisture to reasonable bounds
  grid_df$soil_moisture_at_site <- pmax(0, pmin(0.6, grid_df$soil_moisture_at_site))
  
  # Build features matching model expectations
  features <- data.frame(
    soil_temp_C_mean = month_drivers$soil_temp_C_mean[1],
    air_temp_C_mean = month_drivers$air_temp_C_mean[1],
    soil_moisture_at_site = grid_df$soil_moisture_at_site,
    SI = si_value,
    moisture_x_soilT = grid_df$soil_moisture_at_site * month_drivers$soil_temp_C_mean[1],
    moisture_x_airT = grid_df$soil_moisture_at_site * month_drivers$air_temp_C_mean[1],
    month_sin = sin(2 * pi * month_val / 12),
    month_cos = cos(2 * pi * month_val / 12)
  )
  
  # Predict
  tryCatch({
    pred_asinh <- predict(model, features)$predictions
    pred_flux <- sinh(pred_asinh)
    
    result <- grid_df %>%
      dplyr::select(x, y) %>%
      mutate(
        month = month_val,
        flux_umol_m2_s = pred_flux,
        flux_nmol_m2_s = pred_flux * 1000,
        flux_mg_m2_d = pred_flux * 86400 * 16 * 1e-3,
        moisture_pct = grid_df$soil_moisture_at_site * 100,
        uptake = pred_flux < 0
      )
    
    cat(" done (mean:", round(mean(pred_flux * 1000), 3), "nmol)\n")
    return(result)
    
  }, error = function(e) {
    cat(" failed:", e$message, "\n")
    return(NULL)
  })
}

# =============================================================================
# GENERATE PREDICTIONS
# =============================================================================

cat("\nGenerating monthly soil flux predictions...\n")

# Get usable months
complete_drivers <- DRIVERS %>%
  filter(!is.na(soil_temp_C_mean) & !is.na(air_temp_C_mean))
valid_affine <- MOISTURE_AFFINE_TABLE %>%
  filter(!is.na(alpha_t) & !is.na(beta_t))
usable_months <- intersect(complete_drivers$month, valid_affine$month)

cat("  Usable months:", paste(usable_months, collapse = ", "), "\n\n")

# Generate predictions for each month
prediction_list <- list()
for (m in usable_months) {
  result <- predict_soil_month_extended(m, pred_grid, MOISTURE_AFFINE_TABLE, 
                                        DRIVERS, SoilRF, SI_TABLES)
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
# CREATE MAPS
# =============================================================================

cat("\nCreating soil flux maps...\n")

# Function to create month map
create_month_map <- function(month_val, pred_df) {
  month_data <- pred_df %>% filter(month == month_val)
  
  # Stats
  mean_flux <- mean(month_data$flux_nmol_m2_s)
  pct_uptake <- 100 * mean(month_data$uptake)
  
  # Color limits
  flux_limits <- quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98))
  
  ggplot(month_data, aes(x = x, y = y, fill = flux_nmol_m2_s)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white", 
      high = "red",
      midpoint = 0,
      limits = flux_limits,
      oob = scales::squish,
      name = "nmol/m²/s"
    ) +
    coord_equal() +
    labs(
      title = paste(month.abb[month_val], "- Soil CH₄ Flux"),
      subtitle = paste("Mean:", round(mean_flux, 2), "nmol/m²/s | Uptake:", 
                       round(pct_uptake, 1), "%"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
}

# dplyr::select seasonal months
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
  plot_list <- lapply(season_months, create_month_map, monthly_predictions)
  
  if (length(plot_list) == 4) {
    combined <- (plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]])
  } else if (length(plot_list) == 2) {
    combined <- plot_list[[1]] | plot_list[[2]]
  }
  
  ggsave("../../outputs/figures/soil_flux_extended_seasonal.png", combined, width = 14, height = 10, dpi = 150)
  cat("  Saved: soil_flux_extended_seasonal.png\n")
}

# =============================================================================
# ANNUAL AVERAGE MAP
# =============================================================================

cat("\nCreating annual average map...\n")

annual_avg <- monthly_predictions %>%
  group_by(x, y) %>%
  summarise(
    mean_flux_nmol = mean(flux_nmol_m2_s),
    mean_flux_mg = mean(flux_mg_m2_d),
    pct_uptake = 100 * mean(uptake),
    mean_moisture = mean(moisture_pct),
    .groups = "drop"
  )

p_annual <- ggplot(annual_avg, aes(x = x, y = y, fill = mean_flux_nmol)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "nmol/m²/s"
  ) +
  coord_equal() +
  labs(
    title = "Annual Average Soil CH₄ Flux - Extended Surface",
    subtitle = paste("Mean:", round(mean(annual_avg$mean_flux_nmol), 3), "nmol/m²/s |",
                     "Uptake:", round(mean(annual_avg$pct_uptake), 1), "%"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/soil_flux_extended_annual.png", p_annual, width = 10, height = 8, dpi = 150)
cat("  Saved: soil_flux_extended_annual.png\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== SUMMARY ===\n")

month_summary <- monthly_predictions %>%
  group_by(month) %>%
  summarise(
    mean_flux_nmol = mean(flux_nmol_m2_s),
    median_flux_nmol = median(flux_nmol_m2_s),
    pct_uptake = 100 * mean(uptake),
    .groups = "drop"
  ) %>%
  mutate(month_name = month.abb[month])

print(month_summary)

cat("\nAnnual statistics:\n")
cat("  Mean flux:", round(mean(annual_avg$mean_flux_nmol), 3), "nmol/m²/s\n")
cat("  Uptake frequency:", round(mean(annual_avg$pct_uptake), 1), "%\n")
cat("  Grid points:", nrow(annual_avg), "\n")

# Save data
write.csv(month_summary, "../../outputs/tables/soil_flux_extended_summary.csv", row.names = FALSE)
write.csv(annual_avg, "../../outputs/flux_predictions/soil_flux_extended_annual.csv", row.names = FALSE)

cat("\n=== SOIL FLUX MAPPING COMPLETE ===\n")
cat("Used extended Akima interpolation with river points\n")
cat("Clipped to match final map boundaries\n")


# Save the monthly predictions object for later use
save(monthly_predictions, file = "../../outputs/models/soil_monthly_predictions.RData")