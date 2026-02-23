# ==============================================================================
# ForestGEO Flux Maps
# ==============================================================================
# Purpose: Generates spatial maps of predicted CH4 flux overlaid on the
#   ForestGEO inventory with species-level detail.
#
# Pipeline stage: 4 — Visualization
# Run after: 03_predict_tree_flux.R
#
# Inputs:
#   - RF_MODELS.RData, TAXONOMY_PRIORS.RData (from 02_rf_models)
#   - rf_workflow_input_data_with_2023.RData (from 01_load_and_prep_data)
#
# Outputs:
#   - ForestGEO flux map PNGs (to outputs/figures/)
#   - flux prediction CSVs (to outputs/tables/)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(scales)
library(ranger)

cat("=== CH4 FLUX PREDICTION FOR FORESTGEO TREES ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# STEP 1: LOAD THE EXISTING fg_final DATASET
# =============================================================================

# First run your fg_aligned.R script to create fg_final
source('../../code/07_maps/01_forestgeo_alignment.R')

# Verify fg_final exists
if(!exists("fg_final")) {
  stop("ERROR: fg_final not found. Please run fg_aligned.R first to create the dataset.")
}

cat("Loaded fg_final dataset with", nrow(fg_final), "trees\n")
cat("  Coordinate range:\n")
cat("    Latitude:", round(range(fg_final$Latitude_final), 6), "\n")
cat("    Longitude:", round(range(fg_final$Longitude_final), 6), "\n\n")

# =============================================================================
# STEP 2: LOAD RF MODELS AND SUPPORTING DATA
# =============================================================================

cat("Loading RF models and calibration data...\n")

# Load the saved models and data
load("../../outputs/models/RF_MODELS.RData")  # Contains TreeRF and SoilRF
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")  # Contains all input data
load("../../outputs/models/TAXONOMY_PRIORS.RData")  # Contains taxonomy priors

# Extract components
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY
SI_TABLES <- read.csv("../../outputs/tables/SI_TABLES.csv")
MOISTURE_AFFINE_TABLE <- read.csv("../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")

cat("✓ Models and data loaded\n\n")

# =============================================================================
# STEP 3: PREPARE INVENTORY FROM fg_final
# =============================================================================

# Use fg_final directly as the inventory
INVENTORY <- fg_final %>%
  mutate(
    # Rename to match expected format
    tree_id = as.character(Tag),
    species = Species_Name,
    dbh_m = DBH / 100,  # Convert cm to meters
    latitude = Latitude_final,
    longitude = Longitude_final,
    x = PX,  # Local coordinates for moisture lookup
    y = PY,
    basal_area_m2 = BasalArea_m2,
    dataset = Coord_Source,
    coordinates_estimated = Coordinates_Estimated
  ) %>%
  filter(!is.na(latitude) & !is.na(longitude) & !is.na(dbh_m))

cat("=== INVENTORY PREPARED ===\n")
cat("  Trees in inventory:", nrow(INVENTORY), "\n")
cat("  Species represented:", length(unique(INVENTORY$species)), "\n")
cat("  DBH range:", round(range(INVENTORY$dbh_m * 100), 1), "cm\n")
cat("  Trees with estimated coords:", sum(INVENTORY$coordinates_estimated), "\n\n")

# =============================================================================
# STEP 4: PREDICT CH4 FLUX
# =============================================================================

# Choose month (1-12 or 0 for annual average)
MONTH_TO_MAP <- 7  # July

cat("Predicting CH4 flux for month", MONTH_TO_MAP, "\n")

# Prediction function
predict_tree_flux <- function(inv_df, month_num) {
  
  # Get monthly drivers
  air_temp <- DRIVERS$air_temp_C_mean[month_num]
  soil_temp <- DRIVERS$soil_temp_C_mean[month_num]
  if(is.na(soil_temp)) soil_temp <- air_temp * 0.8
  
  # Get seasonal index
  si_tree <- SI_TABLES %>% 
    filter(group == "tree", month == month_num) %>% 
    pull(SI) %>% 
    {if(length(.) > 0) .[1] else 0}
  
  # Prepare features
  predictions <- inv_df %>%
    mutate(
      # Moisture at tree location
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = MOISTURE_AFFINE_TABLE$alpha_t[month_num] + 
        MOISTURE_AFFINE_TABLE$beta_t[month_num] * moisture_raw,
      soil_moisture_abs = pmax(0, pmin(0.75, soil_moisture_abs)),
      
      # Environmental
      air_temp_C = air_temp,
      soil_temp_C = soil_temp,
      SI_tree = si_tree,
      
      # Cyclic month encoding
      month_sin = sin(2 * pi * month_num / 12),
      month_cos = cos(2 * pi * month_num / 12),
      
      # Chamber type
      chamber_rigid = 1,
      chamber_semirigid = 0,
      
      # Interactions
      moisture_x_airT = soil_moisture_abs * air_temp_C,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      
      # Taxonomy
      genus = sub(" .*", "", species),
      family = case_when(
        genus == "Acer" ~ "Sapindaceae",
        genus == "Betula" ~ "Betulaceae",
        genus %in% c("Fagus", "Quercus") ~ "Fagaceae",
        genus == "Fraxinus" ~ "Oleaceae",
        genus %in% c("Pinus", "Tsuga") ~ "Pinaceae",
        genus == "Carya" ~ "Juglandaceae",
        TRUE ~ "OTHER"
      )
    )
  
  # Add taxonomy priors
  predictions$taxon_prior_asinh <- 0
  
  if(exists("TAXONOMY_PRIORS")) {
    species_priors <- setNames(TAXONOMY_PRIORS$species$med_resid, 
                               TAXONOMY_PRIORS$species$species)
    genus_priors <- setNames(TAXONOMY_PRIORS$genus$med_resid, 
                             TAXONOMY_PRIORS$genus$genus)
    
    for(i in 1:nrow(predictions)) {
      if(predictions$species[i] %in% names(species_priors)) {
        predictions$taxon_prior_asinh[i] <- species_priors[predictions$species[i]]
      } else if(predictions$genus[i] %in% names(genus_priors)) {
        predictions$taxon_prior_asinh[i] <- genus_priors[predictions$genus[i]]
      }
    }
  }
  
  # Build feature matrix
  numeric_features <- predictions %>%
    dplyr::select(dbh_m, air_temp_C, soil_temp_C, soil_moisture_abs, SI_tree,
           moisture_x_airT, moisture_x_soilT, taxon_prior_asinh,
           chamber_rigid, chamber_semirigid, month_sin, month_cos) %>%
    as.matrix()
  
  # Add categorical dummies
  species_dummies <- model.matrix(~ species - 1, data = predictions)
  genus_dummies <- model.matrix(~ genus - 1, data = predictions)
  family_dummies <- model.matrix(~ family - 1, data = predictions)
  
  X_pred <- cbind(numeric_features, species_dummies, genus_dummies, family_dummies)
  
  # Align with model features
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred), 
                           ncol = length(TreeRF$forest$independent.variable.names))
  colnames(X_pred_aligned) <- TreeRF$forest$independent.variable.names
  
  for(col in colnames(X_pred_aligned)) {
    if(col %in% colnames(X_pred)) {
      X_pred_aligned[, col] <- X_pred[, col]
    }
  }
  
  # Predict and back-transform
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  return(sinh(pred_asinh))
}

# Apply predictions
INVENTORY$flux_umol_m2_s <- predict_tree_flux(INVENTORY, MONTH_TO_MAP)
INVENTORY$flux_nmol_m2_s <- INVENTORY$flux_umol_m2_s * 1000

# Add categories
INVENTORY <- INVENTORY %>%
  mutate(
    flux_category = case_when(
      flux_nmol_m2_s < 0 ~ "Sink",
      flux_nmol_m2_s < 1 ~ "Very Low",
      flux_nmol_m2_s < 5 ~ "Low",
      flux_nmol_m2_s < 10 ~ "Moderate",
      flux_nmol_m2_s < 20 ~ "High",
      TRUE ~ "Very High"
    ),
    flux_category = factor(flux_category,
                           levels = c("Sink", "Very Low", "Low", "Moderate", "High", "Very High"))
  )

cat("✓ Flux predictions complete\n")
cat("  Range:", round(range(INVENTORY$flux_nmol_m2_s), 2), "nmol m-2 s-1\n")
cat("  Trees as sinks:", sum(INVENTORY$flux_nmol_m2_s < 0), 
    "(", round(100 * mean(INVENTORY$flux_nmol_m2_s < 0), 1), "%)\n\n")

# =============================================================================
# STEP 5: CREATE VISUALIZATIONS
# =============================================================================

cat("Creating maps...\n")

# Map 1: All trees with flux
p1 <- ggplot(INVENTORY, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = flux_nmol_m2_s, size = dbh_m * 100), alpha = 0.7) +
  scale_color_viridis(
    name = expression(paste("CH"[4], " Flux\n(nmol m"^-2, " s"^-1, ")")),
    option = "turbo",
    trans = "pseudo_log",
    breaks = c(-1, 0, 1, 5, 10, 20)
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.3, 4),
    breaks = c(5, 10, 20, 40, 60)
  ) +
  coord_equal() +
  labs(
    title = paste("ForestGEO CH4 Flux -", month.name[MONTH_TO_MAP]),
    subtitle = paste("N =", nrow(INVENTORY), "trees (geodetically transformed)"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/forestgeo_flux_map.png", p1, width = 12, height = 9, dpi = 300)

# Map 2: By species (top emitters)
species_summary <- INVENTORY %>%
  filter(species != "Unknown") %>%
  group_by(species) %>%
  summarise(
    n = n(),
    mean_flux = mean(flux_nmol_m2_s),
    total_basal = sum(basal_area_m2),
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%
  arrange(desc(mean_flux))

top_species <- head(species_summary$species, 9)

p2 <- INVENTORY %>%
  filter(species %in% top_species) %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point(data = INVENTORY, color = "grey90", size = 0.2, alpha = 0.3) +
  geom_point(aes(color = flux_nmol_m2_s, size = dbh_m * 100), alpha = 0.8) +
  facet_wrap(~ species, ncol = 3) +
  scale_color_viridis(
    name = "Flux",
    option = "plasma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(range = c(0.5, 3)) +
  coord_equal() +
  labs(title = "CH4 Flux by Species") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic", size = 9)
  )

ggsave("../../outputs/figures/forestgeo_flux_by_species.png", p2, width = 12, height = 10, dpi = 300)

# Map 3: Hotspot analysis
library(hexbin)

p3 <- ggplot(INVENTORY, aes(x = longitude, y = latitude)) +
  stat_summary_hex(aes(z = flux_nmol_m2_s), fun = mean, bins = 30, alpha = 0.8) +
  geom_point(data = filter(INVENTORY, flux_nmol_m2_s > quantile(flux_nmol_m2_s, 0.95)),
             aes(size = dbh_m * 100), color = "red", shape = 17, alpha = 0.7) +
  scale_fill_viridis(
    name = expression(paste("Mean Flux\n(nmol m"^-2, " s"^-1, ")")),
    option = "magma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(1, 4)
  ) +
  coord_equal() +
  labs(
    title = "CH4 Flux Hotspots",
    subtitle = "Hexbins = mean flux, Red triangles = top 5% emitters"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/forestgeo_flux_hotspots.png", p3, width = 10, height = 8, dpi = 300)

cat("✓ Maps saved\n\n")

# =============================================================================
# STEP 6: SUMMARY STATISTICS
# =============================================================================

cat("=== SUMMARY STATISTICS ===\n\n")

# Overall
cat("Overall flux statistics:\n")
cat("  Mean:", round(mean(INVENTORY$flux_nmol_m2_s), 2), "nmol m-2 s-1\n")
cat("  Median:", round(median(INVENTORY$flux_nmol_m2_s), 2), "nmol m-2 s-1\n")
cat("  SD:", round(sd(INVENTORY$flux_nmol_m2_s), 2), "nmol m-2 s-1\n\n")

# By species
cat("Top 10 emitting species:\n")
print(head(species_summary, 10))

# By size class
size_summary <- INVENTORY %>%
  mutate(
    size_class = cut(dbh_m * 100, 
                     breaks = c(0, 10, 20, 40, 100),
                     labels = c("<10cm", "10-20cm", "20-40cm", ">40cm"))
  ) %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    mean_flux = round(mean(flux_nmol_m2_s), 2),
    pct_sink = round(100 * mean(flux_nmol_m2_s < 0), 1),
    .groups = "drop"
  )

cat("\nFlux by tree size:\n")
print(size_summary)

# Save results
output <- INVENTORY %>%
  dplyr::select(tree_id, species, latitude, longitude, dbh_m, basal_area_m2,
         flux_umol_m2_s, flux_nmol_m2_s, flux_category,
         transformation_method, coordinates_estimated)

write.csv(output,
          paste0("../../outputs/flux_predictions/forestgeo_flux_predictions_month", MONTH_TO_MAP, ".csv"),
          row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Used fg_final dataset with", nrow(INVENTORY), "trees\n")
cat("All coordinates properly transformed using geodetic methods\n")
cat("Output files saved successfully\n")



















# =============================================================================
# TREE CH4 FLUX MAP USING CORRECTED fg_final DATASET
# Uses the properly cleaned and geodetically transformed data from fg_aligned.R
# Run AFTER running fg_aligned.R to create fg_final
# =============================================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(scales)
library(hexbin)

cat("=== TREE CH4 FLUX PREDICTION USING CORRECTED DATA ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# STEP 1: LOAD THE CORRECTED fg_final DATASET
# =============================================================================

cat("Loading corrected ForestGEO data...\n")

# Run the existing fg_aligned.R script that creates the properly corrected fg_final
# This script handles ALL data cleaning including:
# - Merging fg19, fgplot, fgtag datasets
# - Outlier detection and removal (0.75 × IQR)
# - DBH error corrections
# - Geodetic coordinate transformation
# - Duplicate removal
# - Species name mapping

source('../../code/07_maps/01_forestgeo_alignment.R')

if(!exists("fg_final")) {
  stop("ERROR: fg_final not found. Please ensure fg_aligned.R runs successfully.")
}

cat("✓ Loaded fg_final with", nrow(fg_final), "trees (after outlier removal)\n")
cat("  Geodetic transformation:", sum(fg_final$transformation_method == "Geodetic Transform"), "trees\n")
cat("  Original GPS:", sum(fg_final$transformation_method == "Original GPS"), "trees\n")
cat("  Coordinate ranges:\n")
cat("    Latitude:", round(range(fg_final$Latitude_final), 6), "\n")
cat("    Longitude:", round(range(fg_final$Longitude_final), 6), "\n\n")

# =============================================================================
# STEP 2: LOAD RF MODELS AND SUPPORTING DATA
# =============================================================================

cat("Loading RF models and calibration data...\n")

load("../../outputs/models/RF_MODELS.RData")  # Contains TreeRF and SoilRF
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")  # Contains all input data
load("../../outputs/models/TAXONOMY_PRIORS.RData")  # Contains taxonomy priors

# Extract components
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY
SI_TABLES <- read.csv("../../outputs/tables/SI_TABLES.csv")
MOISTURE_AFFINE_TABLE <- read.csv("../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")

cat("✓ Models loaded\n\n")

# =============================================================================
# STEP 3: PREPARE INVENTORY FROM fg_final
# =============================================================================

cat("Preparing inventory for flux predictions...\n")

INVENTORY <- fg_final %>%
  mutate(
    tree_id = as.character(Tag),
    species = Species_Name,
    species_code = Species_Code,
    latitude = Latitude_final,
    longitude = Longitude_final,
    x = PX,  # Local coordinates for moisture lookup
    y = PY,
    dbh_m = DBH / 100,  # Convert cm to meters
    basal_area_m2 = BasalArea_m2,
    dataset = Coord_Source,
    coordinates_estimated = Coordinates_Estimated,
    transformation_method = transformation_method
  ) %>%
  filter(!is.na(latitude) & !is.na(longitude) & !is.na(dbh_m))

cat("  Trees in inventory:", nrow(INVENTORY), "\n")
cat("  Species:", length(unique(INVENTORY$species[INVENTORY$species != "Unknown"])), "\n")
cat("  DBH range:", round(min(INVENTORY$dbh_m * 100), 1), "-", 
    round(max(INVENTORY$dbh_m * 100), 1), "cm\n\n")

# =============================================================================
# STEP 4: PREDICT CH4 FLUX
# =============================================================================

# Choose month to visualize
MONTH_TO_MAP <- 7  # July

cat("Predicting CH4 flux for month", MONTH_TO_MAP, "(", month.name[MONTH_TO_MAP], ")...\n")

predict_tree_flux <- function(inv_df, month_num) {
  
  # Monthly drivers
  air_temp <- DRIVERS$air_temp_C_mean[month_num]
  soil_temp <- DRIVERS$soil_temp_C_mean[month_num]
  if(is.na(soil_temp)) soil_temp <- air_temp * 0.8
  
  # Seasonal index
  si_tree <- SI_TABLES %>% 
    filter(group == "tree", month == month_num) %>% 
    pull(SI) %>% 
    {if(length(.) > 0) .[1] else 0}
  
  # Build features
  predictions <- inv_df %>%
    mutate(
      # Moisture
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = MOISTURE_AFFINE_TABLE$alpha_t[month_num] + 
        MOISTURE_AFFINE_TABLE$beta_t[month_num] * moisture_raw,
      soil_moisture_abs = pmax(0, pmin(0.75, soil_moisture_abs)),
      
      # Environmental
      air_temp_C = air_temp,
      soil_temp_C = soil_temp,
      SI_tree = si_tree,
      
      # Cyclic month
      month_sin = sin(2 * pi * month_num / 12),
      month_cos = cos(2 * pi * month_num / 12),
      
      # Chamber
      chamber_rigid = 1,
      chamber_semirigid = 0,
      
      # Interactions
      moisture_x_airT = soil_moisture_abs * air_temp_C,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      
      # Taxonomy
      genus = sub(" .*", "", species),
      family = case_when(
        genus == "Acer" ~ "Sapindaceae",
        genus == "Betula" ~ "Betulaceae",
        genus %in% c("Fagus", "Quercus") ~ "Fagaceae",
        genus == "Fraxinus" ~ "Oleaceae",
        genus %in% c("Pinus", "Tsuga") ~ "Pinaceae",
        genus == "Carya" ~ "Juglandaceae",
        genus == "Prunus" ~ "Rosaceae",
        genus == "Sassafras" ~ "Lauraceae",
        genus == "Hamamelis" ~ "Hamamelidaceae",
        genus == "Kalmia" ~ "Ericaceae",
        genus == "Vaccinium" ~ "Ericaceae",
        genus == "Liriodendron" ~ "Magnoliaceae",
        TRUE ~ "OTHER"
      )
    )
  
  # Add taxonomy priors
  predictions$taxon_prior_asinh <- 0
  if(exists("TAXONOMY_PRIORS")) {
    species_priors <- setNames(TAXONOMY_PRIORS$species$med_resid, TAXONOMY_PRIORS$species$species)
    genus_priors <- setNames(TAXONOMY_PRIORS$genus$med_resid, TAXONOMY_PRIORS$genus$genus)
    family_priors <- setNames(TAXONOMY_PRIORS$family$med_resid, TAXONOMY_PRIORS$family$family)
    
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
  }
  
  # Build feature matrix
  numeric_features <- predictions %>%
    dplyr::select(dbh_m, air_temp_C, soil_temp_C, soil_moisture_abs, SI_tree,
           moisture_x_airT, moisture_x_soilT, taxon_prior_asinh,
           chamber_rigid, chamber_semirigid, month_sin, month_cos) %>%
    as.matrix()
  
  # Categorical dummies
  species_dummies <- model.matrix(~ species - 1, data = predictions)
  genus_dummies <- model.matrix(~ genus - 1, data = predictions)
  family_dummies <- model.matrix(~ family - 1, data = predictions)
  
  X_pred <- cbind(numeric_features, species_dummies, genus_dummies, family_dummies)
  
  # Align with model
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred), 
                           ncol = length(TreeRF$forest$independent.variable.names))
  colnames(X_pred_aligned) <- TreeRF$forest$independent.variable.names
  
  for(col in colnames(X_pred_aligned)) {
    if(col %in% colnames(X_pred)) {
      X_pred_aligned[, col] <- X_pred[, col]
    }
  }
  
  # Predict
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  return(sinh(pred_asinh))
}

# Apply predictions
INVENTORY$flux_umol_m2_s <- predict_tree_flux(INVENTORY, MONTH_TO_MAP)
INVENTORY$flux_nmol_m2_s <- INVENTORY$flux_umol_m2_s * 1000

# Categorize flux
INVENTORY <- INVENTORY %>%
  mutate(
    flux_category = case_when(
      flux_nmol_m2_s < 0 ~ "Sink",
      flux_nmol_m2_s < 1 ~ "Very Low",
      flux_nmol_m2_s < 5 ~ "Low", 
      flux_nmol_m2_s < 10 ~ "Moderate",
      flux_nmol_m2_s < 20 ~ "High",
      TRUE ~ "Very High"
    ),
    flux_category = factor(flux_category,
                           levels = c("Sink", "Very Low", "Low", 
                                      "Moderate", "High", "Very High"))
  )

cat("✓ Flux predictions complete\n")
cat("  Range:", round(min(INVENTORY$flux_nmol_m2_s), 2), "-",
    round(max(INVENTORY$flux_nmol_m2_s), 2), "nmol m⁻² s⁻¹\n")
cat("  Trees as sinks:", sum(INVENTORY$flux_nmol_m2_s < 0),
    "(", round(100 * mean(INVENTORY$flux_nmol_m2_s < 0), 1), "%)\n\n")

# =============================================================================
# STEP 5: CREATE VISUALIZATIONS
# =============================================================================

cat("Creating visualizations...\n")

# MAP 1: Overall flux map
p1 <- ggplot(INVENTORY, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = flux_nmol_m2_s, size = dbh_m * 100), alpha = 0.7) +
  scale_color_viridis(
    name = expression(paste("CH"[4], " Flux\n(nmol m"^-2, " s"^-1, ")")),
    option = "turbo",
    trans = "pseudo_log",
    breaks = c(-1, 0, 1, 5, 10, 20)
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.3, 4),
    breaks = c(5, 10, 20, 40, 60)
  ) +
  coord_equal() +
  labs(
    title = paste("ForestGEO CH4 Flux -", month.name[MONTH_TO_MAP]),
    subtitle = paste("N =", nrow(INVENTORY), "trees (corrected & geodetically transformed)"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/fg_flux_map_overall.png", p1, width = 12, height = 9, dpi = 300)

# MAP 2: Species-specific flux (each facet shows ONLY that species)
species_summary <- INVENTORY %>%
  filter(species != "Unknown") %>%
  group_by(species) %>%
  summarise(
    n = n(),
    mean_flux = mean(flux_nmol_m2_s),
    median_flux = median(flux_nmol_m2_s),
    total_basal = sum(basal_area_m2),
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%
  arrange(desc(mean_flux))

top_species <- head(species_summary$species, 12)

p2 <- INVENTORY %>%
  filter(species %in% top_species) %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point(aes(color = flux_nmol_m2_s, size = dbh_m * 100), alpha = 0.8) +
  facet_wrap(~ species, ncol = 3) +
  scale_color_viridis(
    name = expression(paste("CH"[4], " Flux")),
    option = "plasma",
    trans = "pseudo_log",
    breaks = c(0, 1, 5, 10, 20)
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.5, 3),
    breaks = c(10, 20, 40),
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "CH4 Flux by Species",
    subtitle = paste("Top 12 species (n ≥ 10) - Each panel shows only that species")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic", size = 9),
    panel.grid = element_blank()
  )

ggsave("../../outputs/figures/fg_flux_by_species.png", p2, width = 14, height = 12, dpi = 300)

# MAP 3: Alternative species view with context
# Background shows all tree locations very faintly
background_trees <- INVENTORY %>%
  dplyr::select(longitude, latitude) %>%
  expand_grid(facet_species = top_species) %>%
  rename(species = facet_species)

p2b <- INVENTORY %>%
  filter(species %in% top_species) %>%
  ggplot(aes(x = longitude, y = latitude)) +
  # Very subtle context layer
  geom_point(data = background_trees,
             color = "grey95", size = 0.05, alpha = 0.2) +
  # Actual species points
  geom_point(aes(color = flux_nmol_m2_s, size = dbh_m * 100), alpha = 0.9) +
  facet_wrap(~ species, ncol = 3) +
  scale_color_viridis(
    name = expression(paste("CH"[4], " Flux")),
    option = "plasma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  coord_equal() +
  labs(
    title = "CH4 Flux by Species (with spatial context)",
    subtitle = "Colored = species flux, grey = all tree locations"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic", size = 9),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("../../outputs/figures/fg_flux_by_species_context.png", p2b, width = 14, height = 12, dpi = 300)

# MAP 4: Hotspot analysis
p3 <- ggplot(INVENTORY, aes(x = longitude, y = latitude)) +
  stat_summary_hex(aes(z = flux_nmol_m2_s), fun = mean, bins = 30, alpha = 0.8) +
  geom_point(data = filter(INVENTORY, flux_nmol_m2_s > quantile(flux_nmol_m2_s, 0.95)),
             aes(size = dbh_m * 100), color = "red", shape = 17, alpha = 0.7) +
  scale_fill_viridis(
    name = expression(paste("Mean Flux\n(nmol m"^-2, " s"^-1, ")")),
    option = "magma",
    trans = "pseudo_log"
  ) +
  scale_size_continuous(
    name = "Top 5% DBH",
    range = c(1, 4),
    breaks = c(10, 20, 40)
  ) +
  coord_equal() +
  labs(
    title = "CH4 Flux Hotspots",
    subtitle = "Hexbins show mean flux, red triangles = top 5% emitters"
  ) +
  theme_minimal()

ggsave("../../outputs/figures/fg_flux_hotspots.png", p3, width = 10, height = 8, dpi = 300)

# MAP 5: Data source comparison
p4 <- ggplot(INVENTORY, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = dataset, size = dbh_m * 100, alpha = flux_nmol_m2_s)) +
  scale_color_manual(
    name = "Data Source",
    values = c("fg19" = "darkgreen", "fgplot" = "orange", 
               "fgtag" = "purple", "fgtag_missing" = "red", 
               "fgplot_missing" = "pink")
  ) +
  scale_size_continuous(
    name = "DBH (cm)",
    range = c(0.3, 3)
  ) +
  scale_alpha_continuous(
    name = "Flux",
    range = c(0.3, 0.9),
    trans = "pseudo_log"
  ) +
  coord_equal() +
  labs(
    title = "Tree Inventory by Data Source",
    subtitle = paste("Showing", nrow(INVENTORY), "trees from merged datasets")
  ) +
  theme_minimal()

ggsave("../../outputs/figures/fg_flux_by_source.png", p4, width = 12, height = 9, dpi = 300)

cat("✓ Visualizations saved\n\n")

# =============================================================================
# STEP 6: SUMMARY STATISTICS
# =============================================================================

cat("=== SUMMARY STATISTICS ===\n\n")

# Overall
cat("OVERALL:\n")
cat("  Mean flux:", round(mean(INVENTORY$flux_nmol_m2_s), 2), "nmol m⁻² s⁻¹\n")
cat("  Median flux:", round(median(INVENTORY$flux_nmol_m2_s), 2), "nmol m⁻² s⁻¹\n")
cat("  SD:", round(sd(INVENTORY$flux_nmol_m2_s), 2), "\n")
cat("  95% CI:", round(quantile(INVENTORY$flux_nmol_m2_s, 0.025), 2), "-",
    round(quantile(INVENTORY$flux_nmol_m2_s, 0.975), 2), "\n\n")

# Top emitting species
cat("TOP 10 EMITTING SPECIES (mean flux):\n")
print(head(species_summary %>% 
             dplyr::select(species, n, mean_flux, median_flux, total_basal), 10))
cat("\n")

# By size class
size_summary <- INVENTORY %>%
  mutate(
    size_class = cut(dbh_m * 100,
                     breaks = c(0, 10, 20, 40, 100),
                     labels = c("<10cm", "10-20cm", "20-40cm", ">40cm"))
  ) %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    mean_flux = round(mean(flux_nmol_m2_s), 2),
    pct_sink = round(100 * mean(flux_nmol_m2_s < 0), 1),
    .groups = "drop"
  )

cat("BY SIZE CLASS:\n")
print(size_summary)
cat("\n")

# By data source
source_summary <- INVENTORY %>%
  group_by(dataset) %>%
  summarise(
    n = n(),
    mean_flux = round(mean(flux_nmol_m2_s), 2),
    median_flux = round(median(flux_nmol_m2_s), 2),
    .groups = "drop"
  )

cat("BY DATA SOURCE:\n")
print(source_summary)

# Save results
output <- INVENTORY %>%
  dplyr::select(tree_id, species, latitude, longitude, x, y,
         dbh_m, basal_area_m2, flux_umol_m2_s, flux_nmol_m2_s,
         flux_category, dataset, transformation_method, coordinates_estimated)

write.csv(output,
          paste0("../../outputs/flux_predictions/fg_flux_predictions_month", MONTH_TO_MAP, ".csv"),
          row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Used fg_final with", nrow(INVENTORY), "trees (outliers already removed)\n")
cat("All coordinates properly corrected using geodetic transformation\n")
cat("Results saved to fg_flux_predictions_month", MONTH_TO_MAP, ".csv\n")
