# ==============================================================================
# Random Forest CH4 Flux Models
# ==============================================================================
# Purpose: Trains Random Forest models for tree and soil CH4 flux prediction
#   with a species-first deconfounding approach.
#
# Pipeline stage: 3 — Upscaling
# Run after: 01_load_and_prep_data.R
#
# Inputs:
#   - rf_workflow_input_data_with_2023.RData (from 01_load_and_prep_data)
#
# Outputs:
#   - RF_MODELS.RData, TAXONOMY_PRIORS.RData (to outputs/)
#   - MONTHLY_FLUXES.csv, ANNUAL_SUMMARY.csv, SI_TABLES.csv (to outputs/tables/)
#   - diagnostic plots (to outputs/figures/)
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ranger)
library(mgcv)
library(broom)
library(ggplot2)
library(purrr)
library(tibble)

cat("=== RANDOM FOREST CH4 FLUX WORKFLOW - SPECIES-FIRST VERSION (FIXED) ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD INPUT DATA
# =============================================================================

cat("Loading prepared data...\n")
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

# Extract datasets
TREE_JULY <- rf_workflow_data$PLACEHOLDER_TREE_JULY
TREE_YEAR <- rf_workflow_data$PLACEHOLDER_TREE_YEAR
SOIL_YEAR <- rf_workflow_data$PLACEHOLDER_SOIL_YEAR
INVENTORY <- rf_workflow_data$PLACEHOLDER_INVENTORY
MOISTURE_DEC_RASTER <- rf_workflow_data$PLACEHOLDER_MOISTURE_DEC_RASTER
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
TAXONOMY <- rf_workflow_data$PLACEHOLDER_TAXONOMY
PLOT_AREA <- rf_workflow_data$PLACEHOLDER_PLOT_AREA

cat("Note: Using GPS coordinates throughout (x = longitude, y = latitude)\n")

# =============================================================================
# DBH CORRECTION FOR INVENTORY
# =============================================================================

cat("\n*** APPLYING DBH CORRECTIONS ***\n")

INVENTORY <- INVENTORY %>%
  mutate(
    dbh_original = dbh_m,
    dbh_corrected = case_when(
      dbh_m > 3 ~ dbh_m / 100,
      dbh_m >= 0.5 & dbh_m <= 3 ~ dbh_m,
      TRUE ~ dbh_m
    ),
    dbh_final = case_when(
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 100 ~ dbh_corrected / 100,
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 10 ~ dbh_corrected / 10,
      grepl("Betula", species) & (dbh_corrected * 100) > 200 ~ dbh_corrected / 10,
      species == "Pinus strobus" & (dbh_corrected * 100) > 230 ~ dbh_corrected / 10,
      TRUE ~ dbh_corrected
    ),
    dbh_m = dbh_final
  ) %>%
  dplyr::select(tree_id, species, dbh_m, x, y)

n_dbh_corrected <- sum(INVENTORY$dbh_m != rf_workflow_data$PLACEHOLDER_INVENTORY$dbh_m)
cat("  DBH corrections applied to", n_dbh_corrected, "trees\n")
cat("  DBH range (cm):", round(range(INVENTORY$dbh_m * 100), 1), "\n")

# =============================================================================
# UNIT CONVERSION: nmol TO μmol
# =============================================================================

cat("\n*** APPLYING UNIT CONVERSION: nmol to μmol (÷1000) ***\n")

TREE_JULY$stem_flux_umol_m2_s <- TREE_JULY$stem_flux_umol_m2_s / 1000
TREE_YEAR$stem_flux_umol_m2_s <- TREE_YEAR$stem_flux_umol_m2_s / 1000
SOIL_YEAR$soil_flux_umol_m2_s <- SOIL_YEAR$soil_flux_umol_m2_s / 1000

cat("Updated flux ranges (μmol m-2 s-1):\n")
cat("  TREE_JULY: [", round(min(TREE_JULY$stem_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(TREE_JULY$stem_flux_umol_m2_s, na.rm=T), 6), "]\n")
cat("  TREE_YEAR: [", round(min(TREE_YEAR$stem_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(TREE_YEAR$stem_flux_umol_m2_s, na.rm=T), 6), "]\n")
cat("  SOIL_YEAR: [", round(min(SOIL_YEAR$soil_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(SOIL_YEAR$soil_flux_umol_m2_s, na.rm=T), 6), "]\n")

TREE_YEAR <- TREE_YEAR %>%
  left_join(INVENTORY %>% dplyr::select(tree_id, species), by = "tree_id")

cat("\n✓ Data loaded successfully\n\n")

# =============================================================================
# SMART IMPUTATION FUNCTION
# =============================================================================

smart_impute <- function(df, feature_cols) {
  df_imputed <- df
  
  for(col in feature_cols) {
    if(col %in% names(df) && sum(is.na(df[[col]])) > 0) {
      if(col %in% c("air_temp_C", "soil_temp_C", "soil_moisture_abs")) {
        if("month" %in% names(df)) {
          monthly_medians <- df %>%
            group_by(month) %>%
            summarise(med_val = median(!!sym(col), na.rm = TRUE), .groups = "drop")
          
          for(m in unique(df$month)) {
            month_idx <- df$month == m & is.na(df[[col]])
            if(any(month_idx)) {
              month_med <- monthly_medians$med_val[monthly_medians$month == m]
              if(!is.na(month_med)) {
                df_imputed[[col]][month_idx] <- month_med
              }
            }
          }
        }
      }
      
      remaining_na <- is.na(df_imputed[[col]])
      if(any(remaining_na)) {
        overall_median <- median(df[[col]], na.rm = TRUE)
        if(!is.na(overall_median)) {
          df_imputed[[col]][remaining_na] <- overall_median
        }
      }
    }
  }
  
  return(df_imputed)
}

# =============================================================================
# MONTHLY MOISTURE CALIBRATION
# =============================================================================

cat("Building monthly moisture affine calibration...\n")

fit_moisture_affine <- function(month_val, points_df, moisture_lookup_fn) {
  month_data <- points_df %>%
    filter(month == month_val, 
           !is.na(soil_moisture_abs),
           !is.na(x), !is.na(y))
  
  if(nrow(month_data) < 5) {
    return(list(alpha = NA, beta = NA, R2 = NA, n = nrow(month_data)))
  }
  
  month_data$moisture_dec <- moisture_lookup_fn(month_data$x, month_data$y)
  month_data <- month_data %>% filter(!is.na(moisture_dec))
  
  if(nrow(month_data) < 5) {
    return(list(alpha = NA, beta = NA, R2 = NA, n = nrow(month_data)))
  }
  
  fit <- lm(soil_moisture_abs ~ moisture_dec, data = month_data)
  
  return(list(
    alpha = coef(fit)[1],
    beta = coef(fit)[2],
    R2 = summary(fit)$r.squared,
    n = nrow(month_data)
  ))
}

all_moisture_points <- bind_rows(
  TREE_JULY %>% dplyr::select(month, x, y, soil_moisture_abs),
  TREE_YEAR %>% dplyr::select(month, x, y, soil_moisture_abs),
  SOIL_YEAR %>% dplyr::select(month, x, y, soil_moisture_abs)
) %>%
  filter(!is.na(soil_moisture_abs))

MOISTURE_AFFINE_TABLE <- map_df(1:12, function(m) {
  result <- fit_moisture_affine(m, all_moisture_points, moisture_lookup_xy)
  tibble(
    month = m,
    alpha_t = result$alpha,
    beta_t = result$beta,
    R2_t = result$R2,
    n_points = result$n
  )
})

make_Mhat <- function(affine_table, moisture_lookup_fn) {
  function(month_vec, x_vec, y_vec) {
    if(length(month_vec) == 0) return(numeric(0))
    
    moisture_pred <- numeric(length(month_vec))
    
    months_with_calib <- affine_table %>% 
      filter(!is.na(alpha_t)) %>% 
      pull(month)
    
    if(length(months_with_calib) > 0) {
      fallback_moisture <- mean(affine_table$alpha_t[!is.na(affine_table$alpha_t)])
    } else {
      fallback_moisture <- 0.15
    }
    
    unique_months <- unique(month_vec)
    for(m in unique_months) {
      idx <- which(month_vec == m)
      params <- affine_table %>% filter(month == m)
      
      if(nrow(params) == 0 || is.na(params$alpha_t)) {
        moisture_pred[idx] <- fallback_moisture
      } else {
        moisture_dec <- moisture_lookup_fn(x_vec[idx], y_vec[idx])
        pred <- params$alpha_t + params$beta_t * moisture_dec
        pred[pred < 0] <- 0
        pred[pred > 0.6] <- 0.75
        pred[is.na(pred)] <- fallback_moisture
        moisture_pred[idx] <- pred
      }
    }
    
    return(moisture_pred)
  }
}

Mhat <- make_Mhat(MOISTURE_AFFINE_TABLE, moisture_lookup_xy)
cat("✓ Monthly moisture calibration complete\n\n")

# =============================================================================
# DATA PREPARATION WITH OUTLIER REMOVAL
# =============================================================================

cat("Preparing data with chamber type as predictor...\n")

tree_combined <- bind_rows(
  TREE_JULY %>% mutate(chamber_type = "rigid"),
  TREE_YEAR %>% mutate(chamber_type = coalesce(chamber_type, "semirigid"))
) %>%
  filter(!is.na(stem_flux_umol_m2_s))

tree_combined <- tree_combined %>%
  mutate(
    dbh_original = dbh_m,
    dbh_corrected = case_when(
      dbh_m > 3 ~ dbh_m / 100,
      dbh_m >= 0.5 & dbh_m <= 3 ~ dbh_m,
      TRUE ~ dbh_m
    ),
    dbh_final = case_when(
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 100 ~ dbh_corrected / 100,
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 10 ~ dbh_corrected / 10,
      grepl("Betula", species) & (dbh_corrected * 100) > 200 ~ dbh_corrected / 10,
      species == "Pinus strobus" & (dbh_corrected * 100) > 230 ~ dbh_corrected / 10,
      TRUE ~ dbh_corrected
    ),
    dbh_m = dbh_final
  ) %>%
  dplyr::select(-dbh_original, -dbh_corrected, -dbh_final)

tree_combined <- smart_impute(tree_combined, 
                              c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m"))

# Outlier removal for trees
tree_combined$z <- asinh(tree_combined$stem_flux_umol_m2_s)
tree_percentiles <- quantile(tree_combined$z, probs = c(0.01, 0.99), na.rm = TRUE)
outlier_mask_tree <- tree_combined$z >= tree_percentiles[1] & 
  tree_combined$z <= tree_percentiles[2]
n_outliers_tree <- sum(!outlier_mask_tree)

cat("  Removing", n_outliers_tree, "tree flux outliers (", 
    round(100*n_outliers_tree/nrow(tree_combined), 1), "%)\n")

tree_combined <- tree_combined[outlier_mask_tree, ]
tree_combined$stem_flux_corrected <- tree_combined$stem_flux_umol_m2_s 
tree_combined$z_corr <- asinh(tree_combined$stem_flux_corrected)

# Soil data outlier removal
remove_outliers_mad <- function(x, k = 8) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  
  if(mad_val == 0 || is.na(mad_val)) {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    if(iqr == 0) {
      return(!is.na(x))
    }
    
    lower <- q1 - 3 * iqr
    upper <- q3 + 3 * iqr
    return(x >= lower & x <= upper)
  } else {
    lower <- med - k * mad_val
    upper <- med + k * mad_val
    return(x >= lower & x <= upper)
  }
}

SOIL_YEAR <- smart_impute(SOIL_YEAR, 
                          c("air_temp_C", "soil_temp_C", "soil_moisture_abs"))

SOIL_YEAR$z_soil <- asinh(SOIL_YEAR$soil_flux_umol_m2_s)
outlier_mask_soil <- remove_outliers_mad(SOIL_YEAR$z_soil, k = 8)
n_outliers_soil <- sum(!outlier_mask_soil)

cat("  Removing", n_outliers_soil, "soil flux outliers (", 
    round(100*n_outliers_soil/nrow(SOIL_YEAR), 1), "%)\n")

SOIL_YEAR <- SOIL_YEAR[outlier_mask_soil, ]
SOIL_YEAR$z_soil <- NULL

cat("✓ Outlier removal complete\n\n")

# =============================================================================
# EMPIRICAL SEASONAL INDEX
# =============================================================================

cat("Computing empirical seasonal indices...\n")

compute_empirical_SI <- function(df, flux_col, group_name, use_soil_temp = FALSE) {
  df_clean <- df %>%
    filter(!is.na(!!sym(flux_col))) %>%
    mutate(y_asinh = asinh(!!sym(flux_col)))
  
  if(group_name == "tree") {
    X <- df_clean %>%
      dplyr::select(any_of(c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m")))
  } else {
    X <- df_clean %>%
      dplyr::select(any_of(c("soil_temp_C", "air_temp_C", "soil_moisture_abs")))
  }
  
  X <- X[, colSums(is.na(X)) < nrow(X)*0.7, drop = FALSE]
  
  if(ncol(X) > 0 && nrow(X) > 20) {
    rf0 <- ranger(
      x = X,
      y = df_clean$y_asinh,
      num.trees = 200,
      min.node.size = 5,
      mtry = max(1, floor(sqrt(ncol(X)))),
      num.threads = 1,
      keep.inbag = TRUE
    )
    yhat_asinh <- rf0$predictions
  } else {
    yhat_asinh <- rep(mean(df_clean$y_asinh), nrow(df_clean))
  }
  
  df_clean$residual <- df_clean$y_asinh - yhat_asinh
  
  monthly_resid <- df_clean %>%
    group_by(month) %>%
    summarise(
      mean_resid = mean(residual, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  all_months <- tibble(month = 1:12)
  monthly_resid <- all_months %>%
    left_join(monthly_resid, by = "month") %>%
    mutate(mean_resid = ifelse(is.na(mean_resid), 0, mean_resid))
  
  if(sum(!is.na(monthly_resid$mean_resid)) >= 4) {
    extended <- bind_rows(
      monthly_resid %>% mutate(month = month - 12),
      monthly_resid,
      monthly_resid %>% mutate(month = month + 12)
    )
    
    smooth_fit <- smooth.spline(extended$month, extended$mean_resid, df = 4)
    SI_smooth <- predict(smooth_fit, 1:12)$y
  } else {
    SI_smooth <- monthly_resid$mean_resid
  }
  
  return(tibble(
    month = 1:12,
    SI = SI_smooth,
    group = group_name
  ))
}

tree_for_SI <- tree_combined %>%
  dplyr::select(month, stem_flux_corrected, air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m) %>%
  rename(flux = stem_flux_corrected)

SI_tree <- compute_empirical_SI(tree_for_SI, "flux", "tree")

soil_for_SI <- SOIL_YEAR %>%
  rename(flux = soil_flux_umol_m2_s)

SI_soil <- compute_empirical_SI(soil_for_SI, "flux", "soil", use_soil_temp = TRUE)

SI_TABLES <- bind_rows(SI_tree, SI_soil)
cat("✓ Seasonal indices computed\n\n")

# =============================================================================
# FEATURE ENGINEERING WITH SPECIES-FIRST APPROACH
# =============================================================================

cat("Building features for Random Forest models (species-first approach)...\n")

compute_taxon_prior <- function(tree_df, taxonomy_df) {
  df_clean <- tree_df %>%
    filter(!is.na(stem_flux_corrected)) %>%
    mutate(y_asinh = asinh(stem_flux_corrected))
  
  X_env <- df_clean %>%
    dplyr::select(any_of(c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m")))
  
  rf_env <- ranger(
    x = X_env,
    y = df_clean$y_asinh,
    num.trees = 200,
    min.node.size = 5,
    num.threads = 1
  )
  
  df_clean$residual <- df_clean$y_asinh - rf_env$predictions
  
  df_clean <- df_clean %>%
    left_join(taxonomy_df %>% dplyr::select(species, genus, family, order, class),
              by = "species")
  
  med_species <- df_clean %>%
    group_by(species) %>%
    summarise(
      med_resid = ifelse(n() >= 5, median(residual, na.rm = TRUE), NA_real_),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(med_resid))
  
  med_genus <- df_clean %>%
    group_by(genus) %>%
    summarise(
      med_resid = ifelse(n() >= 5, median(residual, na.rm = TRUE), NA_real_),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(med_resid))
  
  med_family <- df_clean %>%
    group_by(family) %>%
    summarise(
      med_resid = ifelse(n() >= 5, median(residual, na.rm = TRUE), NA_real_),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(med_resid))
  
  return(list(
    species = med_species,
    genus = med_genus,
    family = med_family
  ))
}

TAXONOMY_PRIORS <- compute_taxon_prior(tree_combined, TAXONOMY)
cat("✓ Taxonomy priors computed\n")

build_features_tree <- function(df, drivers, Mhat_fn, SI_table, taxonomy, taxon_priors) {
  df <- df %>%
    left_join(drivers %>% dplyr::select(month, air_temp_C_mean, soil_temp_C_mean),
              by = "month")
  
  si_tree <- SI_table %>% filter(group == "tree") %>% dplyr::select(month, SI)
  df <- df %>%
    left_join(si_tree, by = "month")
  
  df$air_temp_C_mean <- as.numeric(df$air_temp_C_mean)
  df$soil_temp_C_mean <- as.numeric(df$soil_temp_C_mean)
  df$SI <- as.numeric(df$SI)
  df$dbh_m <- as.numeric(df$dbh_m)
  
  df$month_sin <- sin(2 * pi * df$month / 12)
  df$month_cos <- cos(2 * pi * df$month / 12)
  
  if("soil_moisture_abs" %in% names(df)) {
    df$soil_moisture_at_tree <- as.numeric(df$soil_moisture_abs)
  } else {
    df$soil_moisture_at_tree <- NA_real_
  }
  
  missing_moisture <- is.na(df$soil_moisture_at_tree)
  if(any(missing_moisture) && !is.null(Mhat_fn)) {
    df$soil_moisture_at_tree[missing_moisture] <- 
      Mhat_fn(df$month[missing_moisture], df$x[missing_moisture], df$y[missing_moisture])
  }
  df$soil_moisture_at_tree <- as.numeric(df$soil_moisture_at_tree)
  
  if("chamber_type" %in% names(df)) {
    df$chamber_rigid <- as.numeric(df$chamber_type == "rigid")
    df$chamber_semirigid <- as.numeric(df$chamber_type == "semirigid")
  } else {
    df$chamber_rigid <- 1
    df$chamber_semirigid <- 0
  }
  
  df <- df %>%
    left_join(taxonomy %>% dplyr::select(species, genus, family), by = "species")
  
  df$taxon_prior_asinh <- 0
  
  for(i in 1:nrow(df)) {
    if(!is.na(df$species[i]) && df$species[i] %in% taxon_priors$species$species) {
      matches <- which(taxon_priors$species$species == df$species[i])
      if(length(matches) == 1) {
        df$taxon_prior_asinh[i] <- taxon_priors$species$med_resid[matches]
      }
    } else if(!is.na(df$genus[i]) && df$genus[i] %in% taxon_priors$genus$genus) {
      matches <- which(taxon_priors$genus$genus == df$genus[i])
      if(length(matches) == 1) {
        df$taxon_prior_asinh[i] <- taxon_priors$genus$med_resid[matches]
      }
    } else if(!is.na(df$family[i]) && df$family[i] %in% taxon_priors$family$family) {
      matches <- which(taxon_priors$family$family == df$family[i])
      if(length(matches) == 1) {
        df$taxon_prior_asinh[i] <- taxon_priors$family$med_resid[matches]
      }
    }
  }
  df$taxon_prior_asinh <- as.numeric(df$taxon_prior_asinh)
  
  species_counts <- table(df$species)
  df$species_clean <- ifelse(species_counts[df$species] >= 10, 
                             df$species, "SPECIES_OTHER")
  
  # CRITICAL: Within-species DBH standardization
  df$dbh_within_z <- ave(
    df$dbh_m, 
    df$species_clean,
    FUN = function(x) {
      if(length(x) >= 5 && sd(x, na.rm=TRUE) > 0) {
        (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
      } else {
        rep(0, length(x))
      }
    }
  )
  
  df$species_factor <- factor(
    df$species_clean,
    levels = c("SPECIES_OTHER", names(sort(table(df$species_clean)[-which(names(table(df$species_clean)) == "SPECIES_OTHER")], 
                                           decreasing = TRUE)))
  )
  
  if("stem_flux_corrected" %in% names(df) && !all(is.na(df$stem_flux_corrected))) {
    df$y_asinh <- asinh(as.numeric(df$stem_flux_corrected))
  } else if("stem_flux_umol_m2_s" %in% names(df) && !all(is.na(df$stem_flux_umol_m2_s))) {
    df$y_asinh <- asinh(as.numeric(df$stem_flux_umol_m2_s))
  } else {
    df$y_asinh <- NA_real_
  }
  
  return(df)
}

build_features_soil <- function(df, drivers, Mhat_fn, SI_table) {
  df <- df %>%
    left_join(drivers %>% dplyr::select(month, air_temp_C_mean, soil_temp_C_mean),
              by = "month")
  
  si_soil <- SI_table %>% filter(group == "soil") %>% dplyr::select(month, SI)
  df <- df %>%
    left_join(si_soil, by = "month")
  
  df$air_temp_C_mean <- as.numeric(df$air_temp_C_mean)
  df$soil_temp_C_mean <- as.numeric(df$soil_temp_C_mean)
  df$SI <- as.numeric(df$SI)
  
  df$month_sin <- sin(2 * pi * df$month / 12)
  df$month_cos <- cos(2 * pi * df$month / 12)
  
  if("soil_moisture_abs" %in% names(df)) {
    df$soil_moisture_at_site <- as.numeric(df$soil_moisture_abs)
  } else {
    df$soil_moisture_at_site <- NA_real_
  }
  
  missing_moisture <- is.na(df$soil_moisture_at_site)
  if(any(missing_moisture) && !is.null(Mhat_fn)) {
    df$soil_moisture_at_site[missing_moisture] <- 
      Mhat_fn(df$month[missing_moisture], df$x[missing_moisture], df$y[missing_moisture])
  }
  df$soil_moisture_at_site <- as.numeric(df$soil_moisture_at_site)
  
  df$moisture_x_soilT <- as.numeric(df$soil_moisture_at_site) * as.numeric(df$soil_temp_C_mean)
  df$moisture_x_airT <- as.numeric(df$soil_moisture_at_site) * as.numeric(df$air_temp_C_mean)
  
  if("soil_flux_umol_m2_s" %in% names(df)) {
    df$y_asinh <- asinh(as.numeric(df$soil_flux_umol_m2_s))
  }
  
  return(df)
}

tree_train <- build_features_tree(tree_combined, DRIVERS, Mhat, SI_TABLES, 
                                  TAXONOMY, TAXONOMY_PRIORS)
soil_train <- build_features_soil(SOIL_YEAR, DRIVERS, Mhat, SI_TABLES)

cat("✓ Features built for training\n\n")

# =============================================================================
# SECTION 5: SPECIES-FIRST RANDOM FOREST MODELS (FULLY FIXED)
# =============================================================================

cat("Training Species-First Random Forest models...\n")

# First identify complete rows - need ALL features to be complete
complete_rows <- !is.na(tree_train$y_asinh) & 
  !is.na(tree_train$soil_moisture_at_tree) & 
  !is.na(tree_train$species_factor)

# Subset the training data to complete cases FIRST
tree_train_complete <- tree_train[complete_rows, ]
y_tree <- tree_train_complete$y_asinh

cat("  Complete cases for training:", nrow(tree_train_complete), "out of", nrow(tree_train), "\n")

# Now build matrices from the complete data
species_dummies_full <- model.matrix(~ species_factor - 1, data = tree_train_complete)

# Create species×moisture interactions manually (more reliable)
species_moisture_full <- species_dummies_full * tree_train_complete$soil_moisture_at_tree
colnames(species_moisture_full) <- paste0(colnames(species_dummies_full), ".soil_moisture_at_tree")

# Verify dimensions
cat("  Species dummies:", nrow(species_dummies_full), "rows,", ncol(species_dummies_full), "columns\n")
cat("  Species×moisture:", nrow(species_moisture_full), "rows,", ncol(species_moisture_full), "columns\n")

# Build feature matrix
X_tree_species_first <- data.frame(
  species_dummies_full,
  species_moisture_full,
  dbh_within_z = tree_train_complete$dbh_within_z,
  air_temp_C_mean = tree_train_complete$air_temp_C_mean,
  soil_temp_C_mean = tree_train_complete$soil_temp_C_mean,
  soil_moisture_at_tree = tree_train_complete$soil_moisture_at_tree,
  SI = tree_train_complete$SI,
  taxon_prior_asinh = tree_train_complete$taxon_prior_asinh,
  chamber_rigid = tree_train_complete$chamber_rigid,
  chamber_semirigid = tree_train_complete$chamber_semirigid,
  month_sin = tree_train_complete$month_sin,
  month_cos = tree_train_complete$month_cos
)

# Remove any columns with zero variance
zero_var_cols <- sapply(X_tree_species_first, function(x) var(x, na.rm=TRUE) == 0)
if(any(zero_var_cols)) {
  cat("  Removing", sum(zero_var_cols), "zero-variance columns\n")
  X_tree_species_first <- X_tree_species_first[, !zero_var_cols]
}

# Store for later use
X_tree <- X_tree_species_first
X_tree_template <- X_tree

cat("  Tree feature matrix:", nrow(X_tree), "obs x", ncol(X_tree), "features\n")
cat("  Using within-species DBH standardization (deconfounded)\n")
cat("  Species×moisture interactions for ALL species including OTHER\n")

# Train TreeRF
TreeRF <- ranger(
  x = as.data.frame(X_tree),
  y = y_tree,
  num.trees = 800,
  min.node.size = 5,
  mtry = floor(sqrt(ncol(X_tree))),
  importance = "impurity",
  num.threads = 1,
  oob.error = TRUE,
  seed = 42
)

cat("\nTreeRF (Species-First) trained:\n")
cat("  OOB R²:", round(TreeRF$r.squared, 3), "\n")
cat("  OOB RMSE (asinh):", round(sqrt(TreeRF$prediction.error), 4), "\n")

# Soil model
soil_features <- c("soil_temp_C_mean", "air_temp_C_mean",
                   "soil_moisture_at_site", "SI",
                   "moisture_x_soilT", "moisture_x_airT",
                   "month_sin", "month_cos")

X_soil <- soil_train[, soil_features, drop = FALSE]

zero_var_soil <- sapply(X_soil, function(x) {
  if(all(is.na(x))) return(TRUE)
  var(x, na.rm = TRUE) == 0
})
X_soil <- X_soil[, !zero_var_soil, drop = FALSE]

y_soil <- soil_train$y_asinh

complete_rows_soil <- !is.na(y_soil)
X_soil <- X_soil[complete_rows_soil, , drop = FALSE]
y_soil <- y_soil[complete_rows_soil]

cat("  Soil feature matrix:", nrow(X_soil), "obs x", ncol(X_soil), "features\n")

X_soil_template <- X_soil

SoilRF <- ranger(
  x = as.data.frame(X_soil),
  y = y_soil,
  num.trees = 800,
  min.node.size = 5,
  mtry = floor(sqrt(ncol(X_soil))),
  importance = "impurity",
  num.threads = 1,
  oob.error = TRUE,
  seed = 42
)

cat("\nSoilRF trained:\n")
cat("  OOB R²:", round(SoilRF$r.squared, 3), "\n")
cat("  OOB RMSE (asinh):", round(sqrt(SoilRF$prediction.error), 4), "\n\n")

## Save training data + feature matrices early (before QC, which may error)
## Add predictions to tree_train_complete for downstream plotting
tree_train_complete$pred_asinh <- TreeRF$predictions
tree_train_complete$pred_flux <- sinh(tree_train_complete$pred_asinh)
tree_train_complete$pred_flux_nmol <- tree_train_complete$pred_flux * 1000
tree_train_complete$obs_flux_nmol <- tree_train_complete$stem_flux_corrected * 1000

## Add predictions to soil data
soil_train_complete <- soil_train[complete_rows_soil, ]
soil_train_complete$pred_asinh <- SoilRF$predictions
soil_train_complete$pred_flux <- sinh(soil_train_complete$pred_asinh)

save(tree_train_complete, X_tree, X_soil,
     soil_train_complete, complete_rows_soil,
     file = "../../outputs/models/TRAINING_DATA.RData")
cat("  Saved training data to outputs/models/TRAINING_DATA.RData\n")

# =============================================================================
# GEOMETRY CALCULATIONS
# =============================================================================

cat("Computing plot geometry...\n")

compute_geometry <- function(inventory_df) {
  inventory_df$S_i <- ifelse(
    inventory_df$species == "Kalmia latifolia",
    pi * inventory_df$dbh_m * 0.75,
    pi * inventory_df$dbh_m * 2
  )
  
  inventory_df$BA_i <- pi * (inventory_df$dbh_m/2)^2
  total_BA <- sum(inventory_df$BA_i)
  A_soil <- PLOT_AREA - total_BA
  
  return(list(
    inventory = inventory_df,
    A_soil = A_soil,
    A_plot = PLOT_AREA
  ))
}

geometry <- compute_geometry(INVENTORY)
cat("  Total basal area:", round(sum(geometry$inventory$BA_i), 1), "m²\n")
cat("  Soil area:", round(geometry$A_soil, 1), "m²\n")
cat("  Soil fraction:", round(geometry$A_soil/geometry$A_plot, 3), "\n\n")

# =============================================================================
# MONTHLY SPATIAL PREDICTIONS (FIXED)
# =============================================================================

cat("\n=== MONTHLY SPATIAL PREDICTIONS (SPECIES-FIRST) ===\n")

# Prepare INVENTORY with within-species DBH
species_counts_train <- table(tree_train$species)
INVENTORY$species_clean <- ifelse(
  INVENTORY$species %in% names(species_counts_train) & species_counts_train[INVENTORY$species] >= 10,
  INVENTORY$species,
  "SPECIES_OTHER"
)

dbh_stats <- tree_train %>%
  group_by(species_clean) %>%
  summarise(
    mean_dbh = mean(dbh_m, na.rm=TRUE),
    sd_dbh = sd(dbh_m, na.rm=TRUE),
    n = n(),
    .groups = "drop"
  )

INVENTORY <- INVENTORY %>%
  left_join(dbh_stats, by = "species_clean") %>%
  mutate(
    dbh_within_z = case_when(
      n >= 5 & sd_dbh > 0 ~ (dbh_m - mean_dbh) / sd_dbh,
      TRUE ~ 0
    )
  )

training_ranges <- list(
  air_temp = quantile(c(tree_combined$air_temp_C, SOIL_YEAR$air_temp_C), 
                      probs = c(0.01, 0.99), na.rm = TRUE),
  soil_temp = quantile(c(tree_combined$soil_temp_C, SOIL_YEAR$soil_temp_C), 
                       probs = c(0.01, 0.99), na.rm = TRUE),
  moisture = quantile(c(tree_combined$soil_moisture_abs, SOIL_YEAR$soil_moisture_abs), 
                      probs = c(0.01, 0.99), na.rm = TRUE),
  dbh = quantile(INVENTORY$dbh_m, probs = c(0.01, 0.99), na.rm = TRUE)
)

clip_to_range <- function(x, range_vals) {
  if (any(is.na(range_vals))) return(x)
  pmax(pmin(x, range_vals[2]), range_vals[1])
}

monthly_results <- tibble(month = 1:12)
tree_fluxes_monthly <- list()

for (t in 1:12) {
  cat("  Month", t, "...")
  
  air_temp_t <- DRIVERS$air_temp_C_mean[t]
  soil_temp_t <- DRIVERS$soil_temp_C_mean[t]
  
  if (is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  SI_tree_row <- SI_TABLES %>% 
    filter(group == "tree", month == t)
  SI_tree_t <- if(nrow(SI_tree_row) > 0) SI_tree_row$SI[1] else 0
  
  inv_predictions <- INVENTORY %>%
    mutate(
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = clip_to_range(
        MOISTURE_AFFINE_TABLE$alpha_t[t] + 
          MOISTURE_AFFINE_TABLE$beta_t[t] * moisture_raw,
        training_ranges$moisture
      ),
      air_temp_C = clip_to_range(air_temp_t, training_ranges$air_temp),
      soil_temp_C = clip_to_range(soil_temp_t, training_ranges$soil_temp),
      SI_tree = SI_tree_t,
      month_sin = sin(2 * pi * t / 12),
      month_cos = cos(2 * pi * t / 12),
      chamber_rigid = 1,
      chamber_semirigid = 0,
      species_factor = factor(
        species_clean,
        levels = levels(tree_train$species_factor)
      ),
      taxon_prior_asinh = 0
    )
  
  # Build prediction features (FIXED)
  species_dummies_pred <- model.matrix(~ species_factor - 1, data = inv_predictions)
  
  # Create species-moisture interactions manually
  species_moisture_pred <- species_dummies_pred * inv_predictions$soil_moisture_abs
  colnames(species_moisture_pred) <- paste0(colnames(species_dummies_pred), ".soil_moisture_at_tree")
  
  X_pred_tree <- data.frame(
    species_dummies_pred,
    species_moisture_pred,
    dbh_within_z = inv_predictions$dbh_within_z,
    air_temp_C_mean = inv_predictions$air_temp_C,
    soil_temp_C_mean = inv_predictions$soil_temp_C,
    soil_moisture_at_tree = inv_predictions$soil_moisture_abs,
    SI = inv_predictions$SI_tree,
    taxon_prior_asinh = inv_predictions$taxon_prior_asinh,
    chamber_rigid = inv_predictions$chamber_rigid,
    chamber_semirigid = inv_predictions$chamber_semirigid,
    month_sin = inv_predictions$month_sin,
    month_cos = inv_predictions$month_cos
  )
  
  # Align columns with training
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred_tree), ncol = ncol(X_tree))
  colnames(X_pred_aligned) <- colnames(X_tree)
  
  for (col_name in colnames(X_tree)) {
    if (col_name %in% colnames(X_pred_tree)) {
      X_pred_aligned[, col_name] <- X_pred_tree[[col_name]]
    }
  }
  
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  pred_flux_umol_m2_s <- sinh(pred_asinh)
  
  S_tree <- ifelse(
    inv_predictions$species == "Kalmia latifolia",
    pi * inv_predictions$dbh_m * 0.75,
    pi * inv_predictions$dbh_m * 2
  )
  
  R_tree <- pred_flux_umol_m2_s * S_tree
  
  tree_fluxes_monthly[[t]] <- list(
    month = t,
    total_R_tree = sum(R_tree, na.rm = TRUE),
    mean_flux = mean(pred_flux_umol_m2_s, na.rm = TRUE),
    median_flux = median(pred_flux_umol_m2_s, na.rm = TRUE),
    n_trees = sum(!is.na(pred_flux_umol_m2_s))
  )
  
  cat(" done. Total emission:", round(sum(R_tree, na.rm = TRUE), 6), "μmol/s\n")
}

tree_results <- bind_rows(tree_fluxes_monthly) %>%
  mutate(
    Phi_tree_umol_m2_s = total_R_tree / geometry$A_plot
  )

cat("  Tree predictions complete.\n")

# Soil predictions
cat("\nPredicting soil fluxes on spatial grid...\n")

x_range <- range(INVENTORY$x, na.rm = TRUE)
y_range <- range(INVENTORY$y, na.rm = TRUE)
grid_res <- 50
x_grid <- seq(x_range[1], x_range[2], length.out = grid_res)
y_grid <- seq(y_range[1], y_range[2], length.out = grid_res)
soil_grid <- expand.grid(x = x_grid, y = y_grid)

lon_dist_m <- diff(x_grid)[1] * 111320 * cos(mean(y_grid) * pi/180)
lat_dist_m <- diff(y_grid)[1] * 111320
pixel_area <- lon_dist_m * lat_dist_m

cat("  Created", nrow(soil_grid), "grid points\n")
cat("  Pixel area:", round(pixel_area, 2), "m²\n")

soil_fluxes_monthly <- list()

for (t in 1:12) {
  cat("  Month", t, "...")
  
  air_temp_t <- DRIVERS$air_temp_C_mean[t]
  soil_temp_t <- DRIVERS$soil_temp_C_mean[t]
  
  if (is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  SI_soil_row <- SI_TABLES %>% 
    filter(group == "soil", month == t)
  SI_soil_t <- if(nrow(SI_soil_row) > 0) SI_soil_row$SI[1] else 0
  
  grid_predictions <- soil_grid %>%
    mutate(
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = clip_to_range(
        MOISTURE_AFFINE_TABLE$alpha_t[t] + 
          MOISTURE_AFFINE_TABLE$beta_t[t] * moisture_raw,
        training_ranges$moisture
      ),
      air_temp_C = clip_to_range(air_temp_t, training_ranges$air_temp),
      soil_temp_C = clip_to_range(soil_temp_t, training_ranges$soil_temp),
      SI_soil = SI_soil_t,
      month_sin = sin(2 * pi * t / 12),
      month_cos = cos(2 * pi * t / 12),
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      moisture_x_airT = soil_moisture_abs * air_temp_C
    )
  
  X_pred_soil <- grid_predictions %>%
    dplyr::select(soil_temp_C, soil_moisture_abs, SI_soil, 
                  moisture_x_soilT, air_temp_C, moisture_x_airT,
                  month_sin, month_cos) %>%
    as.matrix()
  
  X_pred_soil_aligned <- matrix(0, nrow = nrow(X_pred_soil), ncol = ncol(X_soil))
  colnames(X_pred_soil_aligned) <- colnames(X_soil)
  
  for (col_name in colnames(X_soil)) {
    if (col_name %in% colnames(X_pred_soil)) {
      X_pred_soil_aligned[, col_name] <- X_pred_soil[, col_name]
    }
  }
  
  pred_asinh_soil <- predict(SoilRF, X_pred_soil_aligned)$predictions
  pred_flux_soil <- sinh(pred_asinh_soil)
  mean_soil_flux <- mean(pred_flux_soil, na.rm = TRUE)
  
  soil_fluxes_monthly[[t]] <- list(
    month = t,
    mean_flux = mean_soil_flux,
    median_flux = median(pred_flux_soil, na.rm = TRUE),
    n_valid = sum(!is.na(pred_flux_soil))
  )
  
  cat(" done. Mean flux:", round(mean_soil_flux, 6), "μmol/m²/s\n")
}

soil_results <- bind_rows(soil_fluxes_monthly) %>%
  mutate(
    Phi_soil_umol_m2_s = mean_flux * (geometry$A_soil / geometry$A_plot)
  )

cat("  Soil predictions complete.\n")

# Combine results
cat("\nCombining tree and soil predictions...\n")

monthly_results <- monthly_results %>%
  left_join(tree_results %>% dplyr::select(month, Phi_tree_umol_m2_s), by = "month") %>%
  left_join(soil_results %>% dplyr::select(month, Phi_soil_umol_m2_s), by = "month") %>%
  mutate(
    Phi_plot_umol_m2_s = Phi_tree_umol_m2_s + Phi_soil_umol_m2_s
  )

monthly_results <- monthly_results %>%
  mutate(
    Phi_tree_mg_m2_d = Phi_tree_umol_m2_s * 86400 * 16 * 1e-3,
    Phi_soil_mg_m2_d = Phi_soil_umol_m2_s * 86400 * 16 * 1e-3,
    Phi_plot_mg_m2_d = Phi_plot_umol_m2_s * 86400 * 16 * 1e-3
  )

print(monthly_results)

annual_summary <- monthly_results %>%
  summarise(
    annual_tree_mg_m2 = sum(Phi_tree_mg_m2_d) * 30.4,
    annual_soil_mg_m2 = sum(Phi_soil_mg_m2_d) * 30.4,
    annual_plot_mg_m2 = sum(Phi_plot_mg_m2_d) * 30.4
  )

cat("\nAnnual totals (mg CH4 m-2 yr-1):\n")
print(annual_summary)
cat("\nNote: Negative values indicate CH4 uptake (soil consumption)\n")

# Save outputs
cat("\nSaving outputs...\n")
write_csv(monthly_results, "../../outputs/tables/MONTHLY_FLUXES.csv")
write_csv(annual_summary, "../../outputs/tables/ANNUAL_SUMMARY.csv")
write_csv(MOISTURE_AFFINE_TABLE, "../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")
write_csv(SI_TABLES, "../../outputs/tables/SI_TABLES.csv")
save(TAXONOMY_PRIORS, file = "../../outputs/models/TAXONOMY_PRIORS.RData")
save(TreeRF, SoilRF, file = "../../outputs/models/RF_MODELS.RData")

cat("\n=== WORKFLOW COMPLETE - SPECIES-FIRST VERSION (FIXED) ===\n")
cat("Key improvements:\n")
cat("  - Within-species DBH standardization (deconfounded)\n")
cat("  - Species×moisture interactions for ALL species\n")
cat("  - Fixed matrix subsetting issues\n")
cat("  - TreeRF R² =", round(TreeRF$r.squared, 3), "\n")
cat("  - SoilRF R² =", round(SoilRF$r.squared, 3), "\n")





# =============================================================================
# SECTION 11: QUALITY CONTROL PLOTS
# =============================================================================

cat("\nGenerating quality control plots...\n")

# QC1: Chamber type comparison
qc1_data <- tree_train %>%
  mutate(
    chamber = ifelse(chamber_rigid == 1, "rigid", "semirigid"),
    flux_original = sinh(y_asinh)
  )

p1 <- ggplot(qc1_data, aes(x = chamber, y = flux_original)) +
  geom_boxplot(aes(fill = chamber)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Flux Distribution by Chamber Type",
       subtitle = "Chamber type included as RF predictor",
       x = "Chamber Type",
       y = "Flux (μmol m-2 s-1, pseudo-log scale)") +
  theme_minimal()

ggsave("../../outputs/figures/QC_CHAMBER_COMPARISON.png", p1, width = 8, height = 6)

# QC1b: Outlier diagnostics
outlier_summary <- tibble(
  dataset = c("Tree Fluxes", "Soil Fluxes"),
  n_outliers = c(n_outliers_tree, n_outliers_soil),
  pct_outliers = c(100*n_outliers_tree/(n_outliers_tree + nrow(tree_combined)),
                   100*n_outliers_soil/(n_outliers_soil + nrow(SOIL_YEAR))),
  method = "MAD (k=8)"
)

p1b <- ggplot(outlier_summary, aes(x = dataset, y = pct_outliers)) +
  geom_col(fill = "coral") +
  geom_text(aes(label = paste0(round(pct_outliers, 1), "%\n(n=", n_outliers, ")")), 
            vjust = -0.5) +
  labs(title = "Outlier Removal Summary",
       subtitle = "Outliers detected on asinh-transformed scale",
       x = "", y = "Percent of Data Removed") +
  ylim(0, max(outlier_summary$pct_outliers) * 1.2) +
  theme_minimal()

ggsave("../../outputs/figures/QC_OUTLIERS.png", p1b, width = 6, height = 6)

# QC2: Predictions vs observations - Trees
tree_train$pred_asinh <- TreeRF$predictions
tree_train$pred_flux <- sinh(tree_train$pred_asinh)
tree_train$chamber <- ifelse(tree_train$chamber_rigid == 1, "rigid", "semirigid")

p2 <- ggplot(tree_train, aes(x = stem_flux_corrected, y = pred_flux, 
                             color = chamber)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Tree Predictions vs Observations",
       subtitle = "With month as cyclic feature",
       x = "Observed Flux (μmol m-2 s-1)",
       y = "Predicted Flux (μmol m-2 s-1)",
       color = "Chamber Type") +
  scale_color_manual(values = c("rigid" = "blue", "semirigid" = "red")) +
  theme_minimal()

ggsave("../../outputs/figures/QC_PRED_VS_OBS_TREE.png", p2, width = 8, height = 6)

# QC2b: Predictions vs observations - Soil
soil_train$pred_asinh <- SoilRF$predictions
soil_train$pred_flux <- sinh(soil_train$pred_asinh)

p2b <- ggplot(soil_train, aes(x = soil_flux_umol_m2_s, y = pred_flux)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Soil Predictions vs Observations",
       subtitle = paste("R² =", round(SoilRF$r.squared, 3), "- With month as cyclic feature"),
       x = "Observed Flux (μmol m-2 s-1)",
       y = "Predicted Flux (μmol m-2 s-1)") +
  theme_minimal()

ggsave("../../outputs/figures/QC_PRED_VS_OBS_SOIL.png", p2b, width = 8, height = 6)

# QC4: Moisture calibration R² by month
p4 <- ggplot(MOISTURE_AFFINE_TABLE, aes(x = month, y = R2_t)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = ifelse(!is.na(R2_t), round(R2_t, 2), "NA")), vjust = -0.5) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Moisture Calibration Quality by Month",
       x = "Month", y = "R²") +
  ylim(0, 1) +
  theme_minimal()

ggsave("../../outputs/figures/QC_MOISTURE_CALIB.png", p4, width = 8, height = 6)

# QC5: Feature importance - ROBUST FIX
# Debug the structure first
cat("\n  Debugging importance structure:\n")
tree_imp_raw <- importance(TreeRF)
cat("    Class:", class(tree_imp_raw), "\n")
cat("    Length:", length(tree_imp_raw), "\n")
if(is.matrix(tree_imp_raw)) {
  cat("    Dimensions:", dim(tree_imp_raw), "\n")
  cat("    Column names:", colnames(tree_imp_raw), "\n")
}

# Convert to data frame properly
if(is.matrix(tree_imp_raw)) {
  tree_importance <- as.data.frame(tree_imp_raw)
} else {
  # If it's a vector, convert differently
  tree_importance <- data.frame(Importance = as.numeric(tree_imp_raw))
}
tree_importance$feature <- rownames(tree_importance)

# If there's no column named Importance yet, rename the first column
if(!"Importance" %in% names(tree_importance)) {
  names(tree_importance)[1] <- "Importance"
}

# Now sort and dplyr::select top 20
tree_importance <- tree_importance[order(tree_importance$Importance, decreasing = TRUE), ]
tree_importance <- tree_importance[1:min(20, nrow(tree_importance)), ]

p5a <- ggplot(tree_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  labs(title = "Tree Model - Top 20 Feature Importances",
       subtitle = "Now including month_sin and month_cos",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("../../outputs/figures/QC_TREE_IMPORTANCE.png", p5a, width = 10, height = 8)

# Same robust approach for soil
soil_imp_raw <- importance(SoilRF)
if(is.matrix(soil_imp_raw)) {
  soil_importance <- as.data.frame(soil_imp_raw)
} else {
  soil_importance <- data.frame(Importance = as.numeric(soil_imp_raw))
}
soil_importance$feature <- rownames(soil_importance)

if(!"Importance" %in% names(soil_importance)) {
  names(soil_importance)[1] <- "Importance"
}

soil_importance <- soil_importance[order(soil_importance$Importance, decreasing = TRUE), ]

p5b <- ggplot(soil_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Soil Model - Feature Importances",
       subtitle = "Now including month_sin and month_cos",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("../../outputs/figures/QC_SOIL_IMPORTANCE.png", p5b, width = 10, height = 6)

cat("✓ Quality control plots saved\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\nSaving outputs...\n")

# Save main results
write_csv(monthly_results, "../../outputs/tables/MONTHLY_FLUXES.csv")
write_csv(annual_summary, "../../outputs/tables/ANNUAL_SUMMARY.csv")

# Save diagnostic tables
write_csv(MOISTURE_AFFINE_TABLE, "../../outputs/tables/MOISTURE_AFFINE_TABLE.csv")
write_csv(SI_TABLES, "../../outputs/tables/SI_TABLES.csv")

# Save taxonomy priors
save(TAXONOMY_PRIORS, file = "../../outputs/models/TAXONOMY_PRIORS.RData")

# Save models
save(TreeRF, SoilRF, file = "../../outputs/models/RF_MODELS.RData")

# Diagnostics
diagnostics <- list(
  TreeRF_OOB_R2 = TreeRF$r.squared,
  TreeRF_OOB_RMSE = sqrt(TreeRF$prediction.error),
  SoilRF_OOB_R2 = SoilRF$r.squared,
  SoilRF_OOB_RMSE = sqrt(SoilRF$prediction.error),
  Chamber_as_predictor = TRUE,
  Month_as_predictor = TRUE,  # ADDED
  Chamber_beta1 = NA,  # Not applicable with chamber as predictor
  TreeRF_importance = importance(TreeRF),
  SoilRF_importance = importance(SoilRF),
  n_tree_training = nrow(X_tree),
  n_soil_training = nrow(X_soil),
  outliers_removed = list(
    tree_n = n_outliers_tree,
    tree_pct = 100*n_outliers_tree/(n_outliers_tree + nrow(tree_combined)),
    soil_n = n_outliers_soil,
    soil_pct = 100*n_outliers_soil/(n_outliers_soil + nrow(SOIL_YEAR)),
    method = "MAD with k=8"
  ),
  data_completeness = list(
    tree_complete_before_impute = n_before_impute,
    tree_complete_after_impute = n_after_impute,
    soil_complete_before_impute = n_before_soil,
    soil_complete_after_impute = n_after_soil
  ),
  dbh_corrections_applied = n_dbh_corrected,
  coordinate_system = "GPS (x = longitude, y = latitude)"  # ADDED
)

save(diagnostics, file = "../../outputs/models/DIAGNOSTICS.RData")

# Coverage report
coverage_report <- list(
  n_trees_july = length(unique(TREE_JULY$tree_id)),
  n_trees_year = length(unique(TREE_YEAR$tree_id)),
  n_soil_sites = length(unique(SOIL_YEAR$site_id)),
  n_inventory_trees = nrow(INVENTORY),
  soil_area_fraction = geometry$A_soil / geometry$A_plot,
  species_in_training = length(unique(tree_train$species)),
  species_in_inventory = length(unique(INVENTORY$species)),
  plot_area_ha = PLOT_AREA / 10000,
  coordinate_system = "GPS throughout (x = longitude, y = latitude)"  # ADDED
)

save(coverage_report, file = "../../outputs/models/COVERAGE_REPORT.RData")

cat("\n=== WORKFLOW COMPLETE - GPS COORDINATES VERSION ===\n")
cat("All spatial data uses GPS coordinates (x = longitude, y = latitude)\n")




# =============================================================================
# UPDATE PLOTS TO NMOL UNITS
# Add this code after Section 11 (Quality Control Plots) in your workflow
# Or replace the existing QC2 plots section
# =============================================================================

cat("\nGenerating updated plots with nmol units...\n")

# QC2: Predictions vs observations - Trees (in nmol)
tree_train$pred_asinh <- TreeRF$predictions
tree_train$pred_flux <- sinh(tree_train$pred_asinh)
tree_train$chamber <- ifelse(tree_train$chamber_rigid == 1, "rigid", "semirigid")

# Convert to nmol for plotting (multiply by 1000)
tree_train$pred_flux_nmol <- tree_train$pred_flux * 1000
tree_train$obs_flux_nmol <- tree_train$stem_flux_corrected * 1000

p2_nmol <- ggplot(tree_train, aes(x = obs_flux_nmol, y = pred_flux_nmol, 
                                  color = chamber)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Tree Predictions vs Observations",
       subtitle = paste("R² =", round(TreeRF$r.squared, 3), "- With month as cyclic feature"),
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Chamber Type") +
  scale_color_manual(values = c("rigid" = "blue", "semirigid" = "red")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_TREE_nmol.png", p2_nmol, width = 8, height = 6)
cat("  Saved: QC_PRED_VS_OBS_TREE_nmol.png\n")

# QC2b: Predictions vs observations - Soil (in nmol)
soil_train$pred_asinh <- SoilRF$predictions
soil_train$pred_flux <- sinh(soil_train$pred_asinh)

# Convert to nmol for plotting
soil_train$pred_flux_nmol <- soil_train$pred_flux * 1000
soil_train$obs_flux_nmol <- soil_train$soil_flux_umol_m2_s * 1000

p2b_nmol <- ggplot(soil_train, aes(x = obs_flux_nmol, y = pred_flux_nmol)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Soil Predictions vs Observations",
       subtitle = paste("R² =", round(SoilRF$r.squared, 3), "- With month as cyclic feature"),
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")"))) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_SOIL_nmol.png", p2b_nmol, width = 8, height = 6)
cat("  Saved: QC_PRED_VS_OBS_SOIL_nmol.png\n")

# Optional: Create a combined plot showing both tree and soil on same scale
combined_data <- bind_rows(
  tree_train %>% 
    dplyr::select(obs_flux_nmol, pred_flux_nmol) %>%
    mutate(type = "Trees"),
  soil_train %>%
    dplyr::select(obs_flux_nmol, pred_flux_nmol) %>%
    mutate(type = "Soil")
)

p_combined <- ggplot(combined_data, aes(x = obs_flux_nmol, y = pred_flux_nmol, 
                                        color = type)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_color_manual(values = c("Trees" = "forestgreen", "Soil" = "brown")) +
  labs(title = "Random Forest Predictions vs Observations",
       subtitle = "Tree and soil CH₄ fluxes",
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Source") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", alpha = 0.8, color = NA)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_COMBINED_nmol.png", p_combined, width = 8, height = 6)
cat("  Saved: QC_PRED_VS_OBS_COMBINED_nmol.png\n")

# Print summary statistics in nmol units
cat("\n=== FLUX RANGES IN NMOL UNITS ===\n")
cat("\nTree fluxes (nmol m⁻² s⁻¹):\n")
cat("  Observed range:", round(range(tree_train$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Predicted range:", round(range(tree_train$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean observed:", round(mean(tree_train$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean predicted:", round(mean(tree_train$pred_flux_nmol, na.rm = TRUE), 3), "\n")

cat("\nSoil fluxes (nmol m⁻² s⁻¹):\n")
cat("  Observed range:", round(range(soil_train$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Predicted range:", round(range(soil_train$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean observed:", round(mean(soil_train$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean predicted:", round(mean(soil_train$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  % negative (uptake):", round(100 * mean(soil_train$obs_flux_nmol < 0, na.rm = TRUE), 1), "%\n")

cat("\n✓ Plot updates complete - all units now in nmol m⁻² s⁻¹\n")



# =============================================================================
# SECTION 11: QUALITY CONTROL PLOTS FOR SPECIES-FIRST WORKFLOW
# =============================================================================

cat("\nGenerating quality control plots...\n")

# Need to ensure tree_train has predictions for plotting
# Since we used tree_train_complete for training, we need to align predictions
tree_train_complete$pred_asinh <- TreeRF$predictions
tree_train_complete$pred_flux <- sinh(tree_train_complete$pred_asinh)

# QC1: Chamber type comparison
qc1_data <- tree_train_complete %>%
  mutate(
    chamber = ifelse(chamber_rigid == 1, "rigid", "semirigid"),
    flux_original = sinh(y_asinh)
  )

p1 <- ggplot(qc1_data, aes(x = chamber, y = flux_original)) +
  geom_boxplot(aes(fill = chamber)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Flux Distribution by Chamber Type",
       subtitle = "Chamber type included as RF predictor",
       x = "Chamber Type",
       y = "Flux (μmol m-2 s-1, pseudo-log scale)") +
  theme_minimal()

ggsave("../../outputs/figures/QC_CHAMBER_COMPARISON.png", p1, width = 8, height = 6)

# QC1b: Outlier diagnostics
outlier_summary <- tibble(
  dataset = c("Tree Fluxes", "Soil Fluxes"),
  n_outliers = c(n_outliers_tree, n_outliers_soil),
  pct_outliers = c(100*n_outliers_tree/(n_outliers_tree + nrow(tree_combined)),
                   100*n_outliers_soil/(n_outliers_soil + nrow(SOIL_YEAR))),
  method = c("1% tails", "MAD (k=8)")
)

p1b <- ggplot(outlier_summary, aes(x = dataset, y = pct_outliers)) +
  geom_col(fill = "coral") +
  geom_text(aes(label = paste0(round(pct_outliers, 1), "%\n(n=", n_outliers, ")")), 
            vjust = -0.5) +
  labs(title = "Outlier Removal Summary",
       subtitle = "Outliers detected on asinh-transformed scale",
       x = "", y = "Percent of Data Removed") +
  ylim(0, max(outlier_summary$pct_outliers) * 1.2) +
  theme_minimal()

ggsave("../../outputs/figures/QC_OUTLIERS.png", p1b, width = 6, height = 6)

# QC2: Predictions vs observations - Trees (μmol)
tree_train_complete$chamber <- ifelse(tree_train_complete$chamber_rigid == 1, "rigid", "semirigid")

p2 <- ggplot(tree_train_complete, aes(x = stem_flux_corrected, y = pred_flux, 
                                      color = chamber)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Tree Predictions vs Observations (Species-First Model)",
       subtitle = paste("R² =", round(TreeRF$r.squared, 3), "- With deconfounded DBH"),
       x = "Observed Flux (μmol m-2 s-1)",
       y = "Predicted Flux (μmol m-2 s-1)",
       color = "Chamber Type") +
  scale_color_manual(values = c("rigid" = "blue", "semirigid" = "red")) +
  theme_minimal()

ggsave("../../outputs/figures/QC_PRED_VS_OBS_TREE.png", p2, width = 8, height = 6)

# QC2b: Predictions vs observations - Soil
soil_train$pred_asinh <- SoilRF$predictions
soil_train$pred_flux <- sinh(soil_train$pred_asinh)

p2b <- ggplot(soil_train[complete_rows_soil,], aes(x = soil_flux_umol_m2_s, y = pred_flux)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Soil Predictions vs Observations",
       subtitle = paste("R² =", round(SoilRF$r.squared, 3), "- With month as cyclic feature"),
       x = "Observed Flux (μmol m-2 s-1)",
       y = "Predicted Flux (μmol m-2 s-1)") +
  theme_minimal()

ggsave("../../outputs/figures/QC_PRED_VS_OBS_SOIL.png", p2b, width = 8, height = 6)

# QC3: DBH deconfounding visualization
p3 <- ggplot(tree_train_complete, aes(x = dbh_m, y = dbh_within_z, color = species_clean)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "DBH Deconfounding: Raw vs Within-Species Standardized",
       subtitle = "Within-species z-scores remove between-species size differences",
       x = "Raw DBH (m)",
       y = "Within-Species DBH (z-score)",
       color = "Species") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("../../outputs/figures/QC_DBH_DECONFOUNDING.png", p3, width = 10, height = 6)

# QC4: Moisture calibration R² by month
p4 <- ggplot(MOISTURE_AFFINE_TABLE, aes(x = month, y = R2_t)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = ifelse(!is.na(R2_t), round(R2_t, 2), "NA")), vjust = -0.5) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Moisture Calibration Quality by Month",
       x = "Month", y = "R²") +
  ylim(0, 1) +
  theme_minimal()

ggsave("../../outputs/figures/QC_MOISTURE_CALIB.png", p4, width = 8, height = 6)

# QC5: Feature importance
tree_imp_raw <- importance(TreeRF)
tree_importance <- data.frame(
  feature = names(tree_imp_raw),
  Importance = as.numeric(tree_imp_raw)
) %>%
  arrange(desc(Importance)) %>%
  head(20)

# Clean up feature names for display
tree_importance$feature_clean <- gsub("species_factor", "Species: ", tree_importance$feature)
tree_importance$feature_clean <- gsub("\\.soil_moisture_at_tree", " × Moisture", tree_importance$feature_clean)

p5a <- ggplot(tree_importance, aes(x = reorder(feature_clean, Importance), y = Importance)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  labs(title = "Tree Model - Top 20 Feature Importances (Species-First)",
       subtitle = "With deconfounded DBH and species×moisture interactions",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("../../outputs/figures/QC_TREE_IMPORTANCE.png", p5a, width = 10, height = 8)

# Soil importance
soil_imp_raw <- importance(SoilRF)
soil_importance <- data.frame(
  feature = names(soil_imp_raw),
  Importance = as.numeric(soil_imp_raw)
) %>%
  arrange(desc(Importance))

p5b <- ggplot(soil_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Soil Model - Feature Importances",
       subtitle = "With month as cyclic feature",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("../../outputs/figures/QC_SOIL_IMPORTANCE.png", p5b, width = 10, height = 6)

# =============================================================================
# NMOL UNIT PLOTS
# =============================================================================

cat("\nGenerating plots with nmol units...\n")

# Convert to nmol for tree plots
tree_train_complete$pred_flux_nmol <- tree_train_complete$pred_flux * 1000
tree_train_complete$obs_flux_nmol <- tree_train_complete$stem_flux_corrected * 1000

p2_nmol <- ggplot(tree_train_complete, aes(x = obs_flux_nmol, y = pred_flux_nmol, 
                                           color = chamber)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Tree Predictions vs Observations (Species-First Model)",
       subtitle = paste("R² =", round(TreeRF$r.squared, 3), "- With deconfounded DBH"),
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Chamber Type") +
  scale_color_manual(values = c("rigid" = "blue", "semirigid" = "red")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_TREE_nmol.png", p2_nmol, width = 8, height = 6)

# Convert to nmol for soil plots
soil_train_subset <- soil_train[complete_rows_soil,]
soil_train_subset$pred_flux_nmol <- soil_train_subset$pred_flux * 1000
soil_train_subset$obs_flux_nmol <- soil_train_subset$soil_flux_umol_m2_s * 1000

p2b_nmol <- ggplot(soil_train_subset, aes(x = obs_flux_nmol, y = pred_flux_nmol)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  labs(title = "Soil Predictions vs Observations",
       subtitle = paste("R² =", round(SoilRF$r.squared, 3)),
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")"))) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_SOIL_nmol.png", p2b_nmol, width = 8, height = 6)

# Combined plot
combined_data <- bind_rows(
  tree_train_complete %>% 
    dplyr::select(obs_flux_nmol, pred_flux_nmol) %>%
    mutate(type = "Trees"),
  soil_train_subset %>%
    dplyr::select(obs_flux_nmol, pred_flux_nmol) %>%
    mutate(type = "Soil")
)

p_combined <- ggplot(combined_data, aes(x = obs_flux_nmol, y = pred_flux_nmol, 
                                        color = type)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_color_manual(values = c("Trees" = "forestgreen", "Soil" = "brown")) +
  labs(title = "Random Forest Predictions vs Observations (Species-First)",
       subtitle = "Tree and soil CH₄ fluxes with deconfounded DBH",
       x = expression(paste("Observed Flux (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Source") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", alpha = 0.8, color = NA)
  )

ggsave("../../outputs/figures/QC_PRED_VS_OBS_COMBINED_nmol.png", p_combined, width = 8, height = 6)

# Print summary statistics
cat("\n=== FLUX RANGES IN NMOL UNITS ===\n")
cat("\nTree fluxes (nmol m⁻² s⁻¹):\n")
cat("  Observed range:", round(range(tree_train_complete$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Predicted range:", round(range(tree_train_complete$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean observed:", round(mean(tree_train_complete$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean predicted:", round(mean(tree_train_complete$pred_flux_nmol, na.rm = TRUE), 3), "\n")

cat("\nSoil fluxes (nmol m⁻² s⁻¹):\n")
cat("  Observed range:", round(range(soil_train_subset$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Predicted range:", round(range(soil_train_subset$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean observed:", round(mean(soil_train_subset$obs_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  Mean predicted:", round(mean(soil_train_subset$pred_flux_nmol, na.rm = TRUE), 3), "\n")
cat("  % negative (uptake):", round(100 * mean(soil_train_subset$obs_flux_nmol < 0, na.rm = TRUE), 1), "%\n")

cat("\n✓ Quality control plots saved\n")

## (Training data already saved earlier, after model fitting)
