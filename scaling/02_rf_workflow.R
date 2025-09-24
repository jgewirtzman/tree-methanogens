# =============================================================================
# RANDOM FOREST CH4 FLUX WORKFLOW - GPS COORDINATES VERSION
# Includes unit conversion, NA handling, DBH correction, and chamber-as-predictor
# UPDATED: Uses GPS coordinates (lat/lon) consistently throughout
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ranger)  # For Random Forest
library(mgcv)    # For smoothing splines
library(broom)   # For model tidying
library(ggplot2)
library(purrr)
library(tibble)

cat("=== RANDOM FOREST CH4 FLUX WORKFLOW - GPS COORDINATES VERSION ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD INPUT DATA
# =============================================================================

cat("Loading prepared data...\n")
load("rf_workflow_input_data_clean.RData")

# Extract datasets from the loaded list
TREE_JULY <- rf_workflow_data$PLACEHOLDER_TREE_JULY
TREE_YEAR <- rf_workflow_data$PLACEHOLDER_TREE_YEAR
SOIL_YEAR <- rf_workflow_data$PLACEHOLDER_SOIL_YEAR
INVENTORY <- rf_workflow_data$PLACEHOLDER_INVENTORY
MOISTURE_DEC_RASTER <- rf_workflow_data$PLACEHOLDER_MOISTURE_DEC_RASTER
moisture_lookup_xy <- rf_workflow_data$PLACEHOLDER_MOISTURE_LOOKUP_XY  # Note: x = lon, y = lat
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS
TAXONOMY <- rf_workflow_data$PLACEHOLDER_TAXONOMY
PLOT_AREA <- rf_workflow_data$PLACEHOLDER_PLOT_AREA

cat("Note: Using GPS coordinates throughout (x = longitude, y = latitude)\n")

# =============================================================================
# DBH CORRECTION FOR INVENTORY
# =============================================================================

cat("\n*** APPLYING DBH CORRECTIONS ***\n")

# Correct DBH units and outliers
INVENTORY <- INVENTORY %>%
  mutate(
    dbh_original = dbh_m,
    # First correction: fix mixed units
    dbh_corrected = case_when(
      dbh_m > 3 ~ dbh_m / 100,  # Values > 3 are likely cm that weren't converted
      dbh_m >= 0.5 & dbh_m <= 3 ~ dbh_m,  # Values 0.5-3 likely already in meters
      TRUE ~ dbh_m  # Small values < 0.5 already correct
    ),
    # Second correction: fix remaining outliers
    dbh_final = case_when(
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 100 ~ dbh_corrected / 100,
      species == "Kalmia latifolia" & (dbh_corrected * 100) > 10 ~ dbh_corrected / 10,
      grepl("Betula", species) & (dbh_corrected * 100) > 200 ~ dbh_corrected / 10,
      species == "Pinus strobus" & (dbh_corrected * 100) > 230 ~ dbh_corrected / 10,
      TRUE ~ dbh_corrected
    ),
    dbh_m = dbh_final  # Replace original with corrected
  ) %>%
  select(tree_id, species, dbh_m, x, y)  # Keep only needed columns (x = lon, y = lat)

n_dbh_corrected <- sum(INVENTORY$dbh_m != rf_workflow_data$PLACEHOLDER_INVENTORY$dbh_m)
cat("  DBH corrections applied to", n_dbh_corrected, "trees\n")
cat("  DBH range (cm):", round(range(INVENTORY$dbh_m * 100), 1), "\n")

# =============================================================================
# CRITICAL FIX: CONVERT UNITS FROM nmol TO μmol
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
cat("  Soil % negative (CH4 uptake):", 
    round(100 * mean(SOIL_YEAR$soil_flux_umol_m2_s < 0, na.rm=T), 1), "%\n")

# Add species to TREE_YEAR
TREE_YEAR <- TREE_YEAR %>%
  left_join(INVENTORY %>% select(tree_id, species), by = "tree_id")

cat("\n✓ Data loaded successfully\n")
cat("  TREE_JULY:", nrow(TREE_JULY), "observations\n")
cat("  TREE_YEAR:", nrow(TREE_YEAR), "observations\n")
cat("  SOIL_YEAR:", nrow(SOIL_YEAR), "observations\n")
cat("  INVENTORY:", nrow(INVENTORY), "trees\n")
cat("  Plot area:", round(PLOT_AREA/10000, 2), "ha\n\n")

# =============================================================================
# SMART IMPUTATION FUNCTION FOR DATA RETENTION
# =============================================================================

smart_impute <- function(df, feature_cols) {
  df_imputed <- df
  
  for(col in feature_cols) {
    if(col %in% names(df) && sum(is.na(df[[col]])) > 0) {
      # For temporal variables, use monthly median imputation
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
      
      # For remaining NAs, use overall median
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
# SECTION 2: MONTHLY MOISTURE MAP VIA AFFINE CALIBRATION
# =============================================================================

cat("Building monthly moisture affine calibration...\n")

# Function to fit affine calibration for a month
# Note: x = longitude, y = latitude
fit_moisture_affine <- function(month_val, points_df, moisture_lookup_fn) {
  # Filter data for this month
  month_data <- points_df %>%
    filter(month == month_val, 
           !is.na(soil_moisture_abs),
           !is.na(x), !is.na(y))  # x = lon, y = lat
  
  if(nrow(month_data) < 5) {
    return(list(alpha = NA, beta = NA, R2 = NA, n = nrow(month_data)))
  }
  
  # Get December moisture at these GPS points
  month_data$moisture_dec <- moisture_lookup_fn(month_data$x, month_data$y)
  
  # Remove NAs
  month_data <- month_data %>% filter(!is.na(moisture_dec))
  
  if(nrow(month_data) < 5) {
    return(list(alpha = NA, beta = NA, R2 = NA, n = nrow(month_data)))
  }
  
  # Fit linear model: observed ~ December
  fit <- lm(soil_moisture_abs ~ moisture_dec, data = month_data)
  
  return(list(
    alpha = coef(fit)[1],
    beta = coef(fit)[2],
    R2 = summary(fit)$r.squared,
    n = nrow(month_data)
  ))
}

# Combine all moisture observations
all_moisture_points <- bind_rows(
  TREE_JULY %>% select(month, x, y, soil_moisture_abs),
  TREE_YEAR %>% select(month, x, y, soil_moisture_abs),
  SOIL_YEAR %>% select(month, x, y, soil_moisture_abs)
) %>%
  filter(!is.na(soil_moisture_abs))

# Fit affine model for each month
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

print(MOISTURE_AFFINE_TABLE)

# Check moisture calibration availability
cat("\nChecking moisture calibration status:\n")
for(m in 1:12) {
  params <- MOISTURE_AFFINE_TABLE %>% filter(month == m)
  if(is.na(params$alpha_t) || is.na(params$beta_t)) {
    cat("  Month", m, ": NO calibration available (", params$n_points, "points)\n")
  } else {
    cat("  Month", m, ": Calibrated with", params$n_points, "points (R² =", 
        round(params$R2_t, 3), ")\n")
  }
}

# Function to predict moisture at any GPS location for any month (with fallback)
make_Mhat <- function(affine_table, moisture_lookup_fn) {
  function(month_vec, x_vec, y_vec) {  # x = lon, y = lat
    # Handle vectorized input
    if(length(month_vec) == 0) return(numeric(0))
    
    moisture_pred <- numeric(length(month_vec))
    
    # Get mean moisture for months with no calibration
    months_with_calib <- affine_table %>% 
      filter(!is.na(alpha_t)) %>% 
      pull(month)
    
    if(length(months_with_calib) > 0) {
      # Calculate mean moisture from calibrated months as fallback
      fallback_moisture <- mean(affine_table$alpha_t[!is.na(affine_table$alpha_t)])
    } else {
      fallback_moisture <- 0.15  # Default if no calibration exists
    }
    
    # Process each unique month
    unique_months <- unique(month_vec)
    for(m in unique_months) {
      idx <- which(month_vec == m)
      params <- affine_table %>% filter(month == m)
      
      if(nrow(params) == 0 || is.na(params$alpha_t)) {
        # Use fallback for uncalibrated months
        moisture_pred[idx] <- fallback_moisture
        # Only print message once per month
        if(length(idx) > 0 && idx[1] == 1) {
          cat("    Using fallback moisture", round(fallback_moisture, 3), 
              "for month", m, "\n")
        }
      } else {
        moisture_dec <- moisture_lookup_fn(x_vec[idx], y_vec[idx])
        pred <- params$alpha_t + params$beta_t * moisture_dec
        
        # Clip to reasonable bounds
        pred[pred < 0] <- 0
        pred[pred > 0.6] <- 0.75
        pred[is.na(pred)] <- fallback_moisture  # Use fallback for NAs
        
        moisture_pred[idx] <- pred
      }
    }
    
    return(moisture_pred)
  }
}

Mhat <- make_Mhat(MOISTURE_AFFINE_TABLE, moisture_lookup_xy)
cat("✓ Monthly moisture calibration complete\n\n")

# =============================================================================
# SECTION 2.5: CHAMBER AS PREDICTOR (NO CALIBRATION)
# =============================================================================

cat("Preparing data with chamber type as predictor...\n")

# Less aggressive outlier removal function (k=8 instead of k=6)
remove_outliers_mad <- function(x, k = 8) {
  # Remove outliers on transformed scale
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  
  if(mad_val == 0 || is.na(mad_val)) {
    # If MAD is 0, use IQR method instead
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    if(iqr == 0) {
      # If IQR is also 0, keep all non-NA values
      return(!is.na(x))
    }
    
    lower <- q1 - 6 * iqr  # Less aggressive than original
    upper <- q3 + 6 * iqr
    return(x >= lower & x <= upper)
  } else {
    # MAD-based outlier detection
    lower <- med - k * mad_val
    upper <- med + k * mad_val
    return(x >= lower & x <= upper)
  }
}

# Stack tree data with chamber type
tree_combined <- bind_rows(
  TREE_JULY %>% mutate(chamber_type = "rigid"),
  TREE_YEAR %>% mutate(chamber_type = coalesce(chamber_type, "semirigid"))
) %>%
  filter(!is.na(stem_flux_umol_m2_s))

# Apply smart imputation to retain more data
cat("  Applying smart imputation to tree data...\n")
n_before_impute <- sum(complete.cases(tree_combined[, c("air_temp_C", "soil_temp_C", "soil_moisture_abs")]))
tree_combined <- smart_impute(tree_combined, 
                              c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m"))
n_after_impute <- sum(complete.cases(tree_combined[, c("air_temp_C", "soil_temp_C", "soil_moisture_abs")]))
cat("    Complete cases before/after imputation:", n_before_impute, "/", n_after_impute, "\n")

# Transform to asinh scale
tree_combined$z <- asinh(tree_combined$stem_flux_umol_m2_s)

# Remove outliers on transformed scale
outlier_mask_tree <- remove_outliers_mad(tree_combined$z, k = 20)
n_outliers_tree <- sum(!outlier_mask_tree)
cat("  Removing", n_outliers_tree, "tree flux outliers (", 
    round(100*n_outliers_tree/nrow(tree_combined), 1), "%)\n")

tree_combined <- tree_combined[outlier_mask_tree, ]

# No chamber correction - use flux directly
# Convert from nmol to μmol and create corrected column
tree_combined$stem_flux_corrected <- tree_combined$stem_flux_umol_m2_s 
tree_combined$z_corr <- asinh(tree_combined$stem_flux_corrected)  # Recalculate on corrected scale

# For compatibility with downstream code, create dummy chamber calibration variables
beta_1 <- 0  # No offset since chamber is a predictor
beta_1_se <- 0
beta_1_ci <- c(0, 0)
calib_summary <- list(r.squared = NA)  # Placeholder

cat("✓ Chamber type will be used as predictor in RF model\n\n")

# =============================================================================
# SECTION 2.6: SOIL DATA OUTLIER REMOVAL
# =============================================================================

cat("Removing outliers from soil data...\n")

# Apply smart imputation to soil data
n_before_soil <- sum(complete.cases(SOIL_YEAR[, c("air_temp_C", "soil_temp_C", "soil_moisture_abs")]))
SOIL_YEAR <- smart_impute(SOIL_YEAR, 
                          c("air_temp_C", "soil_temp_C", "soil_moisture_abs"))
n_after_soil <- sum(complete.cases(SOIL_YEAR[, c("air_temp_C", "soil_temp_C", "soil_moisture_abs")]))
cat("  Complete cases before/after imputation:", n_before_soil, "/", n_after_soil, "\n")

# Apply outlier removal to soil data
SOIL_YEAR$z_soil <- asinh(SOIL_YEAR$soil_flux_umol_m2_s)
outlier_mask_soil <- remove_outliers_mad(SOIL_YEAR$z_soil, k = 8)
n_outliers_soil <- sum(!outlier_mask_soil)
cat("  Removing", n_outliers_soil, "soil flux outliers (", 
    round(100*n_outliers_soil/nrow(SOIL_YEAR), 1), "%)\n")

SOIL_YEAR <- SOIL_YEAR[outlier_mask_soil, ]
SOIL_YEAR$z_soil <- NULL  # Remove temporary column

cat("✓ Outlier removal complete\n")
cat("  Final tree observations:", nrow(tree_combined), "\n")
cat("  Final soil observations:", nrow(SOIL_YEAR), "\n\n")

# =============================================================================
# SECTION 3: EMPIRICAL SMOOTHED SEASONAL INDEX
# =============================================================================

cat("Computing empirical seasonal indices...\n")

# Function to compute empirical SI
compute_empirical_SI <- function(df, flux_col, group_name, use_soil_temp = FALSE) {
  # Prepare data
  df_clean <- df %>%
    filter(!is.na(!!sym(flux_col))) %>%
    mutate(y_asinh = asinh(!!sym(flux_col)))
  
  # Build basic features (no month/seasonality)
  if(group_name == "tree") {
    # Tree features without seasonality
    X <- df_clean %>%
      select(any_of(c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m")))
  } else {
    # Soil features without seasonality  
    X <- df_clean %>%
      select(any_of(c("soil_temp_C", "air_temp_C", "soil_moisture_abs")))
  }
  
  # Remove columns with too many NAs
  X <- X[, colSums(is.na(X)) < nrow(X)*0.7, drop = FALSE]  # Allow up to 70% NA
  
  # Train baseline RF - with NA handling if using ranger
  if(ncol(X) > 0 && nrow(X) > 20) {
    rf0 <- ranger(
      x = X,
      y = df_clean$y_asinh,
      num.trees = 200,
      min.node.size = 5,
      mtry = max(1, floor(sqrt(ncol(X)))),
      num.threads = 1,  # Enable NA handling
      keep.inbag = TRUE
    )
    
    # Get OOB predictions
    yhat_asinh <- rf0$predictions
  } else {
    # Fallback: use mean
    yhat_asinh <- rep(mean(df_clean$y_asinh), nrow(df_clean))
  }
  
  # Compute residuals
  df_clean$residual <- df_clean$y_asinh - yhat_asinh
  
  # Aggregate by month
  monthly_resid <- df_clean %>%
    group_by(month) %>%
    summarise(
      mean_resid = mean(residual, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Ensure all months present
  all_months <- tibble(month = 1:12)
  monthly_resid <- all_months %>%
    left_join(monthly_resid, by = "month") %>%
    mutate(mean_resid = ifelse(is.na(mean_resid), 0, mean_resid))
  
  # Smooth with cyclic spline
  if(sum(!is.na(monthly_resid$mean_resid)) >= 4) {
    # Duplicate data for cyclic continuity
    extended <- bind_rows(
      monthly_resid %>% mutate(month = month - 12),
      monthly_resid,
      monthly_resid %>% mutate(month = month + 12)
    )
    
    # Fit smooth spline
    smooth_fit <- smooth.spline(extended$month, extended$mean_resid, df = 4)
    
    # Predict for months 1-12
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

# Compute SI for trees (using corrected fluxes)
tree_for_SI <- tree_combined %>%
  select(month, stem_flux_corrected, air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m) %>%
  rename(flux = stem_flux_corrected)

SI_tree <- compute_empirical_SI(tree_for_SI, "flux", "tree")

# Compute SI for soil
soil_for_SI <- SOIL_YEAR %>%
  rename(flux = soil_flux_umol_m2_s)

SI_soil <- compute_empirical_SI(soil_for_SI, "flux", "soil", use_soil_temp = TRUE)

# Combine SI tables
SI_TABLES <- bind_rows(SI_tree, SI_soil)

print(SI_TABLES)
cat("✓ Seasonal indices computed\n\n")

# =============================================================================
# SECTION 4: FEATURE ENGINEERING
# =============================================================================

cat("Building features for Random Forest models...\n")

# 4.3: Compute taxonomy priors first (needed for tree features)
compute_taxon_prior <- function(tree_df, taxonomy_df) {
  # Get environmental-only predictions
  df_clean <- tree_df %>%
    filter(!is.na(stem_flux_corrected)) %>%
    mutate(y_asinh = asinh(stem_flux_corrected))
  
  # Environmental features only
  X_env <- df_clean %>%
    select(any_of(c("air_temp_C", "soil_temp_C", "soil_moisture_abs", "dbh_m")))
  
  # Train environmental-only RF with NA handling
  rf_env <- ranger(
    x = X_env,
    y = df_clean$y_asinh,
    num.trees = 200,
    min.node.size = 5,
    num.threads = 1  # Enable NA handling
  )
  
  # Compute residuals
  df_clean$residual <- df_clean$y_asinh - rf_env$predictions
  
  # Join with taxonomy
  df_clean <- df_clean %>%
    left_join(taxonomy_df %>% select(species, genus, family, order, class),
              by = "species")
  
  # Compute medians at each level (require n>=5)
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
  
  med_order <- df_clean %>%
    group_by(order) %>%
    summarise(
      med_resid = ifelse(n() >= 5, median(residual, na.rm = TRUE), NA_real_),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(med_resid))
  
  med_class <- df_clean %>%
    group_by(class) %>%
    summarise(
      med_resid = ifelse(n() >= 5, median(residual, na.rm = TRUE), NA_real_),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(med_resid))
  
  return(list(
    species = med_species,
    genus = med_genus,
    family = med_family,
    order = med_order,
    class = med_class
  ))
}

# Compute taxonomy priors
TAXONOMY_PRIORS <- compute_taxon_prior(tree_combined, TAXONOMY)
cat("✓ Taxonomy priors computed\n")

# 4.1: Build tree features (ensuring numeric) - UPDATED WITH MONTH
build_features_tree <- function(df, drivers, Mhat_fn, SI_table, taxonomy, taxon_priors) {
  # Join with drivers
  df <- df %>%
    left_join(drivers %>% select(month, air_temp_C_mean, soil_temp_C_mean),
              by = "month")
  
  # Add SI
  si_tree <- SI_table %>% filter(group == "tree") %>% select(month, SI)
  df <- df %>%
    left_join(si_tree, by = "month")
  
  # Ensure numeric conversion for key fields
  df$air_temp_C_mean <- as.numeric(df$air_temp_C_mean)
  df$soil_temp_C_mean <- as.numeric(df$soil_temp_C_mean)
  df$SI <- as.numeric(df$SI)
  df$dbh_m <- as.numeric(df$dbh_m)
  
  # ADD MONTH AS CYCLIC FEATURES
  df$month_sin <- sin(2 * pi * df$month / 12)
  df$month_cos <- cos(2 * pi * df$month / 12)
  
  # Add predicted moisture if needed
  if("soil_moisture_abs" %in% names(df)) {
    # Use observed moisture when available
    df$soil_moisture_at_tree <- as.numeric(df$soil_moisture_abs)
  } else {
    df$soil_moisture_at_tree <- NA_real_
  }
  
  # Fill missing moisture with predicted (x = lon, y = lat)
  missing_moisture <- is.na(df$soil_moisture_at_tree)
  if(any(missing_moisture) && !is.null(Mhat_fn)) {
    df$soil_moisture_at_tree[missing_moisture] <- 
      Mhat_fn(df$month[missing_moisture], df$x[missing_moisture], df$y[missing_moisture])
  }
  df$soil_moisture_at_tree <- as.numeric(df$soil_moisture_at_tree)
  
  # Add chamber type as feature
  if("chamber_type" %in% names(df)) {
    df$chamber_rigid <- as.numeric(df$chamber_type == "rigid")
    df$chamber_semirigid <- as.numeric(df$chamber_type == "semirigid")
  } else {
    # For prediction, default to rigid
    df$chamber_rigid <- 1
    df$chamber_semirigid <- 0
  }
  
  # Join taxonomy info
  df <- df %>%
    left_join(taxonomy %>% select(species, genus, family), by = "species")
  
  # Add taxon prior
  df$taxon_prior_asinh <- 0  # Default
  
  # Look up priors in order: species -> genus -> family -> order -> class
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
  
  # Create one-hot encodings
  # Species (collapse rare)
  species_counts <- table(df$species)
  df$species_clean <- ifelse(species_counts[df$species] >= 10, 
                             df$species, "SPECIES_OTHER")
  
  # Genus (collapse rare)
  genus_counts <- table(df$genus)
  df$genus_clean <- ifelse(genus_counts[df$genus] >= 10,
                           df$genus, "GENUS_OTHER")
  
  # Family (collapse rare)  
  family_counts <- table(df$family)
  df$family_clean <- ifelse(family_counts[df$family] >= 10,
                            df$family, "FAMILY_OTHER")
  
  # Create interaction terms
  df$moisture_x_airT <- as.numeric(df$soil_moisture_at_tree) * as.numeric(df$air_temp_C_mean)
  df$moisture_x_soilT <- as.numeric(df$soil_moisture_at_tree) * as.numeric(df$soil_temp_C_mean)
  
  # Target (only if flux data is present)
  if("stem_flux_corrected" %in% names(df) && !all(is.na(df$stem_flux_corrected))) {
    df$y_asinh <- asinh(as.numeric(df$stem_flux_corrected))
  } else if("stem_flux_umol_m2_s" %in% names(df) && !all(is.na(df$stem_flux_umol_m2_s))) {
    df$y_asinh <- asinh(as.numeric(df$stem_flux_umol_m2_s))
  } else {
    # For prediction on inventory without flux data
    df$y_asinh <- NA_real_
  }
  
  return(df)
}

# 4.2: Build soil features
build_features_soil <- function(df, drivers, Mhat_fn, SI_table) {
  # Join with drivers
  df <- df %>%
    left_join(drivers %>% select(month, air_temp_C_mean, soil_temp_C_mean),
              by = "month")
  
  # Add SI
  si_soil <- SI_table %>% filter(group == "soil") %>% select(month, SI)
  df <- df %>%
    left_join(si_soil, by = "month")
  
  # Ensure numeric
  df$air_temp_C_mean <- as.numeric(df$air_temp_C_mean)
  df$soil_temp_C_mean <- as.numeric(df$soil_temp_C_mean)
  df$SI <- as.numeric(df$SI)
  
  # ADD MONTH AS CYCLIC FEATURES
  df$month_sin <- sin(2 * pi * df$month / 12)
  df$month_cos <- cos(2 * pi * df$month / 12)
  
  # Add predicted moisture if needed
  if("soil_moisture_abs" %in% names(df)) {
    df$soil_moisture_at_site <- as.numeric(df$soil_moisture_abs)
  } else {
    df$soil_moisture_at_site <- NA_real_
  }
  
  # Fill missing moisture (x = lon, y = lat)
  missing_moisture <- is.na(df$soil_moisture_at_site)
  if(any(missing_moisture) && !is.null(Mhat_fn)) {
    df$soil_moisture_at_site[missing_moisture] <- 
      Mhat_fn(df$month[missing_moisture], df$x[missing_moisture], df$y[missing_moisture])
  }
  df$soil_moisture_at_site <- as.numeric(df$soil_moisture_at_site)
  
  # Interaction terms
  df$moisture_x_soilT <- as.numeric(df$soil_moisture_at_site) * as.numeric(df$soil_temp_C_mean)
  df$moisture_x_airT <- as.numeric(df$soil_moisture_at_site) * as.numeric(df$air_temp_C_mean)
  
  # Target
  if("soil_flux_umol_m2_s" %in% names(df)) {
    df$y_asinh <- asinh(as.numeric(df$soil_flux_umol_m2_s))
  }
  
  return(df)
}

# Build training features
tree_train <- build_features_tree(tree_combined, DRIVERS, Mhat, SI_TABLES, 
                                  TAXONOMY, TAXONOMY_PRIORS)
soil_train <- build_features_soil(SOIL_YEAR, DRIVERS, Mhat, SI_TABLES)

cat("✓ Features built for training\n\n")

# =============================================================================
# SECTION 5: RANDOM FOREST MODELS WITH NA HANDLING - UPDATED WITH MONTH
# =============================================================================

cat("Training Random Forest models with NA handling and month features...\n")

# Prepare tree features (including chamber type AND month)
tree_features <- c("dbh_m", "air_temp_C_mean", "soil_temp_C_mean",
                   "soil_moisture_at_tree", "SI", "taxon_prior_asinh",
                   "moisture_x_airT", "moisture_x_soilT",
                   "chamber_rigid", "chamber_semirigid",
                   "month_sin", "month_cos")  # ADDED MONTH FEATURES

# Create feature matrix more carefully
# First extract the base features
X_tree_base <- tree_train[, tree_features, drop = FALSE]

# Handle NAs in categorical variables before creating dummies
tree_train$species_clean[is.na(tree_train$species_clean)] <- "SPECIES_OTHER"
tree_train$genus_clean[is.na(tree_train$genus_clean)] <- "GENUS_OTHER"
tree_train$family_clean[is.na(tree_train$family_clean)] <- "FAMILY_OTHER"

# Create dummy variables ensuring same row alignment
if(length(unique(tree_train$species_clean)) > 1) {
  species_dummies <- model.matrix(~ species_clean - 1, data = tree_train)
} else {
  species_dummies <- matrix(1, nrow = nrow(tree_train), ncol = 1)
  colnames(species_dummies) <- paste0("species_clean", tree_train$species_clean[1])
}

if(length(unique(tree_train$genus_clean)) > 1) {
  genus_dummies <- model.matrix(~ genus_clean - 1, data = tree_train)
} else {
  genus_dummies <- matrix(1, nrow = nrow(tree_train), ncol = 1)
  colnames(genus_dummies) <- paste0("genus_clean", tree_train$genus_clean[1])
}

if(length(unique(tree_train$family_clean)) > 1) {
  family_dummies <- model.matrix(~ family_clean - 1, data = tree_train)
} else {
  family_dummies <- matrix(1, nrow = nrow(tree_train), ncol = 1)
  colnames(family_dummies) <- paste0("family_clean", tree_train$family_clean[1])
}

# Verify dimensions match
cat("  Base features:", nrow(X_tree_base), "rows\n")
cat("  Species dummies:", nrow(species_dummies), "rows\n")
cat("  Genus dummies:", nrow(genus_dummies), "rows\n")
cat("  Family dummies:", nrow(family_dummies), "rows\n")

# Combine all features
X_tree <- cbind(
  as.data.frame(X_tree_base),
  as.data.frame(species_dummies),
  as.data.frame(genus_dummies),
  as.data.frame(family_dummies)
)

# Add species interaction with moisture (only for common species)
common_species <- names(table(tree_train$species_clean))[table(tree_train$species_clean) >= 10]
for(sp in common_species) {
  if(sp != "SPECIES_OTHER") {
    col_name <- paste0("moisture_x_", gsub("[^A-Za-z0-9]", "_", sp))
    X_tree[[col_name]] <- tree_train$soil_moisture_at_tree * 
      (tree_train$species_clean == sp)
  }
}

# Remove columns with zero variance (but keep columns with NAs for ranger)
zero_var <- sapply(X_tree, function(x) {
  if(all(is.na(x))) return(TRUE)
  var(x, na.rm = TRUE) == 0
})
X_tree <- X_tree[, !zero_var, drop = FALSE]

y_tree <- tree_train$y_asinh

# Remove any remaining rows with NA in target
complete_rows <- !is.na(y_tree)
X_tree <- X_tree[complete_rows, , drop = FALSE]
y_tree <- y_tree[complete_rows]

cat("  Tree feature matrix:", nrow(X_tree), "obs x", ncol(X_tree), "features\n")
cat("  Month features included: month_sin, month_cos\n")

# Store the training feature template for prediction
X_tree_template <- X_tree

# Train TreeRF with NA handling
TreeRF <- ranger(
  x = as.data.frame(X_tree),
  y = y_tree,
  num.trees = 800,
  min.node.size = 5,
  mtry = floor(sqrt(ncol(X_tree))),
  importance = "impurity",
  num.threads = 1,  # Required for NA handling
  oob.error = TRUE,
  seed = 42
)

cat("TreeRF trained:\n")
cat("  OOB R²:", round(TreeRF$r.squared, 3), "\n")
cat("  OOB RMSE (asinh):", round(sqrt(TreeRF$prediction.error), 4), "\n")
# Check if RMSE is suspiciously small
if(sqrt(TreeRF$prediction.error) < 0.0001) {
  cat("  Note: RMSE very small - check scale of tree fluxes\n")
}

# Prepare soil features WITH MONTH
soil_features <- c("soil_temp_C_mean", "air_temp_C_mean",
                   "soil_moisture_at_site", "SI",
                   "moisture_x_soilT", "moisture_x_airT",
                   "month_sin", "month_cos")  # ADDED MONTH FEATURES

X_soil <- soil_train[, soil_features, drop = FALSE]

# Remove columns with zero variance
zero_var_soil <- sapply(X_soil, function(x) {
  if(all(is.na(x))) return(TRUE)
  var(x, na.rm = TRUE) == 0
})
X_soil <- X_soil[, !zero_var_soil, drop = FALSE]

y_soil <- soil_train$y_asinh

# Remove rows with NA in target
complete_rows_soil <- !is.na(y_soil)
X_soil <- X_soil[complete_rows_soil, , drop = FALSE]
y_soil <- y_soil[complete_rows_soil]

cat("  Soil feature matrix:", nrow(X_soil), "obs x", ncol(X_soil), "features\n")
cat("  Month features included: month_sin, month_cos\n")

# Store soil template for prediction
X_soil_template <- X_soil

# Train SoilRF with NA handling
SoilRF <- ranger(
  x = as.data.frame(X_soil),
  y = y_soil,
  num.trees = 800,
  min.node.size = 5,
  mtry = floor(sqrt(ncol(X_soil))),
  importance = "impurity",
  num.threads = 1,  # Required for NA handling
  oob.error = TRUE,
  seed = 42
)

cat("SoilRF trained:\n")
cat("  OOB R²:", round(SoilRF$r.squared, 3), "\n")
cat("  OOB RMSE (asinh):", round(sqrt(SoilRF$prediction.error), 4), "\n\n")

# =============================================================================
# SECTION 7: GEOMETRY CALCULATIONS 
# =============================================================================

cat("Computing plot geometry...\n")

# 7.1: Compute geometry
compute_geometry <- function(inventory_df) {
  # Stem lateral area to 2m height (0.75m for Kalmia)
  inventory_df$S_i <- ifelse(
    inventory_df$species == "Kalmia latifolia",
    pi * inventory_df$dbh_m * 0.75,  # 75cm height for Kalmia
    pi * inventory_df$dbh_m * 2      # 2m height for trees
  )
  
  # Basal area
  inventory_df$BA_i <- pi * (inventory_df$dbh_m/2)^2
  
  # Total basal area
  total_BA <- sum(inventory_df$BA_i)
  
  # Soil ground area
  A_soil <- PLOT_AREA - total_BA
  
  return(list(
    inventory = inventory_df,
    A_soil = A_soil,
    A_plot = PLOT_AREA
  ))
}

geometry <- compute_geometry(INVENTORY)
cat("Plot geometry computed:\n")
cat("  Total basal area:", round(sum(geometry$inventory$BA_i), 1), "m²\n")
cat("  Soil area:", round(geometry$A_soil, 1), "m²\n")
cat("  Soil fraction:", round(geometry$A_soil/geometry$A_plot, 3), "\n\n")

# =============================================================================
# SECTION 10: MONTHLY SPATIAL PREDICTIONS - UPDATED WITH MONTH
# =============================================================================

cat("\n=== SECTION 10: MONTHLY SPATIAL PREDICTIONS ===\n")

# Training ranges for clipping (p1, p99)
training_ranges <- list(
  air_temp = quantile(c(tree_combined$air_temp_C, SOIL_YEAR$air_temp_C), 
                      probs = c(0.01, 0.99), na.rm = TRUE),
  soil_temp = quantile(c(tree_combined$soil_temp_C, SOIL_YEAR$soil_temp_C), 
                       probs = c(0.01, 0.99), na.rm = TRUE),
  moisture = quantile(c(tree_combined$soil_moisture_abs, SOIL_YEAR$soil_moisture_abs), 
                      probs = c(0.01, 0.99), na.rm = TRUE),
  dbh = quantile(INVENTORY$dbh_m, probs = c(0.01, 0.99), na.rm = TRUE)
)

cat("Training ranges for clipping:\n")
print(training_ranges)

# Clip function for extrapolation guard
clip_to_range <- function(x, range_vals) {
  if (any(is.na(range_vals))) return(x)
  pmax(pmin(x, range_vals[2]), range_vals[1])
}

# Initialize results storage
monthly_results <- tibble(month = 1:12)

# Extract species priors as named vectors from the dataframes
species_priors <- setNames(
  TAXONOMY_PRIORS$species$med_resid,
  TAXONOMY_PRIORS$species$species
)
genus_priors <- setNames(
  TAXONOMY_PRIORS$genus$med_resid,
  TAXONOMY_PRIORS$genus$genus
)
family_priors <- setNames(
  TAXONOMY_PRIORS$family$med_resid,
  TAXONOMY_PRIORS$family$family
)

# =============================================================================
# 10.1: TREE PREDICTIONS FROM INVENTORY - UPDATED WITH MONTH
# =============================================================================

cat("\nPredicting tree fluxes from inventory (with month features)...\n")

# Storage for monthly tree fluxes
tree_fluxes_monthly <- list()

for (t in 1:12) {
  cat("  Month", t, "...")
  
  # Get monthly drivers
  air_temp_t <- DRIVERS$air_temp_C_mean[t]
  soil_temp_t <- DRIVERS$soil_temp_C_mean[t]
  
  # Handle missing soil temp
  if (is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  # Get seasonal index for trees from SI_TABLES
  SI_tree_row <- SI_TABLES %>% 
    filter(group == "tree", month == t)
  SI_tree_t <- if(nrow(SI_tree_row) > 0) SI_tree_row$SI[1] else 0
  
  # Predict for each inventory tree (x = lon, y = lat)
  inv_predictions <- INVENTORY %>%
    mutate(
      # Look up moisture at tree GPS location using calibrated monthly values
      moisture_raw = moisture_lookup_xy(x, y),
      soil_moisture_abs = clip_to_range(
        MOISTURE_AFFINE_TABLE$alpha_t[t] + 
          MOISTURE_AFFINE_TABLE$beta_t[t] * moisture_raw,
        training_ranges$moisture
      ),
      
      # Clip environmental drivers to training range
      air_temp_C = clip_to_range(air_temp_t, training_ranges$air_temp),
      soil_temp_C = clip_to_range(soil_temp_t, training_ranges$soil_temp),
      dbh_m_clipped = clip_to_range(dbh_m, training_ranges$dbh),
      
      # Add seasonal index
      SI_tree = SI_tree_t,
      
      # ADD MONTH AS CYCLIC FEATURES
      month_sin = sin(2 * pi * t / 12),
      month_cos = cos(2 * pi * t / 12),
      
      # Set chamber type for prediction (default to rigid)
      chamber_rigid = 1,
      chamber_semirigid = 0,
      
      # Create interactions
      moisture_x_airT = soil_moisture_abs * air_temp_C,
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      
      # Extract genus from species name
      genus = sub(" .*", "", species),
      
      # Assign family
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
  
  # Assign taxon priors
  inv_predictions$taxon_prior_asinh <- 0  # Default
  
  # Try species level first
  for (i in 1:nrow(inv_predictions)) {
    sp <- inv_predictions$species[i]
    if (sp %in% names(species_priors)) {
      inv_predictions$taxon_prior_asinh[i] <- species_priors[sp]
    }
  }
  
  # For trees with no species match, try genus
  no_match <- inv_predictions$taxon_prior_asinh == 0
  for (i in which(no_match)) {
    gn <- inv_predictions$genus[i]
    if (gn %in% names(genus_priors)) {
      inv_predictions$taxon_prior_asinh[i] <- genus_priors[gn]
    }
  }
  
  # For trees still with no match, try family
  no_match <- inv_predictions$taxon_prior_asinh == 0
  for (i in which(no_match)) {
    fm <- inv_predictions$family[i]
    if (fm %in% names(family_priors)) {
      inv_predictions$taxon_prior_asinh[i] <- family_priors[fm]
    }
  }
  
  # Build feature matrix - numeric features first (INCLUDING MONTH)
  numeric_features <- inv_predictions %>%
    select(dbh_m_clipped, air_temp_C, soil_temp_C, 
           soil_moisture_abs, SI_tree, moisture_x_airT, 
           moisture_x_soilT, taxon_prior_asinh,
           chamber_rigid, chamber_semirigid,
           month_sin, month_cos) %>%  # ADDED MONTH FEATURES
    as.matrix()
  
  # Add categorical features as one-hot encoding
  species_dummies <- model.matrix(~ species - 1, data = inv_predictions)
  genus_dummies <- model.matrix(~ genus - 1, data = inv_predictions)
  family_dummies <- model.matrix(~ family - 1, data = inv_predictions)
  
  # Combine all features
  X_pred_tree <- cbind(numeric_features, species_dummies, genus_dummies, family_dummies)
  
  # Create aligned matrix with same columns as training
  X_pred_aligned <- matrix(0, nrow = nrow(X_pred_tree), ncol = ncol(X_tree))
  colnames(X_pred_aligned) <- colnames(X_tree)
  
  # Match columns that exist in both
  for (col_name in colnames(X_tree)) {
    if (col_name %in% colnames(X_pred_tree)) {
      X_pred_aligned[, col_name] <- X_pred_tree[, col_name]
    }
  }
  
  # Predict on asinh scale
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  
  # Back-transform
  pred_flux_umol_m2_s <- sinh(pred_asinh)
  
  # Calculate stem surface area (2m height for trees, 0.75m for Kalmia)
  S_tree <- ifelse(
    inv_predictions$species == "Kalmia latifolia",
    pi * inv_predictions$dbh_m_clipped * 0.75,
    pi * inv_predictions$dbh_m_clipped * 2
  )
  
  # Tree-level emission rate
  R_tree <- pred_flux_umol_m2_s * S_tree
  
  # Store results
  tree_fluxes_monthly[[t]] <- list(
    month = t,
    total_R_tree = sum(R_tree, na.rm = TRUE),
    mean_flux = mean(pred_flux_umol_m2_s, na.rm = TRUE),
    median_flux = median(pred_flux_umol_m2_s, na.rm = TRUE),
    n_trees = sum(!is.na(pred_flux_umol_m2_s))
  )
  
  cat(" done. Total emission:", round(sum(R_tree, na.rm = TRUE), 6), "μmol/s\n")
}

# Calculate plot-basis tree flux for each month
tree_results <- bind_rows(tree_fluxes_monthly) %>%
  mutate(
    Phi_tree_umol_m2_s = total_R_tree / geometry$A_plot
  )

cat("  Tree predictions complete.\n")

# =============================================================================
# 10.2: SOIL PREDICTIONS ON SPATIAL GRID - UPDATED WITH MONTH
# =============================================================================

cat("\nPredicting soil fluxes on spatial grid (with month features)...\n")

# Create grid in GPS coordinates (x = lon, y = lat)
x_range <- range(INVENTORY$x, na.rm = TRUE)  # Longitude range
y_range <- range(INVENTORY$y, na.rm = TRUE)  # Latitude range
grid_res <- 50
x_grid <- seq(x_range[1], x_range[2], length.out = grid_res)
y_grid <- seq(y_range[1], y_range[2], length.out = grid_res)
soil_grid <- expand.grid(x = x_grid, y = y_grid)  # x = lon, y = lat

# Calculate pixel area (convert degrees to meters)
lon_dist_m <- diff(x_grid)[1] * 111320 * cos(mean(y_grid) * pi/180)
lat_dist_m <- diff(y_grid)[1] * 111320
pixel_area <- lon_dist_m * lat_dist_m

cat("  Created", nrow(soil_grid), "grid points\n")
cat("  Pixel area:", round(pixel_area, 2), "m²\n")

soil_fluxes_monthly <- list()

for (t in 1:12) {
  cat("  Month", t, "...")
  
  # Get monthly drivers
  air_temp_t <- DRIVERS$air_temp_C_mean[t]
  soil_temp_t <- DRIVERS$soil_temp_C_mean[t]
  
  if (is.na(soil_temp_t)) {
    soil_temp_t <- air_temp_t * 0.8
  }
  
  # Get seasonal index for soil from SI_TABLES
  SI_soil_row <- SI_TABLES %>% 
    filter(group == "soil", month == t)
  SI_soil_t <- if(nrow(SI_soil_row) > 0) SI_soil_row$SI[1] else 0
  
  # Predict for each grid point (x = lon, y = lat)
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
      # ADD MONTH AS CYCLIC FEATURES
      month_sin = sin(2 * pi * t / 12),
      month_cos = cos(2 * pi * t / 12),
      moisture_x_soilT = soil_moisture_abs * soil_temp_C,
      moisture_x_airT = soil_moisture_abs * air_temp_C
    )
  
  # Build feature matrix (INCLUDING MONTH)
  X_pred_soil <- grid_predictions %>%
    select(soil_temp_C, soil_moisture_abs, SI_soil, 
           moisture_x_soilT, air_temp_C, moisture_x_airT,
           month_sin, month_cos) %>%  # ADDED MONTH FEATURES
    as.matrix()
  
  # Create aligned matrix
  X_pred_soil_aligned <- matrix(0, nrow = nrow(X_pred_soil), ncol = ncol(X_soil))
  colnames(X_pred_soil_aligned) <- colnames(X_soil)
  
  # Match columns
  for (col_name in colnames(X_soil)) {
    if (col_name %in% colnames(X_pred_soil)) {
      X_pred_soil_aligned[, col_name] <- X_pred_soil[, col_name]
    }
  }
  
  # Predict
  pred_asinh_soil <- predict(SoilRF, X_pred_soil_aligned)$predictions
  pred_flux_soil <- sinh(pred_asinh_soil)
  mean_soil_flux <- mean(pred_flux_soil, na.rm = TRUE)
  
  # Store results
  soil_fluxes_monthly[[t]] <- list(
    month = t,
    mean_flux = mean_soil_flux,
    median_flux = median(pred_flux_soil, na.rm = TRUE),
    n_valid = sum(!is.na(pred_flux_soil))
  )
  
  cat(" done. Mean flux:", round(mean_soil_flux, 6), "μmol/m²/s\n")
}

# Calculate plot-basis soil flux
soil_results <- bind_rows(soil_fluxes_monthly) %>%
  mutate(
    Phi_soil_umol_m2_s = mean_flux * (geometry$A_soil / geometry$A_plot)
  )

cat("  Soil predictions complete.\n")

# =============================================================================
# 10.3: COMBINE RESULTS
# =============================================================================

cat("\nCombining tree and soil predictions...\n")

monthly_results <- monthly_results %>%
  left_join(tree_results %>% select(month, Phi_tree_umol_m2_s), by = "month") %>%
  left_join(soil_results %>% select(month, Phi_soil_umol_m2_s), by = "month") %>%
  mutate(
    Phi_plot_umol_m2_s = Phi_tree_umol_m2_s + Phi_soil_umol_m2_s
  )

# 7.5: Convert units to mg CH4 m-2 d-1
monthly_results <- monthly_results %>%
  mutate(
    Phi_tree_mg_m2_d = Phi_tree_umol_m2_s * 86400 * 16 * 1e-3,
    Phi_soil_mg_m2_d = Phi_soil_umol_m2_s * 86400 * 16 * 1e-3,
    Phi_plot_mg_m2_d = Phi_plot_umol_m2_s * 86400 * 16 * 1e-3
  )

print(monthly_results)

# Annual totals
annual_summary <- monthly_results %>%
  summarise(
    annual_tree_mg_m2 = sum(Phi_tree_mg_m2_d) * 30.4,  # Approximate days/month
    annual_soil_mg_m2 = sum(Phi_soil_mg_m2_d) * 30.4,
    annual_plot_mg_m2 = sum(Phi_plot_mg_m2_d) * 30.4
  )

cat("\nAnnual totals (mg CH4 m-2 yr-1):\n")
print(annual_summary)
cat("\nNote: Negative values indicate CH4 uptake (soil consumption)\n")

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

ggsave("QC_CHAMBER_COMPARISON.png", p1, width = 8, height = 6)

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

ggsave("QC_OUTLIERS.png", p1b, width = 6, height = 6)

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

ggsave("QC_PRED_VS_OBS_TREE.png", p2, width = 8, height = 6)

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

ggsave("QC_PRED_VS_OBS_SOIL.png", p2b, width = 8, height = 6)

# QC4: Moisture calibration R² by month
p4 <- ggplot(MOISTURE_AFFINE_TABLE, aes(x = month, y = R2_t)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = ifelse(!is.na(R2_t), round(R2_t, 2), "NA")), vjust = -0.5) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Moisture Calibration Quality by Month",
       x = "Month", y = "R²") +
  ylim(0, 1) +
  theme_minimal()

ggsave("QC_MOISTURE_CALIB.png", p4, width = 8, height = 6)

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

# Now sort and select top 20
tree_importance <- tree_importance[order(tree_importance$Importance, decreasing = TRUE), ]
tree_importance <- tree_importance[1:min(20, nrow(tree_importance)), ]

p5a <- ggplot(tree_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  labs(title = "Tree Model - Top 20 Feature Importances",
       subtitle = "Now including month_sin and month_cos",
       x = "Feature", y = "Importance") +
  theme_minimal()

ggsave("QC_TREE_IMPORTANCE.png", p5a, width = 10, height = 8)

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

ggsave("QC_SOIL_IMPORTANCE.png", p5b, width = 10, height = 6)

cat("✓ Quality control plots saved\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\nSaving outputs...\n")

# Save main results
write_csv(monthly_results, "MONTHLY_FLUXES.csv")
write_csv(annual_summary, "ANNUAL_SUMMARY.csv")

# Save diagnostic tables
write_csv(MOISTURE_AFFINE_TABLE, "MOISTURE_AFFINE_TABLE.csv")
write_csv(SI_TABLES, "SI_TABLES.csv")

# Save taxonomy priors
save(TAXONOMY_PRIORS, file = "TAXONOMY_PRIORS.RData")

# Save models
save(TreeRF, SoilRF, file = "RF_MODELS.RData")

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

save(diagnostics, file = "DIAGNOSTICS.RData")

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

save(coverage_report, file = "COVERAGE_REPORT.RData")

cat("\n=== WORKFLOW COMPLETE - GPS COORDINATES VERSION ===\n")
cat("All spatial data uses GPS coordinates (x = longitude, y = latitude)\n")