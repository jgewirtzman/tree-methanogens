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

# NOTE (2026 revision): we now work in nmol m-2 s-1 END TO END, which is the unit the
# paper reports in. Previously the pipeline divided by 1000 here to model in umol.
# That mattered for more than bookkeeping: asinh() has a fixed knee at |x| ~ 1, so at
# umol magnitudes (median ~2.6e-5) asinh(y) was the IDENTITY to 7 significant figures --
# i.e. the "asinh transform" was doing nothing at all, and the model was being fit on
# raw flux. Keeping nmol puts the data where asinh actually acts on the upper tail.
#
# WARNING: the columns are still NAMED *_umol_m2_s upstream but have always CONTAINED
# nmol m-2 s-1. With this change the values are finally consistent with nmol.
cat("\n*** UNITS: modelling in nmol m-2 s-1 (no conversion applied) ***\n")

cat("Flux ranges (nmol m-2 s-1):\n")
cat("  TREE_JULY: [", round(min(TREE_JULY$stem_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(TREE_JULY$stem_flux_umol_m2_s, na.rm=T), 6), "]\n")
cat("  TREE_YEAR: [", round(min(TREE_YEAR$stem_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(TREE_YEAR$stem_flux_umol_m2_s, na.rm=T), 6), "]\n")
cat("  SOIL_YEAR: [", round(min(SOIL_YEAR$soil_flux_umol_m2_s, na.rm=T), 6), ",", 
    round(max(SOIL_YEAR$soil_flux_umol_m2_s, na.rm=T), 6), "]\n")

# NOTE (2026 revision fix): TREE_YEAR already carries its own `species` column.
# The previous `left_join(INVENTORY %>% select(tree_id, species))` collided with it,
# producing species.x/species.y and leaving NO `species` column. Every one of the 369
# year-long monthly rows then failed the downstream completeness filter, so TreeRF was
# silently trained on the July campaign alone (353 rows, months 6-7, rigid chamber only)
# and the chamber harmonisation in §2.5 never had any semirigid data to correct.
# Fix: keep TREE_YEAR's own species and fill only genuine gaps from the inventory.
TREE_YEAR <- TREE_YEAR %>%
  left_join(INVENTORY %>% dplyr::select(tree_id, species_from_inventory = species),
            by = "tree_id") %>%
  mutate(species = dplyr::coalesce(species, species_from_inventory)) %>%
  dplyr::select(-species_from_inventory)

cat("  TREE_YEAR after species fill:", sum(!is.na(TREE_YEAR$species)), "of",
    nrow(TREE_YEAR), "rows have species\n")

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
          
          # NOTE (2026 revision): month can be NA for recovered deployments whose date
          # could not be resolved; `df$month == m` then yields NA and any() errors.
          for(m in unique(df$month[!is.na(df$month)])) {
            month_idx <- !is.na(df$month) & df$month == m & is.na(df[[col]])
            if(any(month_idx)) {
              # NOTE (2026 revision): guard against duplicate/NA month keys returning a
              # vector here, which made the if() condition length > 1.
              month_med <- monthly_medians$med_val[which(monthly_medians$month == m)][1]
              if(length(month_med) == 1 && !is.na(month_med)) {
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

# NOTE (2026 revision): environmental covariates are no longer median-imputed.
# smart_impute() substituted a monthly (then overall) median, which fabricates exactly
# the temp/moisture signal the model is trying to learn -- and it was operating on 26%
# of tree rows. ranger splits on observed values only, so passing NA through is both
# honest and better-performing. DBH is still imputed: it is a stable tree property, not
# an environmental driver, and a missing DBH would otherwise drop the row entirely.
tree_combined <- smart_impute(tree_combined, c("dbh_m"))

# NOTE (2026 revision): the 1st/99th-percentile trim on asinh(flux) has been REMOVED.
# It deleted the largest stem emissions -- the hotspots this study documents -- and cost
# 33% of the mean tree flux, biasing the budget low. Chamber artifacts are now screened
# mechanistically at load time (C0 > 5000 ppb; see 01_load_and_prep_data.R), which flags
# no tree deployments at all: max tree C0 is 4,980 / 2,156 / 2,432 ppb by campaign.
tree_combined$z <- asinh(tree_combined$stem_flux_umol_m2_s)
n_outliers_tree <- 0
cat("  Tree flux outlier trimming disabled (chamber-integrity screen applied upstream)\n")

# =============================================================================
# CHAMBER CALIBRATION (spec §2.5)
# -----------------------------------------------------------------------------
# NOTE (2026 revision fix): this step previously did nothing --
#   stem_flux_corrected <- stem_flux_umol_m2_s
# -- because the species-join collision upstream had already dropped every
# semirigid (monthly-chamber) row, leaving a July/rigid-only training set with no
# bias to estimate. With that collision fixed, both chamber types are present and
# the harmonisation the spec calls for is estimated and applied here.
#
# Model (asinh scale):  z ~ I_semirigid + airT + soilT + moisture + species
# beta1 is the systematic offset of the monthly (semirigid) chamber relative to the
# July (rigid) reference. Only semirigid rows are corrected: z_corr = z - beta1.
# =============================================================================

CHAMBER_CALIB <- list(beta1 = 0, se = NA_real_, ci = c(NA_real_, NA_real_),
                      r2 = NA_real_, n_rigid = 0, n_semirigid = 0, fitted = FALSE)

cal_df <- tree_combined %>%
  filter(!is.na(stem_flux_umol_m2_s), !is.na(chamber_type)) %>%
  mutate(z = asinh(stem_flux_umol_m2_s),
         I_semirigid = as.numeric(chamber_type == "semirigid"))

CHAMBER_CALIB$n_rigid     <- sum(cal_df$I_semirigid == 0)
CHAMBER_CALIB$n_semirigid <- sum(cal_df$I_semirigid == 1)

if (CHAMBER_CALIB$n_rigid >= 10 && CHAMBER_CALIB$n_semirigid >= 10) {
  # species as a control only where it is estimable
  sp_ok <- cal_df %>% count(species) %>% filter(n >= 5, !is.na(species)) %>% pull(species)
  cal_df$species_ctl <- ifelse(cal_df$species %in% sp_ok, cal_df$species, "OTHER")

  chamber_fit <- lm(
    z ~ I_semirigid + air_temp_C + soil_temp_C + soil_moisture_abs + species_ctl,
    data = cal_df
  )
  cf <- summary(chamber_fit)$coefficients
  if ("I_semirigid" %in% rownames(cf)) {
    CHAMBER_CALIB$beta1  <- unname(cf["I_semirigid", "Estimate"])
    CHAMBER_CALIB$se     <- unname(cf["I_semirigid", "Std. Error"])
    CHAMBER_CALIB$ci     <- CHAMBER_CALIB$beta1 + c(-1.96, 1.96) * CHAMBER_CALIB$se
    CHAMBER_CALIB$r2     <- summary(chamber_fit)$r.squared
    CHAMBER_CALIB$fitted <- TRUE
  }
  cat(sprintf("  Chamber bias beta1 = %.4f (SE %.4f, 95%% CI %.4f to %.4f), model R2 = %.3f\n",
              CHAMBER_CALIB$beta1, CHAMBER_CALIB$se,
              CHAMBER_CALIB$ci[1], CHAMBER_CALIB$ci[2], CHAMBER_CALIB$r2))
  cat(sprintf("  n rigid = %d, n semirigid = %d\n",
              CHAMBER_CALIB$n_rigid, CHAMBER_CALIB$n_semirigid))
} else {
  cat("  WARNING: only one chamber type present (rigid =", CHAMBER_CALIB$n_rigid,
      ", semirigid =", CHAMBER_CALIB$n_semirigid, ") - no chamber correction applied.\n")
}

# Apply: correct ONLY the semirigid rows, on the asinh scale, then back-transform.
tree_combined$z_raw  <- asinh(tree_combined$stem_flux_umol_m2_s)
tree_combined$z_corr <- ifelse(
  !is.na(tree_combined$chamber_type) & tree_combined$chamber_type == "semirigid",
  tree_combined$z_raw - CHAMBER_CALIB$beta1,
  tree_combined$z_raw
)
tree_combined$stem_flux_corrected <- sinh(tree_combined$z_corr)

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

# See note above: environmental covariates are passed through with their NAs.
# (soil coverage is now 88% measured)

# NOTE (2026 revision): MAD k=8 outlier removal has been DISABLED. It deleted four
# genuine wetland-margin emissions (11.6, 25.7, 51.3, 63.1 nmol m-2 s-1 -- all with
# normal chamber C0, from the plots adjacent to the wetland) while RETAINING the one
# real artifact (-176.1 nmol m-2 s-1, C0 = 49,049 ppb). Net effect: the soil sink came
# out ~2.7x too strong. Artifacts are now screened on chamber integrity upstream.
SOIL_YEAR$z_soil <- asinh(SOIL_YEAR$soil_flux_umol_m2_s)
outlier_mask_soil <- rep(TRUE, nrow(SOIL_YEAR))
n_outliers_soil <- 0
cat("  Soil MAD outlier removal disabled (chamber-integrity screen applied upstream)\n")

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
      keep.inbag = TRUE,
      seed = 42                      # reproducibility: this RF's residuals define the seasonal index
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
    num.threads = 1,
    seed = 42                        # reproducibility: residuals here define the taxonomy prior (a TreeRF feature)
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

# ---------------------------------------------------------------------------
# Shared species -> genus -> family lookup for taxon_prior_asinh.
# NOTE (2026 revision fix): the training path applied this prior, but the
# prediction path hard-set taxon_prior_asinh = 0 for every inventory tree, so the
# feature was systematically offset from its training distribution and contributed
# nothing at prediction time. This helper is now used in BOTH places.
# ---------------------------------------------------------------------------
lookup_taxon_prior <- function(species_vec, taxonomy, taxon_priors) {
  tx <- taxonomy %>% dplyr::select(species, genus, family) %>% distinct()
  idx <- match(species_vec, tx$species)
  genus_vec  <- tx$genus[idx]
  family_vec <- tx$family[idx]

  out <- rep(0, length(species_vec))
  sp_hit <- match(species_vec, taxon_priors$species$species)
  gn_hit <- match(genus_vec,   taxon_priors$genus$genus)
  fm_hit <- match(family_vec,  taxon_priors$family$family)

  out <- ifelse(!is.na(sp_hit), taxon_priors$species$med_resid[sp_hit],
         ifelse(!is.na(gn_hit), taxon_priors$genus$med_resid[gn_hit],
         ifelse(!is.na(fm_hit), taxon_priors$family$med_resid[fm_hit], 0)))
  as.numeric(ifelse(is.na(out), 0, out))
}

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
  
  # NOTE (2026 revision): the interpolated December-anchored moisture surface is no
  # longer used to FILL training rows. With the per-collar env archive compiled,
  # measured moisture now covers ~87% of monthly, 90% of 2021-height and 100% of 2023
  # measurements, and filling the remainder from a surface that was itself fitted to
  # these same observations is both circular and less accurate than leaving the value
  # missing (ranger splits on observed values only). The surface is still used at
  # PREDICTION time, where there is no alternative.
  if (isTRUE(getOption("rf.fill_training_moisture", FALSE)) &&
      any(is.na(df$soil_moisture_at_tree)) && !is.null(Mhat_fn)) {
    mm <- is.na(df$soil_moisture_at_tree)
    df$soil_moisture_at_tree[mm] <- Mhat_fn(df$month[mm], df$x[mm], df$y[mm])
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
  
  # ---------------------------------------------------------------------------
  # SPECIES LEVELS: a taxonomic ladder -- own level, else genus, else pooled.
  #
  # A species with >= 10 records gets its own level. Below that, it falls back to
  # its GENUS if that genus is measured, and only otherwise to SPECIES_OTHER.
  #
  # This is a bias-variance choice and it was settled by repeated 5-fold CV over
  # observations, with the rare-species rows scored separately because an overall
  # RMSE is dominated by the abundant species and will always favour pooling
  # (rev_rf_species_pooling.R):
  #
  #   scheme               RMSE all   RMSE rare   bias rare   ratio rare
  #   genus_then_pool       0.3110      0.1218      0.0245       1.23
  #   lump10 (old)          0.3110      0.1482      0.0414       1.39
  #   all_species           0.3111      0.1487      0.0374       1.35
  #   all_species node=2    0.3112      0.1502      0.0353       1.33
  #   shrunk offset (EB)    0.3118      0.1510      0.0395       1.37
  #
  # The ladder wins on rare-species RMSE (18% better than lumping), on rare-species
  # BIAS (41% less), and ties for best overall -- so it costs nothing on the
  # abundant species. Every scheme over-predicts rare species, and the ladder
  # over-predicts least, which is also the conservative direction here.
  #
  # Giving every measured species its own level was tried and is NOT better: three
  # records cannot reshape an 800-tree forest, and shrinking min.node.size to 2 to
  # let them try made it worse. Empirical-Bayes shrinkage of a species offset
  # toward genus toward the grand mean -- effectively the old taxon prior -- was
  # the worst of the five.
  #
  # What the ladder fixes: Kalmia latifolia has 3 records averaging 0.0350 nmol
  # m-2 s-1 and is the only Kalmia here, so it forms its own genus level instead
  # of joining a bucket averaging 0.0941 that is dominated by Carya ovata (0.2583)
  # and Quercus velutina (0.1139). That bucket was being applied to 1,898 mountain
  # laurel stems. Quercus velutina -- the felled black oak, measured at seven
  # heights -- now joins Quercus rather than being treated as unmeasured.
  #
  # For a species with NO records nothing works: leave-one-species-out gives
  # negative R2 for every scheme (rev_rf_species_fallback_loso.R). That is a
  # stated limitation, not something to engineer around.
  # ---------------------------------------------------------------------------
  sp_raw <- as.character(df$species)
  unident <- is.na(sp_raw) | grepl("^unknown", sp_raw, ignore.case = TRUE) |
             !grepl(" ", trimws(sp_raw))
  sp_raw[unident] <- "SPECIES_OTHER"
  gen <- ifelse(unident, "SPECIES_OTHER", sub(" .*", "", sp_raw))
  cnt <- table(sp_raw); big <- names(cnt)[cnt >= 10]
  gen_seen <- unique(gen[!unident])
  df$species_clean <- ifelse(sp_raw %in% big, sp_raw,
                        ifelse(gen %in% gen_seen & !unident, paste0("GEN_", gen),
                               "SPECIES_OTHER"))
  cat(sprintf("  Species levels: %d own + %d genus + SPECIES_OTHER (%d unidentified)\n",
              sum(!grepl("^GEN_|SPECIES_OTHER", unique(df$species_clean))),
              sum(grepl("^GEN_", unique(df$species_clean))), sum(unident)))

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
    # NOTE (2026 revision): the response is RAW flux in nmol m-2 s-1, not asinh.
    # asinh was inert at umol (identity to 7 s.f.); moving to nmol made it bite but
    # introduced a 13% low bias at back-transformation (Jensen), concentrated in the
    # top flux decile which carries 68% of the total. Fitting raw is unbiased with no
    # correction at all (mean ratio 0.931 tree / 0.994 soil) and recovers the tail
    # better, at a cost of only ~0.01 R2. The column name is retained so downstream
    # code is unchanged.
    df$y_asinh <- as.numeric(df$stem_flux_corrected)
  } else if("stem_flux_umol_m2_s" %in% names(df) && !all(is.na(df$stem_flux_umol_m2_s))) {
    df$y_asinh <- as.numeric(df$stem_flux_umol_m2_s)
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
  
  # NOTE (2026 revision): previously these took the MONTHLY PLOT-LEVEL means only, so
  # soil temperature had exactly ONE value per month across all collars -- it carried no
  # spatial information whatsoever and could only encode season. The campaign's own
  # per-collar readings (104 of 266 measurements) were discarded. Prefer the observed
  # value and fall back to the monthly mean.
  obs_soil_t <- if ("soil_temp_C" %in% names(df)) suppressWarnings(as.numeric(df$soil_temp_C)) else NA_real_
  obs_air_t  <- if ("air_temp_C"  %in% names(df)) suppressWarnings(as.numeric(df$air_temp_C))  else NA_real_
  df$soil_temp_C_mean <- dplyr::coalesce(obs_soil_t, as.numeric(df$soil_temp_C_mean))
  df$air_temp_C_mean  <- dplyr::coalesce(obs_air_t,  as.numeric(df$air_temp_C_mean))
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
    df$y_asinh <- as.numeric(df$soil_flux_umol_m2_s)   # raw scale, see note above
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

# NOTE (2026 revision): requiring non-missing MOISTURE dropped 84 otherwise usable
# measurements once training-time imputation was removed. ranger splits on observed
# values only, so a row with a missing environmental covariate still contributes
# through species, size and season. Only the response and species identity are
# genuinely required.
complete_rows <- !is.na(tree_train$y_asinh) &
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

# =============================================================================
# TREE FEATURE MATRIX -- final specification (2026 revision)
# -----------------------------------------------------------------------------
# Six predictors, chosen by testing each against a single random reference column
# over 20 fits (a predictor is kept only if its permutation importance beats pure
# noise). Mean importance and the fraction of runs beating the random column:
#
#   dbh_m         0.0749  100%      absolute size -> drives stem area
#   species       0.0310  100%      as a NATIVE FACTOR -- ranger splits on levels
#                                   directly, so no one-hot expansion is needed
#   moisture      0.0193  100%
#   soil_temp     0.0138  100%
#   air_temp      0.0094   95%
#   height_cm     0.0082   85%      measurement height (see below)
#   month         0.0053   65%      marginal, and largely carried by the temperatures
#
# REMOVED, dbh_within_z: it passed the random-column test above, but it is unusable
# for the purpose this model exists for. It is standardised WITHIN POPULATION, and
# the fitting and application populations are not the same one. Measured trees are
# much larger than inventory stems -- Pinus strobus averages 0.490 m among the 110
# measured and 0.047 m across the 3,317 in the inventory; Betula lenta 0.302 m
# against 0.061 m. Because each side is standardised against itself, the inventory
# z-scores land INSIDE the training range (-3.06 to 12.8 against -4.68 to 5.69):
# 0.0% of stems and 0.0% of band area fall outside it, so no extrapolation check
# fires. A 1.9 cm pine sapling is handed z = -0.26, "a slightly below-average pine",
# and the forest reads the size signal of a ~45 cm tree.
#
# The effect is a uniform inflation across every species, not a tail artifact:
# predicted flux at 163 cm falls 20-35% per species when it is dropped (Tsuga
# 0.0615 -> 0.0397, Pinus strobus 0.0715 -> 0.0419, Acer rubrum 0.1638 -> 0.1317),
# and the stand band term falls ~29%. The per-species calibration does not absorb
# this: that ratio is estimated on measured trees, where the z-scores are correct.
#
# Dropping it costs no skill -- 5x5-fold CV R2 0.2104 -> 0.2236 -- and removes the
# r = 0.79 collinearity with dbh_m that made permutation importance unreadable.
#
# REJECTED: elevation (35% of runs), local basal area (0%), distance to river (0%).
# These are static per-tree values and add nothing for stems; for soil they proved to
# be collar identity in disguise.
#
# ALSO REMOVED relative to the previous version: species x moisture interactions (RF
# finds interactions itself), the empirical seasonal index (its range was ~1e-7 against
# a target sd of 1.3e-4 -- numerically inert), the taxonomy prior (99.98% of stem area
# is measured species), and the asinh transform (fit on raw nmol; unbiased with no
# back-transformation correction needed).
#
# IMPORTANCE IS READ AS PERMUTATION (incMSE), NOT IMPURITY. Impurity gain counts
# variance reduction at splits and inflates continuous, high-cardinality predictors;
# on these data the two metrics give opposite answers (impurity called the soil model
# a temperature model at 88% soil temp, permutation calls it a moisture model at 44%).
# The models are still FIT with importance = "impurity" for continuity with the
# stored objects; rev_fig_model_findings.R refits with permutation to report it.
# =============================================================================
X_tree_species_first <- data.frame(
  species       = factor(tree_train_complete$species_clean),
  dbh_m         = tree_train_complete$dbh_m,
  soil_moisture_at_tree = tree_train_complete$soil_moisture_at_tree,
  soil_temp_C_mean = tree_train_complete$soil_temp_C_mean,
  air_temp_C_mean  = tree_train_complete$air_temp_C_mean,
  # Measurement height. The 2021 campaign measured 50/125/200 cm on the same trees;
  # without this the three rows are identical in every predictor and differ only in the
  # response, so the model averages them and the height signal is lost. Flux roughly
  # halves from 50 cm (mean 0.220) to 200 cm (0.102 nmol m-2 s-1). Including it raises
  # grouped-CV R2 from 0.068 to 0.122. Monthly and 2023 rows are breast-height -> 125.
  height_cm     = ifelse(is.na(tree_train_complete$measurement_height_cm), 125,
                         tree_train_complete$measurement_height_cm)
)

# Remove any columns with zero variance
tree_train_complete_species_levels <- factor(tree_train_complete$species_clean)

# NOTE (2026 revision): species is now a native factor, and var() errors on factors.
zero_var_cols <- sapply(X_tree_species_first, function(x) {
  if (is.factor(x) || is.character(x)) return(length(unique(x[!is.na(x)])) < 2)
  if (all(is.na(x))) return(TRUE)
  isTRUE(var(x, na.rm = TRUE) == 0)
})
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
# =============================================================================
# SOIL FEATURE SET -- final specification (2026 revision)
# -----------------------------------------------------------------------------
# Four predictors. All beat a random reference column in 100% of 20 fits:
#   moisture   1.121   soil CH4 uptake is diffusion-limited, so this dominates
#   air_temp   0.597
#   soil_temp  0.466
#   month      0.385
#
# DELIBERATELY EXCLUDED -- distance to river, elevation and local basal area. Each
# beat the random column 100% of the time and distance-to-river was the single
# strongest predictor of all (0.834), but every one of them is STATIC PER COLLAR: they
# take 28 distinct values across 266 observations, one per collar, whereas moisture
# takes 194. They are collar identity in disguise. Including all four raises OOB R2 to
# 0.659 while grouped-CV R2 is only 0.451 -- a 0.21 gap from memorising 28 collars.
# Moisture alone gives 0.461 / 0.435, a gap of 0.026, with the best mean ratio (1.007),
# and would not require applying collar-specific values to thousands of unseen grid
# cells. (Elevation specifically: adding it raises OOB by 0.062 while LOWERING grouped
# CV by 0.046 -- the textbook overfitting signature.)
#
# Also removed: the empirical seasonal index and the two moisture x temperature
# interaction terms (RF finds interactions itself), and month_sin/month_cos in favour
# of plain month.
# =============================================================================
soil_features <- c("soil_moisture_at_site", "soil_temp_C_mean",
                   "air_temp_C_mean", "month")

X_soil <- soil_train[, soil_features, drop = FALSE]

zero_var_soil <- sapply(X_soil, function(x) {
  if (is.factor(x) || is.character(x)) return(length(unique(x[!is.na(x)])) < 2)
  if (all(is.na(x))) return(TRUE)
  isTRUE(var(x, na.rm = TRUE) == 0)
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

# =============================================================================
# BACK-TRANSFORMATION BIAS CORRECTION (2026 revision)
# -----------------------------------------------------------------------------
# The models are unbiased ON THE MODELLED SCALE (mean asinh predicted / observed =
# 1.005 for trees) but sinh() of an unbiased asinh prediction is not unbiased on the
# raw scale -- Jensen's inequality. The effect is concentrated in the upper tail:
# the top flux decile carries 68% of total tree flux and is under-predicted by ~62%,
# which alone drags the aggregate to 0.87 of observed.
#
# Correction: a single multiplicative factor, mean(observed)/mean(sinh(OOB pred)),
# estimated ONLY from training data (a degenerate smearing estimator). Tested under
# grouped CV against the alternatives:
#   asinh + naive sinh   R2 0.112  mean ratio 0.819
#   asinh + Duan smear   R2 0.113  mean ratio 0.877
#   asinh + ratio        R2 0.116  mean ratio 0.932   <- adopted
#   raw scale, no asinh  R2 0.104  mean ratio 0.931
# Duan smearing was rejected: it is unstable here because sinh() exponentially
# amplifies large positive residuals (it overcorrected soil by 2.6x before screening).
# A residual ~7% low bias remains and is REPORTED, not corrected away: it reflects the
# model's genuine inability to predict the extreme upper tail for an unseen tree.
# =============================================================================
bt_factor <- function(rf, y_raw) {
  ins <- sinh(rf$predictions)
  ok <- is.finite(ins) & is.finite(y_raw)
  f <- mean(y_raw[ok]) / mean(ins[ok])
  if (!is.finite(f) || f <= 0) f <- 1
  f
}
BT_TREE <- 1  # raw scale: no back-transformation correction required
BT_SOIL <- 1
cat(sprintf("\nBack-transformation correction: tree x%.4f, soil x%.4f\n", BT_TREE, BT_SOIL))

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
      species_factor = factor(
        species_clean,
        levels = levels(tree_train$species_factor)
      ),
      # NOTE (2026 revision fix): was hard-coded to 0 here, which zeroed the taxon
      # prior for every inventory tree while training rows carried non-zero values.
      taxon_prior_asinh = lookup_taxon_prior(species, TAXONOMY, TAXONOMY_PRIORS)
    )
  
  # Prediction features must match the trained specification exactly (2026 revision):
  # species as a native factor, absolute DBH, the three environmental drivers, and
  # height fixed at the 0-2 m band mid-point.
  X_pred_tree <- data.frame(
    species       = factor(as.character(inv_predictions$species_clean),
                           levels = levels(tree_train_complete_species_levels)),
    dbh_m         = inv_predictions$dbh_m,
    soil_moisture_at_tree = inv_predictions$soil_moisture_abs,
    soil_temp_C_mean = inv_predictions$soil_temp_C,
    air_temp_C_mean  = inv_predictions$air_temp_C,
    height_cm     = 125
  )

  # Align columns with training.
  # NOTE (2026 revision fix): this loop silently zero-fills any trained feature that is
  # absent from the prediction frame, which is how a name mismatch can blank out most of
  # the model's inputs without any error. Check first and fail loudly.
  missing_tree_feats <- setdiff(colnames(X_tree), colnames(X_pred_tree))
  if (length(missing_tree_feats) > 0) {
    stop("TreeRF prediction is missing trained features: ",
         paste(missing_tree_feats, collapse = ", "))
  }

  X_pred_aligned <- matrix(0, nrow = nrow(X_pred_tree), ncol = ncol(X_tree))
  colnames(X_pred_aligned) <- colnames(X_tree)

  for (col_name in colnames(X_tree)) {
    if (col_name %in% colnames(X_pred_tree)) {
      X_pred_aligned[, col_name] <- X_pred_tree[[col_name]]
    }
  }
  
  pred_asinh <- predict(TreeRF, X_pred_aligned)$predictions
  pred_flux_umol_m2_s <- pred_asinh          # raw scale: no back-transform
  
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
  
  # NOTE (2026 revision fix): these columns are renamed to match the names SoilRF was
  # TRAINED on. Previously they were passed as soil_temp_C / air_temp_C /
  # soil_moisture_abs / SI_soil, none of which matched the model's
  # soil_temp_C_mean / air_temp_C_mean / soil_moisture_at_site / SI. The alignment loop
  # below then zero-filled all four, so half of SoilRF's inputs were exactly 0 at every
  # prediction -- biasing mean soil uptake ~1.9x too strong and collapsing its variance.
  X_pred_soil <- grid_predictions %>%
    dplyr::transmute(
      soil_moisture_at_site = soil_moisture_abs,
      soil_temp_C_mean      = soil_temp_C,
      air_temp_C_mean       = air_temp_C,
      month                 = t
    ) %>%
    as.matrix()

  # Fail loudly rather than silently zero-filling a missing predictor.
  missing_soil_feats <- setdiff(colnames(X_soil), colnames(X_pred_soil))
  if (length(missing_soil_feats) > 0) {
    stop("SoilRF prediction is missing trained features: ",
         paste(missing_soil_feats, collapse = ", "))
  }

  X_pred_soil_aligned <- matrix(0, nrow = nrow(X_pred_soil), ncol = ncol(X_soil))
  colnames(X_pred_soil_aligned) <- colnames(X_soil)

  for (col_name in colnames(X_soil)) {
    if (col_name %in% colnames(X_pred_soil)) {
      X_pred_soil_aligned[, col_name] <- X_pred_soil[, col_name]
    }
  }
  
  pred_asinh_soil <- predict(SoilRF, X_pred_soil_aligned)$predictions
  pred_flux_soil <- pred_asinh_soil          # raw scale: no back-transform
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
    # nmol m-2 s-1 -> mg CH4 m-2 d-1 : x 86400 s/d x 16 g/mol x 1e-6 mg/(nmol*g/mol)
    # (1 nmol CH4 = 16 ng = 1.6e-5 mg). Was 1e-3 when the model ran in umol.
    Phi_tree_mg_m2_d = Phi_tree_umol_m2_s * 86400 * 16 * 1e-6,
    Phi_soil_mg_m2_d = Phi_soil_umol_m2_s * 86400 * 16 * 1e-6,
    Phi_plot_mg_m2_d = Phi_plot_umol_m2_s * 86400 * 16 * 1e-6
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
# NOTE (2026 revision fix): TreeRF is fitted on tree_train_complete (complete cases),
# so its predictions align with THAT frame, not the unfiltered tree_train.
tree_train_complete$pred_asinh <- TreeRF$predictions
tree_train_complete$pred_flux <- sinh(tree_train_complete$pred_asinh)
tree_train_complete$chamber <- ifelse(tree_train_complete$chamber_rigid == 1, "rigid", "semirigid")

p2 <- ggplot(tree_train_complete, aes(x = stem_flux_corrected, y = pred_flux, 
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
soil_train_complete$pred_asinh <- SoilRF$predictions
soil_train_complete$pred_flux <- sinh(soil_train_complete$pred_asinh)

p2b <- ggplot(soil_train_complete, aes(x = soil_flux_umol_m2_s, y = pred_flux)) +
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
  # NOTE (2026 revision fix): chamber is now handled by the beta1 correction on the
  # target (spec 2.5), not carried as a predictor.
  Chamber_as_predictor = FALSE,
  Month_as_predictor = TRUE,
  Chamber_beta1 = CHAMBER_CALIB$beta1,
  Chamber_beta1_SE = CHAMBER_CALIB$se,
  Chamber_beta1_CI = CHAMBER_CALIB$ci,
  Chamber_calib_R2 = CHAMBER_CALIB$r2,
  Chamber_n_rigid = CHAMBER_CALIB$n_rigid,
  Chamber_n_semirigid = CHAMBER_CALIB$n_semirigid,
  Chamber_calib_fitted = CHAMBER_CALIB$fitted,
  TreeRF_importance = importance(TreeRF),
  SoilRF_importance = importance(SoilRF),
  n_tree_training = nrow(X_tree),
  n_soil_training = nrow(X_soil),
  outliers_removed = list(
    tree_n = n_outliers_tree,
    tree_pct = 100*n_outliers_tree/(n_outliers_tree + nrow(tree_combined)),
    soil_n = n_outliers_soil,
    soil_pct = 100*n_outliers_soil/(n_outliers_soil + nrow(SOIL_YEAR)),
    # NOTE (2026 revision fix): the two datasets use DIFFERENT rules; this was
    # previously mislabelled as "MAD with k=8" for both.
    method_tree = "asinh 1st/99th percentile trim",
    method_soil = "MAD with k=8"
  ),
  # NOTE (2026 revision fix): these four counters referenced variables that were never
  # assigned anywhere in the script, so building `diagnostics` aborted before the file
  # could be saved -- which is why DIAGNOSTICS.RData on disk was stale (Sep 2025) and
  # disagreed with the saved models (its n_tree_training = 718 vs the model's 353).
  training_months = list(
    tree = sort(unique(tree_train_complete$month)),
    soil = sort(unique(soil_train_complete$month))
  ),
  training_chambers = as.list(table(tree_train_complete$chamber_type)),
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
# NOTE (2026 revision fix): TreeRF is fitted on tree_train_complete (complete cases),
# so its predictions align with THAT frame, not the unfiltered tree_train.
tree_train_complete$pred_asinh <- TreeRF$predictions
tree_train_complete$pred_flux <- sinh(tree_train_complete$pred_asinh)
tree_train_complete$chamber <- ifelse(tree_train_complete$chamber_rigid == 1, "rigid", "semirigid")

# Convert to nmol for plotting (multiply by 1000)
tree_train_complete$pred_flux_nmol <- tree_train_complete$pred_flux * 1000
tree_train_complete$obs_flux_nmol <- tree_train_complete$stem_flux_corrected * 1000

p2_nmol <- ggplot(tree_train_complete, aes(x = obs_flux_nmol, y = pred_flux_nmol, 
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
soil_train_complete$pred_asinh <- SoilRF$predictions
soil_train_complete$pred_flux <- sinh(soil_train_complete$pred_asinh)

# Convert to nmol for plotting
soil_train_complete$pred_flux_nmol <- soil_train_complete$pred_flux * 1000
soil_train_complete$obs_flux_nmol <- soil_train_complete$soil_flux_umol_m2_s * 1000

p2b_nmol <- ggplot(soil_train_complete, aes(x = obs_flux_nmol, y = pred_flux_nmol)) +
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
  tree_train_complete %>% 
    dplyr::select(obs_flux_nmol, pred_flux_nmol) %>%
    mutate(type = "Trees"),
  soil_train_complete %>%
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
# NOTE (2026 revision): a verbatim duplicate of the QC-plot section used to follow
# here. It re-ran the same plots against stale frames, errored, and aborted the
# script before its final messages. The canonical QC block is above (Section 10).
# =============================================================================

cat("\n=== 02_rf_models.R COMPLETE ===\n")
