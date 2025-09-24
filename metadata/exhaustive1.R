# Clean Exhaustive Model Search - Fixed Exclusions
# Load required packages first:
library(randomForest)
library(mgcv)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

prepare_search_data_incomplete <- function(data_subset, target_var, exclude_vars = NULL) {
  model_data <- data_subset %>%
    filter(!is.na(!!sym(target_var))) %>%
    select(-tree_id, -any_of(exclude_vars)) %>%
    mutate(species = if("species_id" %in% names(.)) as.factor(species_id) else NULL) %>%
    select(-any_of("species_id")) %>%
    mutate_if(is.character, as.numeric) %>%
    select(all_of(target_var), everything())
  
  # Remove variables with too many missing values (>70%)
  missing_threshold <- 0.7
  complete_vars <- model_data %>%
    summarise_all(~mean(is.na(.))) %>%
    select_if(~.x <= missing_threshold) %>%
    names()
  
  model_data <- model_data %>%
    select(all_of(complete_vars)) %>%
    filter(!is.na(!!sym(target_var)))
  
  return(model_data)
}

rf_iterative_search_incomplete <- function(model_data, target_var, max_vars = 20) {
  cat("=== RANDOM FOREST ITERATIVE SEARCH ===\n")
  
  predictors <- setdiff(names(model_data), target_var)
  cat("Starting with", length(predictors), "predictors\n")
  
  # Full model for importance ranking
  cat("Building full RF model for variable importance...\n")
  rf_full <- randomForest(
    as.formula(paste(target_var, "~ .")), 
    data = model_data,
    ntree = 500,
    importance = TRUE,
    na.action = na.roughfix
  )
  
  # Get importance ranking
  importance_df <- importance(rf_full) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    select(variable, everything()) %>%
    arrange(desc(`%IncMSE`))
  
  cat("Top 5 variables:", paste(head(importance_df$variable, 5), collapse = ", "), "\n")
  
  # Test models with top k variables
  test_rf_subset <- function(k) {
    if (k > nrow(importance_df)) k <- nrow(importance_df)
    
    top_vars <- importance_df$variable[1:k]
    
    tryCatch({
      set.seed(123)
      train_idx <- sample(1:nrow(model_data), 0.7 * nrow(model_data))
      train_data <- model_data[train_idx, ]
      test_data <- model_data[-train_idx, ]
      
      rf_model <- randomForest(
        as.formula(paste(target_var, "~", paste(top_vars, collapse = " + "))), 
        data = train_data,
        ntree = 500,
        na.action = na.roughfix
      )
      
      if (nrow(test_data) > 0) {
        predictions <- predict(rf_model, test_data)
        actual <- test_data[[target_var]]
        
        complete_cases <- !is.na(predictions) & !is.na(actual)
        if (sum(complete_cases) > 0) {
          test_rmse <- sqrt(mean((predictions[complete_cases] - actual[complete_cases])^2))
          test_r2 <- cor(predictions[complete_cases], actual[complete_cases])^2
        } else {
          test_rmse <- NA
          test_r2 <- NA
        }
      } else {
        test_rmse <- NA
        test_r2 <- NA
      }
      
      tibble(
        n_vars = k,
        variables = paste(top_vars, collapse = ", "),
        oob_mse = tail(rf_model$mse, 1),
        oob_rsq = tail(rf_model$rsq, 1),
        test_rmse = test_rmse,
        test_r2 = test_r2
      )
    }, error = function(e) {
      tibble(
        n_vars = k,
        variables = paste(top_vars, collapse = ", "),
        oob_mse = NA,
        oob_rsq = NA,
        test_rmse = NA,
        test_r2 = NA
      )
    })
  }
  
  # Test different numbers of variables
  k_values <- c(1:min(max_vars, nrow(importance_df)))
  cat("Testing", length(k_values), "different variable combinations...\n")
  
  results <- map_dfr(k_values, test_rf_subset)
  cat("RF search completed!\n")
  
  return(list(
    results = results,
    importance = importance_df,
    full_model = rf_full
  ))
}

gam_specification_search <- function(model_data, target_var, max_vars = 10) {
  cat("=== GAM SPECIFICATION SEARCH ===\n")
  
  predictors <- setdiff(names(model_data), target_var)
  results_list <- list()
  
  # Use RF importance for variable selection
  rf_temp <- randomForest(as.formula(paste(target_var, "~ .")), 
                          data = model_data, importance = TRUE)
  imp_order <- rownames(importance(rf_temp))[order(importance(rf_temp)[,1], decreasing = TRUE)]
  
  # Test different numbers of top variables
  for(i in 1:min(max_vars, length(imp_order))) {
    vars <- imp_order[1:i]
    
    tryCatch({
      formula_str <- paste(target_var, "~", paste(vars, collapse = " + "))
      gam_model <- gam(as.formula(formula_str), data = model_data)
      
      results_list[[length(results_list) + 1]] <- tibble(
        model_type = "linear_gam",
        n_vars = length(vars),
        variables = paste(vars, collapse = ", "),
        aic = AIC(gam_model),
        deviance_explained = summary(gam_model)$dev.expl,
        r_squared = summary(gam_model)$r.sq
      )
    }, error = function(e) NULL)
  }
  
  cat("GAM search completed!\n")
  
  results_df <- bind_rows(results_list) %>%
    arrange(desc(deviance_explained))
  
  return(results_df)
}

run_unified_search <- function(data_combined, target_var, exclude_vars) {
  cat("\n", rep("=", 60), "\n")
  cat("UNIFIED MODEL SEARCH FOR", target_var, "\n")
  cat(rep("=", 60), "\n")
  
  if (length(exclude_vars) > 0) {
    cat("Excluding", length(exclude_vars), "variables\n")
  }
  
  # Prepare data
  model_data <- prepare_search_data_incomplete(data_combined, target_var, exclude_vars)
  
  if (nrow(model_data) < 20) {
    cat("Insufficient data (", nrow(model_data), " observations). Skipping.\n")
    return(NULL)
  }
  
  cat("Model data prepared:", nrow(model_data), "observations,", ncol(model_data)-1, "predictors\n")
  
  # Run searches
  results <- list()
  
  # Random Forest search
  results$rf_iterative <- rf_iterative_search_incomplete(model_data, target_var)
  
  # GAM search on complete cases
  complete_case_data <- model_data %>% na.omit()
  if (nrow(complete_case_data) >= 15) {
    cat("Running GAM search on", nrow(complete_case_data), "complete cases\n")
    results$gam_search <- gam_specification_search(complete_case_data, target_var)
  } else {
    cat("Insufficient complete cases for GAM search\n")
    results$gam_search <- NULL
  }
  
  # Print summaries
  cat("\n=== SEARCH RESULTS SUMMARY ===\n")
  
  if (!is.null(results$rf_iterative$results)) {
    best_rf <- results$rf_iterative$results %>% 
      slice_max(test_r2, n = 1, na_rm = TRUE)
    cat("Best RF model - Test R²:", round(best_rf$test_r2, 4), 
        "with", best_rf$n_vars, "predictors\n")
  }
  
  if (!is.null(results$gam_search)) {
    best_gam <- results$gam_search %>% 
      slice_max(deviance_explained, n = 1)
    cat("Best GAM model - Dev. explained:", round(best_gam$deviance_explained, 4), 
        "with", best_gam$n_vars, "predictors\n")
  }
  
  return(results)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Use the combined_data or merged_final dataset
if (!exists("combined_data")) {
  if (exists("merged_final")) {
    combined_data <- merged_final
  } else {
    stop("No dataset found. Please load your data first.")
  }
}

cat("Using dataset with", nrow(combined_data), "observations\n")

# Define exclusions with simple, explicit variable lists
# Variables to keep for ALL models: environmental, gas, wood properties, Temp_Air_125cm

# All non-probe ddPCR variables to exclude
non_probe_ddpcr <- c(
  "ddpcr_mcra_Inner_loose", "ddpcr_mcra_Inner_strict", 
  "ddpcr_mcra_Outer_loose", "ddpcr_mcra_Outer_strict",
  "ddpcr_mmox_Inner_loose", "ddpcr_mmox_Inner_strict",
  "ddpcr_mmox_Mineral_loose", "ddpcr_mmox_Mineral_strict",
  "ddpcr_mmox_Organic_loose", "ddpcr_mmox_Organic_strict", 
  "ddpcr_mmox_Outer_loose", "ddpcr_mmox_Outer_strict",
  "ddpcr_pmoa_Inner_loose", "ddpcr_pmoa_Inner_strict",
  "ddpcr_pmoa_Mineral_loose", "ddpcr_pmoa_Mineral_strict",
  "ddpcr_pmoa_Organic_loose", "ddpcr_pmoa_Organic_strict",
  "ddpcr_pmoa_Outer_loose", "ddpcr_pmoa_Outer_strict"
)

# Variables to exclude for flux models
vars_to_exclude_flux <- c(
  grep("16s", names(combined_data), value = TRUE),  # All 16S
  grep("Root", names(combined_data), value = TRUE), # All Root measurements
  grep("flux_75cm|flux_50cm|flux_200cm", names(combined_data), value = TRUE), # Non-125cm fluxes
  grep("Temp_Air_75cm|Temp_Air_50cm|Temp_Air_200cm|Temp_Air_RootCrowncm", names(combined_data), value = TRUE), # Non-125cm temps
  non_probe_ddpcr, # All non-probe ddPCR
  "plot"
)

# Variables to exclude for mcrA models - FIXED to avoid excluding target
vars_to_exclude_mcra <- c(
  grep("16s", names(combined_data), value = TRUE),  # All 16S
  grep("Root", names(combined_data), value = TRUE), # All Root measurements
  grep("flux_75cm|flux_50cm|flux_200cm", names(combined_data), value = TRUE), # Non-125cm fluxes
  grep("Temp_Air_75cm|Temp_Air_50cm|Temp_Air_200cm|Temp_Air_RootCrowncm", names(combined_data), value = TRUE), # Non-125cm temps
  non_probe_ddpcr, # All non-probe ddPCR
  grep("strict", names(combined_data), value = TRUE), # All strict variants
  # Exclude OTHER mcrA probe variables (avoid circularity) but NOT the target
  c("ddpcr_mcra_probe_Inner_loose", "ddpcr_mcra_probe_Inner_strict",
    "ddpcr_mcra_probe_Mineral_loose", "ddpcr_mcra_probe_Organic_loose", 
    "ddpcr_mcra_probe_Mineral_strict", "ddpcr_mcra_probe_Organic_strict",
    "ddpcr_mcra_probe_Outer_strict"),
  "plot"
)

cat("Variables excluded for flux models:", length(vars_to_exclude_flux), "\n")
cat("Variables excluded for mcrA models:", length(vars_to_exclude_mcra), "\n")

# Run searches
cat("\n=== SEARCHING FOR CH4 FLUX MODELS ===\n")
flux_results <- run_unified_search(combined_data, "CH4_best.flux_125cm", vars_to_exclude_flux)

cat("\n=== SEARCHING FOR mcrA MODELS ===\n")
# Use ddpcr_mcra_probe_Outer_loose as target (fewest missing values: 81)
mcra_results <- run_unified_search(combined_data, "ddpcr_mcra_probe_Outer_loose", vars_to_exclude_mcra)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n", rep("=", 80), "\n")
cat("FINAL SUMMARY OF RESULTS\n")
cat(rep("=", 80), "\n")

if (!is.null(flux_results)) {
  cat("\nCH4 FLUX MODEL RESULTS:\n")
  if (!is.null(flux_results$rf_iterative$results)) {
    best_flux_rf <- flux_results$rf_iterative$results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
    cat("- Best RF Test R²:", round(best_flux_rf$test_r2, 4), "with", best_flux_rf$n_vars, "variables\n")
    cat("- Top predictors:", flux_results$rf_iterative$importance$variable[1:5] %>% paste(collapse = ", "), "\n")
  }
  if (!is.null(flux_results$gam_search)) {
    best_flux_gam <- flux_results$gam_search %>% slice_max(deviance_explained, n = 1)
    cat("- Best GAM Dev. Explained:", round(best_flux_gam$deviance_explained, 4), "\n")
  }
}

if (!is.null(mcra_results)) {
  cat("\nmcrA MODEL RESULTS:\n")
  if (!is.null(mcra_results$rf_iterative$results)) {
    best_mcra_rf <- mcra_results$rf_iterative$results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
    cat("- Best RF Test R²:", round(best_mcra_rf$test_r2, 4), "with", best_mcra_rf$n_vars, "variables\n")
    cat("- Top predictors:", mcra_results$rf_iterative$importance$variable[1:5] %>% paste(collapse = ", "), "\n")
  }
  if (!is.null(mcra_results$gam_search)) {
    best_mcra_gam <- mcra_results$gam_search %>% slice_max(deviance_explained, n = 1)
    cat("- Best GAM Dev. Explained:", round(best_mcra_gam$deviance_explained, 4), "\n")
  }
}

cat("\n=== SEARCH COMPLETE ===\n")