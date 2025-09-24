# Modified Model Searches
# Version 1: CH4 flux without CO2 flux as predictor
# Version 2: CH4 flux requiring mcrA as predictor

cat("=== RUNNING MODIFIED SEARCHES ===\n")

# =============================================================================
# VERSION 1: CH4 FLUX WITHOUT CO2 FLUX
# =============================================================================

cat("\n=== VERSION 1: CH4 FLUX WITHOUT CO2 FLUX ===\n")

# Modified exclusions for flux models - add CO2 flux to exclusions
vars_to_exclude_flux_no_co2 <- c(
  grep("16s", names(combined_data), value = TRUE),  # All 16S
  grep("Root", names(combined_data), value = TRUE), # All Root measurements
  grep("flux_75cm|flux_50cm|flux_200cm", names(combined_data), value = TRUE), # Non-125cm fluxes
  grep("Temp_Air_75cm|Temp_Air_50cm|Temp_Air_200cm|Temp_Air_RootCrowncm", names(combined_data), value = TRUE), # Non-125cm temps
  non_probe_ddpcr, # All non-probe ddPCR
  "CO2_best.flux_125cm", # EXCLUDE CO2 flux
  "plot"
)

cat("Variables excluded for flux models (no CO2):", length(vars_to_exclude_flux_no_co2), "\n")

flux_results_no_co2 <- run_unified_search(combined_data, "CH4_best.flux_125cm", vars_to_exclude_flux_no_co2)

# =============================================================================
# VERSION 2: CH4 FLUX REQUIRING mcrA
# =============================================================================

cat("\n=== VERSION 2: CH4 FLUX REQUIRING mcrA ===\n")

# Modified function that forces mcrA variables to be included
run_unified_search_force_mcra <- function(data_combined, target_var, exclude_vars, force_vars = NULL) {
  cat("\n", rep("=", 60), "\n")
  cat("UNIFIED MODEL SEARCH FOR", target_var, "(FORCING mcrA INCLUSION)\n")
  cat(rep("=", 60), "\n")
  
  if (length(exclude_vars) > 0) {
    cat("Excluding", length(exclude_vars), "variables\n")
  }
  if (length(force_vars) > 0) {
    cat("Forcing inclusion of:", paste(force_vars, collapse = ", "), "\n")
  }
  
  # Prepare data
  model_data <- prepare_search_data_incomplete(data_combined, target_var, exclude_vars)
  
  # Ensure forced variables are present
  missing_force_vars <- setdiff(force_vars, names(model_data))
  if (length(missing_force_vars) > 0) {
    cat("Warning: Some forced variables not available:", paste(missing_force_vars, collapse = ", "), "\n")
  }
  
  available_force_vars <- intersect(force_vars, names(model_data))
  if (length(available_force_vars) == 0) {
    cat("No forced variables available. Running standard search.\n")
    return(run_unified_search(data_combined, target_var, exclude_vars))
  }
  
  if (nrow(model_data) < 20) {
    cat("Insufficient data (", nrow(model_data), " observations). Skipping.\n")
    return(NULL)
  }
  
  cat("Model data prepared:", nrow(model_data), "observations,", ncol(model_data)-1, "predictors\n")
  cat("Available forced variables:", paste(available_force_vars, collapse = ", "), "\n")
  
  # Run searches
  results <- list()
  
  # Modified RF search that always includes forced variables
  cat("=== RANDOM FOREST SEARCH (FORCING mcrA) ===\n")
  
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
  
  # Test models with forced variables + top k additional variables
  test_rf_subset_forced <- function(k) {
    # Always include forced variables
    other_vars <- setdiff(importance_df$variable, available_force_vars)
    n_additional <- min(k, length(other_vars))
    
    if (n_additional > 0) {
      additional_vars <- other_vars[1:n_additional]
      all_vars <- c(available_force_vars, additional_vars)
    } else {
      all_vars <- available_force_vars
    }
    
    tryCatch({
      set.seed(123)
      train_idx <- sample(1:nrow(model_data), 0.7 * nrow(model_data))
      train_data <- model_data[train_idx, ]
      test_data <- model_data[-train_idx, ]
      
      rf_model <- randomForest(
        as.formula(paste(target_var, "~", paste(all_vars, collapse = " + "))), 
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
        n_vars = length(all_vars),
        n_additional = n_additional,
        variables = paste(all_vars, collapse = ", "),
        forced_vars = paste(available_force_vars, collapse = ", "),
        oob_mse = tail(rf_model$mse, 1),
        oob_rsq = tail(rf_model$rsq, 1),
        test_rmse = test_rmse,
        test_r2 = test_r2
      )
    }, error = function(e) {
      tibble(
        n_vars = length(available_force_vars) + k,
        n_additional = k,
        variables = paste(c(available_force_vars, rep("ERROR", k)), collapse = ", "),
        forced_vars = paste(available_force_vars, collapse = ", "),
        oob_mse = NA,
        oob_rsq = NA,
        test_rmse = NA,
        test_r2 = NA
      )
    })
  }
  
  # Test different numbers of additional variables (beyond forced ones)
  other_predictors <- setdiff(predictors, available_force_vars)
  k_values <- c(0:min(15, length(other_predictors)))
  cat("Testing", length(k_values), "combinations (forced + 0 to", max(k_values), "additional vars)...\n")
  
  results_forced <- map_dfr(k_values, test_rf_subset_forced)
  cat("RF search with forced variables completed!\n")
  
  results$rf_forced <- list(
    results = results_forced,
    importance = importance_df,
    forced_vars = available_force_vars
  )
  
  # GAM search with forced variables
  complete_case_data <- model_data %>% na.omit()
  if (nrow(complete_case_data) >= 15) {
    cat("Running GAM search with forced variables on", nrow(complete_case_data), "complete cases\n")
    
    # Modified GAM that always includes forced variables
    gam_results_list <- list()
    other_predictors_complete <- setdiff(names(complete_case_data), c(target_var, available_force_vars))
    
    if (length(other_predictors_complete) > 0) {
      # Use RF to rank the additional variables
      rf_temp <- randomForest(as.formula(paste(target_var, "~ .")), 
                              data = complete_case_data, importance = TRUE)
      imp_order <- rownames(importance(rf_temp))[order(importance(rf_temp)[,1], decreasing = TRUE)]
      other_vars_ranked <- intersect(imp_order, other_predictors_complete)
      
      # Test forced + different numbers of additional variables
      for(i in 0:min(8, length(other_vars_ranked))) {
        if (i == 0) {
          vars <- available_force_vars
        } else {
          vars <- c(available_force_vars, other_vars_ranked[1:i])
        }
        
        tryCatch({
          formula_str <- paste(target_var, "~", paste(vars, collapse = " + "))
          gam_model <- gam(as.formula(formula_str), data = complete_case_data)
          
          gam_results_list[[length(gam_results_list) + 1]] <- tibble(
            model_type = "forced_gam",
            n_vars = length(vars),
            n_forced = length(available_force_vars),
            variables = paste(vars, collapse = ", "),
            forced_vars = paste(available_force_vars, collapse = ", "),
            aic = AIC(gam_model),
            deviance_explained = summary(gam_model)$dev.expl,
            r_squared = summary(gam_model)$r.sq
          )
        }, error = function(e) NULL)
      }
    } else {
      # Only forced variables available
      tryCatch({
        formula_str <- paste(target_var, "~", paste(available_force_vars, collapse = " + "))
        gam_model <- gam(as.formula(formula_str), data = complete_case_data)
        
        gam_results_list[[1]] <- tibble(
          model_type = "forced_gam",
          n_vars = length(available_force_vars),
          n_forced = length(available_force_vars),
          variables = paste(available_force_vars, collapse = ", "),
          forced_vars = paste(available_force_vars, collapse = ", "),
          aic = AIC(gam_model),
          deviance_explained = summary(gam_model)$dev.expl,
          r_squared = summary(gam_model)$r.sq
        )
      }, error = function(e) NULL)
    }
    
    if (length(gam_results_list) > 0) {
      results$gam_forced <- bind_rows(gam_results_list) %>%
        arrange(desc(deviance_explained))
    } else {
      results$gam_forced <- NULL
    }
    
    cat("GAM search with forced variables completed!\n")
  } else {
    cat("Insufficient complete cases for GAM search\n")
    results$gam_forced <- NULL
  }
  
  # Print summaries
  cat("\n=== FORCED SEARCH RESULTS SUMMARY ===\n")
  
  if (!is.null(results$rf_forced$results)) {
    best_rf <- results$rf_forced$results %>% 
      slice_max(test_r2, n = 1, na_rm = TRUE)
    cat("Best RF model - Test R²:", round(best_rf$test_r2, 4), 
        "with", best_rf$n_vars, "predictors (", best_rf$n_additional, "additional +", 
        length(available_force_vars), "forced)\n")
    cat("Forced variables:", best_rf$forced_vars, "\n")
  }
  
  if (!is.null(results$gam_forced)) {
    best_gam <- results$gam_forced %>% 
      slice_max(deviance_explained, n = 1)
    cat("Best GAM model - Dev. explained:", round(best_gam$deviance_explained, 4), 
        "with", best_gam$n_vars, "predictors (", best_gam$n_forced, "forced)\n")
  }
  
  return(results)
}

# Get available mcrA probe variables
mcra_probe_vars <- grep("mcra_probe.*loose", names(combined_data), value = TRUE)
mcra_probe_vars <- setdiff(mcra_probe_vars, grep("strict", mcra_probe_vars, value = TRUE))

cat("Available mcrA probe variables to force:", paste(mcra_probe_vars, collapse = ", "), "\n")

flux_results_force_mcra <- run_unified_search_force_mcra(
  combined_data, 
  "CH4_best.flux_125cm", 
  vars_to_exclude_flux,
  force_vars = mcra_probe_vars
)

# =============================================================================
# COMPARISON SUMMARY
# =============================================================================

cat("\n", rep("=", 80), "\n")
cat("COMPARISON OF ALL THREE APPROACHES\n")
cat(rep("=", 80), "\n")

# Original results
if (!is.null(flux_results)) {
  best_original_rf <- flux_results$rf_iterative$results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
  cat("\n1. ORIGINAL (with CO2 flux allowed):\n")
  cat("   - Best RF Test R²:", round(best_original_rf$test_r2, 4), "with", best_original_rf$n_vars, "variables\n")
  cat("   - Top predictors:", flux_results$rf_iterative$importance$variable[1:5] %>% paste(collapse = ", "), "\n")
}

# No CO2 results
if (!is.null(flux_results_no_co2)) {
  best_no_co2_rf <- flux_results_no_co2$rf_iterative$results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
  cat("\n2. WITHOUT CO2 FLUX:\n")
  cat("   - Best RF Test R²:", round(best_no_co2_rf$test_r2, 4), "with", best_no_co2_rf$n_vars, "variables\n")
  cat("   - Top predictors:", flux_results_no_co2$rf_iterative$importance$variable[1:5] %>% paste(collapse = ", "), "\n")
}

# Forced mcrA results
if (!is.null(flux_results_force_mcra)) {
  best_forced_rf <- flux_results_force_mcra$rf_forced$results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
  cat("\n3. FORCING mcrA INCLUSION:\n")
  cat("   - Best RF Test R²:", round(best_forced_rf$test_r2, 4), "with", best_forced_rf$n_vars, "variables\n")
  cat("   - Forced variables:", best_forced_rf$forced_vars, "\n")
  cat("   - Top overall predictors:", flux_results_force_mcra$rf_forced$importance$variable[1:5] %>% paste(collapse = ", "), "\n")
}

cat("\n=== COMPARATIVE ANALYSIS COMPLETE ===\n")




# Modified function that requires at least one mcrA variable (not all)
run_unified_search_require_mcra <- function(data_combined, target_var, exclude_vars, mcra_vars) {
  cat("\n", rep("=", 60), "\n")
  cat("UNIFIED MODEL SEARCH FOR", target_var, "(REQUIRING ≥1 mcrA)\n")
  cat(rep("=", 60), "\n")
  
  # Prepare data
  model_data <- prepare_search_data_incomplete(data_combined, target_var, exclude_vars)
  
  # Check which mcrA variables are available
  available_mcra <- intersect(mcra_vars, names(model_data))
  cat("Available mcrA variables:", paste(available_mcra, collapse = ", "), "\n")
  
  if (length(available_mcra) == 0) {
    cat("No mcrA variables available. Running standard search.\n")
    return(run_unified_search(data_combined, target_var, exclude_vars))
  }
  
  cat("Model data prepared:", nrow(model_data), "observations,", ncol(model_data)-1, "predictors\n")
  
  # Get importance ranking from full model
  rf_full <- randomForest(as.formula(paste(target_var, "~ .")), 
                          data = model_data, ntree = 500, importance = TRUE, na.action = na.roughfix)
  
  importance_df <- importance(rf_full) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    arrange(desc(`%IncMSE`))
  
  # Test models ensuring at least one mcrA variable is always included
  test_rf_with_mcra <- function(k_total) {
    other_vars <- setdiff(importance_df$variable, available_mcra)
    
    # For each combination size, test different ways to include mcrA
    results_list <- list()
    
    # Try each mcrA variable as the "required" one, then add k_total-1 other variables
    for (required_mcra in available_mcra) {
      remaining_slots <- k_total - 1
      if (remaining_slots <= 0) {
        # Only the required mcrA variable
        selected_vars <- required_mcra
      } else {
        # Required mcrA + top remaining_slots other variables
        other_candidates <- c(setdiff(available_mcra, required_mcra), other_vars)
        n_others <- min(remaining_slots, length(other_candidates))
        if (n_others > 0) {
          selected_vars <- c(required_mcra, other_candidates[1:n_others])
        } else {
          selected_vars <- required_mcra
        }
      }
      
      tryCatch({
        set.seed(123)
        train_idx <- sample(1:nrow(model_data), 0.7 * nrow(model_data))
        train_data <- model_data[train_idx, ]
        test_data <- model_data[-train_idx, ]
        
        rf_model <- randomForest(
          as.formula(paste(target_var, "~", paste(selected_vars, collapse = " + "))), 
          data = train_data, ntree = 500, na.action = na.roughfix)
        
        if (nrow(test_data) > 0) {
          predictions <- predict(rf_model, test_data)
          actual <- test_data[[target_var]]
          complete_cases <- !is.na(predictions) & !is.na(actual)
          
          if (sum(complete_cases) > 0) {
            test_r2 <- cor(predictions[complete_cases], actual[complete_cases])^2
          } else {
            test_r2 <- NA
          }
        } else {
          test_r2 <- NA
        }
        
        results_list[[length(results_list) + 1]] <- tibble(
          n_vars = length(selected_vars),
          required_mcra = required_mcra,
          variables = paste(selected_vars, collapse = ", "),
          test_r2 = test_r2,
          oob_rsq = tail(rf_model$rsq, 1)
        )
      }, error = function(e) {
        results_list[[length(results_list) + 1]] <- tibble(
          n_vars = k_total,
          required_mcra = required_mcra,
          variables = "ERROR",
          test_r2 = NA,
          oob_rsq = NA
        )
      })
    }
    
    # Return the best result for this k_total
    all_results <- bind_rows(results_list)
    if (nrow(all_results) > 0 && any(!is.na(all_results$test_r2))) {
      best_result <- all_results %>% slice_max(test_r2, n = 1, na_rm = TRUE)
      return(best_result)
    } else {
      return(tibble(n_vars = k_total, required_mcra = NA, variables = NA, test_r2 = NA, oob_rsq = NA))
    }
  }
  
  cat("Testing models with 1 to 15 variables (each requiring ≥1 mcrA)...\n")
  results_df <- map_dfr(1:15, test_rf_with_mcra)
  
  cat("Search requiring mcrA completed!\n")
  
  # Print summary
  best_result <- results_df %>% slice_max(test_r2, n = 1, na_rm = TRUE)
  cat("Best model - Test R²:", round(best_result$test_r2, 4), 
      "with", best_result$n_vars, "variables\n")
  cat("Required mcrA variable:", best_result$required_mcra, "\n")
  
  return(list(
    results = results_df,
    importance = importance_df,
    best_mcra = best_result$required_mcra
  ))
}

# Run the corrected version
flux_results_require_mcra <- run_unified_search_require_mcra(
  combined_data, 
  "CH4_best.flux_125cm", 
  vars_to_exclude_flux,
  mcra_probe_vars
)