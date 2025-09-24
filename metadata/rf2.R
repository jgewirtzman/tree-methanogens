# Random Forest Models with Data Transformations
# Apply log transforms, remove redundant variables, use simplified variable set
# Split by soil moisture with improved threshold

library(randomForest)
library(dplyr)
library(ggplot2)
library(corrplot)
library(GGally)
library(gridExtra)
library(VIM)
library(car)
library(moments)

# Load data if not already loaded
if (!exists("merged_final")) {
  merged_final <- read_csv("merged_tree_dataset_final.csv")
}

cat("=== TRANSFORMED RANDOM FOREST MODELING ===\n")

# =============================================================================
# DATA PREPARATION WITH TRANSFORMATIONS
# =============================================================================

cat("=== DATA PREPARATION AND TRANSFORMATIONS ===\n")

# dplyr::select simplified variable set based on EDA insights
clean_data <- merged_final %>%
  dplyr::select(
    # Tree characteristics
    tree_id, species_id, dbh,
    
    # Gas concentrations (O2 only, skip CO2 and N2O due to redundancy)
    O2_concentration,
    
    # Temperature measurements
    Temp_Air_125cm,
    
    # Soil variables
    VWC_mean, ORP_mean, SoilTemp_mean, OrganicDepth_mean, MineralDepth_mean,
    
    # Wood properties (inner only to reduce redundancy)
    inner_moisture_fresh_percent, inner_density_final,
    
    # ddPCR data - loose variants only (excluding 16S rRNA genes)
    starts_with("ddpcr_mcra"), starts_with("ddpcr_mmox"), starts_with("ddpcr_pmoa"),
    
    # Response variables
    CH4_best.flux_125cm,
    ddpcr_mcra_probe_Inner_loose,
    ddpcr_mcra_probe_Inner_strict
  ) %>%
  # Keep only loose ddPCR variants for predictors
  dplyr::select(-contains("_strict")) %>%
  # Add back the target variables
  mutate(
    ddpcr_mcra_probe_Inner_loose_target = merged_final$ddpcr_mcra_probe_Inner_loose,
    ddpcr_mcra_probe_Inner_strict_target = merged_final$ddpcr_mcra_probe_Inner_strict
  )

cat("dplyr::selected variables:", ncol(clean_data), "\n")

# Apply transformations based on EDA findings
transformed_data <- clean_data %>%
  mutate(
    # Log transform CH4 flux (handle negative values by adding offset)
    CH4_flux_offset = CH4_best.flux_125cm - min(CH4_best.flux_125cm, na.rm = TRUE) + 0.001,
    CH4_flux_log = log(CH4_flux_offset),
    
    # Log transform DBH (was right-skewed)
    dbh_log = log(dbh + 0.001),
    
    # Log transform all ddPCR variables (all were highly skewed)
    across(starts_with("ddpcr_"), ~log(.x + 0.001), .names = "{.col}_log"),
    
    # Log transform organic depth (highly skewed)
    OrganicDepth_log = log(OrganicDepth_mean + 0.001),
    
    # Keep some variables untransformed (approximately normal from EDA)
    VWC_mean = VWC_mean,
    ORP_mean = ORP_mean,
    SoilTemp_mean = SoilTemp_mean,
    Temp_Air_125cm = Temp_Air_125cm,
    MineralDepth_mean = MineralDepth_mean,
    inner_moisture_fresh_percent = inner_moisture_fresh_percent,
    inner_density_final = inner_density_final,
    O2_concentration = O2_concentration
  ) %>%
  # Keep target variables with clear names
  mutate(
    CH4_flux_target = CH4_best.flux_125cm,
    mcra_loose_target = ddpcr_mcra_probe_Inner_loose_target,
    mcra_strict_target = ddpcr_mcra_probe_Inner_strict_target
  ) %>%
  # Remove original versions of transformed variables
  dplyr::select(-CH4_best.flux_125cm, -dbh, -OrganicDepth_mean, -CH4_flux_offset,
         -ddpcr_mcra_probe_Inner_loose_target, -ddpcr_mcra_probe_Inner_strict_target) %>%
  # Keep only log-transformed ddPCR predictors
  dplyr::select(-contains("ddpcr_")) %>%
  rename_with(~str_remove(.x, "_log"), contains("ddpcr_"))

cat("After transformations:", ncol(transformed_data), "variables\n")

# Check VWC distribution for better threshold
vwc_summary <- transformed_data %>%
  filter(!is.na(VWC_mean)) %>%
  summarise(
    min = min(VWC_mean),
    q25 = quantile(VWC_mean, 0.25),
    median = median(VWC_mean),
    q75 = quantile(VWC_mean, 0.75),
    max = max(VWC_mean),
    n_below_20 = sum(VWC_mean < 20),
    n_above_20 = sum(VWC_mean >= 20)
  )

cat("VWC distribution for threshold dplyr::selection:\n")
print(vwc_summary)

# Use median (20%) as threshold for more balanced groups
vwc_threshold <- 20

# Split by soil moisture using median threshold
modeling_data_dry <- transformed_data %>%
  filter(!is.na(VWC_mean), VWC_mean < vwc_threshold)

modeling_data_wet <- transformed_data %>%
  filter(!is.na(VWC_mean), VWC_mean >= vwc_threshold)

cat("Using VWC threshold of", vwc_threshold, "%:\n")
cat("Dry soils (VWC <", vwc_threshold, "%):", nrow(modeling_data_dry), "trees\n")
cat("Wet soils (VWC >=", vwc_threshold, "%):", nrow(modeling_data_wet), "trees\n")

# =============================================================================
# ASSESS TRANSFORMATION EFFECTIVENESS
# =============================================================================

cat("\n=== ASSESSING TRANSFORMATION EFFECTIVENESS ===\n")

# Check skewness before and after transformation for key variables
assess_transform_effectiveness <- function(original, transformed, var_name) {
  orig_clean <- original[!is.na(original)]
  trans_clean <- transformed[!is.na(transformed)]
  
  if (length(orig_clean) < 10 || length(trans_clean) < 10) {
    return(data.frame(Variable = var_name, Original_Skew = NA, Transformed_Skew = NA, Improvement = NA))
  }
  
  orig_skew <- skewness(orig_clean)
  trans_skew <- skewness(trans_clean)
  improvement <- abs(orig_skew) - abs(trans_skew)
  
  return(data.frame(
    Variable = var_name,
    Original_Skew = round(orig_skew, 3),
    Transformed_Skew = round(trans_skew, 3),
    Improvement = round(improvement, 3)
  ))
}

# Test key transformations
transform_effectiveness <- bind_rows(
  assess_transform_effectiveness(clean_data$CH4_best.flux_125cm, transformed_data$CH4_flux_log, "CH4_flux"),
  assess_transform_effectiveness(clean_data$dbh, transformed_data$dbh_log, "dbh"),
  assess_transform_effectiveness(clean_data$ddpcr_mcra_probe_Inner_loose_target, 
                                 transformed_data$ddpcr_mcra_probe_Inner_loose, "mcra_loose")
)

cat("Transformation effectiveness (positive improvement = better):\n")
print(transform_effectiveness)

# =============================================================================
# FUNCTION TO RUN TRANSFORMED MODELS
# =============================================================================

run_transformed_models <- function(data_subset, moisture_label) {
  cat("\n===============================================\n")
  cat("TRANSFORMED MODELS FOR", moisture_label, "SOILS\n")
  cat("===============================================\n")
  
  # MODEL 1: CH4 FLUX PREDICTION (TRANSFORMED)
  cat("\n=== MODEL 1: CH4 FLUX PREDICTION (TRANSFORMED,", moisture_label, ") ===\n")
  
  flux_data <- data_subset %>%
    filter(!is.na(CH4_flux_target)) %>%
    dplyr::select(-mcra_loose_target, -mcra_strict_target, -tree_id, -CH4_flux_target) %>%
    mutate(
      species = as.factor(species_id),
      CH4_flux_log = log(data_subset$CH4_flux_target[!is.na(data_subset$CH4_flux_target)] + 0.001)
    ) %>%
    dplyr::select(-species_id) %>%
    mutate_if(is.character, as.numeric)
  
  cat("Flux modeling dataset (", moisture_label, "):", nrow(flux_data), "trees\n")
  cat("Number of predictors:", ncol(flux_data) - 1, "\n")
  
  if (nrow(flux_data) >= 20) {
    set.seed(123)
    train_indices <- sample(1:nrow(flux_data), 0.7 * nrow(flux_data))
    flux_train <- flux_data[train_indices, ]
    flux_test <- flux_data[-train_indices, ]
    
    cat("Training set:", nrow(flux_train), "trees\n")
    cat("Testing set:", nrow(flux_test), "trees\n")
    
    flux_rf <- randomForest(
      CH4_flux_log ~ ., 
      data = flux_train,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(flux_train) - 1))),
      importance = TRUE,
      na.action = na.roughfix
    )
    
    cat("\nTransformed Flux Model Performance (", moisture_label, "):\n")
    print(flux_rf)
    
    # Test set performance
    if (nrow(flux_test) > 0) {
      flux_pred <- predict(flux_rf, flux_test)
      complete_cases <- !is.na(flux_pred) & !is.na(flux_test$CH4_flux_log)
      
      if (sum(complete_cases) > 0) {
        flux_rmse <- sqrt(mean((flux_pred[complete_cases] - flux_test$CH4_flux_log[complete_cases])^2))
        flux_r2 <- cor(flux_pred[complete_cases], flux_test$CH4_flux_log[complete_cases])^2
        
        cat("Test Set Performance (", sum(complete_cases), "complete predictions):\n")
        cat("RMSE:", round(flux_rmse, 4), "\n")
        cat("R²:", round(flux_r2, 4), "\n")
      }
    }
    
    # Variable importance
    flux_importance <- importance(flux_rf)
    flux_imp_df <- data.frame(
      Variable = rownames(flux_importance),
      IncMSE = flux_importance[, 1],
      IncNodePurity = flux_importance[, 2]
    ) %>%
      arrange(desc(IncMSE))
    
    cat("\nTop 10 Important Variables for Transformed Flux Prediction (", moisture_label, "):\n")
    print(head(flux_imp_df, 10))
    
    # Store for later comparison
    assign(paste0("flux_rf_transformed_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), flux_rf, envir = .GlobalEnv)
    assign(paste0("flux_imp_transformed_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), flux_imp_df, envir = .GlobalEnv)
  } else {
    cat("Insufficient data for flux modeling (", moisture_label, ")\n")
  }
  
  # MODEL 2: mcrA LOOSE PREDICTION (TRANSFORMED)
  cat("\n=== MODEL 2: mcrA ABUNDANCE (LOOSE, TRANSFORMED,", moisture_label, ") ===\n")
  
  mcra_data <- data_subset %>%
    filter(!is.na(mcra_loose_target)) %>%
    dplyr::select(-CH4_flux_target, -mcra_strict_target, -tree_id, -CH4_flux_log) %>%
    # Remove any remaining mcrA predictor columns to avoid circularity
    dplyr::select(-any_of(c("ddpcr_mcra_Inner_loose", "ddpcr_mcra_Outer_loose", 
                     "ddpcr_mcra_probe_Mineral_loose", "ddpcr_mcra_probe_Organic_loose",
                     "ddpcr_mcra_probe_Outer_loose"))) %>%
    mutate(
      species = as.factor(species_id),
      mcra_loose_log = log(mcra_loose_target + 0.001)
    ) %>%
    dplyr::select(-species_id, -mcra_loose_target) %>%
    mutate_if(is.character, as.numeric)
  
  cat("mcrA loose modeling dataset (", moisture_label, "):", nrow(mcra_data), "trees\n")
  
  if (nrow(mcra_data) >= 20) {
    set.seed(123)
    train_indices <- sample(1:nrow(mcra_data), 0.7 * nrow(mcra_data))
    mcra_train <- mcra_data[train_indices, ]
    mcra_test <- mcra_data[-train_indices, ]
    
    mcra_rf <- randomForest(
      mcra_loose_log ~ ., 
      data = mcra_train,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(mcra_train) - 1))),
      importance = TRUE,
      na.action = na.roughfix
    )
    
    cat("\nTransformed mcrA Loose Model Performance (", moisture_label, "):\n")
    print(mcra_rf)
    
    # Test set performance
    if (nrow(mcra_test) > 0) {
      mcra_pred <- predict(mcra_rf, mcra_test)
      complete_cases <- !is.na(mcra_pred) & !is.na(mcra_test$mcra_loose_log)
      
      if (sum(complete_cases) > 0) {
        mcra_rmse <- sqrt(mean((mcra_pred[complete_cases] - mcra_test$mcra_loose_log[complete_cases])^2))
        mcra_r2 <- cor(mcra_pred[complete_cases], mcra_test$mcra_loose_log[complete_cases])^2
        
        cat("Test Set Performance (", sum(complete_cases), "complete predictions):\n")
        cat("RMSE:", round(mcra_rmse, 4), "\n")
        cat("R²:", round(mcra_r2, 4), "\n")
      }
    }
    
    # Variable importance
    mcra_importance <- importance(mcra_rf)
    mcra_imp_df <- data.frame(
      Variable = rownames(mcra_importance),
      IncMSE = mcra_importance[, 1],
      IncNodePurity = mcra_importance[, 2]
    ) %>%
      arrange(desc(IncMSE))
    
    cat("\nTop 10 Important Variables for Transformed mcrA Loose (", moisture_label, "):\n")
    print(head(mcra_imp_df, 10))
    
    # Store model
    assign(paste0("mcra_loose_rf_transformed_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_rf, envir = .GlobalEnv)
    assign(paste0("mcra_loose_imp_transformed_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_imp_df, envir = .GlobalEnv)
  } else {
    cat("Insufficient data for mcrA loose modeling (", moisture_label, ")\n")
  }
}

# Run models for each moisture group
run_transformed_models(modeling_data_dry, "DRY")
run_transformed_models(modeling_data_wet, "WET")

# =============================================================================
# VISUALIZATION WITH TRANSFORMATIONS
# =============================================================================

cat("\n=== CREATING TRANSFORMED MODEL PLOTS ===\n")

# Function to assign colors (same as before)
assign_variable_colors <- function(variables) {
  colors <- character(length(variables))
  
  for (i in 1:length(variables)) {
    var <- variables[i]
    
    # Inside tree or wood (green)
    if (grepl("species|dbh|density|moisture", var) || grepl("Inner", var)) {
      colors[i] <- "#2E8B57"  # Forest green
    }
    # In soil (brown)
    else if (grepl("VWC|ORP|SoilTemp|OrganicDepth|MineralDepth", var) || 
             grepl("Mineral|Organic", var)) {
      colors[i] <- "#8B4513"  # Saddle brown
    }
    # Other (blue)
    else {
      colors[i] <- "#4169E1"  # Royal blue
    }
  }
  
  return(colors)
}

# Create importance plots for transformed models
create_transformed_plot <- function(imp_df, title) {
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    plot_data <- imp_df %>%
      head(15) %>%
      mutate(
        Variable = factor(Variable, levels = rev(Variable)),
        Colors = assign_variable_colors(Variable)
      )
    
    p <- ggplot(plot_data, aes(x = Variable, y = IncMSE)) +
      geom_col(fill = plot_data$Colors, alpha = 0.8, color = "black", size = 0.3) +
      coord_flip() +
      labs(
        title = title,
        x = "Variables",
        y = "Increase in MSE (%)",
        caption = "Green = Tree/Wood, Brown = Soil, Blue = Other"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9),
        plot.caption = element_text(size = 8, color = "gray50")
      )
    
    return(p)
  }
  return(NULL)
}

# Create plots for transformed models
if (exists("flux_imp_transformed_dry")) {
  p1 <- create_transformed_plot(flux_imp_transformed_dry, "Transformed CH4 Flux Predictors (DRY Soils)")
  if (!is.null(p1)) print(p1)
}

if (exists("flux_imp_transformed_wet")) {
  p2 <- create_transformed_plot(flux_imp_transformed_wet, "Transformed CH4 Flux Predictors (WET Soils)")
  if (!is.null(p2)) print(p2)
}

if (exists("mcra_loose_imp_transformed_dry")) {
  p3 <- create_transformed_plot(mcra_loose_imp_transformed_dry, "Transformed mcrA Loose Predictors (DRY Soils)")
  if (!is.null(p3)) print(p3)
}

if (exists("mcra_loose_imp_transformed_wet")) {
  p4 <- create_transformed_plot(mcra_loose_imp_transformed_wet, "Transformed mcrA Loose Predictors (WET Soils)")
  if (!is.null(p4)) print(p4)
}

# Comparison plot for transformed flux models
if (exists("flux_imp_transformed_dry") && exists("flux_imp_transformed_wet")) {
  cat("\n=== TRANSFORMED FLUX MODEL COMPARISON ===\n")
  
  dry_data <- flux_imp_transformed_dry %>%
    head(10) %>%
    mutate(Soil_Type = "DRY (< 20% VWC)")
  
  wet_data <- flux_imp_transformed_wet %>%
    head(10) %>%
    mutate(Soil_Type = "WET (≥ 20% VWC)")
  
  comparison_data <- bind_rows(dry_data, wet_data) %>%
    mutate(
      Colors = assign_variable_colors(Variable),
      Variable = factor(Variable, levels = unique(c(rev(dry_data$Variable), rev(wet_data$Variable))))
    )
  
  comparison_plot <- ggplot(comparison_data, aes(x = Variable, y = IncMSE, fill = Colors)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    scale_fill_identity() +
    facet_wrap(~Soil_Type, scales = "free") +
    coord_flip() +
    labs(
      title = "Transformed CH4 Flux Predictors: DRY vs WET Soils",
      x = "Variables",
      y = "Increase in MSE (%)",
      caption = "Green = Tree/Wood, Brown = Soil, Blue = Other"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 8, color = "gray50")
    )
  
  print(comparison_plot)
}

cat("\n=== TRANSFORMED MODELING COMPLETE ===\n")

# Summary of improvements
cat("\n=== SUMMARY OF IMPROVEMENTS ===\n")
cat("1. Applied log transformations to highly skewed variables\n")
cat("2. Removed redundant variables (strict ddPCR, CO2, N2O, outer moisture)\n")
cat("3. Used median VWC threshold (20%) for better balanced groups\n")
cat("4. Simplified predictor set to reduce multicollinearity\n")
cat("5. Models now predict log-transformed response variables\n")