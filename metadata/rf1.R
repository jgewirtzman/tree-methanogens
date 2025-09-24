# Random Forest Models Split by Soil Moisture
# Models for CH4 flux and mcrA prediction, stratified by VWC < 30% vs VWC >= 30%
# Excludes other flux measurements and 16S rRNA genes from predictors

library(randomForest)
library(dplyr)
library(ggplot2)
library(corrplot)
library(VIM)
library(caret)
library(readr)

# Load data if not already loaded
if (!exists("merged_final")) {
  merged_final <- read_csv("merged_tree_dataset_final.csv")
}

cat("=== RANDOM FOREST MODELING - SOIL MOISTURE SPLIT ===\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

cat("=== PREPARING DATA FOR MODELING ===\n")

# Create modeling dataset with relevant predictors
modeling_data <- merged_final %>%
  dplyr::select(
    # Tree characteristics
    tree_id, species_id, dbh,
    
    # Gas concentrations
    CO2_concentration, CH4_concentration, N2O_concentration, O2_concentration,
    
    # Temperature measurements
    Temp_Air_125cm,
    
    # Soil variables
    VWC_mean, ORP_mean, SoilTemp_mean, OrganicDepth_mean, MineralDepth_mean,
    
    # Wood properties
    outer_moisture_fresh_percent, inner_moisture_fresh_percent, 
    outer_moisture_dry_percent, inner_moisture_dry_percent, 
    outer_density_final, inner_density_final, 
    
    # ddPCR data (excluding 16S rRNA genes)
    starts_with("ddpcr_mcra"), starts_with("ddpcr_mmox"), starts_with("ddpcr_pmoa"),
    
    # Response variables
    CH4_best.flux_125cm,  # Target 1: Flux prediction
    ddpcr_mcra_probe_Inner_loose,  # Target 2: mcrA loose
    ddpcr_mcra_probe_Inner_strict  # Target 3: mcrA strict
  )

cat("Original dataset:", nrow(modeling_data), "trees\n")

# Check data completeness
missing_summary <- modeling_data %>%
  summarise_all(~sum(is.na(.))) %>%
  gather(variable, missing_count) %>%
  mutate(missing_percent = round(missing_count / nrow(modeling_data) * 100, 1)) %>%
  arrange(desc(missing_percent))

cat("Missing data summary (top 10):\n")
print(head(missing_summary, 10))

# =============================================================================
# SOIL MOISTURE STRATIFICATION
# =============================================================================

cat("\n=== SOIL MOISTURE STRATIFICATION ===\n")

# Check VWC distribution
cat("VWC distribution:\n")
vwc_summary <- modeling_data %>%
  filter(!is.na(VWC_mean)) %>%
  summarise(
    min = min(VWC_mean, na.rm = TRUE),
    q25 = quantile(VWC_mean, 0.25, na.rm = TRUE),
    median = median(VWC_mean, na.rm = TRUE),
    q75 = quantile(VWC_mean, 0.75, na.rm = TRUE),
    max = max(VWC_mean, na.rm = TRUE),
    n_below_30 = sum(VWC_mean < 30, na.rm = TRUE),
    n_above_30 = sum(VWC_mean >= 30, na.rm = TRUE)
  )
print(vwc_summary)

# Split modeling data by VWC
modeling_data_dry <- modeling_data %>%
  filter(!is.na(VWC_mean), VWC_mean < 30)

modeling_data_wet <- modeling_data %>%
  filter(!is.na(VWC_mean), VWC_mean >= 30)

cat("Dry soils (VWC < 30%):", nrow(modeling_data_dry), "trees\n")
cat("Wet soils (VWC >= 30%):", nrow(modeling_data_wet), "trees\n")

# =============================================================================
# FUNCTION TO RUN MODELS FOR EACH SOIL MOISTURE GROUP
# =============================================================================

run_models_by_moisture <- function(data_subset, moisture_label) {
  cat("\n===============================================\n")
  cat("MODELS FOR", moisture_label, "SOILS\n")
  cat("===============================================\n")
  
  # MODEL 1: CH4 FLUX PREDICTION
  cat("\n=== MODEL 1: CH4 FLUX PREDICTION (", moisture_label, ") ===\n")
  
  flux_data <- data_subset %>%
    filter(!is.na(CH4_best.flux_125cm)) %>%
    dplyr::select(-ddpcr_mcra_probe_Inner_loose, -ddpcr_mcra_probe_Inner_strict, -tree_id) %>%
    dplyr::select(-contains("_strict")) %>%
    mutate(species = as.factor(species_id)) %>%
    dplyr::select(-species_id) %>%
    mutate_if(is.character, as.numeric)
  
  cat("Flux modeling dataset (", moisture_label, "):", nrow(flux_data), "trees\n")
  
  if (nrow(flux_data) >= 20) {
    set.seed(123)
    train_indices <- sample(1:nrow(flux_data), 0.7 * nrow(flux_data))
    flux_train <- flux_data[train_indices, ]
    flux_test <- flux_data[-train_indices, ]
    
    cat("Training set:", nrow(flux_train), "trees\n")
    cat("Testing set:", nrow(flux_test), "trees\n")
    
    flux_rf <- randomForest(
      CH4_best.flux_125cm ~ ., 
      data = flux_train,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(flux_train) - 1))),
      importance = TRUE,
      na.action = na.roughfix
    )
    
    cat("\nFlux Model Performance (", moisture_label, "):\n")
    print(flux_rf)
    
    # Test set performance
    if (nrow(flux_test) > 0) {
      flux_pred <- predict(flux_rf, flux_test)
      complete_cases <- !is.na(flux_pred) & !is.na(flux_test$CH4_best.flux_125cm)
      
      if (sum(complete_cases) > 0) {
        flux_rmse <- sqrt(mean((flux_pred[complete_cases] - flux_test$CH4_best.flux_125cm[complete_cases])^2))
        flux_r2 <- cor(flux_pred[complete_cases], flux_test$CH4_best.flux_125cm[complete_cases])^2
        
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
    
    cat("\nTop 10 Important Variables for Flux Prediction (", moisture_label, "):\n")
    print(head(flux_imp_df, 10))
    
    # Store flux model for later comparison
    assign(paste0("flux_rf_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), flux_rf, envir = .GlobalEnv)
    assign(paste0("flux_imp_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), flux_imp_df, envir = .GlobalEnv)
  } else {
    cat("Insufficient data for flux modeling (", moisture_label, ")\n")
  }
  
  # MODEL 2: mcrA LOOSE PREDICTION
  cat("\n=== MODEL 2: mcrA ABUNDANCE (LOOSE, ", moisture_label, ") ===\n")
  
  mcra_loose_data <- data_subset %>%
    filter(!is.na(ddpcr_mcra_probe_Inner_loose)) %>%
    dplyr::select(-CH4_best.flux_125cm, -ddpcr_mcra_probe_Inner_strict, -tree_id) %>%
    dplyr::select(-ddpcr_mcra_Inner_loose, -ddpcr_mcra_Inner_strict, -ddpcr_mcra_Outer_loose, 
           -ddpcr_mcra_Outer_strict, -ddpcr_mcra_probe_Mineral_loose, -ddpcr_mcra_probe_Organic_loose,
           -ddpcr_mcra_probe_Organic_strict, -ddpcr_mcra_probe_Outer_loose, -ddpcr_mcra_probe_Outer_strict,
           -ddpcr_mcra_probe_Mineral_strict) %>%
    mutate(species = as.factor(species_id)) %>%
    dplyr::select(-species_id) %>%
    mutate_if(is.character, as.numeric)
  
  cat("mcrA loose modeling dataset (", moisture_label, "):", nrow(mcra_loose_data), "trees\n")
  
  if (nrow(mcra_loose_data) >= 20) {
    set.seed(123)
    train_indices <- sample(1:nrow(mcra_loose_data), 0.7 * nrow(mcra_loose_data))
    mcra_train <- mcra_loose_data[train_indices, ]
    mcra_test <- mcra_loose_data[-train_indices, ]
    
    mcra_loose_rf <- randomForest(
      ddpcr_mcra_probe_Inner_loose ~ ., 
      data = mcra_train,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(mcra_train) - 1))),
      importance = TRUE,
      na.action = na.roughfix
    )
    
    cat("\nmcrA Loose Model Performance (", moisture_label, "):\n")
    print(mcra_loose_rf)
    
    # Variable importance
    mcra_importance <- importance(mcra_loose_rf)
    mcra_imp_df <- data.frame(
      Variable = rownames(mcra_importance),
      IncMSE = mcra_importance[, 1],
      IncNodePurity = mcra_importance[, 2]
    ) %>%
      arrange(desc(IncMSE))
    
    cat("\nTop 10 Important Variables for mcrA Loose (", moisture_label, "):\n")
    print(head(mcra_imp_df, 10))
    
    # Store mcrA loose model
    assign(paste0("mcra_loose_rf_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_loose_rf, envir = .GlobalEnv)
    assign(paste0("mcra_loose_imp_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_imp_df, envir = .GlobalEnv)
  } else {
    cat("Insufficient data for mcrA loose modeling (", moisture_label, ")\n")
  }
  
  # MODEL 3: mcrA STRICT PREDICTION
  cat("\n=== MODEL 3: mcrA ABUNDANCE (STRICT, ", moisture_label, ") ===\n")
  
  mcra_strict_data <- data_subset %>%
    filter(!is.na(ddpcr_mcra_probe_Inner_strict)) %>%
    dplyr::select(-CH4_best.flux_125cm, -ddpcr_mcra_probe_Inner_loose, -tree_id) %>%
    dplyr::select(-ddpcr_mcra_Inner_loose, -ddpcr_mcra_Inner_strict, -ddpcr_mcra_Outer_loose, 
           -ddpcr_mcra_Outer_strict, -ddpcr_mcra_probe_Mineral_loose, -ddpcr_mcra_probe_Organic_loose,
           -ddpcr_mcra_probe_Organic_strict, -ddpcr_mcra_probe_Outer_loose, -ddpcr_mcra_probe_Outer_strict,
           -ddpcr_mcra_probe_Mineral_strict) %>%
    mutate(species = as.factor(species_id)) %>%
    dplyr::select(-species_id) %>%
    mutate_if(is.character, as.numeric)
  
  cat("mcrA strict modeling dataset (", moisture_label, "):", nrow(mcra_strict_data), "trees\n")
  
  if (nrow(mcra_strict_data) >= 20) {
    set.seed(123)
    train_indices <- sample(1:nrow(mcra_strict_data), 0.7 * nrow(mcra_strict_data))
    mcra_strict_train <- mcra_strict_data[train_indices, ]
    mcra_strict_test <- mcra_strict_data[-train_indices, ]
    
    mcra_strict_rf <- randomForest(
      ddpcr_mcra_probe_Inner_strict ~ ., 
      data = mcra_strict_train,
      ntree = 500,
      mtry = max(1, floor(sqrt(ncol(mcra_strict_train) - 1))),
      importance = TRUE,
      na.action = na.roughfix
    )
    
    cat("\nmcrA Strict Model Performance (", moisture_label, "):\n")
    print(mcra_strict_rf)
    
    # Variable importance
    mcra_strict_importance <- importance(mcra_strict_rf)
    mcra_strict_imp_df <- data.frame(
      Variable = rownames(mcra_strict_importance),
      IncMSE = mcra_strict_importance[, 1],
      IncNodePurity = mcra_strict_importance[, 2]
    ) %>%
      arrange(desc(IncMSE))
    
    cat("\nTop 10 Important Variables for mcrA Strict (", moisture_label, "):\n")
    print(head(mcra_strict_imp_df, 10))
    
    # Store mcrA strict model
    assign(paste0("mcra_strict_rf_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_strict_rf, envir = .GlobalEnv)
    assign(paste0("mcra_strict_imp_", gsub("[^a-zA-Z]", "", tolower(moisture_label))), mcra_strict_imp_df, envir = .GlobalEnv)
  } else {
    cat("Insufficient data for mcrA strict modeling (", moisture_label, ")\n")
  }
}

# Run models for each moisture group
run_models_by_moisture(modeling_data_dry, "DRY")
run_models_by_moisture(modeling_data_wet, "WET")

# =============================================================================
# MODEL COMPARISONS
# =============================================================================

cat("\n===============================================\n")
cat("MODEL COMPARISONS\n")
cat("===============================================\n")

# Compare mcrA loose vs strict predictors
cat("\n=== mcrA LOOSE vs STRICT PREDICTOR COMPARISON ===\n")

compare_mcra_predictors <- function(moisture_label) {
  loose_imp_var <- paste0("mcra_loose_imp_", gsub("[^a-zA-Z]", "", tolower(moisture_label)))
  strict_imp_var <- paste0("mcra_strict_imp_", gsub("[^a-zA-Z]", "", tolower(moisture_label)))
  
  if (exists(loose_imp_var) && exists(strict_imp_var)) {
    loose_imp <- get(loose_imp_var)
    strict_imp <- get(strict_imp_var)
    
    cat("\nTop predictors for mcrA (", moisture_label, " soils):\n")
    cat("LOOSE (top 5):\n")
    print(head(loose_imp[, c("Variable", "IncMSE")], 5))
    cat("STRICT (top 5):\n")
    print(head(strict_imp[, c("Variable", "IncMSE")], 5))
    
    # Common predictors
    loose_top10 <- loose_imp$Variable[1:min(10, nrow(loose_imp))]
    strict_top10 <- strict_imp$Variable[1:min(10, nrow(strict_imp))]
    common_predictors <- intersect(loose_top10, strict_top10)
    
    cat("Common top-10 predictors between loose and strict (", moisture_label, "):\n")
    print(common_predictors)
    cat("Number of common predictors:", length(common_predictors), "out of", 
        min(10, length(loose_top10), length(strict_top10)), "\n")
  } else {
    cat("Models not available for", moisture_label, "soils\n")
  }
}

if (nrow(modeling_data_dry) > 0) compare_mcra_predictors("DRY")
if (nrow(modeling_data_wet) > 0) compare_mcra_predictors("WET")

# Compare flux models between moisture groups
cat("\n=== FLUX MODEL COMPARISON: DRY vs WET SOILS ===\n")

if (exists("flux_imp_dry") && exists("flux_imp_wet")) {
  cat("Top 5 flux predictors in DRY soils:\n")
  print(head(flux_imp_dry[, c("Variable", "IncMSE")], 5))
  
  cat("\nTop 5 flux predictors in WET soils:\n")
  print(head(flux_imp_wet[, c("Variable", "IncMSE")], 5))
  
  # Common predictors between moisture groups
  dry_top10 <- flux_imp_dry$Variable[1:min(10, nrow(flux_imp_dry))]
  wet_top10 <- flux_imp_wet$Variable[1:min(10, nrow(flux_imp_wet))]
  common_flux_predictors <- intersect(dry_top10, wet_top10)
  
  cat("\nCommon top-10 flux predictors between DRY and WET soils:\n")
  print(common_flux_predictors)
  cat("Number of common predictors:", length(common_flux_predictors), "out of", 
      min(10, length(dry_top10), length(wet_top10)), "\n")
} else {
  cat("Flux models not available for comparison\n")
}

# Summary of ddPCR predictor importance
cat("\n=== ddPCR PREDICTOR SUMMARY ===\n")

extract_ddpcr_importance <- function(imp_df, model_name) {
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    ddpcr_vars <- imp_df %>%
      filter(grepl("ddpcr_", Variable)) %>%
      head(5)
    
    if (nrow(ddpcr_vars) > 0) {
      cat("Top ddPCR predictors for", model_name, ":\n")
      print(ddpcr_vars[, c("Variable", "IncMSE")])
    }
  }
}

if (exists("flux_imp_dry")) extract_ddpcr_importance(flux_imp_dry, "Flux (DRY)")
if (exists("flux_imp_wet")) extract_ddpcr_importance(flux_imp_wet, "Flux (WET)")
if (exists("mcra_loose_imp_dry")) extract_ddpcr_importance(mcra_loose_imp_dry, "mcrA Loose (DRY)")
if (exists("mcra_loose_imp_wet")) extract_ddpcr_importance(mcra_loose_imp_wet, "mcrA Loose (WET)")
if (exists("mcra_strict_imp_dry")) extract_ddpcr_importance(mcra_strict_imp_dry, "mcrA Strict (DRY)")
if (exists("mcra_strict_imp_wet")) extract_ddpcr_importance(mcra_strict_imp_wet, "mcrA Strict (WET)")

# =============================================================================
# VISUALIZATION WITH COLOR CODING
# =============================================================================

cat("\n=== CREATING VARIABLE IMPORTANCE PLOTS ===\n")

# Function to assign colors based on variable type
assign_variable_colors <- function(variables) {
  colors <- character(length(variables))
  
  for (i in 1:length(variables)) {
    var <- variables[i]
    
    # Inside tree or wood (green)
    if (grepl("species|dbh|density|moisture", var) || grepl("ddpcr_.*Inner", var)) {
      colors[i] <- "#2E8B57"  # Forest green
    }
    # In soil (brown)
    else if (grepl("VWC|ORP|SoilTemp|OrganicDepth|MineralDepth", var) || 
             grepl("ddpcr_.*Mineral|ddpcr_.*Organic", var)) {
      colors[i] <- "#8B4513"  # Saddle brown
    }
    # Other (blue) - gas concentrations, air temperature, outer measurements
    else {
      colors[i] <- "#4169E1"  # Royal blue
    }
  }
  
  return(colors)
}

# Function to create variable importance plots
create_importance_plot <- function(imp_df, title, n_vars = 15) {
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    plot_data <- imp_df %>%
      head(n_vars) %>%
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
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 9),
        plot.caption = element_text(size = 8, color = "gray50")
      )
    
    return(p)
  }
  return(NULL)
}

# Create plots for all models
plots_list <- list()

# Flux models
if (exists("flux_imp_dry")) {
  plots_list[["flux_dry"]] <- create_importance_plot(
    flux_imp_dry, 
    "Variable Importance: CH4 Flux Prediction (DRY Soils, VWC < 30%)"
  )
}

if (exists("flux_imp_wet")) {
  plots_list[["flux_wet"]] <- create_importance_plot(
    flux_imp_wet, 
    "Variable Importance: CH4 Flux Prediction (WET Soils, VWC ≥ 30%)"
  )
}

# mcrA loose models
if (exists("mcra_loose_imp_dry")) {
  plots_list[["mcra_loose_dry"]] <- create_importance_plot(
    mcra_loose_imp_dry, 
    "Variable Importance: mcrA Abundance (LOOSE, DRY Soils)"
  )
}

if (exists("mcra_loose_imp_wet")) {
  plots_list[["mcra_loose_wet"]] <- create_importance_plot(
    mcra_loose_imp_wet, 
    "Variable Importance: mcrA Abundance (LOOSE, WET Soils)"
  )
}

# mcrA strict models
if (exists("mcra_strict_imp_dry")) {
  plots_list[["mcra_strict_dry"]] <- create_importance_plot(
    mcra_strict_imp_dry, 
    "Variable Importance: mcrA Abundance (STRICT, DRY Soils)"
  )
}

if (exists("mcra_strict_imp_wet")) {
  plots_list[["mcra_strict_wet"]] <- create_importance_plot(
    mcra_strict_imp_wet, 
    "Variable Importance: mcrA Abundance (STRICT, WET Soils)"
  )
}

# Display all plots
for (plot_name in names(plots_list)) {
  if (!is.null(plots_list[[plot_name]])) {
    cat("Displaying plot:", plot_name, "\n")
    print(plots_list[[plot_name]])
  }
}

# Create comparison plots for flux models (side by side)
if (exists("flux_imp_dry") && exists("flux_imp_wet")) {
  cat("\n=== FLUX MODEL COMPARISON PLOT ===\n")
  
  # Prepare data for side-by-side comparison
  dry_data <- flux_imp_dry %>%
    head(10) %>%
    mutate(Soil_Type = "DRY (VWC < 30%)")
  
  wet_data <- flux_imp_wet %>%
    head(10) %>%
    mutate(Soil_Type = "WET (VWC ≥ 30%)")
  
  comparison_data <- bind_rows(dry_data, wet_data) %>%
    mutate(
      Colors = assign_variable_colors(Variable),
      Variable = factor(Variable, levels = unique(c(rev(dry_data$Variable), rev(wet_data$Variable))))
    )
  
  flux_comparison_plot <- ggplot(comparison_data, aes(x = Variable, y = IncMSE, fill = Colors)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    scale_fill_identity() +
    facet_wrap(~Soil_Type, scales = "free") +
    coord_flip() +
    labs(
      title = "CH4 Flux Predictors: DRY vs WET Soils Comparison",
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
  
  print(flux_comparison_plot)
}

# Create ddPCR-specific plots
cat("\n=== ddPCR-SPECIFIC IMPORTANCE PLOTS ===\n")

create_ddpcr_plot <- function(imp_df, title) {
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    ddpcr_data <- imp_df %>%
      filter(grepl("ddpcr_", Variable)) %>%
      head(10)
    
    if (nrow(ddpcr_data) > 0) {
      ddpcr_data <- ddpcr_data %>%
        mutate(
          Variable = factor(Variable, levels = rev(Variable)),
          Colors = assign_variable_colors(Variable)
        )
      
      p <- ggplot(ddpcr_data, aes(x = Variable, y = IncMSE)) +
        geom_col(fill = ddpcr_data$Colors, alpha = 0.8, color = "black", size = 0.3) +
        coord_flip() +
        labs(
          title = title,
          x = "ddPCR Variables",
          y = "Increase in MSE (%)",
          caption = "Green = Inner wood, Brown = Mineral/Organic soil, Blue = Outer"
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
  }
  return(NULL)
}

# ddPCR plots for flux models
if (exists("flux_imp_dry")) {
  ddpcr_flux_dry <- create_ddpcr_plot(flux_imp_dry, "ddPCR Importance: CH4 Flux (DRY Soils)")
  if (!is.null(ddpcr_flux_dry)) print(ddpcr_flux_dry)
}

if (exists("flux_imp_wet")) {
  ddpcr_flux_wet <- create_ddpcr_plot(flux_imp_wet, "ddPCR Importance: CH4 Flux (WET Soils)")
  if (!is.null(ddpcr_flux_wet)) print(ddpcr_flux_wet)
}

cat("\n=== MODELING AND VISUALIZATION COMPLETE ===\n")