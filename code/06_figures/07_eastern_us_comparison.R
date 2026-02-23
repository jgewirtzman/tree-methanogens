# ==============================================================================
# Eastern US Forest CH4 Comparison
# ==============================================================================
# Purpose: Compares Yale Myers CH4 data to broader Eastern US forest flux
#   measurements.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - easternforest_clean_wmeta.csv (from data/raw/external/)
# ==============================================================================

# Load required libraries
library(dplyr)
library(broom)
library(ggplot2)
library(car)
library(lme4)
library(lmerTest)
library(gridExtra)

# Read the data

# Check the structure
cat("Data Structure:\n")
str(data)

# Convert categorical variables to factors
data <- data %>%
  mutate(
    gen = as.factor(gen),
    season = as.factor(season),
    bot = as.factor(bot),
    age = as.factor(age)
  )

# Summary of the data
cat("\nData Summary:\n")
summary(data)

# Check for required columns
required_cols <- c("ch4", "co2", "o2", "n2o", "gen", "season", "lat", "bot", "age", "site", "tree")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ============================================================================
# LINEAR MODELS FOR EACH GAS
# ============================================================================

# List of gases to model
gases <- c("ch4", "co2", "o2", "n2o")

# Store all models and results
models <- list()
model_summaries <- list()
anova_results <- list()

for (gas in gases) {
  
  cat("\n", rep("=", 80), "\n")
  cat("MODELING:", toupper(gas), "\n")
  cat(rep("=", 80), "\n\n")
  
  # Build formula
  formula <- as.formula(paste(gas, "~ gen + season + lat + bot + age"))
  
  # Fit linear model
  model <- lm(formula, data = data)
  
  # Store model
  models[[gas]] <- model
  
  # Print summary
  print(summary(model))
  
  # ANOVA Type II (for unbalanced designs)
  cat("\n--- ANOVA Type II Results ---\n")
  anova_ii <- Anova(model, type = "II")
  print(anova_ii)
  anova_results[[gas]] <- anova_ii
  
  # Store tidy summary
  model_summaries[[gas]] <- tidy(model)
  
  # Model diagnostics - DISPLAY
  cat("\n--- Model Diagnostics ---\n")
  par(mfrow = c(2, 2))
  plot(model, main = paste("Diagnostics for", toupper(gas)))
  par(mfrow = c(1, 1))  # Reset plotting layout
  
  # VIF for multicollinearity (exclude intercept)
  cat("\nVariance Inflation Factors (VIF):\n")
  vif_values <- tryCatch({
    vif(model)
  }, error = function(e) {
    cat("Could not calculate VIF:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(vif_values)) {
    print(vif_values)
  }
  
  cat("\n")
}

# ============================================================================
# MIXED EFFECTS MODELS (accounting for site/tree random effects)
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("MIXED EFFECTS MODELS (with random effects for site and tree)\n")
cat(rep("=", 80), "\n\n")

mixed_models <- list()

for (gas in gases) {
  
  cat("\n--- MIXED MODEL FOR:", toupper(gas), "---\n")
  
  # Build formula with random effects
  formula_mixed <- as.formula(paste(gas, "~ gen + season + lat + bot + age + (1|site) + (1|tree)"))
  
  # Fit mixed effects model
  mixed_model <- tryCatch({
    lmer(formula_mixed, data = data, REML = TRUE)
  }, error = function(e) {
    cat("Error fitting mixed model:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(mixed_model)) {
    # Store model
    mixed_models[[gas]] <- mixed_model
    
    # Print summary
    print(summary(mixed_model))
    
    # ANOVA for fixed effects
    cat("\nANOVA for fixed effects:\n")
    print(anova(mixed_model))
  }
  
  cat("\n")
}

# ============================================================================
# COMPILE RESULTS INTO TABLES
# ============================================================================

# Combine all model summaries
all_coefficients <- bind_rows(model_summaries, .id = "gas")

cat("\n--- All Model Coefficients ---\n")
print(all_coefficients)

# Combine ANOVA results
anova_combined <- bind_rows(
  lapply(names(anova_results), function(gas_name) {
    as.data.frame(anova_results[[gas_name]]) %>%
      mutate(gas = gas_name, predictor = rownames(anova_results[[gas_name]]))
  })
)

cat("\n--- Combined ANOVA Results ---\n")
print(anova_combined)

# ============================================================================
# MODEL COMPARISON TABLE
# ============================================================================

model_comparison <- data.frame(
  gas = gases,
  r_squared = sapply(models, function(m) summary(m)$r.squared),
  adj_r_squared = sapply(models, function(m) summary(m)$adj.r.squared),
  rmse = sapply(models, function(m) sqrt(mean(residuals(m)^2))),
  aic = sapply(models, AIC),
  bic = sapply(models, BIC)
)

cat("\n--- Model Comparison ---\n")
print(model_comparison)

# ============================================================================
# EFFECT PLOTS
# ============================================================================

# Function to create and DISPLAY effect plots
create_effect_plots <- function(model, gas_name, data) {
  
  # Get complete cases used in the model
  model_data <- model$model
  
  # Predicted vs Observed (using only data from the model)
  data_pred <- model_data %>%
    mutate(predicted = fitted(model),
           resid = residuals(model))
  
  p1 <- ggplot(data_pred, aes(x = predicted, y = .data[[gas_name]])) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_bw() +
    labs(title = paste("Predicted vs Observed -", toupper(gas_name)),
         x = "Predicted", y = "Observed")
  
  # Effects by genus (using full data for visualization)
  p2 <- ggplot(data, aes(x = gen, y = .data[[gas_name]], fill = gen)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(toupper(gas_name), "by Genus"),
         x = "Genus", y = gas_name) +
    guides(fill = "none")
  
  # Effects by season
  p3 <- ggplot(data, aes(x = season, y = .data[[gas_name]], fill = season)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = paste(toupper(gas_name), "by Season"),
         x = "Season", y = gas_name) +
    guides(fill = "none")
  
  # Effects by latitude
  p4 <- ggplot(data, aes(x = lat, y = .data[[gas_name]])) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    theme_bw() +
    labs(title = paste(toupper(gas_name), "by Latitude"),
         x = "Latitude", y = gas_name)
  
  # Effects by bot (bottom/top of tree)
  p5 <- ggplot(data, aes(x = bot, y = .data[[gas_name]], fill = bot)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = paste(toupper(gas_name), "by Position (Bot)"),
         x = "Position", y = gas_name) +
    guides(fill = "none")
  
  # Effects by age
  p6 <- ggplot(data, aes(x = age, y = .data[[gas_name]], fill = age)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(title = paste(toupper(gas_name), "by Age"),
         x = "Age", y = gas_name) +
    guides(fill = "none")
  
  # Combine and DISPLAY plots
  grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
}

# Create and display effect plots for each gas
cat("\nDisplaying effect plots...\n")
for (gas in gases) {
  cat("\n", rep("-", 80), "\n")
  cat("PLOTS FOR:", toupper(gas), "\n")
  cat(rep("-", 80), "\n")
  tryCatch({
    create_effect_plots(models[[gas]], gas, data)
  }, error = function(e) {
    cat("Error creating plots for", gas, ":", e$message, "\n")
  })
}

# ============================================================================
# SUMMARY STATISTICS BY PREDICTORS
# ============================================================================

summary_by_genus <- data %>%
  group_by(gen) %>%
  summarise(across(all_of(gases), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        n = ~sum(!is.na(.))),
                   .names = "{.col}_{.fn}"))

cat("\n--- Summary by Genus ---\n")
print(summary_by_genus)

summary_by_season <- data %>%
  group_by(season) %>%
  summarise(across(all_of(gases), 
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        n = ~sum(!is.na(.))),
                   .names = "{.col}_{.fn}"))

cat("\n--- Summary by Season ---\n")
print(summary_by_season)

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n")


















data<-read.csv('../../data/raw/external/easternforest_clean_wmeta.csv')

# Load required libraries
library(dplyr)
library(ggplot2)
library(randomForest)
library(gridExtra)


# Convert categorical variables to factors
data <- data %>%
  mutate(
    gen = as.factor(gen),
    season = as.factor(season),
    bot = as.factor(bot),
    age = as.factor(age),
    site = as.factor(site)
  )

cat("Data prepared for Random Forest analysis\n")
cat("Original dataset:", nrow(data), "observations\n\n")

# ============================================================================
# RANDOM FOREST MODELS FOR EACH GAS
# ============================================================================

gases <- c("ch4", "co2", "o2", "n2o")
rf_models <- list()
rf_results <- list()

set.seed(123)

for (gas in gases) {
  
  cat("\n", rep("=", 80), "\n")
  cat("RANDOM FOREST MODEL FOR:", toupper(gas), "\n")
  cat(rep("=", 80), "\n\n")
  
  # Prepare data
  predictor_vars <- c("gen", "season", "lat", "bot", "age", "dbh")
  all_vars <- c(gas, predictor_vars)
  
  model_data <- data %>%
    select(all_of(all_vars)) %>%
    na.omit()
  
  cat("Sample size after removing NAs:", nrow(model_data), "\n")
  cat("Predictors:", paste(predictor_vars, collapse = ", "), "\n\n")
  
  # Fit Random Forest
  cat("Fitting Random Forest (500 trees)...\n")
  
  rf_model <- randomForest(
    x = model_data[, predictor_vars],
    y = model_data[[gas]],
    ntree = 500,
    importance = TRUE,
    mtry = 3,
    keep.forest = TRUE
  )
  
  rf_models[[gas]] <- rf_model
  
  # Calculate performance metrics
  predictions <- predict(rf_model)
  actual <- model_data[[gas]]
  
  # R-squared
  ss_res <- sum((actual - predictions)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  
  # Variance explained (OOB)
  var_explained <- rf_model$rsq[500]
  
  cat("\n--- Model Performance ---\n")
  cat("Out-of-Bag Variance Explained:", round(var_explained * 100, 2), "%\n")
  cat("Out-of-Bag MSE:", formatC(rf_model$mse[500], format = "e", digits = 2), "\n")
  cat("Out-of-Bag RMSE:", formatC(sqrt(rf_model$mse[500]), format = "e", digits = 2), "\n")
  
  # Variable importance
  cat("\n--- Variable Importance (Top predictors) ---\n")
  importance_df <- as.data.frame(importance(rf_model))
  importance_df$Variable <- rownames(importance_df)
  importance_df <- importance_df %>%
    arrange(desc(`%IncMSE`))
  
  print(importance_df[, c("Variable", "%IncMSE", "IncNodePurity")])
  
  # Store results
  rf_results[[gas]] <- list(
    model = rf_model,
    var_explained = var_explained,
    r_squared = r_squared,
    rmse = sqrt(rf_model$mse[500]),
    mse = rf_model$mse[500],
    importance = importance_df,
    n_obs = nrow(model_data)
  )
  
  cat("\n")
}

# ============================================================================
# COMPARISON WITH LINEAR MODELS
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("COMPARING RANDOM FOREST vs LINEAR REGRESSION\n")
cat(rep("=", 80), "\n\n")

# Fit linear models for comparison
lm_results <- list()

for (gas in gases) {
  
  predictor_vars <- c("gen", "season", "lat", "bot", "age", "dbh")
  all_vars <- c(gas, predictor_vars)
  
  model_data <- data %>%
    select(all_of(all_vars)) %>%
    na.omit()
  
  formula <- as.formula(paste(gas, "~", paste(predictor_vars, collapse = " + ")))
  lm_model <- lm(formula, data = model_data)
  
  lm_results[[gas]] <- list(
    r_squared = summary(lm_model)$r.squared,
    adj_r_squared = summary(lm_model)$adj.r.squared,
    rmse = sqrt(mean(residuals(lm_model)^2))
  )
}

# Create comparison table
comparison_df <- data.frame(
  Gas = gases,
  LM_R2 = sapply(gases, function(g) round(lm_results[[g]]$r_squared, 4)),
  LM_RMSE = sapply(gases, function(g) formatC(lm_results[[g]]$rmse, format = "e", digits = 2)),
  RF_VarExp = sapply(gases, function(g) round(rf_results[[g]]$var_explained, 4)),
  RF_RMSE = sapply(gases, function(g) formatC(rf_results[[g]]$rmse, format = "e", digits = 2)),
  Improvement = sapply(gases, function(g) {
    round((rf_results[[g]]$var_explained - lm_results[[g]]$r_squared) * 100, 1)
  }),
  N = sapply(gases, function(g) rf_results[[g]]$n_obs)
)

cat("Model Performance Comparison:\n")
print(comparison_df)
cat("\nNote: Improvement = (RF_VarExp - LM_R2) * 100 (percentage points)\n")

# ============================================================================
# VARIABLE IMPORTANCE PLOTS
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("CREATING VARIABLE IMPORTANCE PLOTS\n")
cat(rep("=", 80), "\n\n")

# Create combined importance plot
all_importance <- bind_rows(
  lapply(gases, function(g) {
    df <- rf_results[[g]]$importance
    df$gas <- g
    df
  })
)

# Plot for each gas
for (gas in gases) {
  
  importance_df <- rf_results[[gas]]$importance %>%
    head(10)  # Top 10 variables
  
  p <- ggplot(importance_df, aes(x = reorder(Variable, `%IncMSE`), y = `%IncMSE`)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste("Variable Importance for", toupper(gas)),
      subtitle = "Random Forest - % Increase in MSE",
      x = NULL,
      y = "% Increase in MSE when variable permuted"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 10)
    )
  
  print(p)
}

# ============================================================================
# PARTIAL DEPENDENCE PLOTS FOR TOP VARIABLES
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("CREATING PARTIAL DEPENDENCE PLOTS\n")
cat(rep("=", 80), "\n\n")

for (gas in gases) {
  
  cat("Creating plots for", toupper(gas), "...\n")
  
  # Get data
  predictor_vars <- c("gen", "season", "lat", "bot", "age", "dbh")
  all_vars <- c(gas, predictor_vars)
  model_data <- data %>%
    select(all_of(all_vars)) %>%
    na.omit()
  
  # Get top 3 numeric/meaningful variables
  top_vars <- rf_results[[gas]]$importance %>%
    filter(Variable %in% c("lat", "dbh", "gen", "season", "age", "bot")) %>%
    head(3) %>%
    pull(Variable)
  
  plot_list <- list()
  
  for (var in top_vars) {
    
    if (var %in% c("lat", "dbh")) {
      # Numeric variable
      var_range <- seq(
        quantile(model_data[[var]], 0.05, na.rm = TRUE),
        quantile(model_data[[var]], 0.95, na.rm = TRUE),
        length.out = 50
      )
      
      pd_values <- sapply(var_range, function(val) {
        temp_data <- model_data
        temp_data[[var]] <- val
        mean(predict(rf_models[[gas]], temp_data))
      })
      
      pd_df <- data.frame(x = var_range, y = pd_values)
      
      p <- ggplot(pd_df, aes(x = x, y = y)) +
        geom_line(color = "steelblue", size = 1.2) +
        theme_bw() +
        labs(
          title = paste("Effect of", var),
          x = var,
          y = paste(gas, "prediction")
        )
      
    } else {
      # Categorical variable
      var_levels <- levels(model_data[[var]])
      
      pd_values <- sapply(var_levels, function(level) {
        temp_data <- model_data
        temp_data[[var]] <- factor(level, levels = var_levels)
        mean(predict(rf_models[[gas]], temp_data))
      })
      
      pd_df <- data.frame(x = var_levels, y = pd_values)
      
      p <- ggplot(pd_df, aes(x = x, y = y)) +
        geom_col(fill = "steelblue", alpha = 0.8) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = paste("Effect of", var),
          x = var,
          y = paste(gas, "prediction")
        )
    }
    
    plot_list[[var]] <- p
  }
  
  if (length(plot_list) >= 2) {
    grid.arrange(grobs = plot_list, ncol = 2)
  }
}

# ============================================================================
# KEY FINDINGS
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("KEY FINDINGS SUMMARY\n")
cat(rep("=", 80), "\n\n")

for (gas in gases) {
  cat("###", toupper(gas), "###\n")
  cat("Random Forest Variance Explained:", round(rf_results[[gas]]$var_explained * 100, 1), "%\n")
  cat("Linear Model R²:", round(lm_results[[gas]]$r_squared * 100, 1), "%\n")
  cat("Improvement:", round((rf_results[[gas]]$var_explained - lm_results[[gas]]$r_squared) * 100, 1), "percentage points\n")
  cat("\nTop 3 Most Important Predictors:\n")
  top3 <- head(rf_results[[gas]]$importance, 3)
  for (i in 1:nrow(top3)) {
    cat(sprintf("  %d. %s (%%IncMSE: %.1f)\n", 
                i, top3$Variable[i], top3$`%IncMSE`[i]))
  }
  cat("\n")
}

cat(rep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n")



















# Load required libraries
library(dplyr)
library(ggplot2)
library(randomForest)
library(gridExtra)


# Convert categorical variables to factors
data <- data %>%
  mutate(
    gen = as.factor(gen),
    season = as.factor(season),
    bot = as.factor(bot),
    age = as.factor(age),
    site = as.factor(site)
  )

cat("Predicting log-transformed CH4 using O2, CO2, and other predictors\n\n")

# ============================================================================
# DATA PREPARATION WITH LOG TRANSFORMATION
# ============================================================================

# Create log-transformed CH4 (add small constant to handle zeros)
data <- data %>%
  mutate(log_ch4 = log(ch4 + 1))

cat("CH4 transformation summary:\n")
cat("Original CH4 - Min:", min(data$ch4, na.rm = TRUE), 
    "Max:", max(data$ch4, na.rm = TRUE), 
    "Median:", median(data$ch4, na.rm = TRUE), "\n")
cat("Log(CH4+1) - Min:", min(data$log_ch4, na.rm = TRUE), 
    "Max:", max(data$log_ch4, na.rm = TRUE), 
    "Median:", median(data$log_ch4, na.rm = TRUE), "\n\n")

# ============================================================================
# MODEL 1: STANDARD PREDICTORS ONLY (baseline)
# ============================================================================

cat(rep("=", 80), "\n")
cat("MODEL 1: Log(CH4) ~ Standard Predictors (baseline)\n")
cat(rep("=", 80), "\n\n")

predictor_vars_1 <- c("gen", "season", "lat", "bot", "age", "dbh")
model_data_1 <- data %>%
  select(all_of(c("log_ch4", predictor_vars_1))) %>%
  na.omit()

cat("Sample size:", nrow(model_data_1), "\n")
cat("Predictors:", paste(predictor_vars_1, collapse = ", "), "\n\n")

set.seed(123)
rf_model_1 <- randomForest(
  x = model_data_1[, predictor_vars_1],
  y = model_data_1$log_ch4,
  ntree = 500,
  importance = TRUE,
  mtry = 3,
  keep.forest = TRUE
)

cat("--- Model 1 Performance ---\n")
cat("Out-of-Bag Variance Explained:", round(rf_model_1$rsq[500] * 100, 2), "%\n")
cat("Out-of-Bag RMSE:", round(sqrt(rf_model_1$mse[500]), 4), "\n\n")

importance_1 <- as.data.frame(importance(rf_model_1))
importance_1$Variable <- rownames(importance_1)
importance_1 <- importance_1 %>% arrange(desc(`%IncMSE`))

cat("Variable Importance:\n")
print(importance_1[, c("Variable", "%IncMSE", "IncNodePurity")])

# ============================================================================
# MODEL 2: ADDING O2 AND CO2 AS PREDICTORS
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("MODEL 2: Log(CH4) ~ Standard Predictors + O2 + CO2\n")
cat(rep("=", 80), "\n\n")

predictor_vars_2 <- c("gen", "season", "lat", "bot", "age", "dbh", "o2", "co2")
model_data_2 <- data %>%
  select(all_of(c("log_ch4", predictor_vars_2))) %>%
  na.omit()

cat("Sample size:", nrow(model_data_2), "\n")
cat("Predictors:", paste(predictor_vars_2, collapse = ", "), "\n\n")

set.seed(123)
rf_model_2 <- randomForest(
  x = model_data_2[, predictor_vars_2],
  y = model_data_2$log_ch4,
  ntree = 500,
  importance = TRUE,
  mtry = 3,
  keep.forest = TRUE
)

cat("--- Model 2 Performance ---\n")
cat("Out-of-Bag Variance Explained:", round(rf_model_2$rsq[500] * 100, 2), "%\n")
cat("Out-of-Bag RMSE:", round(sqrt(rf_model_2$mse[500]), 4), "\n\n")

importance_2 <- as.data.frame(importance(rf_model_2))
importance_2$Variable <- rownames(importance_2)
importance_2 <- importance_2 %>% arrange(desc(`%IncMSE`))

cat("Variable Importance:\n")
print(importance_2[, c("Variable", "%IncMSE", "IncNodePurity")])

# ============================================================================
# MODEL COMPARISON
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("MODEL COMPARISON\n")
cat(rep("=", 80), "\n\n")

comparison <- data.frame(
  Model = c("Model 1: Standard only", "Model 2: + O2 & CO2"),
  Var_Explained = c(
    round(rf_model_1$rsq[500] * 100, 2),
    round(rf_model_2$rsq[500] * 100, 2)
  ),
  RMSE = c(
    round(sqrt(rf_model_1$mse[500]), 4),
    round(sqrt(rf_model_2$mse[500]), 4)
  ),
  N = c(nrow(model_data_1), nrow(model_data_2))
)

print(comparison)

improvement <- rf_model_2$rsq[500] - rf_model_1$rsq[500]
cat("\nImprovement from adding O2 & CO2:", round(improvement * 100, 2), "percentage points\n")

# ============================================================================
# COMPARISON WITH LINEAR MODEL
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("LINEAR MODEL COMPARISON\n")
cat(rep("=", 80), "\n\n")

# Linear model with O2 and CO2
lm_model <- lm(log_ch4 ~ gen + season + lat + bot + age + dbh + o2 + co2, 
               data = model_data_2)

cat("Linear Model Summary:\n")
cat("R-squared:", round(summary(lm_model)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(lm_model)$adj.r.squared, 4), "\n")
cat("RMSE:", round(sqrt(mean(residuals(lm_model)^2)), 4), "\n\n")

cat("Random Forest vs Linear:\n")
cat("RF Var Explained:", round(rf_model_2$rsq[500] * 100, 2), "%\n")
cat("LM R-squared:", round(summary(lm_model)$r.squared * 100, 2), "%\n")
cat("Improvement:", round((rf_model_2$rsq[500] - summary(lm_model)$r.squared) * 100, 2), 
    "percentage points\n")

# ============================================================================
# VISUALIZATIONS
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("CREATING VISUALIZATIONS\n")
cat(rep("=", 80), "\n\n")

# Variable importance plot for Model 2
p1 <- ggplot(importance_2, aes(x = reorder(Variable, `%IncMSE`), y = `%IncMSE`)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Variable Importance for Log(CH4)",
    subtitle = "Model with O2 and CO2 included",
    x = NULL,
    y = "% Increase in MSE"
  ) +
  theme(plot.title = element_text(face = "bold"))

print(p1)

# Partial dependence plots for O2 and CO2
cat("\nCreating partial dependence plots for O2 and CO2...\n")

# O2 partial dependence
o2_range <- seq(
  quantile(model_data_2$o2, 0.05, na.rm = TRUE),
  quantile(model_data_2$o2, 0.95, na.rm = TRUE),
  length.out = 100
)

o2_pd <- sapply(o2_range, function(val) {
  temp_data <- model_data_2
  temp_data$o2 <- val
  mean(predict(rf_model_2, temp_data[, predictor_vars_2]))
})

p2 <- ggplot(data.frame(o2 = o2_range, log_ch4 = o2_pd), 
             aes(x = o2, y = log_ch4)) +
  geom_line(color = "steelblue", size = 1.2) +
  theme_bw() +
  labs(
    title = "Effect of O2 on Log(CH4)",
    x = "O2 concentration",
    y = "Predicted Log(CH4)"
  )

print(p2)

# CO2 partial dependence
co2_range <- seq(
  quantile(model_data_2$co2, 0.05, na.rm = TRUE),
  quantile(model_data_2$co2, 0.95, na.rm = TRUE),
  length.out = 100
)

co2_pd <- sapply(co2_range, function(val) {
  temp_data <- model_data_2
  temp_data$co2 <- val
  mean(predict(rf_model_2, temp_data[, predictor_vars_2]))
})

p3 <- ggplot(data.frame(co2 = co2_range, log_ch4 = co2_pd), 
             aes(x = co2, y = log_ch4)) +
  geom_line(color = "steelblue", size = 1.2) +
  theme_bw() +
  labs(
    title = "Effect of CO2 on Log(CH4)",
    x = "CO2 concentration",
    y = "Predicted Log(CH4)"
  )

print(p3)

# Predicted vs Observed
predictions_2 <- predict(rf_model_2)
pred_df <- data.frame(
  observed = model_data_2$log_ch4,
  predicted = predictions_2,
  residual = model_data_2$log_ch4 - predictions_2
)

p4 <- ggplot(pred_df, aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Predicted vs Observed Log(CH4)",
    subtitle = "Random Forest with O2 & CO2",
    x = "Predicted Log(CH4)",
    y = "Observed Log(CH4)"
  )

print(p4)

# Residuals plot
p5 <- ggplot(pred_df, aes(x = predicted, y = residual)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Residuals vs Predicted",
    x = "Predicted Log(CH4)",
    y = "Residuals"
  )

print(p5)

# ============================================================================
# INTERACTION EFFECTS: O2 x CO2
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("EXPLORING O2 x CO2 INTERACTION\n")
cat(rep("=", 80), "\n\n")

# Create a grid of O2 and CO2 values
o2_grid <- seq(
  quantile(model_data_2$o2, 0.1, na.rm = TRUE),
  quantile(model_data_2$o2, 0.9, na.rm = TRUE),
  length.out = 20
)

co2_grid <- seq(
  quantile(model_data_2$co2, 0.1, na.rm = TRUE),
  quantile(model_data_2$co2, 0.9, na.rm = TRUE),
  length.out = 20
)

interaction_grid <- expand.grid(o2 = o2_grid, co2 = co2_grid)

# Predict for the grid (holding other variables at their medians/modes)
temp_data <- model_data_2[rep(1, nrow(interaction_grid)), ]
temp_data$o2 <- interaction_grid$o2
temp_data$co2 <- interaction_grid$co2

interaction_grid$log_ch4_pred <- predict(rf_model_2, temp_data[, predictor_vars_2])

# Create heatmap
p6 <- ggplot(interaction_grid, aes(x = o2, y = co2, fill = log_ch4_pred)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(interaction_grid$log_ch4_pred),
    name = "Log(CH4)"
  ) +
  theme_bw() +
  labs(
    title = "O2 × CO2 Interaction Effect on Log(CH4)",
    x = "O2 concentration",
    y = "CO2 concentration"
  )

print(p6)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n", rep("=", 80), "\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat(rep("=", 80), "\n\n")

cat("1. LOG TRANSFORMATION:\n")
cat("   - Helps normalize the highly skewed CH4 distribution\n")
cat("   - Makes extreme values less influential\n\n")

cat("2. MODEL PERFORMANCE:\n")
cat("   - Baseline (standard predictors):", round(rf_model_1$rsq[500] * 100, 2), "% variance explained\n")
cat("   - With O2 & CO2:", round(rf_model_2$rsq[500] * 100, 2), "% variance explained\n")
cat("   - Improvement:", round(improvement * 100, 2), "percentage points\n\n")

cat("3. TOP PREDICTORS (with O2 & CO2):\n")
top5 <- head(importance_2, 5)
for (i in 1:nrow(top5)) {
  cat(sprintf("   %d. %s (%%IncMSE: %.1f)\n", 
              i, top5$Variable[i], top5$`%IncMSE`[i]))
}

cat("\n4. BIOLOGICAL INTERPRETATION:\n")
if ("o2" %in% head(importance_2$Variable, 3)) {
  cat("   - O2 is a top predictor, suggesting aerobic/anaerobic conditions matter\n")
}
if ("co2" %in% head(importance_2$Variable, 3)) {
  cat("   - CO2 is a top predictor, indicating respiration-methanogenesis coupling\n")
}

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 80), "\n")

