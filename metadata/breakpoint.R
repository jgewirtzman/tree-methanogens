# Breakpoint Analysis for VWC Thresholds
# Test for data-driven moisture breakpoints rather than arbitrary thresholds

library(dplyr)
library(ggplot2)
# Breakpoint Analysis for VWC Thresholds
# Test for data-driven moisture breakpoints rather than arbitrary thresholds

library(dplyr)
library(ggplot2)
library(segmented)  # Segmented regression
library(strucchange) # Structural change tests
library(gridExtra)

# Load transformed data (assuming from previous analysis)
if (!exists("transformed_data")) {
  stop("Need transformed_data from previous analysis")
}

cat("=== BREAKPOINT ANALYSIS FOR VWC THRESHOLDS ===\n")

# =============================================================================
# 1. VISUAL EXPLORATION OF VWC RELATIONSHIPS
# =============================================================================

cat("\n=== VISUAL EXPLORATION ===\n")

# Create dataset for breakpoint analysis
breakpoint_data <- transformed_data %>%
  filter(!is.na(VWC_mean), !is.na(CH4_flux_target), !is.na(mcra_loose_target)) %>%
  arrange(VWC_mean) %>%
  mutate(
    CH4_flux_log = log(CH4_flux_target - min(CH4_flux_target, na.rm = TRUE) + 0.001),
    mcra_loose_log = log(mcra_loose_target + 0.001)
  )

cat("Complete cases for breakpoint analysis:", nrow(breakpoint_data), "\n")

# Plot relationships to visualize potential breakpoints
p1 <- ggplot(breakpoint_data, aes(x = VWC_mean, y = CH4_flux_log)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "green", alpha = 0.7) +
  labs(title = "CH4 Flux vs VWC", x = "VWC (%)", y = "Log CH4 Flux") +
  theme_minimal()

p2 <- ggplot(breakpoint_data, aes(x = VWC_mean, y = mcra_loose_log)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "green", alpha = 0.7) +
  labs(title = "mcrA Abundance vs VWC", x = "VWC (%)", y = "Log mcrA") +
  theme_minimal()

combined_plot <- grid.arrange(p1, p2, ncol = 2, top = "VWC Relationships (Green line = arbitrary 20% threshold)")
print(combined_plot)

# =============================================================================
# 2. SEGMENTED REGRESSION BREAKPOINT DETECTION
# =============================================================================

cat("\n=== SEGMENTED REGRESSION BREAKPOINT DETECTION ===\n")

# Function to detect breakpoints using segmented regression
detect_breakpoint_segmented <- function(x, y, var_name) {
  # Remove NAs
  complete_cases <- !is.na(x) & !is.na(y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  if (length(x_clean) < 20) {
    return(list(breakpoint = NA, pvalue = NA, method = "Insufficient data"))
  }
  
  tryCatch({
    # Fit linear model first
    lm_model <- lm(y_clean ~ x_clean)
    
    # Try to detect breakpoint (provide initial guess based on data range)
    initial_guess <- median(x_clean)
    seg_model <- segmented(lm_model, seg.Z = ~x_clean, psi = initial_guess)
    
    breakpoint <- seg_model$psi[1, "Est."]
    std_error <- seg_model$psi[1, "St.Err"]
    
    return(list(
      breakpoint = round(breakpoint, 2), 
      std_error = round(std_error, 2),
      method = "Segmented regression",
      model = seg_model
    ))
  }, error = function(e) {
    return(list(breakpoint = NA, std_error = NA, method = paste("Error:", e$message)))
  })
}

# Detect breakpoints
ch4_breakpoint <- detect_breakpoint_segmented(breakpoint_data$VWC_mean, 
                                              breakpoint_data$CH4_flux_log, 
                                              "CH4_flux")

mcra_breakpoint <- detect_breakpoint_segmented(breakpoint_data$VWC_mean, 
                                               breakpoint_data$mcra_loose_log, 
                                               "mcrA_loose")

cat("CH4 Flux breakpoint results:\n")
print(ch4_breakpoint[1:3])

cat("\nmcrA Loose breakpoint results:\n")
print(mcra_breakpoint[1:3])

# =============================================================================
# 3. STRUCTURAL CHANGE TEST
# =============================================================================

cat("\n=== STRUCTURAL CHANGE TESTS ===\n")

# Function for structural change tests
detect_structural_change <- function(x, y, var_name) {
  # Remove NAs
  complete_cases <- !is.na(x) & !is.na(y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  if (length(x_clean) < 20) {
    return(list(breakpoint = NA, pvalue = NA, method = "Insufficient data"))
  }
  
  tryCatch({
    # Create data frame for strucchange
    df <- data.frame(y = y_clean, x = x_clean)
    df <- df[order(df$x), ]
    
    # Structural change test
    ts_data <- ts(df)
    bp_test <- breakpoints(y ~ x, data = df, h = 0.15)  # h = minimum segment size
    
    if (!is.na(bp_test$breakpoints[1])) {
      # Get the x-value corresponding to the breakpoint
      bp_index <- bp_test$breakpoints[1]
      breakpoint_x <- df$x[bp_index]
      
      return(list(
        breakpoint = round(breakpoint_x, 2),
        n_breakpoints = length(bp_test$breakpoints[!is.na(bp_test$breakpoints)]),
        method = "Structural change"
      ))
    } else {
      return(list(breakpoint = NA, n_breakpoints = 0, method = "No structural change detected"))
    }
    
  }, error = function(e) {
    return(list(breakpoint = NA, n_breakpoints = NA, method = paste("Error:", e$message)))
  })
}

# Apply structural change tests
ch4_struc <- detect_structural_change(breakpoint_data$VWC_mean, 
                                      breakpoint_data$CH4_flux_log, 
                                      "CH4_flux")

mcra_struc <- detect_structural_change(breakpoint_data$VWC_mean, 
                                       breakpoint_data$mcra_loose_log, 
                                       "mcrA_loose")

cat("CH4 Flux structural change:\n")
print(ch4_struc[1:3])

cat("\nmcrA Loose structural change:\n")
print(mcra_struc[1:3])

# =============================================================================
# 4. SLIDING WINDOW ANALYSIS (F-STATISTIC OPTIMIZATION)
# =============================================================================

cat("\n=== SLIDING WINDOW ANALYSIS ===\n")

# Function to test multiple thresholds
test_multiple_thresholds <- function(data, response_var, thresholds) {
  results <- data.frame()
  
  for (threshold in thresholds) {
    # Split data
    dry_data <- data %>% filter(VWC_mean < threshold)
    wet_data <- data %>% filter(VWC_mean >= threshold)
    
    # Ensure minimum sample sizes
    if (nrow(dry_data) < 10 || nrow(wet_data) < 10) {
      next
    }
    
    # Calculate means and group differences
    dry_mean <- mean(dry_data[[response_var]], na.rm = TRUE)
    wet_mean <- mean(wet_data[[response_var]], na.rm = TRUE)
    overall_mean <- mean(data[[response_var]], na.rm = TRUE)
    
    # Calculate between-group and within-group variance
    between_ss <- nrow(dry_data) * (dry_mean - overall_mean)^2 + 
      nrow(wet_data) * (wet_mean - overall_mean)^2
    
    within_ss <- sum((dry_data[[response_var]] - dry_mean)^2, na.rm = TRUE) +
      sum((wet_data[[response_var]] - wet_mean)^2, na.rm = TRUE)
    
    # F-statistic (between-group variance / within-group variance)
    df_between <- 1
    df_within <- nrow(dry_data) + nrow(wet_data) - 2
    
    f_stat <- (between_ss / df_between) / (within_ss / df_within)
    
    # Cohen's d (effect size)
    pooled_sd <- sqrt(within_ss / df_within)
    cohens_d <- abs(dry_mean - wet_mean) / pooled_sd
    
    results <- bind_rows(results, data.frame(
      threshold = threshold,
      n_dry = nrow(dry_data),
      n_wet = nrow(wet_data),
      dry_mean = round(dry_mean, 3),
      wet_mean = round(wet_mean, 3),
      f_statistic = round(f_stat, 3),
      cohens_d = round(cohens_d, 3),
      mean_difference = round(abs(dry_mean - wet_mean), 3)
    ))
  }
  
  return(results)
}

# Test thresholds from 12% to 30% VWC
thresholds_to_test <- seq(12, 30, by = 1)

ch4_threshold_results <- test_multiple_thresholds(
  breakpoint_data, "CH4_flux_log", thresholds_to_test
)

mcra_threshold_results <- test_multiple_thresholds(
  breakpoint_data, "mcra_loose_log", thresholds_to_test
)

# Find optimal thresholds
if (nrow(ch4_threshold_results) > 0) {
  optimal_ch4_threshold <- ch4_threshold_results[which.max(ch4_threshold_results$f_statistic), ]
  cat("Optimal CH4 flux threshold (max F-stat):", optimal_ch4_threshold$threshold, "%\n")
  cat("F-statistic:", optimal_ch4_threshold$f_statistic, 
      "| Cohen's d:", optimal_ch4_threshold$cohens_d, "\n")
}

if (nrow(mcra_threshold_results) > 0) {
  optimal_mcra_threshold <- mcra_threshold_results[which.max(mcra_threshold_results$f_statistic), ]
  cat("Optimal mcrA threshold (max F-stat):", optimal_mcra_threshold$threshold, "%\n")
  cat("F-statistic:", optimal_mcra_threshold$f_statistic,
      "| Cohen's d:", optimal_mcra_threshold$cohens_d, "\n")
}

# Plot threshold analysis results
if (nrow(ch4_threshold_results) > 0) {
  p3 <- ggplot(ch4_threshold_results, aes(x = threshold, y = f_statistic)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue") +
    geom_vline(xintercept = 20, linetype = "dashed", color = "green") +
    geom_vline(xintercept = optimal_ch4_threshold$threshold, linetype = "solid", color = "red") +
    labs(title = "CH4 Flux: Threshold Optimization", 
         x = "VWC Threshold (%)", 
         y = "F-statistic (group separation)",
         subtitle = paste("Optimal threshold:", optimal_ch4_threshold$threshold, "% (red line)")) +
    theme_minimal()
  
  print(p3)
}

if (nrow(mcra_threshold_results) > 0) {
  p4 <- ggplot(mcra_threshold_results, aes(x = threshold, y = f_statistic)) +
    geom_line(color = "purple", size = 1) +
    geom_point(color = "purple") +
    geom_vline(xintercept = 20, linetype = "dashed", color = "green") +
    geom_vline(xintercept = optimal_mcra_threshold$threshold, linetype = "solid", color = "red") +
    labs(title = "mcrA Abundance: Threshold Optimization", 
         x = "VWC Threshold (%)", 
         y = "F-statistic (group separation)",
         subtitle = paste("Optimal threshold:", optimal_mcra_threshold$threshold, "% (red line)")) +
    theme_minimal()
  
  print(p4)
}

# =============================================================================
# 5. SIMPLE CHANGE POINT DETECTION (CUMULATIVE SUM)
# =============================================================================

cat("\n=== SIMPLE CHANGE POINT DETECTION ===\n")

# Simple CUSUM-based change point detection
detect_cusum_changepoint <- function(x, y, var_name) {
  # Remove NAs and sort
  complete_cases <- !is.na(x) & !is.na(y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  # Sort by x
  order_idx <- order(x_clean)
  x_sorted <- x_clean[order_idx]
  y_sorted <- y_clean[order_idx]
  
  if (length(x_sorted) < 20) {
    return(list(breakpoint = NA, method = "Insufficient data"))
  }
  
  n <- length(y_sorted)
  overall_mean <- mean(y_sorted)
  
  # Calculate cumulative sum of deviations
  cusum <- cumsum(y_sorted - overall_mean)
  
  # Find maximum absolute CUSUM (most likely change point)
  max_cusum_idx <- which.max(abs(cusum))
  
  # Don't consider endpoints (need at least 15% of data on each side)
  min_idx <- ceiling(0.15 * n)
  max_idx <- floor(0.85 * n)
  
  if (max_cusum_idx >= min_idx && max_cusum_idx <= max_idx) {
    breakpoint_x <- x_sorted[max_cusum_idx]
    return(list(
      breakpoint = round(breakpoint_x, 2),
      cusum_value = round(cusum[max_cusum_idx], 3),
      method = "CUSUM"
    ))
  } else {
    return(list(breakpoint = NA, method = "No valid change point"))
  }
}

# Apply CUSUM detection
ch4_cusum <- detect_cusum_changepoint(breakpoint_data$VWC_mean, 
                                      breakpoint_data$CH4_flux_log, 
                                      "CH4_flux")

mcra_cusum <- detect_cusum_changepoint(breakpoint_data$VWC_mean, 
                                       breakpoint_data$mcra_loose_log, 
                                       "mcrA_loose")

cat("CH4 Flux CUSUM breakpoint:\n")
print(ch4_cusum[1:3])

cat("\nmcrA Loose CUSUM breakpoint:\n")
print(mcra_cusum[1:3])

# =============================================================================
# 6. SUMMARY OF BREAKPOINT ANALYSES
# =============================================================================

cat("\n=== BREAKPOINT ANALYSIS SUMMARY ===\n")

summary_results <- data.frame(
  Variable = c("CH4_flux", "mcrA_loose"),
  Arbitrary_Threshold = c("20%", "20%"),
  Segmented_Regression = c(
    ifelse(is.na(ch4_breakpoint$breakpoint), "Not detected", paste0(ch4_breakpoint$breakpoint, "%")),
    ifelse(is.na(mcra_breakpoint$breakpoint), "Not detected", paste0(mcra_breakpoint$breakpoint, "%"))
  ),
  Structural_Change = c(
    ifelse(is.na(ch4_struc$breakpoint), "Not detected", paste0(ch4_struc$breakpoint, "%")),
    ifelse(is.na(mcra_struc$breakpoint), "Not detected", paste0(mcra_struc$breakpoint, "%"))
  ),
  CUSUM_Method = c(
    ifelse(is.na(ch4_cusum$breakpoint), "Not detected", paste0(ch4_cusum$breakpoint, "%")),
    ifelse(is.na(mcra_cusum$breakpoint), "Not detected", paste0(mcra_cusum$breakpoint, "%"))
  ),
  Optimal_F_Threshold = c(
    ifelse(nrow(ch4_threshold_results) > 0, paste0(optimal_ch4_threshold$threshold, "%"), "Not calculated"),
    ifelse(nrow(mcra_threshold_results) > 0, paste0(optimal_mcra_threshold$threshold, "%"), "Not calculated")
  )
)

print(summary_results)

cat("\n=== RECOMMENDATIONS ===\n")
cat("1. Compare data-driven thresholds with arbitrary 20% threshold\n")
cat("2. If breakpoints differ significantly, rerun models with optimal thresholds\n")
cat("3. Consider continuous VWC modeling instead of binary splits\n")
cat("4. Test non-linear relationships (splines, GAMs)\n")
cat("5. The F-statistic method shows which threshold maximizes group separation\n")
library(segmented)  # Segmented regression
library(strucchange) # Structural change tests
library(gridExtra)

# Load transformed data (assuming from previous analysis)
if (!exists("transformed_data")) {
  stop("Need transformed_data from previous analysis")
}

cat("=== BREAKPOINT ANALYSIS FOR VWC THRESHOLDS ===\n")

# =============================================================================
# 1. VISUAL EXPLORATION OF VWC RELATIONSHIPS
# =============================================================================

cat("\n=== VISUAL EXPLORATION ===\n")

# Create dataset for breakpoint analysis
breakpoint_data <- transformed_data %>%
  filter(!is.na(VWC_mean), !is.na(CH4_flux_target), !is.na(mcra_loose_target)) %>%
  arrange(VWC_mean) %>%
  mutate(
    CH4_flux_log = log(CH4_flux_target - min(CH4_flux_target, na.rm = TRUE) + 0.001),
    mcra_loose_log = log(mcra_loose_target + 0.001)
  )

cat("Complete cases for breakpoint analysis:", nrow(breakpoint_data), "\n")

# Plot relationships to visualize potential breakpoints
p1 <- ggplot(breakpoint_data, aes(x = VWC_mean, y = CH4_flux_log)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "green", alpha = 0.7) +
  labs(title = "CH4 Flux vs VWC", x = "VWC (%)", y = "Log CH4 Flux") +
  theme_minimal()

p2 <- ggplot(breakpoint_data, aes(x = VWC_mean, y = mcra_loose_log)) +
  geom_point(alpha = 0.6, color = "purple") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "green", alpha = 0.7) +
  labs(title = "mcrA Abundance vs VWC", x = "VWC (%)", y = "Log mcrA") +
  theme_minimal()

combined_plot <- grid.arrange(p1, p2, ncol = 2, top = "VWC Relationships (Green line = arbitrary 20% threshold)")
print(combined_plot)

# =============================================================================
# 2. SEGMENTED REGRESSION BREAKPOINT DETECTION
# =============================================================================

cat("\n=== SEGMENTED REGRESSION BREAKPOINT DETECTION ===\n")

# Function to detect breakpoints using segmented regression
detect_breakpoint_segmented <- function(x, y, var_name) {
  # Remove NAs
  complete_cases <- !is.na(x) & !is.na(y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  if (length(x_clean) < 20) {
    return(list(breakpoint = NA, pvalue = NA, method = "Insufficient data"))
  }
  
  tryCatch({
    # Fit linear model first
    lm_model <- lm(y_clean ~ x_clean)
    
    # Try to detect breakpoint
    seg_model <- segmented(lm_model, seg.Z = ~x_clean, npsi = 1)
    
    breakpoint <- seg_model$psi[1, "Est."]
    pvalue <- seg_model$psi[1, "St.Err"]  # Use std error as proxy
    
    return(list(
      breakpoint = round(breakpoint, 2), 
      std_error = round(pvalue, 2),
      method = "Segmented regression",
      model = seg_model
    ))
  }, error = function(e) {
    return(list(breakpoint = NA, pvalue = NA, method = paste("Error:", e$message)))
  })
}

# Detect breakpoints
ch4_breakpoint <- detect_breakpoint_segmented(breakpoint_data$VWC_mean, 
                                              breakpoint_data$CH4_flux_log, 
                                              "CH4_flux")

mcra_breakpoint <- detect_breakpoint_segmented(breakpoint_data$VWC_mean, 
                                               breakpoint_data$mcra_loose_log, 
                                               "mcrA_loose")

cat("CH4 Flux breakpoint results:\n")
print(ch4_breakpoint[1:3])

cat("\nmcrA Loose breakpoint results:\n")
print(mcra_breakpoint[1:3])

# =============================================================================
# 3. BAYESIAN CHANGE POINT DETECTION
# =============================================================================

cat("\n=== BAYESIAN CHANGE POINT DETECTION ===\n")

# Function for Bayesian change point detection
detect_breakpoint_bayesian <- function(x, y, var_name) {
  # Remove NAs and sort by x
  complete_cases <- !is.na(x) & !is.na(y)
  x_clean <- x[complete_cases]
  y_clean <- y[complete_cases]
  
  # Sort by x values
  order_idx <- order(x_clean)
  x_sorted <- x_clean[order_idx]
  y_sorted <- y_clean[order_idx]
  
  if (length(x_sorted) < 20) {
    return(list(breakpoint = NA, probability = NA, method = "Insufficient data"))
  }
  
  tryCatch({
    # Bayesian change point detection
    bcp_result <- bcp(y_sorted, mcmc = 1000, burnin = 100)
    
    # Find most probable change point
    prob_threshold <- 0.5
    change_points <- which(bcp_result$prob.mean > prob_threshold)
    
    if (length(change_points) > 0) {
      # Take the most probable change point
      best_cp_idx <- change_points[which.max(bcp_result$prob.mean[change_points])]
      breakpoint_vwc <- x_sorted[best_cp_idx]
      probability <- bcp_result$prob.mean[best_cp_idx]
      
      return(list(
        breakpoint = round(breakpoint_vwc, 2),
        probability = round(probability, 3),
        method = "Bayesian change point",
        n_changepoints = length(change_points)
      ))
    } else {
      return(list(breakpoint = NA, probability = NA, method = "No significant change points"))
    }
  }, error = function(e) {
    return(list(breakpoint = NA, probability = NA, method = paste("Error:", e$message)))
  })
}

# Apply Bayesian change point detection
ch4_bcp <- detect_breakpoint_bayesian(breakpoint_data$VWC_mean, 
                                      breakpoint_data$CH4_flux_log, 
                                      "CH4_flux")

mcra_bcp <- detect_breakpoint_bayesian(breakpoint_data$VWC_mean, 
                                       breakpoint_data$mcra_loose_log, 
                                       "mcrA_loose")

cat("CH4 Flux Bayesian breakpoint:\n")
print(ch4_bcp[1:3])

cat("\nmcrA Loose Bayesian breakpoint:\n")
print(mcra_bcp[1:3])

# =============================================================================
# 4. SLIDING WINDOW ANALYSIS
# =============================================================================

cat("\n=== SLIDING WINDOW ANALYSIS ===\n")

# Function to test multiple thresholds
test_multiple_thresholds <- function(data, response_var, predictor_vars, thresholds) {
  results <- data.frame()
  
  for (threshold in thresholds) {
    # Split data
    dry_data <- data %>% filter(VWC_mean < threshold)
    wet_data <- data %>% filter(VWC_mean >= threshold)
    
    # Ensure minimum sample sizes
    if (nrow(dry_data) < 15 || nrow(wet_data) < 15) {
      next
    }
    
    # Calculate variance in response variable for each group
    dry_var <- var(dry_data[[response_var]], na.rm = TRUE)
    wet_var <- var(wet_data[[response_var]], na.rm = TRUE)
    total_var <- var(data[[response_var]], na.rm = TRUE)
    
    # Calculate F-statistic for group differences
    between_var <- var(c(rep(mean(dry_data[[response_var]], na.rm = TRUE), nrow(dry_data)),
                         rep(mean(wet_data[[response_var]], na.rm = TRUE), nrow(wet_data))), na.rm = TRUE)
    
    f_stat <- between_var / ((dry_var + wet_var) / 2)
    
    results <- bind_rows(results, data.frame(
      threshold = threshold,
      n_dry = nrow(dry_data),
      n_wet = nrow(wet_data),
      dry_mean = mean(dry_data[[response_var]], na.rm = TRUE),
      wet_mean = mean(wet_data[[response_var]], na.rm = TRUE),
      f_statistic = f_stat,
      between_group_var = between_var
    ))
  }
  
  return(results)
}

# Test thresholds from 10% to 35% VWC
thresholds_to_test <- seq(10, 35, by = 1)

ch4_threshold_results <- test_multiple_thresholds(
  breakpoint_data, "CH4_flux_log", NULL, thresholds_to_test
)

mcra_threshold_results <- test_multiple_thresholds(
  breakpoint_data, "mcra_loose_log", NULL, thresholds_to_test
)

# Find optimal thresholds (maximum F-statistic)
if (nrow(ch4_threshold_results) > 0) {
  optimal_ch4_threshold <- ch4_threshold_results[which.max(ch4_threshold_results$f_statistic), ]
  cat("Optimal CH4 flux threshold (max F-stat):", optimal_ch4_threshold$threshold, "%\n")
  cat("F-statistic:", round(optimal_ch4_threshold$f_statistic, 3), "\n")
}

if (nrow(mcra_threshold_results) > 0) {
  optimal_mcra_threshold <- mcra_threshold_results[which.max(mcra_threshold_results$f_statistic), ]
  cat("Optimal mcrA threshold (max F-stat):", optimal_mcra_threshold$threshold, "%\n")
  cat("F-statistic:", round(optimal_mcra_threshold$f_statistic, 3), "\n")
}

# Plot threshold analysis results
if (nrow(ch4_threshold_results) > 0) {
  p3 <- ggplot(ch4_threshold_results, aes(x = threshold, y = f_statistic)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue") +
    geom_vline(xintercept = 20, linetype = "dashed", color = "green") +
    labs(title = "CH4 Flux: Threshold Optimization", 
         x = "VWC Threshold (%)", 
         y = "F-statistic (group separation)") +
    theme_minimal()
  
  print(p3)
}

if (nrow(mcra_threshold_results) > 0) {
  p4 <- ggplot(mcra_threshold_results, aes(x = threshold, y = f_statistic)) +
    geom_line(color = "purple", size = 1) +
    geom_point(color = "purple") +
    geom_vline(xintercept = 20, linetype = "dashed", color = "green") +
    labs(title = "mcrA Abundance: Threshold Optimization", 
         x = "VWC Threshold (%)", 
         y = "F-statistic (group separation)") +
    theme_minimal()
  
  print(p4)
}

# =============================================================================
# 5. SUMMARY OF BREAKPOINT ANALYSES
# =============================================================================

cat("\n=== BREAKPOINT ANALYSIS SUMMARY ===\n")

summary_results <- data.frame(
  Variable = c("CH4_flux", "mcrA_loose"),
  Arbitrary_Threshold = c(20, 20),
  Segmented_Breakpoint = c(
    ifelse(is.na(ch4_breakpoint$breakpoint), "Not detected", paste(ch4_breakpoint$breakpoint, "%")),
    ifelse(is.na(mcra_breakpoint$breakpoint), "Not detected", paste(mcra_breakpoint$breakpoint, "%"))
  ),
  Bayesian_Breakpoint = c(
    ifelse(is.na(ch4_bcp$breakpoint), "Not detected", paste(ch4_bcp$breakpoint, "%")),
    ifelse(is.na(mcra_bcp$breakpoint), "Not detected", paste(mcra_bcp$breakpoint, "%"))
  ),
  Optimal_F_Threshold = c(
    ifelse(nrow(ch4_threshold_results) > 0, paste(optimal_ch4_threshold$threshold, "%"), "Not calculated"),
    ifelse(nrow(mcra_threshold_results) > 0, paste(optimal_mcra_threshold$threshold, "%"), "Not calculated")
  )
)

print(summary_results)

cat("\n=== RECOMMENDATIONS ===\n")
cat("1. Compare data-driven thresholds with arbitrary 20% threshold\n")
cat("2. If breakpoints differ significantly, rerun models with optimal thresholds\n")
cat("3. Consider continuous VWC modeling instead of binary splits\n")
cat("4. Test non-linear relationships (splines, GAMs)\n")