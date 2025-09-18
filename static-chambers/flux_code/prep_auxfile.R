# Chamber Dimensions Assignment with Weather Data Integration
# This script dynamically loads chamber dimensions from CSV files and assigns them to measurements
# Includes Yale-Myers Forest weather station air temperature integration
# Enhanced with comprehensive air temperature decision hierarchy

library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(purrr)
library(ggplot2)

# Function to assign chamber dimensions by loading data from files
assign_chamber_dimensions <- function(flux_file, surface_area_file, additional_vol_file, simplified_vol_file = NULL, weather_file = NULL) {
  
  cat("=== Loading Data Files ===\n")
  
  # Read the flux data
  flux_data <- read_csv(flux_file, show_col_types = FALSE)
  
  # Read dimension files
  surface_area <- read_csv(surface_area_file, show_col_types = FALSE)
  additional_vol <- read_csv(additional_vol_file, show_col_types = FALSE)
  
  # Read volume file if provided
  if (!is.null(simplified_vol_file) && file.exists(simplified_vol_file)) {
    simplified_vol <- read_csv(simplified_vol_file, show_col_types = FALSE)
  } else {
    simplified_vol <- NULL
    cat("Note: Simplified volume file not provided or not found\n")
  }
  
  # Read weather file if provided
  if (!is.null(weather_file) && file.exists(weather_file)) {
    weather_data <- read_csv(weather_file, show_col_types = FALSE)
    cat("Weather data loaded:", nrow(weather_data), "rows\n")
  } else {
    weather_data <- NULL
    cat("Note: Weather file not provided or not found\n")
  }
  
  cat("Files loaded successfully.\n")
  cat("Original flux data:", nrow(flux_data), "rows\n")
  cat("Surface area data:", nrow(surface_area), "rows\n")
  cat("Additional volume data:", nrow(additional_vol), "rows\n")
  if (!is.null(simplified_vol)) {
    cat("Simplified volume data:", nrow(simplified_vol), "rows\n")
  }
  
  # Check Chamber IDs in the flux data
  unique_chambers <- sort(unique(flux_data$`Chamber ID`[!is.na(flux_data$`Chamber ID`)]))
  cat("Chamber IDs found in flux data:", paste(unique_chambers, collapse = ", "), "\n")
  
  cat("\n=== Processing Surface Area Data ===\n")
  
  # Process surface area data
  surface_area_clean <- surface_area %>%
    mutate(
      chamber_series = case_when(
        grepl("^A", `Chamber ID`, ignore.case = TRUE) ~ "A",
        grepl("^B", `Chamber ID`, ignore.case = TRUE) ~ "B", 
        grepl("^C", `Chamber ID`, ignore.case = TRUE) ~ "C",
        grepl("^D", `Chamber ID`, ignore.case = TRUE) ~ "D",
        TRUE ~ `Chamber ID`
      ),
      surface_area_cm2 = as.numeric(`SA cm2`)
    ) %>%
    select(chamber_series, surface_area_cm2) %>%
    filter(!is.na(surface_area_cm2))
  
  cat("Processed surface area data:\n")
  print(surface_area_clean)
  
  cat("\n=== Processing System Volume Data ===\n")
  
  # Process analyzer volumes
  analyzer_volumes <- additional_vol %>%
    mutate(
      analyzer_name = case_when(
        tolower(instrument) == "lgr_mgga" ~ "LGR",
        tolower(instrument) == "lgr" ~ "LGR",
        tolower(instrument) == "picarro" ~ "Picarro",
        TRUE ~ toupper(instrument)
      ),
      analyzer_cell_volume_cm3 = as.numeric(analyzer_cell),
      tubing_volume_cm3 = as.numeric(tubing),
      # No drierite used - filter volume = 0
      filter_volume_cm3 = 0
    ) %>%
    select(analyzer_name, analyzer_cell_volume_cm3, tubing_volume_cm3, filter_volume_cm3)
  
  cat("Processed analyzer volumes:\n")
  print(analyzer_volumes)
  
  cat("\n=== Processing Weather Data ===\n")
  
  # Process weather data if available
  if (!is.null(weather_data)) {
    # Based on your prep_auxfile.R script, process YMF weather data
    weather_clean <- weather_data %>%
      mutate(
        # Convert timestamp to POSIXct (adjust timezone as needed)
        datetime_weather = as.POSIXct(TIMESTAMP, tz = "UTC"),
        # Use average air temperature for chamber calculations
        air_temp_weather = Tair_Avg,
        # Atmospheric pressure (use default if not available in data)
        pressure_weather = ifelse(!is.na(SVPWPa), 101.325, 101.325)  # Default atmospheric pressure
      ) %>%
      filter(!is.na(datetime_weather) & !is.na(air_temp_weather)) %>%
      select(datetime_weather, air_temp_weather, pressure_weather) %>%
      arrange(datetime_weather)
    
    cat("Weather data processed:", nrow(weather_clean), "records\n")
    cat("Weather date range:", as.character(range(weather_clean$datetime_weather)), "\n")
    cat("Temperature range:", round(range(weather_clean$air_temp_weather, na.rm = TRUE), 1), "°C\n")
  } else {
    weather_clean <- NULL
  }
  
  cat("\n=== Processing Chamber Volume Data ===\n")
  
  # Process chamber volumes if available
  chamber_volumes <- NULL
  if (!is.null(simplified_vol)) {
    # Clean column names (remove newlines and spaces)
    names(simplified_vol) <- gsub("\\n|\\r", "_", names(simplified_vol))
    names(simplified_vol) <- gsub("\\s+", "_", names(simplified_vol))
    
    # Extract chamber volumes
    chamber_volumes <- simplified_vol %>%
      filter(!is.na(`Chamber_Alt_ID`) & !is.na(`Total_Volume_(mL)`)) %>%
      mutate(
        chamber_series = `Chamber_Alt_ID`,
        chamber_volume_cm3 = as.numeric(`Total_Volume_(mL)`)  # mL = cm3
      ) %>%
      filter(!is.na(chamber_volume_cm3) & chamber_series %in% c("A", "B", "C", "D")) %>%
      group_by(chamber_series) %>%
      summarise(
        avg_chamber_volume_cm3 = mean(chamber_volume_cm3, na.rm = TRUE),
        min_volume = min(chamber_volume_cm3, na.rm = TRUE),
        max_volume = max(chamber_volume_cm3, na.rm = TRUE),
        n_measurements = n(),
        .groups = 'drop'
      )
    
    cat("Chamber volumes from simplified_vol data:\n")
    print(chamber_volumes)
  } else {
    cat("No chamber volume data available - will use estimated values\n")
  }
  
  cat("\n=== Creating Chamber Mapping ===\n")
  
  # Create chamber mapping table based on confirmed relationships:
  # 1&2=B, 3&4=C, 5&6=D
  chamber_mapping <- tibble(
    chamber_id = c(1, 2, 3, 4, 5, 6),
    chamber_series = c("B", "B", "C", "C", "D", "D"),
    chamber_pair = c("pair_1", "pair_1", "pair_2", "pair_2", "pair_3", "pair_3")
  )
  
  cat("Chamber ID to series mapping:\n")
  print(chamber_mapping)
  
  # Join mapping with surface areas
  chamber_dimensions <- chamber_mapping %>%
    left_join(surface_area_clean, by = "chamber_series")
  
  # Join with volumes if available
  if (!is.null(chamber_volumes)) {
    chamber_dimensions <- chamber_dimensions %>%
      left_join(chamber_volumes %>% select(chamber_series, avg_chamber_volume_cm3), 
                by = "chamber_series") %>%
      rename(chamber_volume_cm3 = avg_chamber_volume_cm3)
  } else {
    # Use estimated volumes based on relative surface areas
    # B series as baseline, scale others proportionally
    b_surface_area <- surface_area_clean$surface_area_cm2[surface_area_clean$chamber_series == "B"]
    estimated_b_volume <- 590  # From previous observations
    
    chamber_dimensions <- chamber_dimensions %>%
      mutate(
        chamber_volume_cm3 = case_when(
          chamber_series == "B" ~ estimated_b_volume,
          chamber_series == "C" ~ estimated_b_volume * (surface_area_cm2 / b_surface_area),
          chamber_series == "D" ~ estimated_b_volume * (surface_area_cm2 / b_surface_area),
          chamber_series == "A" ~ estimated_b_volume * (surface_area_cm2 / b_surface_area),
          TRUE ~ NA_real_
        )
      )
  }
  
  cat("\nFinal chamber dimensions table:\n")
  print(chamber_dimensions)
  
  cat("\n=== Processing Flux Measurements ===\n")
  
  # Process the flux data
  result <- flux_data %>%
    # Clean up the date format (210806 -> 2021-08-06)
    mutate(
      # Convert 6-digit date format (YYMMDD) to proper date
      date_clean = case_when(
        !is.na(Date) & Date >= 210000 & Date <= 220000 ~ {
          year_part <- floor(Date / 10000) + 2000  # 21 -> 2021
          month_part <- floor((Date %% 10000) / 100)
          day_part <- Date %% 100
          as.Date(paste(year_part, sprintf("%02d", month_part), sprintf("%02d", day_part), sep = "-"))
        },
        TRUE ~ as.Date(NA)
      ),
      
      # Create datetime from date and uncorrected start time for weather matching
      datetime_flux = case_when(
        !is.na(date_clean) & !is.na(`Flux Start Time`) ~ {
          # Convert HHMMSS format to HH:MM:SS
          time_str <- sprintf("%06d", `Flux Start Time`)  # Pad with zeros
          time_formatted <- paste0(substr(time_str, 1, 2), ":", 
                                   substr(time_str, 3, 4), ":", 
                                   substr(time_str, 5, 6))
          as.POSIXct(paste(as.character(date_clean), time_formatted), 
                     format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
        },
        TRUE ~ as.POSIXct(NA)
      ),
      
      # Create unique flux ID
      flux_id = paste(
        format(date_clean, "%Y%m%d"),
        Plot,
        Tree_ID,
        `Measurement height`,
        `Chamber ID`,
        `Flux Start Time`,
        sep = "_"
      ),
      
      # Clean analyzer names - assign all to LGR since all measurements are from LGR3
      analyzer_clean = "LGR",  # Force all to LGR since you confirmed all are LGR3
      
      # Parse start and end times
      start_time = `Flux Start Time`,
      end_time = HMS_End
    ) %>%
    # Join with chamber dimensions
    left_join(chamber_dimensions, by = c("Chamber ID" = "chamber_id")) %>%
    # Join with system volumes
    left_join(analyzer_volumes, by = c("analyzer_clean" = "analyzer_name"))
  
  # Add weather data if available
  if (!is.null(weather_clean)) {
    cat("Merging weather data with flux measurements...\n")
    
    # More efficient approach: use a vectorized function to find closest weather data
    find_closest_weather <- function(flux_datetime, weather_data) {
      if (is.na(flux_datetime)) {
        return(data.frame(air_temp_weather = NA, pressure_weather = NA, time_diff_hours = NA))
      }
      
      time_diffs <- abs(as.numeric(weather_data$datetime_weather - flux_datetime))
      closest_idx <- which.min(time_diffs)
      closest_time_diff_hours <- min(time_diffs) / 3600
      
      # Only use weather data if within 6 hours (not 24)
      if (closest_time_diff_hours <= 6) {
        return(data.frame(
          air_temp_weather = weather_data$air_temp_weather[closest_idx],
          pressure_weather = weather_data$pressure_weather[closest_idx],
          time_diff_hours = closest_time_diff_hours
        ))
      } else {
        return(data.frame(air_temp_weather = NA, pressure_weather = NA, time_diff_hours = NA))
      }
    }
    
    # Apply weather matching to all rows
    weather_matches <- map_dfr(result$datetime_flux, ~find_closest_weather(.x, weather_clean))
    
    # Add weather data to results
    result <- result %>%
      bind_cols(weather_matches) %>%
      rename(
        met_air_temp = air_temp_weather,  # Rename to be clear it's met station data
        met_pressure = pressure_weather,
        met_time_diff_hours = time_diff_hours
      )
    
    # Report weather data usage
    weather_matched <- sum(!is.na(result$met_air_temp))
    
    cat("Weather data integration summary:\n")
    cat("  Met station temperatures matched (within 6 hours):", weather_matched, "of", nrow(result), "\n")
    if (weather_matched > 0) {
      cat("  Met temperature range:", round(range(result$met_air_temp, na.rm = TRUE), 1), "°C\n")
      cat("  Average time difference:", round(mean(result$met_time_diff_hours, na.rm = TRUE), 1), "hours\n")
    }
    
  } else {
    # No weather data - add empty columns
    result <- result %>%
      mutate(
        met_air_temp = NA_real_,
        met_pressure = NA_real_,
        met_time_diff_hours = NA_real_
      )
  }
  
  # Continue with volume calculations
  result <- result %>%
    # Calculate total system volume
    mutate(
      # Total system volume = chamber + tubing + analyzer cell + filter
      total_system_volume_cm3 = case_when(
        !is.na(chamber_volume_cm3) & !is.na(tubing_volume_cm3) & 
          !is.na(analyzer_cell_volume_cm3) & !is.na(filter_volume_cm3) ~
          chamber_volume_cm3 + tubing_volume_cm3 + analyzer_cell_volume_cm3 + filter_volume_cm3,
        TRUE ~ NA_real_
      ),
      
      # Convert to liters
      total_system_volume_L = total_system_volume_cm3 / 1000
    ) %>%
    # Select and rename columns for consistency
    select(
      # Original identifiers
      date = date_clean,
      datetime = datetime_flux,
      personnel = Personnel,
      plot = Plot,
      tree_id = Tree_ID,
      species = Tree_sp,
      measurement_height = `Measurement height`,
      chamber_id = `Chamber ID`,
      temp_stem = Temp_stem,
      temp_air = Temp_air,  # Keep original field air temperature
      dbh = DBH,
      start_time,
      end_time,
      analyzer = analyzer_clean,
      gas_analyser_original = Gas_analyser,
      
      # Generated IDs
      flux_id,
      
      # Chamber characteristics
      chamber_series,
      chamber_pair,
      surface_area_cm2,
      chamber_volume_cm3,
      
      # System volumes
      analyzer_cell_volume_cm3,
      tubing_volume_cm3,
      filter_volume_cm3,
      total_system_volume_cm3,
      total_system_volume_L,
      
      # Met station data (separate columns, no fallbacks)
      met_air_temp,
      met_pressure,
      met_time_diff_hours
    )
  
  cat("\n=== Creating Air Temperature Final Column ===\n")
  
  # Sort data by datetime for temporal operations
  result <- result %>% arrange(datetime)
  
  # Step 1: Create stem-air temperature regression for predictions
  cat("Building stem-air temperature regression...\n")
  
  # Get data with both stem and air temperatures
  regression_data <- result %>%
    filter(!is.na(temp_stem) & !is.na(temp_air))
  
  stem_air_model <- NULL
  if (nrow(regression_data) >= 3) {  # Need at least 3 points for regression
    stem_air_model <- lm(temp_air ~ temp_stem, data = regression_data)
    cat("Stem-air regression built with", nrow(regression_data), "observations\n")
    cat("Model: Air_temp =", round(coef(stem_air_model)[1], 2), "+", 
        round(coef(stem_air_model)[2], 2), "* Stem_temp\n")
    cat("R-squared:", round(summary(stem_air_model)$r.squared, 3), "\n")
  } else {
    cat("Insufficient data for stem-air regression (need ≥3 paired observations)\n")
  }
  
  # Function to find observations within 1 hour window
  find_within_hour <- function(target_datetime, all_datetimes, all_values) {
    if (is.na(target_datetime)) return(NA_real_)
    
    time_diffs <- abs(as.numeric(all_datetimes - target_datetime)) / 3600  # Convert to hours
    within_hour_idx <- which(time_diffs <= 1 & !is.na(all_values))
    
    if (length(within_hour_idx) > 0) {
      return(mean(all_values[within_hour_idx], na.rm = TRUE))
    } else {
      return(NA_real_)
    }
  }
  
  # Function to predict air temp from stem temp using regression
  predict_air_from_stem <- function(stem_temp, model) {
    if (is.na(stem_temp) || is.null(model)) return(NA_real_)
    tryCatch({
      predict(model, newdata = data.frame(temp_stem = stem_temp))
    }, error = function(e) NA_real_)
  }
  
  # Apply the decision hierarchy
  cat("Applying temperature decision hierarchy...\n")
  
  result <- result %>%
    mutate(
      # Step 1: Use field-measured air temp if available
      air_temp_final = case_when(
        !is.na(temp_air) ~ temp_air,
        TRUE ~ NA_real_
      ),
      temp_source = case_when(
        !is.na(temp_air) ~ "field_air",
        TRUE ~ NA_character_
      )
    )
  
  # Step 2: For missing values, try mean of field air temps within 1 hour
  missing_idx <- which(is.na(result$air_temp_final))
  cat("Step 2: Finding field air temps within 1 hour for", length(missing_idx), "missing values...\n")
  
  for (i in missing_idx) {
    hourly_mean <- find_within_hour(result$datetime[i], result$datetime, result$temp_air)
    if (!is.na(hourly_mean)) {
      result$air_temp_final[i] <- hourly_mean
      result$temp_source[i] <- "field_air_1hr_mean"
    }
  }
  
  # Step 3: For remaining missing values, use stem temp regression
  still_missing_idx <- which(is.na(result$air_temp_final))
  cat("Step 3: Using stem temp regression for", length(still_missing_idx), "remaining missing values...\n")
  
  if (!is.null(stem_air_model)) {
    for (i in still_missing_idx) {
      if (!is.na(result$temp_stem[i])) {
        predicted_temp <- predict_air_from_stem(result$temp_stem[i], stem_air_model)
        if (!is.na(predicted_temp)) {
          result$air_temp_final[i] <- predicted_temp
          result$temp_source[i] <- "stem_regression"
        }
      }
    }
  }
  
  # Step 4: For remaining missing values, try mean of stem regression values within 1 hour
  still_missing_idx <- which(is.na(result$air_temp_final))
  cat("Step 4: Finding stem regression temps within 1 hour for", length(still_missing_idx), "missing values...\n")
  
  if (!is.null(stem_air_model)) {
    # Calculate all possible stem regression predictions
    stem_predicted_temps <- map_dbl(result$temp_stem, ~predict_air_from_stem(.x, stem_air_model))
    
    for (i in still_missing_idx) {
      hourly_mean_stem <- find_within_hour(result$datetime[i], result$datetime, stem_predicted_temps)
      if (!is.na(hourly_mean_stem)) {
        result$air_temp_final[i] <- hourly_mean_stem
        result$temp_source[i] <- "stem_regression_1hr_mean"
      }
    }
  }
  
  # Step 5: For remaining missing values, use met station air temp
  final_missing_idx <- which(is.na(result$air_temp_final))
  cat("Step 5: Using met station air temp for", length(final_missing_idx), "final missing values...\n")
  
  for (i in final_missing_idx) {
    if (!is.na(result$met_air_temp[i])) {
      result$air_temp_final[i] <- result$met_air_temp[i]
      result$temp_source[i] <- "met_station"
    }
  }
  
  # Final summary of temperature sources
  cat("\n=== Air Temperature Final Summary ===\n")
  source_summary <- table(result$temp_source, useNA = "ifany")
  print(source_summary)
  
  final_missing <- sum(is.na(result$air_temp_final))
  cat("Final missing air temperatures:", final_missing, "of", nrow(result), "\n")
  
  if (final_missing > 0) {
    cat("Rows still missing air temperature:\n")
    missing_rows <- result[is.na(result$air_temp_final), c("flux_id", "temp_air", "temp_stem", "met_air_temp")]
    print(head(missing_rows, 10))
  }
  
  # Temperature range by source
  cat("\nTemperature ranges by source:\n")
  for (source in names(source_summary)) {
    if (!is.na(source)) {
      source_temps <- result$air_temp_final[result$temp_source == source & !is.na(result$temp_source)]
      if (length(source_temps) > 0) {
        cat(source, ":", round(range(source_temps, na.rm = TRUE), 1), "°C (n =", length(source_temps), ")\n")
      }
    }
  }
  
  cat("\n=== Results Summary ===\n")
  cat("Processed measurements:", nrow(result), "\n")
  
  # Check completeness
  complete_dimensions <- sum(!is.na(result$surface_area_cm2) & !is.na(result$total_system_volume_L))
  cat("Measurements with complete dimensions:", complete_dimensions, "of", nrow(result), "\n")
  
  missing_chambers <- sum(is.na(result$surface_area_cm2))
  if (missing_chambers > 0) {
    cat("Measurements missing chamber dimensions:", missing_chambers, "\n")
    missing_chamber_ids <- unique(result$chamber_id[is.na(result$surface_area_cm2)])
    cat("Missing Chamber IDs:", paste(missing_chamber_ids[!is.na(missing_chamber_ids)], collapse = ", "), "\n")
  }
  
  missing_analyzers <- sum(is.na(result$analyzer_cell_volume_cm3))
  if (missing_analyzers > 0) {
    cat("Measurements missing analyzer dimensions:", missing_analyzers, "\n")
    missing_analyzer_types <- unique(result$analyzer[is.na(result$analyzer_cell_volume_cm3)])
    cat("Missing analyzer types:", paste(missing_analyzer_types[!is.na(missing_analyzer_types)], collapse = ", "), "\n")
  }
  
  # Distribution summaries
  cat("\nChamber distribution:\n")
  chamber_counts <- table(result$chamber_id, useNA = "ifany")
  print(chamber_counts)
  
  cat("\nAnalyzer distribution:\n")
  analyzer_counts <- table(result$analyzer, useNA = "ifany")
  print(analyzer_counts)
  
  cat("\nDate distribution:\n")
  date_counts <- table(format(result$date, "%Y-%m-%d"), useNA = "ifany")
  print(head(date_counts, 10))
  
  # Show examples by chamber series
  cat("\n=== Sample Results by Chamber Series ===\n")
  for (series in c("B", "C", "D")) {
    series_data <- result[result$chamber_series == series & !is.na(result$chamber_series), ]
    if (nrow(series_data) > 0) {
      cat("\n", series, "series chambers (first 2 examples):\n")
      sample_cols <- c("flux_id", "chamber_id", "chamber_series", "surface_area_cm2", 
                       "total_system_volume_L", "temp_air", "air_temp_final", "temp_source")
      print(head(series_data[, sample_cols], 2))
    }
  }
  
  return(result)
}

# Execute the function with error handling
tryCatch({
  cat("Starting chamber dimension assignment with comprehensive air temperature hierarchy...\n")
  
  # File paths - update these to match your actual file locations
  base_path <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/dims"
  
  result <- assign_chamber_dimensions(
    flux_file = '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/field/Gewirtzman_datasheets_fixing_dates.csv',
    surface_area_file = file.path(base_path, "surface_area.csv"),
    additional_vol_file = file.path(base_path, "additional_vol.csv"),
    simplified_vol_file = file.path(base_path, "simplified_volume.csv"),  # Optional
    weather_file = "/Users/jongewirtzman/YSE Dropbox/Jonathan Gewirtzman/Yale-Myers Weather Station Data/ymf_clean_sorted.csv"  # YMF weather data
  )
  
  # Save results
  write_csv(result, "flux_data_with_dimensions_weather_and_final_airtemp.csv")
  cat("\nResults saved to: flux_data_with_dimensions_weather_and_final_airtemp.csv\n")
  
  # Create visualizations
  cat("\nCreating visualizations...\n")
  
  # Plot surface area by chamber series
  p1 <- ggplot(result, aes(x = factor(chamber_series), y = surface_area_cm2, color = factor(chamber_id))) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
    labs(title = "Surface Area by Chamber Series",
         x = "Chamber Series", 
         y = "Surface Area (cm²)",
         color = "Chamber ID") +
    theme_minimal()
  
  print(p1)
  
  # Plot total volume by chamber series
  p2 <- ggplot(result, aes(x = factor(chamber_series), y = total_system_volume_L, color = factor(chamber_id))) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
    labs(title = "Total System Volume by Chamber Series",
         x = "Chamber Series",
         y = "Total System Volume (L)",
         color = "Chamber ID") +
    theme_minimal()
  
  print(p2)
  
  # Plot field measured air temp vs. station measured air temp, colored by plot
  temp_comparison_data <- result %>%
    filter(!is.na(temp_air) & !is.na(met_air_temp))
  
  if (nrow(temp_comparison_data) > 0) {
    p3 <- ggplot(temp_comparison_data, aes(y = temp_air, x = met_air_temp, color = plot)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      labs(title = "Field vs. Weather Station Air Temperature",
           subtitle = paste("Comparison for", nrow(temp_comparison_data), "measurements with both field and station data"),
           y = "Field Measured Air Temperature (°C)",
           x = "Weather Station Air Temperature (°C)",
           color = "Plot") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      coord_equal()
    
    print(p3)
    
    # Calculate correlation and bias statistics
    correlation <- cor(temp_comparison_data$temp_air, temp_comparison_data$met_air_temp, use = "complete.obs")
    bias <- mean(temp_comparison_data$temp_air - temp_comparison_data$met_air_temp, na.rm = TRUE)
    rmse <- sqrt(mean((temp_comparison_data$temp_air - temp_comparison_data$met_air_temp)^2, na.rm = TRUE))
    
    cat("\n=== Temperature Comparison Statistics ===\n")
    cat("Number of paired measurements:", nrow(temp_comparison_data), "\n")
    cat("Correlation (r):", round(correlation, 3), "\n")
    cat("Mean bias (Field - Station):", round(bias, 2), "°C\n")
    cat("RMSE:", round(rmse, 2), "°C\n")
    cat("Field temp range:", round(range(temp_comparison_data$temp_air, na.rm = TRUE), 1), "°C\n")
    cat("Station temp range:", round(range(temp_comparison_data$met_air_temp, na.rm = TRUE), 1), "°C\n")
    
  } else {
    cat("\nNo measurements with both field and station air temperature data for comparison plot.\n")
  }
  
  # Plot air temperature final by source
  if (!all(is.na(result$air_temp_final))) {
    p4 <- ggplot(result, aes(x = temp_source, y = air_temp_final, fill = temp_source)) +
      geom_boxplot(alpha = 0.7) +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.4, size = 0.8) +
      labs(title = "Final Air Temperature by Data Source",
           x = "Temperature Data Source",
           y = "Final Air Temperature (°C)",
           fill = "Source") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none")
    
    print(p4)
  }
  
  # Show stem vs air temperature relationship if regression was built
  stem_air_data <- result %>% filter(!is.na(temp_stem) & !is.na(temp_air))
  if (nrow(stem_air_data) > 0) {
    p5 <- ggplot(stem_air_data, aes(x = temp_stem, y = temp_air)) +
      geom_point(alpha = 0.6, color = "darkblue") +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = "Stem vs. Air Temperature Relationship",
           subtitle = paste("Used for regression predictions (n =", nrow(stem_air_data), ")"),
           x = "Stem Temperature (°C)",
           y = "Field Air Temperature (°C)") +
      theme_minimal()
    
    print(p5)
  }
  
  cat("\nStep 1 completed successfully!\n")
  cat("Chamber dimensions loaded dynamically from CSV files\n")
  cat("Weather station air temperature integrated (6-hour window)\n")
  cat("Air temperature final column created with comprehensive decision hierarchy\n")
  cat("Next step: Create auxfile format for goFlux\n")
  
}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
  cat("Please check your input files and file paths.\n")
  print(traceback())
})


# Create goFlux-compatible auxfile
cat("\n=== Creating goFlux-compatible auxfile ===\n")

auxfile_goflux <- result %>%
  filter(!is.na(air_temp_final) & !is.na(total_system_volume_L) & !is.na(surface_area_cm2)) %>%
  mutate(
    # Required goFlux columns
    UniqueID = flux_id,
    start.time = datetime,
    Area = surface_area_cm2,  # cm²
    Vtot = total_system_volume_L,  # L
    Tcham = air_temp_final,  # °C
    Pcham = case_when(
      !is.na(met_pressure) ~ met_pressure,
      TRUE ~ 101.325  # Default atmospheric pressure in kPa
    ),
    
    # Optional but useful columns
    obs.length = case_when(
      !is.na(end_time) & !is.na(start_time) ~ 
        as.numeric(difftime(end_time, start_time, units = "secs")),
      TRUE ~ 180  # Default 3 minutes if not available
    ),
    
    # Additional metadata (optional)
    plot = plot,
    tree_id = tree_id,
    species = species,
    chamber_id = chamber_id,
    measurement_height = measurement_height,
    temp_source = temp_source
  ) %>%
  select(
    # Core goFlux requirements
    UniqueID, start.time, Area, Vtot, Tcham, Pcham,
    # Optional
    obs.length,
    # Metadata
    plot, tree_id, species, chamber_id, measurement_height, temp_source
  ) %>%
  arrange(start.time)

# Validate auxfile
cat("goFlux auxfile validation:\n")
cat("  Total measurements:", nrow(auxfile_goflux), "\n")
cat("  Complete measurements (all required columns):", 
    sum(complete.cases(auxfile_goflux[, c("UniqueID", "start.time", "Area", "Vtot", "Tcham", "Pcham")])), "\n")

# Check for required columns
required_cols <- c("UniqueID", "start.time", "Area", "Vtot", "Tcham", "Pcham")
missing_cols <- required_cols[!required_cols %in% names(auxfile_goflux)]
if (length(missing_cols) > 0) {
  cat("  WARNING: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
} else {
  cat("  All required columns present ✓\n")
}

# Check data ranges
cat("  Temperature range:", round(range(auxfile_goflux$Tcham, na.rm = TRUE), 1), "°C\n")
cat("  Pressure range:", round(range(auxfile_goflux$Pcham, na.rm = TRUE), 1), "kPa\n")
cat("  Area range:", round(range(auxfile_goflux$Area, na.rm = TRUE), 1), "cm²\n")
cat("  Volume range:", round(range(auxfile_goflux$Vtot, na.rm = TRUE), 3), "L\n")

# Save auxfile
write_csv(auxfile_goflux, "goflux_auxfile.csv")
cat("\ngoFlux-compatible auxfile saved to: goflux_auxfile.csv\n")

# Display sample
cat("\nSample auxfile (first 3 rows):\n")
print(head(auxfile_goflux, 3))


p3 <- ggplot(temp_comparison_data, aes(y = temp_stem, x = temp_air, color = plot)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Field vs. Weather Station Air Temperature",
       color = "Plot") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_equal()+
  geom_smooth(method="lm")

print(p3)
