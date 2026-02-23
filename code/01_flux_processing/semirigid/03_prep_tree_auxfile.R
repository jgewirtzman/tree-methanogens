# ==============================================================================
# Prepare Tree Flux Auxiliary File
# ==============================================================================
# Purpose: Prepares goFlux auxiliary metadata file for tree flux measurements,
#   adding chamber geometry and weather station data.
#
# Pipeline stage: 01 Flux Processing (semirigid)
# Run after: 02_join_flux_geometry.R
#
# Inputs:
#   - flux_with_geometry_fixed.csv (from 02_join_flux_geometry.R)
#   - ymf_clean_sorted.csv (from data/raw/weather/)
#
# Outputs:
#   - auxfile_goFlux_with_weather.csv
# ==============================================================================

library(dplyr)
library(readr)
library(lubridate)

# Read the CSV file
flux_data <- read_csv("flux_with_geometry_fixed.csv")

# Read the weather data
ymf_met <- read.csv("../../../data/raw/weather/ymf_clean_sorted.csv")

# Prepare weather data
cat("=== PREPARING WEATHER DATA ===\n")
weather_clean <- ymf_met %>%
  mutate(
    # Convert timestamp to POSIXct
    TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "UTC"),
    # Convert temperature from Celsius to match goFlux expectations
    Tcham_weather = Tair_Avg,  # Using average air temperature
    # Convert pressure - note: need to check if SVPWPa is atmospheric pressure
    # SVPWPa appears to be saturation vapor pressure, not atmospheric pressure
    # We'll calculate atmospheric pressure from other variables or use a default
    Pcham_weather = 101.325  # Default for now - may need site-specific correction
  ) %>%
  select(TIMESTAMP, Tcham_weather, Pcham_weather) %>%
  arrange(TIMESTAMP)

cat("Weather data range:", as.character(range(weather_clean$TIMESTAMP)), "\n")
cat("Temperature range:", round(range(weather_clean$Tcham_weather, na.rm = TRUE), 1), "°C\n")

# Check for problematic time values first
cat("\n=== INVESTIGATING TIME FORMAT ISSUES ===\n")
problematic_times <- flux_data %>%
  mutate(numeric_time = as.numeric(`Time of sampling`)) %>%
  filter(is.na(numeric_time) & !is.na(`Time of sampling`)) %>%
  select(Date, `Plot Tag`, `Plot Letter`, `Time of sampling`, `File name`)

if(nrow(problematic_times) > 0) {
  cat("Found", nrow(problematic_times), "rows with non-numeric time values:\n")
  print(problematic_times)
} else {
  cat("No problematic time values found.\n")
}

# Show examples of time format variations
cat("\nExamples of time formats in your data:\n")
time_examples <- flux_data %>%
  filter(!is.na(`Time of sampling`)) %>%
  select(`Time of sampling`) %>%
  distinct() %>%
  head(10)
print(time_examples)

# Check for potential duplicates before processing
potential_duplicates <- flux_data %>%
  mutate(UniqueID_base = paste(Date, `Plot Tag`, `Plot Letter`, Chamber_ID, sep = "_")) %>%
  group_by(UniqueID_base) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1)

if(nrow(potential_duplicates) > 0) {
  cat("\nFound", nrow(potential_duplicates), "sets of duplicate measurements:\n")
  print(potential_duplicates)
}

cat("\n=== CONVERTING TO AUXFILE FORMAT WITH WEATHER DATA ===\n")

# Convert to auxfile format for goFlux
auxfile <- flux_data %>%
  # Create UniqueID by combining relevant identifiers
  mutate(
    # Convert time, handling both decimal and integer HHMM formats
    numeric_time = as.numeric(`Time of sampling`),
    
    # Create a base unique identifier
    UniqueID_base = paste(Date, `Plot Tag`, `Plot Letter`, Chamber_ID, sep = "_"),
    
    # Make truly unique IDs by adding sequence numbers for duplicates
    UniqueID = ave(UniqueID_base, UniqueID_base, FUN = function(x) {
      if(length(x) == 1) {
        return(x)
      } else {
        return(paste(x, seq_along(x), sep = "_"))
      }
    }),
    
    # IMPROVED DATETIME HANDLING FOR HHMM FORMAT (with or without decimals)
    # First calculate hours and minutes for all rows
    time_int = floor(numeric_time),  # Remove decimal part
    hours = floor(time_int / 100),
    minutes = time_int %% 100,
    
    # Convert HHMM time format (e.g., 1357.0 or 1357) to HH:MM:SS
    formatted_time = case_when(
      is.na(numeric_time) ~ NA_character_,  # Keep NA for missing values
      hours >= 0 & hours <= 23 & minutes >= 0 & minutes <= 59 ~ 
        sprintf("%02d:%02d:00", hours, minutes),  # Valid time
      TRUE ~ NA_character_  # Keep NA for invalid times - don't create fake times
    ),
    
    # Combine Date and formatted time - ensuring datetime is preserved
    datetime_string = case_when(
      !is.na(formatted_time) ~ paste(Date, formatted_time),
      TRUE ~ NA_character_  # Keep NA if no valid time
    ),
    start.time = as.POSIXct(datetime_string, 
                            format = "%Y-%m-%d %H:%M:%S", 
                            tz = "UTC"),
    
    # Convert surface area from cm2 (already correct unit)
    Area = surface_area_cm2,
    
    # Calculate total system volume
    # Chamber volume (from your data) + tubing volume + system volume
    # Tubing: 1/8" ID × 12 ft = π × (1/16)² × 12 × 12 = π × (0.0625)² × 144 in³
    # Convert to cm³: × 16.387 
    tubing_volume_cm3 = pi * (1/16)^2 * 12 * 12 * 16.387, # ≈ 18.6 cm³
    system_volume_cm3 = 70, # Given system volume
    
    # Total volume = chamber + tubing + system, convert to L
    Vtot = (volume_cm3 + tubing_volume_cm3 + system_volume_cm3) / 1000
  ) %>%
  # Remove rows with missing critical data before merging
  filter(!is.na(start.time), !is.na(Area), !is.na(Vtot)) %>%
  # Merge with weather data using nearest timestamp
  # This finds the closest weather measurement to each flux measurement
  left_join(
    weather_clean,
    by = c("start.time" = "TIMESTAMP")
  ) %>%
  # If no exact match, find nearest weather data
  mutate(
    # For unmatched rows, we'll use a rolling join approach
    Tcham = case_when(
      !is.na(Tcham_weather) ~ Tcham_weather,
      TRUE ~ NA_real_
    ),
    Pcham = case_when(
      !is.na(Pcham_weather) ~ Pcham_weather,
      TRUE ~ NA_real_
    )
  )

# For rows without weather matches, find nearest weather data
unmatched_rows <- which(is.na(auxfile$Tcham))

if(length(unmatched_rows) > 0) {
  cat("Finding nearest weather data for", length(unmatched_rows), "unmatched measurements...\n")
  
  for(i in unmatched_rows) {
    measurement_time <- auxfile$start.time[i]
    
    # Find closest weather observation
    time_diffs <- abs(as.numeric(weather_clean$TIMESTAMP - measurement_time))
    closest_idx <- which.min(time_diffs)
    closest_time_diff <- min(time_diffs) / 3600  # Convert to hours
    
    # Only use weather data if it's within 24 hours
    if(closest_time_diff <= 24) {
      auxfile$Tcham[i] <- weather_clean$Tcham_weather[closest_idx]
      auxfile$Pcham[i] <- weather_clean$Pcham_weather[closest_idx]
    }
  }
}

# Final processing - NO DEFAULTS for missing times
auxfile <- auxfile %>%
  # Add defaults only for temperature and pressure, not time
  mutate(
    Tcham = case_when(
      !is.na(Tcham) ~ Tcham,
      TRUE ~ 15.0  # Default if no weather data within 24 hours
    ),
    Pcham = case_when(
      !is.na(Pcham) ~ Pcham,
      TRUE ~ 101.325  # Default if no weather data within 24 hours
    )
  ) %>%
  # Select columns needed for goFlux auxfile plus original columns to preserve
  select(UniqueID, start.time, Area, Vtot, Tcham, Pcham, 
         Date, `Plot Tag`, `Plot Letter`, Chamber_ID, `File name`, Notes) %>%
  # Arrange by start time to maintain chronological order
  arrange(start.time)

# DATETIME VERIFICATION
cat("\n=== DATETIME VERIFICATION ===\n")
cat("Sample of start.time values (showing full datetime):\n")
print(head(auxfile$start.time, 10))
cat("\nClass of start.time column:", class(auxfile$start.time), "\n")
cat("Timezone:", attr(auxfile$start.time, "tzone"), "\n")

# Check if any times are missing (would show NA)
missing_times <- sum(is.na(auxfile$start.time))
cat("Number of measurements with missing start.time (these will be excluded):", missing_times, "\n")

# Display the first few rows to check the format
cat("\nFirst 6 rows of auxfile with weather data:\n")
print(auxfile)

# Check the structure
cat("\nStructure of auxfile:\n")
str(auxfile)

# SAVE WITH PROPER DATETIME FORMAT
# Save the auxfile ensuring datetime format is preserved
write.table(auxfile, "auxfile_goFlux_with_weather.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Also save as CSV for easier viewing
write_csv(auxfile, "auxfile_goFlux_with_weather.csv")

# ADDITIONAL: Create a human-readable version showing the datetime formatting
auxfile_readable <- auxfile %>%
  mutate(
    start.time_formatted = format(start.time, "%Y-%m-%d %H:%M:%S %Z")
  )

write_csv(auxfile_readable, "auxfile_goFlux_with_weather_formatted_datetime.csv")

# Print summary information
cat("\n=== SUMMARY ===\n")
cat("Auxfile created with", nrow(auxfile), "measurements\n")
cat("Date range:", as.character(range(auxfile$start.time, na.rm = TRUE)), "\n")
cat("Time range:", 
    format(min(auxfile$start.time, na.rm = TRUE), "%H:%M:%S"), "to",
    format(max(auxfile$start.time, na.rm = TRUE), "%H:%M:%S"), "\n")
cat("Temperature range:", round(range(auxfile$Tcham, na.rm = TRUE), 1), "°C\n")
cat("Files saved as:\n")
cat("  - auxfile_goFlux_with_weather.txt (main file for goFlux)\n")
cat("  - auxfile_goFlux_with_weather.csv (for viewing)\n")
cat("  - auxfile_goFlux_with_weather_formatted_datetime.csv (with readable datetime)\n")

# Display comprehensive quality checks
cat("\n=== DATA QUALITY CHECKS ===\n")
cat("Original flux rows:", nrow(flux_data), "\n")
cat("Final rows:", nrow(auxfile), "\n")
cat("Rows removed:", nrow(flux_data) - nrow(auxfile), "\n")
cat("Missing start.time:", sum(is.na(auxfile$start.time)), "\n")
cat("Missing Area:", sum(is.na(auxfile$Area)), "\n")
cat("Missing Vtot:", sum(is.na(auxfile$Vtot)), "\n")
cat("Missing Tcham:", sum(is.na(auxfile$Tcham)), "\n")
cat("Missing Pcham:", sum(is.na(auxfile$Pcham)), "\n")
cat("Duplicate UniqueIDs:", sum(duplicated(auxfile$UniqueID)), "\n")

# Show time distribution
cat("\n=== TIME DISTRIBUTION ===\n")
time_summary <- auxfile %>%
  mutate(hour = hour(start.time)) %>%
  count(hour, sort = TRUE)
cat("Measurements by hour of day:\n")
print(time_summary)

# Weather data usage summary
weather_matched <- sum(!is.na(auxfile$Tcham) & auxfile$Tcham != 15.0)
default_temp <- sum(auxfile$Tcham == 15.0)

cat("\n=== WEATHER DATA USAGE ===\n")
cat("Measurements with weather data:", weather_matched, "\n")
cat("Measurements using default temp (15°C):", default_temp, "\n")
cat("Weather data coverage:", round(100 * weather_matched / nrow(auxfile), 1), "%\n")

# Show which rows were removed and why (if any)
removed_rows <- flux_data %>%
  mutate(
    numeric_time = as.numeric(`Time of sampling`),
    time_int = floor(numeric_time),
    hours = floor(time_int / 100),
    minutes = time_int %% 100,
    formatted_time = case_when(
      is.na(numeric_time) ~ NA_character_,
      hours >= 0 & hours <= 23 & minutes >= 0 & minutes <= 59 ~ 
        sprintf("%02d:%02d:00", hours, minutes),
      TRUE ~ NA_character_  # Keep NA for invalid times
    ),
    datetime_string = case_when(
      !is.na(formatted_time) ~ paste(Date, formatted_time),
      TRUE ~ NA_character_
    ),
    start.time = as.POSIXct(datetime_string, 
                            format = "%Y-%m-%d %H:%M:%S", 
                            tz = "UTC"),
    Area = surface_area_cm2,
    tubing_volume_cm3 = pi * (1/16)^2 * 12 * 12 * 16.387,
    system_volume_cm3 = 70,
    Vtot = (volume_cm3 + tubing_volume_cm3 + system_volume_cm3) / 1000
  ) %>%
  filter(is.na(start.time) | is.na(Area) | is.na(Vtot)) %>%
  select(Date, `Plot Tag`, `Plot Letter`, `Time of sampling`, surface_area_cm2, volume_cm3)

if(nrow(removed_rows) > 0) {
  cat("\nRows that were removed due to missing critical data:\n")
  print(removed_rows)
} else {
  cat("\nNo rows were removed - all data successfully converted!\n")
}

cat("\n=== AUXFILE READY FOR goFlux! ===\n")
cat("Your auxfile contains all required columns with proper datetime formatting:\n")
cat("- UniqueID: Unique identifier for each measurement\n")
cat("- start.time: Full datetime in POSIXct format (YYYY-MM-DD HH:MM:SS UTC)\n")
cat("- Area: Chamber surface area in cm²\n")
cat("- Vtot: Total system volume in L (chamber + tubing + system)\n")
cat("- Tcham: Chamber temperature in °C (from weather data where available)\n")
cat("- Pcham: Chamber pressure in kPa\n")
cat("\nTime parsing handles both formats: 1357.0 and 1357 → 13:57:00\n")
cat("Weather data integrated from Yale-Myers weather station!\n")
cat("Original columns preserved: Date, Plot Tag, Plot Letter, Chamber_ID, File name, Notes\n")
cat("The start.time column preserves both date AND time information!\n")
