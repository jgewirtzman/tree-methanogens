# ==============================================================================
# Prepare Soil Flux Auxiliary File
# ==============================================================================
# Purpose: Prepares goFlux auxiliary metadata file for soil flux measurements.
#
# Pipeline stage: 01 Flux Processing (semirigid)
# Run after: None (independent of tree pipeline)
#
# Inputs:
#   - soilflux_total.csv (from iPad data)
#   - ymf_clean_sorted.csv (from data/raw/weather/)
#
# Outputs:
#   - auxfile_goFlux_soilflux_with_weather.csv
# ==============================================================================

library(dplyr)
library(readr)
library(lubridate)

# Read the CSV file with explicit column types to prevent auto-parsing
flux_data <- read_csv("../../../data/raw/lgr/semirigid_2020-2021/ipad_data/Cleaned data/soilflux_total.csv",
                      col_types = cols(
                        Date = col_character(),
                        Site = col_character(),
                        `Plot Tag` = col_character(),  # Force as character
                        `Plot letter` = col_character(),
                        `Person sampling` = col_character(),
                        `Time of sampling` = col_character(),  # Keep as character initially
                        `File name` = col_character(),
                        Notes = col_character()
                      ))

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
    # Convert pressure - using default for now (may need site-specific correction)
    Pcham_weather = 101.325  # Default atmospheric pressure in kPa
  ) %>%
  select(TIMESTAMP, Tcham_weather, Pcham_weather) %>%
  arrange(TIMESTAMP)

cat("Weather data range:", as.character(range(weather_clean$TIMESTAMP)), "\n")
cat("Temperature range:", round(range(weather_clean$Tcham_weather, na.rm = TRUE), 1), "°C\n")

cat("=== SOIL FLUX CHAMBER SPECIFICATIONS ===\n")
cat("Chamber surface area: 507.7 cm² (from 10.01\" interior diameter)\n")
cat("Chamber volume: 17.75 L (cap + collar headspace)\n")
cat("Tubing: 1/8\" ID × 12 ft\n")
cat("System volume: 70 cm³\n\n")

# Define chamber geometry constants
CHAMBER_SURFACE_AREA_CM2 <- 507.7  # Based on 25.43 cm interior diameter
CHAMBER_VOLUME_L <- 17.75  # Cap (14.42 L) + collar headspace (3.33 L)

# Calculate system volumes
tubing_volume_cm3 <- pi * (1/16)^2 * 12 * 12 * 16.387  # ≈ 18.6 cm³
system_volume_cm3 <- 70  # Given analyzer volume
total_system_volume_cm3 <- CHAMBER_VOLUME_L * 1000 + tubing_volume_cm3 + system_volume_cm3

cat("=== VOLUME CALCULATIONS ===\n")
cat("Chamber volume:", CHAMBER_VOLUME_L, "L\n")
cat("Tubing volume:", round(tubing_volume_cm3, 1), "cm³\n")
cat("System volume:", system_volume_cm3, "cm³\n")
cat("Total system volume:", round(total_system_volume_cm3, 1), "cm³ =", 
    round(total_system_volume_cm3/1000, 3), "L\n\n")

# Check for problematic data
cat("=== DATA QUALITY CHECKS ===\n")
cat("Total rows in dataset:", nrow(flux_data), "\n")

# Check for missing file names
missing_files <- sum(is.na(flux_data$`File name`))
cat("Missing file names:", missing_files, "\n")

# Check time format variations
cat("Time format examples:\n")
time_examples <- flux_data %>%
  filter(!is.na(`Time of sampling`)) %>%
  select(`Time of sampling`) %>%
  distinct() %>%
  head(15)
print(time_examples)

# Convert date format from M/D/YY to YYYY-MM-DD
flux_data <- flux_data %>%
  mutate(
    # Convert date format (handles M/D/YY format)
    Date_parsed = mdy(Date),  # This handles formats like 4/24/21
    Date_formatted = format(Date_parsed, "%Y-%m-%d")
  )

cat("\nDate conversion check:\n")
date_check <- flux_data %>%
  select(Date, Date_parsed, Date_formatted) %>%
  head(5)
print(date_check)

# Create auxfile for goFlux
cat("\n=== CREATING AUXFILE FOR GOFLUX ===\n")

auxfile <- flux_data %>%
  # Filter out rows with missing file names (can't process without gas data)
  filter(!is.na(`File name`)) %>%
  
  # Create unique identifier and process datetime
  mutate(
    # Handle time conversion from HHMM format
    numeric_time = as.numeric(`Time of sampling`),
    time_int = floor(numeric_time),
    hours = floor(time_int / 100),
    minutes = time_int %% 100,
    
    # Format time as HH:MM:SS
    formatted_time = case_when(
      hours >= 0 & hours <= 23 & minutes >= 0 & minutes <= 59 ~ 
        sprintf("%02d:%02d:00", hours, minutes),
      TRUE ~ NA_character_
    ),
    
    # Create datetime string and convert to POSIXct
    datetime_string = case_when(
      !is.na(formatted_time) ~ paste(Date_formatted, formatted_time),
      TRUE ~ NA_character_
    ),
    start.time = as.POSIXct(datetime_string, 
                            format = "%Y-%m-%d %H:%M:%S", 
                            tz = "UTC"),
    
    # Create unique ID with clear format
    # Format: MEASUREDATE_PLOT-LETTER_FILE (e.g., "20200615_2Jan-U_13")
    UniqueID = paste0(gsub("-", "", Date_formatted), "_", 
                      `Plot Tag`, "-", `Plot letter`, "_", 
                      `File name`),
    
    # Set chamber geometry (same for all measurements)
    Area = CHAMBER_SURFACE_AREA_CM2,  # cm²
    Vtot = total_system_volume_cm3 / 1000  # Convert to L
  ) %>%
  
  # Remove rows with missing critical data before merging with weather
  filter(!is.na(start.time), !is.na(Area), !is.na(Vtot)) %>%
  
  # Merge with weather data using nearest timestamp approach
  left_join(
    weather_clean,
    by = c("start.time" = "TIMESTAMP")
  ) %>%
  
  # Handle weather data matching
  mutate(
    Tcham = case_when(
      !is.na(Tcham_weather) ~ Tcham_weather,
      TRUE ~ NA_real_
    ),
    Pcham = case_when(
      !is.na(Pcham_weather) ~ Pcham_weather,
      TRUE ~ NA_real_
    )
  )

# For rows without exact weather matches, find nearest weather data
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

# Apply defaults for any remaining missing weather data
auxfile <- auxfile %>%
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
  
  # Select columns needed for goFlux
  select(UniqueID, start.time, Area, Vtot, Tcham, Pcham,
         Date = Date_formatted, Site, `Plot Tag`, `Plot letter`, 
         `Person sampling`, `File name`, Notes) %>%
  
  # Remove any rows with missing start.time (already filtered above, but safety check)
  filter(!is.na(start.time)) %>%
  
  # Arrange chronologically
  arrange(start.time)

# DATETIME AND WEATHER VERIFICATION
cat("\n=== DATETIME AND WEATHER VERIFICATION ===\n")
cat("Sample of start.time values:\n")
print(head(auxfile$start.time, 10))
cat("Class of start.time column:", class(auxfile$start.time), "\n")
cat("Timezone:", attr(auxfile$start.time, "tzone"), "\n")

# Weather data integration summary
weather_matched <- sum(!is.na(auxfile$Tcham) & auxfile$Tcham != 15.0)
default_temp <- sum(auxfile$Tcham == 15.0)

cat("\n=== WEATHER DATA INTEGRATION ===\n")
cat("Measurements with weather station data:", weather_matched, "\n")
cat("Measurements using default temp (15°C):", default_temp, "\n")
cat("Weather data coverage:", round(100 * weather_matched / nrow(auxfile), 1), "%\n")

# Display results
cat("Auxfile created successfully!\n")
cat("Final dataset contains", nrow(auxfile), "measurements\n")
cat("Date range:", as.character(range(auxfile$start.time)), "\n")

# Show first few rows
cat("\nFirst 6 rows of auxfile:\n")
print(head(auxfile))

# Check for any remaining issues
cat("\n=== FINAL QUALITY CHECKS ===\n")
cat("Missing start.time:", sum(is.na(auxfile$start.time)), "\n")
cat("Missing Area:", sum(is.na(auxfile$Area)), "\n")
cat("Missing Vtot:", sum(is.na(auxfile$Vtot)), "\n")
cat("Duplicate UniqueIDs:", sum(duplicated(auxfile$UniqueID)), "\n")

# Show measurement distribution by date
cat("\nMeasurements by date:\n")
date_summary <- auxfile %>%
  count(Date, sort = TRUE)
print(date_summary)

# Show time distribution
cat("\nMeasurements by hour:\n")
time_summary <- auxfile %>%
  mutate(hour = hour(start.time)) %>%
  count(hour, sort = TRUE)
print(time_summary)

# Save the auxfile with weather integration
write.table(auxfile, "auxfile_goFlux_soilflux_with_weather.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write_csv(auxfile, "auxfile_goFlux_soilflux_with_weather.csv")

# Create a version with readable datetime formatting
auxfile_readable <- auxfile %>%
  mutate(start.time_formatted = format(start.time, "%Y-%m-%d %H:%M:%S %Z"))

write_csv(auxfile_readable, "auxfile_goFlux_soilflux_with_weather_formatted.csv")

cat("\n=== FILES SAVED ===\n")
cat("Main auxfile: auxfile_goFlux_soilflux_with_weather.txt\n")
cat("CSV version: auxfile_goFlux_soilflux_with_weather.csv\n")
cat("Readable version: auxfile_goFlux_soilflux_with_weather_formatted.csv\n")

cat("\n=== CHAMBER SPECIFICATIONS APPLIED ===\n")
cat("Surface Area (Area):", CHAMBER_SURFACE_AREA_CM2, "cm²\n")
cat("Total Volume (Vtot):", round(auxfile$Vtot[1], 3), "L\n")
cat("  - Chamber volume:", CHAMBER_VOLUME_L, "L\n")
cat("  - Tubing volume:", round(tubing_volume_cm3, 1), "cm³\n")
cat("  - System volume:", system_volume_cm3, "cm³\n")

cat("\n=== READY FOR GOFLUX! ===\n")
cat("Your auxfile contains:\n")
cat("- UniqueID: Unique identifier for each measurement\n")
cat("- start.time: POSIXct datetime (YYYY-MM-DD HH:MM:SS UTC)\n")
cat("- Area: Chamber surface area (507.7 cm²)\n")
cat("- Vtot: Total system volume (~17.77 L)\n")
cat("- Tcham: Chamber temperature from Yale-Myers weather station (°C)\n")
cat("- Pcham: Chamber pressure from weather data or default (kPa)\n")
cat("- Original metadata preserved for reference\n")

# Show any excluded data
excluded_count <- nrow(flux_data) - nrow(auxfile)
if(excluded_count > 0) {
  cat("\nNOTE:", excluded_count, "rows excluded due to missing file names\n")
  excluded_data <- flux_data %>%
    filter(is.na(`File name`)) %>%
    select(Date, `Plot Tag`, `Plot letter`, `Time of sampling`)
  print(head(excluded_data))
}
