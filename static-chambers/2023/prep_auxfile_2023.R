# Prep auxfile for YMF 2023 data - Simple version with correct column names
library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(purrr)

# Function to create goFlux-compatible auxfile from YMF 2023 data
create_ymf2023_auxfile <- function(data_path, 
                                   surface_area_file,
                                   additional_vol_file,
                                   simplified_vol_file = NULL) {
  
  cat("=== Loading YMF 2023 Data ===\n")
  
  # Read the data
  ymf23 <- read_csv(data_path, show_col_types = FALSE)
  
  cat("Loaded", nrow(ymf23), "measurements\n")
  
  # Load dimension files
  cat("\n=== Loading Dimension Files ===\n")
  surface_area <- read_csv(surface_area_file, show_col_types = FALSE)
  additional_vol <- read_csv(additional_vol_file, show_col_types = FALSE)
  
  if (!is.null(simplified_vol_file) && file.exists(simplified_vol_file)) {
    simplified_vol <- read_csv(simplified_vol_file, show_col_types = FALSE)
  } else {
    simplified_vol <- NULL
  }
  
  # Process surface area data
  cat("\n=== Processing Surface Area Data ===\n")
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
  
  print(surface_area_clean)
  
  # Process analyzer volumes
  cat("\n=== Processing System Volume Data ===\n")
  analyzer_volumes <- additional_vol %>%
    mutate(
      analyzer_name = case_when(
        tolower(instrument) %in% c("lgr_mgga", "lgr", "lgr3") ~ "LGR",
        tolower(instrument) == "picarro" ~ "Picarro",
        TRUE ~ toupper(instrument)
      ),
      analyzer_cell_volume_cm3 = as.numeric(analyzer_cell),
      tubing_volume_cm3 = as.numeric(tubing),
      filter_volume_cm3 = 0  # No drierite used
    ) %>%
    select(analyzer_name, analyzer_cell_volume_cm3, tubing_volume_cm3, filter_volume_cm3)
  
  print(analyzer_volumes)
  
  # Process chamber volumes if available
  cat("\n=== Processing Chamber Volume Data ===\n")
  if (!is.null(simplified_vol)) {
    names(simplified_vol) <- gsub("\\n|\\r", "_", names(simplified_vol))
    names(simplified_vol) <- gsub("\\s+", "_", names(simplified_vol))
    
    chamber_volumes <- simplified_vol %>%
      filter(!is.na(`Chamber_Alt_ID`) & !is.na(`Total_Volume_(mL)`)) %>%
      mutate(
        chamber_series = `Chamber_Alt_ID`,
        chamber_volume_cm3 = as.numeric(`Total_Volume_(mL)`)
      ) %>%
      filter(!is.na(chamber_volume_cm3) & chamber_series %in% c("A", "B", "C", "D")) %>%
      group_by(chamber_series) %>%
      summarise(
        avg_chamber_volume_cm3 = mean(chamber_volume_cm3, na.rm = TRUE),
        .groups = 'drop'
      )
    
    print(chamber_volumes)
  } else {
    chamber_volumes <- NULL
    cat("No chamber volume data - using estimates\n")
  }
  
  # Create chamber mapping for 2023 data
  cat("\n=== Processing Chamber Mapping ===\n")
  
  # Get unique chamber IDs from the data
  unique_chambers <- unique(ymf23$`Chamber ID`[!is.na(ymf23$`Chamber ID`)])
  cat("Chamber IDs found:", paste(unique_chambers, collapse = ", "), "\n")
  
  # Parse chamber IDs (e.g., D8, C3, C4)
  chamber_mapping <- data.frame(
    chamber_id_original = unique_chambers,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      # Extract letter from ID, or assign D to "6"
      chamber_series = case_when(
        chamber_id_original == "6" ~ "D",  # Chamber 6 is D series
        TRUE ~ gsub("[0-9]", "", chamber_id_original)  # Extract letter from others
      )
    )
  
  # Join with surface areas
  chamber_dimensions <- chamber_mapping %>%
    left_join(surface_area_clean, by = "chamber_series")
  
  # Join with volumes
  if (!is.null(chamber_volumes)) {
    chamber_dimensions <- chamber_dimensions %>%
      left_join(chamber_volumes, by = "chamber_series") %>%
      rename(chamber_volume_cm3 = avg_chamber_volume_cm3)
  } else {
    # Use estimated volumes based on relative surface areas
    b_surface_area <- surface_area_clean$surface_area_cm2[surface_area_clean$chamber_series == "B"]
    estimated_b_volume <- 590
    
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
  
  print(chamber_dimensions)
  
  # Process the 2023 data
  cat("\n=== Processing YMF 2023 Data ===\n")
  
  ymf23_processed <- ymf23 %>%
    mutate(
      # Parse date (M/D/YY format)
      date_parts = str_split(Date, "/"),
      month = as.numeric(map_chr(date_parts, 1)),
      day = as.numeric(map_chr(date_parts, 2)),
      year = as.numeric(map_chr(date_parts, 3)) + 2000,
      date_clean = as.Date(paste(year, sprintf("%02d", month), sprintf("%02d", day), sep = "-")),
      
      # Clean start time
      start_time_raw = as.character(`Flux Start (System) Time`),
      start_time_clean = case_when(
        grepl("PM|AM", start_time_raw) ~ gsub(" PM| AM", "", start_time_raw),
        TRUE ~ start_time_raw
      ),
      
      # Create datetime
      datetime_start = case_when(
        !is.na(start_time_clean) & start_time_clean != "<NA>" & start_time_clean != "" & !is.na(date_clean) ~ {
          as.POSIXct(paste(as.character(date_clean), start_time_clean), 
                     format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
        },
        TRUE ~ as.POSIXct(NA)
      ),
      
      # Clean end time
      end_time_raw = as.character(`Flux End (System) Time`),
      end_time_clean = case_when(
        !is.na(end_time_raw) & end_time_raw != "<NA>" ~ {
          case_when(
            end_time_raw == "21:31:15" ~ "12:35:15",  # Fix typo
            TRUE ~ gsub(" PM| AM", "", end_time_raw)
          )
        },
        TRUE ~ NA_character_
      ),
      
      # Create end datetime
      datetime_end = case_when(
        !is.na(end_time_clean) & !is.na(date_clean) ~ {
          as.POSIXct(paste(as.character(date_clean), end_time_clean), 
                     format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
        },
        TRUE ~ as.POSIXct(NA)
      ),
      
      # Calculate observation length
      obs_length_seconds = as.numeric(difftime(datetime_end, datetime_start, units = "secs")),
      
      # Create UniqueID
      UniqueID = paste(
        format(date_clean, "%Y%m%d"),
        `Tree Tag`,
        gsub(":", "", start_time_clean),
        sep = "_"
      ),
      
      # Temperature conversions
      air_temp_F = as.numeric(`Air Temp (°F)`),
      air_temp_C = (air_temp_F - 32) * 5/9,
      stem_temp_C = as.numeric(`Stem Temp (°C)`),
      
      # Calculate mean VWC from the three VWC columns
      # Handle "Saturated" as probe maximum (~60% for most probes)
      vwc_1_raw = `VWC 1 (%)`,
      vwc_2_raw = `VWC 2 (%)`,
      vwc_3_raw = `VWC 3 (%)`,
      
      vwc_1 = case_when(
        vwc_1_raw == "Saturated" ~ 60,  # Typical probe maximum
        TRUE ~ as.numeric(vwc_1_raw)
      ),
      vwc_2 = case_when(
        vwc_2_raw == "Saturated" ~ 60,
        TRUE ~ as.numeric(vwc_2_raw)
      ),
      vwc_3 = case_when(
        vwc_3_raw == "Saturated" ~ 60,
        TRUE ~ as.numeric(vwc_3_raw)
      ),
      vwc_mean = rowMeans(cbind(vwc_1, vwc_2, vwc_3), na.rm = TRUE),
      
      # Clean analyzer name
      analyzer_clean = "LGR"  # All 2023 data uses LGR3
    ) %>%
    # Join with chamber dimensions
    left_join(chamber_dimensions, 
              by = c("Chamber ID" = "chamber_id_original")) %>%
    # Join with system volumes
    left_join(analyzer_volumes, by = c("analyzer_clean" = "analyzer_name"))
  
  # Calculate total system volume
  ymf23_processed <- ymf23_processed %>%
    mutate(
      total_system_volume_cm3 = case_when(
        !is.na(chamber_volume_cm3) & !is.na(tubing_volume_cm3) & 
          !is.na(analyzer_cell_volume_cm3) & !is.na(filter_volume_cm3) ~
          chamber_volume_cm3 + tubing_volume_cm3 + analyzer_cell_volume_cm3 + filter_volume_cm3,
        TRUE ~ NA_real_
      ),
      total_system_volume_L = total_system_volume_cm3 / 1000
    )
  
  # Create final auxfile with all columns
  cat("\n=== Creating goFlux Auxfile ===\n")
  
  auxfile <- ymf23_processed %>%
    mutate(
      # Core goFlux columns
      start.time = datetime_start,
      Area = surface_area_cm2,
      Vtot = total_system_volume_L,
      Tcham = case_when(
        !is.na(air_temp_C) ~ air_temp_C,
        !is.na(stem_temp_C) ~ stem_temp_C,
        TRUE ~ 20
      ),
      Pcham = 101.325,
      obs.length = case_when(
        !is.na(obs_length_seconds) & obs_length_seconds > 0 ~ obs_length_seconds,
        TRUE ~ 600
      )
    ) %>%
    filter(!is.na(datetime_start))  # Keep only valid measurements
  
  # Validation
  cat("\nValidation:\n")
  cat("  Total measurements:", nrow(auxfile), "\n")
  
  complete_required <- auxfile %>%
    filter(!is.na(UniqueID) & !is.na(start.time) & !is.na(Area) & 
             !is.na(Vtot) & !is.na(Tcham) & !is.na(Pcham))
  cat("  Complete measurements:", nrow(complete_required), "\n")
  
  # Check for anomalies
  suspicious_lengths <- auxfile %>%
    filter(!is.na(obs.length) & (obs.length < 180 | obs.length > 900))
  
  if (nrow(suspicious_lengths) > 0) {
    cat("\nWARNING: Found", nrow(suspicious_lengths), "measurements with unusual durations:\n")
    print(suspicious_lengths %>% 
            select(UniqueID, `Tree Tag`, start.time, obs.length) %>%
            head(5))
  }
  
  # Summary statistics
  cat("\nSummary:\n")
  cat("  Temperature range:", round(range(auxfile$Tcham, na.rm = TRUE), 1), "°C\n")
  cat("  Area range:", round(range(auxfile$Area, na.rm = TRUE), 1), "cm²\n")
  cat("  Volume range:", round(range(auxfile$Vtot, na.rm = TRUE), 3), "L\n")
  cat("  Mean VWC range:", round(range(auxfile$vwc_mean, na.rm = TRUE), 1), "%\n")
  
  # Chamber distribution
  cat("\nChamber distribution:\n")
  print(table(auxfile$`Chamber ID`, useNA = "ifany"))
  
  # Species distribution
  cat("\nSpecies distribution:\n")
  print(table(auxfile$`Species Code`, useNA = "ifany"))
  
  return(auxfile)
}

# Execute the function
cat("Creating auxfile for YMF 2023 data...\n\n")

# Update these paths to your file locations
base_path <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers"

# Create the auxfile
auxfile_2023 <- create_ymf2023_auxfile(
  data_path = file.path(base_path, "2023/Compiled YMF Data 2023 - Sheet1.csv"),
  surface_area_file = file.path(base_path, "dims/surface_area.csv"),
  additional_vol_file = file.path(base_path, "dims/additional_vol.csv"),
  simplified_vol_file = file.path(base_path, "dims/simplified_volume.csv")
)

# Save the auxfile
write_csv(auxfile_2023, "ymf2023_goflux_auxfile.csv")
cat("\n=== Auxfile saved to: ymf2023_goflux_auxfile.csv ===\n")

# Display sample
cat("\nSample of auxfile (first 3 rows, key columns):\n")
print(auxfile_2023 %>% 
        select(UniqueID, start.time, Area, Vtot, Tcham, obs.length, `Chamber ID`, `Tree Tag`) %>%
        head(3))