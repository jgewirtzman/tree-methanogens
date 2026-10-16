# ==============================================================================
# 01_load_and_prep_data.R -- random-forest workflow inputs from the canonical tables
# ------------------------------------------------------------------------------
# Reads data/tables/ only (built earlier in stage B by code/tables/) plus code/lib:
#   species, site, collars, inventory_stems, trees, flux_stem, flux_soil,
#   soil_readings, weather, moisture_survey, stream_points.
#
# Builds the PLACEHOLDER_* objects 02_rf_models.R reads:
#   TREE_JULY   2021 multi-height campaign + 2023 cross-species campaign
#   TREE_YEAR   monthly survey 2020-21 (semi-rigid chambers)
#   SOIL_YEAR   monthly soil collars 2020-21
#   MOISTURE_DEC_RASTER / MOISTURE_LOOKUP_XY   December 2020 moisture surface
#   DRIVERS     monthly air and soil temperature
#   TAXONOMY, PLOT_AREA, geodetic functions and constants
#
# Stem rows are flux_stem deployments (measurement_id), one row each. Trees are
# trees.csv's: tree_id, species_code and dbh_cm, and the position of the tree's
# census stem (or its GPS point), so one physical tree has one id, species, DBH and
# position across campaigns.
# Environmental covariates are the per-deployment values of flux_stem / flux_soil
# (air_temp_C, soil_temp_C, soil_vwc; their env_source column says where each came
# from). Run from code/05_model/ (pipeline workdir = script).
#
# "# KNOWN:" marks current behaviour that is a candidate for revision.
#
# Output: data/processed/integrated/rf_workflow_input_data_with_2023.RData
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(akima)
})

cat("=== RF WORKFLOW DATA IMPORT (data/tables) ===\n")
cat("Coordinate System: WGS84 (x=longitude, y=latitude)\n")
cat("Started at:", format(Sys.time()), "\n\n")

TAB <- "../../data/tables"
source("../lib/geometry.R")         # geodetic_transform(), geodetic_inverse()
source("../lib/species_levels.R")   # MODEL_SPECIES

# Every column as text. Numbers: the tables hold write_csv's round-trip decimal
# text, which a correctly rounded parser (jsonlite, i.e. C strtod) reads back to the
# double that was written; R's and readr's parsers misread some long decimals in the
# last digits.
rd <- function(t) read_csv(file.path(TAB, paste0(t, ".csv")), col_types = cols(.default = col_character()),
                           na = "", trim_ws = FALSE, progress = FALSE)
num <- function(x) {
  out <- rep(NA_real_, length(x))
  js <- !is.na(x) & grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x)
  if (any(js)) out[js] <- as.double(jsonlite::parse_json(paste0("[", paste(x[js], collapse = ","), "]"),
                                                         simplifyVector = TRUE))
  oth <- !is.na(x) & !js
  if (any(oth)) out[oth] <- as.numeric(x[oth])   # Inf, NaN
  out
}
wall <- function(x) as.POSIXct(substr(sub("T", " ", x), 1, 19), tz = "UTC")   # local wall-clock digits

SPECIES  <- rd("species")
SITE     <- rd("site")
COLLARS  <- rd("collars")
STEMS    <- rd("inventory_stems")
TREES    <- rd("trees")
FLUX     <- rd("flux_stem")
FLUX_S   <- rd("flux_soil")
SOILR    <- rd("soil_readings")
WEATHER  <- rd("weather")
MOIST    <- rd("moisture_survey")
STREAM   <- rd("stream_points")
stopifnot(nrow(SITE) == 1)
# KNOWN: CH4 flux (flux_stem, flux_soil) and the 2023 field-sheet covariates are
#   readr's parse of the processed flux files, as the model has always read them;
#   readr misreads some long decimals (up to ~1e-13 relative).

# =============================================================================
# 1. TAXONOMY
# =============================================================================
TAXONOMY <- data.frame(
  species_code = MODEL_SPECIES,
  species = SPECIES$species[match(MODEL_SPECIES, SPECIES$species_code)],
  genus   = SPECIES$genus[match(MODEL_SPECIES, SPECIES$species_code)],
  family  = SPECIES$family[match(MODEL_SPECIES, SPECIES$species_code)],
  order   = SPECIES$order[match(MODEL_SPECIES, SPECIES$species_code)],
  class   = SPECIES$class[match(MODEL_SPECIES, SPECIES$species_code)])
stopifnot(!anyNA(TAXONOMY))
species_mapping <- setNames(TAXONOMY$species, TAXONOMY$species_code)
cat("Taxonomy:", nrow(TAXONOMY), "species\n")

# =============================================================================
# 2. GEODETIC CONSTANTS (site.csv)
# =============================================================================
ref_lat <- num(SITE$ref_latitude)
ref_lon <- num(SITE$ref_longitude)
rotation_angle <- num(SITE$rotation_rad)
earth_radius <- num(SITE$earth_radius_m)
lat1_rad <- ref_lat * pi/180
lon1_rad <- ref_lon * pi/180
cat("Reference point:", ref_lat, "N,", ref_lon, "W; rotation", rotation_angle * 180/pi, "degrees\n")

# =============================================================================
# 3. MEASURED TREES (trees.csv): species, DBH, position
# =============================================================================
# Every flux row takes its tree's trees.csv attributes through flux_stem tree_id:
# species_code and the position of the census stem the tree is linked to
# (in_inventory: inventory_stems PX, PY through geodetic_transform), else its GPS
# point; trees with neither have no position.
TREE_ATTR <- local({
  st <- STEMS[match(TREES$in_inventory, STEMS$stem_id), ]
  ll <- geodetic_transform(num(st$PX), num(st$PY), lat1_rad, lon1_rad, rotation_angle, earth_radius)
  data.frame(tree_id = TREES$tree_id, species_code = TREES$species_code,
             species = unname(species_mapping[TREES$species_code]),
             x = coalesce(ll$longitude, num(TREES$longitude)),
             y = coalesce(ll$latitude, num(TREES$latitude)))
})
# the census positions are those trees.csv records (to its 7-decimal rounding)
stopifnot(!anyDuplicated(TREE_ATTR$tree_id),
          all(abs(TREE_ATTR$x - num(TREES$longitude)) < 1e-6, na.rm = TRUE),
          all(is.na(TREE_ATTR$x) == is.na(TREES$longitude)))
# DBH of a deployment: the diameter recorded at that measurement where there is one
# (2023: the field sheet, flux_stem dbh_recorded_cm; monthly survey: the survey's tree
# list; 2021: the 2021 survey consensus -- for those two, trees.dbh_cm is that record
# when its dbh_source says so), else the tree's trees.dbh_cm.
DBH_RECORD <- c(monthly_2020_2021 = "^monthly survey tree list", `2021_multiheight` = "^2021 survey DBH")
deployment_dbh_m <- function(flux) {
  i <- match(flux$tree_id, TREES$tree_id)
  own <- ifelse(flux$campaign %in% names(DBH_RECORD) &
                  mapply(grepl, coalesce(DBH_RECORD[flux$campaign], "^$"), coalesce(TREES$dbh_source[i], "")),
                num(TREES$dbh_cm[i]), NA_real_)
  rec <- coalesce(num(flux$dbh_recorded_cm), own)
  cat(sprintf("  DBH: %d of %d deployments at the diameter recorded with them; %d from trees.csv; %d none\n",
              sum(!is.na(rec)), nrow(flux), sum(is.na(rec) & !is.na(TREES$dbh_cm[i])), sum(is.na(rec) & is.na(TREES$dbh_cm[i]))))
  coalesce(rec, num(TREES$dbh_cm[i])) / 100
}
tree_attr <- function(tree_id) {
  a <- TREE_ATTR[match(tree_id, TREE_ATTR$tree_id), ]
  if (anyNA(a$species)) stop("trees without a model species: ",
                             paste(unique(tree_id[is.na(a$species)]), collapse = ", "), call. = FALSE)
  a
}
PLOT_AREA <- num(SITE$plot_side_m)^2

# =============================================================================
# 4. WEATHER (monthly air-temperature drivers)
# =============================================================================
# KNOWN: weather timestamps are the logger clock as recorded (weather.csv
#   timestamp_logger; about UTC-7 before the 2021-08-03 reset, UTC-4 after).
weather_clean <- data.frame(TIMESTAMP = as.POSIXct(sub("T", " ", WEATHER$timestamp_logger), tz = "UTC"),
                            air_temp_C = num(WEATHER$air_temp_avg_C)) %>%
  filter(!is.na(air_temp_C)) %>%
  arrange(TIMESTAMP)
cat("Weather:", nrow(weather_clean), "records,", format(min(weather_clean$TIMESTAMP), "%Y-%m-%d"),
    "to", format(max(weather_clean$TIMESTAMP), "%Y-%m-%d"), "\n\n")

# =============================================================================
# 5. 2021 MULTI-HEIGHT CAMPAIGN (rigid chambers; 50, 125, 200 cm; Kalmia 75 cm)
# =============================================================================
# One row per deployment, dated at its own logged start. Kalmia is read at 75 cm,
# else at the root crown: a root-crown deployment on a tree that also has a 75 cm
# one is left out.
# Soil temperature and VWC: flux_stem's (the tree mean of the 2021 survey soil
# readings under the flux file's spelling of the tree), else the mean of the tree's
# soil_readings by tree_id, which adds the readings filed under a second spelling
# of the same tree (WP6, WP8).
# KNOWN: 75 cm and root-crown deployments carry no measurement height
#   (02_rf_models.R reads a missing height as 125 cm).
sr21 <- SOILR %>% filter(!is.na(tree_id)) %>% group_by(tree_id) %>%
  summarise(soil_t = mean(num(soil_temp_C), na.rm = TRUE), soil_v = mean(num(soil_vwc_pct), na.rm = TRUE) / 100,
            .groups = "drop") %>%
  mutate(across(c(soil_t, soil_v), ~ ifelse(is.nan(.x), NA_real_, .x)))
h21 <- FLUX %>% filter(campaign == "2021_multiheight")
rc <- h21$height_label %in% "root crown"
stopifnot(all(rc | num(h21$height_cm) %in% c(50, 75, 125, 200)))
h21 <- h21[!(rc & h21$tree_id %in% h21$tree_id[num(h21$height_cm) %in% 75]), ]
a21 <- tree_attr(h21$tree_id)
TREE_JULY_2021 <- h21 %>%
  transmute(
    measurement_id, tree_id,
    species = a21$species, species_code = a21$species_code,
    Date = wall(start_time), month = month(Date), year = year(Date),
    stem_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = coalesce(num(soil_temp_C), sr21$soil_t[match(tree_id, sr21$tree_id)]),
    soil_moisture_abs = coalesce(num(soil_vwc), sr21$soil_v[match(tree_id, sr21$tree_id)]),
    dbh_m = deployment_dbh_m(h21), chamber_type = "rigid",
    measurement_height_cm = ifelse(num(height_cm) %in% c(50, 125, 200), num(height_cm), NA_real_),
    x = a21$x, y = a21$y) %>%
  filter(!is.na(stem_flux_umol_m2_s))
stopifnot(!anyNA(TREE_JULY_2021$Date))
cat(sprintf("2021 campaign: %d deployments on %d trees (50 cm %d, 125 cm %d, 200 cm %d, 75 cm or root crown %d)\n",
            nrow(TREE_JULY_2021), n_distinct(TREE_JULY_2021$tree_id),
            sum(TREE_JULY_2021$measurement_height_cm %in% 50), sum(TREE_JULY_2021$measurement_height_cm %in% 125),
            sum(TREE_JULY_2021$measurement_height_cm %in% 200), sum(is.na(TREE_JULY_2021$measurement_height_cm))))
cat(sprintf("  dates %s to %s, %d field days\n", format(min(TREE_JULY_2021$Date), "%Y-%m-%d"),
            format(max(TREE_JULY_2021$Date), "%Y-%m-%d"), length(unique(as.Date(TREE_JULY_2021$Date)))))

# =============================================================================
# 6. 2023 CROSS-SPECIES CAMPAIGN (rigid chambers, 125 cm)
# =============================================================================
y23 <- FLUX %>% filter(campaign == "2023_cross_species")
a23 <- tree_attr(y23$tree_id)
TREE_JULY_2023 <- y23 %>%
  transmute(
    measurement_id, tree_id,
    species = a23$species, species_code = a23$species_code,
    Date = as.POSIXct(date, tz = "UTC"), month = month(Date), year = year(Date),
    stem_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    # KNOWN: where the sheet has no air temperature, flux_stem air_temp_C is the
    #   tower reading nearest 00:00 of the date.
    air_temp_C = num(air_temp_C),
    soil_temp_C = num(soil_temp_C), soil_moisture_abs = num(soil_vwc),
    dbh_m = deployment_dbh_m(y23), chamber_type = "rigid", x = a23$x, y = a23$y) %>%
  filter(!is.na(stem_flux_umol_m2_s))
cat(sprintf("2023 campaign: %d deployments on %d trees (%d without a position)\n", nrow(TREE_JULY_2023),
            n_distinct(TREE_JULY_2023$tree_id), sum(is.na(TREE_JULY_2023$x))))

# =============================================================================
# 7. TREE_JULY = 2021 + 2023
# =============================================================================
TREE_JULY <- bind_rows(TREE_JULY_2021, TREE_JULY_2023) %>% arrange(Date)
n_2021 <- sum(TREE_JULY$year == 2021, na.rm = TRUE)
n_2023 <- sum(TREE_JULY$year == 2023, na.rm = TRUE)
cat(sprintf("TREE_JULY: %d rows from %d trees (2021 %d, 2023 %d); with coordinates %d\n\n", nrow(TREE_JULY),
            length(unique(TREE_JULY$tree_id)), n_2021, n_2023, sum(!is.na(TREE_JULY$x) & !is.na(TREE_JULY$y))))

# =============================================================================
# 8. TREE_YEAR: monthly survey 2020-21 (semi-rigid chambers, 125 cm)
# =============================================================================
# KNOWN: the 45 rescued untagged deployments are dated at 12:00 of their date (their
#   flux_stem air_temp_C is the tower reading nearest that noon), the others at the
#   logged start.
# KNOWN: soil temperature and VWC are flux_stem's: the collar readings averaged by
#   date and landscape position, nearest date within 3 days.
mon <- FLUX %>% filter(campaign == "monthly_2020_2021")
am <- tree_attr(mon$tree_id)
TREE_YEAR <- mon %>%
  transmute(
    measurement_id, tree_id,
    species = am$species, species_code = am$species_code,
    Date = if_else(untagged_rescue == "TRUE", as.POSIXct(paste(date, "12:00:00"), tz = "UTC"), wall(start_time)),
    month = month(Date), year = year(Date),
    stem_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C), soil_temp_C = num(soil_temp_C), soil_moisture_abs = num(soil_vwc),
    chamber_type = "semirigid", dbh_m = deployment_dbh_m(mon), x = am$x, y = am$y) %>%
  filter(!is.na(stem_flux_umol_m2_s))
stopifnot(!anyDuplicated(c(TREE_JULY$measurement_id, TREE_YEAR$measurement_id)))
cat(sprintf("TREE_YEAR: %d deployments from %d trees; soil moisture for %d; %d without a position\n\n",
            nrow(TREE_YEAR), length(unique(TREE_YEAR$tree_id)), sum(!is.na(TREE_YEAR$soil_moisture_abs)),
            sum(is.na(TREE_YEAR$x) | is.na(TREE_YEAR$y))))

# =============================================================================
# 9. SOIL_YEAR: monthly soil collars 2020-21
# =============================================================================
# Deployments that failed the closure screen (flux_soil qc_pass) and the two collars
# outside the stand (collars in_stand: WS_1-1, WS_1-2) are left out.
# KNOWN: air temperature is flux_soil's: the tower reading nearest 00:00 of the date
#   (else chamber temperature), not the deployment time.
# KNOWN: soil temperature and VWC are flux_soil's: the collar's own reading that
#   date, else the plot-tag mean that date (soilmoisture_total.csv).
in_stand_collar <- COLLARS$collar_id[COLLARS$in_stand == "TRUE"]
SOIL_YEAR <- FLUX_S %>%
  filter(qc_pass == "TRUE") %>%
  transmute(
    site_id = collar_id,
    Date = as.POSIXct(as.Date(date), tz = "UTC"),
    month = month(Date),
    year = year(Date),
    soil_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = num(soil_temp_C),
    soil_moisture_abs = num(soil_vwc)) %>%
  filter(!is.na(soil_flux_umol_m2_s))
n_before_stand <- nrow(SOIL_YEAR)
SOIL_YEAR <- SOIL_YEAR %>% filter(site_id %in% in_stand_collar)
SOIL_YEAR <- SOIL_YEAR %>%
  mutate(x = num(COLLARS$longitude[match(site_id, COLLARS$collar_id)]),
         y = num(COLLARS$latitude[match(site_id, COLLARS$collar_id)])) %>%
  select(site_id, Date, month, year, soil_flux_umol_m2_s, air_temp_C,
         soil_temp_C, soil_moisture_abs, x, y)
stopifnot(!anyNA(SOIL_YEAR$x), !anyNA(SOIL_YEAR$y))
cat(sprintf("SOIL_YEAR: %d deployments on %d collars (%d on out-of-stand collars dropped)\n\n",
            nrow(SOIL_YEAR), length(unique(SOIL_YEAR$site_id)), n_before_stand - nrow(SOIL_YEAR)))

# =============================================================================
# 10. DECEMBER 2020 MOISTURE RASTER
# =============================================================================
# Survey points plus the stream-channel flags at saturation (VWC = 1), interpolated
# (akima, linear) onto a 100 x 100 lon/lat grid.
combined_moisture <- rbind(
  data.frame(Longitude = num(MOIST$longitude), Latitude = num(MOIST$latitude), VWC = num(MOIST$soil_vwc_pct) / 100),
  data.frame(Longitude = num(STREAM$longitude), Latitude = num(STREAM$latitude), VWC = 1.0))
lon_range <- range(combined_moisture$Longitude)
lat_range <- range(combined_moisture$Latitude)
MOISTURE_DEC_RASTER <- interp(
  x = combined_moisture$Longitude,
  y = combined_moisture$Latitude,
  z = combined_moisture$VWC,
  xo = seq(lon_range[1], lon_range[2], length = 100),
  yo = seq(lat_range[1], lat_range[2], length = 100),
  duplicate = "mean"
)

# Lookup function for moisture at any x,y coordinate
moisture_lookup_xy <- function(x, y) {
  moisture_values <- numeric(length(x))
  for(i in 1:length(x)) {
    if(!is.na(x[i]) && !is.na(y[i])) {
      lon_idx <- which.min(abs(MOISTURE_DEC_RASTER$x - x[i]))
      lat_idx <- which.min(abs(MOISTURE_DEC_RASTER$y - y[i]))
      moisture_values[i] <- MOISTURE_DEC_RASTER$z[lon_idx, lat_idx]
    } else {
      moisture_values[i] <- NA
    }
  }
  return(moisture_values)
}
cat("December moisture raster:", length(MOISTURE_DEC_RASTER$x), "x", length(MOISTURE_DEC_RASTER$y), "\n")

# =============================================================================
# 11. MONTHLY DRIVERS
# =============================================================================
# KNOWN: air by calendar month of the logger-clock timestamp, all years pooled; soil
#   temperature by month from the training rows themselves (TREE_JULY, TREE_YEAR,
#   SOIL_YEAR).
air_monthly <- weather_clean %>%
  mutate(month = month(TIMESTAMP)) %>%
  group_by(month) %>%
  summarise(
    air_temp_C_mean = mean(air_temp_C, na.rm = TRUE),
    air_temp_C_sd = sd(air_temp_C, na.rm = TRUE),
    .groups = "drop"
  )
soil_monthly <- bind_rows(
  TREE_JULY %>% select(month, soil_temp_C),
  TREE_YEAR %>% select(month, soil_temp_C),
  SOIL_YEAR %>% select(month, soil_temp_C)
) %>%
  filter(!is.na(soil_temp_C)) %>%
  group_by(month) %>%
  summarise(
    soil_temp_C_mean = mean(soil_temp_C, na.rm = TRUE),
    soil_temp_C_sd = sd(soil_temp_C, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(
    soil_temp_C_mean = ifelse(is.nan(soil_temp_C_mean), NA_real_, soil_temp_C_mean),
    soil_temp_C_sd = ifelse(is.nan(soil_temp_C_sd), NA_real_, soil_temp_C_sd)
  )
DRIVERS <- tibble(month = 1:12) %>%
  left_join(air_monthly, by = "month") %>%
  left_join(soil_monthly, by = "month")
cat("Monthly drivers: air", sum(!is.na(DRIVERS$air_temp_C_mean)), "months, soil",
    sum(!is.na(DRIVERS$soil_temp_C_mean)), "months\n")

# =============================================================================
# 12. CHECKS AND SAVE
# =============================================================================
need <- list(TREE_JULY = c("measurement_id", "tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
             TREE_YEAR = c("measurement_id", "tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
             SOIL_YEAR = c("site_id", "Date", "soil_flux_umol_m2_s", "x", "y"))
for (n in names(need)) {
  miss <- setdiff(need[[n]], names(get(n)))
  if (nrow(get(n)) == 0 || length(miss)) stop(n, " is empty or lacks: ", paste(miss, collapse = ", "))
}

rf_workflow_data <- list(
  PLACEHOLDER_TREE_JULY = TREE_JULY,
  PLACEHOLDER_TREE_YEAR = TREE_YEAR,
  PLACEHOLDER_SOIL_YEAR = SOIL_YEAR,
  PLACEHOLDER_MOISTURE_DEC_RASTER = MOISTURE_DEC_RASTER,
  PLACEHOLDER_MOISTURE_LOOKUP_XY = moisture_lookup_xy,
  PLACEHOLDER_DRIVERS = DRIVERS,
  PLACEHOLDER_TAXONOMY = TAXONOMY,
  PLACEHOLDER_PLOT_AREA = PLOT_AREA,
  geodetic_transform = geodetic_transform,
  geodetic_inverse = geodetic_inverse,
  ref_lat = ref_lat,
  ref_lon = ref_lon,
  rotation_angle = rotation_angle,
  earth_radius = earth_radius,
  metadata = list(
    coordinate_system = "WGS84 geodetic",
    x_variable = "longitude (decimal degrees)",
    y_variable = "latitude (decimal degrees)",
    creation_date = Sys.time(),
    data_years = unique(c(TREE_JULY$year, TREE_YEAR$year, SOIL_YEAR$year)),
    n_species = length(unique(c(TREE_JULY$species, TREE_YEAR$species))),
    script_version = "2024-12 with 2023 data and geodetic coordinates"
  )
)
save(rf_workflow_data, file = "../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")
cat("\nSaved data/processed/integrated/rf_workflow_input_data_with_2023.RData\n")
cat("Completed at:", format(Sys.time()), "\n")
