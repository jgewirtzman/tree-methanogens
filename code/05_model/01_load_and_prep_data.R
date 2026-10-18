# ==============================================================================
# 01_load_and_prep_data.R -- random-forest workflow inputs from the canonical tables
# ------------------------------------------------------------------------------
# Reads data/tables/ only (built earlier in stage B by code/tables/) plus code/lib:
#   species, site, collars, inventory_stems, trees, flux_stem, flux_soil,
#   weather.
#
# Builds the list rf_workflow_data, read by 02_rf_models.R and the stage-C scripts:
#   TREE_JULY   2021 multi-height campaign + 2023 cross-species campaign
#   TREE_YEAR   monthly survey 2020-21 (semi-rigid chambers)
#   SOIL_YEAR   monthly soil collars 2020-21
#   DRIVERS     monthly air and soil temperature climatologies (stand predictions)
# Flux columns (stem_flux_nmol_m2_s, soil_flux_nmol_m2_s) are CH4 in nmol m-2 s-1.
#
# Stem rows are flux_stem deployments (measurement_id), one row each. Trees are
# trees.csv's: tree_id, species_code and dbh_cm, and the position of the tree's
# census stem (or its GPS point), so one physical tree has one id, species, DBH and
# position across campaigns.
# Conditions at the measurement: air temperature is flux_stem / flux_soil air_temp_C,
# the value measured at the deployment (carried, not a predictor of the forests; the
# with-air forests are the sensitivity 05_model/15_rf_air_temperature_test.R). Soil
# temperature and VWC are attached from data/tables/soil_readings.csv by
# code/lib/soil_fill.R attach_soil() -- the one soil rule -- with soil_match recording
# how each was found ("none": no reading within 3 days; the value stays NA). The
# monthly 24-h tower climatologies (section 10) drive the stand predictions only.
# Run from code/05_model/ (workdir = script).
#
# "# KNOWN:" marks current behaviour that is a candidate for revision.
#
# Output: data/processed/integrated/rf_workflow_input_data_with_2023.RData
# ==============================================================================

source("../lib/project_library.R")   # the project package library (renv.lock), however this script is started
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
})

cat("=== RF WORKFLOW DATA IMPORT (data/tables) ===\n")
cat("Coordinate System: WGS84 (x=longitude, y=latitude)\n")
cat("Started at:", format(Sys.time()), "\n\n")

TAB <- "../../data/tables"
source("../lib/geometry.R")         # geodetic_transform(), geodetic_inverse()
source("../lib/species_levels.R")   # MODEL_SPECIES
source("../lib/soil_fill.R")        # attach_soil(): soil temperature and VWC from soil_readings

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
WEATHER  <- rd("weather")
stopifnot(nrow(SITE) == 1)
# KNOWN: CH4 flux (flux_stem, flux_soil) and the 2023 field-sheet covariates are
#   readr's parse of the processed flux files, as the model has always read them;
#   readr misreads some long decimals (up to ~1e-13 relative).

# =============================================================================
# 1. MODEL SPECIES (species.csv)
# =============================================================================
species_mapping <- setNames(SPECIES$species[match(MODEL_SPECIES, SPECIES$species_code)], MODEL_SPECIES)
stopifnot(!anyNA(species_mapping))
cat("Model species:", length(species_mapping), "\n")

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
# list; 2021: the 2021 datasheets or tree DBH list -- for those two, trees.dbh_cm is that
# record when its dbh_source says so), else the tree's trees.dbh_cm. Measured DBH only:
# trees with no recorded DBH (six 2021 trees; trees.csv dbh_cm NA) keep dbh_m NA, and
# the random forest takes the missing value natively (ranger 0.17), so their
# deployments stay in training. dbh_source is carried per row.
DBH_RECORD <- c(monthly_2020_2021 = "^monthly survey tree list", `2021_multiheight` = "^2021 (flux datasheet|tree DBH list)")
deployment_dbh_m <- function(flux) {
  i <- match(flux$tree_id, TREES$tree_id)
  own <- ifelse(flux$campaign %in% names(DBH_RECORD) &
                  mapply(grepl, coalesce(DBH_RECORD[flux$campaign], "^$"), coalesce(TREES$dbh_source[i], "")),
                num(TREES$dbh_cm[i]), NA_real_)
  rec <- coalesce(num(flux$dbh_recorded_cm), own)
  cat(sprintf("  DBH: %d of %d deployments at the diameter recorded with them; %d from trees.csv; %d none\n",
              sum(!is.na(rec)), nrow(flux), sum(is.na(rec) & !is.na(TREES$dbh_cm[i])), sum(is.na(rec) & is.na(TREES$dbh_cm[i]))))
  out <- coalesce(rec, num(TREES$dbh_cm[i])) / 100
  # a missing DBH must be one trees.csv records as never measured
  stopifnot(all(TREES$dbh_source[i][is.na(out)] %in% "none: no DBH measurement recorded"))
  if (anyNA(out)) cat(sprintf("  DBH not measured (dbh_m NA, kept): %s\n", paste(unique(flux$tree_id[is.na(out)]), collapse = ", ")))
  out
}
deployment_dbh_source <- function(flux) {
  i <- match(flux$tree_id, TREES$tree_id)
  ifelse(!is.na(num(flux$dbh_recorded_cm)), "flux_stem dbh_recorded_cm (2023 field sheet)", paste("trees.csv:", TREES$dbh_source[i]))
}
# soil temperature and VWC (fraction) at a deployment: code/lib/soil_fill.R
soil_at <- function(flux, unit) {
  a <- attach_soil(flux[, c(if (unit == "tree") "tree_id" else "collar_id", "start_time")],
                   if (unit == "tree") "tree_id" else "collar_id", "start_time", unit)
  data.frame(soil_temp_C = a$soil_temp_C, soil_moisture_abs = a$soil_vwc_pct / 100, soil_match = a$soil_match)
}
tree_attr <- function(tree_id) {
  a <- TREE_ATTR[match(tree_id, TREE_ATTR$tree_id), ]
  if (anyNA(a$species)) stop("trees without a model species: ",
                             paste(unique(tree_id[is.na(a$species)]), collapse = ", "), call. = FALSE)
  a
}

# =============================================================================
# 4. TOWER RECORD (the monthly 24-h air-temperature climatology, section 10)
# =============================================================================
# Hourly records by local true time (weather.csv timestamp: the end of the hour; the
# local wall-clock digits, as wall()). Training rows take the air temperature measured
# at the deployment (flux_stem / flux_soil air_temp_C), not a value from this record.
weather_clean <- data.frame(TIMESTAMP = wall(WEATHER$timestamp),
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
# Soil temperature and VWC: attach_soil(unit = "tree") -- the tree's 2021 survey holes
# read within 3 days of the deployment.
# KNOWN: 75 cm and root-crown deployments carry no measurement height
#   (02_rf_models.R reads a missing height as 125 cm).
h21 <- FLUX %>% filter(campaign == "2021_multiheight")
rc <- h21$height_label %in% "root crown"
stopifnot(all(rc | num(h21$height_cm) %in% c(50, 75, 125, 200)))
h21 <- h21[!(rc & h21$tree_id %in% h21$tree_id[num(h21$height_cm) %in% 75]), ]
a21 <- tree_attr(h21$tree_id)
s21 <- soil_at(h21, "tree")
TREE_JULY_2021 <- h21 %>%
  transmute(
    measurement_id, tree_id, in_rf_training = in_rf_training == "TRUE",
    species = a21$species, species_code = a21$species_code,
    Date = wall(start_time), month = month(Date), year = year(Date),
    stem_flux_nmol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = s21$soil_temp_C, soil_moisture_abs = s21$soil_moisture_abs, soil_match = s21$soil_match,
    dbh_m = deployment_dbh_m(h21), dbh_source = deployment_dbh_source(h21), chamber_type = "rigid",
    measurement_height_cm = ifelse(num(height_cm) %in% c(50, 125, 200), num(height_cm), NA_real_),
    x = a21$x, y = a21$y) %>%
  filter(!is.na(stem_flux_nmol_m2_s))
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
s23 <- soil_at(y23, "tree")      # the 2023 field-sheet readings (soil_readings), same date
TREE_JULY_2023 <- y23 %>%
  transmute(
    measurement_id, tree_id, in_rf_training = in_rf_training == "TRUE",
    species = a23$species, species_code = a23$species_code,
    Date = as.POSIXct(date, tz = "UTC"), month = month(Date), year = year(Date),
    stem_flux_nmol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = s23$soil_temp_C, soil_moisture_abs = s23$soil_moisture_abs, soil_match = s23$soil_match,
    dbh_m = deployment_dbh_m(y23), dbh_source = deployment_dbh_source(y23), chamber_type = "rigid", x = a23$x, y = a23$y) %>%
  filter(!is.na(stem_flux_nmol_m2_s))
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
# KNOWN: the 45 rescued untagged deployments are dated at 12:00 of their date, the
#   others at the logged start.
# Soil temperature and VWC: attach_soil(unit = "tree") -- the collar readings averaged
#   over the tree's landscape position (trees.csv), that date or within 3 days.
mon <- FLUX %>% filter(campaign == "monthly_2020_2021")
am <- tree_attr(mon$tree_id)
smo <- soil_at(mon, "tree")
TREE_YEAR <- mon %>%
  transmute(
    measurement_id, tree_id, in_rf_training = in_rf_training == "TRUE",
    species = am$species, species_code = am$species_code,
    Date = if_else(untagged_rescue == "TRUE", as.POSIXct(paste(date, "12:00:00"), tz = "UTC"), wall(start_time)),
    month = month(Date), year = year(Date),
    stem_flux_nmol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = smo$soil_temp_C, soil_moisture_abs = smo$soil_moisture_abs, soil_match = smo$soil_match,
    chamber_type = "semirigid", dbh_m = deployment_dbh_m(mon), dbh_source = deployment_dbh_source(mon), x = am$x, y = am$y) %>%
  filter(!is.na(stem_flux_nmol_m2_s))
stopifnot(!anyDuplicated(c(TREE_JULY$measurement_id, TREE_YEAR$measurement_id)))
cat(sprintf("TREE_YEAR: %d deployments from %d trees; soil moisture for %d; %d without a position\n\n",
            nrow(TREE_YEAR), length(unique(TREE_YEAR$tree_id)), sum(!is.na(TREE_YEAR$soil_moisture_abs)),
            sum(is.na(TREE_YEAR$x) | is.na(TREE_YEAR$y))))

# =============================================================================
# 9. SOIL_YEAR: monthly soil collars 2020-21
# =============================================================================
# Deployments that failed the closure screen (flux_soil qc_pass) and the two collars
# outside the stand (collars in_stand: WS_1-1, WS_1-2) are left out.
# Air temperature: flux_soil air_temp_C, the tower hour containing the logged start.
# Soil temperature and VWC: attach_soil(unit = "collar") (soil_match says which step).
in_stand_collar <- COLLARS$collar_id[COLLARS$in_stand == "TRUE"]
FLUX_S <- FLUX_S %>% filter(qc_pass == "TRUE")
ssy <- soil_at(FLUX_S, "collar")
SOIL_YEAR <- FLUX_S %>%
  transmute(
    site_id = collar_id,
    Date = as.POSIXct(as.Date(date), tz = "UTC"),
    month = month(Date),
    year = year(Date),
    soil_flux_nmol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = ssy$soil_temp_C,
    soil_moisture_abs = ssy$soil_moisture_abs, soil_match = ssy$soil_match) %>%
  filter(!is.na(soil_flux_nmol_m2_s))
n_before_stand <- nrow(SOIL_YEAR)
SOIL_YEAR <- SOIL_YEAR %>% filter(site_id %in% in_stand_collar)
SOIL_YEAR <- SOIL_YEAR %>%
  mutate(x = num(COLLARS$longitude[match(site_id, COLLARS$collar_id)]),
         y = num(COLLARS$latitude[match(site_id, COLLARS$collar_id)])) %>%
  select(site_id, Date, month, year, soil_flux_nmol_m2_s, air_temp_C,
         soil_temp_C, soil_moisture_abs, soil_match, x, y)
stopifnot(!anyNA(SOIL_YEAR$x), !anyNA(SOIL_YEAR$y))
cat(sprintf("SOIL_YEAR: %d deployments on %d collars (%d on out-of-stand collars dropped)\n\n",
            nrow(SOIL_YEAR), length(unique(SOIL_YEAR$site_id)), n_before_stand - nrow(SOIL_YEAR)))

# =============================================================================
# 10. MONTHLY DRIVERS
# =============================================================================
# The stand-prediction drivers (not training inputs): 24-h climatologies over the
# climatology years (CLIM_YEARS; the years of the moisture climatology,
# 04_drivers/02_moisture_climatology.R).
# Air: the mean of the tower's hourly air temperature by calendar month of the local
#   true time, CLIM_YEARS.
# Soil temperature: the monthly climatology of code/04_drivers/00_soil_temp_climatology.R
#   (soil_readings date means regressed on trailing tower air, applied to CLIM_YEARS).
CLIM_YEARS <- 2019:2022
air_monthly <- weather_clean %>%
  filter(year(TIMESTAMP) %in% CLIM_YEARS) %>%
  mutate(month = month(TIMESTAMP)) %>%
  group_by(month) %>%
  summarise(
    air_temp_C_mean = mean(air_temp_C),
    air_temp_C_sd = sd(air_temp_C),
    .groups = "drop"
  )
STCLIM <- read.csv("../../outputs/tables/soil_temp_climatology_monthly.csv")
soil_monthly <- data.frame(month = STCLIM$mo, soil_temp_C_mean = STCLIM$soil_temp_C)
stopifnot(setequal(soil_monthly$month, 1:12), !anyNA(soil_monthly$soil_temp_C_mean))
DRIVERS <- tibble(month = 1:12) %>%
  left_join(air_monthly, by = "month") %>%
  left_join(soil_monthly, by = "month")
cat("Monthly drivers: air", sum(!is.na(DRIVERS$air_temp_C_mean)), "months, soil",
    sum(!is.na(DRIVERS$soil_temp_C_mean)), "months (soil_temp_climatology_monthly.csv)\n")
cat("Soil readings attached (soil_match):\n")
print(table(table = rep(c("TREE_JULY", "TREE_YEAR", "SOIL_YEAR"), c(nrow(TREE_JULY), nrow(TREE_YEAR), nrow(SOIL_YEAR))),
            soil_match = c(TREE_JULY$soil_match, TREE_YEAR$soil_match, SOIL_YEAR$soil_match)))

# =============================================================================
# 11. CHECKS AND SAVE
# =============================================================================
need <- list(TREE_JULY = c("measurement_id", "tree_id", "species", "Date", "stem_flux_nmol_m2_s", "x", "y"),
             TREE_YEAR = c("measurement_id", "tree_id", "species", "Date", "stem_flux_nmol_m2_s", "x", "y"),
             SOIL_YEAR = c("site_id", "Date", "soil_flux_nmol_m2_s", "x", "y"))
for (n in names(need)) {
  miss <- setdiff(need[[n]], names(get(n)))
  if (nrow(get(n)) == 0 || length(miss)) stop(n, " is empty or lacks: ", paste(miss, collapse = ", "))
}

rf_workflow_data <- list(
  TREE_JULY = TREE_JULY,
  TREE_YEAR = TREE_YEAR,
  SOIL_YEAR = SOIL_YEAR,
  DRIVERS = DRIVERS
)
save(rf_workflow_data, file = "../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")
cat("\nSaved data/processed/integrated/rf_workflow_input_data_with_2023.RData\n")
cat("Completed at:", format(Sys.time()), "\n")
