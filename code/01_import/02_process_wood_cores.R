# ==============================================================================
# 02_process_wood_cores.R -- wood-core sections: density and moisture by radial position
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/01_import/).
#
# Each core (5.15 mm diameter) was cut into outer (sapwood), middle and inner
# (heartwood) sections. Per section: volume = pi r^2 x length; density = dry mass /
# volume (g cm-3); moisture on a fresh-mass and a dry-mass basis (%).
#
# Imputation: in seven listed trees (problematic_trees) a section whose dry mass exceeds
# its fresh mass, whose dry-basis moisture is > 300 %, or whose density is outside
# 0.1-2.0 g cm-3 has its moisture and density replaced by the species median of the
# remaining sections at that position. No other value is changed.
#
# Inputs:  data/raw/tree_cores/Tree Core Sectioning Data.csv
# Outputs (data/processed/tree_cores/):
#   tree_core_filled_complete.csv    every core: raw columns, original values, and the
#                                    imputed *_percent / *_density_final columns
#   tree_core_data_filled.csv        the same, with the imputed density under the
#                                    original names (*_density_g_cm3)
#   tree_core_filled_results.csv     the imputed values only, cores with any moisture
#   tree_core_filled_statistics.csv  mean and SD by species and section
# ==============================================================================

library(dplyr)
library(readr)

OUT <- "../../data/processed/tree_cores"
data <- read_csv("../../data/raw/tree_cores/Tree Core Sectioning Data.csv")

core_diameter_mm <- 5.15
core_radius_cm <- (core_diameter_mm / 2) / 10
calculate_volume <- function(length_cm) pi * core_radius_cm^2 * length_cm

processed_data <- data %>%
  rename(
    outer_length_cm = `Outer Section Length (cm)`,
    outer_fresh_mass_g = `OUTER SECTION FRESH MASS (g)`,
    outer_dry_mass_g = `Post-Freeze Dry OUTER SECTION Sample Mass (g)`,
    inner_length_cm = `Inner Section Length (cm)`,
    inner_fresh_mass_g = `INNER SECTION FRESH MASS (g)`,
    inner_dry_mass_g = `Post-Freeze Dry INNER SECTION Sample Mass (g)`,
    middle_length_mm = `Middle Section Length (mm)`,
    middle_fresh_mass_g = `MIDDLE SECTION FRESH MASS (g)`,
    middle_dry_mass_g = `Post-Freeze Dry MIDDLE SECTION Sample Mass (g)`,
    tree_id = `Tree Core ID`,
    species = `Tree Core Species Code`
  ) %>%
  mutate(
    outer_length_cm = as.numeric(outer_length_cm),
    inner_length_cm = as.numeric(inner_length_cm),
    middle_length_mm = as.numeric(middle_length_mm),
    outer_fresh_mass_g = as.numeric(outer_fresh_mass_g),
    outer_dry_mass_g = as.numeric(outer_dry_mass_g),
    inner_fresh_mass_g = as.numeric(inner_fresh_mass_g),
    inner_dry_mass_g = as.numeric(inner_dry_mass_g),
    middle_fresh_mass_g = as.numeric(middle_fresh_mass_g),
    middle_dry_mass_g = as.numeric(middle_dry_mass_g),
    middle_length_cm = middle_length_mm / 10
  )

complete_data <- processed_data %>%
  mutate(
    outer_volume_cm3 = ifelse(!is.na(outer_length_cm) & outer_length_cm > 0,
                              calculate_volume(outer_length_cm), NA),
    inner_volume_cm3 = ifelse(!is.na(inner_length_cm) & inner_length_cm > 0,
                              calculate_volume(inner_length_cm), NA),
    middle_volume_cm3 = ifelse(!is.na(middle_length_cm) & middle_length_cm > 0,
                               calculate_volume(middle_length_cm), NA),

    outer_moisture_fresh = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_fresh_mass_g > 0,
                                  ((outer_fresh_mass_g - outer_dry_mass_g) / outer_fresh_mass_g) * 100, NA),
    inner_moisture_fresh = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_fresh_mass_g > 0,
                                  ((inner_fresh_mass_g - inner_dry_mass_g) / inner_fresh_mass_g) * 100, NA),
    middle_moisture_fresh = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_fresh_mass_g > 0,
                                   ((middle_fresh_mass_g - middle_dry_mass_g) / middle_fresh_mass_g) * 100, NA),

    outer_moisture_dry = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_dry_mass_g > 0,
                                ((outer_fresh_mass_g - outer_dry_mass_g) / outer_dry_mass_g) * 100, NA),
    inner_moisture_dry = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_dry_mass_g > 0,
                                ((inner_fresh_mass_g - inner_dry_mass_g) / inner_dry_mass_g) * 100, NA),
    middle_moisture_dry = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_dry_mass_g > 0,
                                 ((middle_fresh_mass_g - middle_dry_mass_g) / middle_dry_mass_g) * 100, NA),

    outer_density_g_cm3 = ifelse(!is.na(outer_dry_mass_g) & !is.na(outer_volume_cm3) & outer_volume_cm3 > 0,
                                 outer_dry_mass_g / outer_volume_cm3, NA),
    inner_density_g_cm3 = ifelse(!is.na(inner_dry_mass_g) & !is.na(inner_volume_cm3) & inner_volume_cm3 > 0,
                                 inner_dry_mass_g / inner_volume_cm3, NA),
    middle_density_g_cm3 = ifelse(!is.na(middle_dry_mass_g) & !is.na(middle_volume_cm3) & middle_volume_cm3 > 0,
                                  middle_dry_mass_g / middle_volume_cm3, NA)
  )

# ---- imputation of the flagged sections (see header) ------------------------------
problematic_trees <- c("AB6'", "H104", "HX2", "RM1", "rms1", "WA14", "YB15")
SECTIONS <- c("outer", "inner", "middle")

is_problematic <- function(tree_id, fresh_mass, dry_mass, moisture_dry, density) {
  dry_exceeds_fresh  <- !is.na(fresh_mass) & !is.na(dry_mass) & dry_mass > fresh_mass
  extreme_moisture   <- !is.na(moisture_dry) & moisture_dry > 300
  impossible_density <- !is.na(density) & (density > 2.0 | density < 0.1)
  tree_id %in% problematic_trees & (dry_exceeds_fresh | extreme_moisture | impossible_density)
}
section_flag <- function(d, s)
  is_problematic(d$tree_id, d[[paste0(s, "_fresh_mass_g")]], d[[paste0(s, "_dry_mass_g")]],
                 d[[paste0(s, "_moisture_dry")]], d[[paste0(s, "_density_g_cm3")]])

# species medians from the unflagged sections
clean_data <- complete_data
for (s in SECTIONS) {
  bad <- section_flag(complete_data, s)
  for (v in c("moisture_fresh", "moisture_dry", "density_g_cm3"))
    clean_data[[paste0(s, "_", v, "_clean")]] <- ifelse(bad, NA, complete_data[[paste0(s, "_", v)]])
}
medians_by_species <- clean_data %>%
  group_by(species) %>%
  summarise(across(ends_with("_clean"), ~ median(.x, na.rm = TRUE),
                   .names = "{sub('_clean$', '_median', .col)}"), .groups = "drop")

final_data <- complete_data %>% left_join(medians_by_species, by = "species")
for (s in SECTIONS) {
  bad <- section_flag(final_data, s)
  p <- function(v) paste0(s, "_", v)
  final_data[[p("moisture_fresh_percent")]] <- ifelse(bad, final_data[[p("moisture_fresh_median")]], final_data[[p("moisture_fresh")]])
  final_data[[p("moisture_dry_percent")]]   <- ifelse(bad, final_data[[p("moisture_dry_median")]],   final_data[[p("moisture_dry")]])
  final_data[[p("density_final")]]          <- ifelse(bad, final_data[[p("density_g_cm3_median")]],  final_data[[p("density_g_cm3")]])
}
final_data <- final_data %>% select(-contains("median"), -contains("clean"))

cat("Cores:", nrow(final_data), "| trees checked for imputation:", paste(problematic_trees, collapse = ", "), "\n")
cat("\nImputed sections (original -> species median):\n")
for (s in SECTIONS) {
  bad <- section_flag(final_data, s)
  if (any(bad)) print(as.data.frame(final_data[bad, c("tree_id", "species", paste0(s, c(
    "_moisture_fresh", "_moisture_fresh_percent", "_moisture_dry", "_moisture_dry_percent",
    "_density_g_cm3", "_density_final")))]), row.names = FALSE)
}

# ---- outputs ------------------------------------------------------------------------
data_filled <- final_data %>%
  select(-outer_density_g_cm3, -inner_density_g_cm3, -middle_density_g_cm3) %>%
  rename(outer_density_g_cm3 = outer_density_final,
         inner_density_g_cm3 = inner_density_final,
         middle_density_g_cm3 = middle_density_final)
write_csv(data_filled, file.path(OUT, "tree_core_data_filled.csv"))

results_summary <- final_data %>%
  select(tree_id, species,
         outer_moisture_fresh_percent, outer_moisture_dry_percent, outer_density_final,
         inner_moisture_fresh_percent, inner_moisture_dry_percent, inner_density_final,
         middle_moisture_fresh_percent, middle_moisture_dry_percent, middle_density_final) %>%
  filter(!is.na(outer_moisture_fresh_percent) | !is.na(inner_moisture_fresh_percent) | !is.na(middle_moisture_fresh_percent))

section_stats <- function(s, label) {
  d <- final_data
  d$dens <- d[[paste0(s, "_density_final")]]
  d$mf <- d[[paste0(s, "_moisture_fresh_percent")]]
  d$md <- d[[paste0(s, "_moisture_dry_percent")]]
  d %>%
    filter(!is.na(dens)) %>%
    group_by(species) %>%
    summarise(
      section = label,
      n_samples = n(),
      mean_density = round(mean(dens, na.rm = TRUE), 3),
      sd_density = round(sd(dens, na.rm = TRUE), 3),
      mean_moisture_fresh = round(mean(mf, na.rm = TRUE), 1),
      sd_moisture_fresh = round(sd(mf, na.rm = TRUE), 1),
      mean_moisture_dry = round(mean(md, na.rm = TRUE), 1),
      sd_moisture_dry = round(sd(md, na.rm = TRUE), 1),
      .groups = "drop"
    )
}
all_stats <- bind_rows(section_stats("outer", "Outer (Sapwood)"),
                       section_stats("inner", "Inner (Heartwood)"),
                       section_stats("middle", "Middle")) %>%
  arrange(species, section)
cat("\nDensity and moisture by species and section:\n"); print(all_stats, n = Inf)

write_csv(results_summary, file.path(OUT, "tree_core_filled_results.csv"))
write_csv(final_data, file.path(OUT, "tree_core_filled_complete.csv"))
write_csv(all_stats, file.path(OUT, "tree_core_filled_statistics.csv"))

all_densities <- c(final_data$outer_density_final, final_data$inner_density_final, final_data$middle_density_final)
cat(sprintf("\nDensity range after imputation: %.3f-%.3f g cm-3\n",
            min(all_densities, na.rm = TRUE), max(all_densities, na.rm = TRUE)))
cat("Wrote", OUT, ": tree_core_filled_complete.csv, tree_core_data_filled.csv,",
    "tree_core_filled_results.csv, tree_core_filled_statistics.csv\n")
