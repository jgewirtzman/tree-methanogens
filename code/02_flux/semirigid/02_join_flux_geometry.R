# ==============================================================================
# 02_join_flux_geometry.R -- chamber area and volume for each monthly stem deployment
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/).
#
# Each deployment's geometry comes from the chamber used that day (treeflux_total.xlsx
# Chamber) and the stem's diameter, by the rule in code/lib/stem_chamber_geometry.R:
# the chamber sheet's row for that stem x month x chamber where there is one, else the
# formula on the chamber's dimensions (chambers 7-9: mean of chambers 1-6). The source
# of each deployment's geometry is written with it (geometry_source, chamber_dims,
# D_stem_cm) and carried to data/tables/flux_stem.csv.
#
# Untagged stems (non-numeric Plot Tag) get no geometry here: the rescue
# (code/02_flux/rescue/01_untagged_auxfile.R) handles them with the same rule.
#
# Inputs:  data/raw/field_data/ipad_data/treeflux_total.xlsx (deployments, chamber used)
#          data/raw/flux_chamber_corrections/Chamber vol corr/Tree_fluxes_*.xlsx
#          data/raw/inventory/spatial_data/YM_trees_measured.csv.csv (D_stem of stems on
#          no chamber sheet)
# Outputs: data/processed/flux/flux_with_geometry_fixed.csv (tagged stems only)
#          outputs/audit/semirigid_chamber_dims.csv (the dimensions used per chamber)
# ==============================================================================
source("../../lib/project_library.R")   # the project package library (renv.lock), however this script is started
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(readr) })
ROOT <- "../../.."
source(file.path(ROOT, "code/lib/stem_chamber_geometry.R"))

flux_data <- read_excel(file.path(ROOT, "data/raw/field_data/ipad_data/treeflux_total.xlsx"), sheet = "Sheet1")
flux_clean <- flux_data %>%
  mutate(Sample_ID = as.character(`Plot Tag`), Chamber_ID = as.character(Chamber)) %>%
  filter(!is.na(Sample_ID) & !is.na(Chamber_ID))

rows <- semirigid_sheet_rows(ROOT)
dims <- semirigid_chamber_dims(rows)
dir.create(file.path(ROOT, "outputs/audit"), showWarnings = FALSE, recursive = TRUE)
write_csv(dims, file.path(ROOT, "outputs/audit/semirigid_chamber_dims.csv"))
cat("Chamber dimensions used where a deployment has no sheet row:\n"); print(dims, row.names = FALSE)

ym <- format(as.Date(flux_clean$Date), "%Y-%m")
ch <- sub("\\.0$", "", trimws(flux_clean$Chamber_ID))
month_chambers <- split(ch, ym)
tagged <- grepl("^[0-9]+(\\.0)?$", trimws(flux_clean$Sample_ID))
tl <- read.csv(file.path(ROOT, "data/raw/inventory/spatial_data/YM_trees_measured.csv.csv"), stringsAsFactors = FALSE)
smp <- sub("\\.0$", "", trimws(flux_clean$Sample_ID))
D_list <- suppressWarnings(as.numeric(tl$D_stem[match(smp, as.character(tl$Sample))]))

geo <- semirigid_deployment_geometry(smp, ym, ch, D_list, rows, dims, month_chambers)
geo[!tagged, ] <- NA
final_joined <- bind_cols(flux_clean, geo) %>%
  mutate(surface_area_cm2 = Sc_cm2, volume_cm3 = Vc_cm3, has_geometry = !is.na(surface_area_cm2)) %>%
  select(-Sc_cm2, -Vc_cm3)

cat("\nGeometry by source (tagged stems):\n")
print(as.data.frame(final_joined %>% filter(tagged) %>% count(geometry_source)), row.names = FALSE)
cat("\nBy month and source:\n")
print(table(ym[tagged], final_joined$geometry_source[tagged], useNA = "ifany"))
cat("\nTagged deployments with no geometry (stem on no chamber sheet and not in the tree list):\n")
print(as.data.frame(final_joined %>% filter(tagged, !has_geometry) %>% select(Date, Sample_ID, Chamber_ID)))
cat(sprintf("\nSheet rows whose chamber differs from the flux sheet's (geometry computed for the flux sheet's chamber): %d\n",
            sum(tagged & grepl("^computed", final_joined$geometry_source) & ym %in% SEMIRIGID_SHEET_MONTHS &
                  paste(smp, ym) %in% paste(rows$sample, rows$month))))
final_joined <- final_joined[tagged, ]          # untagged stems: code/02_flux/rescue/
write_csv(final_joined, file.path(ROOT, "data/processed/flux/flux_with_geometry_fixed.csv"))
cat("Wrote data/processed/flux/flux_with_geometry_fixed.csv:", nrow(final_joined), "rows\n")
