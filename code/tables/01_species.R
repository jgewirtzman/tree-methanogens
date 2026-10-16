#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 01_species.R -- data/tables/species.csv: one row per species code used anywhere
# ------------------------------------------------------------------------------
# Codes come from the stand census, the four flux campaigns, internal gas, ddPCR
# and the trait table. Binomials are those of the species_mapping vectors in the
# analysis scripts (all agree code by code; 08_dbh_by_species_campaign.R reports
# BEPO and CALA at genus level as uncertain field identifications) and of
# code/lib/trait_species.csv; family/order/class are those the model's TAXONOMY used,
# extended to the census-only genera. The model takes its TAXONOMY from this table
# (codes MODEL_SPECIES, code/lib/species_levels.R). "TSLA" (2021 aux file)
# is a typo for TSCA and is not a species (code/lib/tree_ids.R).
#
# Checks: every code in the sources and in MODEL_SPECIES has a row; binomial agrees
# with trait_species.csv; gymnosperm agrees with geometry.R.
# Writes data/tables/species.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/tree_ids.R"); source("code/lib/geometry.R"); source("code/lib/species_levels.R")
OUT <- "data/tables"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SP <- tribble(
  ~species_code, ~species,                  ~common_name,              ~family,          ~order,         ~class,
  "ACPE",        "Acer pensylvanicum",      "striped maple",           "Sapindaceae",    "Sapindales",   "Magnoliopsida",
  "ACRU",        "Acer rubrum",             "red maple",               "Sapindaceae",    "Sapindales",   "Magnoliopsida",
  "ACSA",        "Acer saccharum",          "sugar maple",             "Sapindaceae",    "Sapindales",   "Magnoliopsida",
  "AMSP",        "Amelanchier sp.",         "serviceberry",            "Rosaceae",       "Rosales",      "Magnoliopsida",
  "BEAL",        "Betula alleghaniensis",   "yellow birch",            "Betulaceae",     "Fagales",      "Magnoliopsida",
  "BELE",        "Betula lenta",            "black birch",             "Betulaceae",     "Fagales",      "Magnoliopsida",
  "BEPA",        "Betula papyrifera",       "paper birch",             "Betulaceae",     "Fagales",      "Magnoliopsida",
  "BEPO",        "Betula populifolia",      "gray birch",              "Betulaceae",     "Fagales",      "Magnoliopsida",
  "CALA",        "Carya laciniosa",         "shellbark hickory",       "Juglandaceae",   "Fagales",      "Magnoliopsida",
  "CAOV",        "Carya ovata",             "shagbark hickory",        "Juglandaceae",   "Fagales",      "Magnoliopsida",
  "FAGR",        "Fagus grandifolia",       "American beech",          "Fagaceae",       "Fagales",      "Magnoliopsida",
  "FRAM",        "Fraxinus americana",      "white ash",               "Oleaceae",       "Lamiales",     "Magnoliopsida",
  "HAVI",        "Hamamelis virginiana",    "American witch-hazel",    "Hamamelidaceae", "Saxifragales", "Magnoliopsida",
  "KALA",        "Kalmia latifolia",        "mountain laurel",         "Ericaceae",      "Ericales",     "Magnoliopsida",
  "LITU",        "Liriodendron tulipifera", "tuliptree",               "Magnoliaceae",   "Magnoliales",  "Magnoliopsida",
  "PIST",        "Pinus strobus",           "eastern white pine",      "Pinaceae",       "Pinales",      "Pinopsida",
  "PRSE",        "Prunus serotina",         "black cherry",            "Rosaceae",       "Rosales",      "Magnoliopsida",
  "QUAL",        "Quercus alba",            "white oak",               "Fagaceae",       "Fagales",      "Magnoliopsida",
  "QUCO",        "Quercus coccinea",        "scarlet oak",             "Fagaceae",       "Fagales",      "Magnoliopsida",
  "QURU",        "Quercus rubra",           "northern red oak",        "Fagaceae",       "Fagales",      "Magnoliopsida",
  "QUVE",        "Quercus velutina",        "black oak",               "Fagaceae",       "Fagales",      "Magnoliopsida",
  "SAAL",        "Sassafras albidum",       "sassafras",               "Lauraceae",      "Laurales",     "Magnoliopsida",
  "TSCA",        "Tsuga canadensis",        "eastern hemlock",         "Pinaceae",       "Pinales",      "Pinopsida",
  "VACO",        "Vaccinium corymbosum",    "highbush blueberry",      "Ericaceae",      "Ericales",     "Magnoliopsida") %>%
  mutate(genus = sub(" .*", "", species), gymnosperm = class == "Pinopsida") %>%
  select(species_code, species, common_name, genus, family, order, class, gymnosperm)
stopifnot(!anyDuplicated(SP$species_code), !anyDuplicated(SP$species))

# ---- every code in use has a row ---------------------------------------------------
codes <- function(x) unique(na.omit(fix_species_code(x)))
rd <- function(f, ...) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, ...)
used <- list(
  census_2019   = codes(rd("data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv", fileEncoding = "UTF-8-BOM")$Species_Code),
  census_bytag  = codes(rd("data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_bytag.csv", fileEncoding = "UTF-8-BOM")$`Species (4 char)`),
  monthly       = codes(rd("data/raw/inventory/spatial_data/YM_trees_measured.csv.csv", fileEncoding = "UTF-8-BOM")$Species),
  flux_2021     = codes(rd("data/processed/flux/tree_flux_2021_multiheight.csv")$species),
  flux_2023     = codes(rd("data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 - Sheet1.csv")$`Species Code`),
  gas_2021      = codes(rd("data/processed/internal_gas/sample_data_only.csv")$Species.ID),
  ddpcr         = codes(rd("data/raw/ddpcr/ddPCR_meta_all_data.csv")$species),
  traits        = codes(rd("code/lib/trait_species.csv")$spcode))
miss <- setdiff(unlist(used), SP$species_code)
if (length(miss)) stop("species codes with no row in species.csv: ", paste(miss, collapse = ", "))

# ---- agreement with the existing definitions ---------------------------------------
if (length(setdiff(MODEL_SPECIES, SP$species_code)))
  stop("MODEL_SPECIES codes with no row in species.csv: ", paste(setdiff(MODEL_SPECIES, SP$species_code), collapse = ", "))
TR <- rd("code/lib/trait_species.csv")
stopifnot(all(SP$species[match(TR$spcode, SP$species_code)] == TR$species),
          all(SP$gymnosperm[match(TR$spcode, SP$species_code)] == (TR$gymnosperm == 1)),
          setequal(SP$species[SP$gymnosperm], GYMNOSPERMS))

write_csv(SP, file.path(OUT, "species.csv"), na = "")
cat(sprintf("species.csv: %d codes (%d in the model TAXONOMY, %d with traits)\n", nrow(SP),
            length(MODEL_SPECIES), nrow(TR)))
for (k in names(used)) cat(sprintf("  %-13s %2d codes\n", k, length(used[[k]])))
