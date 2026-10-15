#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 10_species_traits.R -- species trait table from archived TRY and GWDD records
# (Methods S5; read by code/09_tables_stats/25_plant-traits.R and Fig S25)
# ------------------------------------------------------------------------------
# Species: code/lib/trait_species.csv (16 YMF species; binomials as in TRY
# AccSpeciesName and GWDD `species`; gymnosperm flag).
# Inputs (Zenodo drop-in; see data/raw/external/README.md):
#   data/raw/external/gwdd/gwdd_v2.2_records_16spp.csv    GWDD v2.2 (with bark)
#   data/raw/external/try/try_request47607_records.csv    TRY request 47607
#
# GWDD (Fischer et al. 2026): species mean of `wsg` over the species' records.
#   wood density: type_tissue missing or not "bark"; bark density: type_tissue
#   "bark". Species-level matches only (all 16 have wood records).
# TRY (Kattge et al. 2020): 11 traits. Records are dropped, in this order:
#   Q1 duplicates (OrigObsDataID non-empty);
#   Q2 outliers (ErrorRisk > 4; records without ErrorRisk kept);
#   Q3 one quantity per trait: 6 maximum rooting depth only (no "Minimum");
#      59 maximum age (DataID 205 or 3686, or 190 with ValueKind "Maximum");
#      896 no FRED "median" diameter; 1158 single-exponential k (DataID 3083);
#      7 no "Original term" rows (dataset 162 repeats its harmonised type);
#   Q4 values republished from another dataset: (61, 662), (61, 68), (6, 662),
#      (407, 313) as (TraitID, DatasetID).
#   Value: StdValue for 407, 6, 896, 1781, 59; numeric OrigValueStr for 4079,
#   1260, 617, 1158, 61 (0-5 scores for 61; USDA text classes drop out);
#   mycorrhiza (7) coded from OrigValueStr: non-/no/NM -> NA, ecto/EcM/EM -> 1,
#   AM/arbuscular/VAM/endomycorrhiza -> 0, other -> NA.
#   Species value = mean over datasets of the within-dataset means.
# Rooting depth is built and archived but fails R1 in 25_plant-traits.R (6 of 10).
# No wood pH: no species-level source covers the analysis species.
#
# Writes outputs/data/species_traits.csv        (one row per species)
#        outputs/data/species_trait_sources.csv (one row per column x species x
#          contributing dataset, with record IDs: TRY ObsDataID, GWDD id)
#        outputs/audit/species_traits.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))

sp <- read_csv("code/lib/trait_species.csv", show_col_types = FALSE)
stopifnot(nrow(sp) == 16, !anyDuplicated(sp$species))

# ---- GWDD -----------------------------------------------------------------------
gw <- read_csv("data/raw/external/gwdd/gwdd_v2.2_records_16spp.csv", show_col_types = FALSE,
               col_types = cols(id = "c", wsg = "d", type_tissue = "c", .default = "c")) %>%
  filter(species %in% sp$species, is.finite(wsg))
gw <- bind_rows(
  gw %>% filter(is.na(type_tissue) | !grepl("^bark$", type_tissue, ignore.case = TRUE)) %>% mutate(column = "wood_density_gcm3"),
  gw %>% filter(!is.na(type_tissue) & grepl("bark", type_tissue, ignore.case = TRUE)) %>% mutate(column = "bark_density_gcm3")) %>%
  mutate(spcode = sp$spcode[match(species, sp$species)])
gw_val <- gw %>% group_by(column, spcode) %>% summarise(value = mean(wsg), .groups = "drop")
stopifnot(all(sp$spcode %in% gw_val$spcode[gw_val$column == "wood_density_gcm3"]))
gw_src <- gw %>% group_by(column, spcode, dataset = source_short) %>%
  summarise(reference = paste(unique(source_long), collapse = " | "), n_records = n(),
            dataset_value = mean(wsg), record_ids = paste(id, collapse = ";"), .groups = "drop") %>%
  mutate(source = "GWDD", dataset_id = NA_character_)

# ---- TRY ------------------------------------------------------------------------
try_cols <- c(`4079` = "try_bark_volume_rel_wood", `1260` = "try_wood_CN_ratio", `407` = "try_stem_C",
  `617` = "try_bark_C", `1158` = "try_CWD_stem_decomp_rate_k", `61` = "try_waterlogging_tolerance",
  `6` = "try_rooting_depth", `896` = "try_fine_root_diameter", `1781` = "try_fine_root_tissue_density",
  `7` = "try_ectomycorrhizal", `59` = "try_plant_longevity")
r <- read_csv("data/raw/external/try/try_request47607_records.csv", col_types = cols(.default = "c"),
              na = "") %>%
  filter(AccSpeciesName %in% sp$species, TraitID %in% names(try_cols)) %>%
  mutate(column = unname(try_cols[TraitID]), spcode = sp$spcode[match(AccSpeciesName, sp$species)],
         TraitID = as.integer(TraitID), DatasetID = as.integer(DatasetID), DataID = as.integer(DataID),
         ErrorRisk = as.numeric(ErrorRisk))
o <- r$OrigValueStr
r <- r %>% mutate(excluded = case_when(
  !is.na(OrigObsDataID) & OrigObsDataID != "" ~ "Q1 duplicate",
  !is.na(ErrorRisk) & ErrorRisk > 4 ~ "Q2 ErrorRisk > 4",
  (TraitID == 6    & grepl("Minimum", DataName)) |
  (TraitID == 59   & !(DataID %in% c(205, 3686) | (DataID == 190 & ValueKindName %in% "Maximum"))) |
  (TraitID == 896  & grepl("median", DataName)) |
  (TraitID == 1158 & DataID != 3083) |
  (TraitID == 7    & grepl("^Original term", DataName)) ~ "Q3 other quantity",
  paste(TraitID, DatasetID) %in% c("61 662", "61 68", "6 662", "407 313") ~ "Q4 republished",
  .default = NA_character_),
  value = case_when(
    TraitID %in% c(407, 6, 896, 1781, 59) ~ suppressWarnings(as.numeric(StdValue)),
    TraitID %in% c(4079, 1260, 617, 1158, 61) ~ suppressWarnings(as.numeric(o)),
    TraitID == 7 & grepl("^non|^no$|^NM$", o, ignore.case = TRUE) ~ NA_real_,
    TraitID == 7 & grepl("^ecto|^ecm$|^em$", o, ignore.case = TRUE) ~ 1,
    TraitID == 7 & grepl("^AM$|arbuscular|VAM|endomycorrhiza", o, ignore.case = TRUE) ~ 0,
    .default = NA_real_),
  excluded = ifelse(is.na(excluded) & !is.finite(value), "no usable value", excluded))
kept <- r %>% filter(is.na(excluded))
try_src <- kept %>% group_by(column, spcode, dataset_id = as.character(DatasetID)) %>%
  summarise(dataset = first(Dataset), reference = paste(unique(Reference), collapse = " | "),
            n_records = n(), dataset_value = mean(value), record_ids = paste(ObsDataID, collapse = ";"),
            .groups = "drop") %>% mutate(source = "TRY")
try_val <- try_src %>% group_by(column, spcode) %>% summarise(value = mean(dataset_value), .groups = "drop")

# ---- assemble ---------------------------------------------------------------------
cols <- c("wood_density_gcm3", "bark_density_gcm3", unname(try_cols))
vals <- bind_rows(gw_val, try_val)
out <- sp %>% select(spcode, species, gymnosperm)
for (cn in cols) out[[cn]] <- with(vals[vals$column == cn, ], value[match(out$spcode, spcode)])
src <- bind_rows(gw_src, try_src) %>%
  select(column, spcode, source, dataset_id, dataset, reference, n_records, dataset_value, record_ids) %>%
  arrange(match(column, cols), spcode, source, as.numeric(dataset_id), dataset)

# every value has at least one source row, and every source row a value
have <- out %>% pivot_longer(all_of(cols), names_to = "column") %>% filter(!is.na(value))
miss <- anti_join(have, src, by = c("column", "spcode"))
if (nrow(miss)) stop("values without a source row: ", paste(miss$column, miss$spcode, collapse = "; "))
stopifnot(nrow(anti_join(distinct(src, column, spcode), have, by = c("column", "spcode"))) == 0,
          all(src$n_records > 0))

write_csv(out, out_path("species_traits.csv"))
write_csv(src, out_path("species_trait_sources.csv"))

sink(out_path("species_traits.txt")); op <- options(width = 160)
cat(sprintf("Species traits: %d species; GWDD %d records; TRY %d records, %d kept\n\n",
            nrow(sp), nrow(gw), nrow(r), nrow(kept)))
cat("TRY records by trait and first exclusion rule:\n")
print(as.data.frame(r %>% mutate(excluded = replace_na(excluded, "kept")) %>% count(TraitID, column, excluded) %>%
  pivot_wider(names_from = excluded, values_from = n, values_fill = 0) %>%
  select(TraitID, column, kept, any_of(c("Q1 duplicate", "Q2 ErrorRisk > 4", "Q3 other quantity", "Q4 republished")), everything())), row.names = FALSE)
cat("\nContributing TRY datasets (species of 16):\n")
print(as.data.frame(src %>% filter(source == "TRY") %>% group_by(column, dataset_id, dataset) %>%
  summarise(species = n(), records = sum(n_records), .groups = "drop") %>%
  mutate(dataset = substr(dataset, 1, 60))), row.names = FALSE)
cat("\nGWDD sources and records per column:\n")
print(as.data.frame(src %>% filter(source == "GWDD") %>% group_by(column) %>%
  summarise(species = n_distinct(spcode), sources = n_distinct(dataset), records = sum(n_records))), row.names = FALSE)
cat("\nSpecies values:\n")
print(as.data.frame(out %>% select(-species) %>% mutate(across(all_of(cols), ~ signif(., 4)))), row.names = FALSE)
options(op)
sink()
cat(readLines(out_path("species_traits.txt")), sep = "\n")
