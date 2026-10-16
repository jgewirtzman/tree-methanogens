#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 20_samples.R -- samples.csv (data/tables/): one row per physical DNA sample
# ------------------------------------------------------------------------------
# 2021 survey  wood core sections (inner, outer) and soil horizons (organic,
#              mineral) of each tree, from the four extraction guides
#              (data/raw/ddpcr/Data Files/{inner,outer,organic,mineral}_guide.csv:
#              lab ID, mass added to the tube, extraction plate and well).
#              Guide wells marked "Blank Well? = Yes" held no sample (sample missing,
#              insufficient mass, sectioned twice) and are not samples -- except three
#              that have a 16S library (SM15Organic, YB12Mineral, YB15Inner: kept,
#              well_marked_blank = TRUE, no mass). Samples sequenced but absent from the
#              guides (the RM301 soils, one bag labelled "301" without A/B; soil between
#              RO7 and RO8; WPBerenice mineral) come from the 16S library map.
# 2022 felled oak  tissues by height and the ground samples around it
#              (data/raw/field_data/black_oak/bo_extraction_mass.csv).
#
# SAMPLE IDs
#   2021        NCBI isolation_source form: the tree in NCBI spelling (its tree_id;
#               code/lib/tree_ids.R with code/lib/tree_name_aliases.csv) + compartment,
#               e.g. AB10Inner, AB3slash5Mineral, BB10primeOuter, H8Inner (16S map:
#               HM8Inner), ROBOWSERInner (guide: BOWSER; 16S map: ROBowser). Where one tree has two
#               samples of one compartment (YB9, WA102, WA14: two lab IDs), the lab
#               ID is appended (YB9Inner44c), as in the 16S library names and the
#               facility qPCR table.
#   felled oak  the sample ID of the extraction sheet without hyphens
#               (QUVE50SAP, QUVEBRANCH4BARK), as in the 16S library names.
#               The RM301 soils (one bag labelled "301") are RM301AMineral and
#               RM301AOrganic, following the alias RM301 -> RM301A. Soil sampled
#               between RO7 and RO8 (16S map RO7/8) is RO7slash8Mineral, tree_id NA.
#
# Writes data/tables/samples.csv, outputs/audit/samples_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/tree_mentions.R")          # ncbi_tree_id(), S21, OAK, FELLED_OAK_ID
OUT <- "data/tables"
NO_ALIAS <- data.frame(spelling = character(0), tree_id = character(0))
RO78 <- "RO7slash8"                          # soil sampled between RO7 and RO8: no single tree
TREES <- read.csv(file.path(OUT, "trees.csv"), stringsAsFactors = FALSE)

# ==============================================================================
# 1. 2021 SURVEY: EXTRACTION GUIDES
# ==============================================================================
G <- bind_rows(lapply(c("inner", "outer", "organic", "mineral"), function(k)
  read.csv(sprintf("data/raw/ddpcr/Data Files/%s_guide.csv", k), skip = 1, check.names = FALSE,
           colClasses = "character", na.strings = "") %>% mutate(guide = k)))
G <- G %>%
  mutate(label = gsub("_x000D_|\r", "", `Inner Core Sample ID`),
         line1 = trimws(map_chr(strsplit(label, "\n"), 1)),
         line2 = trimws(map_chr(strsplit(label, "\n"), 2)),
         lab_id = trimws(sub(".*Lab ID:\\s*", "", label)),
         field_name = sub("^[A-Z]{4}\\s*_\\s*", "", line1),
         species_recorded = substr(line1, 1, 4),
         material = tolower(sub("_.*", "", line2)),
         position = tolower(sub(".*_", "", line2)),
         blank = grepl("yes", `Blank Well?`, ignore.case = TRUE),
         extraction_mass_mg = suppressWarnings(as.numeric(`Sample Mass Added to Tube (mg)`)),
         extraction_plate = as.integer(sub("Plate\\s*", "", `Extraction Plate ID`)),
         extraction_well = trimws(`Plate Position`),
         notes = pmap_chr(list(Notes, `Sample Notes`), function(a, b) {
           x <- trimws(c(a, b)); x <- x[!is.na(x) & x != ""]; if (length(x)) paste(x, collapse = "; ") else NA_character_ }))
stopifnot(nrow(G) == 632, !anyDuplicated(G$lab_id), all(G$position == G$guide),
          !anyDuplicated(paste(G$extraction_plate, G$extraction_well)))
M <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", check.names = FALSE, colClasses = "character")
G$sequenced <- G$lab_id %in% M$lab_id
blank_wells <- G %>% filter(blank, !sequenced)
G <- G %>% filter(!blank | sequenced)

POS <- c(inner = "Inner", outer = "Outer", organic = "Organic", mineral = "Mineral")
S1 <- G %>%
  mutate(tree_id = ncbi_tree_id(field_name), sample_id = paste0(tree_id, POS[position])) %>%
  group_by(sample_id) %>% mutate(sample_id = if (n() > 1) paste0(sample_id, lab_id) else sample_id) %>% ungroup() %>%
  transmute(sample_id, tree_id, campaign = S21, material, position, height_cm = NA_real_, replicate = NA_character_,
            extraction_mass_mg, extraction_plate, extraction_well, lab_id, well_marked_blank = blank, notes,
            source = "extraction guide")

# ---- sequenced samples that are not in the guides ------------------------------------
# 16S library map entries naming a tree and a compartment that no guide sample matches.
extra <- M %>% filter(core_type %in% POS) %>%
  mutate(sample_id = paste0(coalesce(ncbi_tree_id(seq_id), ncbi_tree_id(seq_id, NO_ALIAS)), core_type)) %>%
  filter(lab_id %in% c("MISC", ""), !sample_id %in% S1$sample_id,
         !ForwardFastqFile %in% "206_B01_16S_S3_R1_001.fastq") %>%     # library-map patch (RO104 Inner), see 23_16s.R
  transmute(sample_id, tree_id = na_if(ncbi_tree_id(seq_id), RO78), campaign = S21,
            material = ifelse(core_type %in% c("Inner", "Outer"), "wood", "soil"), position = tolower(core_type),
            height_cm = NA_real_, replicate = NA_character_, extraction_mass_mg = NA_real_,
            extraction_plate = NA_integer_, extraction_well = NA_character_, lab_id = NA_character_,
            well_marked_blank = NA, notes = ifelse(seq_id == "RO7/8", "soil sampled between RO7 and RO8; in the 16S library map only",
                           "in the 16S library map only (the guide wells for this tree and compartment are blank)"),
            source = "16S library map")

# ==============================================================================
# 2. 2022 FELLED OAK
# ==============================================================================
B <- read.csv("data/raw/field_data/black_oak/bo_extraction_mass.csv", check.names = FALSE, stringsAsFactors = FALSE)
bo_lab <- read.csv("data/raw/ddpcr/black_oak_mcrA.csv", check.names = FALSE, stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8-BOM") %>% filter(`Sample ID` != "")
TISSUE <- tribble(~Component, ~material, ~position,
  "Sapwood", "wood", "sapwood",        "Heartwood", "wood", "heartwood",   "Bark", "bark", "bark",
  "Branch", "wood", "branch",          "Coarse Roots", "root", "coarse root", "Fine Roots", "root", "fine root",
  "Foliage", "foliage", "foliage",     "Heart Rot", "wood", "heart rot",   "Leaf Litter", "litter", "litter",
  "Mineral Soil", "soil", "mineral",   "Organic Soil", "soil", "organic")
S2 <- B %>% left_join(TISSUE, by = "Component") %>%
  mutate(position = case_when(grepl("-Bark$", `Sample Name`) ~ "branch bark",
                              grepl("Branch.*Interior", `Sample Name`) ~ "branch interior",
                              grepl("Epiderm", `Sample Name`) ~ "coarse root epidermis",
                              grepl("Coarse.*Interior", `Sample Name`) ~ "coarse root interior",
                              TRUE ~ position),
         material = ifelse(position == "branch bark", "bark", material),
         replicate = str_extract(`Sample Name`, "[0-9]+(\\+[0-9]+)?(?=(-[A-Za-z]+)?$)")) %>%
  transmute(sample_id = gsub("-", "", `Sample ID`), tree_id = FELLED_OAK_ID, campaign = OAK, material, position,
            height_cm = suppressWarnings(as.numeric(`Height (cm)`)), replicate,
            extraction_mass_mg = `Sample Mass Added to Tube (mg)`, extraction_plate = NA_integer_,
            extraction_well = NA_character_, lab_id = bo_lab$`Sample No`[match(`Sample ID`, bo_lab$`Sample ID`)],
            well_marked_blank = NA,
            notes = ifelse(Component == "Sapwood", "sapwood cores 1 and 2 combined when grinding", NA_character_),
            source = "felled-oak extraction sheet")
stopifnot(nrow(S2) == 60, !anyNA(S2$material), all(bo_lab$`Sample ID` %in% B$`Sample ID`))

# ==============================================================================
# 3. ASSEMBLE AND CHECK
# ==============================================================================
SAMPLES <- bind_rows(S1, extra, S2) %>% arrange(campaign, sample_id)
if (anyDuplicated(SAMPLES$sample_id)) stop("duplicate sample_id: ", paste(SAMPLES$sample_id[duplicated(SAMPLES$sample_id)], collapse = ", "))
stopifnot(!any(grepl("[\r\n ]", SAMPLES$sample_id)))
unlinked <- SAMPLES %>% filter(is.na(tree_id) | !tree_id %in% TREES$tree_id)
SAMPLES$tree_id[!SAMPLES$tree_id %in% TREES$tree_id] <- NA
write_csv(SAMPLES, file.path(OUT, "samples.csv"), na = "")

sink(out_path("samples_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("samples %d (2021 survey %d from the guides + %d from the 16S map only; felled oak %d)\n",
            nrow(SAMPLES), nrow(S1), nrow(extra), nrow(S2)))
print(as.data.frame(count(SAMPLES, campaign, material, position)), row.names = FALSE)
cat(sprintf("\nGuide wells marked blank (no sample; not in samples.csv): %d\n", nrow(blank_wells)))
cat("Guide wells marked blank that have a 16S library (kept, well_marked_blank = TRUE):\n")
print(as.data.frame(SAMPLES %>% filter(well_marked_blank %in% TRUE) %>% select(sample_id, lab_id, extraction_plate, extraction_well, notes)), row.names = FALSE)
print(as.data.frame(count(blank_wells, position)), row.names = FALSE)
cat("\nSamples added from the 16S library map:\n"); print(as.data.frame(extra %>% select(sample_id, tree_id)), row.names = FALSE)
cat("\nOne tree, two samples of one compartment (lab ID appended):\n")
print(as.data.frame(SAMPLES %>% filter(grepl("[0-9]+[a-e]$", sample_id) & campaign == S21) %>% select(sample_id, tree_id, lab_id, extraction_mass_mg, notes)), row.names = FALSE)
cat("\n2021 samples whose guide or 16S spelling differs from the tree_id spelling:\n")
al <- bind_rows(G %>% transmute(spelling = field_name, src = "guide"), M %>% transmute(spelling = seq_id, src = "16S map")) %>%
  distinct() %>% mutate(no_alias = ncbi_tree_id(spelling, NO_ALIAS), tree_id = ncbi_tree_id(spelling)) %>%
  filter(!is.na(tree_id), no_alias != tree_id | spelling != tree_id, spelling != "MISC") %>% arrange(tree_id)
print(as.data.frame(al), row.names = FALSE)
cat(sprintf("\nSamples not linked to a tree: %d\n", nrow(unlinked)))
if (nrow(unlinked)) print(as.data.frame(unlinked %>% select(sample_id, tree_id, source)), row.names = FALSE)
cat(sprintf("Samples without an extraction mass: %d\n", sum(is.na(SAMPLES$extraction_mass_mg))))
print(as.data.frame(SAMPLES %>% filter(is.na(extraction_mass_mg)) %>% select(sample_id, source, notes)), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("samples_build.txt", "audit")), sep = "\n")
