#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 22_qpcr.R -- qpcr.csv (data/tables/): sequencing-facility qPCR, one row per
# sample x target
# ------------------------------------------------------------------------------
# 16S  data/raw/16s/16s_w_metadata.csv (archived as "16S metadata"; it is the
#      facility's 16S qPCR plate table: Sample ID, Well, 16S_per_ul), 2021 survey
#      and felled-oak samples.
# ITS  2021 survey: ITS_per_ul of data/raw/external/tree-microbiome/
#      tree_data_methanogen_group.csv (one row per 16S library; 33 reported as
#      "Undetected"). Felled oak: data/raw/external/tree-microbiome/bo_its_load.csv.
# Values are copies per µL of DNA extract as reported. "Undetected" is stored as
# copies NA with undetected = TRUE (zeros reported as 0 stay 0).
# Facility sample names are matched to samples.sample_id by the tree part in NCBI
# spelling (HM8 -> H8, ROBowser -> ROBOWSER, RM301 -> RM301A) + compartment + any
# lab-ID suffix; felled-oak names directly. Facility rows that are not samples
# (empty wells, muck, aquatic and other non-survey material) are left out and
# listed in the audit.
# Writes data/tables/qpcr.csv, outputs/audit/qpcr_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/tree_ids.R")
OUT <- "data/tables"
SAMPLES <- read.csv(file.path(OUT, "samples.csv"), stringsAsFactors = FALSE)

# facility / library name -> sample_id
# (names already spelled out -- AB3slash5, BB10prime, WA103star -- are spelled back first)
to_sample_id <- function(x) {
  m <- str_match(x, "^(.*?)(Inner|Outer|Mineral|Organic)([0-9]+[a-e])?$")
  tree <- gsub("star", "*", gsub("prime", "'", gsub("slash", "/", m[, 2])))
  tid <- ifelse(tree == "RO7/8", "RO7slash8", ncbi_tree_id(tree))          # soil between RO7 and RO8
  y <- ifelse(is.na(m[, 1]), x, paste0(tid, m[, 3], coalesce(m[, 4], "")))
  ifelse(y %in% SAMPLES$sample_id, y, NA_character_)
}

# ---- 16S ----------------------------------------------------------------------------
Q <- read.csv("data/raw/16s/16s_w_metadata.csv", check.names = FALSE, colClasses = "character") %>%
  filter(`Sample ID` != "None") %>%
  transmute(name = `Sample ID`, facility_well = Well, target = "16S", raw = `16S_per_ul`, source = "16s_w_metadata.csv")

# ---- ITS ----------------------------------------------------------------------------
T <- read.csv("data/raw/external/tree-microbiome/tree_data_methanogen_group.csv", check.names = FALSE, colClasses = "character")
names(T)[1] <- "library"
I21 <- T %>% distinct(library, ITS_per_ul, X16S_per_ul) %>%
  transmute(name = sub("\\.16S\\..*$", "", library), facility_well = NA_character_, target = "ITS", raw = ITS_per_ul,
            x16 = X16S_per_ul, source = "tree_data_methanogen_group.csv")
stopifnot(!anyDuplicated(I21$name))
IBO <- read.csv("data/raw/external/tree-microbiome/bo_its_load.csv", check.names = FALSE, colClasses = "character") %>%
  transmute(name = `Sample ID`, facility_well = Well, target = "ITS", raw = ITS_per_ul, source = "bo_its_load.csv")

A <- bind_rows(Q, I21 %>% select(-x16), IBO) %>%
  filter(!is.na(raw), raw != "") %>%
  mutate(sample_id = to_sample_id(name),
         undetected = grepl("^undetected$", raw, ignore.case = TRUE),
         copies_per_uL_extract = suppressWarnings(as.numeric(raw)))
stopifnot(all(!is.na(A$copies_per_uL_extract) | A$undetected))
unlinked <- A %>% filter(is.na(sample_id))
QPCR <- A %>% filter(!is.na(sample_id)) %>%
  transmute(sample_id, target, copies_per_uL_extract, undetected, facility_well, source) %>%
  arrange(target, sample_id)
k <- paste(QPCR$sample_id, QPCR$target)
if (anyDuplicated(k)) { print(QPCR[k %in% k[duplicated(k)], ]); stop("qpcr: sample x target not unique") }
write_csv(QPCR, file.path(OUT, "qpcr.csv"), na = "")

# ---- check: the 16S values carried in tree_data_methanogen_group.csv are the facility's
chk <- I21 %>% filter(!is.na(x16), x16 != "") %>% mutate(sample_id = to_sample_id(name)) %>%
  left_join(QPCR %>% filter(target == "16S") %>% select(sample_id, q16 = copies_per_uL_extract), by = "sample_id")
bad16 <- chk %>% filter(is.na(q16) | abs(as.numeric(x16) - q16) > 1e-9)

sink(out_path("qpcr_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("qpcr %d rows\n", nrow(QPCR)))
print(as.data.frame(QPCR %>% left_join(SAMPLES %>% select(sample_id, campaign), by = "sample_id") %>%
                      count(target, campaign, undetected)), row.names = FALSE)
cat(sprintf("\n16S values in tree_data_methanogen_group.csv vs the facility table: %d of %d agree\n", nrow(chk) - nrow(bad16), nrow(chk)))
if (nrow(bad16)) print(as.data.frame(bad16), row.names = FALSE)
cat(sprintf("\nFacility rows not linked to a sample (left out): %d\n", nrow(unlinked)))
print(as.data.frame(unlinked %>% select(name, target, raw, source)), row.names = FALSE)
cat("\nSamples with a 16S library but no 16S qPCR value:\n")
L16 <- I21 %>% mutate(sample_id = to_sample_id(name)) %>% filter(!is.na(sample_id), !sample_id %in% QPCR$sample_id[QPCR$target == "16S"])
print(L16$sample_id)
options(op); sink()
cat(readLines(out_path("qpcr_build.txt", "audit")), sep = "\n")
