#!/usr/bin/env Rscript
# ==============================================================================
# 33_isotopes.R -- isotopes.csv (data/tables/)
# ------------------------------------------------------------------------------
# One row per Picarro run (408): every data/raw/internal_gas/picarro/*_results.csv
# file except the *_merged_results.csv copies (code/lib/isotope_samples.R,
# picarro_runs()). Values are the run means: d13CH4 = HR_Delta_iCH4_Raw_mean,
# CH4 = HR_12CH4_dry_mean, d13CO2 = Delta_Raw_iCO2_mean, CO2 = 12CO2_mean.
#
# sample_type, redone, in_whole_tree_set and exclusion_reason exactly as the
# archived table (code/zenodo/01_compile_datasets.R); the whole-tree set is
# isotope_whole_tree_samples() of code/lib/isotope_samples.R (the one definition:
# whole-tree single-borehole vials, no atmosphere / room air / incubations /
# heartwood-sapwood pairs / standards, the re-run replaces the original,
# CH4 >= 1.5 ppm).
# tree_id: tree_tags.csv (SampleName spellings, both campaigns), for every run of
# a tree sample -- whole-tree vials (2021 survey) and heartwood / sapwood vials
# (2024); NA for atmosphere, room air, incubations and standards.
# run_id: <results file stamp>_<autosampler position>.
#
# Checks: N_WHOLE_TREE runs in the whole-tree set (the manuscript's n), each with a tree_id.
#
# Writes data/tables/isotopes.csv, outputs/audit/isotopes_build.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/isotope_samples.R")
OUT <- "data/tables"
DIR <- "data/raw/internal_gas/picarro"
N_WHOLE_TREE <- 128   # whole-tree set: manuscript n = 128 trees

f <- list.files(DIR, pattern = "_results.csv$", full.names = TRUE)
f <- f[!grepl("merged", f)]
raw <- map_dfr(f, ~ read_csv(.x, show_col_types = FALSE) %>% mutate(run_file = basename(.x)))
stopifnot(isTRUE(all.equal(as.data.frame(select(raw, -run_file)), as.data.frame(picarro_runs()))))
keep <- isotope_whole_tree_samples(raw %>% select(-run_file))
stopifnot(nrow(keep) == N_WHOLE_TREE)

TT <- read.csv(file.path(OUT, "tree_tags.csv"), colClasses = "character", na.strings = "") %>%
  filter(source_file == "data/raw/internal_gas/picarro/*_results.csv", source_column == "SampleName") %>%
  distinct(tag = trimws(tag_recorded), tree_id, campaign, .keep_all = TRUE)
stopifnot(all(is.na(TT$species_recorded)), !anyDuplicated(TT$tag))

ISO <- raw %>% transmute(
  run_id = paste0(sub("_results.csv$", "", run_file), "_", Position), run_file,
  run_date = as.Date(as.character(Rundate), "%Y%m%d"),
  run_time = sub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1:\\2:\\3", sprintf("%06d", as.integer(Runtime))),
  position = as.integer(Position),
  sample_name = SampleName, d13CH4_permil = HR_Delta_iCH4_Raw_mean, CH4_ppm = HR_12CH4_dry_mean,
  d13CO2_permil = Delta_Raw_iCO2_mean, CO2_ppm = `12CO2_mean`) %>%
  mutate(sample_type = case_when(
           grepl("^Amb|room[[:space:]]*air", sample_name, ignore.case = TRUE) ~ "atmosphere or room air",
           grepl("^V", sample_name) ~ "incubation",
           grepl("[HS]$", sample_name) ~ "heartwood/sapwood pair",
           sample_name %in% ISO_STANDARDS ~ "standard",
           TRUE ~ "whole tree"),
         redone = sample_name %in% sub("_Redo$", "", sample_name[grepl("_Redo$", sample_name)]),
         in_whole_tree_set = sample_name %in% keep$SampleName,
         exclusion_reason = case_when(
           in_whole_tree_set ~ NA_character_,
           sample_type != "whole tree" ~ "not a whole-tree sample",
           redone ~ "replaced by its re-run",
           is.na(d13CH4_permil) ~ "no delta13C value",
           CH4_ppm < ISO_CH4_FLOOR ~ "CH4 below 1.5 ppm",
           TRUE ~ "delta13CH4 outside -115 to +25 permil"),
         tree_id = TT$tree_id[match(trimws(sample_name), TT$tag)],
         campaign = TT$campaign[match(trimws(sample_name), TT$tag)]) %>%
  relocate(tree_id, campaign, .after = sample_name)
cat(sprintf("isotopes: %d Picarro runs, %d in the whole-tree set\n", nrow(ISO), sum(ISO$in_whole_tree_set)))
stopifnot(!anyDuplicated(ISO$run_id), sum(ISO$in_whole_tree_set) == N_WHOLE_TREE,
          !anyNA(ISO$tree_id[ISO$in_whole_tree_set]))
# the set's tree_id (canonical) is the tree named by isotope_whole_tree_samples()
stopifnot(all(ISO$tree_id[match(keep$SampleName, ISO$sample_name)] ==
                TT$tree_id[match(keep$tree_id, TT$tag)]))

sink(out_path("isotopes_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("isotopes: %d runs in %d files; whole-tree set %d runs (%d trees)\n", nrow(ISO), n_distinct(ISO$run_file),
            sum(ISO$in_whole_tree_set), n_distinct(ISO$tree_id[ISO$in_whole_tree_set])))
print(table(ISO$sample_type, ISO$in_whole_tree_set, dnn = c("sample_type", "in_whole_tree_set")))
cat("\nexclusion reasons:\n"); print(table(ISO$exclusion_reason, useNA = "ifany"))
cat(sprintf("\nruns with a tree_id: %d (%s)\n", sum(!is.na(ISO$tree_id)), paste(names(table(ISO$campaign)), table(ISO$campaign), collapse = ", ")))
options(op); sink()
cat(readLines(out_path("isotopes_build.txt", "audit")), sep = "\n")

write_csv(ISO, file.path(OUT, "isotopes.csv"), na = "")
