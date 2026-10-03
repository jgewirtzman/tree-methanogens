#!/usr/bin/env Rscript
# ==============================================================================
# 04_tag_corrections.R -- correct misread tree tags in the stored flux tables
# ------------------------------------------------------------------------------
# Each row of code/lib/tag_corrections.csv names a file, a measurement date, the tag
# and species as recorded, and the corrected tag, with the evidence. A row is
# applied only where tag, species and date all match, so the step is idempotent and
# cannot touch another measurement; the recorded tag is kept in tag_as_recorded.
# Run from the repo root after 03_precision_and_detection.R.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })
fix <- read.csv("code/lib/tag_corrections.csv", stringsAsFactors = FALSE, colClasses = "character")
for (f in unique(fix$file)) {
  d <- read_csv(f, show_col_types = FALSE, col_types = cols(.default = col_guess(), `Tree Tag` = col_character()))
  if (!"tag_as_recorded" %in% names(d)) d$tag_as_recorded <- d$`Tree Tag`
  date_col <- if ("date_clean" %in% names(d)) "date_clean" else "Date"
  for (i in which(fix$file == f)) {
    hit <- d$tag_as_recorded == fix$tag_recorded[i] & d$`Species Code` == fix$species_recorded[i] &
           as.character(as.Date(d[[date_col]])) == fix$date[i]
    hit[is.na(hit)] <- FALSE
    d$`Tree Tag`[hit] <- fix$tag_corrected[i]
    cat(sprintf("  %s: %s %s on %s -> %s (%d row)\n", basename(f), fix$tag_recorded[i], fix$species_recorded[i],
                fix$date[i], fix$tag_corrected[i], sum(hit)))
    stopifnot(sum(hit) == 1)
  }
  write_csv(d, f)
}
