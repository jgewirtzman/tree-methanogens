#!/usr/bin/env Rscript
# ==============================================================================
# harvest_dictionary.R -- seed the column dictionary from the hand-written
# README_data.md. Run ONCE to bootstrap code/zenodo/column_dictionary.csv;
# after that the dictionary is the source of truth and is edited directly.
#
#   Rscript code/zenodo/harvest_dictionary.R
# ==============================================================================
suppressMessages({library(dplyr); library(readr)})

rd  <- readLines("data/compiled/README_data.md", warn = FALSE)
hdr <- grep("^## [0-9]+[.] `", rd)
ds  <- unlist(regmatches(rd[hdr], regexpr("[a-zA-Z0-9_]+[.]csv", rd[hdr])))

rows <- list()
for (i in seq_along(hdr)) {
  from <- hdr[i]
  to   <- if (i < length(hdr)) hdr[i + 1] - 1L else length(rd)
  blk  <- rd[from:to]
  tbl  <- grep("^\\|", blk, value = TRUE)
  tbl  <- tbl[!grepl("^\\|[-: |]+\\|$", tbl)]              # drop the rule row
  tbl  <- tbl[!grepl("^\\|\\s*Column\\s*\\|", tbl)]        # drop the header row
  for (ln in tbl) {
    cell <- trimws(strsplit(sub("^\\|", "", sub("\\|$", "", ln)), "\\|")[[1]])
    if (length(cell) < 4) next
    # Some rows document several columns at once ("year, month, day, hour").
    # Split them so the dictionary is keyed one row per actual column.
    for (cn in trimws(strsplit(cell[1], ",")[[1]])) {
      if (!nzchar(cn)) next
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds[i], column = gsub("`", "", cn), type = cell[2],
        unit = cell[3], description = cell[4], stringsAsFactors = FALSE)
    }
  }
}
D <- bind_rows(rows) %>% filter(nzchar(column)) %>% distinct(dataset, column, .keep_all = TRUE)

dir.create("code/zenodo", showWarnings = FALSE, recursive = TRUE)
write_csv(D, "code/zenodo/column_dictionary.csv")
cat(sprintf("harvested %d column definitions across %d datasets -> code/zenodo/column_dictionary.csv\n",
            nrow(D), length(unique(D$dataset))))
cat("datasets:", paste(sort(unique(D$dataset)), collapse = ", "), "\n")
