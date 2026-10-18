#!/usr/bin/env Rscript
# ==============================================================================
# 03_generate_data_readme.R -- write data/compiled/README_data.md: the canonical tables
# (data/tables/, archived as they are) documented from their column dictionaries,
# code/tables/dictionary*.csv, plus the list of data/compiled/results/.
#
# Generated, so the dictionary is the only thing to maintain, and the check below
# means a new column cannot ship undocumented.
#
#   Rscript code/zenodo/03_generate_data_readme.R          # write the README
#   Rscript code/zenodo/03_generate_data_readme.R --check  # report gaps, write nothing
#
# Exits non-zero if any column of a table has no dictionary entry, or if the dictionary
# describes a table or column that no longer exists.
# ==============================================================================
suppressMessages({library(dplyr); library(readr)})

TABLES  <- "data/tables"
OUT     <- "data/compiled/README_data.md"
check_only <- "--check" %in% commandArgs(trailingOnly = TRUE)

stopifnot(dir.exists(TABLES))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
D <- bind_rows(lapply(list.files("code/tables", "^dictionary.*[.]csv$", full.names = TRUE),
                      function(f) read_csv(f, col_types = cols(.default = "c")))) %>%
  transmute(dataset = paste0(table, ".csv"), column, type, unit, description) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "", as.character(.x))))

files <- sort(basename(list.files(TABLES, "[.]csv$")))

# Datasets whose columns are sample identifiers rather than fixed fields; a
# per-column dictionary is meaningless for these, so they are described once.
WIDE <- c(asv_counts_16S.csv =
            "One row per ASV of the main 16S run, one column per sequenced library. Values are read counts. Column names are `library_id` in samples_16S.csv (which gives the sample); the first column, `feature_id`, is the ASV identifier used in asv_taxonomy_16S.csv.",
          asv_counts_16S_felled_oak.csv =
            "As asv_counts_16S.csv, for the separately denoised felled-oak run (its `feature_id`s are not the same sequences as the main run's); taxonomy in asv_taxonomy_16S_felled_oak.csv.")

# ---- gather the tables' actual shape --------------------------------------
shape <- lapply(files, function(f) {
  p <- file.path(TABLES, f)
  h <- read.csv(p, nrows = 1, check.names = FALSE, stringsAsFactors = FALSE)
  n <- tryCatch(length(count.fields(p, sep = ",", quote = "\"")) - 1L,
                error = function(e) NA_integer_)
  list(file = f, cols = names(h), nrow = n,
       types = vapply(h, function(x) class(x)[1], character(1)))
})
names(shape) <- files

# ---- coverage check ---------------------------------------------------------
missing <- list(); orphan <- list()
for (f in files) {
  if (f %in% names(WIDE)) next
  have <- D$column[D$dataset == f]
  want <- shape[[f]]$cols
  m <- setdiff(want, have); o <- setdiff(have, want)
  if (length(m)) missing[[f]] <- m
  if (length(o)) orphan[[f]]  <- o
}
dead <- setdiff(unique(D$dataset), c(files, names(WIDE)))

cat("=== dictionary coverage ===\n")
cat(sprintf("  tables: %d | dictionary entries: %d\n", length(files), nrow(D)))
if (length(missing)) {
  cat("\n  COLUMNS WITH NO DICTIONARY ENTRY:\n")
  for (f in names(missing))
    cat(sprintf("    %-34s %s\n", f, paste(head(missing[[f]], 12), collapse = ", ")))
}
if (length(orphan)) {
  cat("\n  dictionary entries for columns no longer present:\n")
  for (f in names(orphan))
    cat(sprintf("    %-34s %s\n", f, paste(head(orphan[[f]], 12), collapse = ", ")))
}
if (length(dead)) cat("\n  dictionary describes tables not in data/tables:\n    ",
                      paste(dead, collapse = ", "), "\n")

nbad <- length(missing) + length(orphan) + length(dead)
if (check_only) {
  cat(sprintf("\n%s\n", if (nbad) "GAPS FOUND" else "dictionary and tables agree"))
  quit(status = if (nbad) 1L else 0L)
}
if (length(missing))
  stop("refusing to write a README with undocumented columns; add them to ",
       "code/tables/dictionary*.csv or run with --check to list them", call. = FALSE)

# ---- emit -------------------------------------------------------------------
fmt_n <- function(x) formatC(x, big.mark = ",", format = "d")
L <- c(
"# Data tables",
"",
"The canonical tables of the Yale-Myers Forest tree methane project, in `data/tables/`:",
"one CSV per kind of observation, with documented columns and units, for reanalysis and",
"meta-analysis.",
"",
"**Citation**: Gewirtzman, J. (2026). jgewirtzman/tree-methanogens: Tree Microbiomes",
"and Methane Exchange in Upland Forests — Code & Data. Zenodo.",
"https://doi.org/10.5281/zenodo.18779715",
"",
"**License**: CC-BY 4.0",
"",
"**Scope**: the canonical tables of the analysis, one copy of each data item, one table",
"per kind of observation, each row kept and flagged rather than deleted. Every analysis",
"in the code repository reads these tables. Results computed from them (stand budget,",
"scaling grid, model skill, the isotope summary and the paper's tables) are in",
"`data/compiled/results/`. Raw instrument files and field sheets are in `data/raw/`;",
"`data/processed/` holds the stage-A products the tables are built from.",
"",
"*This file is generated by `code/zenodo/03_generate_data_readme.R` from",
"`code/tables/dictionary*.csv`. Edit the dictionaries, not this file.*",
"",
"---",
"",
"## How the datasets relate",
"",
"**Keys.** `species.csv` (species_code), `site.csv`, `collars.csv` (collar_id),",
"`trees.csv` (tree_id: one physical tree; species, DBH and its source, landscape",
"position, dead status) and `tree_tags.csv` (every spelling of a tree name in the source",
"files and the tree_id it refers to). `inventory_stems.csv` is the whole stand from the",
"2018-2019 censuses and the by-tag survey, `in_stand` marking the budget basis.",
"",
sprintf("**CH4 flux.** `flux_stem.csv` (%s deployments) and `flux_soil.csv` (%s) hold every", fmt_n(shape[["flux_stem.csv"]]$nrow), fmt_n(shape[["flux_soil.csv"]]$nrow)),
"chamber deployment once, from all campaigns, with the goFlux fit diagnostics, the",
"minimum detectable flux and detection status, and the flags `qc_pass` (procedural",
"failure only), `dead_stem`, `untagged_rescue`, `in_rf_training` and `exclusion_reason`.",
"Filter on `in_rf_training` to reproduce the random forests; use all rows to describe",
"the measurements. `measurement_id` is the goFlux UniqueID.",
"",
"**Environment.** `weather.csv` (hourly station record in true time, screened),",
"`soil_readings.csv` (soil temperature, moisture and horizons at trees and collars),",
"`moisture_survey.csv` and `stream_points.csv` (the December 2020 moisture survey).",
"",
"**Samples and genes.** `samples.csv` is one row per DNA sample (tree, compartment,",
"extraction mass); `ddpcr.csv` holds every ddPCR run (copies per gram, QC flags),",
"`qpcr.csv` the facility 16S and ITS qPCR, `samples_16S.csv` the sequenced libraries,",
"and `asv_counts_16S*.csv` / `asv_taxonomy_16S*.csv` the ASV tables of the main and",
"felled-oak runs. `picrust_pathways.csv` holds the PICRUSt2 MetaCyc pathway abundance",
"per 16S library (organelle ASVs left out; also without the ASVs predicted to carry mcrA or",
"pmoA) and `picrust_pathway_info.csv` the pathway names.",
"",
"**Wood, gas and isotopes.** `wood_cores.csv` (density and moisture by radial",
"section), `internal_gas.csv` (one row per gas sample) and `isotopes.csv` (every Picarro",
"run, with `in_whole_tree_set` marking the samples analysed).",
"",
"**Plant traits.** `species_traits.csv` holds one value per species and trait (16",
"species), built from TRY (request 47607) and GWDD v2.2 records; `species_trait_sources.csv`",
"lists, for every value, the contributing TRY datasets or GWDD sources with their record",
"IDs. Cite TRY and each contributing dataset's reference when reusing these values.",
"",
"---",
"")
for (i in seq_along(files)) {
  f <- files[i]; sh <- shape[[f]]
  L <- c(L, sprintf("## %d. `%s`", i, f), "",
         sprintf("%s rows x %d columns.", fmt_n(sh$nrow), length(sh$cols)), "")
  if (f %in% names(WIDE)) {
    L <- c(L, WIDE[[f]], "", "---", "")
    next
  }
  dd <- D[D$dataset == f, ]
  L <- c(L, "| Column | Type | Unit | Description |",
            "|--------|------|------|-------------|")
  for (cn in sh$cols) {
    r <- dd[dd$column == cn, ][1, ]
    L <- c(L, sprintf("| %s | %s | %s | %s |", cn, r$type, r$unit, r$description))
  }
  L <- c(L, "", "---", "")
}
RES <- sort(list.files("data/compiled/results", "[.]csv$"))
if (length(RES)) L <- c(L, "## Results (`data/compiled/results/`)", "",
  "Copied from the pipeline's outputs without recomputation; each is the source of the",
  "numbers the paper quotes. Generated by `code/zenodo/02_compile_results.R`.", "",
  paste0("- `", RES, "`"), "")
writeLines(L, OUT)
cat(sprintf("\nwrote %s (%d tables, %d documented columns)\n",
            OUT, length(files), sum(vapply(shape, function(s) length(s$cols), 1L))))
