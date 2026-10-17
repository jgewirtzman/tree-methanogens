#!/usr/bin/env Rscript
# ==============================================================================
# validate_tables.R -- check every table in data/tables/ against the dictionary
# ------------------------------------------------------------------------------
# Rules come from code/tables/dictionary*.csv (table, column, type, unit, key,
# range):
#   every column of every table is in the dictionary, and every dictionary row
#     names a column that exists;
#   type parses: character, integer, numeric, logical (TRUE/FALSE), date (ISO);
#   numeric and integer columns have a unit (count for counts, identifier for codes);
#   key PK: the PK columns (with any PKopt columns) are unique together; PK
#     columns are never missing. NN: never missing.
#   key FK:<table>.<column>: every non-missing value is present in that column;
#   range "min;max" for numbers and dates; lat_bbox / lon_bbox: within the
#     plot's bounding box (site.csv, plot corners) plus 100 m.
# Also: a tree's in_inventory stem belongs to that tree (inventory_stems.tree_id);
# picrust_pathways libraries and pathways exist in samples_16S and picrust_pathway_info.
# Reads only; stops with a list of failures.
#   Rscript code/tables/validate_tables.R
# ==============================================================================
DIR <- "data/tables"
# one dictionary file per group of tables: dictionary.csv plus dictionary_*.csv
D <- do.call(rbind, lapply(sort(Sys.glob("code/tables/dictionary*.csv")), utils::read.csv,
                           stringsAsFactors = FALSE, na.strings = ""))
FAILS <- character(0)
fail <- function(...) FAILS <<- c(FAILS, sprintf(...))

files <- list.files(DIR, "\\.csv$")
tabs <- sub("\\.csv$", "", files)
for (t in setdiff(unique(D$table), tabs)) fail("%s: in the dictionary but no %s/%s.csv", t, DIR, t)
for (t in setdiff(tabs, D$table)) fail("%s: table has no dictionary entries", t)
X <- setNames(lapply(intersect(tabs, D$table), function(t)
  utils::read.csv(file.path(DIR, paste0(t, ".csv")), colClasses = "character", na.strings = "",
                  check.names = FALSE)), intersect(tabs, D$table))

# plot bounding box from the site row (flat-earth offsets of the four corners)
S <- X$site
bbox <- local({
  if (is.null(S) || nrow(S) != 1) return(NULL)
  la0 <- as.numeric(S$ref_latitude); lo0 <- as.numeric(S$ref_longitude)
  th <- as.numeric(S$rotation_rad); L <- as.numeric(S$plot_side_m); buf <- 100
  px <- c(-buf, L + buf, -buf, L + buf); py <- c(-buf, -buf, L + buf, L + buf)
  e <- px * cos(th) + py * sin(th); n <- -px * sin(th) + py * cos(th)
  lat <- la0 + n / 111320; lon <- lo0 + e / (111320 * cos(la0 * pi / 180))
  list(lat = range(lat), lon = range(lon))
})

parse <- function(x, type) switch(type,
  character = x,
  integer   = { v <- suppressWarnings(as.numeric(x)); ifelse(is.finite(v) & v == round(v), v, NA) },
  numeric   = suppressWarnings(as.numeric(x)),
  logical   = ifelse(x %in% c("TRUE", "FALSE"), x == "TRUE", NA),
  date      = { v <- as.Date(x, format = "%Y-%m-%d"); ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", x), v, NA) },
  stop("unknown type ", type))

for (t in names(X)) {
  x <- X[[t]]; d <- D[D$table == t, ]
  # wide matrices: a row with column "*" stands for every column not listed by name;
  # key FKnames:<table>.<column> requires those column names to be values of that column
  w <- d[d$column %in% "*", ]
  if (nrow(w)) {
    wide <- setdiff(names(x), d$column)
    k <- ifelse(is.na(w$key[1]), "", w$key[1])
    if (startsWith(k, "FKnames:")) {
      ref <- strsplit(sub("^FKnames:", "", k), ".", fixed = TRUE)[[1]]
      miss <- setdiff(wide, X[[ref[1]]][[ref[2]]])
      if (length(miss)) fail("%s: %d columns not in %s.%s (e.g. %s)", t, length(miss), ref[1], ref[2], paste(head(miss, 3), collapse = ", "))
    }
    w$key <- NA
    d <- rbind(d[!d$column %in% "*", ], w[rep(1, length(wide)), , drop = FALSE])
    if (length(wide)) d$column[(nrow(d) - length(wide) + 1):nrow(d)] <- wide
  }
  for (cn in setdiff(names(x), d$column)) fail("%s.%s: column not in the dictionary", t, cn)
  for (cn in setdiff(d$column, names(x))) fail("%s.%s: in the dictionary but not in the table", t, cn)
  d <- d[d$column %in% names(x), ]
  for (i in seq_len(nrow(d))) {
    cn <- d$column[i]; ty <- d$type[i]; v <- x[[cn]]; p <- parse(v, ty)
    bad <- !is.na(v) & is.na(p)
    if (any(bad)) fail("%s.%s: %d values do not parse as %s (e.g. '%s')", t, cn, sum(bad), ty, v[bad][1])
    if (ty %in% c("integer", "numeric") && is.na(d$unit[i])) fail("%s.%s: numeric column without a unit", t, cn)
    k <- ifelse(is.na(d$key[i]), "", d$key[i])
    if (k %in% c("PK", "NN") && anyNA(v)) fail("%s.%s: %d missing values (%s)", t, cn, sum(is.na(v)), k)
    if (startsWith(k, "FK:")) {
      ref <- strsplit(sub("^FK:", "", k), ".", fixed = TRUE)[[1]]
      if (is.null(X[[ref[1]]])) fail("%s.%s: references missing table %s", t, cn, ref[1])
      else {
        miss <- unique(v[!is.na(v) & !v %in% X[[ref[1]]][[ref[2]]]])
        if (length(miss)) fail("%s.%s: %d values not in %s.%s (e.g. %s)", t, cn, length(miss), ref[1], ref[2],
                               paste(head(miss, 5), collapse = ", "))
      }
    }
    r <- d$range[i]
    if (!is.na(r)) {
      lim <- if (r == "lat_bbox") bbox$lat else if (r == "lon_bbox") bbox$lon else
        if (ty == "date") as.numeric(as.Date(strsplit(r, ";")[[1]])) else as.numeric(strsplit(r, ";")[[1]])
      out <- !is.na(p) & (as.numeric(p) < lim[1] | as.numeric(p) > lim[2])
      if (any(out)) fail("%s.%s: %d values outside %s (e.g. %s)", t, cn, sum(out), r, v[out][1])
    }
  }
  pk <- d$column[d$key %in% "PK"]; pko <- d$column[d$key %in% "PKopt"]
  if (!length(pk)) fail("%s: no primary key in the dictionary", t)
  else if (anyDuplicated(do.call(paste, c(x[c(pk, pko)], sep = "\r"))))
    fail("%s: primary key (%s) not unique: %d duplicates", t, paste(c(pk, pko), collapse = ", "),
         sum(duplicated(do.call(paste, c(x[c(pk, pko)], sep = "\r")))))
}

# a tree's census stem belongs to that tree
if (!is.null(X$trees) && !is.null(X$inventory_stems)) {
  tr <- X$trees[!is.na(X$trees$in_inventory), ]
  own <- X$inventory_stems$tree_id[match(tr$in_inventory, X$inventory_stems$stem_id)]
  bad <- which(!is.na(own) & own != tr$tree_id)
  if (length(bad)) fail("trees.in_inventory: %d stems belong to another tree_id (e.g. %s -> %s)", length(bad),
                        tr$tree_id[bad[1]], own[bad[1]])
}

# PICRUSt2 pathways (two-column key, so the references are checked here): libraries
# of the main 16S run; pathways described in picrust_pathway_info
if (!is.null(X$picrust_pathways)) {
  pp <- X$picrust_pathways
  L <- X$samples_16S$library_id[X$samples_16S$run %in% "main"]
  if (any(!pp$library_id %in% L)) fail("picrust_pathways.library_id: %d values not a main-run samples_16S.library_id",
                                       sum(!unique(pp$library_id) %in% L))
  if (any(!pp$pathway %in% X$picrust_pathway_info$pathway)) fail("picrust_pathways.pathway: not in picrust_pathway_info")
}

cat(sprintf("validate_tables: %d tables, %d columns checked\n", length(X), sum(D$table %in% names(X))))
for (t in names(X)) cat(sprintf("  %-16s %6d rows\n", t, nrow(X[[t]])))
if (length(FAILS)) {
  cat("FAILED:\n", paste0("  ", FAILS, collapse = "\n"), "\n", sep = "")
  stop(length(FAILS), " table check(s) failed", call. = FALSE)
}
cat("all checks passed\n")
