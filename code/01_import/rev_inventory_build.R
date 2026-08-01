#!/usr/bin/env Rscript
# ==============================================================================
# rev_inventory_build.R  --  THE stem inventory, rebuilt from raw.
# ------------------------------------------------------------------------------
# Run from the repo root:  Rscript code/01_import/rev_inventory_build.R
# Reads  data/raw/inventory/  (Zenodo drop-in)
# Writes outputs/tables/inventory_stems.csv   -- one row per stem, all sources
#
# This replaces a workflow in which the mapped tree list existed ONLY inside
# outputs/tables/tree_flux_predictions.csv, an untracked file that
# rev_predict_tree_flux_current.R read as its own input and overwrote as its own
# output, and that rev_budget_canonical.R also wrote with a different stem set.
# Running those two scripts in the wrong order destroyed the list irrecoverably.
# Everything now derives from raw files by a deterministic script instead.
#
# THREE CORRECTIONS TO THE PREVIOUS INVENTORY
#
# 1. THE bytag SURVEY WAS BEING DISCARDED, AND ITS STEMS GIVEN RANDOM POSITIONS.
#    deprecated/metadata_scripts/fg_aligned.R selected only lat.1/lon.1 from the
#    bytag table and never its PX/PY columns, so those stems arrived with PX = NA
#    and were filled by `runif(n(), 0, 200)` -- unseeded. The geodetic transform
#    of those random values then OVERWROTE the coordinates. (lat.1/lon.1 would not
#    have helped: that column is a single constant, the plot origin, repeated.)
#    The real PX/PY are in the bytag file and are recovered here. They place every
#    one of those stems in quadrats 506-908, i.e. the block PX/PY 100-200 -- a
#    hectare that fg19 censused at 5.3 m2/ha basal area against 43.1 m2/ha over
#    the rest of the plot. They are not duplicates: median distance to the nearest
#    fg19 stem is 4.61 m, against 0.54 m for fg19-to-fg19.
#
# 2. TAG IDs ARE NOT COMPARABLE ACROSS SOURCES. fg19 tags run 10000+, bytag runs
#    34-998. The ~25 numeric collisions are DIFFERENT TREES -- the species
#    disagree on 9 of 25 and the diameters differ by 3x to 386x. De-duplicating by
#    tag across sources silently deletes real stems, so we de-duplicate strictly
#    within each source and key on (source, tag).
#
# 3. DIAMETER OUTLIERS ARE DECIMAL TYPOS AND ARE REPAIRED, NOT DROPPED.
#    Units, established three ways: 229 census stems also appear in the flux
#    training data, where dbh_m is metres, and 223 of them give a ratio of
#    exactly 100.000 -- so fg19 DBH is CENTIMETRES. bytag is millimetres as
#    labelled (median 2.25 cm, against 2.30 for fg19 and 2.40 for the shapefile
#    copy). Only 8 stems in 8,014 exceed a defensible 100 cm; dividing each by 10
#    until it falls in range lands every one between 19.6 and 43.2 cm. This single
#    rule replaces the previous fx() stack of species-specific thresholds, which
#    was also over-correcting (fg19 basal area 135.48 m2 here vs 133.83 under fx).
#
# STEMS WITHOUT COORDINATES ARE KEPT. 51 stems (0.64%, 1.24% of basal area) have
# no usable PX/PY. They contribute flux regardless of whether we can place them,
# so they are retained and flagged `located = FALSE`. They are assigned to the
# stand: the notch was never surveyed, so no stem record can have originated
# there, and an unlocated record therefore comes from censused ground.
# ==============================================================================
suppressPackageStartupMessages({library(dplyr)})
source("code/lib/rev_geometry.R")

RAW <- "data/raw/inventory"
OUT <- "outputs/tables/inventory_stems.csv"
MAX_DBH_CM <- 100      # largest defensible stem at this site; see header note 3

# ---- species code -> name, from the workflow taxonomy ------------------------
e <- new.env(); load("data/processed/integrated/rf_workflow_input_data_with_2023.RData", envir = e)
lut <- e$rf_workflow_data$PLACEHOLDER_INVENTORY %>%
  distinct(species_code, species) %>% filter(!is.na(species_code)) %>%
  mutate(species_code = trimws(species_code))

# ---- decimal-typo repair -----------------------------------------------------
#' Shift the decimal left until the diameter is defensible. Returns the repaired
#' value and the number of shifts applied, so the repair is auditable.
repair_dbh <- function(cm) {
  shifts <- rep(0L, length(cm)); ok <- is.finite(cm)
  while (any(ok & cm > MAX_DBH_CM, na.rm = TRUE)) {
    i <- ok & cm > MAX_DBH_CM & !is.na(cm)
    cm[i] <- cm[i] / 10; shifts[i] <- shifts[i] + 1L
  }
  list(cm = cm, shifts = shifts)
}

num <- function(x) suppressWarnings(as.numeric(x))

fg19 <- read.csv(file.path(RAW, "ForestGEO_data2021UPDATE_6_21_DW_2019.csv"), stringsAsFactors = FALSE)
btag <- read.csv(file.path(RAW, "ForestGEO_data2021UPDATE_6_21_DW_bytag.csv"), stringsAsFactors = FALSE)

A <- fg19 %>% transmute(
  source = "fg19", tag = Tag,
  stem   = if ("Stem_Tag" %in% names(fg19)) as.character(Stem_Tag) else NA_character_,
  species_code = trimws(Species_Code),
  dbh_cm_raw = num(DBH),                     # centimetres
  PX = num(PX), PY = num(PY)) %>%
  distinct(tag, stem, .keep_all = TRUE)

B <- btag %>% transmute(
  source = "bytag", tag = Tag.ID, stem = NA_character_,
  species_code = trimws(Species..4.char.),
  dbh_cm_raw = num(Diam..mm.) / 10,          # millimetres -> centimetres
  PX = num(PX), PY = num(PY)) %>%
  distinct(tag, .keep_all = TRUE)

INV <- bind_rows(A, B) %>% filter(is.finite(dbh_cm_raw), dbh_cm_raw > 0)

rp <- repair_dbh(INV$dbh_cm_raw)
INV <- INV %>% mutate(
  dbh_cm       = rp$cm,
  dbh_shifts   = rp$shifts,
  dbh_m        = dbh_cm / 100,
  basal_area_m2 = pi * (dbh_m / 2)^2,
  species      = lut$species[match(species_code, lut$species_code)],
  located      = is.finite(PX) & is.finite(PY) & PX != 0 & PY != 0,
  in_notch     = located & in_notch(PX, PY),
  in_square    = located & in_square(PX, PY),
  # budget set: censused ground. Unlocated stems join it (see header).
  in_stand     = ifelse(located, in_square & !in_notch(PX, PY), TRUE))

stopifnot(max(INV$dbh_cm) <= MAX_DBH_CM, all(INV$dbh_m > 0))

# ---- verify the uncensused quadrats against the raw census ------------------
# UNCENSUSED_QUADRATS in rev_geometry.R is a constant so the geometry is stable
# and auditable, but it must keep matching the raw tables. Re-derive it here from
# the quadrat identifiers and stop on any drift: a quadrat appearing in the census
# that we are excluding would silently delete real ground.
local({
  qn <- function(v) suppressWarnings(as.integer(v))
  present <- sort(unique(na.omit(c(qn(fg19$Quadrat), qn(btag$X...Quadrat)))))
  grid <- expand.grid(row = 0:9, col = 0:9)
  grid$q <- grid$row*100 + grid$col
  got <- sort(grid$q[!(grid$q %in% present)])
  cat(sprintf("\n-- uncensused quadrats --\n  %d of 100 quadrats present in the raw census; %d absent\n",
              sum(grid$q %in% present), length(got)))
  cat(sprintf("  absent: %s\n", paste(got, collapse = ", ")))
  if (!identical(as.integer(got), as.integer(sort(UNCENSUSED_QUADRATS)))) {
    cat(sprintf("  CONSTANT: %s\n", paste(sort(UNCENSUSED_QUADRATS), collapse = ", ")))
    stop("uncensused quadrats have drifted from UNCENSUSED_QUADRATS in rev_geometry.R")
  }
  cat("  matches UNCENSUSED_QUADRATS in rev_geometry.R\n")
  ins <- sum(INV$located & in_gap(INV$PX, INV$PY))
  cat(sprintf("  located stems falling inside them: %d (excluded with the quadrats)\n", ins))
})

# ---- report ------------------------------------------------------------------
cat("================ INVENTORY REBUILD ================\n")
cat(sprintf("stem records: %d   (fg19 %d, bytag %d)\n", nrow(INV),
            sum(INV$source == "fg19"), sum(INV$source == "bytag")))
cat(sprintf("total basal area: %.2f m2\n\n", sum(INV$basal_area_m2)))

cat("-- decimal-typo repair --\n")
cat(sprintf("repaired %d of %d stems (%.3f%%); max DBH after repair %.1f cm\n",
            sum(INV$dbh_shifts > 0), nrow(INV), 100*mean(INV$dbh_shifts > 0), max(INV$dbh_cm)))
print(INV %>% filter(dbh_shifts > 0) %>%
      transmute(source, tag, species, raw_cm = dbh_cm_raw, shifts = dbh_shifts,
                fixed_cm = round(dbh_cm, 2)) %>% arrange(-raw_cm))

cat("\n-- coordinates --\n")
f <- function(sel, lab) cat(sprintf("  %-34s %5d stems (%6.2f%%)  %8.3f m2 BA (%6.2f%%)\n",
  lab, sum(sel), 100*mean(sel), sum(INV$basal_area_m2[sel]),
  100*sum(INV$basal_area_m2[sel])/sum(INV$basal_area_m2)))
f(INV$located,  "located")
f(!INV$located, "no coordinate (kept, in stand)")
f(INV$located & INV$in_notch, "in the unsurveyed notch (excluded)")
f(INV$located & !INV$in_square, "outside the square (excluded)")
f(INV$in_stand, "IN STAND -> budget set")

cat("\n-- stand totals --\n")
S <- INV %>% filter(in_stand)
cat(sprintf("  stems           %d\n", nrow(S)))
cat(sprintf("  basal area      %.2f m2  ->  %.1f m2 ha-1 over %.2f ha\n",
            sum(S$basal_area_m2), sum(S$basal_area_m2)/(STAND_AREA_M2/1e4), STAND_AREA_M2/1e4))
cat(sprintf("  stem density    %.0f ha-1\n", nrow(S)/(STAND_AREA_M2/1e4)))

cat("\n-- composition (stand, top 10 by basal area) --\n")
print(S %>% group_by(species) %>%
      summarise(stems = n(), BA = round(sum(basal_area_m2), 2), .groups = "drop") %>%
      mutate(pct_BA = round(100*BA/sum(BA), 2)) %>% arrange(-BA) %>% head(10) %>% as.data.frame())

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(INV, OUT, row.names = FALSE)
cat(sprintf("\nwritten: %s  (%d rows)\n", OUT, nrow(INV)))

# ---- provenance record -------------------------------------------------------
# The built inventory lives under outputs/ and is therefore gitignored (the repo
# tracks code, not regenerable products). Traceability comes from this record:
# which raw files were read, their checksums, and exactly which rows were altered.
PROV <- "outputs/audit/inventory_provenance.txt"
dir.create(dirname(PROV), showWarnings = FALSE, recursive = TRUE)
srcs <- file.path(RAW, c("ForestGEO_data2021UPDATE_6_21_DW_2019.csv",
                         "ForestGEO_data2021UPDATE_6_21_DW_bytag.csv"))
md5 <- tryCatch(tools::md5sum(srcs), error = function(e) setNames(rep(NA, length(srcs)), srcs))
con <- file(PROV, "w")
wr <- function(...) cat(sprintf(...), file = con)
wr("INVENTORY PROVENANCE\ngenerated by code/01_import/rev_inventory_build.R\n")
wr("built: %s\nR: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R.version.string)
wr("RAW INPUTS (md5)\n")
for (s in srcs) wr("  %-58s %s\n", basename(s), md5[[s]])
wr("\nGEOMETRY\n  stand = %d x %d m minus a %d-cell staircase gap (%d m cells)\n",
   PLOT_SIDE_M, PLOT_SIDE_M, GAP_AREA_M2/GAP_CELL_M^2, GAP_CELL_M)
wr("  gap by PY band (PX at which it starts): %s\n",
   paste(ifelse(is.na(GAP_PX_MIN_BY_PY), "-", GAP_PX_MIN_BY_PY), collapse = " "))
wr("  stand area %d m2 (%.2f ha); nominal %d m2; gap %d m2\n",
   STAND_AREA_M2, STAND_AREA_M2/1e4, PLOT_AREA_M2, GAP_AREA_M2)
wr("\nRECORD COUNTS\n  fg19 %d | bytag %d | total %d\n",
   sum(INV$source == "fg19"), sum(INV$source == "bytag"), nrow(INV))
wr("  de-duplicated WITHIN source only (tag namespaces are disjoint)\n")
wr("\nDBH\n  units: fg19 cm, bytag mm; repaired if > %d cm by /10 until in range\n", MAX_DBH_CM)
wr("  repaired %d stems:\n", sum(INV$dbh_shifts > 0))
for (i in which(INV$dbh_shifts > 0))
  wr("    %-6s tag %-6s %-24s %9.1f -> %6.2f cm (/10 x%d)\n",
     INV$source[i], INV$tag[i], ifelse(is.na(INV$species[i]), "<unknown>", INV$species[i]),
     INV$dbh_cm_raw[i], INV$dbh_cm[i], INV$dbh_shifts[i])
wr("\nCOORDINATES\n  located %d | unlocated %d (kept, assigned to stand)\n",
   sum(INV$located), sum(!INV$located))
wr("  excluded: %d in the unsurveyed gap, %d outside square\n",
   sum(INV$located & INV$in_notch), sum(INV$located & !INV$in_square))
wr("\nBUDGET SET (in_stand)\n  %d stems | %.2f m2 basal area | %.1f m2 ha-1 | %.0f stems ha-1\n",
   nrow(S), sum(S$basal_area_m2), sum(S$basal_area_m2)/(STAND_AREA_M2/1e4),
   nrow(S)/(STAND_AREA_M2/1e4))
close(con)
cat(sprintf("written: %s\n", PROV))
