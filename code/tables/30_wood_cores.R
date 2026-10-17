#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 30_wood_cores.R -- wood_cores.csv (data/tables/)
# ------------------------------------------------------------------------------
# One row per core section: tree x core x radial position (outer, middle, inner),
# from the 2022 core-sectioning sheet of the 2021 survey cores.
#
# SOURCE SHEET. data/raw/tree_cores/Tree Core Sectioning Data.csv (wide: one row
# per core, 158 rows, outer + middle + inner sections) -- the sheet
# code/01_import/02_process_wood_cores.R reads. The other sheet,
# Tree_Core_Sectioning_Data.csv (long, 316 rows), is a reshape of the outer and
# inner sections of the same 158 cores (no middle sections) and is not used here:
#   - its core_type is swapped: rows labelled Inner hold the outer section (their
#     Ground Sample ID reads "<tree> Outer") and vice versa, for all 316 rows;
#   - two fresh masses differ, both typos corrected in the wide sheet only
#     (WAX1 outer 0.3349 -> 0.8349 g; AB6' outer 6.3458 -> 0.6346 g; the long
#     sheet's Percent_Change matches its own typo).
#   (code/lib/build_phyloseq.R used to read it, merged on core_type; it now reads the tables.)
#
# CORES. Three trees were cored twice (WA14, WA102, YB9: two rows each in the
# sheet). core_id = <tree_id>-<n>, n = 1, 2 in sheet order; key core_id + position. Analyses average a tree's cores
# (mean with na.rm; code/lib/tables_read.R survey_2021_trees()).
#
# VALUES as in 02_process_wood_cores.R:
#   section volume = pi * (0.515 cm / 2)^2 * length (5.15 mm corer); middle
#     lengths are recorded in mm (converted); length 0 or missing -> no volume.
#   density_g_cm3      = dry mass / volume
#   moisture_fresh_pct = (fresh - dry) / fresh * 100
#   moisture_dry_pct   = (fresh - dry) / dry * 100
#   filled (gap-filling rule): a section of one of the seven trees listed in
#     02_process_wood_cores.R (AB6', H104, HX2, RM1, rms1, WA14, YB15 -- matched on
#     the sheet's Tree Core ID) is filled when dry mass
#     > fresh mass, moisture_dry > 300 %, or density outside 0.1-2.0 g cm-3. Its
#     density and both moistures are then replaced by the median of that position
#     over the other cores of the same species (sheet species code), problem
#     sections excluded. Masses and length stay as recorded.
# Rows: every section with any recorded length, mass or derived value.
#
# (The check against merged_tree_dataset_final.csv was removed 2026-10-16 when that
# file was retired.)
#
# Writes data/tables/wood_cores.csv, outputs/audit/wood_cores_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
OUT <- "data/tables"
SRC <- "data/raw/tree_cores/Tree Core Sectioning Data.csv"
CORE_RADIUS_CM <- (5.15 / 2) / 10
PROBLEM_TREES <- c("AB6'", "H104", "HX2", "RM1", "rms1", "WA14", "YB15")   # 02_process_wood_cores.R

raw <- read_csv(SRC, show_col_types = FALSE, name_repair = "minimal")
stopifnot(nrow(raw) == 158)
num <- function(x) suppressWarnings(as.numeric(x))
col <- function(nm) raw[[which(names(raw) == nm)[1]]]
raw$core <- ave(seq_len(nrow(raw)), raw$`Tree Core ID`, FUN = seq_along)

pos_cols <- list(
  outer  = c(len = "Outer Section Length (cm)", fresh = "OUTER SECTION FRESH MASS (g)",
             dry = "Post-Freeze Dry OUTER SECTION Sample Mass (g)", note = "OUTER SECTION Processing Issues/Notes"),
  middle = c(len = "Middle Section Length (mm)", fresh = "MIDDLE SECTION FRESH MASS (g)",
             dry = "Post-Freeze Dry MIDDLE SECTION Sample Mass (g)", note = "MIDDLE SECTION Processing Issues/Notes"),
  inner  = c(len = "Inner Section Length (cm)", fresh = "INNER SECTION FRESH MASS (g)",
             dry = "Post-Freeze Dry INNER SECTION Sample Mass (g)", note = "INNER SECTION Processing Issues/Notes"))
S <- bind_rows(lapply(names(pos_cols), function(p) {
  cc <- pos_cols[[p]]
  tibble(tag = raw$`Tree Core ID`, species = raw$`Tree Core Species Code`, core = raw$core,
         date = as.Date(raw$Date, "%m/%d/%y"), position = p,
         section_length_cm = num(col(cc[["len"]])) / ifelse(p == "middle", 10, 1),
         fresh_mass_g = num(col(cc[["fresh"]])), dry_mass_g = num(col(cc[["dry"]])),
         core_note = raw$`Tree Core Notes`, sec_note = col(cc[["note"]]))
}))

# ---- derived values (02_process_wood_cores.R) -----------------------------------
S <- S %>% mutate(
  volume = ifelse(!is.na(section_length_cm) & section_length_cm > 0, pi * CORE_RADIUS_CM^2 * section_length_cm, NA),
  m_fresh = ifelse(!is.na(fresh_mass_g) & !is.na(dry_mass_g) & fresh_mass_g > 0, (fresh_mass_g - dry_mass_g) / fresh_mass_g * 100, NA),
  m_dry = ifelse(!is.na(fresh_mass_g) & !is.na(dry_mass_g) & dry_mass_g > 0, (fresh_mass_g - dry_mass_g) / dry_mass_g * 100, NA),
  dens = ifelse(!is.na(dry_mass_g) & !is.na(volume) & volume > 0, dry_mass_g / volume, NA),
  filled = tag %in% PROBLEM_TREES &
    ((!is.na(fresh_mass_g) & !is.na(dry_mass_g) & dry_mass_g > fresh_mass_g) |
       (!is.na(m_dry) & m_dry > 300) | (!is.na(dens) & (dens > 2.0 | dens < 0.1))))
med <- S %>% group_by(species, position) %>%
  summarise(md_fresh = median(m_fresh[!filled], na.rm = TRUE), md_dry = median(m_dry[!filled], na.rm = TRUE),
            md_dens = median(dens[!filled], na.rm = TRUE), .groups = "drop")
S <- S %>% left_join(med, by = c("species", "position")) %>% mutate(
  density_g_cm3 = ifelse(filled, md_dens, dens),
  moisture_fresh_pct = ifelse(filled, md_fresh, m_fresh),
  moisture_dry_pct = ifelse(filled, md_dry, m_dry))

# ---- tree_id from tree_tags ------------------------------------------------------
TT <- read.csv(file.path(OUT, "tree_tags.csv"), colClasses = "character", na.strings = "") %>%
  filter(source_file == SRC, source_column == "Tree Core ID")
TT <- TT %>% distinct(tag = trimws(tag_recorded), tree_id, .keep_all = TRUE)   # spellings as read.csv keeps them
stopifnot(all(is.na(TT$species_recorded)), !anyDuplicated(TT$tag))
S$tree_id <- TT$tree_id[match(S$tag, TT$tag)]
if (anyNA(S$tree_id)) stop("no tree_tags entry for ", paste(unique(S$tag[is.na(S$tree_id)]), collapse = ", "))

clean_note <- function(...) { x <- trimws(c(...)); x <- unique(x[!is.na(x) & x != ""]); if (length(x)) paste(x, collapse = "; ") else NA_character_ }
WC <- S %>%
  filter(!is.na(section_length_cm) | !is.na(fresh_mass_g) | !is.na(dry_mass_g) | !is.na(density_g_cm3) |
           !is.na(moisture_fresh_pct)) %>%
  rowwise() %>% mutate(notes = clean_note(core_note, sec_note)) %>% ungroup() %>%
  transmute(core_id = paste0(tree_id, "-", core), tree_id, position, sectioning_date = date, section_length_cm, fresh_mass_g, dry_mass_g,
            density_g_cm3, moisture_fresh_pct, moisture_dry_pct, filled, notes) %>%
  arrange(tree_id, core_id, match(position, c("outer", "middle", "inner")))
stopifnot(!anyDuplicated(WC[c("core_id", "position")]))

sink(out_path("wood_cores_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("wood_cores: %d sections from %d cores of %d trees (%s)\n", nrow(WC), n_distinct(WC$core_id),
            n_distinct(WC$tree_id), paste(names(table(WC$position)), table(WC$position), collapse = ", ")))
cat(sprintf("trees cored twice: %s\n", paste(unique(WC$tree_id[grepl("-2$", WC$core_id)]), collapse = ", ")))
cat(sprintf("filled sections: %d\n", sum(WC$filled)))
print(as.data.frame(WC %>% filter(filled) %>% select(core_id, position, density_g_cm3, moisture_fresh_pct, moisture_dry_pct)), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("wood_cores_build.txt", "audit")), sep = "\n")

write_csv(WC, file.path(OUT, "wood_cores.csv"), na = "")
