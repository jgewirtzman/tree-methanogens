#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 21_ddpcr.R -- ddpcr.csv (data/tables/): one row per sample x gene x assay x plate
# ------------------------------------------------------------------------------
# Survey (2021 samples): the QX200 exports in data/raw/ddpcr/Data Files/Raw Data/,
# one file per ddPCR plate (run) and threshold ("loose" / "strict" in the file name;
# the paper uses loose throughout). A well is the extraction well of the same
# position on the extraction plate named in "Sample description 1" (as in
# 01_import/05_ddpcr_import_and_align.R). Exceptions handled here:
#   * the probe-mcrA plate-8 run also re-ran 11 wells of extraction plate 1
#     ("Plate_1_4EFGH" in ddPCR wells E06-H06 = extraction wells E4-H4;
#     "Plate_1_5ABCDEFG" in A07-G07 = A5-G5; 10 hold samples). Its strict export is
#     "..._Plate1Redos_..._strict"; its loose export has the plate-8 wells only.
#   * the probe-mcrA plate-1 file ("combined", loose only; written the day after
#     the plate-8 run) carries, for 9 of those 10 wells, the plate-8 re-run (same
#     accepted droplets). Those loose values go on the plate-8 run's row (one row per
#     sample x run); the 10th (YB15Mineral: different droplets) keeps a plate-1 row
#     (loose) and a plate-8 row (strict only).
# Standard wells and wells whose extraction well holds no sample (blank guide wells)
# are not samples and are left out. qc_pass = FALSE where accepted droplets < 5,000
# or the software made no call, or the extraction guide marks the well blank (three
# wells that nonetheless have a 16S library; samples.well_marked_blank) -- the
# import's filters; those rows are kept.
# Felled oak: data/raw/ddpcr/black_oak_mcrA.csv (one export of one probe-mcrA plate;
# its threshold is not recorded -- stored in the loose columns, as fig07 uses it).
#
# copies_per_g = code/lib/ddpcr_constants.R (Conc x 25/2.5 x elution 75 wood / 100
# soil x 800/250 / mass), loose Conc, extraction mass of samples.csv. Mass basis:
# freeze-dried, cryo-ground wood (dry); sieved, frozen soil (fresh) -- as stated in
# the paper. below_lod: loose Conc < 3 copies per 20 µL well (0.15 copies/µL).
# saturated: every accepted droplet positive at the loose threshold (QuantaSoft then
# reports 1,000,000 copies/µL; six 16S-bacteria / EvaGreen-mcrA wells).
#
# Checks (stop on failure): every analysed row of processed_ddpcr_data.csv is
# reproduced (Conc, droplets, positives); compartment medians equal
# 09_tables_stats/06_copies-per-gram.R's. (The check against the tree-level means of
# merged_tree_dataset_final.csv was removed 2026-10-16 when that file was retired.)
# Writes data/tables/ddpcr.csv, outputs/audit/ddpcr_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/ddpcr_constants.R")
OUT <- "data/tables"
LOD_COPIES_PER_WELL <- 3; MIN_DROPLETS <- 5000
SAMPLES <- read.csv(file.path(OUT, "samples.csv"), stringsAsFactors = FALSE)

# ==============================================================================
# 1. SURVEY RUNS
# ==============================================================================
files <- list.files("data/raw/ddpcr/Data Files/Raw Data", "\\.csv$", full.names = TRUE)
files <- files[!grepl("Icon|\\._|DS_Store|\\.tmp", basename(files))]
R <- bind_rows(lapply(files, function(f) {
  d <- read.csv(f, check.names = FALSE, colClasses = "character")
  tibble(file = basename(f), ddpcr_well = d$Well, sd1 = d$`Sample description 1`, sd2 = d$`Sample description 2`,
         target = d$Target, supermix = d$Supermix, conc = d$`Conc(copies/µL)`, accepted = d$`Accepted Droplets`,
         positives = d$Positives)
}))
R <- R %>%
  mutate(threshold = str_extract(file, "(loose|strict)(?=\\.csv$)"),
         stamp = str_extract(file, "[0-9]{8}_[0-9]{6}_[0-9]{3}"),
         gene = recode(tolower(target), mcra_probe = "mcrA", mcra = "mcrA", pmoa = "pmoA", mmox = "mmoX",
                       `16s_arc` = "16S_archaea", `16s_bact` = "16S_bacteria"),
         assay = ifelse(grepl("Probes", supermix), "probe", ifelse(grepl("EvaGreen", supermix), "EvaGreen", NA)),
         conc = suppressWarnings(as.numeric(conc)), accepted = as.integer(accepted), positives = as.integer(positives),
         row = substr(ddpcr_well, 1, 1), col = as.integer(substr(ddpcr_well, 2, 3)),
         # extraction plate and well; plate-1 re-runs carry plate, column and rows in sd1
         redo = str_match(sd1, "^Plate_([0-9])_([0-9]+)([A-H]+)$"),
         extraction_plate = as.integer(ifelse(!is.na(redo[, 2]), redo[, 2], str_extract(sd1, "[0-9]+$"))),
         extraction_well = ifelse(!is.na(redo[, 2]), paste0(row, redo[, 3]), paste0(row, col)),
         rerun = !is.na(redo[, 2]))
stopifnot(!anyNA(R$threshold), !anyNA(R$gene), !anyNA(R$assay),
          all(mapply(grepl, R$row[R$rerun], R$redo[R$rerun, 4])))
# one plate_id per run: the loose file's stem (the plate-8 strict file has another name)
runs <- R %>% distinct(file, stamp, threshold) %>% group_by(stamp) %>%
  summarise(plate_id = sub("_(loose|strict)\\.csv$", "", if (any(threshold == "loose")) file[threshold == "loose"][1] else file[1]),
            .groups = "drop")
R <- R %>% left_join(runs, by = "stamp")

R$std <- grepl("std", R$sd2, ignore.case = TRUE)
hit <- match(paste(R$extraction_plate, R$extraction_well),
             paste(SAMPLES$extraction_plate, SAMPLES$extraction_well)[SAMPLES$campaign == "2021_survey"])
R$sample_id <- SAMPLES$sample_id[SAMPLES$campaign == "2021_survey"][hit]
dropped <- R %>% filter(std | is.na(sample_id)) %>% mutate(why = ifelse(std, "standard", "no sample in the extraction well"))
R <- R %>% filter(!std, !is.na(sample_id))

wide <- function(d) d %>%
  select(sample_id, gene, assay, plate_id, ddpcr_well, rerun, threshold, accepted, positives, conc) %>%
  pivot_wider(names_from = threshold, values_from = c(accepted, positives, conc))
W <- wide(R)
for (cn in c("accepted_loose", "accepted_strict", "positives_strict", "conc_strict")) if (!cn %in% names(W)) W[[cn]] <- NA
# loose and strict exports of one run hold the same droplets
same <- with(W, is.na(accepted_loose) | is.na(accepted_strict) | accepted_loose == accepted_strict)
stopifnot(all(same))
# plate-1 "combined" loose values that are the plate-8 re-run -> onto the re-run row
W$loose_from <- NA_character_; drop <- integer(0)
for (i in which(W$rerun)) {
  j <- which(W$sample_id == W$sample_id[i] & W$gene == W$gene[i] & W$assay == W$assay[i] &
             grepl("combined", W$plate_id) & W$accepted_loose %in% W$accepted_strict[i])
  if (length(j) == 1) {
    for (cn in c("accepted_loose", "positives_loose", "conc_loose")) W[[cn]][i] <- W[[cn]][j]
    W$loose_from[i] <- W$plate_id[j]; drop <- c(drop, j)
  }
}
W <- W[-drop, ]

# ==============================================================================
# 2. FELLED OAK
# ==============================================================================
B <- read.csv("data/raw/ddpcr/black_oak_mcrA.csv", check.names = FALSE, colClasses = "character",
              fileEncoding = "UTF-8-BOM") %>% filter(`Sample ID` != "")
stopifnot(all(B$Target == "mcrA_probe"), all(gsub("-", "", B$`Sample ID`) %in% SAMPLES$sample_id))
WB <- tibble(sample_id = gsub("-", "", B$`Sample ID`), gene = "mcrA", assay = "probe", plate_id = "black_oak_mcrA",
             ddpcr_well = B$Well, rerun = FALSE, loose_from = NA_character_, accepted_loose = as.integer(B$`Accepted Droplets`),
             positives_loose = as.integer(B$Positives), conc_loose = as.numeric(B$`Conc(copies/µL)`))

# ==============================================================================
# 3. ASSEMBLE
# ==============================================================================
S <- SAMPLES %>% select(sample_id, campaign, material, position, extraction_mass_mg, well_marked_blank)
DD <- bind_rows(W, WB) %>% left_join(S, by = "sample_id") %>%
  mutate(accepted_droplets = coalesce(accepted_loose, accepted_strict),
         qc_pass = !is.na(accepted_droplets) & accepted_droplets >= MIN_DROPLETS & (!is.na(conc_loose) | !is.na(conc_strict)) &
                   !well_marked_blank %in% TRUE,
         qc_reason = case_when(well_marked_blank %in% TRUE ~ "extraction guide marks the well blank",
                               is.na(conc_loose) & is.na(conc_strict) ~ "no call",
                               accepted_droplets < MIN_DROPLETS ~ "accepted droplets < 5000",
                               TRUE ~ NA_character_),
         mat = ifelse(material == "soil", "Soil", "Wood"),
         copies_per_g = ifelse(qc_pass, ddpcr_copies_per_g(conc_loose, extraction_mass_mg, mat), NA),
         mass_basis = case_when(material == "soil" ~ "fresh", material == "wood" ~ "dry", TRUE ~ NA_character_),
         below_lod = ifelse(qc_pass, conc_loose * DDPCR_WELL_UL < LOD_COPIES_PER_WELL, NA))
stopifnot(all(DD$material %in% c("wood", "soil")), !anyNA(DD$mass_basis))
DDPCR <- DD %>%
  transmute(sample_id, gene, assay, plate_id, ddpcr_well, accepted_droplets,
            positive_droplets = positives_loose, copies_per_uL_rxn = conc_loose,
            positive_droplets_strict = positives_strict, copies_per_uL_rxn_strict = conc_strict,
            copies_per_g, mass_basis, below_lod, saturated = positives_loose == accepted_droplets, qc_pass, qc_reason,
            note = case_when(rerun & !is.na(loose_from) ~ "plate-1 well re-run on the plate-8 run; loose values from the plate-1 combined file (same droplets)",
                             rerun ~ "plate-1 well re-run on the plate-8 run; strict export only (the plate-1 file's loose value has other droplets)",
                             plate_id == "black_oak_mcrA" ~ "felled oak; single export, threshold not recorded",
                             TRUE ~ NA_character_)) %>%
  arrange(gene, assay, sample_id, plate_id)
k <- paste(DDPCR$sample_id, DDPCR$gene, DDPCR$assay, DDPCR$plate_id)
if (anyDuplicated(k)) stop("ddpcr key not unique: ", paste(head(k[duplicated(k)]), collapse = ", "))
write_csv(DDPCR, file.path(OUT, "ddpcr.csv"), na = "")

# ==============================================================================
# 4. RECONCILIATION
# ==============================================================================
# (a) every analysed row of processed_ddpcr_data.csv (the source of everything downstream)
P <- read.csv("data/processed/molecular/processed_ddpcr_data.csv", check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(extraction_plate = as.integer(sub("Plate ", "", `Extraction.Plate.ID`)), extraction_well = `Plate.Position`,
            gene = recode(Target, mcra_probe = "mcrA", mcra = "mcrA", pmoa = "pmoA", mmox = "mmoX",
                          `16s_arc` = "16S_archaea", `16s_bact` = "16S_bacteria"),
            assay = ifelse(Target %in% c("mcra_probe"), "probe", "EvaGreen"), threshold = analysis_type,
            conc = `Conc.copies.µL.`, accepted = Accepted.Droplets, positives = Positives,
            mass = `Sample.Mass.Added.to.Tube..mg.`, material) %>%
  mutate(sample_id = SAMPLES$sample_id[match(paste(extraction_plate, extraction_well),
                                             paste(SAMPLES$extraction_plate, SAMPLES$extraction_well))])
stopifnot(!anyNA(P$sample_id))
L <- bind_rows(DD %>% filter(campaign == "2021_survey") %>% transmute(sample_id, gene, assay, threshold = "loose", conc = conc_loose, positives = positives_loose, accepted = accepted_droplets, qc_pass),
               DD %>% filter(!rerun, campaign == "2021_survey") %>% transmute(sample_id, gene, assay, threshold = "strict", conc = conc_strict, positives = positives_strict, accepted = accepted_droplets, qc_pass)) %>%
  filter(!is.na(positives))
cmpP <- P %>% left_join(L, by = c("sample_id", "gene", "assay", "threshold"), suffix = c("", ".t"))
okP <- with(cmpP, !is.na(positives.t) & abs(conc - conc.t) < 1e-6 * pmax(1, conc) & accepted == accepted.t & positives == positives.t)
if (!all(okP)) { print(head(cmpP[!okP, ])); stop(sum(!okP), " processed_ddpcr_data rows not reproduced") }
# rows the import drops are exactly the qc failures (plus the plate-1 re-runs)
extraP <- L %>% anti_join(P, by = c("sample_id", "gene", "assay", "threshold"))
if (any(extraP$qc_pass)) { print(extraP %>% filter(qc_pass)); stop("rows with qc_pass missing from processed_ddpcr_data.csv") }
massP <- P %>% distinct(sample_id, mass) %>% left_join(S, by = "sample_id")
stopifnot(all(abs(massP$mass - massP$extraction_mass_mg) < 1e-9 | (is.na(massP$mass) & is.na(massP$extraction_mass_mg))))

# (c) 06_copies-per-gram.R compartment medians (probe mcrA, loose, real mass)
C6 <- DD %>% filter(gene == "mcrA", assay == "probe", !is.na(conc_loose), campaign == "2021_survey", qc_pass, !is.na(extraction_mass_mg)) %>%
  group_by(position) %>% summarise(n = n(), median_copies_per_g = median(copies_per_g), .groups = "drop")
f6 <- "outputs/data/copies_per_g_compartment.csv"
if (file.exists(f6)) {
  O6 <- read.csv(f6) %>% mutate(position = tolower(sub("Soil_", "", recode(compartment, Heartwood = "inner", Sapwood = "outer"))))
  c6 <- left_join(C6, O6, by = "position")
  stopifnot(all(c6$n.x == c6$n.y), all(abs(c6$median_copies_per_g - c6$median_real) <= 1e-6 * c6$median_real))
}

# ==============================================================================
# 5. AUDIT
# ==============================================================================
sink(out_path("ddpcr_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("ddpcr %d rows (%d samples); survey %d, felled oak %d\n", nrow(DDPCR), n_distinct(DDPCR$sample_id),
            sum(DDPCR$plate_id != "black_oak_mcrA"), sum(DDPCR$plate_id == "black_oak_mcrA")))
print(as.data.frame(DDPCR %>% count(gene, assay, mass_basis, qc_pass)), row.names = FALSE)
cat(sprintf("\nRaw survey wells: %d (loose + strict exports); left out: %d standards, %d wells with no sample in the extraction well\n",
            nrow(R) + nrow(dropped), sum(dropped$why == "standard"), sum(dropped$why != "standard")))
print(as.data.frame(dropped %>% filter(why != "standard") %>% count(gene, assay, threshold)), row.names = FALSE)
cat("\nPlate-1 probe-mcrA wells re-run on the plate-8 run:\n")
print(as.data.frame(DDPCR %>% filter(grepl("re-run", note)) %>% select(sample_id, ddpcr_well, accepted_droplets, positive_droplets, positive_droplets_strict, copies_per_uL_rxn, copies_per_uL_rxn_strict, note)), row.names = FALSE)
cat("\nqc failures (kept, qc_pass = FALSE):\n")
print(as.data.frame(DDPCR %>% filter(!qc_pass) %>% select(sample_id, gene, assay, plate_id, ddpcr_well, accepted_droplets, qc_reason)), row.names = FALSE)
cat(sprintf("\nRows with loose values only: %d; strict only: %d\n", sum(!is.na(DDPCR$copies_per_uL_rxn) & is.na(DDPCR$copies_per_uL_rxn_strict)),
            sum(is.na(DDPCR$copies_per_uL_rxn) & !is.na(DDPCR$copies_per_uL_rxn_strict))))
cat(sprintf("Rows without copies_per_g although qc_pass (no extraction mass): %d\n", sum(DDPCR$qc_pass & !is.na(DDPCR$copies_per_uL_rxn) & is.na(DDPCR$copies_per_g))))
cat("\nRECONCILIATION\n")
cat(sprintf("  processed_ddpcr_data.csv: %d of %d rows reproduced exactly (Conc, accepted, positives; masses identical)\n", sum(okP), nrow(P)))
cat("  06_copies-per-gram.R medians (probe mcrA, loose):\n"); print(as.data.frame(C6), row.names = FALSE)
cat(sprintf("  %s\n", if (file.exists(f6)) "match outputs/data/copies_per_g_compartment.csv" else "outputs/data/copies_per_g_compartment.csv not found: not compared"))
cat("\nbelow_lod (qc_pass rows) by gene x assay x position:\n")
print(as.data.frame(DDPCR %>% filter(qc_pass, !is.na(below_lod)) %>% left_join(S, by = "sample_id") %>%
                      group_by(gene, assay, position) %>% summarise(n = n(), below_lod = sum(below_lod), .groups = "drop")), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("ddpcr_build.txt", "audit")), sep = "\n")
