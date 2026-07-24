#!/usr/bin/env Rscript
# ==============================================================================
# rev_stat_campaign_counts.R  —  REPRODUCIBLE campaign / sample-size accounting
# ------------------------------------------------------------------------------
# Generates every N in the study-design table (code/revision/notes/study_design_
# table.md) directly from the archived data, so the counts are reproducible and
# there is a single, documented definition for each.
#
# FLUX-COUNTING RULE (definitive): count every chamber deployment that produced
# a CH4 flux value (non-NA best.flux). **No quality filtering** — QC exclusions
# are a separate, downstream analysis concern, not part of the totals. The same
# rule is applied identically to every campaign (stem and soil).
#
# Run from repo root:  Rscript code/revision/rev_stat_campaign_counts.R
# Writes: outputs/revision/campaign_counts.csv  and  campaign_counts.txt
# ==============================================================================
suppressMessages({library(dplyr)})
options(warn = -1, stringsAsFactors = FALSE)
nn  <- function(x) sum(!is.na(x))                 # n non-NA
uq  <- function(x) length(unique(x[!is.na(x) & x != ""]))
P   <- function(...) cat(sprintf(...), "\n")
rows <- list()  # accumulate table rows

flux_dir <- "data/processed/flux"
comp     <- "data/compiled"

# ------------------------------------------------------------------ 1. FLUX ----
# Stem — three campaigns, each in its own file (different structures by design).
srt <- read.csv(file.path(flux_dir, "semirigid_tree_final_complete_dataset.csv"))       # monthly 2020-21
ch4 <- read.csv(file.path(flux_dir, "CH4_best_flux_lgr_results.csv"))                    # height 2021 (flux)
aux <- read.csv(file.path(flux_dir, "goflux_auxfile.csv"))                               # height 2021 (metadata)
y23 <- read.csv(file.path(flux_dir, "methanogen_tree_flux_complete_dataset.csv"))        # cross-species 2023
srs <- read.csv(file.path(flux_dir, "semirigid_tree_final_complete_dataset_soil.csv"))   # monthly soil

n_monthly <- nn(srt$CH4_best.flux.x)
n_height  <- nn(ch4$best.flux)
n_2023    <- nn(y23$CH4_best.flux)
n_soil    <- nn(srs$CH4_best.flux)
stem_total <- n_monthly + n_height + n_2023

# trees / species per campaign (keep-all: over rows with a flux value)
hd <- merge(ch4, aux[, c("UniqueID", "tree_id", "species", "plot")], by = "UniqueID", all.x = TRUE)
hd <- hd[!is.na(hd$best.flux), ]
mon_trees <- uq(srt$Plot.Tag[!is.na(srt$CH4_best.flux.x)])
h_trees   <- length(unique(paste(hd$plot, hd$tree_id)))
h_spp     <- uq(hd$species)
t23_ok    <- !is.na(y23$CH4_best.flux)
t23_trees <- uq(y23$Tree.Tag[t23_ok]);  t23_spp <- uq(y23$Species.Code[t23_ok])

# ---------------------------------------------------------------- 2. TREES -----
# Cross-campaign unique-tree union. Tag systems differ (monthly & 2023 numeric;
# 2021 uses plot/quadrat codes), so we normalise via the ID crosswalk where it
# maps, then union. Reported as a reconciled union + the naive union for context.
map_f <- "data/processed/tree_data/tree_id_comprehensive_mapping.csv"
norm  <- function(x) toupper(trimws(as.character(x)))
mon_id <- norm(srt$Plot.Tag[!is.na(srt$CH4_best.flux.x)])
h_id   <- norm(paste(hd$plot, hd$tree_id))
t23_id <- norm(y23$Tree.Tag[t23_ok])
naive_union <- length(unique(c(mon_id, h_id, t23_id)))
recon_union <- NA
if (file.exists(map_f)) {
  M <- read.csv(map_f, check.names = FALSE)
  vcols <- grep("^name_in_|^variant_|^name_variants$|primary_id|Tree_ID_normalized", names(M), value = TRUE)
  lut <- c()  # variant -> primary
  prim <- if ("primary_id" %in% names(M)) M$primary_id else M$Tree_ID_normalized
  for (cc in vcols) { v <- norm(M[[cc]]); ok <- v != "" & !is.na(v); lut[v[ok]] <- prim[ok] }
  canon <- function(ids) ifelse(ids %in% names(lut), lut[ids], ids)
  recon_union <- length(unique(canon(c(mon_id, h_id, t23_id))))
}

# ---------------------------------------------------------------- 3. ddPCR -----
dd <- read.csv(file.path(comp, "ddpcr_gene_abundances.csv"))
dds <- dd[!duplicated(dd$sample_id), ]
dd_wood <- sum(dds$material == "Wood");  dd_soil <- sum(dds$material == "Soil")
dd_total <- nrow(dds)
bo_mcra <- tryCatch(nrow(read.csv("data/raw/ddpcr/black_oak_mcrA.csv", check.names = FALSE)), error = function(e) NA)

# ----------------------------------------------------------------- 4. 16S ------
s16 <- read.csv(file.path(comp, "sample_metadata_16S.csv"), check.names = FALSE)
s16_wood <- sum(s16$Material == "Wood"); s16_soil <- sum(s16$Material == "Soil")
s16_oak  <- sum(s16$Material == "QUVE"); s16_total <- nrow(s16)

# -------------------------------------------------- 5. INTERNAL GAS / ISOTOPES -
n_gas <- tryCatch(nrow(read.csv("data/processed/internal_gas/sample_data_only.csv", check.names = FALSE)), error = function(e) NA)
n_iso <- tryCatch(nrow(read.csv("data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv", check.names = FALSE)), error = function(e) NA)

# ------------------------------------------------------------- 6. BLACK OAK ----
bo_flux <- tryCatch({b <- read.csv("data/compiled/black_oak_experiment.csv"); nrow(b)}, error = function(e) NA)

# -------------------------------------------------------------- 7. INVENTORY ---
inv <- read.csv(file.path(comp, "forest_inventory.csv"))
inv_live <- sum(inv$status == "LI", na.rm = TRUE)
inv_spp  <- uq(inv$species_code)
inv_area_ha <- round(diff(range(inv$x_m, na.rm = TRUE)) * diff(range(inv$y_m, na.rm = TRUE)) / 1e4, 2)

# ------------------------------------------------------------------- REPORT ----
sink("outputs/revision/campaign_counts.txt")
cat("CAMPAIGN / SAMPLE-SIZE ACCOUNTING  (rule: keep all fluxes, non-NA best.flux, no QC)\n")
cat("Generated by code/revision/rev_stat_campaign_counts.R\n")
cat(strrep("=", 78), "\n\n")
cat("FLUX (stem — count = deployments with a CH4 value):\n")
P("  Monthly survey  2020-21 : %4d stem fluxes | %d trees",           n_monthly, mon_trees)
P("  Height survey   2021    : %4d stem fluxes | %d trees | %d species (3 heights: 50/125/200 cm)", n_height, h_trees, h_spp)
P("  Cross-species   2023    : %4d stem fluxes | %d trees | %d species (breast height)", n_2023, t23_trees, t23_spp)
P("  ------------------------------------------------")
P("  STEM TOTAL              : %4d", stem_total)
P("  Monthly SOIL flux       : %4d", n_soil)
cat("\nUNIQUE TREES (flux, across campaigns):\n")
P("  naive union (normalised tags)      : %d", naive_union)
P("  reconciled via ID crosswalk        : %s", ifelse(is.na(recon_union), "n/a", recon_union))
cat("\nMOLECULAR / COMMUNITY (summer-2021 survey unless noted):\n")
P("  ddPCR unique samples    : %d  (wood %d / soil %d)", dd_total, dd_wood, dd_soil)
P("  16S samples             : %d  (wood %d / soil %d / black-oak %d)", s16_total, s16_wood, s16_soil, s16_oak)
P("  Internal gas / isotopes : %s / %s", n_gas, n_iso)
cat("\nFELLED BLACK OAK (Oct 2022, analysed separately — NOT in totals above):\n")
P("  stem fluxes (7 hts+seam): %s", bo_flux)
P("  ddPCR (mcrA only)       : %s measurements", bo_mcra)
P("  16S tissue samples      : %d", s16_oak)
cat("\nINVENTORY / UPSCALING:\n")
P("  live stems              : %d | %d species | extent-based area ~%s ha (CONFIRM plot area)", inv_live, inv_spp, inv_area_ha)
sink()

# machine-readable table
tab <- data.frame(
  campaign = c("Monthly survey","Height + molecular","Felled black oak","Cross-species survey","Inventory"),
  period   = c("Jun 2020-May 2021","Summer 2021","Oct 2022","Summer 2023","-"),
  trees    = c(mon_trees, h_trees, 1, t23_trees, inv_live),
  species  = c(NA, h_spp, 1, t23_spp, inv_spp),
  stem_flux= c(n_monthly, n_height, bo_flux, n_2023, NA),
  soil_flux= c(n_soil, NA, NA, NA, NA),
  ddpcr    = c(NA, dd_total, bo_mcra, NA, NA),
  s16      = c(NA, s16_wood + s16_soil, s16_oak, NA, NA),
  gas_iso  = c(NA, paste(n_gas, n_iso, sep="/"), "yes", NA, NA)
)
write.csv(tab, "outputs/revision/campaign_counts.csv", row.names = FALSE)
cat(readLines("outputs/revision/campaign_counts.txt"), sep = "\n")
cat("\n\nWrote outputs/revision/campaign_counts.{txt,csv}\n")
