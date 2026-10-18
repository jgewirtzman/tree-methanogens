#!/usr/bin/env Rscript
# ==============================================================================
# 02_campaign_counts.R  —  REPRODUCIBLE campaign / sample-size accounting
# ------------------------------------------------------------------------------
# Regenerates every N in the study-design table from the archived data, under
# one consistent, documented definition per quantity.
#
# FLUX RULE (definitive): count every chamber deployment that produced a CH4
# flux value (non-NA best.flux). **No quality filtering** — QC exclusions are a
# separate, downstream analysis concern. Applied identically to all campaigns.
#
# SOURCE: the canonical tables in data/tables/ (flux_stem, flux_soil, trees, ddpcr,
# samples, samples_16S, internal_gas, isotopes, wood_cores, inventory_stems), so every
# count is of the same rows every other script reads. Not the model's training export
# (outputs/data/flux_measurements_tree.csv): that is complete-cased on the model
# covariates and is the wrong source for how many measurements were taken.
# Trees are counted by tree_id (trees.csv: one physical tree; tags reused by two trees are
# disambiguated as <tag>_a / _b), species by trees.csv species_code; no crosswalk.
#
# NOTE: the 2021 survey includes black-oak individuals (tree_ids BO1/BO2/ROBOWSER,
# species QUVE) — these are SURVEY trees and stay in the survey totals. They are
# distinct from the single 2022 FELLED black oak, which is a separate section.
#
# Run from repo root:  Rscript code/09_tables_stats/02_campaign_counts.R
# Writes: outputs/audit/campaign_counts.txt  (+ outputs/data/campaign_counts.csv)
# ==============================================================================
suppressMessages(library(dplyr)); options(warn = -1, stringsAsFactors = FALSE)
nn <- function(x) sum(!is.na(x)); uq <- function(x) length(unique(x[!is.na(x) & x != ""]))
P  <- function(...) cat(sprintf(...), "\n")
source("code/lib/tables_read.R")
trees <- trees_table()
pos_of <- function(id) trees$landscape_position[match(id, trees$tree_id)]

# ============================================================== A. STEM FLUX ===
fs  <- read_table("flux_stem")                       # every deployment has a CH4 flux (NN)
mon <- fs[fs$campaign == "monthly_2020_2021", ]       # incl. recovered untagged and dead stems
hd  <- fs[fs$campaign == "2021_multiheight", ]
y23 <- fs[fs$campaign == "2023_cross_species", ]
oak <- fs[fs$campaign == "2022_felled_oak", ]
srs <- read_table("flux_soil")                       # monthly soil, all deployments (no QC)

n_monthly <- nn(mon$CH4_flux_nmol_m2_s); n_soil <- nn(srs$CH4_flux_nmol_m2_s)
n_height  <- nn(hd$CH4_flux_nmol_m2_s);  n_2023 <- nn(y23$CH4_flux_nmol_m2_s)
stem_total <- n_monthly + n_height + n_2023

h_lab <- ifelse(is.na(hd$height_cm), hd$height_label, as.character(hd$height_cm))
h_by <- sort(table(h_lab), decreasing = TRUE)
mon_trees <- uq(mon$tree_id)
mon_pos   <- table(pos_of(mon$tree_id))                # U/I/WD/WS (trees.csv)
mon_species <- tree_species(unique(mon$tree_id))
mon_spp <- uq(mon_species); mon_unres <- sum(is.na(mon_species))
h_trees <- uq(hd$tree_id); h_spp <- uq(tree_species(unique(hd$tree_id)))
t23_trees <- uq(y23$tree_id); t23_spp <- uq(tree_species(unique(y23$tree_id)))
union_trees <- uq(c(mon$tree_id, hd$tree_id, y23$tree_id))

# ============================================================= B. COVARIATES ===
# per deployment: air temperature in flux_stem; DBH per tree in trees.csv (measured values
# only; dbh_source says which record). Soil temperature and moisture are attached by
# code/lib/soil_fill.R (attach_soil) where an analysis needs them.
cov_line <- function(d) sprintf("air temperature %d of %d deployments; DBH for %d of %d trees",
  nn(d$air_temp_C), nrow(d), sum(!is.na(trees$dbh_cm[match(unique(d$tree_id), trees$tree_id)])), uq(d$tree_id))

# ============================================== C. MOLECULAR / COMMUNITY / GAS ==
smp <- read_table("samples")
dd  <- read_table("ddpcr"); dd$campaign <- smp$campaign[match(dd$sample_id, smp$sample_id)]
dd$target <- paste(dd$gene, dd$assay)
# runs that passed QC (ddpcr.qc_pass: >= 5000 accepted droplets, a call, well not marked
# blank); the 3 samples whose only runs fail are counted separately
dsv_all <- dd[dd$campaign == "2021_survey", ]
dsv <- dsv_all[dsv_all$qc_pass, ]
mat_of <- function(id) smp$material[match(id, smp$sample_id)]
dds <- unique(dsv$sample_id)
dd_surv <- length(dds); dd_w <- sum(mat_of(dds) == "wood"); dd_s <- sum(mat_of(dds) == "soil")
dd_failed <- uq(dsv_all$sample_id) - dd_surv
dd_genes <- sort(unique(dsv$target))
# one ddPCR reaction = a sample x target (re-runs of a sample count once); loose and strict
# thresholds are two analyses of one reaction
dd_by_target <- sort(table(dsv$target[!duplicated(paste(dsv$sample_id, dsv$target))]), decreasing = TRUE)
dd_reactions <- sum(dd_by_target)
dd_datapoints <- nn(dsv$copies_per_uL_rxn) + nn(dsv$copies_per_uL_rxn_strict)
s16 <- read_table("samples_16S")
s16$campaign <- smp$campaign[match(s16$sample_id, smp$sample_id)]
s16$position <- smp$position[match(s16$sample_id, smp$sample_id)]
s16_w <- sum(s16$campaign %in% "2021_survey" & s16$position %in% c("inner", "outer"))
s16_s <- sum(s16$campaign %in% "2021_survey" & s16$position %in% c("organic", "mineral"))
s16_none <- sum(s16$run == "main" & s16$library_type != "sample")   # mock community + other libraries
s16_oak <- sum(s16$campaign %in% "2022_felled_oak")
s16_oak_all <- sum(s16$run == "felled_oak")
gas <- read_table("internal_gas")
n_gas <- sum(gas$campaign == "2021_survey")
iso <- read_table("isotopes")
iso_samp <- sum(iso$in_whole_tree_set)       # whole-tree d13C set (code/lib/isotope_samples.R rule)

# ================================================= D. FELLED BLACK OAK (2022) ==
bo_flux <- nrow(oak)
bo_gas  <- sum(gas$campaign == "2022_felled_oak")
dbo <- dd[dd$campaign %in% "2022_felled_oak" & dd$qc_pass, ]
bo_ddpcr <- uq(dbo$sample_id)                    # felled-tree ddPCR samples (mcrA only)
bo_ext <- sum(smp$campaign == "2022_felled_oak") # extractions
bo_tiss <- table(smp$position[smp$campaign == "2022_felled_oak"])

# ============================================== E. WOOD PROPERTIES + INVENTORY ==
wc <- read_table("wood_cores")
wood_dens  <- uq(wc$tree_id[!is.na(wc$density_g_cm3)])
wood_moist <- uq(wc$tree_id[!is.na(wc$moisture_dry_pct)])
# inventory_stems.csv is the whole stand: the 2019 plot census merged with the
# supplementary by-tag survey. Individuals are physical trees (tree_id; tags are not
# unique). Census status is a census-only attribute that the tables do not carry (the
# stand inventory is live stems by design); the count is taken from the 2019 census
# sheet itself, as recorded.
inv <- read_table("inventory_stems")
inv_stems <- nrow(inv); inv_dbh <- sum(!is.na(inv$dbh_cm)); inv_indiv <- uq(inv$tree_id)
inv_spp <- uq(inv$species_code)
inv_area <- round(diff(range(inv$PX,na.rm=TRUE))*diff(range(inv$PY,na.rm=TRUE))/1e4, 2)
inv_census   <- sum(inv$source == "fg19", na.rm = TRUE)
inv_bytag    <- sum(inv$source == "bytag", na.rm = TRUE)
census <- read.csv("data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv", check.names = FALSE, fileEncoding = "UTF-8-BOM")
inv_live     <- sum(census$Status == "LI", na.rm = TRUE); census_rows <- nrow(census)
inv_in_stand <- sum(inv$in_stand, na.rm = TRUE)
inv_spp_stand <- uq(inv$species_code[inv$in_stand %in% TRUE])   # Table 1 reports the stand (the budget basis)

# ================================================================== REPORT =====
sink("outputs/audit/campaign_counts.txt")
cat("CAMPAIGN / SAMPLE-SIZE ACCOUNTING\n")
cat("Rule: keep all fluxes (every deployment in flux_stem/flux_soil, no QC). Black oak = separate section.\n")
cat("Source: code/09_tables_stats/02_campaign_counts.R\n"); cat(strrep("=",78),"\n\n")

cat("A. STEM FLUX (deployments with a CH4 value)\n")
P("  Monthly 2020-21 (breast ht; positions U/I/WD/WS): %d fluxes | %d trees | %d species (%d tags unresolved)", n_monthly, mon_trees, mon_spp, mon_unres)
P("      positions (fluxes): %s", paste(sprintf("%s=%d",names(mon_pos),mon_pos), collapse="  "))
P("  Height 2021 (3 heights; upland-focused, in plot): %d fluxes | %d trees | %d spp", n_height, h_trees, h_spp)
P("      by height (cm): %s", paste(sprintf("%s=%d",names(h_by),h_by), collapse="  "))
P("  Cross-species 2023 (breast ht; spatial in plot):  %d fluxes | %d trees | %d spp", n_2023, t23_trees, t23_spp)
P("  ---------------------------------------------")
P("  Survey stem: %d   + felled oak %d   = GRAND STEM TOTAL %d", stem_total, bo_flux, stem_total + bo_flux)
P("  Soil (monthly only): %d", n_soil)
P("  Unique flux trees: survey %d  + felled oak 1  = %d  (canonical tree_id)", union_trees, union_trees + 1)

cat("\nB. COVARIATES per flux measurement (VWC / DBH / temp)\n")
P("  Monthly 2020-21 : %s", cov_line(mon))
P("  Height 2021     : %s", cov_line(hd))
P("  Cross-species 23: %s", cov_line(y23))

cat("\nC. MOLECULAR / COMMUNITY / GAS (summer-2021 survey; incl. survey black oaks bo1/bo2/bowser)\n")
P("  ddPCR survey samples: %d  (wood %d / soil %d; + %d run but failing QC)  targets: %s", dd_surv, dd_w, dd_s, dd_failed, paste(dd_genes, collapse=", "))
P("  16S survey samples  : %d  (wood inner/outer %d / soil org/min %d) -- paired to ddPCR", s16_w+s16_s, s16_w, s16_s)
P("      (+%d mock-community and other libraries of the main run — excluded from the paired survey count)", s16_none)
P("  Internal gas (survey): %d      Isotopes (Picarro): %d runs, %d whole-tree samples", n_gas, nrow(iso), iso_samp)

cat("\nD. FELLED BLACK OAK (single tree, Oct 2022) — separate individual from the 2021 survey black oaks\n")
P("  Stem flux    : %d (heights 0.5,1.25,2,4,6,8,10 m + basal seam)", bo_flux)
P("  Internal gas : %d samples (one per height). No isotopes for this tree.", bo_gas)
P("  ddPCR        : %d samples (mcrA only)", bo_ddpcr)
P("  16S tissues  : %d biological (heartwood/sapwood + bark/branch/roots/foliage/litter/rot/soil); %d libraries in the felled-oak run", s16_oak, s16_oak_all)
P("  Extractions  : %d", bo_ext)

cat("\nMOLECULAR ASSAY COUNT (ddPCR reactions = sample x target)\n")
P("  survey by target: %s", paste(sprintf("%s=%d",names(dd_by_target),dd_by_target), collapse="  "))
P("  survey ddPCR reactions: %d (%d data points incl. loose+strict) + felled oak %d = %d reactions",
  dd_reactions, dd_datapoints, nrow(dbo), dd_reactions + nrow(dbo))
P("  + 16S libraries %d  => ~%d molecular assays total", (s16_w+s16_s)+s16_oak, dd_reactions + nrow(dbo) + (s16_w+s16_s)+s16_oak)

cat("\nGRAND TOTALS (survey + felled oak):\n")
P("  Stem fluxes: %d   Soil fluxes: %d   Trees: %d", stem_total + bo_flux, n_soil, union_trees + 1)
P("  ddPCR samples: %d   16S: %d   Internal gas: %d   Isotopes: %d", dd_surv + bo_ddpcr, (s16_w+s16_s) + s16_oak, n_gas + bo_gas, iso_samp)

cat("\nE. WOOD PROPERTIES (2021 cores)\n")
P("  Wood density: %d trees | wood moisture: %d trees (inner/outer/middle)", wood_dens, wood_moist)
P("  Internal gas measured for CH4, CO2, N2O, O2 (the 157 survey samples)")

cat("\nF. INVENTORY / UPSCALING (all censused stems)\n")
P("  Total stems: %d (%d with DBH) | %d individuals (tree_id) | %d species", inv_stems, inv_dbh, inv_indiv, inv_spp)
P("      sources: %d from the 2019 census, %d from the by-tag survey", inv_census, inv_bytag)
P("      in the censused stand (the budget basis): %d", inv_in_stand)
P("      live (LI) in the 2019 census sheet: %d of %d rows (%d census stems in the inventory); by-tag stems carry no census status",
  inv_live, census_rows, inv_census)
P("      extent-based area ~%s ha (CONFIRM censused plot area)", inv_area)

cat("\nNOTE: covariates (VWC, DBH, air/soil/stem temp) are recorded per flux measurement / per tree\n")
cat("      as attributes — a Methods sentence, not a separate sample-size stream.\n")
sink()
cat(readLines("outputs/audit/campaign_counts.txt"), sep="\n")

# machine-readable
write.csv(data.frame(
  campaign=c("Monthly survey","Height+molecular","Felled black oak","Cross-species","Inventory"),
  period=c("Jun2020-May2021","Summer2021","Oct2022","Summer2023","-"),
  location=c("hillslope U/I/wetland-margin","upland, in plot","stand (1 tree)","spatial, in plot","census plot"),
  heights=c("breast","50/125/200cm","0.5-10m+seam","breast","-"),
  trees=c(mon_trees,h_trees,1,t23_trees,inv_in_stand), species=c(mon_spp,h_spp,1,t23_spp,inv_spp_stand),
  stem_flux=c(n_monthly,n_height,bo_flux,n_2023,NA), soil_flux=c(n_soil,NA,NA,NA,NA),
  ddpcr=c(NA,dd_surv,bo_ddpcr,NA,NA), s16=c(NA,s16_w+s16_s,s16_oak,NA,NA),
  gas=c(NA,n_gas,bo_gas,NA,NA), isotopes=c(NA,iso_samp,NA,NA,NA),
  covariates=c("air T/DBH; soil by attach_soil","air T/DBH; soil holes; wood","tree","air T/DBH; soil by attach_soil",NA)
), "outputs/data/campaign_counts.csv", row.names=FALSE)
cat("\n\nWrote outputs/audit/campaign_counts.txt, outputs/data/campaign_counts.csv\n")
