#!/usr/bin/env Rscript
# ==============================================================================
# rev_stat_campaign_counts.R  —  REPRODUCIBLE campaign / sample-size accounting
# ------------------------------------------------------------------------------
# Regenerates every N in the study-design table from the archived data, under
# one consistent, documented definition per quantity.
#
# FLUX RULE (definitive): count every chamber deployment that produced a CH4
# flux value (non-NA best.flux). **No quality filtering** — QC exclusions are a
# separate, downstream analysis concern. Applied identically to all campaigns.
#
# NOTE: the 2021 survey includes black-oak individuals (tree_ids bo1/bo2/bowser,
# species QUVE) — these are SURVEY trees and stay in the survey totals. They are
# distinct from the single 2022 FELLED black oak, which is a separate section.
#
# Run from repo root:  Rscript code/revision/rev_stat_campaign_counts.R
# Writes: outputs/revision/campaign_counts.txt  (+ .csv table)
# ==============================================================================
suppressMessages(library(dplyr)); options(warn = -1, stringsAsFactors = FALSE)
nn <- function(x) sum(!is.na(x)); uq <- function(x) length(unique(x[!is.na(x) & x != ""]))
P  <- function(...) cat(sprintf(...), "\n")
fd <- "data/processed/flux"; comp <- "data/compiled"

# ============================================================== A. STEM FLUX ===
srt <- read.csv(file.path(fd, "semirigid_tree_final_complete_dataset.csv"))       # monthly 2020-21
srs <- read.csv(file.path(fd, "semirigid_tree_final_complete_dataset_soil.csv"))  # monthly soil
ch4 <- read.csv(file.path(fd, "CH4_best_flux_lgr_results.csv"))                   # 2021 height flux
aux <- read.csv(file.path(fd, "goflux_auxfile.csv"))                             # 2021 height metadata
y23 <- read.csv(file.path(fd, "methanogen_tree_flux_complete_dataset.csv"))       # 2023 cross-species

n_monthly <- nn(srt$CH4_best.flux.x); n_soil <- nn(srs$CH4_best.flux)
n_height  <- nn(ch4$best.flux);       n_2023 <- nn(y23$CH4_best.flux)
stem_total <- n_monthly + n_height + n_2023

# 2021 height flux by measurement height
hd <- merge(ch4, aux[, c("UniqueID","measurement_height","tree_id","species","plot")], by="UniqueID", all.x=TRUE)
hd <- hd[!is.na(hd$best.flux), ]
h_by <- sort(table(hd$measurement_height), decreasing = TRUE)
mon_trees <- uq(srt$Plot.Tag[!is.na(srt$CH4_best.flux.x)])
mon_pos   <- table(srt$Plot.Letter[!is.na(srt$CH4_best.flux.x)])   # U/I/WD/WS
# monthly species: not in the flux file; resolve from the field tree list (fallback: 2023 map)
mon_tags <- unique(as.character(srt$Plot.Tag[!is.na(srt$CH4_best.flux.x)]))
ymt <- read.csv("data/raw/field_data/static_chamber_field/YM_trees_measured.csv.csv", check.names=FALSE)
sp_lut <- setNames(as.character(ymt$Species), as.character(ymt$Label))
sp23   <- setNames(as.character(y23$Species.Code), as.character(y23$Tree.Tag))
mon_species <- ifelse(mon_tags %in% names(sp_lut), sp_lut[mon_tags],
                      ifelse(mon_tags %in% names(sp23), sp23[mon_tags], NA))
mon_spp <- uq(mon_species); mon_unres <- sum(is.na(mon_species))
h_trees <- length(unique(paste(hd$plot, hd$tree_id))); h_spp <- uq(hd$species)
t23_ok <- !is.na(y23$CH4_best.flux); t23_trees <- uq(y23$Tree.Tag[t23_ok]); t23_spp <- uq(y23$Species.Code[t23_ok])

# unique trees across campaigns, reconciled via the ID crosswalk.
# NOTE (verified): tree_id_comprehensive_mapping.csv is 2021-only — it collapses the
# 2021 survey trees' name variants (SM6, BO1, AB 3/5 ...) but does NOT contain the
# numeric monthly/2023 tags, so it cannot bridge campaigns. It is applied to the 2021
# IDs; monthly & 2023 numeric tags reconcile to each other directly.
norm <- function(x) toupper(trimws(as.character(x)))
canon <- function(x) norm(x)
map_f <- "data/processed/tree_data/tree_id_comprehensive_mapping.csv"
if (file.exists(map_f)) {
  M <- read.csv(map_f, check.names = FALSE)
  vcols <- grep("^name_in_|^variant_|^name_variants$|primary_id|Tree_ID_normalized", names(M), value = TRUE)
  prim <- if ("primary_id" %in% names(M)) M$primary_id else M$Tree_ID_normalized
  lut <- c(); for (cc in vcols) { v <- norm(M[[cc]]); ok <- v != "" & !is.na(v); lut[v[ok]] <- norm(prim[ok]) }
  canon <- function(x) { x <- norm(x); ifelse(x %in% names(lut), lut[x], x) }
}
ids <- unique(c(canon(srt$Plot.Tag[!is.na(srt$CH4_best.flux.x)]), canon(hd$tree_id), canon(y23$Tree.Tag[t23_ok])))
union_trees <- length(ids)

# ============================================================= B. COVARIATES ===
cov_has <- function(f, pat){ h <- names(read.csv(f, nrows=1, check.names=FALSE)); any(grepl(pat, h, ignore.case=TRUE)) }
cov_monthly <- cov_has(file.path(fd,"semirigid_tree_final_complete_dataset.csv"), "vwc|dbh|diam|temp|moist")
cov_2023    <- cov_has(file.path(fd,"methanogen_tree_flux_complete_dataset.csv"), "vwc")
cov_2021    <- cov_has("data/processed/integrated/merged_tree_dataset_final.csv", "VWC")

# ============================================== C. MOLECULAR / COMMUNITY / GAS ==
dd  <- read.csv(file.path(comp,"ddpcr_gene_abundances.csv")); dds <- dd[!duplicated(dd$sample_id), ]
# survey ddPCR = all wood/soil cores (includes the 2021 survey black oaks bo1/bo2/bowser)
dd_surv <- nrow(dds); dd_w <- sum(dds$material=="Wood"); dd_s <- sum(dds$material=="Soil")
dd_genes <- sort(unique(dd$target_gene))
# 16S: pair to the ddPCR sampling by counting core-typed samples (inner/outer wood,
# organic/mineral soil). The file also holds ~"None"/control samples that inflate a raw total.
s16 <- read.csv(file.path(comp,"sample_metadata_16S.csv"), check.names=FALSE)
wood_ct <- c("Inner","Outer"); soil_ct <- c("Organic","Mineral")
s16_w <- sum(s16$Material=="Wood" & s16$core_type %in% wood_ct)
s16_s <- sum(s16$Material=="Soil" & s16$core_type %in% soil_ct)
s16_none <- sum(s16$Material %in% c("Wood","Soil") & !(s16$core_type %in% c(wood_ct, soil_ct)))
# 2022 felled-tree tissues: real biological samples exclude core_type None/empty
s16_oak_all <- sum(s16$Material=="QUVE")
s16_oak <- sum(s16$Material=="QUVE" & !(s16$core_type %in% c("None","empty","")) & !is.na(s16$core_type))
n_gas <- nn(read.csv("data/processed/internal_gas/sample_data_only.csv", check.names=FALSE)[[1]])
iso <- read.csv("data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv", check.names=FALSE)
iso_stcol <- grep("Sample.Type|Sample Type", names(iso), value=TRUE)[1]
iso_samp <- if(!is.na(iso_stcol)) sum(grepl("Sample", iso[[iso_stcol]])) else nrow(iso)

# ================================================= D. FELLED BLACK OAK (2022) ==
bof <- read.csv("data/compiled/black_oak_experiment.csv"); bo_flux <- nrow(bof)
bog <- read.csv("data/raw/field_data/black_oak/quve_gc_data.csv", check.names=FALSE)
bog_st <- grep("Sample Type|Sample.Type", names(bog), value=TRUE)[1]
bo_gas <- if(!is.na(bog_st)) sum(grepl("Internal Gas", bog[[bog_st]])) else NA
bom <- read.csv("data/raw/ddpcr/black_oak_mcrA.csv", check.names=FALSE)
bo_sid <- grep("Sample ID|Sample Name", names(bom), value=TRUE)[1]
bo_ddpcr <- length(unique(bom[[bo_sid]]))                       # total felled-tree ddPCR samples (mcrA)
bo_mcra_h <- length(unique(bom[[grep("Height", names(bom))[1]]]))
bo_ext <- nrow(read.csv("data/processed/molecular/black_oak/bo_extraction_mass.csv", check.names=FALSE))
bo_tiss <- table(s16$core_type[s16$Material=="QUVE"])

# ============================================== E. WOOD PROPERTIES + INVENTORY ==
tp <- read.csv(file.path(comp,"tree_properties.csv"), check.names=FALSE)
wood_dens <- sum(apply(!is.na(tp[,grep("density_final",names(tp))]),1,any))
wood_moist<- sum(apply(!is.na(tp[,grep("moisture_dry_percent",names(tp))]),1,any))
inv <- read.csv(file.path(comp,"forest_inventory.csv"))
inv_stems <- nrow(inv); inv_dbh <- sum(!is.na(inv$dbh_mm)); inv_indiv <- uq(inv$tag)
inv_live <- sum(inv$status=="LI", na.rm=TRUE); inv_spp <- uq(inv$species_code)
inv_area <- round(diff(range(inv$x_m,na.rm=TRUE))*diff(range(inv$y_m,na.rm=TRUE))/1e4, 2)

# ================================================================== REPORT =====
sink("outputs/revision/campaign_counts.txt")
cat("CAMPAIGN / SAMPLE-SIZE ACCOUNTING\n")
cat("Rule: keep all fluxes (non-NA best.flux, no QC). Black oak = separate section.\n")
cat("Source: code/revision/rev_stat_campaign_counts.R\n"); cat(strrep("=",78),"\n\n")

cat("A. STEM FLUX (deployments with a CH4 value)\n")
P("  Monthly 2020-21 (breast ht; positions U/I/WD/WS): %d fluxes | %d trees | %d species (%d tags unresolved)", n_monthly, mon_trees, mon_spp, mon_unres)
P("      positions (fluxes): %s", paste(sprintf("%s=%d",names(mon_pos),mon_pos), collapse="  "))
P("  Height 2021 (3 heights; upland-focused, in plot): %d fluxes | %d trees | %d spp", n_height, h_trees, h_spp)
P("      by height (cm): %s", paste(sprintf("%s=%d",names(h_by),h_by), collapse="  "))
P("  Cross-species 2023 (breast ht; spatial in plot):  %d fluxes | %d trees | %d spp", n_2023, t23_trees, t23_spp)
P("  ---------------------------------------------")
P("  Survey stem: %d   + felled oak %d   = GRAND STEM TOTAL %d", stem_total, bo_flux, stem_total + bo_flux)
P("  Soil (monthly only): %d", n_soil)
P("  Unique flux trees: survey %d  + felled oak 1  = %d  (crosswalk-reconciled; 2021-only crosswalk)", union_trees, union_trees + 1)

cat("\nB. COVARIATES per flux measurement (VWC / DBH / temp)\n")
P("  Monthly 2020-21 : merged in the Fig 1 workflow (06_soil_tree_timeseries.R) from soil-moisture/DBH/weather; not in the raw flux file")
P("  Height 2021     : in merged 2021 table (VWC, DBH, air/soil temp, wood moisture+density): %s", cov_2021)
P("  Cross-species 23: inline in flux file (DBH, VWC x3, air/soil/stem temp): %s", cov_2023)

cat("\nC. MOLECULAR / COMMUNITY / GAS (summer-2021 survey; incl. survey black oaks bo1/bo2/bowser)\n")
P("  ddPCR survey samples: %d  (wood %d / soil %d)  genes: %s", dd_surv, dd_w, dd_s, paste(dd_genes, collapse=", "))
P("  16S survey samples  : %d  (wood inner/outer %d / soil org/min %d) -- paired to ddPCR", s16_w+s16_s, s16_w, s16_s)
P("      (+%d 'None'/uncharacterised + control samples in the 16S file — excluded from the paired survey count)", s16_none)
P("  Internal gas (survey): %d      Isotopes (Picarro): %d runs, %d samples (rest standards)", n_gas, nrow(iso), iso_samp)

cat("\nD. FELLED BLACK OAK (single tree, Oct 2022) — separate individual from the 2021 survey black oaks\n")
P("  Stem flux    : %d (heights 0.5,1.25,2,4,6,8,10 m + basal seam)", bo_flux)
P("  Internal gas : %d samples (multiple heights; file also holds standards/blanks). No isotopes for this tree.", bo_gas)
P("  ddPCR        : %d samples (mcrA only)", bo_ddpcr)
P("  16S tissues  : %d biological (inner/outer + bark/branch/roots/foliage/litter/rot/soil); %d total incl. %d None/controls", s16_oak, s16_oak_all, s16_oak_all - s16_oak)
P("  Extractions  : %d", bo_ext)

cat("\nGRAND TOTALS (survey + felled oak):\n")
P("  Stem fluxes: %d   Soil fluxes: %d   Trees: %d", stem_total + bo_flux, n_soil, union_trees + 1)
P("  ddPCR: %d   16S: %d   Internal gas: %d   Isotopes: %d", dd_surv + bo_ddpcr, (s16_w+s16_s) + s16_oak, n_gas + bo_gas, iso_samp)

cat("\nE. WOOD PROPERTIES (2021 cores)\n")
P("  Wood density: %d trees | wood moisture: %d trees (inner/outer/middle)", wood_dens, wood_moist)
P("  Internal gas measured for CH4, CO2, N2O, O2 (the 157 survey samples)")

cat("\nF. INVENTORY / UPSCALING (all censused stems)\n")
P("  Total stems: %d (%d with DBH) | %d individuals | %d species", inv_stems, inv_dbh, inv_indiv, inv_spp)
P("      of which live (LI): %d; status not recorded: %d", inv_live, inv_stems - inv_live)
P("      extent-based area ~%s ha (CONFIRM censused plot area)", inv_area)

cat("\nNOTE: covariates (VWC, DBH, air/soil/stem temp) are recorded per flux measurement / per tree\n")
cat("      as attributes — a Methods sentence, not a separate sample-size stream.\n")
sink()
cat(readLines("outputs/revision/campaign_counts.txt"), sep="\n")

# machine-readable
write.csv(data.frame(
  campaign=c("Monthly survey","Height+molecular","Felled black oak","Cross-species","Inventory"),
  period=c("Jun2020-May2021","Summer2021","Oct2022","Summer2023","-"),
  location=c("hillslope U/I/wetland-margin","upland, in plot","stand (1 tree)","spatial, in plot","census plot"),
  heights=c("breast","50/125/200cm","0.5-10m+seam","breast","-"),
  trees=c(mon_trees,h_trees,1,t23_trees,inv_stems), species=c(mon_spp,h_spp,1,t23_spp,inv_spp),
  stem_flux=c(n_monthly,n_height,bo_flux,n_2023,NA), soil_flux=c(n_soil,NA,NA,NA,NA),
  ddpcr=c(NA,dd_surv,bo_ddpcr,NA,NA), s16=c(NA,s16_w+s16_s,s16_oak,NA,NA),
  gas=c(NA,n_gas,bo_gas,NA,NA), isotopes=c(NA,iso_samp,NA,NA,NA),
  covariates=c("merged in Fig1","VWC/DBH/temp+wood","tree","VWC/DBH/temp",NA)
), "outputs/revision/campaign_counts.csv", row.names=FALSE)
cat("\n\nWrote outputs/revision/campaign_counts.{txt,csv}\n")
