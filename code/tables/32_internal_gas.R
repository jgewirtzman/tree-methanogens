#!/usr/bin/env Rscript
# ==============================================================================
# 32_internal_gas.R -- internal_gas.csv (data/tables/)
# ------------------------------------------------------------------------------
# One row per stem gas sample analysed by GC, three campaigns:
#   2021_survey             whole-stem samples (bark to pith), GC run 2022-09/10;
#                           calibrated by code/01_import/03_process_internal_gas.R
#                           (data/processed/internal_gas/sample_data_only.csv)
#   2024_heartwood_sapwood  paired heartwood / sapwood samples (+ 3 atmosphere
#                           vials, no tree), GC run 2024-10-10/11; calibrated by
#                           code/01_import/04_process_internal_gas_2024.R
#                           (data/processed/internal_gas/stem_gas_2024_calibrated.csv);
#                           10 vials of its sample list have no GC result (rows
#                           kept, values and analysis_date NA)
#   2022_felled_oak         trunk gas at 7 heights along the felled black oak,
#                           GC run 2022-10-07/08 (data/raw/field_data/black_oak/
#                           quve_gc_data.csv); calibrated HERE, moved from
#                           fig07_decay-methanogenesis.R and 20_manuscript_statistics.R
#                           (identical code in both): CH4 linear through N2, SB1-SB3
#                           for peak area <= 1000, else linear through all GHG
#                           standards (N2, SB1-SB6, outdoor-air series); O2 linear
#                           through the low O2 standards (outdoor air 2-4, O2
#                           standards 1-4) up to their largest area, else through
#                           those plus outdoor air 5-7 and O2 standard 5.
#                           CO2 and N2O are not calibrated for this run (NA). The
#                           48 bark / sapwood / heartwood vials in the same file
#                           are tissue incubations, not internal gas, and are not
#                           carried.
# Values are as calibrated: nothing clamped or substituted. LOD flags exist for
# the 2021 run only (3 x SD of SB1 replicate predictions). CH4_above_top_standard:
# peak area above the highest certified standard of the run (2021 SB6, value held
# at SB6 = lower bound; 2024 SB5; felled oak SB6).
# sample_id: the GC lab ID (2021), sample name (2024), oak_<lab no> (felled oak).
# Two 2021 samples carry lab ID "BLANK" (trees BOWSER, H301): BLANK_<lab no>.
# tree_id: tree_tags.csv; the 2024 names "RO??" and "BB11/4?" map to
# RO_unidentified_2024 and BB11_or_BB4_2024. species_recorded: the species code on
# the 2021 sample sheet (Species.ID), as written; NA for the other campaigns.
#
# Checks: the felled-oak values equal those of the fig07 / 20_manuscript_statistics.R code, run verbatim
# (to 1e-12, relative: the two fits differ only in floating-point summation order).
#
# Writes data/tables/internal_gas.csv, outputs/audit/internal_gas_build.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages(library(tidyverse))
OUT <- "data/tables"
G21 <- "data/processed/internal_gas/sample_data_only.csv"
G24 <- "data/processed/internal_gas/stem_gas_2024_calibrated.csv"
OAKDIR <- "data/raw/field_data/black_oak"
TT_ALL <- read.csv(file.path(OUT, "tree_tags.csv"), colClasses = "character", na.strings = "")
tags_for <- function(file, column) {
  t <- TT_ALL %>% filter(source_file == file, source_column == column) %>% distinct(tag = trimws(tag_recorded), tree_id, .keep_all = TRUE)
  stopifnot(all(is.na(t$species_recorded)), !anyDuplicated(t$tag)); t
}
mdy <- function(x) as.Date(x, "%m/%d/%y")

# ---- 2021 survey -------------------------------------------------------------------
a <- read.csv(G21, stringsAsFactors = FALSE)
t21 <- tags_for(G21, "Tree.ID")
g21 <- a %>% transmute(
  sample_id = ifelse(Lab.ID == "BLANK", paste0("BLANK_", Lab.No), Lab.ID),
  tree_id = t21$tree_id[match(trimws(Tree.ID), t21$tag)], species_recorded = Species.ID, campaign = "2021_survey",
  tissue = "whole stem (bark to pith)", height_m = NA_real_, analysis_date = as.Date(FID.Date),
  CH4_ppm = CH4_concentration, CO2_ppm = CO2_concentration, N2O_ppm = N2O_concentration, O2_ppm = O2_concentration,
  CH4_below_lod, CO2_below_lod, N2O_below_lod, CH4_above_top_standard = CH4_above_SB6)
stopifnot(!anyNA(g21$tree_id))

# ---- 2024 heartwood / sapwood ------------------------------------------------------------
b <- read.csv(G24, stringsAsFactors = FALSE)
top24 <- max(b$CH4_Area[b$Sample.ID %in% c("SB1", "SB3", "SB4", "SB5")], na.rm = TRUE)   # the run's certified tanks
t24 <- tags_for(G24, "Tree.No")
g24 <- b %>% filter(Sample.Type == "Sample") %>% transmute(
  sample_id = Sample.ID, tree_id = t24$tree_id[match(trimws(Tree.No), t24$tag)], campaign = "2024_heartwood_sapwood",
  tissue = recode(Tissue, H = "heartwood", S = "sapwood", .missing = "atmosphere"), height_m = NA_real_,
  analysis_date = mdy(FID_Date), CH4_ppm = CH4_calibrated_ppm, CO2_ppm = CO2_calibrated_ppm,
  N2O_ppm = NA_real_, O2_ppm = NA_real_, CH4_below_lod = NA, CO2_below_lod = NA, N2O_below_lod = NA,
  CH4_above_top_standard = CH4_Area > top24)
stopifnot(identical(is.na(g24$tree_id), g24$tissue == "atmosphere"))

# ---- 2022 felled oak: calibration ------------------------------------------------------------
rd <- function(f) { d <- read.csv(file.path(OAKDIR, f), fileEncoding = "UTF-8-BOM"); colnames(d) <- make.names(colnames(d)); d }
GC <- rd("quve_gc_data.csv"); O2S <- rd("o2_standards.csv"); GHGS <- rd("ghg_standards.csv")
ppm <- function(x) as.numeric(gsub(",", "", x))
GC$O2.Area[GC$Sample.Type == "N2"] <- NA
o2d <- O2S %>% left_join(GC, by = c("Sample" = "Sample.Type")) %>% mutate(O2.Area = as.numeric(O2.Area), O2_ppm = ppm(X.O2...ppm.))
O2_LOW <- c("Outdoor Air 2", "Outdoor Air 3", "Outdoor Air 4", paste("Oxygen Standard", 1:4))
O2_HIGH <- c("Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7", "Oxygen Standard 5")
o2_lo <- lm(O2_ppm ~ O2.Area, o2d %>% filter(Sample %in% O2_LOW))
o2_hi <- lm(O2_ppm ~ O2.Area, o2d %>% filter(Sample %in% c(O2_LOW, O2_HIGH)))
o2_break <- max(o2d$O2.Area[o2d$Sample %in% O2_LOW], na.rm = TRUE)
ghd <- GHGS %>% left_join(GC, by = c("Sample" = "Sample.Type")) %>%
  mutate(CH4.Area = as.numeric(CH4.Area), CH4_ppm = ppm(X.CH4...ppm.)) %>% filter(!is.na(CH4.Area), !is.na(CH4_ppm))
ch4_lo <- lm(CH4_ppm ~ CH4.Area, ghd %>% filter(Sample %in% c("N2", "SB1", "SB2", "SB3")))
ch4_hi <- lm(CH4_ppm ~ CH4.Area, ghd)
CH4_BREAK_AREA <- 1000
top_oak <- max(ghd$CH4.Area[ghd$Sample == "SB6"])
oak <- GC %>% filter(!is.na(Lab.ID), Tree.Tissue == "Trunk Gas") %>% transmute(
  sample_id = paste0("oak_", Lab.ID), tree_id = "QUVE_felled_2022", campaign = "2022_felled_oak",
  tissue = "whole stem (trunk gas)", height_m = as.numeric(Tree.Height), analysis_date = mdy(FID.Date),
  CH4_ppm = ifelse(CH4.Area <= CH4_BREAK_AREA, predict(ch4_lo, data.frame(CH4.Area = CH4.Area)),
                   predict(ch4_hi, data.frame(CH4.Area = CH4.Area))),
  CO2_ppm = NA_real_, N2O_ppm = NA_real_,
  O2_ppm = ifelse(O2.Area <= o2_break, predict(o2_lo, data.frame(O2.Area = O2.Area)),
                  predict(o2_hi, data.frame(O2.Area = O2.Area))),
  CH4_below_lod = NA, CO2_below_lod = NA, N2O_below_lod = NA, CH4_above_top_standard = CH4.Area > top_oak)

IG <- bind_rows(g21, g24, oak) %>% mutate(across(c(CH4_ppm, CO2_ppm, N2O_ppm, O2_ppm), as.numeric))
stopifnot(!anyDuplicated(IG$sample_id))

# ==============================================================================
# CHECKS
# ==============================================================================
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
near <- function(a, b) same(a, b) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(a)))   # floating point
# (b) felled oak vs the fig07 code, run verbatim
fig07 <- local({
  O2s<-read.csv("data/raw/field_data/black_oak/o2_standards.csv",fileEncoding="UTF-8-BOM")
  GHGs<-read.csv("data/raw/field_data/black_oak/ghg_standards.csv",fileEncoding="UTF-8-BOM")
  GC<-read.csv("data/raw/field_data/black_oak/quve_gc_data.csv",fileEncoding="UTF-8-BOM")
  colnames(GC)<-make.names(colnames(GC));colnames(O2s)<-make.names(colnames(O2s));colnames(GHGs)<-make.names(colnames(GHGs))
  GC$O2.Area[GC$Sample.Type=="N2"]<-NA
  O2d<-O2s%>%left_join(GC,by=c("Sample"="Sample.Type"));O2d$O2.Area<-as.numeric(O2d$O2.Area);O2d$O2_ppm<-as.numeric(gsub(",","",O2d$X.O2...ppm.))
  lo<-c("Outdoor Air 2","Outdoor Air 3","Outdoor Air 4","Oxygen Standard 1","Oxygen Standard 2","Oxygen Standard 3","Oxygen Standard 4")
  hi<-c("Outdoor Air 5","Outdoor Air 6","Outdoor Air 7","Oxygen Standard 5")
  lod<-O2d%>%filter(Sample%in%lo);hid<-bind_rows(lod,O2d%>%filter(Sample%in%hi))
  O2lo<-lm(O2_ppm~O2.Area,lod);O2hi<-lm(O2_ppm~O2.Area,hid);mxO2<-max(lod$O2.Area,na.rm=TRUE)
  GHGd<-GHGs%>%left_join(GC,by=c("Sample"="Sample.Type"))
  for(g in c("N2O","CH4","CO2")){GHGd[[paste0(g,".Area")]]<-as.numeric(GHGd[[paste0(g,".Area")]]);GHGd[[paste0(g,"_ppm")]]<-as.numeric(gsub(",","",GHGd[[paste0("X.",g,"...ppm.")]]))}
  GHGd<-GHGd%>%filter(!is.na(CH4.Area)&!is.na(CH4_ppm))
  lord<-GHGd%>%filter(Sample%in%c("N2","SB1","SB2","SB3"))
  CH4lo<-lm(CH4_ppm~CH4.Area,lord);CH4hi<-lm(CH4_ppm~CH4.Area,GHGd)
  GC<-GC%>%mutate(O2_concentration=if_else(O2.Area<=mxO2,predict(O2lo,data.frame(O2.Area=O2.Area)),predict(O2hi,data.frame(O2.Area=O2.Area))),
    CH4_concentration=if_else(CH4.Area<=1000,predict(CH4lo,data.frame(CH4.Area=CH4.Area)),predict(CH4hi,data.frame(CH4.Area=CH4.Area))))
  int_gas<-GC%>%filter(!is.na(Lab.ID),Tree.Tissue=="Trunk Gas");int_gas$Tree.Height<-as.numeric(int_gas$Tree.Height)
  int_gas
})
stopifnot(nrow(fig07) == nrow(oak))
obad <- c(CH4 = sum(!near(fig07$CH4_concentration, oak$CH4_ppm)), O2 = sum(!near(fig07$O2_concentration, oak$O2_ppm)),
          height = sum(fig07$Tree.Height != oak$height_m))

sink(out_path("internal_gas_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("internal_gas: %d samples (%s)\n", nrow(IG), paste(names(table(IG$campaign)), table(IG$campaign), collapse = ", ")))
cat(sprintf("trees: %d; samples without a tree (2024 atmosphere vials): %d\n", n_distinct(IG$tree_id, na.rm = TRUE), sum(is.na(IG$tree_id))))
cat(sprintf("CH4 above the top standard: %s\n", paste(IG$sample_id[IG$CH4_above_top_standard %in% TRUE], collapse = ", ")))
cat("\n(b) felled oak vs fig07 / 20_manuscript_statistics.R calibration -- values differing:\n"); print(obad)
cat("\nFelled-oak trunk gas:\n"); print(as.data.frame(oak %>% select(sample_id, height_m, analysis_date, CH4_ppm, O2_ppm)), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("internal_gas_build.txt", "audit")), sep = "\n")
if (any(obad > 0)) stop("internal_gas does not reproduce today's values (outputs/audit/internal_gas_build.txt)")

write_csv(IG, file.path(OUT, "internal_gas.csv"), na = "")
