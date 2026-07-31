#!/usr/bin/env Rscript
# ==============================================================================
# rev_compile_zenodo_semirigid.R   (REVISION override of compile_zenodo_datasets §1)
# Rebuilds the public Zenodo semirigid_chamber_flux.csv to INCLUDE the 7 recovered
# untagged/dead-snag monthly trees, so the archived flux dataset matches the
# revised manuscript (500 trees / 1,200 stem fluxes). Tagged tree + soil blocks
# are reproduced verbatim from compile_zenodo_datasets.R; the untagged block is
# added from the goFlux output (untagged_monthly_fluxes.csv), carrying real
# geometry/model/quality and a live/dead note. Original compiler left untouched.
# Run this AFTER compile_zenodo_datasets.R (it overwrites data/compiled/semirigid_chamber_flux.csv).
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse)); options(warn=-1)
out_dir<-"data/compiled"; dir.create(out_dir,showWarnings=FALSE)

# --- tagged tree stems (verbatim schema from compile_zenodo_datasets.R) ---
tree_flux<-read_csv("data/processed/flux/semirigid_tree_final_complete_dataset.csv",show_col_types=FALSE)
tree_clean<-tree_flux %>% transmute(
  unique_id=UniqueID, date=as.Date(Date), tree_tag=as.character(`Plot Tag`),
  landscape_position=`Plot Letter`, chamber_id=Chamber_ID, chamber_area_m2=Area,
  chamber_vol_L=Vtot, chamber_temp_C=Tcham, CH4_best_flux_nmol_m2_s=CH4_best.flux.y,
  CH4_model=CH4_model.y, CH4_quality=CH4_quality.check.y, CH4_LM_r2=CH4_LM.r2.y,
  CH4_LM_pval=CH4_LM.p.val.y, CO2_best_flux_nmol_m2_s=CO2_best.flux, CO2_model=CO2_model,
  CO2_quality=CO2_quality.check, measurement_type="tree_stem", notes=Notes)

# --- recovered untagged/dead-snag tree stems (goFlux output) ---
u<-read_csv("outputs/revision/untagged_monthly_fluxes.csv",show_col_types=FALSE)
untag_clean<-u %>% transmute(
  unique_id=UniqueID,
  date=as.Date(sub(".*_(\\d{8})(_.*)?$","\\1",UniqueID),format="%Y%m%d"),
  tree_tag=tree_id,                                        # per-tree label (measurements share it)
  landscape_position=dplyr::recode(Plot_Type, U="U", I="I", W="WS"),
  chamber_id=NA_character_, chamber_area_m2=Area, chamber_vol_L=Vtot, chamber_temp_C=NA_real_,
  CH4_best_flux_nmol_m2_s=CH4_best.flux, CH4_model=CH4_model, CH4_quality=CH4_quality,
  CH4_LM_r2=CH4_LM.r2, CH4_LM_pval=NA_real_, CO2_best_flux_nmol_m2_s=CO2_best.flux,
  CO2_model=NA_character_, CO2_quality=NA_character_, measurement_type="tree_stem",
  notes=paste0("recovered untagged monthly tree; ", ifelse(dead,"dead snag","live")))

# --- soil (verbatim) ---
soil_flux<-read_csv("data/processed/flux/semirigid_tree_final_complete_dataset_soil_CORRECTED.csv",show_col_types=FALSE)
soil_clean<-soil_flux %>% transmute(
  unique_id=UniqueID, date=as.Date(Date), tree_tag=as.character(`Plot Tag`),
  landscape_position=`Plot letter`, chamber_id=NA_character_, chamber_area_m2=Area,
  chamber_vol_L=Vtot, chamber_temp_C=Tcham, CH4_best_flux_nmol_m2_s=CH4_best.flux,
  CH4_model=CH4_model, CH4_quality=CH4_quality.check, CH4_LM_r2=CH4_LM.r2,
  CH4_LM_pval=CH4_LM.p.val, CO2_best_flux_nmol_m2_s=CO2_best.flux, CO2_model=CO2_model,
  CO2_quality=CO2_quality.check, measurement_type="soil", notes=Notes)

semirigid<-bind_rows(tree_clean, untag_clean, soil_clean)
write_csv(semirigid, file.path(out_dir,"semirigid_chamber_flux.csv"))
cat(sprintf("Zenodo semirigid_chamber_flux.csv rebuilt: %d tree + %d untagged + %d soil = %d rows\n",
    nrow(tree_clean), nrow(untag_clean), nrow(soil_clean), nrow(semirigid)))
cat(sprintf("  tree_stem rows: %d  (untagged tagged as 'recovered untagged monthly tree' in notes)\n",
    sum(semirigid$measurement_type=="tree_stem")))
