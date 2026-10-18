#!/usr/bin/env Rscript
# ==============================================================================
# 02_untagged_fluxes.R -- fluxes of the untagged monthly stems from their hand-picked windows
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root. Runs after 01_untagged_auxfile.R.
#
# The hand-picked closure windows (untagged_manID.rds, code/tools/click_untagged_windows.R)
# are matched to the auxfile by start time (a stable, species-independent key). Each
# window takes its measurement ID, species, live/dead status, chamber geometry (Area,
# Vtot) and chamber temperature from the auxfile, not the values saved at picking: the
# picking-time IDs carry the species label of the time (e.g. UNTAG_Upland_BELE_20200826_1_1
# for the dead red maple UNTAG_Upl_ACRU_dead). goFlux and best.flux are then re-run
# (package defaults).
# QC (per measurement): linear or HM r2 >= 0.7 and CH4 flux within -100 to 200 nmol m-2 s-1.
#
# Inputs:  data/raw/flux_windows/untagged/untagged_manID.rds
#          data/processed/flux/untagged_rescue/untagged_auxfile.csv
# Outputs (data/processed/flux/untagged_rescue/):
#          untagged_monthly_fluxes.csv  per measurement
#          untagged_monthly_trees.csv   per stem (median and mean of the QC'd CH4 fluxes)
# ==============================================================================
suppressMessages({library(goFlux);library(dplyr)})
manID<-readRDS("data/raw/flux_windows/untagged/untagged_manID.rds")
aux<-read.csv("data/processed/flux/untagged_rescue/untagged_auxfile.csv",check.names=FALSE)
aux$key<-aux$start.time_formatted
stopifnot(!any(duplicated(aux$key)))                      # start.time is a unique measurement key
manID$key<-format(as.POSIXct(manID$start.time_formatted,tz="UTC"),"%Y-%m-%d %H:%M:%S")
mp<-aux %>% select(key, UniqueID_c=UniqueID, Area_c=Area, Vtot_c=Vtot, Tcham_c=Tcham, site, species, dead, Sample, Dstem)
# IDs, geometry and chamber temperature from the auxfile (see header)
manID<-manID %>% select(-any_of(c("Area","Vtot","Tcham","site","species","Sample","Dstem"))) %>%
  left_join(mp,by="key") %>% rename(Area=Area_c, Vtot=Vtot_c, Tcham=Tcham_c)   # geometry and chamber temperature from the auxfile
stopifnot(!anyNA(manID$Tcham), !anyNA(manID$UniqueID_c),
          all(tapply(manID$UniqueID_c, manID$UniqueID, function(x) length(unique(x)) == 1)),
          all(tapply(manID$UniqueID, manID$UniqueID_c, function(x) length(unique(x)) == 1)))
cat("measurement IDs re-labelled to the corrected species:", sum(manID$UniqueID != manID$UniqueID_c & !duplicated(manID$UniqueID)), "\n")
manID$UniqueID<-manID$UniqueID_c; manID$UniqueID_c<-NULL
cat("re-keyed", length(unique(manID$key)), "measurements to notes-based species/geometry\n")

ch4<-goFlux(manID,"CH4dry_ppb"); ch4b<-best.flux(ch4); co2b<-best.flux(goFlux(manID,"CO2dry_ppm"))
meta<-manID %>% distinct(UniqueID,key,site,species,dead,Sample,Dstem,Area,Vtot)
f<-meta %>%
  left_join(ch4 %>% transmute(UniqueID,CH4_LM.r2=LM.r2,CH4_HM.r2=HM.r2),by="UniqueID") %>%
  left_join(ch4b %>% transmute(UniqueID,CH4_best.flux=best.flux,CH4_model=model,CH4_quality=quality.check,
                                 CH4_LM.p.val=LM.p.val,CH4_flux.term=flux.term,CH4_nb.obs=nb.obs,
                                 CH4_MDF=MDF,CH4_prec=prec),by="UniqueID") %>%   # precision terms for detection
  left_join(co2b %>% transmute(UniqueID,CO2_best.flux=best.flux),by="UniqueID") %>%
  mutate(Plot_Type=case_when(site=="Upland"~"U",site=="Intermediate"~"I",site=="Wetland"~"W"),
         tree_id=paste0("UNTAG_",substr(site,1,3),"_",species,ifelse(dead,"_dead","")),
         QC=(suppressWarnings(as.numeric(CH4_LM.r2))>=0.7 | suppressWarnings(as.numeric(CH4_HM.r2))>=0.7) & CH4_best.flux>=-100 & CH4_best.flux<=200)
write.csv(f %>% select(-key),"data/processed/flux/untagged_rescue/untagged_monthly_fluxes.csv",row.names=FALSE)
trees<-f %>% group_by(tree_id,site,Plot_Type,species,dead,dbh=Dstem) %>%
  summarise(n_meas=n(),n_qc=sum(QC,na.rm=TRUE),
            CH4_median_qc=round(median(CH4_best.flux[QC],na.rm=TRUE),3),
            CH4_mean_qc=round(mean(CH4_best.flux[QC],na.rm=TRUE),3),.groups="drop") %>% arrange(Plot_Type,species)
write.csv(trees,"data/processed/flux/untagged_rescue/untagged_monthly_trees.csv",row.names=FALSE)
cat(sprintf("\n=== %d untagged monthly stems (%d dead + %d live) ===\n",
    nrow(trees),sum(trees$dead),sum(!trees$dead)))
print(as.data.frame(trees))
cat(sprintf("\nQC'd measurements: %d of %d\n",sum(f$QC,na.rm=TRUE),nrow(f)))
