#!/usr/bin/env Rscript
# ==============================================================================
# rev_rescue_untagged_finalize.R
# Attaches the NOTES-based species + verified live/dead + corrected geometry to
# the clicked measurements (untagged_manID.rds), re-keying on start.time (stable,
# species-independent), overwrites the chamber geometry in manID, and recomputes
# goFlux so fluxes use the correct Sc. Writes:
#   outputs/revision/untagged_monthly_fluxes.csv   (per measurement, canonical)
#   outputs/revision/untagged_monthly_trees.csv    (per tree; 7 trees, 6 dead + 1 live)
# Prereq: rev_rescue_untagged_prep.R (writes untagged_auxfile.csv with species/dead/geometry).
# ==============================================================================
suppressMessages({library(goFlux);library(dplyr)}); options(warn=-1)
manID<-readRDS("outputs/revision/untagged_manID.rds")
aux<-read.csv("outputs/revision/untagged_auxfile.csv",check.names=FALSE)
aux$key<-aux$start.time_formatted
stopifnot(!any(duplicated(aux$key)))                      # start.time is a unique measurement key
manID$key<-format(as.POSIXct(manID$start.time_formatted,tz="UTC"),"%Y-%m-%d %H:%M:%S")
mp<-aux %>% select(key, Area_c=Area, Vtot_c=Vtot, site, species, dead, Sample, Dstem)
manID<-manID %>% select(-any_of(c("Area","Vtot","site","species","Sample","Dstem"))) %>%
  left_join(mp,by="key") %>% rename(Area=Area_c, Vtot=Vtot_c)     # overwrite geometry w/ corrected
cat("re-keyed", length(unique(manID$key)), "measurements to notes-based species/geometry\n")

ch4<-goFlux(manID,"CH4dry_ppb"); ch4b<-best.flux(ch4); co2b<-best.flux(goFlux(manID,"CO2dry_ppm"))
meta<-manID %>% distinct(UniqueID,key,site,species,dead,Sample,Dstem,Area,Vtot)
f<-meta %>%
  left_join(ch4 %>% transmute(UniqueID,CH4_LM.r2=LM.r2,CH4_HM.r2=HM.r2),by="UniqueID") %>%
  left_join(ch4b %>% transmute(UniqueID,CH4_best.flux=best.flux,CH4_model=model,CH4_quality=quality.check),by="UniqueID") %>%
  left_join(co2b %>% transmute(UniqueID,CO2_best.flux=best.flux),by="UniqueID") %>%
  mutate(Plot_Type=case_when(site=="Upland"~"U",site=="Intermediate"~"I",site=="Wetland"~"W"),
         tree_id=paste0("UNTAG_",substr(site,1,3),"_",species,ifelse(dead,"_dead","")),
         QC=(suppressWarnings(as.numeric(CH4_LM.r2))>=0.7 | suppressWarnings(as.numeric(CH4_HM.r2))>=0.7) & CH4_best.flux>=-100 & CH4_best.flux<=200)
write.csv(f %>% select(-key),"outputs/revision/untagged_monthly_fluxes.csv",row.names=FALSE)
trees<-f %>% group_by(tree_id,site,Plot_Type,species,dead,dbh=Dstem) %>%
  summarise(n_meas=n(),n_qc=sum(QC,na.rm=TRUE),
            CH4_median_qc=round(median(CH4_best.flux[QC],na.rm=TRUE),3),
            CH4_mean_qc=round(mean(CH4_best.flux[QC],na.rm=TRUE),3),.groups="drop") %>% arrange(Plot_Type,species)
write.csv(trees,"outputs/revision/untagged_monthly_trees.csv",row.names=FALSE)
cat(sprintf("\n=== %d untagged monthly trees (%d dead + %d live) -> fold monthly 41 to %d ===\n",
    nrow(trees),sum(trees$dead),sum(!trees$dead),41+nrow(trees)))
print(as.data.frame(trees))
cat(sprintf("\nQC'd measurements: %d of %d\n",sum(f$QC,na.rm=TRUE),nrow(f)))
