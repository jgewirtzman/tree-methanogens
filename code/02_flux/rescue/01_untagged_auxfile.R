#!/usr/bin/env Rscript
# ==============================================================================
# 01_untagged_auxfile.R  ("rescue" — monthly untagged tree fluxes)
# The 2020-21 monthly campaign fluxed ~7 untagged/dead-snag trees, but they were
# dropped at 02_join_flux_geometry (non-numeric Plot Tag -> no geometry join),
# BEFORE flux calc. Their measurements live in treeflux_total.xlsx and their
# chamber geometry lives in the Tree_fluxes_*.xlsx sheets (labelled by site+species).
# This builds a goFlux auxfile for the 49 untagged measurements so the standard
# click+goFlux workflow (code/tools/click_untagged_windows.R) can compute their fluxes.
# Geometry: the rule of code/lib/stem_chamber_geometry.R (the chamber sheet's row for
# that stem x month x chamber, else the formula on the chamber used that day), as for
# the tagged stems. (Until 2026-10-16: the sheet's uncorrected Sc/Vc of any row of the
# same position x species, falling back to another month and an arbitrary row.)
# Writes data/processed/flux/untagged_rescue/untagged_auxfile.csv
# ==============================================================================
suppressMessages({library(readxl);library(dplyr);library(tidyr)}); options(warn=-1)
num<-function(x) suppressWarnings(as.numeric(x))

## ---- untagged stems on the chamber sheets (labelled by position + species) ----
source("code/lib/stem_chamber_geometry.R")
rows<-semirigid_sheet_rows(".")
dims<-semirigid_chamber_dims(rows)
geo<-rows[grepl("untag|no tag|dead|upland|wetland|intermediate",rows$sample,ignore.case=TRUE),]
geo$site<-case_when(grepl("upland|no tag",geo$sample,ignore.case=TRUE)~"Upland",
  grepl("wetland",geo$sample,ignore.case=TRUE)~"Wetland", grepl("intermediate",geo$sample,ignore.case=TRUE)~"Intermediate")
geo$sp<-toupper(sub(".*?(ACRU|BEAL|BELE|FRAM|TSCA|BEPA|QURU).*","\\1",geo$sample,ignore.case=TRUE))
geo$sp[grepl("no tag",geo$sample,ignore.case=TRUE)]<-"BELE"
geo<-geo %>% filter(!is.na(site))
geo$label<-geo$sample; geo$sample<-paste(geo$site,geo$sp)            # one untagged stem per position x species
stopifnot(all(tapply(geo$D,geo$sample,function(d) diff(range(d))==0)))   # one stem diameter per stem

## ---- untagged measurements (treeflux_total); species from the field NOTES ----
d<-suppressMessages(read_excel("data/raw/field_data/ipad_data/treeflux_total.xlsx",sheet="Sheet1"))
d$PT<-as.character(d[["Plot Tag"]]); d$PL<-as.character(d[["Plot Letter"]]); d$CH<-as.character(d$Chamber); d$NO<-as.character(d$Notes)
spof<-function(s){ s<-toupper(s); if(grepl("TUSCAN|TSCA|TSUC",s))"TSCA" else if(grepl("BEAL",s))"BEAL" else if(grepl("BELE",s))"BELE" else if(grepl("ACRU",s))"ACRU" else if(grepl("FRAM",s))"FRAM" else NA_character_ }
u<-d[grepl("[a-zA-Z]",d$PT)&!grepl("^[0-9]+$",trimws(d$PT)),] %>%
  transmute(Date=as.Date(Date), time=num(`Time of sampling`), CH, PL, NO=paste(PT,NO),
            site=case_when(PL=="U"~"Upland",PL=="I"~"Intermediate",grepl("^W",PL)~"Wetland"),
            month=format(Date,"%Y-%m"), sp=vapply(NO,spof,character(1)),
            chamber=ifelse(CH %in% c("na","NA","") | is.na(CH),NA,sub("[.]0$","",CH)))
u<-u[!is.na(u$site)&!is.na(u$time),]
# 2 rows have no species note: Upland ch1 -> BELE; Intermediate -> BELE (TSCA is the noted one that day)
u$sp[is.na(u$sp) & u$site=="Upland"]<-"BELE"; u$sp[is.na(u$sp) & u$site=="Intermediate"]<-"BELE"

## ---- geometry: the chamber used that day (code/lib/stem_chamber_geometry.R) ----
mc<-split(sub("[.]0$","",as.character(d$CH)), format(as.Date(d$Date),"%Y-%m"))
m<-cbind(u, semirigid_deployment_geometry(paste(u$site,u$sp), u$month, ifelse(is.na(u$chamber),"na",u$chamber),
                                          NA_real_, geo, dims, mc))
stopifnot(!anyNA(m$Sc_cm2))
m$Sample<-geo$label[match(paste(m$site,m$sp),geo$sample)]; m$Dstem<-m$D_stem_cm
# dead = marked dead in ANY source (verified); only Upland BELE is never marked dead
m$dead <- !(m$site=="Upland" & m$sp=="BELE")

## ---- build goFlux auxfile (mirror 03_prep_tree_auxfile derivations) ----
# Tcham: station air temperature of the hour containing the start, in true time
# (code/lib/weather_time.R, the rule of code/tables/06_weather.R)
source("code/lib/weather_time.R")
wx<-read_station_air("data/raw/weather/ymf_clean_sorted.csv")
hh<-floor(m$time/100); mm<-floor(m$time)%%100
m$start.time<-as.POSIXct(paste(m$Date,sprintf("%02d:%02d:00",hh,mm)),format="%Y-%m-%d %H:%M:%S",tz="UTC")
m$Area<-m$Sc_cm2
source("code/lib/chamber_constants.R")  # ANALYZER_VOLUME_CM3 (was a literal 70)
m$Vtot<-(m$Vc_cm3 + pi*(1/16)^2*12*12*16.387 + ANALYZER_VOLUME_CM3)/1000
m$Tcham<-station_air_at(m$start.time,wx); m$Pcham<-101.325
m$UniqueID<-paste("UNTAG",m$site,m$sp,format(m$Date,"%Y%m%d"),ifelse(is.na(m$chamber),"na",m$chamber),sep="_")
m$UniqueID<-ave(m$UniqueID,m$UniqueID,FUN=function(x) if(length(x)==1) x else paste(x,seq_along(x),sep="_"))
aux<-m %>% filter(!is.na(start.time),!is.na(Area),!is.na(Vtot)) %>%
  transmute(UniqueID, start.time, start.time_formatted=format(start.time,"%Y-%m-%d %H:%M:%S"),
            Area, Vtot, Tcham, Pcham, site, species=sp, dead, Sample, Dstem, obs.length=600,
            geometry_source, chamber_dims, chamber_volume_cm3=Vc_cm3)
dir.create("data/processed/flux/untagged_rescue", showWarnings = FALSE, recursive = TRUE)
write.csv(aux,"data/processed/flux/untagged_rescue/untagged_auxfile.csv",row.names=FALSE)
cat("untagged measurements:",nrow(u)," -> auxfile rows (with geometry+time):",nrow(aux),"\n")
cat("trees covered:\n"); print(as.data.frame(aux %>% count(site,species,name="n_meas")))
cat("geometry by source:\n"); print(as.data.frame(aux %>% count(geometry_source)))
cat("Sc range (cm2):",paste(round(range(aux$Area,na.rm=TRUE)),collapse="-"),"| wrote data/processed/flux/untagged_rescue/untagged_auxfile.csv\n")
