#!/usr/bin/env Rscript
# ==============================================================================
# rev_rescue_untagged_prep.R  (revision "rescue" — monthly untagged tree fluxes)
# The 2020-21 monthly campaign fluxed ~7 untagged/dead-snag trees, but they were
# dropped at 02_join_flux_geometry (non-numeric Plot Tag -> no geometry join),
# BEFORE flux calc. Their measurements live in treeflux_total.xlsx and their
# chamber geometry lives in the Tree_fluxes_*.xlsx sheets (labelled by site+species).
# This builds a goFlux auxfile for the 49 untagged measurements so the standard
# click+goFlux workflow (rev_rescue_untagged_click.R) can compute their fluxes.
# Writes outputs/revision/untagged_auxfile.csv
# ==============================================================================
suppressMessages({library(readxl);library(dplyr);library(tidyr)}); options(warn=-1)
num<-function(x) suppressWarnings(as.numeric(x))
dir<-"data/raw/flux_chamber_corrections/Chamber vol corr"
mons<-c(june="2020-06",july="2020-07",august="2020-08",september="2020-09",october="2020-10",december="2020-12")

## ---- untagged geometry per tree x month x chamber (from Tree_fluxes sheets) ----
geo<-bind_rows(lapply(names(mons),function(m){
  p<-file.path(dir,paste0("Tree_fluxes_",m,".xlsx")); if(!file.exists(p))return(NULL)
  x<-suppressMessages(read_excel(p)); x$Sample<-as.character(x$Sample)
  u<-x[grepl("untag|no tag|dead|upland|wetland|intermediate",x$Sample,ignore.case=TRUE),,drop=FALSE]; if(!nrow(u))return(NULL)
  Dst<-num(u[["D stem (cm)"]]); Dext<-num(u[["D exterior (cm)"]]); H<-num(u[["Chamber height (cm)"]]); K<-num(u[["K"]]); Vw<-num(u[["Two V wedges (cm3)"]])
  tibble(month=mons[m], Sample=u$Sample, chamber=as.character(u[["Chamber used"]]),
         Sc_cm2=num(u[["Sc (cm2)"]]), Vc_cm3=coalesce(num(u$Vc), K*(pi/4)*(Dext^2-Dst^2)*H-Vw), Dstem=Dst)
}))
geo$site<-case_when(grepl("upland|no tag",geo$Sample,ignore.case=TRUE)~"Upland",
  grepl("wetland",geo$Sample,ignore.case=TRUE)~"Wetland", grepl("intermediate",geo$Sample,ignore.case=TRUE)~"Intermediate")
geo$sp<-toupper(sub(".*?(ACRU|BEAL|BELE|FRAM|TSCA|BEPA|QURU).*","\\1",geo$Sample,ignore.case=TRUE))
geo$sp[grepl("no tag",geo$Sample,ignore.case=TRUE)]<-"BELE"
geo<-geo %>% filter(chamber!="XXXX",!is.na(site))

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

## ---- attach geometry for the (now known) site x species x month x chamber tree ----
assign_geo<-function(r){
  cand<-geo[geo$site==r$site & geo$sp==r$sp,]
  if(!nrow(cand)) return(geo[1,][NA,])
  c2<-cand[cand$month==r$month,]; if(nrow(c2)) cand<-c2
  if(!is.na(r$chamber) && any(cand$chamber==r$chamber)) cand<-cand[cand$chamber==r$chamber,]
  cand[1,]
}
rows<-lapply(seq_len(nrow(u)),function(i){a<-assign_geo(u[i,]); cbind(u[i,],a[,c("Sample","Sc_cm2","Vc_cm3","Dstem")])})
m<-bind_rows(rows)
# dead = marked dead in ANY source (verified); only Upland BELE is never marked dead
m$dead <- !(m$site=="Upland" & m$sp=="BELE")

## ---- build goFlux auxfile (mirror 03_prep_tree_auxfile derivations) ----
wx<-read.csv("data/raw/weather/ymf_clean_sorted.csv"); wx$TS<-as.POSIXct(wx$TIMESTAMP,tz="UTC")
nearestT<-function(t){ i<-which.min(abs(as.numeric(wx$TS-t))); if(length(i)) num(wx$Tair_Avg[i]) else NA }
hh<-floor(m$time/100); mm<-floor(m$time)%%100
m$start.time<-as.POSIXct(paste(m$Date,sprintf("%02d:%02d:00",hh,mm)),format="%Y-%m-%d %H:%M:%S",tz="UTC")
m$Area<-m$Sc_cm2
m$Vtot<-(m$Vc_cm3 + pi*(1/16)^2*12*12*16.387 + 70)/1000
m$Tcham<-vapply(m$start.time,nearestT,numeric(1)); m$Pcham<-101.325
m$UniqueID<-paste("UNTAG",m$site,m$sp,format(m$Date,"%Y%m%d"),ifelse(is.na(m$chamber),"na",m$chamber),sep="_")
m$UniqueID<-ave(m$UniqueID,m$UniqueID,FUN=function(x) if(length(x)==1) x else paste(x,seq_along(x),sep="_"))
aux<-m %>% filter(!is.na(start.time),!is.na(Area),!is.na(Vtot)) %>%
  transmute(UniqueID, start.time, start.time_formatted=format(start.time,"%Y-%m-%d %H:%M:%S"),
            Area, Vtot, Tcham, Pcham, site, species=sp, dead, Sample, Dstem, obs.length=600)
write.csv(aux,"outputs/revision/untagged_auxfile.csv",row.names=FALSE)
cat("untagged measurements:",nrow(u)," -> auxfile rows (with geometry+time):",nrow(aux),"\n")
cat("trees covered:\n"); print(as.data.frame(aux %>% count(site,species,name="n_meas")))
cat("Sc range (cm2):",paste(round(range(aux$Area,na.rm=TRUE)),collapse="-"),"| wrote outputs/revision/untagged_auxfile.csv\n")
