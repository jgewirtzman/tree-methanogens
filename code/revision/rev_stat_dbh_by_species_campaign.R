#!/usr/bin/env Rscript
# ==============================================================================
# rev_stat_dbh_by_species_campaign.R  (SI Table — R3 L286)
# "Number of trees of each species and their DBH (mean +/- sd) in each location."
# Tree SETS mirror rev_stat_campaign_counts EXACTLY (unique trees with valid CH4
# flux), so per-campaign totals == campaign_counts (41 monthly / 150 height / 332
# in 2023). DBH is joined onto those canonical trees.
#   * 2020-2021 monthly  -> unique Plot.Tag w/ flux; grouped by moisture class
#                           (Plot.Letter U/I/WD/WS -> Upland/Intermediate/Wetland);
#                           species+DBH from YM_trees_measured (D_stem = breast-ht dia, cm)
#   * 2021 summer height -> unique (plot,tree_id) w/ flux; DBH from merged 2021 table;
#                           ungrouped (adjacent to inventory plot; no destructive sampling in plot)
#   * 2023 cross-species -> unique Tree.Tag w/ flux; ungrouped (in inventory plot)
# Writes outputs/revision/dbh_by_species_campaign.csv (+ .png).
# ==============================================================================
suppressPackageStartupMessages({library(tidyverse);library(gridExtra);library(grid)}); options(warn=-1)
num<-function(x) suppressWarnings(as.numeric(x)); uq<-function(x) length(unique(x[!is.na(x)]))
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",BEPO="Betula sp.",CALA="Carya sp.",
  CAOV="Carya ovata",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",KALA="Kalmia latifolia",
  PIST="Pinus strobus",PRSE="Prunus serotina",QUAL="Quercus alba",QURU="Quercus rubra",
  QUVE="Quercus velutina",SAAL="Sassafras albidum",TSCA="Tsuga canadensis")
# BEPO/CALA are in-plot 2023 IDs absent from the census -> reported at genus level (uncertain sp.)
# confirmed data-entry code corrections (not real distinct species)
corr <- c(TSLA="TSCA")   # TSLA typo for Tsuga canadensis (2021 aux); BEPO/CALA left pending Jon
lab<-function(code){code<-toupper(trimws(code)); code<-ifelse(code %in% names(corr), corr[code], code)
  ifelse(code %in% names(sp_map), sp_map[code], code)}
fd<-"data/processed/flux"

## 2020-2021 monthly — unique Plot.Tag with valid flux (== campaign_counts 41)
srt<-read.csv(file.path(fd,"semirigid_tree_final_complete_dataset.csv"))
mon<-data.frame(tag=as.character(srt$Plot.Tag), pos=as.character(srt$Plot.Letter),
                ok=!is.na(srt$CH4_best.flux.x), stringsAsFactors=FALSE)
mon<-mon[mon$ok,]; mon<-mon[!duplicated(mon$tag),]
ymt<-read.csv("data/raw/field_data/static_chamber_field/YM_trees_measured.csv.csv",check.names=FALSE)
sp_lut<-setNames(as.character(ymt$Species),as.character(ymt$Label)); dbh_lut<-setNames(num(ymt$D_stem),as.character(ymt$Label))
y23<-read.csv(file.path(fd,"methanogen_tree_flux_complete_dataset.csv"),check.names=FALSE)
sp23<-setNames(as.character(y23$`Species Code`),as.character(y23$`Tree Tag`))
mon$code<-ifelse(mon$tag %in% names(sp_lut), sp_lut[mon$tag], ifelse(mon$tag %in% names(sp23), sp23[mon$tag], NA))
mon$dbh<-dbh_lut[mon$tag]
# gap-fill: a few monthly trees lack a YM D_stem but are the SAME physical tree
# re-measured in 2023 (shared numeric tag) -> use the 2023 DBH (static over ~2 yr)
dbh23<-tapply(num(y23$`DBH (cm)`), as.character(y23$`Tree Tag`), function(v) mean(v,na.rm=TRUE))
need<-is.na(mon$dbh) & mon$tag %in% names(dbh23)
mon$dbh[need]<-dbh23[mon$tag[need]]
mon$location<-dplyr::recode(mon$pos, U="Upland", I="Intermediate", WD="Wetland", WS="Wetland", .default=mon$pos)
t1<-transmute(mon, campaign="2020-2021 monthly", location, code, dbh)

## 2021 summer height — unique (plot,tree_id) with valid flux (== campaign_counts 150)
ch4<-read.csv(file.path(fd,"CH4_best_flux_lgr_results.csv")); aux<-read.csv(file.path(fd,"goflux_auxfile.csv"))
hd<-merge(ch4, aux[,c("UniqueID","measurement_height","tree_id","species","plot")], by="UniqueID", all.x=TRUE)
hd<-hd[!is.na(hd$best.flux),]; hd<-hd[!duplicated(paste(hd$plot,hd$tree_id)),]
mg<-read.csv("data/processed/integrated/merged_tree_dataset_final.csv",check.names=FALSE)
# reconcile 2021 tree_id variants (aux SM1.. vs merged ab1/sm6/bo1..) via the ID crosswalk
norm<-function(x) toupper(trimws(as.character(x)))
M<-read.csv("data/processed/tree_data/tree_id_comprehensive_mapping.csv",check.names=FALSE)
vcols<-grep("^name_in_|^variant_|^name_variants$|primary_id|Tree_ID_normalized",names(M),value=TRUE)
prim<-if("primary_id" %in% names(M)) M$primary_id else M$Tree_ID_normalized
lut<-c(); for(cc in vcols){v<-norm(M[[cc]]); okk<-v!="" & !is.na(v); lut[v[okk]]<-norm(prim[okk])}
canon<-function(x){x<-norm(x); ifelse(x %in% names(lut), lut[x], x)}
dbh_by_cid<-tapply(num(mg$dbh), canon(mg$tree_id), function(v) if(all(is.na(v))) NA_real_ else mean(v,na.rm=TRUE))
hd$dbh<-as.numeric(dbh_by_cid[canon(hd$tree_id)])
t2<-transmute(hd, campaign="2021 summer (height)", location="adjacent to inventory plot", code=species, dbh)

## 2023 cross-species — unique tree with valid flux. Untagged trees share a non-unique
## "untagged"/"untagged, forked" tag -> key those by UniqueID so each is a distinct tree.
o23<-!is.na(y23$CH4_best.flux)
tg23<-trimws(as.character(y23$`Tree Tag`)); tk23<-ifelse(grepl("^[0-9]+$",tg23), tg23, as.character(y23$UniqueID))
t23<-y23[o23,]; t23$tk<-tk23[o23]; t23<-t23[!duplicated(t23$tk),]
t3<-transmute(t23, campaign="2023 cross-species", location="in inventory plot", code=`Species Code`, dbh=num(`DBH (cm)`))

allt<-bind_rows(t1,t2,t3) %>% filter(!is.na(code), toupper(trimws(code))!="NA", code!="") %>% mutate(Species=lab(code))
cat(sprintf("tree-set check vs campaign_counts:  monthly %d (want 41)  |  height %d (want 150)  |  2023 %d (want 336, untagged-corrected)\n",
    nrow(t1), nrow(t2), nrow(t3)))
cat("DBH available:", sprintf("monthly %d/%d, height %d/%d, 2023 %d/%d\n",
    sum(!is.na(t1$dbh)),nrow(t1), sum(!is.na(t2$dbh)),nrow(t2), sum(!is.na(t3$dbh)),nrow(t3)))

tab<-allt %>% group_by(Campaign=campaign,Location=location,Species) %>%
  summarise(n=n(), m=mean(dbh,na.rm=TRUE), s=sd(dbh,na.rm=TRUE), .groups="drop") %>%
  mutate(`DBH cm (mean +/- sd)`=ifelse(is.nan(m),"-",ifelse(is.na(s)|n<2, sprintf("%.1f",m), sprintf("%.1f +/- %.1f",m,s)))) %>%
  select(Campaign,Location,Species,n,`DBH cm (mean +/- sd)`) %>% arrange(Campaign,Location,Species)
write.csv(tab,"outputs/revision/dbh_by_species_campaign.csv",row.names=FALSE)
g<-tableGrob(tab,rows=NULL,theme=ttheme_minimal(base_size=8,colhead=list(fg_params=list(fontface="bold"))))
ggplot2::ggsave("outputs/revision/dbh_by_species_campaign.png",g,width=8,height=0.22*nrow(tab)+0.6,dpi=200,bg="white",limitsize=FALSE)
cat("wrote dbh_by_species_campaign.{csv,png}; rows:",nrow(tab),"\n")
