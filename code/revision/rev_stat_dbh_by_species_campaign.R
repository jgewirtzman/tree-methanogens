#!/usr/bin/env Rscript
# ==============================================================================
# rev_stat_dbh_by_species_campaign.R  (SI Table — R3 L286)
# "A table of the number of trees of each species and their DBH (mean +/- sd) in
# each location." Built per CAMPAIGN from the canonical campaign datasets (same
# sources as rev_stat_campaign_counts), NOT the pooled tree_properties.csv:
#   * 2020-2021 monthly  -> grouped by moisture class (Upland / Intermediate / Wetland)
#   * 2021 summer height -> ungrouped (trees immediately ADJACENT to the inventory
#                           plot; destructive sampling avoided inside the plot)
#   * 2023 cross-species -> ungrouped (all trees IN the inventory plot)
# Writes outputs/revision/dbh_by_species_campaign.csv (+ .png render).
# ==============================================================================
suppressPackageStartupMessages({library(tidyverse);library(gridExtra);library(grid)}); options(warn=-1)
num<-function(x) suppressWarnings(as.numeric(x))
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",BEPO="Betula populifolia",CAOV="Carya ovata",
  FAGR="Fagus grandifolia",FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",
  PRSE="Prunus serotina",QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",
  SAAL="Sassafras albidum",TSCA="Tsuga canadensis")
lab<-function(code) ifelse(code %in% names(sp_map), sp_map[code], code)  # unmapped -> keep code

fd<-"data/processed/flux"
# 2020-2021 monthly (moisture class = Site)
ym<-read.csv("data/raw/field_data/static_chamber_field/YM_trees_measured.csv.csv",check.names=FALSE)
t1<-tibble(campaign="2020-2021 monthly", location=ym$Site, code=toupper(trimws(ym$Species)), dbh=num(ym$D_stem))
# 2021 summer height (unique trees; adjacent to inventory plot)
m<-read.csv("data/processed/integrated/merged_tree_dataset_final.csv",check.names=FALSE)
t2<-tibble(campaign="2021 summer (height)", location="adjacent to inventory plot",
           code=toupper(trimws(m$species_id)), dbh=num(m$dbh))
# 2023 cross-species (unique trees; in inventory plot)
y<-read.csv(file.path(fd,"methanogen_tree_flux_complete_dataset.csv"),check.names=FALSE)
y<-y[!duplicated(y[["Tree Tag"]]),]
t3<-tibble(campaign="2023 cross-species", location="in inventory plot",
           code=toupper(trimws(y[["Species Code"]])), dbh=num(y[["DBH (cm)"]]))

all<-bind_rows(t1,t2,t3) %>% filter(!is.na(code),code!="",code!="NA",!is.na(dbh)) %>% mutate(Species=lab(code))
tab<-all %>% group_by(Campaign=campaign,Location=location,Species) %>%
  summarise(n=n(), m=mean(dbh), s=sd(dbh), .groups="drop") %>%
  mutate(`DBH cm (mean +/- sd)`=ifelse(is.na(s)|n<2, sprintf("%.1f",m), sprintf("%.1f +/- %.1f",m,s))) %>%
  select(Campaign,Location,Species,n,`DBH cm (mean +/- sd)`) %>%
  arrange(Campaign,Location,Species)
write.csv(tab,"outputs/revision/dbh_by_species_campaign.csv",row.names=FALSE)
unmapped<-setdiff(unique(all$code[!all$code %in% names(sp_map)]),"")
cat("rows:",nrow(tab)," | unmapped species codes (kept as-is):",ifelse(length(unmapped),paste(unmapped,collapse=", "),"none"),"\n")
# render
g<-tableGrob(tab,rows=NULL,theme=ttheme_minimal(base_size=8,colhead=list(fg_params=list(fontface="bold"))))
ggplot2::ggsave("outputs/revision/dbh_by_species_campaign.png",g,width=8,height=0.22*nrow(tab)+0.6,dpi=200,bg="white",limitsize=FALSE)
cat("wrote dbh_by_species_campaign.{csv,png}\n"); print(as.data.frame(tab),row.names=FALSE)
