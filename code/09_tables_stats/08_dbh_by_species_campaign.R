#!/usr/bin/env Rscript
# ==============================================================================
# 08_dbh_by_species_campaign.R  (SI Table)
# "Number of trees of each species and their DBH (mean +/- sd) in each location."
# Tree sets mirror 02_campaign_counts.R: every tree (tree_id) with a deployment in
# data/tables/flux_stem.csv for the campaign. Species: trees.csv species_code, named from
# species.csv. Status: Dead where any of the tree's deployments in that campaign is on a
# dead stem (flux_stem.dead_stem). DBH, from the canonical records only:
#   * 2020-2021 monthly  -> trees.csv dbh_cm (monthly survey tree list); grouped by
#                           landscape position (trees.csv: U/I/WD/WS -> Upland/Intermediate/Wetland)
#   * 2021 summer height -> trees.csv dbh_cm (2021 datasheets / DBH list); ungrouped
#   * 2023 cross-species -> the diameter written on the 2023 field sheet for the deployment
#                           (flux_stem.dbh_recorded_cm; mean per tree); ungrouped
#   trees.csv carries measured DBH only; trees without a measurement (dbh_cm NA) are left
#   out of the means; the count of trees is unaffected.
# Writes outputs/data/dbh_by_species_campaign.csv (+ .png).
# ==============================================================================
suppressPackageStartupMessages({library(tidyverse);library(gridExtra);library(grid)}); options(warn=-1)
source("code/lib/tables_read.R")
trees <- trees_table()
measured_dbh <- function(id) trees$dbh_cm[match(id, trees$tree_id)]
per_tree <- function(campaign) flux_stem_campaign(campaign) %>% group_by(tree_id) %>%
  summarise(dead = any(dead_stem), dbh_rec = if (all(is.na(dbh_recorded_cm))) NA_real_ else mean(dbh_recorded_cm, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(code = tree_species(tree_id), status = ifelse(dead, "Dead", "Live"))

## 2020-2021 monthly (tagged and recovered untagged stems)
t1 <- per_tree("monthly_2020_2021") %>%
  transmute(campaign = "2020-2021 monthly",
            location = dplyr::recode(trees$landscape_position[match(tree_id, trees$tree_id)],
                                     U = "Upland", I = "Intermediate", WD = "Wetland", WS = "Wetland"),
            code, dbh = measured_dbh(tree_id), status)
## 2021 summer height
t2 <- per_tree("2021_multiheight") %>%
  transmute(campaign = "2021 summer (height)", location = "adjacent to inventory plot", code,
            dbh = measured_dbh(tree_id), status)
## 2023 cross-species
t3 <- per_tree("2023_cross_species") %>%
  transmute(campaign = "2023 cross-species", location = "in inventory plot", code, dbh = dbh_rec, status)

allt <- bind_rows(t1, t2, t3) %>% mutate(Species = species_name(code))
stopifnot(!anyNA(allt$Species))
cat(sprintf("trees per campaign (as campaign_counts):  monthly %d  |  height %d  |  2023 %d\n",
    nrow(t1), nrow(t2), nrow(t3)))
cat("DBH available:", sprintf("monthly %d/%d, height %d/%d, 2023 %d/%d\n",
    sum(!is.na(t1$dbh)),nrow(t1), sum(!is.na(t2$dbh)),nrow(t2), sum(!is.na(t3$dbh)),nrow(t3)))

tab<-allt %>% group_by(Campaign=campaign,Location=location,Species,Status=status) %>%
  summarise(n=n(), m=mean(dbh,na.rm=TRUE), s=sd(dbh,na.rm=TRUE), .groups="drop") %>%
  mutate(`DBH cm (mean +/- sd)`=ifelse(is.nan(m),"-",ifelse(is.na(s)|n<2, sprintf("%.1f",m), sprintf("%.1f +/- %.1f",m,s)))) %>%
  select(Campaign,Location,Species,Status,n,`DBH cm (mean +/- sd)`) %>% arrange(Campaign,Location,Species,Status)
write.csv(tab,"outputs/data/dbh_by_species_campaign.csv",row.names=FALSE)
g<-tableGrob(tab,rows=NULL,theme=ttheme_minimal(base_size=8,colhead=list(fg_params=list(fontface="bold"))))
ggplot2::ggsave("outputs/figures/generated/dbh_by_species_campaign.png",g,width=8,height=0.22*nrow(tab)+0.6,dpi=200,bg="white",limitsize=FALSE)
cat("wrote dbh_by_species_campaign.{csv,png}; rows:",nrow(tab),"\n")
