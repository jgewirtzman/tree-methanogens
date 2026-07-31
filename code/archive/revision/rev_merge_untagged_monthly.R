#!/usr/bin/env Rscript
# ==============================================================================
# rev_merge_untagged_monthly.R
# Folds the recovered untagged/dead-snag monthly measurements into the monthly
# tree-flux dataset, producing a REVISION merged file (original left untouched):
#   data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv
# Schema-matched to the original (Plot.Tag, Plot.Letter, Date, CH4_best.flux.x,
# CH4_LM.r2.x, CH4_HM.r2.x, CO2_LM.r2/HM.r2) + a new `dead` flag. The tagged trees
# 4156 & 5091 (flagged dead in field notes) get dead=TRUE too.
# Revision Fig 1 + manuscript-stat recompute read this file.
# ==============================================================================
suppressMessages(library(dplyr)); options(warn=-1); num<-function(x)suppressWarnings(as.numeric(x))
srt<-read.csv("data/processed/flux/semirigid_tree_final_complete_dataset.csv")   # check.names=TRUE
# dead tagged trees from field notes
suppressMessages(library(readxl))
tf<-read_excel("data/raw/field_data/ipad_data/treeflux_total.xlsx",sheet="Sheet1"); tf$PT<-trimws(as.character(tf[["Plot Tag"]])); tf$NO<-as.character(tf$Notes)
dead_tags<-unique(tf$PT[grepl("^[0-9]+$",tf$PT)&grepl("dead|snag",tf$NO,ignore.case=TRUE)])
srt$dead <- trimws(as.character(srt$Plot.Tag)) %in% dead_tags
srt$Plot.Tag <- as.character(srt$Plot.Tag)                       # allow string untagged IDs

# untagged recovered measurements -> srt schema
u<-read.csv("outputs/revision/untagged_monthly_fluxes.csv",check.names=FALSE)
u$Date<-as.Date(sub(".*_(\\d{8})(_.*)?$","\\1",u$UniqueID),format="%Y%m%d")
n<-nrow(u); untag<-srt[rep(NA_integer_,n),,drop=FALSE]; rownames(untag)<-NULL   # typed NA clone
untag$Plot.Tag        <- u$tree_id                              # per-TREE id (measurements share it)
untag$Plot.Letter     <- dplyr::recode(u$Plot_Type, U="U", I="I", W="WS")
untag$Date            <- as.character(u$Date)
untag$CH4_best.flux.x <- u$CH4_best.flux ; untag$CH4_best.flux.y <- u$CH4_best.flux
untag$CH4_LM.r2.x     <- num(u$CH4_LM.r2); untag$CH4_HM.r2.x     <- num(u$CH4_HM.r2)
untag$CH4_LM.r2.y     <- num(u$CH4_LM.r2); untag$CH4_HM.r2.y     <- num(u$CH4_HM.r2)
untag$dead            <- u$dead

merged<-bind_rows(srt, untag)
write.csv(merged,"data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",row.names=FALSE)
cat("merged monthly dataset:",nrow(srt),"tagged +",nrow(untag),"untagged =",nrow(merged),"rows\n")
cat("dead measurements:",sum(merged$dead,na.rm=TRUE)," | dead trees:",length(unique(merged$Plot.Tag[merged$dead])),"\n")
cat("wrote data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv\n")
