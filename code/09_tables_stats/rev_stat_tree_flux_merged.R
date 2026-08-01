#!/usr/bin/env Rscript
# ==============================================================================
# rev_stat_tree_flux_merged.R   (REVISION recompute of manuscript_statistics §1)
# Recomputes the main-text tree-flux statistics on the MERGED monthly dataset
# (tagged + 7 recovered untagged/dead-snag trees), using the SAME QC filter as
# manuscript_statistics.R:
#     (CO2_LM.r2|CO2_HM.r2|CH4_LM.r2.x|CH4_HM.r2.x) >= 0.7  &  CH4 in [-100,200]
# Reports ORIGINAL vs MERGED side by side so the revised numbers are traceable,
# and adds the dead-vs-live breakdown (new to the revision). Does NOT edit the
# original pipeline script. Writes outputs/revision/tree_flux_merged_stats.txt.
# ==============================================================================
suppressPackageStartupMessages(library(dplyr)); options(warn=-1)
num<-function(x) suppressWarnings(as.numeric(x)); NMOL_TO_UG<-57.744
qc<-function(d) d %>% filter((num(CO2_LM.r2)>=0.7|num(CO2_HM.r2)>=0.7|num(CH4_LM.r2.x)>=0.7|num(CH4_HM.r2.x)>=0.7),
                             num(CH4_best.flux.x)>=-100, num(CH4_best.flux.x)<=200) %>%
  mutate(pos=case_when(Plot.Letter%in%c("WD","WS")~"W", Plot.Letter=="I"~"I", Plot.Letter=="U"~"U", TRUE~Plot.Letter),
         f=num(CH4_best.flux.x))
fd<-"data/processed/flux"
orig<-qc(read.csv(file.path(fd,"semirigid_tree_final_complete_dataset.csv")))
merg<-qc(read.csv(file.path(fd,"semirigid_tree_final_complete_dataset_with_untagged.csv")))
if(!"dead"%in%names(merg)) merg$dead<-FALSE; merg$dead[is.na(merg$dead)]<-FALSE

out<-file("outputs/audit/tree_flux_merged_stats.txt","w"); P<-function(...) {s<-sprintf(...); cat(s,"\n"); writeLines(s,out)}
P("=== TREE-FLUX STATS: manuscript-statistics.R §1 recomputed on MERGED monthly ===")
P("QC filter identical to main text. ORIGINAL = tagged only; MERGED = + 7 untagged/dead-snag trees.\n")
P("%-38s %10s %10s %8s","metric","original","merged","delta")
P("%-38s %10d %10d %+8d","QC tree-flux measurements",nrow(orig),nrow(merg),nrow(merg)-nrow(orig))
P("%-38s %10d %10d %+8d","unique trees (QC)",length(unique(orig$Plot.Tag)),length(unique(merg$Plot.Tag)),
  length(unique(merg$Plot.Tag))-length(unique(orig$Plot.Tag)))

P("\n-- Tree CH4 flux by landscape position (merged; ug CH4 m-2 hr-1 in [] ) --")
for(p in c("U","I","W")){
  do<-orig%>%filter(pos==p); dm<-merg%>%filter(pos==p)
  lab<-c(U="Upland",I="Intermediate",W="Transitional wetland")[p]
  P("  %-20s n: %d->%d | mean %.3f -> %.3f nmol  [%.1f -> %.1f ug] | %%pos %.1f -> %.1f",
    lab, nrow(do), nrow(dm), mean(do$f), mean(dm$f), mean(do$f)*NMOL_TO_UG, mean(dm$f)*NMOL_TO_UG,
    100*sum(do$f>0)/nrow(do), 100*sum(dm$f>0)/nrow(dm))
}

P("\n-- Dead vs live (merged, QC; measurement-level means) --")
for(s in c(FALSE,TRUE)){
  d<-merg%>%filter(dead==s)
  P("  %-5s n=%d (%d trees) | mean %.3f nmol [%.1f ug] | median %.4f | %%pos %.1f",
    ifelse(s,"DEAD","LIVE"), nrow(d), length(unique(d$Plot.Tag)), mean(d$f), mean(d$f)*NMOL_TO_UG,
    median(d$f), 100*sum(d$f>0)/nrow(d))
}
P("\n(NOTE: formal dead-vs-live inference is tree-level, controlling for species/plot;")
P(" see dead-vs-live analysis. Above is descriptive only. All fold-in adds are pre-QC 45 meas / 7 trees.)")
close(out); cat("\nwrote outputs/revision/tree_flux_merged_stats.txt\n")
