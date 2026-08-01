#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS_black-oak-cross-sections.R  (SI)
# Field photo plate: felled black-oak stem cross-sections by height (base -> top).
# Documents the visible decay cone (basal rot cavity + discolouration near the
# base, sound wood by ~6 m) that underlies the within-tree profiles in Fig 7 (the
# mcrA/CH4/flux peak at 4-6 m is the leading edge of this cone). Photos are field
# shots (ruler for scale); centre-cropped to square and labelled by height.
# Input photos: outputs/figures/black_oak_cross_section_photos/{height}.jpg
# Writes outputs/figures/generated/figS_black_oak_cross_sections.png
# ==============================================================================
suppressMessages({library(jpeg);library(grid);library(gridExtra)});options(warn=-1)
dir<-"outputs/figures/black_oak_cross_section_photos"
ord<-c("50cm","125cm","2m","4m","6m","8m","10m")
labs<-c("0.5 m","1.25 m","2 m","4 m","6 m","8 m","10 m")
grobs<-Map(function(f,l){
  img<-readJPEG(file.path(dir,paste0(f,".jpg")))
  H<-dim(img)[1];W<-dim(img)[2];s<-min(H,W);x0<-floor((W-s)/2);y0<-floor((H-s)/2)
  crop<-img[(y0+1):(y0+s),(x0+1):(x0+s),,drop=FALSE]
  arrangeGrob(rasterGrob(crop,interpolate=TRUE),
    top=textGrob(l,gp=gpar(fontsize=17,fontface="bold"),vjust=.3),padding=unit(2,"mm"))
},ord,labs)
# title/caption intentionally omitted -> lives in the manuscript figure caption
plate<-arrangeGrob(grobs=grobs,nrow=1,padding=unit(3,"mm"))
ggplot2::ggsave("outputs/figures/generated/figS_black_oak_cross_sections.png",plate,width=15.5,height=3.1,dpi=220,bg="white")
cat("wrote outputs/figures/generated/figS_black_oak_cross_sections.png\n")
