#!/usr/bin/env Rscript
# ==============================================================================
# fig07_decay-methanogenesis.R  (candidate expanded MAIN Fig 7 — DRAFT)
# Folds the deterioration + fungal-decay story into Fig 7:
#   TOP ROW  (population-level, a-b):
#     a  stem decay (bark loss) -> CH4 flux hump   [2023 survey]
#     b  fungal colonization (ITS load) -> mcrA hump, wood vs soil  [2021 survey]
#   BOTTOM ROW (within the felled black oak, c-f) — vertical profiles vs height:
#     c  internal CH4    d  mcrA (copies/g)    e  CH4 flux    f  fungal ITS load
#   Methanogens/gas/flux peak at 4-6 m — the leading edge of the ~4 m decay cone
#   (cavity <1 m, staining 1-2 m, clean >6 m) = intermediate/early decay, matching
#   the two population humps.
# Matches rev_fig07 felled-oak aesthetics. Data: canonical tables (code/lib/tables_read.R):
#   felled oak  internal_gas.csv (calibrated trunk gas by height), ddpcr.csv + samples.csv
#               (mcrA copies/g), qpcr.csv (ITS copies per uL extract, converted with each
#               sample's own extraction mass), flux_stem.csv (2022_felled_oak)
#   2023 survey flux_stem.csv (2023_cross_species; bark-loss score; species from trees.csv)
#   2021 survey ddpcr.csv + qpcr.csv + samples.csv per DNA sample (dna_samples())
# Writes outputs/figures/generated/fig7_decay_methanogenesis.png
# ==============================================================================
suppressMessages({library(tidyverse);library(patchwork);library(lme4);library(lmerTest)})
source("code/lib/ddpcr_constants.R"); source("code/lib/tables_read.R")
options(warn=-1)
asinh10<-function(x) asinh(x/0.1)/log(10)
# marginal R2 for a mixed model (fixed-effects variance / total)
r2m<-function(m){vf<-var(predict(m,re.form=NA));vc<-as.data.frame(VarCorr(m))
  vr<-sum(vc$vcov[vc$grp!="Residual"]);ve<-vc$vcov[vc$grp=="Residual"];as.numeric(vf/(vf+vr+ve))}
pfmt<-function(p) ifelse(p<0.001,"p < 0.001",sprintf("p = %.3f",p))
# panel stats show R2 + p only; dAIC (hump vs linear) reported in text
stat_lab<-function(R2,p) sprintf("R\u00b2 = %.2f\n%s",R2,pfmt(p))
hlab<-function(b) vapply(b,function(v) format(v,drop0trailing=TRUE,trim=TRUE),"")   # 10, 8 ... 1.25, 0.5
shared_theme <- theme_classic(base_size=12)+theme(axis.title=element_text(size=12),axis.text=element_text(size=10),
  legend.position="bottom",legend.title=element_text(size=9.5),legend.text=element_text(size=9),
  legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,0,0,0),plot.margin=margin(5,7,5,5))

## ===================== BOTTOM ROW: felled-oak vertical profiles =====================
# felled-oak trunk gas (internal_gas.csv; GC calibration in code/tables/32_internal_gas.R)
int_gas<-read_table("internal_gas")%>%filter(campaign=="2022_felled_oak")%>%
  transmute(Tree.Height=height_m,CH4_concentration=CH4_ppm,O2_concentration=O2_ppm)
hb<-sort(unique(int_gas$Tree.Height))
# felled-oak DNA samples: mcrA copies/g (dry wood) and ITS per uL extract, with height and tissue
bo<-dna_samples("2022_felled_oak")
mcra<-bo%>%filter(core_type%in%c("Heartwood","Sapwood"),!is.na(mcra_probe_g))%>%
  transmute(`Sample ID`=sample_id,Height_cm=height_cm,Component=core_type,mcrA_g=mcra_probe_g)
flux<-flux_stem_campaign("2022_felled_oak")%>%transmute(Height_m=height_cm/100,CH4_best.flux=CH4_flux_nmol_m2_s)
fdf<-flux%>%filter(!is.na(Height_m))%>%select(height=Height_m,flux=CH4_best.flux)%>%
  left_join(int_gas%>%group_by(Tree.Height)%>%summarize(ch4=mean(CH4_concentration,na.rm=TRUE),.groups="drop"),by=c("height"="Tree.Height"))
# felled-oak ITS load (heartwood), by height -> copies/g dry (same basis as mcrA panel d),
# each sample with its own extraction mass
itw<-bo%>%filter(core_type=="Heartwood",!is.na(ITS_per_ul))%>%
  transmute(height=height_cm/100,tissue="Heartwood",ITS_g=extract_copies_per_g(ITS_per_ul,extraction_mass_mg,"Wood"))

col_ht<-c(Heartwood="#a6611a",Sapwood="#dfc27d")          # Fig 4 heartwood / sapwood
col_ht_line<-c(Heartwood="#a6611a",Sapwood="#c49a3e")     # sapwood line a shade darker so it reads on white
pc<-ggplot(int_gas,aes(Tree.Height,CH4_concentration))+geom_smooth(se=FALSE,color="black")+
  geom_point(aes(fill=O2_concentration/1e6*100),size=3,shape=21,color="black",stroke=.6,alpha=.85)+shared_theme+
  scale_fill_gradient(low="#c6dbef",high="#08306b",name=expression(O[2]~"(%)"))+coord_flip()+
  ylab(expression(CH[4]~"(ppm)"))+xlab("Height (m)")+scale_x_continuous(breaks=hb,labels=hlab,minor_breaks=NULL)+
  guides(fill=guide_colorbar(barwidth=7,barheight=.5,title.position="top"))
pd<-ggplot(mcra,aes(Height_cm/100,mcrA_g))+geom_smooth(se=FALSE,aes(color=Component))+
  geom_jitter(size=3,shape=21,aes(fill=Component),color="black",stroke=.6,alpha=.85)+shared_theme+
  scale_color_manual(values=col_ht_line)+scale_fill_manual(values=col_ht)+coord_flip()+
  ylab(expression(italic(mcrA)~"(copies g"^-1*")"))+xlab("")+scale_x_continuous(breaks=hb,labels=hlab,minor_breaks=NULL)+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
# flux points in the same purple as panel a (CH4-ppm colour dropped: redundant with panel c)
pe<-ggplot(fdf,aes(height,flux))+geom_smooth(se=FALSE,color="black")+
  geom_point(size=3,shape=21,fill="#756BB1",color="black",stroke=.6,alpha=.85)+shared_theme+coord_flip()+xlab("")+
  ylab(expression(CH[4]~"flux (nmol m"^-2*" s"^-1*")"))+scale_x_continuous(breaks=hb,labels=hlab,minor_breaks=NULL)
# linear axis: makes the single 2 m heartwood spike read as the true peak (log/arcsinh
# under-weight it and the LOESS wiggles toward the ~1e4 points)
pf<-ggplot(itw,aes(height,ITS_g))+geom_smooth(se=FALSE,color="black",span=1)+
  geom_jitter(width=.15,height=0,size=3,shape=21,fill="#a6611a",color="black",stroke=.6,alpha=.85)+shared_theme+
  coord_flip()+scale_y_continuous(labels=function(b) format(b/1e9,drop0trailing=TRUE,trim=TRUE))+xlab("")+
  ylab(expression("fungal ITS ("*10^9~"copies g"^-1*")"))+theme(plot.margin=margin(5,16,5,5))+scale_x_continuous(breaks=hb,labels=hlab,minor_breaks=NULL)

## ===================== TOP ROW: population humps =====================
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-flux_stem_campaign("2023_cross_species")
y<-y[!is.na(y$CH4_flux_nmol_m2_s),];y$fx<-asinh10(y$CH4_flux_nmol_m2_s);y$sp<-as.factor(y$species_code);y$x<-ord(y$bark_missing);y<-y[!is.na(y$x),]
mq<-lmer(fx~x+I(x^2)+(1|sp),data=y);co<-summary(mq)$coefficients;pqa<-co["I(x^2)","Pr(>|t|)"]
# AIC compares fixed effects, so from ML refits (REML AICs are not comparable)
mla<-lmer(fx~x+(1|sp),data=y);dAICa<-as.numeric(AIC(update(mla,REML=FALSE))-AIC(update(mq,REML=FALSE)));R2a<-r2m(mq)
gA<-data.frame(x=seq(1,4,.05));Xa<-model.matrix(~x+I(x^2),gA);ba<-fixef(mq);Va<-as.matrix(vcov(mq));gA$fit<-as.vector(Xa%*%ba);gA$se<-sqrt(diag(Xa%*%Va%*%t(Xa)))
obrk<-c(0,0.05,0.2,1)
# stem-flux purple, matching Fig 2 (avoids implying 'wood' via brown)
pa<-ggplot(y,aes(x,fx))+geom_hline(yintercept=0,linetype=3,colour="grey70")+
  geom_jitter(width=.12,alpha=.15,size=1.1,colour="#756BB1")+
  geom_ribbon(data=gA,aes(x,ymin=fit-1.96*se,ymax=fit+1.96*se),alpha=.2,fill="#756BB1",inherit.aes=FALSE)+
  geom_line(data=gA,aes(x,fit),linewidth=1,colour="#54278f",inherit.aes=FALSE)+
  annotate("text",Inf,asinh10(0.9),hjust=1.05,vjust=1,size=3.2,lineheight=.9,label=stat_lab(R2a,pqa))+
  scale_x_continuous(breaks=1:4,labels=c("healthy","moderate","severe","dead"))+scale_y_continuous(breaks=asinh10(obrk),labels=obrk)+
  coord_cartesian(ylim=asinh10(c(-0.3,1)))+shared_theme+labs(x="bark loss (decay)",y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))

# 2021 survey, per DNA sample: mcrA copies/g (ddpcr.csv, probe assay, loose threshold);
# ITS: facility qPCR copies per uL of extract x elution / mass (extract_copies_per_g)
o<-dna_samples("2021_survey")
o$mcra_g<-o$mcra_probe_g
o$ITS_g<-extract_copies_per_g(o$ITS_per_ul,o$extraction_mass_mg,o$material)
d<-o[is.finite(o$mcra_g)&is.finite(o$ITS_g)&o$mcra_g>0&o$ITS_g>0,];d$sp<-as.factor(d$species_code)
d$X<-log10(d$ITS_g);d$Y<-log10(d$mcra_g)
fm<-function(mat,q){s<-d[d$material==mat,];s$xc<-s$X-mean(s$X);mm<-if(q)lmer(Y~xc+I(xc^2)+(1|sp),s) else lmer(Y~X+(1|sp),s)
  g<-data.frame(X=seq(min(s$X),max(s$X),length=60));g$xc<-g$X-mean(s$X);Xm<-if(q)model.matrix(~xc+I(xc^2),g) else model.matrix(~X,g)
  b<-fixef(mm);V<-as.matrix(vcov(mm));g$fit<-as.vector(Xm%*%b);g$se<-sqrt(diag(Xm%*%V%*%t(Xm)));g$material<-mat;g}
gw<-fm("Wood",TRUE);gs<-fm("Soil",FALSE)
# wood stats: quadratic (hump) vs linear
dw<-d[d$material=="Wood",];dw$xc<-dw$X-mean(dw$X)
mqb<-lmer(Y~xc+I(xc^2)+(1|sp),dw);mlb<-lmer(Y~xc+(1|sp),dw)
pqb<-summary(mqb)$coefficients[3,"Pr(>|t|)"];dAICb<-as.numeric(AIC(update(mlb,REML=FALSE))-AIC(update(mqb,REML=FALSE)));R2b<-r2m(mqb)
# paper BrBG convention: wood = brown, soil = teal (matches Fig 4/5 compartment palette)
col_ws<-c(Wood="#a6611a",Soil="#018571")
pb<-ggplot(d,aes(X,Y,colour=material))+geom_point(alpha=.4,size=1.2)+
  geom_ribbon(data=rbind(gw,gs),aes(X,ymin=fit-1.96*se,ymax=fit+1.96*se,fill=material),alpha=.15,colour=NA,inherit.aes=FALSE)+
  geom_line(data=gw,aes(X,fit),colour=col_ws["Wood"],linewidth=1,inherit.aes=FALSE)+
  geom_line(data=gs,aes(X,fit),colour=col_ws["Soil"],linewidth=1,inherit.aes=FALSE)+
  scale_colour_manual(values=col_ws,name=NULL)+scale_fill_manual(values=col_ws,guide="none")+
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))+
  annotate("text",Inf,Inf,hjust=1.05,vjust=1.2,size=3.2,lineheight=.9,label=paste0("wood:\n",stat_lab(R2b,pqb)))+
  shared_theme+
  scale_x_continuous(labels=function(b) parse(text=paste0("10^",b)))+scale_y_continuous(labels=function(b) parse(text=paste0("10^",b)))+
  labs(x=expression(fungal~ITS~load~(copies~g^-1)),y=expression(italic(mcrA)~(copies~g^-1)))

fig<-(pa|pb)/(pc|pd|pe|pf)+plot_layout(heights=c(1,1.15))+
  plot_annotation(tag_levels="a",tag_prefix="(",tag_suffix=")",
    theme=theme(plot.tag=element_text(size=11,face="bold"))) & theme(plot.tag=element_text(size=11,face="bold"))
ggsave("outputs/figures/generated/fig7_decay_methanogenesis.png",fig,width=9.8,height=7,dpi=300,bg="white")   # near print size
cat(sprintf("bark-flux: R2=%.2f p=%.4f dAIC=%.1f | wood mcrA-ITS: R2=%.2f p=%.4f dAIC=%.1f\n",R2a,pqa,dAICa,R2b,pqb,dAICb))
cat("wrote outputs/figures/generated/fig7_decay_methanogenesis.png\n")
