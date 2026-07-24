#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig07_decay-methanogenesis.R  (candidate expanded MAIN Fig 7 — DRAFT)
# Folds the deterioration + fungal-decay story into Fig 7:
#   TOP ROW  (population-level, a-b):
#     a  stem decay (bark loss) -> CH4 flux hump   [2023 survey]
#     b  fungal colonization (ITS load) -> mcrA hump, wood vs soil  [2021 survey]
#   BOTTOM ROW (within the felled black oak, c-f) — vertical profiles vs height:
#     c  internal CH4    d  mcrA (copies/g)    e  CH4 flux    f  fungal ITS load
#   Methanogens/gas/flux peak at 4-6 m — the leading edge of the ~4 m decay cone
#   (cavity <1 m, staining 1-2 m, clean >6 m) = intermediate/early decay, matching
#   the two population humps.
# Matches rev_fig07 felled-oak aesthetics. Fungal ITS from sibling tree-microbiome.
# Writes outputs/revision/fig7_decay_methanogenesis.png
# ==============================================================================
suppressMessages({library(tidyverse);library(patchwork);library(lme4);library(lmerTest)})
options(warn=-1)
TM<-"/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-microbiome"
asinh10<-function(x) asinh(x/0.1)/log(10)
# marginal R2 for a mixed model (fixed-effects variance / total)
r2m<-function(m){vf<-var(predict(m,re.form=NA));vc<-as.data.frame(VarCorr(m))
  vr<-sum(vc$vcov[vc$grp!="Residual"]);ve<-vc$vcov[vc$grp=="Residual"];as.numeric(vf/(vf+vr+ve))}
pfmt<-function(p) ifelse(p<0.001,"p < 0.001",sprintf("p = %.3f",p))
stat_lab<-function(R2,p,dAIC) sprintf("R2 = %.2f\n%s\ndAIC = %.1f",R2,pfmt(p),dAIC)
shared_theme <- theme_classic(base_size=12)+theme(axis.title=element_text(size=12),axis.text=element_text(size=10),
  legend.position="bottom",legend.title=element_text(size=9),legend.text=element_text(size=8),
  legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,0,0,0),plot.margin=margin(5,7,5,5))

## ===================== BOTTOM ROW: felled-oak vertical profiles =====================
O2s<-read.csv("data/raw/field_data/black_oak/o2_standards.csv",fileEncoding="UTF-8-BOM")
GHGs<-read.csv("data/raw/field_data/black_oak/ghg_standards.csv",fileEncoding="UTF-8-BOM")
GC<-read.csv("data/raw/field_data/black_oak/quve_gc_data.csv",fileEncoding="UTF-8-BOM")
colnames(GC)<-make.names(colnames(GC));colnames(O2s)<-make.names(colnames(O2s));colnames(GHGs)<-make.names(colnames(GHGs))
GC$O2.Area[GC$Sample.Type=="N2"]<-NA
O2d<-O2s%>%left_join(GC,by=c("Sample"="Sample.Type"));O2d$O2.Area<-as.numeric(O2d$O2.Area);O2d$O2_ppm<-as.numeric(gsub(",","",O2d$X.O2...ppm.))
lo<-c("Outdoor Air 2","Outdoor Air 3","Outdoor Air 4","Oxygen Standard 1","Oxygen Standard 2","Oxygen Standard 3","Oxygen Standard 4")
hi<-c("Outdoor Air 5","Outdoor Air 6","Outdoor Air 7","Oxygen Standard 5")
lod<-O2d%>%filter(Sample%in%lo);hid<-bind_rows(lod,O2d%>%filter(Sample%in%hi))
O2lo<-lm(O2_ppm~O2.Area,lod);O2hi<-lm(O2_ppm~O2.Area,hid);mxO2<-max(lod$O2.Area,na.rm=TRUE)
GHGd<-GHGs%>%left_join(GC,by=c("Sample"="Sample.Type"))
for(g in c("N2O","CH4","CO2")){GHGd[[paste0(g,".Area")]]<-as.numeric(GHGd[[paste0(g,".Area")]]);GHGd[[paste0(g,"_ppm")]]<-as.numeric(gsub(",","",GHGd[[paste0("X.",g,"...ppm.")]]))}
GHGd<-GHGd%>%filter(!is.na(CH4.Area)&!is.na(CH4_ppm))
lord<-GHGd%>%filter(Sample%in%c("N2","SB1","SB2","SB3"))
CH4lo<-lm(CH4_ppm~CH4.Area,lord);CH4hi<-lm(CH4_ppm~CH4.Area,GHGd)
GC<-GC%>%mutate(O2_concentration=if_else(O2.Area<=mxO2,predict(O2lo,data.frame(O2.Area=O2.Area)),predict(O2hi,data.frame(O2.Area=O2.Area))),
  CH4_concentration=if_else(CH4.Area<=1000,predict(CH4lo,data.frame(CH4.Area=CH4.Area)),predict(CH4hi,data.frame(CH4.Area=CH4.Area))))
int_gas<-GC%>%filter(!is.na(Lab.ID),Tree.Tissue=="Trunk Gas");int_gas$Tree.Height<-as.numeric(int_gas$Tree.Height)
hb<-sort(unique(int_gas$Tree.Height))
mcra<-readr::read_csv("data/raw/ddpcr/black_oak_mcrA.csv",show_col_types=FALSE)
cc<-grep("Conc",colnames(mcra),value=TRUE)[1];hc<-grep("Height",colnames(mcra),value=TRUE)[1]
mcra<-mcra%>%rename(Height_cm=!!hc,Conc=!!cc)
bm<-readr::read_csv("data/processed/molecular/black_oak/bo_extraction_mass.csv",show_col_types=FALSE)
mcra<-mcra%>%left_join(bm%>%transmute(`Sample ID`,mass_mg=`Sample Mass Added to Tube (mg)`),by="Sample ID")%>%
  mutate(mcrA_g=Conc*75/(mass_mg/1000))%>%filter(Component%in%c("Heartwood","Sapwood"))
flux<-read.csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv")
flux$Height_m<-suppressWarnings(as.numeric(flux$Height_m))
fdf<-flux%>%filter(!is.na(Height_m))%>%select(height=Height_m,flux=CH4_best.flux)%>%
  left_join(int_gas%>%group_by(Tree.Height)%>%summarize(ch4=mean(CH4_concentration,na.rm=TRUE),.groups="drop"),by=c("height"="Tree.Height"))
# felled-oak ITS load (heartwood), by height -> copies/g dry (same basis as mcrA panel d)
its<-read.csv(file.path(TM,"Other-Metadata/ITS_w_metadata.csv"),check.names=FALSE);its$ITS<-suppressWarnings(as.numeric(its$ITS_per_ul))
iid<-toupper(its[[1]]);itw<-its[grepl("^QUVE[0-9]+HEART",iid)&!is.na(its$ITS),]
itw$height<-as.numeric(sub("^QUVE([0-9]+).*","\\1",toupper(itw[[1]])))/100;itw$tissue<-"Heartwood"
# heartwood extraction masses (mean per height) from bo_extraction_mass -> copies/g
bm2<-bm%>%filter(Component=="Heartwood",!is.na(`Sample Mass Added to Tube (mg)`))%>%
  group_by(h_cm=`Height (cm)`)%>%summarize(mass_mg=mean(`Sample Mass Added to Tube (mg)`),.groups="drop")
itw$h_cm<-round(itw$height*100);itw<-itw%>%left_join(bm2,by="h_cm")
itw$ITS_g<-itw$ITS*75000/itw$mass_mg

col_ht<-c(Heartwood="#a6611a",Sapwood="#1f78b4")
pc<-ggplot(int_gas,aes(Tree.Height,CH4_concentration))+geom_smooth(se=FALSE,color="black")+
  geom_point(aes(fill=O2_concentration/1e6*100),size=3,shape=21,color="black",stroke=.6,alpha=.85)+shared_theme+
  scale_fill_distiller(palette="RdYlBu",direction=-1,name=expression(O[2]~"(%)"))+coord_flip()+
  ylab(expression(CH[4]~"(ppm)"))+xlab("Height (m)")+scale_x_continuous(breaks=hb,minor_breaks=NULL)+
  guides(fill=guide_colorbar(barwidth=7,barheight=.5,title.position="top"))
pd<-ggplot(mcra,aes(Height_cm/100,mcrA_g))+geom_smooth(se=FALSE,aes(color=Component))+
  geom_jitter(size=3,shape=21,aes(fill=Component),color="black",stroke=.6,alpha=.85)+shared_theme+
  scale_color_manual(values=col_ht)+scale_fill_manual(values=col_ht)+coord_flip()+
  ylab(expression(italic(mcrA)~"(copies g"^-1*")"))+xlab("")+scale_x_continuous(breaks=hb,minor_breaks=NULL)+
  guides(fill=guide_legend(title=NULL),color=guide_legend(title=NULL))
pe<-ggplot(fdf,aes(height,flux,fill=ch4))+geom_smooth(se=FALSE,color="black")+
  geom_point(size=3,shape=21,color="black",stroke=.6,alpha=.85)+shared_theme+
  scale_fill_viridis_c(option="E",name=expression(CH[4]~"(ppm)"))+coord_flip()+xlab("")+
  ylab(expression(CH[4]~"flux (nmol m"^-2*" s"^-1*")"))+scale_x_continuous(breaks=hb,minor_breaks=NULL)+
  guides(fill=guide_colorbar(barwidth=7,barheight=.5,title.position="top"))
pf<-ggplot(itw,aes(height,ITS_g+1))+geom_smooth(se=FALSE,color="black",span=1)+
  geom_jitter(width=.15,size=3,shape=21,fill="#a6611a",color="black",stroke=.6,alpha=.85)+shared_theme+
  coord_flip()+scale_y_log10(breaks=c(1,1e2+1,1e4+1,1e6+1,1e8+1),labels=c("0","1e2","1e4","1e6","1e8"))+xlab("")+
  ylab(expression("fungal ITS load (copies g"^-1*")"))+scale_x_continuous(breaks=hb,minor_breaks=NULL)

## ===================== TOP ROW: population humps =====================
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE);names(y)<-make.names(names(y))
y<-y[!is.na(y$CH4_best.flux),];y$fx<-asinh10(y$CH4_best.flux);y$sp<-as.factor(y$Species);y$x<-ord(y$Bark.Missing);y<-y[!is.na(y$x),]
mq<-lmer(fx~x+I(x^2)+(1|sp),data=y);co<-summary(mq)$coefficients;pqa<-co["I(x^2)","Pr(>|t|)"]
mla<-lmer(fx~x+(1|sp),data=y);dAICa<-as.numeric(AIC(mla)-AIC(mq));R2a<-r2m(mq)
gA<-data.frame(x=seq(1,4,.05));Xa<-model.matrix(~x+I(x^2),gA);ba<-fixef(mq);Va<-as.matrix(vcov(mq));gA$fit<-as.vector(Xa%*%ba);gA$se<-sqrt(diag(Xa%*%Va%*%t(Xa)))
obrk<-c(0,0.05,0.2,1)
# stem-flux purple, matching Fig 2 (avoids implying 'wood' via brown)
pa<-ggplot(y,aes(x,fx))+geom_hline(yintercept=0,linetype=3,colour="grey70")+
  geom_jitter(width=.12,alpha=.15,size=1.1,colour="#756BB1")+
  geom_ribbon(data=gA,aes(x,ymin=fit-1.96*se,ymax=fit+1.96*se),alpha=.2,fill="#756BB1",inherit.aes=FALSE)+
  geom_line(data=gA,aes(x,fit),linewidth=1,colour="#54278f",inherit.aes=FALSE)+
  annotate("text",1,asinh10(0.9),hjust=0,vjust=1,size=2.9,lineheight=.9,label=stat_lab(R2a,pqa,dAICa))+
  scale_x_continuous(breaks=1:4,labels=c("healthy","moderate","severe","dead"))+scale_y_continuous(breaks=asinh10(obrk),labels=obrk)+
  coord_cartesian(ylim=asinh10(c(-0.3,1)))+shared_theme+labs(x="bark loss (decay)",y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))

o<-read.csv("data/processed/molecular/tree_data_methanogen_group.csv",check.names=FALSE);o<-o[,names(o)!=""]
o$mass<-suppressWarnings(as.numeric(o$`Sample.Mass.Added.to.Tube..mg..x`))
# copies/g dry = copies/uL * 75 uL elution / mass_mg * 1000  (same basis as mcrA panel d)
o$mcra_g<-suppressWarnings(as.numeric(o$mcra_probe_loose))*75000/o$mass
o$ITS_g<-suppressWarnings(as.numeric(o$ITS_per_ul))*75000/o$mass
d<-o[is.finite(o$mcra_g)&is.finite(o$ITS_g)&o$mcra_g>0&o$ITS_g>0&o$material%in%c("Wood","Soil"),];d$sp<-as.factor(o$species.x[as.integer(rownames(d))])
d$X<-log10(d$ITS_g);d$Y<-log10(d$mcra_g)
fm<-function(mat,q){s<-d[d$material==mat,];s$xc<-s$X-mean(s$X);mm<-if(q)lmer(Y~xc+I(xc^2)+(1|sp),s) else lmer(Y~X+(1|sp),s)
  g<-data.frame(X=seq(min(s$X),max(s$X),length=60));g$xc<-g$X-mean(s$X);Xm<-if(q)model.matrix(~xc+I(xc^2),g) else model.matrix(~X,g)
  b<-fixef(mm);V<-as.matrix(vcov(mm));g$fit<-as.vector(Xm%*%b);g$se<-sqrt(diag(Xm%*%V%*%t(Xm)));g$material<-mat;g}
gw<-fm("Wood",TRUE);gs<-fm("Soil",FALSE)
# wood stats: quadratic (hump) vs linear
dw<-d[d$material=="Wood",];dw$xc<-dw$X-mean(dw$X)
mqb<-lmer(Y~xc+I(xc^2)+(1|sp),dw);mlb<-lmer(Y~xc+(1|sp),dw)
pqb<-summary(mqb)$coefficients[3,"Pr(>|t|)"];dAICb<-as.numeric(AIC(mlb)-AIC(mqb));R2b<-r2m(mqb)
# paper BrBG convention: wood = brown, soil = teal (matches Fig 4/5 compartment palette)
col_ws<-c(Wood="#a6611a",Soil="#018571")
pb<-ggplot(d,aes(X,Y,colour=material))+geom_point(alpha=.4,size=1.2)+
  geom_ribbon(data=rbind(gw,gs),aes(X,ymin=fit-1.96*se,ymax=fit+1.96*se,fill=material),alpha=.15,colour=NA,inherit.aes=FALSE)+
  geom_line(data=gw,aes(X,fit),colour=col_ws["Wood"],linewidth=1,inherit.aes=FALSE)+
  geom_line(data=gs,aes(X,fit),colour=col_ws["Soil"],linewidth=1,inherit.aes=FALSE)+
  scale_colour_manual(values=col_ws,name=NULL)+scale_fill_manual(values=col_ws,guide="none")+
  annotate("text",-Inf,Inf,hjust=-.05,vjust=1.2,size=2.9,lineheight=.9,label=paste0("wood:\n",stat_lab(R2b,pqb,dAICb)))+
  shared_theme+theme(legend.position="right")+
  labs(x=expression(log[10]~fungal~ITS~load~(copies~g^-1)),y=expression(log[10]~italic(mcrA)~(copies~g^-1)))

fig<-(pa|pb)/(pc|pd|pe|pf)+plot_layout(heights=c(1,1.15))+
  plot_annotation(tag_levels="a",tag_prefix="(",tag_suffix=")",
    title="Methane production peaks at moderate wood decay",
    theme=theme(plot.tag=element_text(size=11,face="bold"),plot.title=element_text(face="bold",size=14)))
ggsave("outputs/revision/fig7_decay_methanogenesis.png",fig,width=13,height=9,dpi=250)
cat(sprintf("bark-flux hump p=%.4f | wood mcrA-ITS hump p=%.4f\n",pqa,pqb))
cat("wrote outputs/revision/fig7_decay_methanogenesis.png\n")
