#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_decay-methanogenesis.R   (candidate MAIN figure — DRAFT)
# "Methane production peaks at moderate wood decay." Three convergent lines:
#   A  Stem decay (bark loss) -> CH4 flux hump          [2023 flux survey]
#   B  Fungal colonization (ITS load) -> methanogen (mcrA) hump, WOOD only; soil flat
#                                                        [2021 molecular survey]
#   C  Wood fungal community is decay/saprotroph-dominated (~0 mycorrhizal) vs soil
#      mycorrhizal-dominated -> specificity contrast     [ITS FUNGuild composition]
# Fungal ITS data come from the sibling repo (tree-microbiome); path documented below.
# Reproducible; writes outputs/revision/fig_decay_methanogenesis.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest)})
options(warn=-1,stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
TM<-"/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-microbiome"  # sibling ITS repo
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
predci<-function(m,g,q){X<-if(q)model.matrix(~x+I(x^2),g) else model.matrix(~x,g);b<-fixef(m);V<-as.matrix(vcov(m))
  g$fit<-as.vector(X%*%b);g$se<-sqrt(diag(X%*%V%*%t(X)));g$lo<-g$fit-1.96*g$se;g$hi<-g$fit+1.96*g$se;g}

## ---- A: bark decay -> flux (2023 survey) ----
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species); y$x<-ord(y$Bark.Missing); y<-y[!is.na(y$x),]
mqA<-lmer(fx~x+I(x^2)+(1|sp),data=y); pqA<-summary(mqA)$coefficients["I(x^2)","Pr(>|t|)"]
gA<-predci(mqA,data.frame(x=seq(1,4,.05)),TRUE); obrk<-c(-0.2,0,0.05,0.2,1)
pA<-ggplot(y,aes(x,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_jitter(width=.12,alpha=.14,size=1,colour="#b8732b")+
  geom_ribbon(data=gA,aes(x,ymin=lo,ymax=hi),alpha=.18,fill="#8c510a",inherit.aes=FALSE)+
  geom_line(data=gA,aes(x,fit),linewidth=1,colour="#5e3a06",inherit.aes=FALSE)+
  annotate("text",1,asinh10(0.9),hjust=0,size=3,label=paste("hump",ifelse(pqA<0.001,"p<0.001",sprintf("p=%.3f",pqA))))+
  scale_x_continuous(breaks=1:4,labels=c("healthy","mod","sev","dead"))+
  scale_y_continuous(breaks=asinh10(obrk),labels=obrk)+coord_cartesian(ylim=asinh10(c(-0.3,1)))+
  labs(title="A  Stem decay -> CH4 flux",subtitle="2023 flux survey (n=282)",
       x="bark loss (decay)",y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))+th

## ---- molecular join (2021): mcrA, ITS load, guilds ----
o<-read.csv("data/processed/molecular/tree_data_methanogen_group.csv",check.names=FALSE); o<-o[,names(o)!=""]
o$mcra<-suppressWarnings(as.numeric(o$mcra_probe_loose)); o$ITS<-suppressWarnings(as.numeric(o$ITS_per_ul))
o$key<-paste0(o$seq_id,o$core_type)
g<-read.csv(file.path(TM,"metabolisms/ITS_metabolisms_weighted.csv"),check.names=FALSE)
names(g)[1]<-"gid"; g$key<-sub("[.]ITS[.]S[0-9]+$","",g$gid)
gk<-c("Saprotroph","Symbiotroph","Pathotroph")   # mutually-exclusive FUNGuild trophic modes
m<-merge(o[,c("key","core_type","material","species.x","mcra","ITS")],g[,c("key",gk)],by="key")
names(m)<-make.names(names(m)); m$sp<-as.factor(m$species.x)

## ---- B: ITS load -> mcrA, wood hump vs soil flat ----
d<-m[is.finite(m$mcra)&is.finite(m$ITS)&m$mcra>0&m$ITS>0 & m$material%in%c("Wood","Soil"),]
d$x<-log10(d$ITS); d$y<-log10(d$mcra)
fitmat<-function(mat,q){s<-d[d$material==mat,]; s$xc<-s$x-mean(s$x)
  mm<-if(q) lmer(y~xc+I(xc^2)+(1|sp),data=s) else lmer(y~x+(1|sp),data=s)
  g<-data.frame(x=seq(min(s$x),max(s$x),length=80)); g$xc<-g$x-mean(s$x)
  if(q){X<-model.matrix(~xc+I(xc^2),g)}else{X<-model.matrix(~x,g)}
  b<-fixef(mm);V<-as.matrix(vcov(mm)); g$fit<-as.vector(X%*%b);g$se<-sqrt(diag(X%*%V%*%t(X)))
  g$lo<-g$fit-1.96*g$se;g$hi<-g$fit+1.96*g$se;g$material<-mat;g}
gw<-fitmat("Wood",TRUE); gs<-fitmat("Soil",FALSE)
pqB<-summary(lmer(y~I(x-mean(x))+I((x-mean(x))^2)+(1|sp),data=d[d$material=="Wood",]))$coefficients[3,"Pr(>|t|)"]
cols<-c(Wood="#2c7a2c",Soil="#8c6d3f")
pB<-ggplot(d,aes(x,y,colour=material))+
  geom_point(alpha=.4,size=1.3)+
  geom_ribbon(data=rbind(gw,gs),aes(x,ymin=lo,ymax=hi,fill=material),alpha=.15,colour=NA,inherit.aes=FALSE)+
  geom_line(data=gw,aes(x,fit),colour=cols["Wood"],linewidth=1,inherit.aes=FALSE)+
  geom_line(data=gs,aes(x,fit),colour=cols["Soil"],linewidth=1,linetype=1,inherit.aes=FALSE)+
  scale_colour_manual(values=cols,name=NULL)+scale_fill_manual(values=cols,guide="none")+
  annotate("text",-Inf,Inf,hjust=-.05,vjust=1.3,size=3,label=sprintf("wood hump p=%.4f\nsoil: flat",pqB))+
  labs(title="B  Fungal load -> methanogens",subtitle="2021 molecular survey",
       x=expression(log[10]~fungal~ITS~load),y=expression(log[10]~italic(mcrA)~load))+
  th+theme(legend.position=c(.99,.01),legend.justification=c(1,0))

## ---- C: guild composition by compartment ----
library(tidyr)
cc<-aggregate(cbind(Saprotroph,Symbiotroph,Pathotroph)~core_type,data=m,mean)
cc$Other<-pmax(0,100-rowSums(cc[,-1])); cl<-pivot_longer(cc,-core_type,names_to="guild",values_to="pct")
cl$core_type<-factor(cl$core_type,levels=c("Inner","Outer","Organic","Mineral"),
                     labels=c("heartwood","sapwood","organic","mineral"))
cl$guild<-factor(cl$guild,levels=c("Saprotroph","Pathotroph","Symbiotroph","Other"))
gcol<-c(Saprotroph="#8c510a",Pathotroph="#d8a35b",Symbiotroph="#2c7a2c",Other="grey80")
pC<-ggplot(cl,aes(core_type,pct,fill=guild))+geom_col(width=.7)+
  scale_fill_manual(values=gcol,name="fungal guild")+
  annotate("segment",x=.5,xend=2.5,y=104,yend=104,colour="#2c7a2c",linewidth=1)+
  annotate("segment",x=2.5,xend=4.5,y=104,yend=104,colour="#8c6d3f",linewidth=1)+
  annotate("text",1.5,108,label="wood",size=3,colour="#2c7a2c")+annotate("text",3.5,108,label="soil",size=3,colour="#8c6d3f")+
  labs(title="C  Fungal community by compartment",subtitle="wood = saprotroph (decay); soil = symbiotroph (mycorrhizal)",
       x=NULL,y="% of community")+coord_cartesian(ylim=c(0,112))+th+theme(axis.text.x=element_text(angle=20,hjust=1))

fig<-(pA|pB|pC)+plot_layout(widths=c(1,1.1,1))+
  plot_annotation(title="Methane production peaks at moderate wood decay",
    theme=theme(plot.title=element_text(face="bold",size=14)))
ggsave("outputs/revision/fig_decay_methanogenesis.png",fig,width=15,height=5,dpi=200)
cat(sprintf("A bark-flux hump p=%.4f | B wood mcrA-ITS hump p=%.4f\n",pqA,pqB))
cat("wrote outputs/revision/fig_decay_methanogenesis.png\n")
