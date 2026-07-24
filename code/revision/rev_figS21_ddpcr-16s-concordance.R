#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS21_ddpcr-16s-concordance.R  (SI)
# ddPCR functional-gene detection vs 16S TAXONOMY-based abundance of the same
# groups (methanogen families / methanotroph ASVs; SAME curated definitions as
# Fig 5, not FAPROTAX). Split by WOOD vs SOIL. ddPCR is far more sensitive, so:
# many samples ddPCR+/16S=0, almost none ddPCR=0/16S+. The low 16S methanogen
# detection is compartment-structured (16S sees methanogens in heartwood, misses
# them in soil/sapwood) and matches the manuscript ("methanogens in heart/sapwood,
# not soil") -> this VALIDATES using ddPCR for quantification. Non-detects shown
# at a labelled "0"; black line = quantitative fit (both detect) with p. Full read
# depth (fairer for a detection comparison; Fig 5 rarefies for composition).
# Self-contained (our data/raw/16s/). Writes outputs/revision/figS21_...png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork)}); options(warn=-1,stringsAsFactors=FALSE)
num<-function(x) suppressWarnings(as.numeric(x))
otu<-read.delim("data/raw/16s/OTU_table.txt",header=TRUE,row.names=1,check.names=FALSE)
taxc<-intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),names(otu))
tax<-otu[,taxc]; cnt<-as.matrix(sapply(otu[,setdiff(names(otu),taxc),drop=FALSE],num)); rownames(cnt)<-rownames(otu); tot<-colSums(cnt,na.rm=TRUE)
mgfam<-c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae","Methanocellaceae",
         "Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae","Methanomethyliaceae","Methanocorpusculaceae")
mg_asv<-rownames(otu)[tax$Family %in% mgfam]
source("code/00_harmonization/load_methanotroph_definitions.R")
dp<-if(file.exists("data/processed/molecular/methanotroph_definitions_revised.csv"))"data/processed/molecular/methanotroph_definitions_revised.csv" else "data/compiled/methanotroph_definitions.csv"
tax$mt<-classify_methanotrophs(tax,load_methanotroph_defs(dp),include_conditional=TRUE)
mt_asv<-rownames(otu)[tax$mt %in% c("Known","Putative")]
s16<-data.frame(key=sub("[.]16S[.]S[0-9]*$","",colnames(cnt)),
                mg16=100*colSums(cnt[mg_asv,,drop=FALSE],na.rm=TRUE)/tot,
                mt16=100*colSums(cnt[mt_asv,,drop=FALSE],na.rm=TRUE)/tot)
o<-read.csv("data/processed/molecular/tree_data_methanogen_group.csv",check.names=FALSE); o<-o[,names(o)!=""]
o$mcra<-num(o$mcra_probe_loose); o$pmoa<-num(o$pmoa_loose); o$mmox<-num(o$mmox_loose); o$key<-paste0(o$seq_id,o$core_type)
m<-merge(o[,c("key","material","mcra","pmoa","mmox")],s16,by="key")
cols<-c("both detect"="#2c7fb8","ddPCR+ / 16S = 0"="#d95f02","ddPCR = 0 / 16S+"="#7570b3","both = 0"="grey65")
set.seed(1)
panel<-function(gene,func,mat,glab,flab){
  s<-m[m$material==mat & is.finite(m[[gene]])&is.finite(m[[func]]),]; s$g<-s[[gene]]; s$f<-s[[func]]
  s$cat<-factor(with(s,ifelse(g>0&f>0,"both detect",ifelse(g>0&f==0,"ddPCR+ / 16S = 0",
             ifelse(g==0&f>0,"ddPCR = 0 / 16S+","both = 0")))),levels=names(cols))
  fx<-function(v){p<-v[v>0]; fl<-if(length(p))10^(floor(log10(min(p)))-1) else 1e-3; list(v=ifelse(v==0, fl*10^runif(length(v),-.2,.2), v),fl=fl,top=if(length(p))max(p) else 1)}
  gg<-fx(s$g); ff<-fx(s$f); s$gp<-gg$v; s$fp<-ff$v
  br<-function(fl,top){b<-10^(seq(floor(log10(fl)),ceiling(log10(top)))); list(b=b,l=ifelse(b==fl,"0",parse(text=paste0("10^",round(log10(b))))))}
  bg<-br(gg$fl,gg$top); bf<-br(ff$fl,ff$top); pos<-s[s$g>0&s$f>0,]; fit<-NULL; rp<-""
  if(nrow(pos)>=4){ct<-suppressWarnings(cor.test(log10(pos$g),log10(pos$f)))
    cf<-coef(lm(log10(g)~log10(f),data=pos)); fit<-data.frame(f=10^seq(log10(min(pos$f)),log10(max(pos$f)),length=50)); fit$g<-10^(cf[1]+cf[2]*log10(fit$f))
    rp<-sprintf("r=%.2f, p%s",ct$estimate,ifelse(ct$p.value<0.001,"<0.001",sprintf("=%.3f",ct$p.value)))}
  tab<-table(s$cat)
  lab<-sprintf("both: %d  %s\nddPCR+/16S=0: %d ; 16S+/ddPCR=0: %d\ndetect ddPCR %.0f%% / 16S %.0f%%",
     tab["both detect"], rp, tab["ddPCR+ / 16S = 0"], tab["ddPCR = 0 / 16S+"], 100*mean(s$g>0),100*mean(s$f>0))
  p<-ggplot(s,aes(fp,gp,colour=cat))+
    annotate("rect",xmin=ff$fl/3,xmax=ff$fl*3,ymin=-Inf,ymax=Inf,fill="grey92",alpha=.6)+
    annotate("rect",ymin=gg$fl/3,ymax=gg$fl*3,xmin=-Inf,xmax=Inf,fill="grey92",alpha=.6)+
    geom_point(alpha=.6,size=1.3)+scale_colour_manual(values=cols,name=NULL,drop=FALSE)+
    scale_x_log10(breaks=bf$b,labels=bf$l)+scale_y_log10(breaks=bg$b,labels=bg$l)+
    annotate("text",ff$fl,gg$top,label=lab,hjust=0,vjust=1,size=2.5)+
    labs(title=glab,x=flab,y=bquote(.(gene)~"ddPCR (copies, 0=nd)"))+
    theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position="bottom")+
    guides(colour=guide_legend(override.aes=list(size=3,alpha=1),nrow=1))
  if(!is.null(fit)) p<-p+geom_line(data=fit,aes(f,g),colour="black",linewidth=.5,inherit.aes=FALSE)
  p
}
g<-list(c("mcra","mg16","mcrA","methanogen"),c("pmoa","mt16","pmoA","methanotroph"),c("mmox","mt16","mmoX","methanotroph"))
mk<-function(mat) lapply(g,function(x) panel(x[1],x[2],mat,sprintf("%s vs 16S %s (%s)",x[3],x[4],mat),sprintf("16S %s taxa (%%, 0=nd)",x[4])))
W<-mk("Wood"); S<-mk("Soil")
fig<-(W[[1]]|W[[2]]|W[[3]])/(S[[1]]|S[[2]]|S[[3]])+plot_layout(guides="collect")+
  plot_annotation(title="ddPCR detects methane-cycling genes where 16S taxa are absent (wood vs soil)",
    subtitle="16S methanogens are seen mainly in heartwood, ~none in soil (matches manuscript); ddPCR detects them throughout. taxonomy-based (Fig 5 defs).",
    theme=theme(plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=9),legend.position="bottom"))
ggsave("outputs/revision/figS21_ddpcr_16s_concordance.png",fig,width=13.5,height=9,dpi=200)
cat("wrote outputs/revision/figS21_ddpcr_16s_concordance.png\n")
