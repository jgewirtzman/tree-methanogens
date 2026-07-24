suppressMessages({library(microeco);library(tidyverse)});options(warn=-1)
BO<-"/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-microbiome/Black Oak/ITS"
otu<-read.delim(file.path(BO,"OTU_table.txt"),row.names=1,check.names=FALSE)
otu[]<-lapply(otu,function(c)suppressWarnings(as.numeric(as.character(c))))
otu<-otu[,colSums(!is.na(otu))>0,drop=FALSE];otu[is.na(otu)]<-0  # drop non-count (e.g. taxonomy) cols
txr<-read.delim(file.path(BO,"taxonomy_table.txt"),check.names=FALSE);rownames(txr)<-txr[[1]]
ranks<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
parts<-strsplit(as.character(txr$Taxon),";\\s*")
taxm<-t(sapply(parts,function(p){p<-trimws(p);length(p)<-7;p}))
taxm[is.na(taxm)|taxm%in%c("none","unidentified","","k__","p__","c__","o__","f__","g__","s__")]<-""
tax_df<-as.data.frame(taxm,stringsAsFactors=FALSE);colnames(tax_df)<-ranks;rownames(tax_df)<-rownames(txr)
common<-intersect(rownames(otu),rownames(tax_df))
otu<-otu[common,,drop=FALSE];tax_df<-tax_df[common,,drop=FALSE]
samp<-data.frame(SampleID=colnames(otu));rownames(samp)<-colnames(otu)
tax_df<-tidy_taxonomy(tax_df)
mt<-suppressMessages(microtable$new(otu_table=otu,tax_table=tax_df,sample_table=samp))
tf<-trans_func$new(mt)
invisible(capture.output(tf$cal_spe_func(fungi_database="FUNGuild")))
invisible(capture.output(tf$cal_spe_func_perc(abundance_weighted=TRUE)))
perc<-as.data.frame(tf$res_spe_func_perc)
cat("guild columns available:\n");print(grep("aprotroph|Patho|Symbio|Myco",colnames(perc),value=TRUE))
perc$sample<-rownames(perc)
q<-perc[grepl("^QUVE[0-9]+(HEART|SAP)",toupper(gsub("[.].*","",perc$sample))),]
qid<-toupper(gsub("[.].*","",q$sample))
q$height<-as.numeric(sub("^QUVE([0-9]+).*","\\1",qid))/100
q$tissue<-ifelse(grepl("SAP",qid),"Sapwood","Heartwood")
out<-q%>%transmute(sample,height,tissue,
  Saprotroph=Saprotroph,`Wood Saprotroph`=`Wood Saprotroph`,Pathotroph=Pathotroph,Symbiotroph=Symbiotroph)%>%arrange(tissue,height)
cat("\nper-QUVE-sample guild % (abundance-weighted):\n");print(as.data.frame(out),digits=3)
write.csv(out,"outputs/revision/bo_quve_funguild_perc.csv",row.names=FALSE)
cat("\nwrote outputs/revision/bo_quve_funguild_perc.csv\n")
