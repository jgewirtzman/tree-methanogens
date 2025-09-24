library(plyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
theme_set(theme_bw(base_size = 12))

#### Import guide csvs ####
inner_key=read.csv("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/inner_guide.csv", skip=1, header = TRUE, stringsAsFactors = TRUE)
outer_key=read.csv("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/outer_guide.csv", skip=1, header = TRUE, stringsAsFactors = TRUE)
mineral_key=read.csv("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/mineral_guide.csv", skip=1, header = TRUE, stringsAsFactors = TRUE)
organic_key=read.csv("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/organic_guide.csv", skip=1, header = TRUE, stringsAsFactors = TRUE)


inner_key=rbind(inner_key, outer_key, mineral_key, organic_key)

#Extract 4 char species ID
inner_key$species=substr(inner_key$Inner.Core.Sample.ID, 1, 4)

#Extract core inner/outer identifier, make unique identifier for samples (Plate#WellLocation, e.g. 6A1)
inner_key$core_type=sub("\nL.*","",inner_key$Inner.Core.Sample.ID)
inner_key$core_type=sub(".*\n","",inner_key$core_type)
inner_key$material=sub("_.*","",inner_key$core_type)
inner_key$core_type=sub(".*_","",inner_key$core_type)
inner_key$lab_id=sub(".*: ","",inner_key$Inner.Core.Sample.ID)
inner_key$plate_identifier=paste(sub("Plate ","",inner_key$Extraction.Plate.ID), inner_key$Plate.Position, sep="")

#### Import ddPCR files ####

data_files = paste(("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/data/"),list.files("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/data"), sep="")

#Check DNA or RNA, and strict or loose

ddpcr_data=data.frame()
for (i in 1:length(data_files)){
  print(paste("Imported file #:",i))
  if(str_detect(data_files[i], "loose")){        # check file name for loose/strict
    data_file=read.csv(data_files[i])            # read in file
    data_file=data_file[,1:23]                   # cutoff junk
    data_file$Sample.description.4="loose"       # add into column
    ddpcr_data=rbind(ddpcr_data, data_file)      # append to datafile
  }
  else if(str_detect(data_files[i], "strict")){
    data_file=read.csv(data_files[i])
    data_file=data_file[,1:23]
    data_file$Sample.description.4="strict"
    ddpcr_data=rbind(ddpcr_data, data_file)
  }
  else{
    data_file=read.csv(data_files[i])
    data_file$Sample.description.4="none"
    data_file=data_file[,1:23]
    ddpcr_data=rbind(ddpcr_data, data_file)
  }
}

# Unique plate identifier based on plate # and well, need to format first
ddpcr_data$ï..Well=paste(substr(ddpcr_data$ï..Well, 0,1),as.numeric(substr(ddpcr_data$ï..Well, 2,4)), sep="")
ddpcr_data$Sample.description.1=sub("_"," ",ddpcr_data$Sample.description.1)
ddpcr_data$plate_identifier=paste(sub("Plate ","",ddpcr_data$Sample.description.1),ddpcr_data$ï..Well, sep="")

# Merge metadta and ddPCR data
merged_all=merge(inner_key, ddpcr_data, by = "plate_identifier")

# Remove stds and blank wells
merge_no_blanks=subset(merged_all, grepl(merged_all$Sample.description.2,pattern = "std", ignore.case=TRUE)==FALSE)
merge_no_blanks=subset(merge_no_blanks, grepl(merge_no_blanks$Blank.Well., pattern = "yes", ignore.case = TRUE)==FALSE)

#Redo
redo_samples=subset(merge_no_blanks, Accepted.Droplets<5000)

#Standardize target formatting
merge_no_blanks$Target=tolower(merge_no_blanks$Target)

#### Quick Plots ####
mcra_only=subset(merge_no_blanks, Target=="mmox" & Sample.description.4=="strict")

species_bar<-ggplot(data=mcra_only, aes(x=species, y=log10(1+mcra_only$Copies.20ÂµLWell), col=mcra_only$Sample.description.4)) +
  geom_point(size=0.9, shape=21)+
 # geom_boxplot()+
  stat_summary(size=1.1)+
  ggtitle(paste("mcrA DNA (FAM Probe)"))+ 
  facet_grid(.~mcra_only$core_type, scales = "free", space="free")+
#  ylim(0,4.75)+
  ylab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  xlab("Species ID")+ 
  scale_color_manual(name="Cutoff Threshold \n (Mean SE)", values=wes_palette("Darjeeling1"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(species_bar)


#### Transpose variables ####

# Need to have data cleaned before transpose; metadata must be lost
ddpcr_data_clean=subset(ddpcr_data, grepl(ddpcr_data$Sample.description.2,pattern = "std", ignore.case = TRUE)==FALSE)
ddpcr_data_clean$Target=tolower(ddpcr_data_clean$Target)
ddpcr_data_clean=subset(ddpcr_data_clean, ddpcr_data_clean$Accepted.Droplets>7500)
inner_key_clean=subset(inner_key, grepl(inner_key$Blank.Well., pattern = "yes", ignore.case = TRUE)==FALSE)
ddpcr_trans=inner_key_clean
ddpcr_trans$seq_id=sub(".*_","",sub("\n.*","",ddpcr_trans$Inner.Core.Sample.ID))


# Assign 
var_trans=function(target, threshold, name){
  mcra_probe_loose=subset(ddpcr_data_clean, Target==target & Sample.description.4==threshold)   # Speciy gene and cutoff thresh
  mcra_probe_loose=mcra_probe_loose[,c(16,24)]                                                  # Extraxt just copies and key
  colnames(mcra_probe_loose)=c(name,"plate_identifier")                                         # Name copies column gene name
  return(merge(ddpcr_trans, mcra_probe_loose, by="plate_identifier", all.x=TRUE))               # Merge by key
}

ddpcr_trans=var_trans("mcra_probe", "loose", "mcra_probe_loose")
ddpcr_trans=var_trans("mcra_probe", "strict", "mcra_probe_strict")
ddpcr_trans=var_trans("mcra", "loose", "mcra_eva_loose")
ddpcr_trans=var_trans("mcra", "strict", "mcra_eva_strict")
ddpcr_trans=var_trans("16s_bact", "loose", "bact_loose")
ddpcr_trans=var_trans("16s_bact", "strict", "bact_strict")
ddpcr_trans=var_trans("16s_arc", "loose", "arc_loose")
ddpcr_trans=var_trans("16s_arc", "strict", "arc_strict")
ddpcr_trans=var_trans("mmox", "loose", "mmox_loose")
ddpcr_trans=var_trans("mmox", "strict", "mmox_strict")
ddpcr_trans=var_trans("pmoa", "loose", "pmoa_loose")
ddpcr_trans=var_trans("pmoa", "strict", "pmoa_strict")

#### Export CSV ####


ddpcr_trans$Raw=paste(sub("H","HM",ddpcr_trans$seq_id), ddpcr_trans$core_type, sep="") # hemlock is HM in seq but H in key
#write.csv(ddpcr_trans,"ddPCR_tree_transposed_data.csv")

#### Import Soil Data ####

soil_data=read.csv("/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/Soils_data.csv")

for(i in 1:dim(soil_data)){  # Fill in missing labels
  if(soil_data$Tree_ID[i]==""){
    soil_data$Tree_ID[i]=soil_data$Tree_ID[i-1]
  }
}

for(i in 1:dim(soil_data)){  # Fill in missing labels
  if(soil_data$Tree_sp[i]==""){
    soil_data$Tree_sp[i]=soil_data$Tree_sp[i-1]
  }
}

library(dplyr)

soil_collapse=soil_data %>%  # Collapse cardinal directions into mean
  group_by(Tree_ID) %>%
  summarise(mean_vwc=mean(c(as.numeric(VWC1....),as.numeric(VWC2....),as.numeric(VWC3....)), na.rm=TRUE),
            mean_orp=mean(as.numeric(ORP..mV.), na.rm=TRUE),
            mean_st=mean(as.numeric(Soil.Temp), na.rm=TRUE))

colnames(soil_collapse)[1]="seq_id"
write.csv(soil_collapse,"sd_needs_clean.csv") # write for manual corrections
soil_col_clean=read.csv("C:/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/sd_clean.csv")


ddpcr_meta=merge(ddpcr_trans,soil_col_clean, by="seq_id", all.x=TRUE) # Lost some samples?



#### Soil graphs ####

vwc<-ggplot(data=ddpcr_meta,aes(x=log10(ddpcr_meta$mcra_probe_loose+mmox_strict+pmoa_strict+1), 
                                y=mean_vwc)) +
  geom_point(data=ddpcr_meta, aes(col=species), size=3)+
  ggtitle(paste("DNA Abundances vs. VWC"))+ 
  stat_cor()+
  geom_smooth(method="lm")+
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Mean~VWC))+ 
  xlab(expression(Log[10](Sum~All~Targets~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  scale_color_viridis_c(name="VWC %")
print(vwc)

st<-ggplot(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+mmox_strict+pmoa_strict+1), 
                              y=mean_st)) +
  geom_point(data=ddpcr_meta, aes(col=species), size=3)+
  ggtitle(paste("DNA Abundances vs. Soil Temperature"))+ 
  stat_cor()+
  geom_smooth(method="lm")+
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  scale_color_continuous(name="°C")
print(st)

orp<-ggplot(data=ddpcr_meta) +
  geom_point(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+1), 
                                  y=log10(pmoa_loose+1), col=mean_orp), size=3)+
  ggtitle(paste("DNA Abundances vs. ORP"))+ 
  
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_gradient(name="ORP", low="green", high="red")
print(orp)

library(RColorBrewer)
coul <- c(brewer.pal(9, "Set3"),brewer.pal(12, "Paired"),rev(brewer.pal(8, "Accent")),brewer.pal(8, "Dark2"),
          brewer.pal(12, "Paired"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))



species1<-ggplot(data=ddpcr_meta, aes(x=log10((ddpcr_meta$mcra_probe_loose*37.5)/ddpcr_meta$Sample.Mass.Added.to.Tube..mg.+1), 
                                      y=log10(((pmoa_strict+mmox_strict)*37.5)/ddpcr_meta$Sample.Mass.Added.to.Tube..mg.+1), col=core_type)) +
  geom_point(data=ddpcr_meta, size=4, alpha=0.9)+
  ggtitle(paste("DNA Abundances vs. Species"))+ 
  
 # facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  xlim(0,4.9)+
  ylim(0,3.7)+
  ylab(expression(Log[10](pmoA+mmoX~(Copies~mg^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~mg^-1))))+ 
  scale_shape_discrete(name="Mean")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Core Type", values=coul[c(4:15,17:25)])+
  stat_mean(data=ddpcr_meta, aes(shape=core_type), size=4, col="black")
#  theme(legend.position = "none")
print(species1)

# Add histograms to the margins
species1_with_hist <- ggExtra::ggMarginal(species1, type="histogram", margins="both", groupFill = TRUE, groupColour = TRUE)
print(species1_with_hist)


species<-ggplot(data=ddpcr_meta) +
  geom_point(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+1), 
                                  y=log10(pmoa_loose+1), col=species), size=3)+
  ggtitle(paste("DNA Abundances vs. Species"))+ 
  
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[c(1,3:15,17:20)])+
    stat_ellipse(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+1), 
                                     y=log10(pmoa_loose+1), col=species), size=2, alpha=0.75)
print(species)

cowplot::plot_grid(vwc, st, orp, nrow=1)
cowplot::plot_grid(species1, species, nrow=1)



#### Flux data import ####

flux_data_raw=read.csv("C:/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/flux_data_TEMP_mod.csv")

flux_50 = subset(flux_data_raw, Measurement.height==50)
flux_50=as.data.frame(cbind(flux_50$seq_id, flux_50$CH4_lin_flux, flux_50$CO2_lin_flux))
colnames(flux_50)=c("seq_id","ch4_50","co2_50")

flux_125= subset(flux_data_raw, Measurement.height==125)
flux_125=as.data.frame(cbind(flux_125$seq_id, flux_125$CH4_lin_flux, flux_125$CO2_lin_flux))
colnames(flux_125)=c("seq_id","ch4_125","co2_125")

flux_200= subset(flux_data_raw, Measurement.height==200)
flux_200=as.data.frame(cbind(flux_200$seq_id, flux_200$CH4_lin_flux, flux_200$CO2_lin_flux))
colnames(flux_200)=c("seq_id","ch4_200","co2_200")

fluxes_split=merge(flux_50, flux_125, by="seq_id", all.x = TRUE )
fluxes_split=merge(fluxes_split, flux_200, by="seq_id", all.x = TRUE)

ddpcr_meta_flux=merge(ddpcr_meta, fluxes_split, by="seq_id", all.x=TRUE)


# Export csv again
ddpcr_meta_flux$Raw=paste(sub("H","HM",ddpcr_meta_flux$seq_id), ddpcr_meta_flux$core_type, sep="") # hemlock is HM in seq but H in key
#write.csv(ddpcr_all,"ddPCR_meta_transposed_data.csv")


#### Flux graphs ####

ddpcr_meta_flux=subset(ddpcr_meta_flux, species!="KALA"&
species!="SAAL" &
  species!="QUAL" &
  species!="QUVE" &
  species!="CAOV" &
  species!="PRSE")

f50<-ggplot(data=ddpcr_meta_flux,aes(x=log10(ddpcr_meta_flux$mcra_probe_loose+1), 
                                                             y=log10(as.numeric(ddpcr_meta_flux$ch4_50)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("DNA Abundances vs. Flux at 50cm"))+ 
  facet_wrap(ddpcr_meta_flux$core_type~species, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(f50)

f125<-ggplot(data=ddpcr_meta_flux,aes(x=log10(ddpcr_meta_flux$mcra_probe_loose+1), 
                                     y=log10(as.numeric(ddpcr_meta_flux$ch4_125)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("DNA Abundances vs. Flux at 125cm"))+ 
  facet_grid(ddpcr_meta_flux$core_type~.)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(f125)

f200<-ggplot(data=ddpcr_meta_flux,aes(x=log10(ddpcr_meta_flux$mcra_probe_loose+1), 
                                      y=log10(as.numeric(ddpcr_meta_flux$ch4_200)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("DNA Abundances vs. Flux at 200cm"))+ 
  facet_grid(ddpcr_meta_flux$core_type~.)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(f200)

cowplot::plot_grid(f50, f125, f200, nrow=1)




ddpcr_meta_flux_justsoil=subset(ddpcr_all, material=="Soil")
ddpcr_meta_flux_justsoil=subset(ddpcr_meta_flux_justsoil, species!="KALA"&
                         species!="SAAL" &
                         species!="QUAL" &
                         species!="QUVE" &
                         species!="CAOV" &
                         species!="PRSE")




int_soil<-ggplot(data=ddpcr_all,aes(x=log10(ddpcr_all$mcra_probe_loose+1), 
                                      y=log10(as.numeric(ddpcr_all$CO2_int)))) +
  geom_point(size=4, aes(col=species))+
  ggtitle(paste("Soil mcrA Abundances vs. Tree Internal CO2"))+ 
  facet_wrap(ddpcr_all$core_type~species, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Conc)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(int_soil)


soil_trophs<-ggplot(data=ddpcr_meta_flux_justsoil,aes(x=log10(ddpcr_meta_flux_justsoil$mmox_loose+1), 
                                                  y=log10(as.numeric(ddpcr_meta_flux_justsoil$mcra_probe_loose+1)))) +
  geom_point(size=4, aes(col=mean_vwc))+
  ggtitle(paste("Soil mmoX Abundances vs. Soil pmoA Abundances"))+ 
 # facet_wrap(ddpcr_meta_flux_justsoil$core_type~., ncol=10)+
  geom_vline(xintercept=mean(log10(ddpcr_meta_flux_justsoil$mmox_loose+1), na.rm=TRUE))+
  geom_hline(yintercept=mean(log10(ddpcr_meta_flux_justsoil$pmoa_loose+1), na.rm=TRUE))+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mmoX~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_viridis_c(name="Mean VWC")
 # scale_color_manual(name="Mean VWC", values=coul[10:30])
print(soil_trophs)


mmox_bar<-ggplot(data=ddpcr_meta_flux, aes(y=log10(as.numeric(ddpcr_meta_flux$mcra_probe_loose*37.5)/ddpcr_meta_flux$Sample.Mass.Added.to.Tube..mg.+1), x=species)) +
  geom_jitter(shape=21, aes(fill=ddpcr_meta_flux$species), size=3)+
  stat_summary(size=1) +
  ggtitle(paste("mcrA Abundances by Species"))+ 
  facet_wrap(core_type~., scales = "free_x", nrow=1)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](mcrA~(Copies~mg^-1))))+ 
  xlab("Species")+ 
 # scale_fill_manual(values=coul[23:32])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none")
print(mmox_bar)









flux_soil<-ggplot(data=ddpcr_meta_flux_justsoil,aes(x=log10(ddpcr_meta_flux_justsoil$mcra_probe_loose+1), 
                                                        y=log10(as.numeric(ddpcr_meta_flux_justsoil$ch4_50)))) +
  geom_point(size=4, aes(col=species))+
  ggtitle(paste("Soil mcrA Abundances vs. Tree Flux at 50cm"))+ 
  facet_wrap(ddpcr_meta_flux_justsoil$core_type~species, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm", col="black")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(flux_soil)

flux_soil_125<-ggplot(data=ddpcr_meta_flux_justsoil,aes(x=log10(ddpcr_meta_flux_justsoil$mcra_probe_loose+1), 
                                                    y=log10(as.numeric(ddpcr_meta_flux_justsoil$ch4_125)))) +
  geom_point(size=4, aes(col=species))+
  ggtitle(paste("Soil mcrA Abundances vs. Tree Flux at 125cm"))+ 
  facet_wrap(ddpcr_meta_flux_justsoil$core_type~species, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm", col="black")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mmoX~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(flux_soil_125)

flux_soil_200<-ggplot(data=ddpcr_meta_flux_justsoil,aes(x=log10(ddpcr_meta_flux_justsoil$mcra_probe_loose+1), 
                                                        y=log10(as.numeric(ddpcr_meta_flux_justsoil$ch4_200)))) +
  geom_point(size=4, aes(col=species))+
  ggtitle(paste("Soil mcrA Abundances vs. Tree Flux at 200cm"))+ 
  facet_wrap(ddpcr_meta_flux_justsoil$core_type~species, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm", col="black")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~Flux)))+ 
  xlab(expression(Log[10](mmoX~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(flux_soil_200)

cowplot::plot_grid(flux_soil, flux_soil_125, flux_soil_200, ncol=1)


ddpcr_meta_flux_justsoil=subset(ddpcr_all, material=="Soil")
ddpcr_meta_flux_justwood=subset(ddpcr_all, material=="Wood")

mmox_remerge=merge(ddpcr_meta_flux_justsoil, ddpcr_meta_flux_justwood, by="seq_id")
mmox_remerge=subset(mmox_remerge, species.x!="KALA"&
                                  species.x!="SAAL" &
                                  species.x!="QUAL" &
                                  species.x!="QUVE" &
                                  species.x!="CAOV" &
                                  species.x!="PRSE")



flux_soil_mmox<-ggplot(data=mmox_remerge,aes(x=log10(mcra_probe_loose.x+1), 
                                                        y=log10(mcra_probe_loose.y+1))) +
  geom_point(size=4, aes(col=species.x))+
  ggtitle(paste("Soil mmoX vs. Tree mmoX"))+ 
  facet_wrap(core_type.y~core_type.x*species.x, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](mmoX~(Copies~rxn^-1))~Tree))+ 
  xlab(expression(Log[10](mmoX~(Copies~rxn^-1))~Soil))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(flux_soil_mmox)

flux_soil_pmoa<-ggplot(data=mmox_remerge,aes(x=log10(pmoa_strict.x+1), 
                                             y=log10(pmoa_loose.y+1))) +
  geom_point(size=4, aes(col=species.x))+
  ggtitle(paste("Soil pmoA vs. Tree pmoA"))+ 
  facet_wrap(core_type.y~core_type.x*species.x, scales="free", ncol=10)+
  stat_cor()+
  geom_smooth(method="lm", se=FALSE, col="black")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))~Tree))+ 
  xlab(expression(Log[10](pmoA~(Copies~rxn^-1))~Soil))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(flux_soil_pmoa)



#### Import Internal Gas ####
internal_gas=read.csv("C:/Users/Wyatt/Desktop/Peat Bog/Tree/ddPCR/int_conc_prelim.csv")
internal_gas=subset(internal_gas, Tree.ID!="BLANK")        # Remove blank
internal_gas=subset(internal_gas, is.na(Tree.ID)==FALSE)   # Remove NA
internal_gas=subset(internal_gas, grepl("Amb",internal_gas$Tree.ID)==FALSE) # Remove ambient
internal_gas=subset(internal_gas, Tree.ID!="CHECK") # Remove ambient
internal_gas_clean=as.data.frame(cbind(internal_gas$Tree.ID, internal_gas$CH4.ppm, internal_gas$N2O.ppm, internal_gas$O2.ppm, internal_gas$CO2.ppm))
colnames(internal_gas_clean)=c("seq_id","CH4_int","N2O_int","O2_int","CO2_int")

ddpcr_all=merge(ddpcr_meta_flux,internal_gas_clean, by="seq_id", all.x=TRUE)

#### Internal gas graphs ####
ddpcr_all=subset(ddpcr_all, species!="KALA"&
                   species!="SAAL" &
                   species!="QUAL" &
                   species!="QUVE" &
                   species!="CAOV" &
                   species!="PRSE")

fint<-ggplot(data=ddpcr_all,aes(x=log10(as.numeric(ddpcr_all$mcra_probe_loose+1)), 
                                      y=log10(as.numeric(ddpcr_all$CH4_int)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("DNA Abundances vs. Internal Gas Concentration"))+ 
  facet_grid(ddpcr_all$core_type~., scales="free")+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(fint)


int_box=ggplot(data=ddpcr_all, aes(y=log10(as.numeric(ddpcr_all$CH4_int)), x=ddpcr_all$species, fill=species))+
  geom_jitter()+
  ggtitle("Internal Gas Concentrations")+
  geom_boxplot(alpha=0.8)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Tree~Species))
int_box
                 
          
fint_50<-ggplot(data=ddpcr_all,aes(x=log10(as.numeric(ddpcr_all$ch4_50)), 
                                y=log10(as.numeric(ddpcr_all$CH4_int)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("Flux vs. Internal Gas Concentration: 50cm"))+ 
  # facet_grid(ddpcr_all$core_type~.)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Log[10](CH[4]~Flux)))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(fint_50)

fint_125<-ggplot(data=ddpcr_all_inner,aes(x=log10(as.numeric(ddpcr_all_inner$ch4_125)), 
                                   y=log10(as.numeric(ddpcr_all_inner$CH4_int)), col=species)) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("Flux vs. Internal Gas Concentration: 125cm"))+ 
 # facet_wrap(.~species)+
  stat_cor()+
  geom_smooth(method="lm", se=FALSE)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Log[10](CH[4]~Flux)))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(fint_125)

ddpcr_all_inner=subset(ddpcr_all, core_type=="Inner")

fint_200<-ggplot(data=ddpcr_all_inner,aes(x=log10(as.numeric(ddpcr_all_inner$ch4_200)), 
                                   y=log10(as.numeric(ddpcr_all_inner$CH4_int)))) +
  geom_point(size=3, aes(col=species))+
  ggtitle(paste("Flux vs. Internal Gas Concentration: 200cm"))+ 
  facet_wrap(.~species)+
  stat_cor()+
  geom_smooth(method="lm")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Log[10](CH[4]~Flux)))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Species", values=coul[10:30])
print(fint_200)

cowplot::plot_grid(fint_50, fint_125, fint_200, ncol=3)
                           


vwc<-ggplot(data=ddpcr_all,aes(x=log10(as.numeric(ddpcr_all$ch4_125)), 
                               y=log10(as.numeric(ddpcr_all$CH4_int)))) +
  geom_point(data=ddpcr_all, aes(col=as.numeric(O2_int)), size=3)+
  ggtitle(paste("Flux vs. VWC"))+ 
  #  ylim(0,4.75)+
  ylab(expression(Log[10](CH[4]~ppm)))+ 
  xlab(expression(Log[10](CH[4]~Flux)))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_viridis_c(name="VWC %")
print(vwc)

st<-ggplot(data=ddpcr_all,aes(x=log10(as.numeric(ddpcr_all$ch4_125)), 
                              y=log10(as.numeric(ddpcr_all$CH4_int)))) +
  geom_point(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+1), 
                                  y=log10(pmoa_loose+1), col=mean_st), size=3)+
  ggtitle(paste("DNA Abundances vs. Soil Temperature"))+ 
  
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_continuous(name="°C")
print(st)

orp<-ggplot(data=ddpcr_all,aes(x=log10(as.numeric(ddpcr_all$ch4_125)), 
                               y=log10(as.numeric(ddpcr_all$CH4_int)))) +
  geom_point(data=ddpcr_meta, aes(x=log10(ddpcr_meta$mcra_probe_loose+1), 
                                  y=log10(pmoa_loose+1), col=mean_orp), size=3)+
  ggtitle(paste("DNA Abundances vs. ORP"))+ 
  
  facet_grid(.~ddpcr_meta$core_type)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1))))+ 
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  #  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_gradient(name="ORP", low="green", high="red")
print(orp)


#### More graphs ####

merge_no_blanks=subset(merge_no_blanks, merge_no_blanks$Accepted.Droplets>7500)
merge_no_blanks=subset(merge_no_blanks, Target!="16s_bact" & Target!="16s_arc")
#write.csv(merge_no_blanks,"ddPCR_tree_full_data.csv")
merge_no_blanks=subset(merge_no_blanks, species!="KALA"&
                              species!="SAAL" &
                              species!="QUAL" &
                              species!="QUVE" &
                              species!="CAOV" &
                              species!="PRSE")

#merge_no_blanks=subset(merge_no_blanks, species=="QURU"|
#                         species=="QUAL" |
#                         species=="QUVE" )




species_bar<-ggplot(data=merge_no_blanks, aes(x=species, y=log10(1+Copies.20ÂµLWell), col=Target)) +
  stat_summary(size=1.2)+
  ggtitle(paste("DNA Abundances"))+ 
  facet_wrap(.~core_type, scales = "free_x", ncol=2)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](Gene~Copies)))+ 
  xlab("Species ID")+ 
  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(species_bar)

species_bar<-ggplot(data=merge_no_blanks, aes(x=plate_identifier, y=log10(1+Copies.20ÂµLWell), col=Target, shape=core_type)) +
  stat_summary(size=0.7)+
  ggtitle(paste("DNA Abundances"))+ 
  facet_wrap(species~., scales = "free_x", nrow = 2)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](Gene~Copies)))+ 
  xlab("Sample")+ 
  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(species_bar)


ddpcr_trans=subset(ddpcr_trans, species!="KALA"&
                              species!="SAAL" &
                              species!="QUAL" &
                              species!="QUVE" &
                              species!="CAOV" &
                              species!="PRSE")

middle_bands<-ggplot(data=ddpcr_trans) +
  stat_summary(data=ddpcr_trans, aes(x=species, y=log10(mcra_probe_loose+1), col="loose"), size=1, na.rm = TRUE)+
  stat_summary(data=ddpcr_trans, aes(x=species, y=log10(mcra_probe_strict+1), col="strict"), size=1, na.rm = TRUE)+
  ggtitle(paste("DNA Abundances"))+ 
  facet_grid(.~ddpcr_trans$core_type, scales = "free", space="free")+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](Gene~(Copies~g~wood^-1))))+ 
  xlab("Species ID")+ 
  scale_color_discrete(name="Target")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(middle_bands)














ddpcr_trans_mcra_all=subset(ddpcr_trans, is.na(ddpcr_trans$mcra_eva_loose)==FALSE)
ddpcr_trans_mcra_all=subset(ddpcr_trans_mcra_all, species!="KALA"&
                                species!="SAAL" &
                                species!="QUAL" &
                                species!="QUVE" &
                                species!="CAOV" &
                                species!="PRSE")

species_bar<-ggplot(data=ddpcr_trans_mcra_all, aes(x=plate_identifier, y=log10(1+(mcra_probe_strict+mcra_probe_loose)/2), col="Probe")) +
  geom_point(size=3, aes(shape=core_type),alpha=0.6)+
  geom_point(aes(x=plate_identifier, y=log10(1+(mcra_eva_strict+mcra_eva_loose)/2), col="EvaGreen", shape=core_type), size=3, alpha=0.6)+
  # geom_boxplot()+
 # stat_summary(size=1.1)+
  ggtitle(paste("mcrA DNA"))+ 
  facet_wrap(.~ddpcr_trans_mcra_all$species, scales = "free_x", nrow=2)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  xlab("ID")+ 
  scale_color_manual(name="Fluorophore", values=wes_palette("Darjeeling1")[c(2,1)])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(species_bar)

species_bar<-ggplot(data=ddpcr_trans_mcra_all, aes(x=species, y=log10(1+(mcra_probe_strict)), col="Probe Strict")) +
  stat_summary(size=1)+
  stat_summary(aes(x=species, y=log10(1+(mcra_probe_loose)), col="Probe Loose"), size=1, shape=18, alpha=0.5)+
  stat_summary(aes(x=species, y=log10(1+(mcra_eva_loose)), col="EvaGreen Loose"), size=1, shape=18, alpha=0.5)+
  stat_summary(aes(x=species, y=log10(1+(mcra_eva_strict)), col="EvaGreen Strict"), size=1, alpha=0.5)+
  # geom_boxplot()+
  # stat_summary(size=1.1)+
  ggtitle(paste("mcrA DNA"))+ 
  facet_wrap(.~ddpcr_trans_mcra_all$core_type, scales = "free_x", nrow=1)+
  #  ylim(0,4.75)+
  ylab(expression(Log[10](mcrA~(Copies~rxn^-1))))+ 
  xlab("Species ID")+ 
  scale_color_manual(name="Cutoff Threshold \n (Mean SE)", values=wes_palette("Darjeeling1")[c(2,5,1,3)])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(species_bar)



species_bar<-ggplot(data=ddpcr_trans_mcra_all, aes(x=log10(1+mcra_eva_loose), y=log10(1+mcra_probe_loose))) +
  geom_point(size=3, shape=21, aes(fill=ddpcr_trans_mcra_all$species))+
  stat_cor() +
  geom_smooth(method="lm", formula = y~x, col="black")+
  ggtitle(paste("mcrA DNA Corr"))+ 
  facet_wrap(core_type~ddpcr_trans_mcra_all$species, scales = "free", nrow=2)+
  #  ylim(0,4.75)+
  ylab("Probe")+ 
  xlab("EvaGreen")+ 
 # scale_color_manual(name="Cutoff Threshold \n (Mean SE)", values=wes_palette("Darjeeling1"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none")
print(species_bar)





species_bar<-ggplot(data=ddpcr_trans, aes(x=log10(), y=log10(1+mmox_strict))) +
  geom_point(size=3, shape=21, aes(fill=ddpcr_trans$species))+
  stat_cor() +
  geom_smooth(method="lm", formula = y~x, col="black")+
  ggtitle(paste("mcrA DNA Corr"))+ 
  facet_wrap(.~material, scales = "free", nrow=2)+
  #  ylim(0,4.75)+
  ylab("Probe")+ 
  xlab("EvaGreen")+ 
  # scale_color_manual(name="Cutoff Threshold \n (Mean SE)", values=wes_palette("Darjeeling1"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position = "none")
print(species_bar)

#### LMER tests ####
ddpcr_inner= subset(ddpcr_all, core_type=="Inner")
ddpcr_outer = subset(ddpcr_all, core_type=="Outer")
ddpcr_org = subset(ddpcr_all, core_type=="Organic")
ddpcr_min = subset(ddpcr_all, core_type=="Mineral")


model=lm(log10(as.numeric(CH4_int))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_inner)
anova(model)

model=lm(log10(as.numeric(ch4_50))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_inner)
anova(model)



model=lm(log10(as.numeric(CH4_int))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_outer)
anova(model)

model=lm(log10(as.numeric(ch4_50))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_outer)
anova(model)



model=lm(log10(as.numeric(CH4_int))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_org)
anova(model)

model=lm(log10(as.numeric(ch4_50))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_org)
anova(model)



model=lm(log10(as.numeric(CH4_int))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_min)
anova(model)

model=lm(log10(as.numeric(ch4_50))~log10(mcra_probe_loose+1)+mean_vwc+log10(mmox_strict+1)+log10(pmoa_strict+1)+species, data=ddpcr_min)
anova(model)





lm_mod=blmer(scale(as.numeric(CO2_int)) ~ scale(mcra_probe_loose) + (1|species)+(0+scale(mcra_probe_loose)|species), data=ddpcr_inner)
lm_mod
anova(lm_mod)
parameters::p_value(lm_mod)

subset(ddpcr_inner, is.na(mcra_probe_loose)==FALSE & is.na(CO2_int)==FALSE) %>% 
  # save predicted values
  mutate(pred_ch4 = fitted(lm_mod)) %>% 
  # graph
  ggplot(aes(x=scale(mcra_probe_loose), y=pred_ch4, group=species, color=species)) + theme_classic() +
  geom_line(size=2) 


