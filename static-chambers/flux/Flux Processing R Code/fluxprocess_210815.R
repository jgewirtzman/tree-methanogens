
# Linear and exponential fluxes calculated from LGR measurements    (created by jbarba at 09262018)
## This script contents data compilation from different LGR files and fluxes processing using linear and exponential fits

# Load packages
library(data.table)
library(lubridate)
library(doBy)
library(dplyr)
library(padr)
library(zoo)
library(sqldf)
library(flux)
library(ggplot2)
library(cowplot)



# Merging LGR files -----------------------------------------------------------

## With this chunk of script we can make one data table from a list of many. It reads all files from a folder and their Tree subfolders and merge them.
# It does not read compressed files

# LGR sometimes saves files in .zip. The first part of this chunk is for decompress all .zip files from a given folder
# List all .zip files including sub-folders
list_of_zip_files <- list.files(path = "/Users/jongewirtzman/Google Drive/YMF Methane/Data/2021 raw flux data", recursive=TRUE, 
                                pattern="\\.zip$", full.names=TRUE)

# Decompress all .zip files from a given folder and subfolders
sapply(list_of_zip_files, function(i) unzip(i, exdir=gsub("\\.zip$", "", i))) #This makes a copy of the uncompressed folder in the same folder as .zip was

# List all txt files, merge them in one single file and create a variable with the name of each original file 
# List all txt files
list_of_txt_files <- list.files(path = "/Users/jongewirtzman/Google Drive/YMF Methane/Data/2021 raw flux data", recursive = TRUE,
                                pattern = "\\.txt$", full.names = T)  #with full.name=F the function save the name of each file instead of the name of each path. This is useful for the idcol in the next section 

# Read all the files and create a Path column to store filenames
LGR_data <- rbindlist(sapply(list_of_txt_files, fread, simplify = FALSE),
                      use.names = TRUE, idcol = "Path", fill=T)

LGR_data<- as.data.frame(LGR_data) 

LGR_data<-LGR_data[!duplicated(LGR_data$date_time), ]

#plot(LGR_data$`[CH4]d_ppm`~LGR_data$date_time, ylim=c(0,5))



###Write Files###
#write.table(LGR_data,file="/Users/jongewirtzman/Google Drive/YMF Methane/Data/2021 raw flux data/LGR_data.txt",sep="\t")   
#write.csv(LGR_data,file="/Users/jongewirtzman/Google Drive/YMF Methane/Data/2021 raw flux data/LGR_data.csv")   



## LGR_field_notes
# Create new variables for meteo, diameter, collar height

setwd("/Users/jongewirtzman/Google Drive/YMF Methane/Data")                        
LGR_fn <- read.table("Gewirtzman_datasheets_fixing_dates.csv", fill=TRUE,header=TRUE,sep=",",na.strings=c("NA","#N/A!","#N/A","#NA!"))    


LGR_fn$Tree_sp[which(LGR_fn$Tree_sp=="TSLA")]<-"TSCA"


LGR_fn$Date_full<-paste("20", LGR_fn$Date, sep="")
LGR_fn$date.parsed<-ymd(LGR_fn$Date_full)

LGR_fn$starttime<-ymd_hms(paste(LGR_fn$date.parsed, LGR_fn$HMS_Start))

LGR_fn$endtime<-ymd_hms(paste(LGR_fn$date.parsed, LGR_fn$HMS_End))


LGR_fn$fluxid<-paste(LGR_fn$Tree_ID, LGR_fn$Measurement.height, sep="_")




LGR_fn.labels<-LGR_fn[,21:23]
#LGR_data.subset<-LGR_data[20000:25000,]
#lgrds<-LGR_data.subset[,c(1,40)]

#LGR_fn.labels$start.time.final<-as_datetime(LGR_fn.labels$start.time.final)
#LGR_fn.labels$end.time.final<-as_datetime(LGR_fn.labels$end.time.final)

lgrfn<-LGR_fn.labels


names(lgrfn)<-c('start', 'end', 'fluxid')

#LGR_data_key<-LGR_data[,c(1:34, 36)]
LGR_data_key <- LGR_data[, !duplicated(colnames(LGR_data))]


str(LGR_data)


result = sqldf("select * from LGR_fn
                left join LGR_data_key
                on LGR_data_key.date_time between LGR_fn.starttime and LGR_fn.endtime")


str(result)
result1<-result[order(result$date_time), ]

fluxes<-unique(result1$fluxid)


####Plot###
par(mfrow=(c(3,2)))
for(i in 1:length(fluxes)){
  plot(result1$`[CH4]d_ppm`[which(result1$fluxid==fluxes[i])]~
         result1$date_time[which(result1$fluxid==fluxes[i])],
       xlab="Time", ylab="[CH4] Dry", main=fluxes[i])
 # abline(lm(
  #  result1$`[CH4]d_ppm`[which(result1$fluxid==fluxes[i])]~
  #    result1$date_time[which(result1$fluxid==fluxes[i])]
  #))
  
  plot(result1$`[CO2]d_ppm`[which(result1$fluxid==fluxes[i])]~
         result1$date_time[which(result1$fluxid==fluxes[i])],
       xlab="Time", ylab="[CO2] Dry", main=fluxes[i])
}



#dev.off()



#time in min
result1$tmin<-(result1$date_time-result1$starttime)/60

#convert names
result1$ch4dppm<-result1$`[CH4]d_ppm`


######TESTING ONLY
#result1<-result1[which(result1$fluxid=="AB5_200"),]
#############


#remove unneeded columns#
#result1a<-result1[,c('fluxid', 'ch4dppm', 'tmin', 'vol', 'Temp_air', 'area')]
result1a<-result1


result1a$tmin<-as.numeric(result1a$tmin)
rownames(result1a)<-NULL

#write.csv(result1, "result1.csv")

#result1a<-read.csv("result1_new.csv")

#result1a<-result1

#result1a$CH4Code<-c(rep(0, nrow(result1a)))
#result1a$ch4ppb<-result1a$ch4dppm*1000

#fluxes.chopped<-chop(dat=result1a, factors=c('tree', "height"), 
#                     nmes = c("tree", "height"))

#x<-summary(lm(result1a$ch4ppb~result1a$tmin))
#plot(result1a$ch4ppb~result1a$tmin)
#abline(lm(result1a$ch4ppb~result1a$tmin))


sum(is.na(result1a$tmin))
result1a<-result1a[which(is.na(result1a$tmin)==F),]
result1a<-result1a[which(is.na(result1a$tmin)==F),]


fluxes<-fluxes[fluxes %in% result1a$fluxid]


#n_obs<-c()
#for(i in 1:length(fluxes)){
#  n_obs[i]<-(nrow(result1a[which(result1a$fluxid==fluxes[i]),]))
#}

#result1a<-result1a[which(n_obs>3),]

result1a$tsec<-result1a$tmin*60


ch4_slopes<-c()
co2_slopes<-c()
ch4_r2<-c()
co2_r2<-c()
for(i in c(1:15, 17:411)){
  ch4_slopes[i]<-
    coefficients(
      summary(lm(result1a$`[CH4]d_ppm`[which(result1a$fluxid==fluxes[i])]~
                   result1a$tsec[which(result1a$fluxid==fluxes[i])]))
    )[2,1]
  
  co2_slopes[i]<-
    coefficients(
      summary(lm(result1a$`[CO2]d_ppm`[which(result1a$fluxid==fluxes[i])]~
                   result1a$tsec[which(result1a$fluxid==fluxes[i])]))
    )[2,1]
  
  ch4_r2[i]<-
    summary(lm(result1a$`[CH4]d_ppm`[which(result1a$fluxid==fluxes[i])]~
                 result1a$tsec[which(result1a$fluxid==fluxes[i])]))$r.squared
  
  co2_r2[i]<-
    summary(lm(result1a$`[CO2]d_ppm`[which(result1a$fluxid==fluxes[i])]~
                 result1a$tsec[which(result1a$fluxid==fluxes[i])]))$r.squared
  
  
}

outputs<-cbind(fluxes, ch4_slopes, co2_slopes, ch4_r2, co2_r2)
outputs<-as.data.frame(outputs)
names(outputs)[1]<-"fluxid"



#plot(ch4_r2~co2_r2)
#plot(co2_slopes~ch4_slopes)


output2<-left_join(outputs, LGR_fn, by="fluxid")
#summary(lm(output2$ch4_slopes~output2$Measurement.height+output2$Tree_sp))



str(output2)
output2$ch4_slopes<-as.numeric(output2$ch4_slopes)
output2$Measurement.height<-as.numeric(output2$Measurement.height)
output2$Tree_sp<-as.factor(output2$Tree_sp)
output2<-output2[which(is.na(output2$Measurement.height)==F),]
output2$Measurement.height[which(output2$Measurement.height==75)]<-50



######calcualte fluxes from slopes######
chamber_area<-c()
chamber_vol<-c()

unique(output2$Chamber.ID)
output2$Chamber.ID[which(is.na(output2$Chamber.ID)==T)]<-6

#in m3
chamber_vol[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-0.3*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002
chamber_vol[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-0.75*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002
chamber_vol[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-2.3*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002

#in m2
chamber_area[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-14.16*9.69*0.0001
chamber_area[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-20.04*13.6*0.0001
chamber_area[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-28.07*19.3*0.0001 


AtmPress<-c(rep(101.325, nrow(output2))) # If Atmospheric Pressure is not available at the measuring time, we can put a constant value (kg m-2 s-1)
#pressure in kPa
glc<-c(rep(0.00831447, nrow(output2))) # ideal gas law constant in kg m2 micromol-1 k-1


output2$chamber_vol<-chamber_vol
output2$chamber_area<-chamber_area
output2$glc<-glc
output2$pressure<-AtmPress


output2$CH4_lin_flux<-
  output2$ch4_slopes*(output2$chamber_vol/output2$chamber_area)*
  (output2$pressure/(output2$glc*(output2$Temp_air+273)))*1000 #nmol m-2 s-1


output2$CO2_lin_flux<-
  as.numeric(output2$co2_slopes)*(output2$chamber_vol/output2$chamber_area)*
  (output2$pressure/(output2$glc*(output2$Temp_air+273))) #micromol m-2 s-1



####




#####Plots 

output2<-output2[-which(output2$fluxid=="BO1_200"),]



summary(lme(CH4_lin_flux~Measurement.height, random=~1|Tree_ID, 
    data=output3, na.action = na.omit))



###

setwd("/Users/jongewirtzman/Google Drive/YMF Methane/figs")


###


p1<-ggplot(data=output2, aes(x=Tree_sp, y=CH4_lin_flux,
                         color=Tree_sp))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  geom_jitter()+
  xlab("Tree Species")+ylab("CH4 Flux, nmol m-2 s-1")+theme_cowplot()+ theme(legend.position = "none")+
  #geom_hline(yintercept=mean(output2$CH4_lin_flux, na.rm=T), linetype="dashed")+
  #geom_hline(yintercept=-2, linetype="dashed")+
  geom_hline(yintercept=0)#+
  #geom_hline(yintercept=-2+mean(output2$CH4_lin_flux, na.rm=T), linetype="dotted")
ggsave("fig1.png", p1)





p2<-ggplot(data=output2, aes(x=Tree_sp, y=CH4_lin_flux,
                             color=Tree_sp))+
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean")+
  geom_jitter()+
  geom_hline(yintercept=mean(output2$CH4_lin_flux, na.rm=T), linetype="dashed")+
  geom_hline(yintercept=-2, linetype="dashed")+theme_cowplot()+theme(legend.position = "none")+
  xlab("Tree Species")+ylab("CH4 Flux, nmol m-2 s-1")+
  geom_hline(yintercept=0)#+
ggsave("fig2.png", p2)



p3<-
ggplot(data=output2, aes(x=as.factor(Measurement.height), y=CH4_lin_flux))+
  #geom_point(aes(col=Tree_sp))+
  geom_jitter(aes(col=Tree_sp))+theme_cowplot()+
  #geom_violin()+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()
ggsave("fig3.png", p3)




p4<-
ggplot(data=output2, aes(x=Measurement.height, y=CH4_lin_flux, fill=as.factor(Measurement.height)))+
  geom_point()+
  #geom_smooth()+
  geom_boxplot()+
  facet_wrap(vars(Tree_sp))+
  #ylim(-2, 10)+
  coord_flip()+
  xlab("Stem Height")+ylab("CH4 Flux, nmol m-2 s-1")+theme_cowplot()+theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")
ggsave("fig4.png", p4)

  

output3<-output2[which(output2$Tree_sp!="BEAL" & output2$Tree_sp!="FRAM"),]


p5<-
ggplot(data=output3, aes(x=as.factor(Measurement.height), y=CH4_lin_flux))+
  #geom_point(aes(col=Tree_sp))+
  geom_jitter(aes(col=Tree_sp))+
  #geom_violin()+
  geom_boxplot(outlier.shape = NA)+theme_cowplot()+theme(legend.position = "none")+
  xlab("Stem Height")+ylab("CH4 Flux, nmol m-2 s-1")+
  coord_flip()
  ggsave("fig5.png", p5)


p6<-
  ggplot(data=output3, aes(x=Measurement.height, y=CH4_lin_flux, 
                         fill=as.factor(Measurement.height)))+
  geom_jitter()+
  #geom_smooth()+
  geom_boxplot()+
  facet_wrap(vars(Tree_sp))+
  coord_flip()+theme_cowplot()+theme(legend.position="none")+
  geom_hline(yintercept=0, linetype="dashed")+
  xlab("Stem Height")+ylab("CH4 Flux, nmol m-2 s-1")+
  ylim(-.01, 0.3)
ggsave("fig6.png", p6)


p7<-
ggplot(data=output2, aes(x=as.numeric(DBH), y=CH4_lin_flux, 
                         color=as.factor(Tree_sp)))+
  geom_jitter()+
  geom_smooth(method="lm")+theme_cowplot()+
  xlab("Stem Height")+ylab("CH4 Flux, nmol m-2 s-1")
ggsave("fig7.png", p7)


p8<-
ggplot(data=output3, aes(x=as.numeric(DBH), y=CH4_lin_flux, 
                         color=as.factor(Tree_sp)))+
  geom_jitter()+theme_cowplot()+
  xlab("DBH (cm)")+ylab("CH4 Flux, nmol m-2 s-1")+theme(legend.position="none")
  #geom_smooth(method="lm")
ggsave("fig8.png", p8)


p9<-
ggplot(data=output3, aes(x=as.factor(Measurement.height), y=CH4_lin_flux))+
  #geom_point(aes(col=Tree_sp))+
  geom_violin()+
  geom_jitter(aes(col=Tree_sp))+theme_cowplot()+
  coord_flip()+
  geom_hline(yintercept=mean(output3$CH4_lin_flux, na.rm=T), linetype="dashed")
ggsave("fig9.png", p9)



###
  
  
  soils<-read.csv("/Users/jongewirtzman/Downloads/Soils Data.csv")
  names(soils)[1]<-"Tree_ID"
  
  alldata<-left_join(output2, soils)
  
  
  ###
  
  p10<-
ggplot(data=alldata, aes(x=Average.of.VWC..Mean., y=CH4_lin_flux,
                         color=Tree_sp))+theme_cowplot()+
  geom_jitter()+
  ylim(-0.1, 1)
  #geom_smooth(method="lm")
  ggsave("fig10.png", p10)
  

p11<-
ggplot(data=alldata, aes(x=Average.of.ORP..mV., y=CH4_lin_flux,
                         color=Tree_sp))+theme_cowplot()+
  geom_jitter()

ggsave("fig11.png", p11)


p12<-
ggplot(data=alldata, aes(x=CO2_lin_flux, y=CH4_lin_flux))+
  geom_jitter(aes(color=Tree_sp))+
  xlim(0, 7.5)+theme_cowplot()+
  geom_smooth(method="lm")
ggsave("fig12.png", p12)

  
p13<-
ggplot(data=alldata2, aes(x=CO2_lin_flux, y=CH4_lin_flux,
                         color=Tree_sp))+
  geom_jitter()+theme_cowplot()+
  geom_smooth(method="lm")
ggsave("fig13.png", p13)



  alldata2<-alldata[which(output2$Tree_sp!="BEAL" & output2$Tree_sp!="FRAM"),]

  
  p14<-
  ggplot(data=alldata2, aes(x=Average.of.VWC..Mean., y=CH4_lin_flux
                           ))+
    geom_jitter(aes(color=Tree_sp))+theme_cowplot()+
    xlab("Soil Moisture (%VWC)")+ylab("CH4 Flux, nmol m-2 s-1")+theme(legend.position="none")
    #geom_smooth(method="lm")
  ggsave("fig14.png", p14)
  
  p142<-
    ggplot(data=alldata2, aes(x=Average.of.Soil.Temp, y=CH4_lin_flux
    ))+
    geom_jitter(aes(color=Tree_sp))+theme_cowplot()+
    xlab("Soil Temperature (C)")+ylab("CH4 Flux, nmol m-2 s-1")+theme(legend.position="none")#+
  #  geom_smooth(method="lm")
  ggsave("fig142.png", p142)
  

  p15<-
  ggplot(data=alldata2, aes(x=Average.of.ORP..mV., y=CH4_lin_flux))+
    geom_jitter(aes(color=Tree_sp))+
    xlab("Redox Potential (mV)")+ylab("CH4 Flux, nmol m-2 s-1")+theme_cowplot()+theme(legend.position="none")
    #    geom_smooth(method="lm")
    #ylim(-0.01, 0.1)
  ggsave("fig15.png", p15)
  
  p122<-
    ggplot(data=alldata2, aes(x=CO2_lin_flux, y=CH4_lin_flux))+
    geom_jitter(aes(color=Tree_sp))+
    xlim(0, 7.5)+
    geom_smooth(method="lm")+theme_cowplot()
  ggsave("fig122.png", p12)
  
  
  
  
  p16<-
  ggplot(data=alldata, aes(x=Tree_sp, y=CH4_lin_flux))+
    geom_boxplot()+theme_cowplot()
  ggsave("fig16.png", p16)
  
  
  p17<-
    ggplot(data=alldata2, aes(x=Tree_sp, y=CH4_lin_flux))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter()+theme_cowplot()
  ggsave("fig17.png", p17)
  
  
  p18<-
  ggplot(data=alldata2, aes(x=as.factor(Measurement.height), y=CH4_lin_flux))+
    geom_boxplot()+
    geom_jitter()+theme_cowplot()
  ggsave("fig18.png", p18)
  
  
  
  plot(alldata$CO2_lin_flux~alldata$chamber_vol)
  
  
  
  ####


wood<-c()
wood[which(output2$Tree_sp=="TSCA"|output2$Tree_sp=="PIST")]<-"Softwood"
wood[which(output2$Tree_sp!="TSCA"&output2$Tree_sp!="PIST")]<-"Hardwood"
output2$wood<-wood

p19<-
ggplot(data=output2, aes(x=wood, y=CH4_lin_flux))+
  geom_jitter(aes(color=Tree_sp))+
  geom_boxplot(outlier.shape = NA)+xlab("")+ylab("CH4 Flux, nmol m-2 s-1")+theme_cowplot()
ggsave("fig19.png", p19)

p192<-
  ggplot(data=output3, aes(Temp_stem, y=CH4_lin_flux))+
  geom_point(aes(color=Tree_sp))+
  xlab("Stem Temperature (C)")+ylab("CH4 Flux, nmol m-2 s-1")+theme_cowplot()+theme(legend.position="none")
ggsave("fig192.png", p192)


p20<-
ggplot(data=alldata2, aes(x=wood, y=CH4_lin_flux))+
  geom_jitter(aes(color=Tree_sp))+
  xlab("")+ylab("CH4 Flux, nmol m-2 s-1")+
  geom_boxplot(outlier.shape = NA, fill=alpha(0.25))+theme_cowplot()+theme(legend.position="none")
ggsave("fig20.png", p20)


p21<-
ggplot(data=alldata2, aes(x=Temp_air, y=CH4_lin_flux))+
  geom_point()+theme_cowplot()
ggsave("fig21.png", p21)


p22<-
ggplot(data=alldata, aes(x=Average.of.Organic.Sample.Depth..cm., y=CH4_lin_flux))+
  geom_point()+theme_cowplot()
ggsave("fig22.png", p22)




lm<-lm(alldata$CH4_lin_flux~
         alldata$CO2_lin_flux+alldata$Average.of.VWC..Mean.+
         alldata$Average.of.ORP..mV.+alldata$Average.of.Organic.Sample.Depth..cm.+
         alldata$Tree_sp+alldata$Measurement.height+alldata$Average.of.Soil.Temp+
         alldata$Temp_stem+alldata$Temp_air+as.numeric(alldata$DBH))
summary(lm)

lm2<-lm(alldata$CH4_lin_flux~
         alldata$CO2_lin_flux+alldata$Average.of.VWC..Mean.+
         alldata$Average.of.ORP..mV.+alldata$Average.of.Organic.Sample.Depth..cm.+
         alldata$Measurement.height+alldata$Average.of.Soil.Temp+
         alldata$Temp_stem+alldata$Temp_air+as.numeric(alldata$DBH))
summary(lm2)


lm3<-lm(alldata$CH4_lin_flux~
          alldata$Average.of.Soil.Temp*alldata$Tree_sp)
summary(lm3)



mean(output2$CH4_lin_flux[which(output2$Tree_sp=="PIST")])
mean(output2$CH4_lin_flux[which(output2$Tree_sp=="TSCA")])
mean(output2$CH4_lin_flux[which(output2$Tree_sp=="QURU")])
mean(output2$CH4_lin_flux[which(output2$Tree_sp=="BELE")], na.rm=T)
mean(output2$CH4_lin_flux[which(output2$Tree_sp=="ACRU")])



####

library(ggrepel)


ef<-read.csv("/Users/jongewirtzman/Dropbox (Yale_FES)/1 Eastern US internal CH4/From Covey 9_21/easternforest_clean_wmeta.csv")

str(ef$ch4)

mean(ef$ch4[which(ef$spcode=="PIST")])
#length(ef$ch4[which(ef$spcode=="PIST")])

mean(ef$ch4[which(ef$spcode=="TSCA")])
#length(ef$ch4[which(ef$spcode=="TSCA")])

mean(ef$ch4[which(ef$spcode=="QURU")])
#length(ef$ch4[which(ef$spcode=="QURU")])

mean(ef$ch4[which(ef$spcode=="BELE")])
#length(ef$ch4[which(ef$spcode=="BELE")])

mean(ef$ch4[which(ef$spcode=="ACRU")])
#length(ef$ch4[which(ef$spcode=="ACRU")])

ef$ch4

figx<-
ggplot(data=output3,
       aes(y=CH4_lin_flux, x=Measurement.height))+
  geom_point()+
  facet_wrap(vars(Tree_ID))+
  coord_flip()+theme_cowplot()
ggsave("figx.png", figx)


  
spmeans<-
  ef %>%
  group_by(spcode) %>%
  summarise(meanch4 = mean(ch4), 
            medianch4=median(ch4),
            meano2=mean(o2, na.rm=T),
            mediano2=median(o2, na.rm = T),
            n = n())

spmean_flux<-
  output2%>%
  group_by(Tree_sp) %>%
  summarise(meanch4_flux=mean(CH4_lin_flux, na.rm=T),
            medianch4_flux=median(CH4_lin_flux, na.rm=T),
            n=n())

names(spmean_flux)[1]<-"spcode"
spmean_flux$spcode<-as.character(spmean_flux$spcode)

int_c<-merge(spmean_flux, spmeans, by="spcode")



species<-unique(output2$Tree_sp)
ef<-ef[which(ef$spcode %in% species),]

dev.off()
plot(int_c$meanch4_flux~int_c$meanch4)

plot(int_c$medianch4_flux~int_c$medianch4)

p23<-
ggplot(data=ef, aes(x=spcode, y=log(ch4), color=spcode))+
  geom_jitter()+
  geom_boxplot()+theme_cowplot()
ggsave("fig23.png", p23)

p24<-
ggplot(data=ef, aes(x=spcode, y=o2, color=spcode))+
  geom_jitter()+
  geom_boxplot()+theme_cowplot()
ggsave("fig24.png", p24)

p25<-
ggplot(data=ef, aes(x=o2, y=log(ch4), color=spcode))+
  geom_jitter()+
  geom_smooth(method="lm")+theme_cowplot()
ggsave("fig25.png", p25)

p26<-
ggplot(data=ef, aes(x=o2, y=log(ch4)))+
  geom_jitter()+
  geom_smooth(method="lm")+theme_cowplot()
ggsave("fig26.png", p26)


p27<-
ggplot(data=int_c, aes(x=medianch4, y=medianch4_flux,
                       label=spcode))+
  geom_point(aes(color=spcode))+
  geom_smooth(method="lm")+
  geom_text_repel(aes(label = spcode),
                  segment.color = 'grey50')+theme_cowplot()
ggsave("fig27.png", p27)




p28<-
ggplot(data=int_c, aes(x=meanch4_flux, y=medianch4_flux,
                       label=spcode))+
  geom_point(aes(color=spcode))+
  geom_smooth(method="lm")+
  geom_text_repel(aes(label = spcode),
                  segment.color = 'grey50')+theme_cowplot()
ggsave("fig28.png", p28)

p28<-
  ggplot(data=int_c, aes(y=meanch4_flux, x=mediano2/10000,
                         label=spcode))+
  geom_point(aes(color=spcode))+
  geom_smooth(method="lm")+
  geom_text_repel(aes(label = spcode),
                  segment.color = 'grey50')+
  xlab("Species Median Trunk Oxygen (%)")+ylab("Species Median CH4 Flux, nmol m-2 s-1")+theme_cowplot()+theme(legend.position="none")
  #geom_vline(xintercept=20.95)
ggsave("fig28.png", p28)



###

library(ggmap)

map<-get_stamenmap(bbox = c(left = -70,
                       bottom = 33,
                       right = -93,
                       top = 45),
              maptype = "terrain", 
              crop = FALSE,
              zoom = 7)
m1<-
ggmap(map)+
  geom_point(data=ef, aes(x=long, y=lat2))
ggsave("m1.png", m1)


map2<-get_stamenmap(bbox = c(left = -72.08342248359584-5,
                             bottom = 41.95937699237806-5,
                             right = -72.08342248359584+5,
                             top = 41.95937699237806+5),
                    maptype = "terrain", 
                    crop = FALSE,
                    zoom = 7)

ymflat<-c(41.95937699237806)
ymflong<-c(-72.08342248359584)
label<-c("Yale-Myers Forest")
mappt<-as.data.frame(cbind(ymflat, ymflong, label))
mappt[,1:2]<-as.numeric(as.character(mappt[,1:2]))

m2<-
ggmap(map2)+
  geom_point(data=mappt, aes(x=ymflong, y=ymflat), color="yellow")+
  geom_label_repel(data=mappt, aes(x=ymflong, y=ymflat,
                       label = "Yale-Myers Forest"),
                  segment.color = 'grey50')
ggsave("m2.png", m2)

latlong<-strsplit(ef$latlong, ",")
ef$long<-as.numeric(substr(lapply(latlong, `[[`, 1), start=2, stop=100))
ef$lat2<-as.numeric(substr(lapply(latlong, `[[`, 2), start=1, stop=4))
min(ef$long)

length(which(output2$ch4_slopes>0))/length(output2$ch4_slopes)

library(ggpubr)
pg<-
ggarrange(ncol=3,nrow=2,
          p14,p15,p142,p20,p8,p192)
ggsave("pg.png", pg)

length(output2$fluxid)
length(which(output2$ch4_slopes>0))
