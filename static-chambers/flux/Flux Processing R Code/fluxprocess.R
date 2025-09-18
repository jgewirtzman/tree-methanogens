
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
par(mfrow=(c(3,4)))
for(i in 1:length(fluxes)){
  plot(result1$`[CH4]d_ppm`[which(result1$fluxid==fluxes[i])]~
         result1$date_time[which(result1$fluxid==fluxes[i])],
       xlab="Time", ylab="[CH4] Dry", main=fluxes[i])
#  abline(lm(
#    result1$`[CH4]d_ppm`[which(result1$fluxid==fluxes[i])]~
#      result1$date_time[which(result1$fluxid==fluxes[i])]
#  ))
}
dev.off()


####VOLUME AND AREA#####
result1$area<-c(rep(1.5, nrow(result1)))
result1$vol<-c(rep(0.75, nrow(result1)))



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

write.csv(result1, "result1.csv")

#result1a<-read.csv("result1_new.csv")

result1a<-result1

result1a$CH4Code<-c(rep(0, nrow(result1a)))
result1a$ch4ppb<-result1a$ch4dppm*1000

#fluxes.chopped<-chop(dat=result1a, factors=c('tree', "height"), 
#                     nmes = c("tree", "height"))

x<-summary(lm(result1a$ch4ppb~result1a$tmin))
plot(result1a$ch4ppb~result1a$tmin)
abline(lm(result1a$ch4ppb~result1a$tmin))


lgrfn$ch4_slope<-c()
lgrfn$ch4_r2<-c()
lgrfn$ch4_pval<-c()

lgrfn$co2_slope<-c()
lgrfn$co2_r2<-c()
lgrfn$co2_pval<-c()


sum(is.na(result1a$tmin))
result1a<-result1a[which(is.na(result1a$tmin)==F),]
result1a<-result1a[which(is.na(result1a$tmin)==F),]


fluxes<-fluxes[fluxes %in% result1a$fluxid]


#n_obs<-c()
#for(i in 1:length(fluxes)){
#  n_obs[i]<-(nrow(result1a[which(result1a$fluxid==fluxes[i]),]))
#}

#result1a<-result1a[which(n_obs>3),]



plot(ch4_r2~co2_r2)
plot(co2_slopes~ch4_slopes)


output2<-left_join(outputs, LGR_fn, by="fluxid")
summary(lm(output2$ch4_slopes~output2$Measurement.height+output2$Tree_sp))



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
chamber_vol[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-0.3*0.001
chamber_vol[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-0.75*0.001
chamber_vol[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-2.3*0.001

#in m2
chamber_area[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-14.16*9.69*6.31*0.0001
chamber_area[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-20.04*13.6*6.31*0.0001
chamber_area[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-28.07*19.3*7.87*0.0001 


AtmPress<-c(rep(101.325, nrow(output2))) # If Atmospheric Pressure is not available at the measuring time, we can put a constant value (kg m-2 s-1)
glc<-c(rep(0.00831447, nrow(output2))) # ideal gas law constant in kg m2 micromol-1 k-1

output2$chamber_vol<-chamber_vol
output2$chamber_area<-chamber_area
output2$glc<-glc
output2$pressure<-AtmPress


output2$CH4_lin_flux<-
  output2$ch4_slopes*(output2$chamber_vol/output2$chamber_area)*
  (output2$pressure/(output2$glc*(output2$Temp_air+273)))*1000 #nmol m-2 s-1



#####Plots 


library(ggplot2)


ggplot(data=output2, aes(x=Tree_sp, y=CH4_lin_flux))+
  geom_jitter()

ggplot(data=output2, aes(x=as.factor(Measurement.height), y=CH4_lin_flux))+
  geom_point(aes(col=Tree_sp))+
  #geom_jitter(aes(col=Tree_sp))+
  geom_violin()+
  coord_flip()

  facet_wrap(vars(Tree_sp))
  #ylim(-2, 10)+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")

ggplot(data=output2, aes(x=Measurement.height, y=CH4_lin_flux, fill=as.factor(Measurement.height)))+
  geom_point()+
  #geom_smooth()+
  geom_boxplot()+
  facet_wrap(vars(Tree_sp))+
  #ylim(-2, 10)+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")
  

output3<-output2[which(output2$Tree_sp!="BEAL" & output2$Tree_sp!="FRAM"),]
ggplot(data=output3, aes(x=Measurement.height, y=CH4_lin_flux, 
                         fill=as.factor(Measurement.height)))+
  geom_jitter()+
  #geom_smooth()+
  geom_boxplot()+
  facet_wrap(vars(Tree_sp))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-500, 1500)



mean(output2$CH4_lin_flux, na.rm=T)
mean(output3$CH4_lin_flux, na.rm=T)








fluxes[16]
check<-result1a[which(result1a$fluxid=="BB9_125"),]
View(check)



for(i in 1:length(fluxes)){
  plot(result1a$`[CH4]d_ppm`[which(result1a$fluxid==fluxes[i])]
       ~result1a$tmin[which(result1a$fluxid==fluxes[i])])
  abline(lm(result1a$`[CH4]d_ppm`[which(result1a$fluxid==fluxes[i])]~
              result1a$tmin[which(result1a$fluxid==fluxes[i])]))
  






x<-summary(lm(result1a$ch4ppb~result1a$tmin))

x$coefficients[1]
x$coefficients[4]
x$r.squared
