# Linear and exponential fluxes calculated from LGR measurements    (modified from scrip by josep barba 09262018)
#Modified by Jon Gewirtzman over 2022-2023
## This script performs data compilation from different LGR files and fluxes processing using linear fits

# Load packages -----------------------------------------------------------

#function to load packages, after first installing if not yet installed
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[, 1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

# Load packages
usePackage("data.table")
usePackage("lubridate")
usePackage("doBy")
usePackage("dplyr")
usePackage("padr")
usePackage("zoo")
usePackage("sqldf")
usePackage("flux")
usePackage("ggplot2")
usePackage("cowplot")
usePackage("data.table")


# Merging LGR files -----------------------------------------------------------

#takes all LGR files of interest and stitches them into a single chronologically-ordered dataframe/CSV

# data is stored on the LGR instrument in dated folders, e.g. "2023-01-27".
#these folders should be dragged directly, as-is, to a single folder.
#make a back-up copy before running this code, as it will modify folder structure.
#the name of the folder does not matter. 
#no need to unzip subfolders, etc.; this script will do that

###

# Set file locations -----------------------------------------------------------

#set working directory
#this should be the folder that contains all LGR files, within dated folders
wd_folder<-"/Users/jongewirtzman/Downloads/LGR_files_25Jan_12feb (copy)"
setwd(wd_folder)

#set output directory
#this is where your final text file/CSV will be saved
output_location<-"/Users/jongewirtzman/Downloads/LGR_data_compiled"

# Merging LGR files -----------------------------------------------------------

## With this chunk of script we can make one data table from a list of many. It reads all files from a folder and their subfolders and merge them.
# It does not read compressed files
# LGR sometimes saves files in .zip. The first part of this chunk is for decompress all .zip files from a given folder

# List all .zip files including sub-folders
list_of_zip_files <- list.files(path = wd_folder, recursive=TRUE, 
                                pattern="\\.zip$", full.names=TRUE)

# Decompress all .zip files from a given folder and subfolders
sapply(list_of_zip_files, function(i) unzip(i, exdir=gsub("\\.zip$", "", i))) #This makes a copy of the uncompressed folder in the same folder as .zip was

# List all txt files, merge them in one single file and create a variable with the name of each original file 
# List all txt files
list_of_txt_files <- list.files(path = wd_folder, recursive = TRUE,
                                pattern = "\\.txt$", full.names = T)  #with full.name=F the function save the name of each file instead of the name of each path. This is useful for the idcol in the next section 

# Read all the files and create a Path column to store filenames
LGR_data <- rbindlist(sapply(list_of_txt_files, fread, simplify = FALSE),
                      use.names = TRUE, idcol = "Path", fill=T)

#convert table to dataframe format
LGR_data<- as.data.frame(LGR_data) 

#confirm data types are correct
str(LGR_data)

# Exporting LGR data -----------------------------------------------------------

###Write Files###
#Optionally writes compiled data as TXT or CSV file
#write.table(LGR_data,file=paste(output_location, "/LGR_data_", Sys.time(), ".txt", sep=""))   
write.csv(LGR_data,file=paste(output_location, "/LGR_data_", Sys.time(), ".csv", sep=""))


# Load field notes -----------------------------------------------------------


# Create new variables for field data
#This must include a start and end time
#there must also be sufficient information to create a unique ID for each flux observation (see below)
#the code in sections below will use those unique flux ID keys in storing/calculating each flux

setwd("/Users/jongewirtzman/Google Drive/YMF Methane/Data")                

#note "LGR_fn" here means "field notes" not "function"... not the best naming, my apologies

LGR_fn <- read.table("Gewirtzman_datasheets_fixing_dates.csv", fill=TRUE,header=TRUE,sep=",",na.strings=c("NA","#N/A!","#N/A","#NA!"))    


#These lines fix the date formats, as read from the input field notes
#depending on the formats of dates in your field notes these lines may need to be altered

LGR_fn$Date_full<-paste("20", LGR_fn$Date, sep="")
LGR_fn$date.parsed<-ymd(LGR_fn$Date_full)

LGR_fn$starttime<-ymd_hms(paste(LGR_fn$date.parsed, LGR_fn$HMS_Start))

LGR_fn$endtime<-ymd_hms(paste(LGR_fn$date.parsed, LGR_fn$HMS_End))


#These lines defines a unique flux ID from my YMF field notes
#you will need to come up with a system for assigning unique IDs based on your specific dataset

LGR_fn$fluxid<-paste(LGR_fn$Tree_ID, LGR_fn$Measurement.height, sep="_")


#Subset a new dataframe with just the flux IDs and start/end times (may need to be modified to fit your field data columns)
#here you will again need to select the correct columns based on your field note format-
#in this case, it is columns 21-23 in my field notes

LGR_fn.labels<-LGR_fn[,21:23]
lgrfn<-LGR_fn.labels
names(lgrfn)<-c('start', 'end', 'fluxid')

#Subset LGR data by removing duplicated columns
#LGR_data_key<-LGR_data[,c(1:34, 36)]
LGR_data_key <- LGR_data[, !duplicated(colnames(LGR_data))]

#Check structure and remaining columns in LGR data
str(LGR_data)


# Select observations for each flux measurement -----------------------------------------------------------

#these lines assign the unique flux ID to all rows of data that belong to a particular flux
#i.e., if a flux "Tree_1A" starts at 0:01:00 and goes to 0:03:00, it will add a column noting "Tree_1A" to all rows of LGR data between those two times
#then, in following sections, we can easily select rows of data based on this key ID

#Use SQL package to merge flux ID (from LGR.fn) with LGR data (from LGR files)
#Appends the ID with all measurement points falling between the flux's  associated start and end time
result = sqldf("select * from LGR_fn
                left join LGR_data_key
                on LGR_data_key.date_time between LGR_fn.starttime and LGR_fn.endtime")


#check structure of result
str(result)

#make new dataframe where all data is ordered chronologically
#for some reason I chose the horrendously uninformative name "result1" (shame on you, past self)
#note that "result1" is the dataframe used for flux calclations below
result1<-result[order(result$date_time), ]

#make new list of all unique fliuxes to be calculated
fluxes<-unique(result1$fluxid)


# QA Plots  -----------------------------------------------------------
#This makes a plot of each flux, by plotting all points from start time to end time
#The loop will break if there are any issues where start times/end times are wrong-- e.g. if a decimal is wrong and there's a negative time length etc.
#future versions of the code could definitely be altered to allow the loop to keep running and just flag the issues... 

#This should be used to check that each subset looks like a linear flux and has correct start and end times
#If manual changes to start and end time are needed, change them in the field notes file, and run script again up to here
#Can move on past this step when all plots look like linear (or expected/correct duration) fluxes
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


# Flux calculation  -----------------------------------------------------------


#Convert time to minutes
#time in min
result1$tmin<-(result1$date_time-result1$starttime)/60

#Convert name to be compatible with other functions
#convert names
result1$ch4dppm<-result1$`[CH4]d_ppm`

#remove unneeded columns#optional
#result1a<-result1[,c('fluxid', 'ch4dppm', 'tmin', 'vol', 'Temp_air', 'area')]
result1a<-result1

#fix formatting if needed
result1a$tmin<-as.numeric(result1a$tmin)
rownames(result1a)<-NULL

#save results as CSV if needed
#write.csv(result1, "result1.csv")
#result1a<-read.csv("result1_new.csv")


#remove lines with issues with timestamps
#this is a temporary fix-- you should have all times correct in the code by now...
sum(is.na(result1a$tmin))
result1a<-result1a[which(is.na(result1a$tmin)==F),]
result1a<-result1a[which(is.na(result1a$tmin)==F),]

#make list of all fluxes that are still in the dataframe
fluxes<-fluxes[fluxes %in% result1a$fluxid]

#convert time to seconds
result1a$tsec<-result1a$tmin*60


###CALCULATE FLUXES###


##First we calculate slopes for each measurement; these are used below to calculate fluxes
#make vectors to hold calculated outputs
ch4_slopes<-c()
co2_slopes<-c()
ch4_r2<-c()
co2_r2<-c()

#Run loop that calculates slope and r-squared value for each flux
#This is a simple linear regression for each flux; y=dry-corrected ch4 ppm; x=time (seconds)
#Eventually this should also calculate non-linear (exponential) fluxes but I haven't coded that yet
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

#make a dataframe from the calculated outputs
outputs<-cbind(fluxes, ch4_slopes, co2_slopes, ch4_r2, co2_r2)
outputs<-as.data.frame(outputs)
names(outputs)[1]<-"fluxid"

output2<-left_join(outputs, LGR_fn, by="fluxid")
#summary(lm(output2$ch4_slopes~output2$Measurement.height+output2$Tree_sp))

#some manual corrections for my field notes... ignore these in your own version of the code
#just leaving them here so you can demonstrate running the code on my own input files
str(output2)
output2$ch4_slopes<-as.numeric(output2$ch4_slopes)
output2$Measurement.height<-as.numeric(output2$Measurement.height)
output2$Tree_sp<-as.factor(output2$Tree_sp)
output2<-output2[which(is.na(output2$Measurement.height)==F),]
output2$Measurement.height[which(output2$Measurement.height==75)]<-50


######Calcualte fluxes from slopes######

#this is relevant if different chambers are being used
#if a single chamber was used for all fluxes, single area and volume parameters can be added here

#make vectors to store chamber area and volume for each flux
chamber_area<-c()
chamber_vol<-c()

#check all chamber IDs
unique(output2$Chamber.ID)

#manual correction for my field notes-- ignore
output2$Chamber.ID[which(is.na(output2$Chamber.ID)==T)]<-6


##Enter the volume for each chamber

#Here, chambers 1/2, 3/4, 5/6 are the same size as each other
#Volume is in units of m3
#Chamber volumes must have been previously calculated empirically to enter here
chamber_vol[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-0.3*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002
chamber_vol[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-0.75*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002
chamber_vol[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-2.3*0.001+ pi*0.00175^2*(1.62+1.70) + 0.0002

#Here, enter surface area enclosed by each chamber
#Same caveats apply as for volume
#in m2
chamber_area[which(output2$Chamber.ID==1|output2$Chamber.ID==2)]<-14.16*9.69*0.0001
chamber_area[which(output2$Chamber.ID==3|output2$Chamber.ID==4)]<-20.04*13.6*0.0001
chamber_area[which(output2$Chamber.ID==5|output2$Chamber.ID==6)]<-28.07*19.3*0.0001 


#Enter atmospheric pressure
#Can either use standard atmospheric pressure, or empircally measured pressure if available
#If so, should be a column in field notes rather than a variable here
AtmPress<-c(rep(101.325, nrow(output2))) # If Atmospheric Pressure is not available at the measuring time, we can put a constant value (kg m-2 s-1)
#pressure in kPa

#Enter gas law constant
glc<-c(rep(0.00831447, nrow(output2))) # ideal gas law constant in kg m2 micromol-1 k-1


#append columns to dataframe
output2$chamber_vol<-chamber_vol
output2$chamber_area<-chamber_area
output2$glc<-glc
output2$pressure<-AtmPress


#Calculate fluxes corrected for volume, area, and pressure

#CH4 flux
output2$CH4_lin_flux<-
  output2$ch4_slopes*(output2$chamber_vol/output2$chamber_area)*
  (output2$pressure/(output2$glc*(output2$Temp_air+273)))*1000 #nmol m-2 s-1

#CO2 flux
output2$CO2_lin_flux<-
  as.numeric(output2$co2_slopes)*(output2$chamber_vol/output2$chamber_area)*
  (output2$pressure/(output2$glc*(output2$Temp_air+273))) #micromol m-2 s-1

###