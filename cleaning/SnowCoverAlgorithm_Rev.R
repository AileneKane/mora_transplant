##############################################
##############################################

# TITLE: Mean Temprature algorithm


# CONTACT: Ailene Ettinger 


# DESCRIPTION:  The purpose of this code is to summarize temperature data from many temperature sensors, to yield Growing Degree Days (GDD), with a 5C threshold, from the soil surface.

# To acquire GDD, the sensor must not be covered in snow. For the algorithm to consider a soil temperature sensor covered by snow, two criteria must be met:

# 1) the diurnal ground temperature range on that day does not exceed a certain value 
# (Recommended: 1 deg C).

# 2) the temperature at that time step does not exceed a maximum value (Recommended: in theory, 
# 0 deg C, but some records flatline as high as 2 deg C, likely due to bias in the sensor).

# In Mark's original algorithm, there was a third criterion: the temperature is constant for a user-
# specified number of timesteps. For Mark's original algorithm to score a sensor as "snow 
# covered" on a given date, two of the three criteria must be met. In Kevin's modified algorithm 
# (this script), the algorithm only considers the first two criteria and both of those criteria must 
# be met for the algorithm to score a sensor as "snow covered" on a given date. Kevin did not use 
# the third criterion because his dataset contained temperature sensors with different data logging 
# intervals.

# You can try the recommended vaules for the algorithm criteria, but you should always check the 
# output to make sure these parameters were appropriate. The idiosyncrasies of any one study require
# the algorithm to be modified to each study. This script plots soil temperature and snow algorithm
# output on the same graph for each sensor to help you do these checks.

# This script is an example of how to use the snow cover algorithm to estimate snow cover on a 
# daily basis for two example temperature sensors (ibuttons deployed at Mount Rainier National Park). 
# You can modify the script to analyze your own data. However, the .csv data files containing your soil 
# temperature measurements must have the exact same format as .csv files I provide as an example.


# REQUIRED DATA FILES:

##############################################
##############################################
setwd("~/Dropbox/Documents/Work/uw/Research/Mount Rainier/2014ManuscriptTransplant/microclimate")

#system.time({ #start the clock to record how long the code takes to run

## READ IN THE DATA FILES FOR EACH TEMPERATURE SENSOR

### !!! LIST YOUR FILES HERE !!! ###
path<-"/Users/aileneettinger/Dropbox/Documents/Work/uw/Research/Mount Rainier/2014ManuscriptTransplant/microclimate/CombinedFiles_2010_to_2013"
files <- c("NISQ-GC-1ALL.csv","NISQ-GC-2ALL.csv","NISQ-GC-3ALL.csv","NISQ-GC-4ALL.csv","NISQ-GC-5ALL.csv","NISQ-GX-1ALL.csv","NISQ-GX-2ALL.csv","NISQ-GX-3ALL.csv","NISQ-GX-4ALL.csv","NISQ-GX-5ALL.csv","NISQ-NC-1ALL.csv","NISQ-NC-2ALL.csv","NISQ-NC-3ALL.csv","NISQ-NC-4ALL.csv","NISQ-NC-5ALL.csv","NISQ-NX-1ALL.csv","NISQ-NX-2ALL.csv","NISQ-NX-3ALL.csv","NISQ-NX-4ALL.csv","NISQ-NX-5ALL.csv","HIGH-GC-1ALL.csv","HIGH-GC-2ALL.csv","HIGH-GC-3ALL.csv","HIGH-GC-4ALL.csv","HIGH-GC-5ALL.csv","HIGH-GX-1ALL.csv","HIGH-GX-2ALL.csv","HIGH-GX-3ALL.csv","HIGH-GX-4ALL.csv","HIGH-GX-5ALL.csv","HIGH-NC-1ALL.csv","HIGH-NC-2ALL.csv","HIGH-NC-3ALL.csv","HIGH-NC-4ALL.csv","HIGH-NC-5ALL.csv","HIGH-NX-1ALL.csv","HIGH-NX-2ALL.csv","HIGH-NX-3ALL.csv","HIGH-NX-4ALL.csv","HIGH-NX-5ALL.csv","PARA-GC-1ALL.csv","PARA-GC-2ALL.csv","PARA-GC-3ALL.csv","PARA-GC-4ALL.csv","PARA-GC-5ALL.csv","PARA-GX-1ALL.csv","PARA-GX-2ALL.csv","PARA-GX-3ALL.csv","PARA-GX-4ALL.csv","PARA-GX-5ALL.csv","PARA-NC-1ALL.csv","PARA-NC-2ALL.csv","PARA-NC-3ALL.csv","PARA-NC-4ALL.csv","PARA-NC-5ALL.csv","PARA-NX-1ALL.csv","PARA-NX-2ALL.csv","PARA-NX-3ALL.csv","PARA-NX-4ALL.csv","PARA-NX-5ALL.csv","AE10-GC-1ALL.csv","AE10-GC-2ALL.csv","AE10-GC-3ALL.csv","AE10-GC-4ALL.csv","AE10-GC-5ALL.csv","AE10-GX-1ALL.csv","AE10-GX-2ALL.csv","AE10-GX-3ALL.csv","AE10-GX-4ALL.csv","AE10-GX-5ALL.csv","AE10-NC-1ALL.csv","AE10-NC-2ALL.csv","AE10-NC-3ALL.csv","AE10-NC-4ALL.csv","AE10-NC-5ALL.csv","AE10-NX-1ALL.csv","AE10-NX-2ALL.csv","AE10-NX-3ALL.csv","AE10-NX-4ALL.csv","AE10-NX-5ALL.csv","AM16-GC-1ALL.csv","AM16-GC-2ALL.csv","AM16-GC-3ALL.csv","AM16-GC-4ALL.csv","AM16-GC-5ALL.csv","AM16-GX-1ALL.csv","AM16-GX-2ALL.csv","AM16-GX-3ALL.csv","AM16-GX-4ALL.csv","AM16-GX-5ALL.csv","AM16-NC-1ALL.csv","AM16-NC-2ALL.csv","AM16-NC-3ALL.csv","AM16-NC-4ALL.csv","AM16-NC-5ALL.csv","AM16-NX-1ALL.csv","AM16-NX-2ALL.csv","AM16-NX-3ALL.csv","AM16-NX-4ALL.csv","AM16-NX-5ALL.csv","AV06-GC-1ALL.csv","AV06-GC-2ALL.csv","AV06-GC-3ALL.csv","AV06-GC-4ALL.csv","AV06-GC-5ALL.csv","AV06-GX-1ALL.csv","AV06-GX-2ALL.csv","AV06-GX-3ALL.csv","AV06-GX-4ALL.csv","AV06-GX-5ALL.csv","AV06-NC-1ALL.csv","AV06-NC-2ALL.csv","AV06-NC-3ALL.csv","AV06-NC-4ALL.csv","AV06-NC-5ALL.csv","AV06-NX-1ALL.csv","AV06-NX-2ALL.csv","AV06-NX-3ALL.csv","AV06-NX-4ALL.csv","AV06-NX-5ALL.csv","TO04-GC-1ALL.csv","TO04-GC-2ALL.csv","TO04-GC-3ALL.csv","TO04-GC-4ALL.csv","TO04-GC-5ALL.csv","TO04-GX-1ALL.csv","TO04-GX-2ALL.csv","TO04-GX-3ALL.csv","TO04-GX-4ALL.csv","TO04-GX-5ALL.csv","TO04-NC-1ALL.csv","TO04-NC-2ALL.csv","TO04-NC-3ALL.csv","TO04-NC-4ALL.csv","TO04-NC-5ALL.csv","TO04-NX-1ALL.csv","TO04-NX-2ALL.csv","TO04-NX-3ALL.csv","TO04-NX-4ALL.csv","TO04-NX-5ALL.csv") # the data files for each temperature sensor to be analyzed 
#files <- c("NISQ-GC-2ALL.csv")
#"NISQ-NC-5ALL.csv","HIGH-GX-5ALL.csv", "HIGH-NC-1ALL.csv",and "HIGH-GX-1ALL.csv"removed for now because it is missing temperature data for many dates
#NISQ files,"NISQ-GC-1ALL.csv","NISQ-GC-2ALL.csv","NISQ-GC-3ALL.csv","NISQ-GC-4ALL.csv","NISQ-GC-5ALL.csv","NISQ-GX-1ALL.csv","NISQ-GX-2ALL.csv","NISQ-GX-3ALL.csv","NISQ-GX-4ALL.csv","NISQ-GX-5ALL.csv","NISQ-NC-1ALL.csv","NISQ-NC-2ALL.csv","NISQ-NC-3ALL.csv","NISQ-NC-4ALL.csv","NISQ-NX-1ALL.csv","NISQ-NX-2ALL.csv","NISQ-NX-3ALL.csv","NISQ-NX-4ALL.csv","NISQ-NX-5ALL.csv",
#HIGH files: "HIGH-GC-1ALL.csv","HIGH-GC-2ALL.csv","HIGH-GC-3ALL.csv","HIGH-GC-4ALL.csv","HIGH-GC-5ALL.csv","HIGH-GX-2ALL.csv","HIGH-GX-3ALL.csv","HIGH-GX-4ALL.csv","HIGH-NC-2ALL.csv","HIGH-NC-3ALL.csv","HIGH-NC-4ALL.csv","HIGH-NC-5ALL.csv","HIGH-NX-1ALL.csv","HIGH-NX-2ALL.csv","HIGH-NX-3ALL.csv","HIGH-NX-4ALL.csv","HIGH-NX-5ALL.csv","PARA-GC-1ALL.csv","PARA-GC-2ALL.csv","PARA-GC-3ALL.csv","PARA-GC-4ALL.csv","PARA-GC-5ALL.csv","PARA-GX-1ALL.csv","PARA-GX-2ALL.csv","PARA-GX-3ALL.csv","PARA-GX-4ALL.csv","PARA-GX-5ALL.csv","PARA-NC-1ALL.csv","PARA-NC-2ALL.csv","PARA-NC-3ALL.csv","PARA-NC-4ALL.csv","PARA-NC-5ALL.csv","PARA-NX-1ALL.csv","PARA-NX-2ALL.csv","PARA-NX-3ALL.csv","PARA-NX-4ALL.csv","PARA-NX-5ALL.csv","AE10-GC-1ALL.csv","AE10-GC-2ALL.csv","AE10-GC-3ALL.csv","AE10-GC-4ALL.csv","AE10-GC-5ALL.csv","AE10-GX-1ALL.csv","AE10-GX-2ALL.csv","AE10-GX-3ALL.csv","AE10-GX-4ALL.csv","AE10-GX-5ALL.csv","AE10-NC-1ALL.csv","AE10-NC-2ALL.csv","AE10-NC-3ALL.csv","AE10-NC-4ALL.csv","AE10-NC-5ALL.csv","AE10-NX-1ALL.csv","AE10-NX-2ALL.csv","AE10-NX-3ALL.csv","AE10-NX-4ALL.csv","AE10-NX-5ALL.csv","AM16-GC-1ALL.csv","AM16-GC-2ALL.csv","AM16-GC-3ALL.csv","AM16-GC-4ALL.csv","AM16-GC-5ALL.csv","AM16-GX-1ALL.csv","AM16-GX-2ALL.csv","AM16-GX-3ALL.csv","AM16-GX-4ALL.csv","AM16-GX-5ALL.csv","AM16-NC-1ALL.csv","AM16-NC-2ALL.csv","AM16-NC-3ALL.csv","AM16-NC-4ALL.csv","AM16-NC-5ALL.csv","AM16-NX-1ALL.csv","AM16-NX-2ALL.csv","AM16-NX-3ALL.csv","AM16-NX-4ALL.csv","AM16-NX-5ALL.csv","AV06-GC-1ALL.csv","AV06-GC-2ALL.csv","AV06-GC-3ALL.csv","AV06-GC-4ALL.csv","AV06-GC-5ALL.csv","AV06-GX-1ALL.csv","AV06-GX-2ALL.csv","AV06-GX-3ALL.csv","AV06-GX-4ALL.csv","AV06-GX-5ALL.csv","AV06-NC-1ALL.csv","AV06-NC-2ALL.csv","AV06-NC-3ALL.csv","AV06-NC-4ALL.csv","AV06-NC-5ALL.csv","AV06-NX-1ALL.csv","AV06-NX-2ALL.csv","AV06-NX-3ALL.csv","AV06-NX-4ALL.csv","AV06-NX-5ALL.csv","TO04-GC-1ALL.csv","TO04-GC-2ALL.csv","TO04-GC-3ALL.csv","TO04-GC-4ALL.csv","TO04-GC-5ALL.csv","TO04-GX-1ALL.csv","TO04-GX-2ALL.csv","TO04-GX-3ALL.csv","TO04-GX-4ALL.csv","TO04-GX-5ALL.csv","TO04-NC-1ALL.csv","TO04-NC-2ALL.csv",

#,,,
snow_appearance_dateyr1 <- c() 
snow_disappearance_dateyr1 <- c()
snow_cover_durationyr1 <- c() # in days, counts days with snowcover
snow_cover_durationyr1_cont<-c()#in days, assumes continual coverage from first appearance date to last appearance date
snow_appearance_dateyr2 <- c() 
snow_disappearance_dateyr2 <- c()
snow_cover_durationyr2 <- c() # in days
snow_cover_durationyr2_cont<-c()#in days, assumes continual coverage from first appearance date to last appearance date

year1<-c()
year2<-c()
startdateyr1<-c()
enddateyr1<-c()
startdateyr2<-c()
enddateyr2<-c()
cumna.tempyr1<-c()
cumna.tempyr2<-c()
cumnumrowsyr1<-c()
cumnumrowsyr2<-c()
for(k in 1:length(files)){
  ##READING IN DATA FROM ONE FILE###############################################
    file <- file.path(path, files[k])
    d <- read.csv(file)
    d$Temp_C<-as.numeric(d$Temp_C)
    d$Light_Lux<-as.numeric(d$Light_Lux)
    d$date<-as.Date(d$DateTime, format = "%m/%d/%Y") 
    if(length(which(d$Temp_C>100))>1){d$Temp_C[which(d$Temp_C>100)]<-NA}#remove extreme temperature values

    Date <- as.factor(as.Date(d$date)) #pulls out the date data for each row
    dates.unique <- levels(Date) #all the unique dates with temperature data
    
    ###############################################
    
    ##CRITERION 1
    
  RangeThresh <- 1 #the threshold temperature range in degrees C (i.e. if the range of soil temperatures on a given day exceed RangeThresh, snow probably absent)
    # Calculate the daily temperature range for each day: 
  ##CRITERION 2
  MaxThresh <- 2 # the threshold maximum temperature (i.e. if the max temperature on a given day exceeded MaxThresh, snow probably absent)
  temp_min<-aggregate(x=subset(d, select=c("date","Temp_C")), by=list(d$date), FUN=min,na.rm=T)#daily min temperature
  temp_max<-aggregate(x=subset(d, select=c("date","Temp_C")), by=list(d$date), FUN=max,na.rm=T)
  temp_mean<-aggregate(x=subset(d, select=c("date","Temp_C")), by=list(d$date), FUN=mean,na.rm=T)
  na_count <-function(y) {length(which(is.na(y)))} 
  temp_nas<-aggregate(x=subset(d, select=c("date","Temp_C")), by=list(d$date), FUN=na_count)
  numrows<-aggregate(d$Temp_C, by=list(d$date), FUN=length)
  allmicro<-cbind(temp_min,temp_max$Temp_C,temp_mean$Temp_C,temp_nas$Temp_C,numrows[,2])
  allmicro<-allmicro[,-1]
  colnames(allmicro)<-c("date","mint","maxt","meant","numnast","numrows")
  allmicro$temp_range<-allmicro$maxt-allmicro$mint
  allmicro$range_logical <- NA# will store the binary code for daily temperature range for each day
  allmicro[which(allmicro$temp_range <= RangeThresh),]$range_logical<-1
  allmicro[which(allmicro$temp_range > RangeThresh),]$range_logical<-0
# Convert the daily temperature mean into a binary code: 
# 1 = daily temperature max did not exceed threshold (MaxThresh) AND range was not greater than 1 degree celsius, which equals 2 for snow
# 0 = daily temperature mean did  exceed threshold (MaxThresh)
allmicro$snow_logical<-0
allmicro[which(allmicro$meant<=MaxThresh & allmicro$range_logical==1),]$snow_logical<-1
allmicro$year<-paste("20",substr(as.Date(allmicro$date),3,4), sep="")

#summarize to annual values
yr1startdate<-"0010-11-01"
yr1enddate<-"0011-10-31"
yr1jan1date<-"0011-01-01"
yr2startdate<-"0011-11-01"
yr2enddate<-"0012-09-15"
yr2jan1date<-"0012-01-01"
yr1gstdate1<-"0011-03-01"# March 1 
yr2gstdate1<-"0012-03-01"# March 1 
yr1gstdate2<-"0011-08-31"# August 31
yr2gstdate2<-"0012-08-31"# August 31
allmicro$styear<-1#assign a study year (1 or 2)
allmicro[which(allmicro$date>=yr2startdate),]$styear<-2
yr1firstrow<-which(allmicro$date==yr1startdate)#start with Nov 1, 2010
yr1lastrow<-which(allmicro$date==yr1enddate)#end with Oct 31 ,2011
yr2firstrow<-which(allmicro$date==yr2startdate)#start with Nov 1, 2011
yr2lastrow<-max(which(allmicro$date<=yr2enddate))#end with sept 15 ,2012 or earlier- whatever the max date is before this (data only goes this far)
firstgddrow<-which(allmicro$date=="0011-01-01")#for growing degree days 
year1[k]<-allmicro[which(allmicro$date==yr1enddate),]$year
year2[k]<-allmicro[yr2lastrow,]$year
startdateyr1[k]<-as.character(allmicro[yr1firstrow,]$date)
enddateyr1[k]<-as.character(allmicro[yr1lastrow,]$date)
startdateyr2[k]<-as.character(allmicro[yr2firstrow,]$date)
enddateyr2[k]<-as.character(allmicro[yr2lastrow,]$date)
cumna.tempyr1[k]<-sum(allmicro[which(allmicro$date>=yr1startdate & allmicro$date<=yr1enddate),]$numnast)#n
cumna.tempyr2[k]<-sum(allmicro[which(allmicro$date>=yr2startdate & allmicro$date<=yr2enddate),]$numnast)#n
cumnumrowsyr1[k]<-sum(allmicro[which(allmicro$date>=yr1startdate & allmicro$date<=yr1enddate),]$numrows)#number of total rows (i.e. measurement points/samples) in year 2- should be ~8-12 per day
cumnumrowsyr2[k]<-sum(allmicro[which(allmicro$date>=yr2startdate & allmicro$date<=yr2enddate),]$numrows)#number of total rows (i.e. measurement points/samples) in year 2- should be ~8-12 per day
snow_presentyr1 <- allmicro[which(allmicro$snow_logical==1 & allmicro$date>=yr1startdate & allmicro$date<=yr1enddate),]$date # the dates where snow is present
snow_presentyr2 <- allmicro[which(allmicro$snow_logical==1 & allmicro$date>=yr2startdate & allmicro$date<=yr2enddate),]$date # the dates where snow is present

snow_appearance_dateyr1[k] <- as.character(snow_presentyr1[1]) # first day when snow covered sensor
if(length(snow_presentyr1)>0){
  snow_disappearance_dateyr1[k] <- as.character(snow_presentyr1[length(snow_presentyr1)])
  }
else{snow_disappearance_dateyr1[k] <-NA}# last day when snow covered sensor
snow_cover_durationyr1[k] <- sum(allmicro[which(allmicro$date>=yr1startdate & allmicro$date<=yr1enddate),]$snow_logical) #"snow cover duration" the total number of days with snow cover
if(!is.na(snow_appearance_dateyr1[k])){
  snow_cover_durationyr1_cont[k] <- as.Date(snow_disappearance_dateyr1[k])-as.Date(snow_appearance_dateyr1[k]) #"snow cover duration" the total number of days with snow cover
} else{
  snow_cover_durationyr1_cont[k] <-NA
}
  snow_appearance_dateyr2[k] <- as.character(snow_presentyr2[1]) # first day when snow covered sensor
snow_disappearance_dateyr2[k] <- as.character(snow_presentyr2[length(snow_presentyr2)]) # last day when snow covered sensor
snow_cover_durationyr2[k] <- sum(allmicro[which(allmicro$date>=yr2startdate & allmicro$date<=yr2enddate),]$snow_logical) #"snow cover duration" the total number of days with snow cover
if(!is.na(snow_appearance_dateyr2[k])){
  snow_cover_durationyr2_cont[k] <- as.Date(snow_disappearance_dateyr2[k])-as.Date(snow_appearance_dateyr2[k]) #"snow cover duration" the total number of days with snow cover
} else{
  snow_cover_durationyr2_cont[k] <-NA
}
}
################################################################################################
output1 <- data.frame(files,year1,snow_appearance_dateyr1,snow_disappearance_dateyr1, snow_cover_durationyr1,snow_cover_durationyr1_cont,startdateyr1,enddateyr1,cumna.tempyr1,cumnumrowsyr1)
output2 <- data.frame(files,year2,snow_appearance_dateyr2,snow_disappearance_dateyr2, snow_cover_durationyr2,snow_cover_durationyr2_cont,startdateyr2,enddateyr2,cumna.tempyr2,cumnumrowsyr2)


colnames(output2)<-colnames(output1)
output<-rbind(output1,output2)
colnames(output)<-c("filesname","year","snowapdate","snowdidate","snow_dur","snow_dur_cont","yrstartdate","yrenddate","cumna.temp","cumnumrows" )
write.csv(output,file="snowmicroclim20102013.csv", row.names=FALSE)


datesTest <- as.Date(allmicro$date)
Date<-allmicro$date
DateTmp <- as.Date(Date)
DateClip <- c()
TempClip <- c()
for (i in 1:length(DateTmp)) {
  if (DateTmp[i] >= min(datesTest) & DateTmp[i] <= max(datesTest) )
  { 
    DateClip <- c(DateClip, as.character(DateTmp[i])) 
    TempClip <- c(TempClip, Temp[i])
  }
}





date.plot <- as.Date(DateClip)

quartz(15,9)

#Left Y axis
par(mar=c(4,4,3,4))
plot(date.plot, TempClip, main=files[k], axes=F, xlab="", ylab="", pch=20, col="red")
axis(2, col="red", col.axis="red", col.ticks="red")
leftY <- paste("Temperature (?C)")
text(par("usr")[1] - 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj=.5, 
     leftY, srt=90, xpd=TRUE, col="red")

#Right Y axis
par(new=TRUE)
plot(allmicro$date, allmicro$GDD_logical, pch=20, cex=1.5,col="blue", ylim=c(0,1.2), axes=FALSE, xaxt="n",yaxt="n",xlab="",ylab="")
axis(4, col="blue", col.axis="blue", col.ticks="blue", yaxp=c(0,1,1), labels=F)
rightY <- paste ("GDD")
text(par("usr")[2] + 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj=.5,
     rightY, srt = 270, xpd = TRUE, col="blue")
text(par("usr")[2] + 15, 0, 0, srt=270, xpd=T, col="blue")
text(par("usr")[2] + 15, 1, 1, srt=270, xpd=T, col="blue")

#X axis
axis.Date(1, date.plot, 
          at=seq(from=min(date.plot), to=max(date.plot), by="months"))

}
################################################################################################
output1 <- data.frame(files,year1,GDD_totalyr1,GDD_totalnyr1, Light_Meanyr1,Light_GDDyr1,GST_Meanyr1,MATyr1)
output2 <- data.frame(files,year2,GDD_totalyr2,GDD_totalnyr2, Light_Meanyr2,Light_GDDyr2,GST_Meanyr2,MATyr2)
colnames(output2)<-colnames(output1)
output<-rbind(output1,output2)
write.csv(output,file="ALLmicroclim20102013.csv")
#plot some correlations
plot(output$Light_GDDyr1,output$Light_Meanyr1, xlim=c(0,10000),ylim=c(0,5000))
## PLOT SOIL TEMPERATURE AND THE GDD ALGORITH OUTPUT TO MAKE SURE OUTPUT IS REASONABLE

datesTest <- as.Date(allmicro$date)
DateTmp <- as.Date(Date)
DateClip <- c()
TempClip <- c()
for (i in 1:length(DateTmp)) {
	if (DateTmp[i] >= min(datesTest) & DateTmp[i] <= max(datesTest) )
		{ 
		DateClip <- c(DateClip, as.character(DateTmp[i])) 
		TempClip <- c(TempClip, Temp[i])
		}
}



date.plot <- as.Date(DateClip)

quartz(15,9)

#Left Y axis
par(mar=c(4,4,3,4))
plot(date.plot, TempClip, main=files[k], axes=F, xlab="", ylab="", pch=20, col="red")
axis(2, col="red", col.axis="red", col.ticks="red")
leftY <- paste("Temperature (?C)")
text(par("usr")[1] - 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj=.5, 
	leftY, srt=90, xpd=TRUE, col="red")

#Right Y axis
par(new=TRUE)
plot(allmicro$date, allmicro$GDD_logical, pch=20, cex=1.5,col="blue", ylim=c(0,1.2), axes=FALSE, xaxt="n",yaxt="n",xlab="",ylab="")
axis(4, col="blue", col.axis="blue", col.ticks="blue", yaxp=c(0,1,1), labels=F)
rightY <- paste ("GDD")
text(par("usr")[2] + 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj=.5,
	rightY, srt = 270, xpd = TRUE, col="blue")
text(par("usr")[2] + 15, 0, 0, srt=270, xpd=T, col="blue")
text(par("usr")[2] + 15, 1, 1, srt=270, xpd=T, col="blue")

#X axis
axis.Date(1, date.plot, 
at=seq(from=min(date.plot), to=max(date.plot), by="months"))


#ALL STands:"NISQ-GC-1ALL.csv","NISQ-GC-2ALL.csv","NISQ-GC-3ALL.csv","NISQ-GC-4ALL.csv","NISQ-GC-5ALL.csv","NISQ-GX-1ALL.csv","NISQ-GX-2ALL.csv","NISQ-GX-3ALL.csv","NISQ-GX-4ALL.csv","NISQ-GX-5ALL.csv","NISQ-NC-1ALL.csv","NISQ-NC-2ALL.csv","NISQ-NC-3ALL.csv","NISQ-NC-4ALL.csv","NISQ-NC-5ALL.csv","NISQ-NX-1ALL.csv","NISQ-NX-2ALL.csv","NISQ-NX-3ALL.csv","NISQ-NX-4ALL.csv","NISQ-NX-5ALL.csv","HIGH-GC-1ALL.csv","HIGH-GC-2ALL.csv","HIGH-GC-3ALL.csv","HIGH-GC-4ALL.csv","HIGH-GC-5ALL.csv","HIGH-GX-1ALL.csv","HIGH-GX-2ALL.csv","HIGH-GX-3ALL.csv","HIGH-GX-4ALL.csv","HIGH-GX-5ALL.csv","HIGH-NC-1ALL.csv","HIGH-NC-2ALL.csv","HIGH-NC-3ALL.csv","HIGH-NC-4ALL.csv","HIGH-NC-5ALL.csv","HIGH-NX-1ALL.csv","HIGH-NX-2ALL.csv","HIGH-NX-3ALL.csv","HIGH-NX-4ALL.csv","HIGH-NX-5ALL.csv","PARA-GC-1ALL.csv","PARA-GC-2ALL.csv","PARA-GC-3ALL.csv","PARA-GC-4ALL.csv","PARA-GC-5ALL.csv","PARA-GX-1ALL.csv","PARA-GX-2ALL.csv","PARA-GX-3ALL.csv","PARA-GX-4ALL.csv","PARA-GX-5ALL.csv","PARA-NC-1ALL.csv","PARA-NC-2ALL.csv","PARA-NC-3ALL.csv","PARA-NC-4ALL.csv","PARA-NC-5ALL.csv","PARA-NX-1ALL.csv","PARA-NX-2ALL.csv","PARA-NX-3ALL.csv","PARA-NX-4ALL.csv","PARA-NX-5ALL.csv","AE10-GC-1ALL.csv","AE10-GC-2ALL.csv","AE10-GC-3ALL.csv","AE10-GC-4ALL.csv","AE10-GC-5ALL.csv","AE10-GX-1ALL.csv","AE10-GX-2ALL.csv","AE10-GX-3ALL.csv","AE10-GX-4ALL.csv","AE10-GX-5ALL.csv","AE10-NC-1ALL.csv","AE10-NC-2ALL.csv","AE10-NC-3ALL.csv","AE10-NC-4ALL.csv","AE10-NC-5ALL.csv","AE10-NX-1ALL.csv","AE10-NX-2ALL.csv","AE10-NX-3ALL.csv","AE10-NX-4ALL.csv","AE10-NX-5ALL.csv","AM16-GC-1ALL.csv","AM16-GC-2ALL.csv","AM16-GC-3ALL.csv","AM16-GC-4ALL.csv","AM16-GC-5ALL.csv","AM16-GX-1ALL.csv","AM16-GX-2ALL.csv","AM16-GX-3ALL.csv","AM16-GX-4ALL.csv","AM16-GX-5ALL.csv","AM16-NC-1ALL.csv","AM16-NC-2ALL.csv","AM16-NC-3ALL.csv","AM16-NC-4ALL.csv","AM16-NC-5ALL.csv","AM16-NX-1ALL.csv","AM16-NX-2ALL.csv","AM16-NX-3ALL.csv","AM16-NX-4ALL.csv","AM16-NX-5ALL.csv","AV06-GC-1ALL.csv","AV06-GC-2ALL.csv","AV06-GC-3ALL.csv","AV06-GC-4ALL.csv","AV06-GC-5ALL.csv","AV06-GX-1ALL.csv","AV06-GX-2ALL.csv","AV06-GX-3ALL.csv","AV06-GX-4ALL.csv","AV06-GX-5ALL.csv","AV06-NC-1ALL.csv","AV06-NC-2ALL.csv","AV06-NC-3ALL.csv","AV06-NC-4ALL.csv","AV06-NC-5ALL.csv","AV06-NX-1ALL.csv","AV06-NX-2ALL.csv","AV06-NX-3ALL.csv","AV06-NX-4ALL.csv","AV06-NX-5ALL.csv","TO04-GC-1ALL.csv","TO04-GC-2ALL.csv","TO04-GC-3ALL.csv","TO04-GC-4ALL.csv","TO04-GC-5ALL.csv","TO04-GX-1ALL.csv","TO04-GX-2ALL.csv","TO04-GX-3ALL.csv","TO04-GX-4ALL.csv","TO04-GX-5ALL.csv","TO04-NC-1ALL.csv","TO04-NC-2ALL.csv","TO04-NC-3ALL.csv","TO04-NC-4ALL.csv","TO04-NC-5ALL.csv","TO04-NX-1ALL.csv","TO04-NX-2ALL.csv","TO04-NX-3ALL.csv","TO04-NX-4ALL.csv","TO04-NX-5ALL.csv"