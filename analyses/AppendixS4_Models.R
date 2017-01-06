#R Code for GCB Submission
#Competition and facilitation may lead to asymmetric range shift dynamics with climate change
#By Ailene Ettinger & Janneke HilleRisLambers
#Code for Germination, Survival and Growth Data
#Code written by Ailene Ettinger
#Completed January 6, 2016
setwd("~/git/mora_transplant")
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library(dplyr)
library(plyr)
library(survival)
library(car)
library(lme4)
library(bbmle)
library(boot)
library(raster)
library(AICcmodavg)
library(glmmADMB)
library(RColorBrewer)

#Begin with transplant data: get the survival times, growth variables, etc
transdat<-read.csv("data/2013TransplantStatusHeight(October).csv", header=TRUE)#csv file has been sorted so that all dates appear together, all status columns are together, etc
dim(transdat)#3959 rows(=individuals),  38 columns
transdat$PlantedStand2<-as.numeric(transdat$PlantedStand)
transdat$PlantedStand<-as.factor(transdat$PlantedStand)
transdat$OriginStand<-as.factor(transdat$OriginStand)
transdat$Block<-as.factor(transdat$Block)
transdat$UniqueID<-as.factor(transdat$UniqueID)
transdat$Date1<-as.Date(transdat$Date1,format='%m/%d/%y')
transdat$Date2<-as.Date(transdat$Date2,format='%m/%d/%y')
transdat$Date3<-as.Date(transdat$Date3,format='%m/%d/%y')
transdat$Date4<-as.Date(transdat$Date4,format='%m/%d/%y')
transdat$Date5<-as.Date(transdat$Date5,format='%m/%d/%y')
transdat[which(transdat$Date3=="2025-06-14"),]$Date3<-"0012-06-14"#fix typo
##add columns for relative growth rate and annual relative growth rate
transdat$rgr=NA
transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$rgr=(transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$Height2-transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$Initial.Height)/transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$Initial.Height#rgr for plants that only survived 1 year
transdat[which(is.na(transdat$HeightDate5)&transdat$HeightDate4>0),]$rgr=(transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$HeightDate4-transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$Initial.Height)/transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$Initial.Height#rgr for plants that survived 2 years
transdat[which(transdat$HeightDate5>0),]$rgr=(transdat[which(transdat$HeightDate5>0),]$HeightDate5-transdat[which(transdat$HeightDate5>0),]$Initial.Height)/transdat[which(transdat$HeightDate5>0),]$Initial.Height#rgr for plants that survived 3 years
#also try just height increment/yr
transdat$hi=NA
transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$hi=(transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$Height2-transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$Initial.Height)#height incr for plants that only survived 1 year
transdat[which(is.na(transdat$HeightDate5)&transdat$HeightDate4>0),]$hi=(transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$HeightDate4-transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$Initial.Height)#height incr. for plants that survived 2 years
transdat[which(transdat$HeightDate5>0),]$hi=(transdat[which(transdat$HeightDate5>0),]$HeightDate5-transdat[which(transdat$HeightDate5>0),]$Initial.Height)#height incr for plants that survived 3 years
#create column for annual rgr and annual hi
transdat$yrs=NA
transdat$annrgr=NA     
transdat$annhi=NA 
#first, identify # of yrs over which growth meaured
transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),]$yrs=rep(1,times=dim(transdat[which(is.na(transdat$HeightDate4)& transdat$Height2>0),])[1])
transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),]$yrs=rep(2,times=dim(transdat[which(is.na(transdat$HeightDate5)& transdat$HeightDate4>0),])[1])
transdat[which(transdat$HeightDate5>0),]$yrs=rep(3,times=dim(transdat[which(transdat$HeightDate5>0),])[1]) 
transdat$finalrgr=(transdat$RtCrnHeightDate5-transdat$Initial.Height)/transdat$Initial.Height
transdat$annrgr=transdat$rgr/transdat$yrs     
transdat$annhi=transdat$hi/transdat$yrs
#Add column for survival data  (0=alive, 1=dead)
transdat$Death=NA
transdat[which(transdat$StatusDate5==0),]$Death=1
transdat[which(transdat$StatusDate5==1),]$Death=0
#Add column for time to death in days
transdat$StartDate=as.Date("9/1/2010",format='%m/%d/%Y')
transdat$StartDateDate5=transdat$Date5-transdat$StartDate
transdat$StartDateDate4=transdat$Date4-transdat$StartDate
transdat$StartDateDate3=transdat$Date3-transdat$StartDate
transdat$StartDateDate2=transdat$Date2-transdat$StartDate
transdat$StartDateDate1=transdat$Date1-transdat$StartDate
DaysDeath<-c()
alltime1<-c()
alltime2<-c()
for (i in 1:dim(transdat)[1]){
  if(transdat$StatusDate5[i]==1)inddd<-transdat$StartDateDate5[i]
  if(transdat$StatusDate5[i]==0 & transdat$StatusDate4[i]==1) inddd<-transdat$StartDateDate5[i]
  if(transdat$StatusDate3[i]==1 & transdat$StatusDate4[i]==0) inddd<-transdat$StartDateDate4[i]
  if(transdat$StatusDate2[i]==1 & transdat$StatusDate3[i]==0) inddd<-transdat$StartDateDate3[i]
  if(transdat$StatusDate1[i]==1 & transdat$StatusDate2[i]==0) inddd<-transdat$StartDateDate2[i]
  if(transdat$StatusDate1[i]==0) inddd<-transdat$StartDateDate1[i]
  DaysDeath<-c(DaysDeath,inddd)
  #for interval censored data, need to specificy that we know it died between which intervals
  #If alive on final census, time2=NA, time1=#ofday on last census
  if(transdat$StatusDate5[i]==1)time1<-transdat$StartDateDate5[i]
  if(transdat$StatusDate5[i]==1)time2<-NA
  #If dead on final census, but alive on 4th census, time2=#ofdaydate5, time1=#ofdaydate4
  if(transdat$StatusDate5[i]==0 & transdat$StatusDate4[i]==1) time2<-transdat$StartDateDate5[i]
  if(transdat$StatusDate5[i]==0 & transdat$StatusDate4[i]==1) time1<-transdat$StartDateDate4[i]
  #If dead on 4th census, but alive on 3rd census, time2=#ofdaydate4, time1=#ofdaysondate3
  if(transdat$StatusDate3[i]==1 & transdat$StatusDate4[i]==0) time2<-transdat$StartDateDate4[i]
  if(transdat$StatusDate3[i]==1 & transdat$StatusDate4[i]==0) time1<-transdat$StartDateDate3[i]
  #If dead on 3rd census, but alive on 2nd census, time2=#ofdaydate3, time1=#ofdaysondate2
  if(transdat$StatusDate2[i]==1 & transdat$StatusDate3[i]==0) time2<-transdat$StartDateDate3[i]
  if(transdat$StatusDate2[i]==1 & transdat$StatusDate3[i]==0) time1<-transdat$StartDateDate2[i]
  #If dead on 2nd census, but alive on 1st census, time2=#ofdaydate2, time1=#ofdaysondate1
  if(transdat$StatusDate1[i]==1 & transdat$StatusDate2[i]==0) time2<-transdat$StartDateDate2[i]
  if(transdat$StatusDate1[i]==1 & transdat$StatusDate2[i]==0) time1<-transdat$StartDateDate1[i]
  #If dead on 1st census, time2=#ofdaydate1, time1=14 (all plants checked after 2 weeks and still alive)
  if(transdat$StatusDate1[i]==0) time2<-transdat$StartDateDate1[i]
  if(transdat$StatusDate1[i]==0) time1<-14
  alltime1<-c(alltime1,time1)
  alltime2<-c(alltime2,time2)
}
transdat$DaysDeath=DaysDeath
transdat$time1=as.numeric(alltime1)
transdat$time2=alltime2
tsmedat<-transdat[transdat$Species=="TSME",]
tshedat<-transdat[transdat$Species=="TSHE",]
abamdat<-transdat[transdat$Species=="ABAM",]
tsmedat$PlantedStand=factor(tsmedat$PlantedStand)
tshedat$PlantedStand=factor(tshedat$PlantedStand)
abamdat$PlantedStand=factor(abamdat$PlantedStand)
tsmedat$OriginStand=factor(tsmedat$OriginStand)
tshedat$OriginStand=factor(tshedat$OriginStand)
abamdat$OriginStand=factor(abamdat$OriginStand)
tsmedat$Block=factor(tsmedat$Block)
tshedat$Block=factor(tshedat$Block)
abamdat$Block=factor(abamdat$Block)
#First fit models for survival with all categorial data
#ABAM
constmod.abam<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=abamdat)
#TSME
constmod.tsme<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tsmedat)
#TSHE
constmod.tshe<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tshedat)
Anova(constmod.tsme, test.statistic="LR", type="III")
Anova(constmod.abam, test.statistic="LR", type="III")
Anova(constmod.tshe,test.statistic="LR", type="III")

#Next models for growth
consthimod.abam<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=abamdat)
#TSME
consthimod.tsme<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tsmedat)
#TSHE
consthimod.tshe<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tshedat)
Anova(consthimod.tsme, type="III")
Anova(consthimod.abam, type="III")
Anova(consthimod.tshe, type="III")

#Now add microclimate data, including GDD, snowduration, and mean light for each block
climdat<-read.csv("analyses/clim_gf.csv", header=TRUE)
#Combine the climate data with the transplant data columns that we want to use for survival and growth analysis
transdatsub<-subset(transdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi"))
climdat$Canopy<-"CompAbsent"
climdat[which(climdat$canopy=="N"),]$Canopy<-"CompPresent"
climdat$Understory<-"CompAbsent"
climdat[which(climdat$understory=="C"),]$Understory<-"CompPresent"
#take average climate across both years, for model comparison
climdatsub<-subset(climdat,select=c("stand","block","Canopy","Understory","GDD_total","GDD_totaln","GST_Mean","Light_Mean","Light_GDD","MAT","snow_dur","snow_dur_cont"))
climdatsub2<-aggregate(climdatsub, by=list(climdatsub$stand,climdatsub$block,climdatsub$Canopy,climdatsub$Understory), FUN=mean,na.rm=F)
climdatsub3<-subset(climdatsub2,select=c("Group.1","Group.2","Group.3","Group.4","GDD_total","GDD_totaln","GST_Mean","Light_Mean","Light_GDD","MAT","snow_dur","snow_dur_cont"))
colnames(climdatsub3)[1:4]<-c("Stand","Block","Canopy","Understory")
climdatsub3$Block<-as.character(climdatsub3$Block)
transdatsub$Block<-as.character(transdatsub$Block)
colnames(climdatsub3)[7]<-"meanGST"

alldat <- join(transdatsub, climdatsub3, by=c("Block","Canopy","Understory"), match="all")
alldat<-subset(alldat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi","Canopy","Understory","GDD_total","GDD_totaln","meanGST","Light_Mean","Light_GDD","MAT","snow_dur","snow_dur_cont"))
#Ok, so now I have climate data and survival/growth data in the same place. 
#separate out species
tsmedat<-alldat[alldat$Species=="TSME",]
tshedat<-alldat[alldat$Species=="TSHE",]
abamdat<-alldat[alldat$Species=="ABAM",]
tsmedat$PlantedStand=factor(tsmedat$PlantedStand)
tshedat$PlantedStand=factor(tshedat$PlantedStand)
abamdat$PlantedStand=factor(abamdat$PlantedStand)
tsmedat$OriginStand=factor(tsmedat$OriginStand)
tshedat$OriginStand=factor(tshedat$OriginStand)
abamdat$OriginStand=factor(abamdat$OriginStand)
tsmedat$Block=factor(tsmedat$Block)
tshedat$Block=factor(tshedat$Block)
abamdat$Block=factor(abamdat$Block)
tsmedat$snow_dur=as.numeric(tsmedat$snow_dur)
tshedat$snow_dur=as.numeric(tshedat$snow_dur)
abamdat$snow_dur=as.numeric(abamdat$snow_dur)
tsmedat$GDD_total=as.numeric(tsmedat$GDD_total)
tshedat$GDD_total=as.numeric(tshedat$GDD_total)
abamdat$GDD_total=as.numeric(abamdat$GDD_total)
tsmedat$meanGST=as.numeric(tsmedat$meanGST)
tshedat$meanGST=as.numeric(tshedat$meanGST)
abamdat$meanGST=as.numeric(abamdat$meanGST)
tsmedat$Light_Mean=as.numeric(tsmedat$Light_Mean)
tshedat$Light_Mean=as.numeric(tshedat$Light_Mean)
abamdat$Light_Mean=as.numeric(abamdat$Light_Mean)
tsmedat$snow_duration_cent=scale(as.numeric(tsmedat$snow_dur))
tshedat$snow_cover_duration_cent=scale(as.numeric(tshedat$snow_dur))
abamdat$snow_cover_duration_cent=scale(as.numeric(abamdat$snow_dur))
tsmedat$GDD_total_cent=scale(as.numeric(tsmedat$GDD_total))
tshedat$GDD_total_cent=scale(as.numeric(tshedat$GDD_total))
abamdat$GDD_total_cent=scale(as.numeric(abamdat$GDD_total))
tsmedat$GDD_totaln_cent=scale(as.numeric(tsmedat$GDD_totaln))
tshedat$GDD_totaln_cent=scale(as.numeric(tshedat$GDD_totaln))
abamdat$GDD_totaln_cent=scale(as.numeric(abamdat$GDD_totaln))

tsmedat$meanGST_cent=scale(as.numeric(tsmedat$meanGST))
tshedat$meanGST_cent=scale(as.numeric(tshedat$meanGST))
abamdat$meanGST_cent=scale(as.numeric(abamdat$meanGST))
tsmedat$Light_Mean_cent=scale(as.numeric(tsmedat$Light_Mean))
tshedat$Light_Mean_cent=scale(as.numeric(tshedat$Light_Mean))
abamdat$Light_Mean_cent=scale(as.numeric(abamdat$Light_Mean))
tsmedat$Light_GDD_cent=scale(as.numeric(tsmedat$Light_GDD))
tshedat$Light_GDD_cent=scale(as.numeric(tshedat$Light_GDD))
abamdat$Light_GDD_cent=scale(as.numeric(abamdat$Light_GDD))
tsmedat$MAT_cent=scale(as.numeric(tsmedat$MAT))
tshedat$MAT_cent=scale(as.numeric(tshedat$MAT))
abamdat$MAT_cent=scale(as.numeric(abamdat$MAT))
tsmedat$snow_dur_cent=scale(as.numeric(tsmedat$snow_dur))
tshedat$snow_dur_cent=scale(as.numeric(tshedat$snow_dur))
abamdat$snow_dur_cent=scale(as.numeric(abamdat$snow_dur))
tsmedat$snow_dur_cont_cent=scale(as.numeric(tsmedat$snow_dur_cont))
tshedat$snow_dur_cont_cent=scale(as.numeric(tshedat$snow_dur_cont))
abamdat$snow_dur_cont_cent=scale(as.numeric(abamdat$snow_dur_cont))
#Now fit models with microclimate data
tsmesurvdat<-subset(tsmedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#tsmesurvdat<- tsmesurvdat[apply(tsmesurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshesurvdat<-subset(tshedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#tshesurvdat<- tshesurvdat[apply(tshesurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamsurvdat<-subset(abamdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#abamsurvdat<- abamsurvdat[apply(abamsurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
#ABAM
constmod.abam<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=abamsurvdat)
#TSME
constmod.tsme<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tsmesurvdat)
#TSHE
constmod.tshe<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tshesurvdat)
#climate models
#MAT
matmod.abam<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=abamsurvdat)
matmod.tsme<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
matmod.tshe<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=tshesurvdat)
#Growing season length
gddmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gddmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gddmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=tshesurvdat)
gddmod2.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gddmod2.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gddmod2.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=tshesurvdat)
#Growing Season Temperature
gstmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gstmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gstmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=tshesurvdat)
#Snow duration
sdmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+snow_dur_cent:OriginStand, dist="lognormal", data=abamsurvdat)
sdmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+snow_dur_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
sdmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+snow_dur_cent:OriginStand, dist="lognormal", data=tshesurvdat)
sdmod2.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+snow_dur_cont_cent:OriginStand, dist="lognormal", data=abamsurvdat)
sdmod2.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+snow_dur_cont_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
sdmod2.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+snow_dur_cont_cent:OriginStand, dist="lognormal", data=tshesurvdat)
#light models
lightmod.abam<-survreg(Surv(time1,time2, type="interval2")~Light_Mean_cent+Light_Mean_cent:OriginStand, dist="lognormal", data=abamsurvdat)
lightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~Light_Mean_cent+Light_Mean_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
lightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~Light_Mean_cent+Light_Mean_cent:OriginStand, dist="lognormal", data=tshesurvdat)
lightmod2.abam<-survreg(Surv(time1,time2, type="interval2")~Light_GDD_cent+Light_GDD_cent:OriginStand, dist="lognormal", data=abamsurvdat)
lightmod2.tsme<-survreg(Surv(time1,time2, type="interval2")~Light_GDD_cent+Light_GDD_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
lightmod2.tshe<-survreg(Surv(time1,time2, type="interval2")~Light_GDD_cent+Light_GDD_cent:OriginStand, dist="lognormal", data=tshesurvdat)
# Light and GDD
constgddlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constgddlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constgddlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constgdd2lightmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constgdd2lightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constgdd2lightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constgddlight2mod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constgddlight2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constgddlight2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
constgdd2light2mod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constgdd2light2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constgdd2light2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
#Light and snow duration models
constsdlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constsdlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constsdlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constsdlight2mod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constsdlight2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constsdlight2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
constsd2lightmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constsd2lightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constsd2lightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constsd2light2mod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constsd2light2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constsd2light2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
#Light and GST models
constgstlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constgstlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constgstlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constgstlight2mod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constgstlight2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constgstlight2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
#Light and MAT models
constmatlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=abamsurvdat)
constmatlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tsmesurvdat)
constmatlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent, dist="lognormal", data=tshesurvdat)
constmatlight2mod.abam<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=abamsurvdat)
constmatlight2mod.tsme<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tsmesurvdat)
constmatlight2mod.tshe<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent, dist="lognormal", data=tshesurvdat)
#nul models
nullmod.abam<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=abamsurvdat)
nullmod.tsme<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=tsmesurvdat)
nullmod.tshe<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=tshesurvdat)
#model comparison to find best-fit model
AICtab(nullmod.abam,constmod.abam, matmod.abam,gddmod.abam,gddmod2.abam,sdmod.abam,sdmod2.abam,gstmod.abam,lightmod.abam,lightmod2.abam,constgddlightmod.abam,constgddlight2mod.abam,constgdd2lightmod.abam,constgdd2light2mod.abam,constgstlightmod.abam,constgstlight2mod.abam,constsdlightmod.abam,constsd2lightmod.abam,constsdlight2mod.abam,constsd2light2mod.abam,constmatlightmod.abam,constmatlight2mod.abam)#lowest aic for constmod.abam- 
AICtab(nullmod.tsme,constmod.tsme, matmod.tsme,gddmod.tsme,gddmod2.tsme,sdmod.tsme,sdmod2.tsme,gstmod.tsme,lightmod.tsme,lightmod2.tsme,constgddlightmod.tsme,constgddlight2mod.tsme,constgdd2lightmod.tsme,constgdd2light2mod.tsme,constgstlightmod.tsme,constgstlight2mod.tsme,constsdlightmod.tsme,constsd2lightmod.tsme,constsdlight2mod.tsme,constsd2light2mod.tsme,constmatlightmod.tsme,constmatlight2mod.tsme)#lowest aic for constmod.tsme- #
AICtab(nullmod.tshe,constmod.tshe, matmod.tshe,gddmod.tshe,gddmod2.tshe,sdmod.tshe,sdmod2.tshe,gstmod.tshe,lightmod.tshe,lightmod2.tshe,constgddlightmod.tshe,constgddlight2mod.tshe,constgdd2lightmod.tshe,constgdd2light2mod.tshe,constgstlightmod.tshe,constgstlight2mod.tshe,constsdlightmod.tshe,constsd2lightmod.tshe,constsdlight2mod.tshe,constsd2light2mod.tshe,constmatlightmod.tshe,constmatlight2mod.tshe)#lowest aic for constmod.tshe- 

#Now growth data and microclimate data
tsmegrdat<-subset(tsmedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#tsmegrdat<- tsmegrdat[apply(tsmegrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshegrdat<-subset(tshedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#tshegrdat<- tshegrdat[apply(tshegrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamgrdat<-subset(abamdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","GDD_total_cent","GDD_totaln_cent","meanGST_cent","Light_Mean_cent","Light_GDD_cent","MAT_cent","snow_dur_cent","snow_dur_cont_cent"))
#abamgrdat<- abamgrdat[apply(abamgrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na

#ABAM
consthimod.abam<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=abamgrdat)
#TSME
consthimod.tsme<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tsmegrdat)
#TSHE
consthimod.tshe<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tshegrdat)
Anova(consthimod.tsme, type="III")
Anova(consthimod.abam, type="III")
Anova(consthimod.tshe, type="III")
#climate models
#MAT
mathimod.abam<-lmer(annhi~OriginStand+MAT_cent+OriginStand:MAT_cent+(1|Block), REML=FALSE, data=abamgrdat)
mathimod.tsme<-lmer(annhi~OriginStand+MAT_cent+OriginStand:MAT_cent+(1|Block), REML=FALSE, data=tsmegrdat)
mathimod.tshe<-lmer(annhi~OriginStand+MAT_cent+OriginStand:MAT_cent+(1|Block), REML=FALSE, data=tshedat)

#Growing season length
gddhimod.abam<-lmer(annhi~OriginStand+GDD_total_cent+OriginStand:GDD_total_cent+(1|Block), REML=FALSE, data=abamgrdat)
gddhimod.tsme<-lmer(annhi~OriginStand+GDD_total_cent+OriginStand:GDD_total_cent+(1|Block), REML=FALSE, data=tsmegrdat)
gddhimod.tshe<-lmer(annhi~OriginStand+GDD_total_cent+OriginStand:GDD_total_cent+(1|Block), REML=FALSE, data=tshegrdat)
gddhimod2.abam<-lmer(annhi~OriginStand+GDD_totaln_cent+OriginStand:GDD_totaln_cent+(1|Block), REML=FALSE, data=abamgrdat)
gddhimod2.tsme<-lmer(annhi~OriginStand+GDD_totaln_cent+OriginStand:GDD_totaln_cent+(1|Block), REML=FALSE, data=tsmegrdat)
gddhimod2.tshe<-lmer(annhi~OriginStand+GDD_totaln_cent+OriginStand:GDD_totaln_cent+(1|Block), REML=FALSE, data=tshegrdat)

#GST
gsthimod.abam<-lmer(annhi~OriginStand+meanGST_cent+OriginStand:meanGST_cent+(1|Block), REML=FALSE, data=abamgrdat)
gsthimod.tsme<-lmer(annhi~OriginStand+meanGST_cent+OriginStand:meanGST_cent+(1|Block), REML=FALSE, data=tsmegrdat)
gsthimod.tshe<-lmer(annhi~OriginStand+meanGST_cent+OriginStand:meanGST_cent+(1|Block), REML=FALSE, data=tshedat)

#Snow
sdhimod.abam<-lmer(annhi~OriginStand+snow_dur_cent+OriginStand:snow_dur_cent+(1|Block), REML=FALSE, data=abamgrdat)
sdhimod.tsme<-lmer(annhi~OriginStand+snow_dur_cent+OriginStand:snow_dur_cent+(1|Block), REML=FALSE, data=tsmegrdat)
sdhimod.tshe<-lmer(annhi~OriginStand+snow_dur_cent+OriginStand:snow_dur_cent+(1|Block), REML=FALSE, data=tshegrdat)
sdhimod2.abam<-lmer(annhi~OriginStand+snow_dur_cent+OriginStand:snow_dur_cent+(1|Block), REML=FALSE, data=abamgrdat)
sdhimod2.tsme<-lmer(annhi~OriginStand+snow_dur_cont_cent+OriginStand:snow_dur_cont_cent+(1|Block), REML=FALSE, data=tsmegrdat)
sdhimod2.tshe<-lmer(annhi~OriginStand+snow_dur_cont_cent+OriginStand:snow_dur_cont_cent+(1|Block), REML=FALSE, data=tshegrdat)

#light models
lighthimod.abam<-lmer(annhi~OriginStand+Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
lighthimod.tsme<-lmer(annhi~OriginStand+Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
lighthimod.tshe<-lmer(annhi~OriginStand+Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
lighthimod2.abam<-lmer(annhi~OriginStand+Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
lighthimod2.tsme<-lmer(annhi~OriginStand+Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
lighthimod2.tshe<-lmer(annhi~OriginStand+Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)

#GDD Models AND Light models
constgddlighthimod.abam<-lmer(annhi~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgddlighthimod.tsme<-lmer(annhi~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgddlighthimod.tshe<-lmer(annhi~GDD_total_cent+OriginStand+Light_Mean_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constgddlight2himod.abam<-lmer(annhi~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgddlight2himod.tsme<-lmer(annhi~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgddlight2himod.tshe<-lmer(annhi~GDD_total_cent+OriginStand+Light_GDD_cent+GDD_total_cent:OriginStand+GDD_total_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
constgdd2lighthimod.abam<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgdd2lighthimod.tsme<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgdd2lighthimod.tshe<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_Mean_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constgdd2light2himod.abam<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgdd2light2himod.tsme<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgdd2light2himod.tshe<-lmer(annhi~GDD_totaln_cent+OriginStand+Light_GDD_cent+GDD_totaln_cent:OriginStand+GDD_totaln_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
#Snow AND Light models
constsdlighthimod.abam<-lmer(annhi~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constsdlighthimod.tsme<-lmer(annhi~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constsdlighthimod.tshe<-lmer(annhi~snow_dur_cent+OriginStand+Light_Mean_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constsdlight2himod.abam<-lmer(annhi~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constsdlight2himod.tsme<-lmer(annhi~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constsdlight2himod.tshe<-lmer(annhi~snow_dur_cent+OriginStand+Light_GDD_cent+snow_dur_cent:OriginStand+snow_dur_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
constsd2lighthimod.abam<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constsd2lighthimod.tsme<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constsd2lighthimod.tshe<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_Mean_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constsd2light2himod.abam<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constsd2light2himod.tsme<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constsd2light2himod.tshe<-lmer(annhi~snow_dur_cont_cent+OriginStand+Light_GDD_cent+snow_dur_cont_cent:OriginStand+snow_dur_cont_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
#GST and light models
constgstlighthimod.abam<-lmer(annhi~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgstlighthimod.tsme<-lmer(annhi~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgstlighthimod.tshe<-lmer(annhi~meanGST_cent+OriginStand+Light_Mean_cent+meanGST_cent:OriginStand+meanGST_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constgstlight2himod.abam<-lmer(annhi~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constgstlight2himod.tsme<-lmer(annhi~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constgstlight2himod.tshe<-lmer(annhi~meanGST_cent+OriginStand+Light_GDD_cent+meanGST_cent:OriginStand+meanGST_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
#MAT and light models
constmatlighthimod.abam<-lmer(annhi~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=abamgrdat)
constmatlighthimod.tsme<-lmer(annhi~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constmatlighthimod.tshe<-lmer(annhi~MAT_cent+OriginStand+Light_Mean_cent+MAT_cent:OriginStand+MAT_cent:Light_Mean_cent+OriginStand:Light_Mean_cent+(1|Block), REML=FALSE, data=tshegrdat)
constmatlight2himod.abam<-lmer(annhi~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=abamgrdat)
constmatlight2himod.tsme<-lmer(annhi~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tsmegrdat)
constmatlight2himod.tshe<-lmer(annhi~MAT_cent+OriginStand+Light_GDD_cent+MAT_cent:OriginStand+MAT_cent:Light_GDD_cent+OriginStand:Light_GDD_cent+(1|Block), REML=FALSE, data=tshegrdat)
nullhimod.abam<-lmer(annhi~1+(1|Block), REML=FALSE, data=abamgrdat)
nullhimod.tsme<-lmer(annhi~1+(1|Block), REML=FALSE, data=tsmegrdat)
nullhimod.tshe<-lmer(annhi~1+(1|Block), REML=FALSE, data=tshegrdat)

##########################
####Germination Models####
##########################
germdat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
germdat[which(germdat$TotalGerms>50 &germdat$SeedsAdded==50),]$TotalGerms=50
##First, select only rows in which seeds were added
addat<-germdat[germdat$SeedsAdded>0,]
#select out 2011 data for analysis, since sample size was much higher in this year (i.e. more seeds were added
dat2011<-addat[addat$Year=="2011",]
colnames(dat2011)[2]<-"Stand"
#divide data by species
abamgermdat<-dat2011[dat2011$SpPlant=="ABAM",]
tsmegermdat<-dat2011[dat2011$SpPlant=="TSME",]
tshegermdat<-dat2011[dat2011$SpPlant=="TSHE",]

abamgermdat$Stand=factor(abamgermdat$Stand)
tsmegermdat$Stand=factor(tsmegermdat$Stand)
tshegermdat$Stand=factor(tshegermdat$Stand)
abamgermdat$Origin=factor(abamgermdat$Origin)
tsmegermdat$Origin=factor(tsmegermdat$Origin)
tshegermdat$Origin=factor(tshegermdat$Origin)
tsmegermdat$Block=factor(tsmegermdat$Block)
tshegermdat$Block=factor(tshegermdat$Block)
abamgermdat$Block=factor(abamgermdat$Block)

#ABAM
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)

#TSME
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)

#Tshe
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

#fit categorical models
germod.abam<-glmmadmb(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=abamgermdat, family="binomial")
germod.tsme<-glmmadmb(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="binomial")
germod.tshe<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="binomial")

#Add in moss data:
grcov<-read.csv("data/transplantcover.csv", header=T )
colnames(dat2011)[2]<-"Elev"
germdat2<-join(dat2011, grcov, by=c("Block","Plot","Canopy"), match="all")
abamgermdat$Moss_cent<-scale(abamgermdat$Moss)
germod.abam.moss<-glmmadmb(abamy~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
summary(germod.abam.moss)
germod.tshe<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="binomial")
Anova(germod.tshe,type="III")
tshegermdat$Moss_cent<-scale(tshegermdat$Moss)
germod.tshe.moss<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss+(1|Block), data=tshegermdat, family="binomial")
Anova(germod.tshe.moss,type="III")
colnames(climdatsub)[3]<-"Stand"
#now add in microclimate data
allgerdat <- join(dat2011, climdatsub, by=c("Stand","Block","Canopy","Understory"), match="first")
allgerdat<-subset(allgerdat,select=c("SpPlant","Block","Stand","Origin","Canopy","Understory","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean","Moss","TotalGerms","TotalFails"))
#climate models
abamgermdat<-allgerdat[allgerdat$SpPlant=="ABAM",]
tsmegermdat<-allgerdat[allgerdat$SpPlant=="TSME",]
tshegermdat<-allgerdat[allgerdat$SpPlant=="TSHE",]

abamgermdat$Stand=factor(abamgermdat$Stand)
tsmegermdat$Stand=factor(tsmegermdat$Stand)
tshegermdat$Stand=factor(tshegermdat$Stand)
abamgermdat$Origin=factor(abamgermdat$Origin)
tsmegermdat$Origin=factor(tsmegermdat$Origin)
tshegermdat$Origin=factor(tshegermdat$Origin)
tsmegermdat$TotalGerms=as.numeric(tsmegermdat$TotalGerms)
tshegermdat$TotalGerms=as.numeric(tshegermdat$TotalGerms)
abamgermdat$TotalGerms=as.numeric(abamgermdat$TotalGerms)
tsmegermdat$Block=factor(tsmegermdat$Block)
tshegermdat$Block=factor(tshegermdat$Block)
abamgermdat$Block=factor(abamgermdat$Block)
tsmegermdat$Moss=as.numeric(tsmegermdat$Moss)
tshegermdat$Moss=as.numeric(tshegermdat$Moss)
abamgermdat$Moss=as.numeric(abamgermdat$Moss)
#ABAM
#make sure that datasets are the same, regardless of explanatory variables
abamgermdat<- abamgermdat[apply(abamgermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)

#TSME
#make sure that datasets are the same, regardless of explanatory variables
tsmegermdat<- tsmegermdat[apply(tsmegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)

#Tshe
tshegermdat<- tshegermdat[apply(tshegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)
#Fit categorical models again, with reduced dataset
tsmegermdat$Moss_cent=scale(as.numeric(tsmegermdat$Moss))
tshegermdat$Moss_cent=scale(as.numeric(tshegermdat$Moss))
abamgermdat$Moss_cent=scale(as.numeric(abamgermdat$Moss))

constgermod.tsme<-glmmadmb(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="binomial")
constgermod.abam<-glmmadmb(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constgermod.tshe<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

tsmegermdat$snow_cover_duration=as.numeric(tsmegermdat$snow_cover_duration)
tshegermdat$snow_cover_duration=as.numeric(tshegermdat$snow_cover_duration)
abamgermdat$snow_cover_duration=as.numeric(abamgermdat$snow_cover_duration)
tsmegermdat$GDD_total=as.numeric(tsmegermdat$GDD_total)
tshegermdat$GDD_total=as.numeric(tshegermdat$GDD_total)
abamgermdat$GDD_total=as.numeric(abamgermdat$GDD_total)
tsmegermdat$meanGST=as.numeric(tsmegermdat$meanGST)
tshegermdat$meanGST=as.numeric(tshegermdat$meanGST)
abamgermdat$meanGST=as.numeric(abamgermdat$meanGST)
tsmegermdat$Light_Mean=as.numeric(tsmegermdat$Light_Mean)
tshegermdat$Light_Mean=as.numeric(tshegermdat$Light_Mean)
abamgermdat$Light_Mean=as.numeric(abamgermdat$Light_Mean)
tsmegermdat$snow_cover_duration_cent=scale(as.numeric(tsmegermdat$snow_cover_duration))
tshegermdat$snow_cover_duration_cent=scale(as.numeric(tshegermdat$snow_cover_duration))
abamgermdat$snow_cover_duration_cent=scale(as.numeric(abamgermdat$snow_cover_duration))
tsmegermdat$GDD_total_cent=scale(as.numeric(tsmegermdat$GDD_total))
tshegermdat$GDD_total_cent=scale(as.numeric(tshegermdat$GDD_total))
abamgermdat$GDD_total_cent=scale(as.numeric(abamgermdat$GDD_total))
tsmegermdat$meanGST_cent=scale(as.numeric(tsmegermdat$meanGST))
tshegermdat$meanGST_cent=scale(as.numeric(tshegermdat$meanGST))
abamgermdat$meanGST_cent=scale(as.numeric(abamgermdat$meanGST))
tsmegermdat$Light_Mean_cent=scale(as.numeric(tsmegermdat$Light_Mean))
tshegermdat$Light_Mean_cent=scale(as.numeric(tshegermdat$Light_Mean))
abamgermdat$Light_Mean_cent=scale(as.numeric(abamgermdat$Light_Mean))
#Now climate models
gddgermod.abam<-glmmadmb(abamy~Origin+GDD_total_cent+Origin:GDD_total_cent+Moss_cent+(1|Block),data=abamgermdat, family="binomial")
gddgermod.tsme<-glmmadmb(tsmey~Origin+GDD_total_cent+Origin:GDD_total_cent+(1|Block),data=tsmegermdat, family="binomial")
gddgermod.tshe<-glmmadmb(tshey~Origin+GDD_total_cent+Origin:GDD_total_cent+Moss_cent+(1|Block),data=tshegermdat, family="binomial")
gstgermod.abam<-glmmadmb(abamy~Origin+meanGST_cent+Origin:meanGST_cent+Moss_cent+(1|Block),data=abamgermdat, family="binomial")
gstgermod.tsme<-glmmadmb(tsmey~Origin+meanGST_cent+Origin:meanGST_cent+(1|Block), data=tsmegermdat, family="binomial")
gstgermod.tshe<-glmmadmb(tshey~Origin+meanGST_cent+Origin:meanGST_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
sdgermod.abam<-glmmadmb(abamy~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
sdgermod.tsme<-glmmadmb(tsmey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+(1|Block), data=tsmegermdat, family="binomial")
sdgermod.tshe<-glmmadmb(tshey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

#light models
lightgermod.abam<-glmmadmb(abamy~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
lightgermod.tsme<-glmmadmb(tsmey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tsmegermdat, family="binomial")
lightgermod.tshe<-glmmadmb(tshey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

#Try replacing canopy and understory with "light" AND planted stand with GDD- this might tell us how much competition is due to light
constgddlightgermod.abam<-glmmadmb(abamy~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constgddlightgermod.tsme<-glmmadmb(tsmey~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
constgddlightgermod.tshe<-glmmadmb(tshey~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

constsdlightgermod.abam<-glmmadmb(abamy~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constsdlightgermod.tsme<-glmmadmb(tsmey~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
constsdlightgermod.tshe<-glmmadmb(tshey~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
#Try replacing canopy and understory with "light" AND planted stand with meanGST- this might tell us how much competition is due to light
constgstlightgermod.abam<-glmmadmb(abamy~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constgstlightgermod.tsme<-glmmadmb(tsmey~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
constgstlightgermod.tshe<-glmmadmb(tshey~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
#compare models
AICctab(constgermod.abam,sdgermod.abam,lightgermod.abam,gstgermod.abam,gddgermod.abam,constgstlightgermod.abam,constgddlightgermod.abam,constsdlightgermod.abam)#lowest aic for constgddhimod.abam
AICctab(constgermod.tsme,sdgermod.tsme,lightgermod.tsme,gstgermod.tsme,gddgermod.tsme,constgstlightgermod.tsme,constgddlightgermod.tsme,constsdlightgermod.tsme)#lowest aic for constsdhimod.tsme 
AICctab(constgermod.tshe,sdgermod.tshe,lightgermod.tshe,gstgermod.tshe,gddgermod.tshe,constgstlightgermod.tshe,constgddlightgermod.tshe,constsdlightgermod.tshe)#lowest aic for constgddhimod.tshe 

