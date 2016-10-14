#Analyses suggested by reviewers for GCB Transplant Study
#Started by Ailene Ettinger on October 6, 2016
setwd("~/git/mora_transplant")
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library(dplyr)
library(plyr)
library(survival)
library(car)
library(lme4)
library(bbmle)

#Begin with transplant data: get the survival times, growth variables, etc
transdat<-read.csv("data/2013TransplantStatusHeight(October).csv", header=TRUE)#csv file has been sorted so that all dates appear together, all status columns are together, etc
dim(transdat)#3959 rows(=individuals),  39 columns
microclim<-read.csv("data/AllStands_clim1.csv", header=T)#this is just data from 2012; 
dim(microclim)#136 rows, 15 columns
transdat$PlantedStand2<-as.numeric(transdat$PlantedStand)
transdat$PlantedStand<-as.factor(transdat$PlantedStand)
transdat$OriginStand<-as.factor(transdat$OriginStand)
transdat$Block<-as.factor(transdat$Block)
transdat$UniqueID<-as.factor(transdat$UniqueID)
transdat$Date1<-as.Date(transdat$Date1,format='%m/%d/%Y')
transdat$Date2<-as.Date(transdat$Date2,format='%m/%d/%Y')
transdat$Date3<-as.Date(transdat$Date3,format='%m/%d/%Y')
transdat$Date4<-as.Date(transdat$Date4,format='%m/%d/%Y')
transdat$Date5<-as.Date(transdat$Date5,format='%m/%d/%Y')
transdat$CompAmt.1<-as.factor(transdat$CompAmt.1)
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
#Now add microclimate data, including GDD, snowduration, and mean light for each block
climdat<-read.csv("data/AllStands_clim1.csv", header=TRUE)
climdat$Block<-as.factor(climdat$Block)
#Combine the climate data with the transplant data columns that we want to use for survival and growth analysis
transdatsub<-subset(transdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi"))
climdatsub<-subset(climdat,select=c("Stand","Block","Elevation_m","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
colnames(climdatsub)[3]<-"PlantedStand"
alldat <- join(transdatsub, climdatsub, by=c("PlantedStand","Block","Canopy","Understory"), match="first")
alldat<-subset(alldat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))

###remove any rows that don't contain all the data, so that we can use AICc to compare models

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
#Next step is to modify the models so that the explanatory variables are GDD and light instead of range position/treatment and 
#Previous models:
#in previous preliminary analyses (available upon request), I used model selection to identify that the lognormal distirbution is best-fit for these data
#make sure that datasets are the same, regardless of explanatory variables
tsmesurvdat<-subset(tsmedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
tsmesurvdat<- tsmesurvdat[apply(tsmesurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshesurvdat<-subset(tshedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
tshesurvdat<- tshesurvdat[apply(tshesurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamsurvdat<-subset(abamdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
abamsurvdat<- abamsurvdat[apply(abamsurvdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
#ABAM
constmod.abam<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=abamsurvdat)
#TSME
constmod.tsme<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tsmesurvdat)
#TSHE
constmod.tshe<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tshesurvdat)
Anova(constmod.tsme, test.statistic="LR", type="III")
Anova(constmod.abam, test.statistic="LR", type="III")
Anova(constmod.tshe,test.statistic="LR", type="III")
#For new models, use only light and gdd
#does not converge with gdd or with snow duration AND light- "ran out of iterations"- perhaps because of discrete nature of the data? its continuous data, but we only have it in chunks- it can be fit on its known, but not as an interaction with light
#climmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST*Light_Mean*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#climmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST*Light_Mean*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#climmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST*Light_Mean*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#climcanmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST*Canopy*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#climcanmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST*Canopy*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#climcanmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST*Canopy*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#climundmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST*Understory*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#climundmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST*Understory*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#climundmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST*Understory*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#gddmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#gddmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#gddmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#lightmod.abam<-survreg(Surv(time1,time2, type="interval2")~Light_Mean*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#lightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~Light_Mean*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#lightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~Light_Mean*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#gstmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#gstmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#gstmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#sdmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration*OriginStand, dist="lognormal", data=abamdat,maxiter=1000)
#sdmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration*OriginStand, dist="lognormal", data=tshedat,maxiter=1000)
#sdmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration*OriginStand, dist="lognormal", data=tsmedat,maxiter=1000)
#Try replacing canopy and understory with "light"- this might tell us how much competition is due to light
constlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=abamsurvdat)
constlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tsmesurvdat)
constlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tshesurvdat)
#try replacing planted stand with GDD
constgddmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory, dist="lognormal", data=abamsurvdat)
constgddmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory, dist="lognormal", data=tsmesurvdat)
constgddmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory, dist="lognormal", data=tshesurvdat)
#try replacing planted stand with snow duration
constsdmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory, dist="lognormal", data=abamsurvdat)
constsdmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory, dist="lognormal", data=tsmesurvdat)
constsdmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory, dist="lognormal", data=tshesurvdat)

##The below do not converge
#Try replacing canopy and understory with "light" AND planted stand with GDD- this might tell us how much competition is due to light
#constgddlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=abamdat)
#constgddlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tsmedat)
#constgddlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tshedat)
#Try replacing canopy and understory with "light" AND planted stand with sd- this might tell us how much competition is due to light
#constsdlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=abamdat)
#constsdlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tsmedat)
#constsdlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tshedat)
#Try replacing canopy and understory with "light" AND planted stand with meanGST- this might tell us how much competition is due to light
constgddlightmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=abamsurvdat)
constgddlightmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tsmesurvdat)
constgddlightmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean, dist="lognormal", data=tshesurvdat)


#AIC(climmod.abam,climcanmod.abam,climundmod.abam,constmod.abam, gddmod.abam,lightmod.abam,gstmod.abam,sdmod.abam,constlightmod.abam,constgddmod.abam,constsdmod.abam,constgddlightmod.abam)#lowest aic for constlightmod.abam- suggests that there are other important differences across elevations besides GDD/SD/temperature; but that light is main driver of competition
#AIC(climmod.tsme,climcanmod.tsme,climundmod.tsme,constmod.tsme,gddmod.tsme,lightmod.tsme,gstmod.tsme,sdmod.tsme,constlightmod.tsme,constgddmod.tsme,constsdmod.tsme,constgddlightmod.tsme)#lowest aic for constlightmod.tsme
#AIC(climmod.tshe,climcanmod.tshe,climundmod.tshe,constmod.tshe,gddmod.tshe,lightmod.tshe,gstmod.tshe,sdmod.tshe,constlightmod.tshe,constgddmod.tshe,constsdmod.tshe,constgddlightmod.tshe)#lowest aic for constmod- other differences across elevations not explained by light and climate variables- light is not important driver of competition (other resources)
AIC(constmod.abam, constlightmod.abam,constgddmod.abam,constsdmod.abam,constgddlightmod.abam)#lowest aic for constmod.abam- 
AIC(constmod.tsme,constlightmod.tsme,constgddmod.tsme,constsdmod.tsme,constgddlightmod.tsme)#lowest aic for constmod.tsme
AIC(constmod.tshe,constlightmod.tshe,constgddmod.tshe,constsdmod.tshe,constgddlightmod.tshe)#lowest aic for constmod- other differences across elevations not explained by light and climate variables- light is not important driver of competition (other resources)

#Look at growth data
#make sure that datasets are the same, regardless of explanatory variables
tsmegrdat<-subset(tsmedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
tsmegrdat<- tsmegrdat[apply(tsmegrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshegrdat<-subset(tshedat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
tshegrdat<- tshegrdat[apply(tshegrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamgrdat<-subset(abamdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","annhi","Canopy","Understory","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
abamgrdat<- abamgrdat[apply(abamgrdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na

#ABAM
consthimod.abam<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=abamgrdat)
#TSME
consthimod.tsme<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tsmegrdat)
#TSHE
consthimod.tshe<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block), data=tshegrdat)
Anova(consthimod.tsme, type="III")
Anova(consthimod.abam, type="III")
Anova(consthimod.tshe, type="III")
#Try replacing canopy and understory with "light"- this might tell us how much competition is due to light
constlighthimod.abam<-lmer(annhi~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamgrdat)
constlighthimod.tsme<-lmer(annhi~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmegrdat)
constlighthimod.tshe<-lmer(annhi~PlantedStand+OriginStand+Light_Mean+PlantedStand:OriginStand+PlantedStand:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshegrdat)
#try replacing planted stand with GDD
constgddhimod.abam<-lmer(annhi~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamgrdat)
constgddhimod.tsme<-lmer(annhi~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmegrdat)
constgddhimod.tshe<-lmer(annhi~GDD_total+OriginStand+Canopy+Understory+GDD_total:OriginStand+GDD_total:Canopy+GDD_total:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+GDD_total:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshegrdat)
#try replacing planted stand with snow duration
constsdhimod.abam<-lmer(annhi~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamgrdat)
constsdhimod.tsme<-lmer(annhi~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmegrdat)
constsdhimod.tshe<-lmer(annhi~snow_cover_duration+OriginStand+Canopy+Understory+snow_cover_duration:OriginStand+snow_cover_duration:Canopy+snow_cover_duration:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+snow_cover_duration:Canopy:Understory+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshegrdat)

##The below do not converge
#Try replacing canopy and understory with "light" AND planted stand with GDD- this might tell us how much competition is due to light
#constgddlighthimod.abam<-lmer(annhi~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamdat)
#constgddlighthimod.tsme<-lmer(annhi~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmedat)
#constgddlighthimod.tshe<-lmer(annhi~GDD_total+OriginStand+Light_Mean+GDD_total:OriginStand+GDD_total:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshedat)
#Try replacing canopy and understory with "light" AND planted stand with sd- this might tell us how much competition is due to light
#constsdlighthimod.abam<-lmer(annhi~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamdat)
#constsdlighthimod.tsme<-lmer(annhi~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmedat)
#constsdlighthimod.tshe<-lmer(annhi~snow_cover_duration+OriginStand+Light_Mean+snow_cover_duration:OriginStand+snow_cover_duration:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshedat)
#Try replacing canopy and understory with "light" AND planted stand with meanGST- this might tell us how much competition is due to light
constgddlighthimod.abam<-lmer(annhi~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=abamgrdat)
constgddlighthimod.tsme<-lmer(annhi~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tsmegrdat)
constgddlighthimod.tshe<-lmer(annhi~meanGST+OriginStand+Light_Mean+meanGST:OriginStand+meanGST:Light_Mean+OriginStand:Light_Mean+(1|Block), REML=FALSE, contrasts=c(unordered="contr.sum", ordered="contr.poly"), data=tshegrdat)


AICctab(consthimod.abam,constlighthimod.abam,constgddhimod.abam,constsdhimod.abam)#lowest aic for constgddhimod.abam
AICctab(consthimod.tsme,constlighthimod.tsme,constgddhimod.tsme,constsdhimod.tsme)#lowest aic for constsdhimod.tsme 
AICctab(consthimod.tshe,constlighthimod.tshe,constgddhimod.tshe,constsdhimod.tshe)#lowest aic for constgddhimod.tshe 

##Greenhouse germination
ggermdat<-read.csv("data/AllGreenhouseGerm.csv", head=T)
ggermdat$grate<-ggermdat$SeedsGerm/ggermdat$SeedsSown
ggermdat$SeedsFail<-ggermdat$SeedsSown-ggermdat$SeedsGerm
ggermdat$yrsstore<-ggermdat$YearSown-ggermdat$YearCollected# number of years stored before sowing
tsmeggerm<-ggermdat[ggermdat$Species=="TSME",]
tsheggerm<-ggermdat[ggermdat$Species=="TSHE",]
abamggerm<-ggermdat[ggermdat$Species=="ABAM",]
gmod.tsme<-glm(cbind(tsmeggerm$SeedsGerm,tsmeggerm$SeedsFail)~SourceStand+as.factor(YearCollected)+yrsstore, data=tsmeggerm, family=binomial)
summary(gmod.tsme)
Anova(gmod.tsme)
gmod.tshe<-glm(cbind(tsheggerm$SeedsGerm,tsheggerm$SeedsFail)~SourceStand+as.factor(YearCollected)+yrsstore, data=tsheggerm, family=binomial)
summary(gmod.tshe)
Anova(gmod.tshe)
gmod.abam<-glm(cbind(abamggerm$SeedsGerm,abamggerm$SeedsFail)~SourceStand+as.factor(YearCollected)+yrsstore, data=abamggerm, family=binomial)
summary(gmod.abam)
Anova(gmod.abam)
#years stored has negative effect on germination for TSHE but no other species
ggermdat09<-ggermdat[which(ggermdat$YearCollected==2009),]
ggermdat10<-ggermdat[which(ggermdat$YearCollected==2010),]

gmean09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), mean,na.rm=T)
gmean10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), mean,na.rm=T)
gsd09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), sd,na.rm=T)
gsd10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), sd,na.rm=T)
gn09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), length)
gn10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), length)
gse09<-gsd09/sqrt(gn09)
gse10<-gsd10/sqrt(gn10)
#Figure of germination in the greenhouse, for supplemental
quartz(height=7,width=4)
par(mfcol=c(3,1),mai=c(.6,.7,.2,.1), omi=c(.7,.01,.2,.1))
X<-c(1,2,3,4,5)
#TSME
plot(gmean09[3,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("a",side=3,line=1, adj=-.1)
for (i in 1:5){
  arrows(X[i],gmean09[3,i,]-gse09[3,i,],X[i],gmean09[3,i,]+gse09[3,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[3,i,1]-gse10[3,i,1],X[i],gmean10[3,i,1]+gse10[3,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[3,i,2]-gse10[3,i,2],X[i],gmean10[3,i,2]+gse10[3,i,2],length=.05,angle=90,code=0)}
points(gmean09[3,,]~X,pch=21,cex=1.5,bg=c("black","gray","white"))
points(gmean10[3,,1]~X,pch=22,cex=1.5,bg=c("black","gray","white"))
points(gmean10[3,,2]~X,pch=23,cex=1.5,bg=c("black","gray","white"))
legend(1,.6,legend=c("2009","2010-1", "2010-2"),bty="n",pch=c(21,22,23),pt.bg=c("black"),angle=45,cex=1.1, pt.cex=1.5)

##ABAM
plot(gmean09[1,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("b",side=3,line=1, adj=-.1)
for (i in 1:5){
  arrows(X[i],gmean09[1,i,]-gse09[1,i,],X[i],gmean09[1,i,]+gse09[1,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[1,i,1]-gse10[1,i,1],X[i],gmean10[1,i,1]+gse10[1,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[1,i,2]-gse10[1,i,2],X[i],gmean10[1,i,2]+gse10[1,i,2],length=.05,angle=90,code=0)}
points(gmean09[1,,]~X,pch=21,cex=1.5,bg=c("black","gray","gray","white"))
points(gmean10[1,,1]~X,pch=22,cex=1.5,bg=c("black","gray","gray","white"))
points(gmean10[1,,2]~X,pch=23,cex=1.5,bg=c("black","gray","gray","white"))
mtext("Proportion Seeds Germinating", side=2, line=3)

##TSHE
plot(gmean09[2,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
for (i in 1:5){
  arrows(X[i],gmean09[2,i,]-gse09[2,i,],X[i],gmean09[2,i,]+gse09[2,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[2,i,1]-gse10[2,i,1],X[i],gmean10[2,i,1]+gse10[2,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[2,i,2]-gse10[2,i,2],X[i],gmean10[2,i,2]+gse10[2,i,2],length=.05,angle=90,code=0)}
points(gmean09[2,,]~X,pch=21,cex=1.5,bg=c("black","gray","gray","white"))
points(gmean10[2,,1]~X,pch=22,cex=1.5,bg=c("black","gray","gray","white"))
points(gmean10[2,,2]~X,pch=23,cex=1.5,bg=c("black","gray","gray","white"))
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("c",side=3,line=1, adj=-.1)
mtext ("Source elevation (m)", line=-13, cex=.9)
axis(1, at = c(1,2,3,4,5), labels = colnames(gmean09), tick = FALSE, cex.axis=1.1, line=-.5)

####"Seed bed" analysis for Reviewer 2. 
grcov<-read.csv("data/perccover_MORAplots.csv", header=T )
grcov2<-subset(grcov,select=c("Name","Stand","Block","Gap.","Plot","Litter","Moss","Wood","Feces","Rush.1","Sand.Soil","Rock.1","Root","Bare.1"))
###do some cleaning of the cover data
grcov2$Litter<-as.numeric(grcov2$Litter)
grcov2[which(grcov2$Moss=="5 ground, 5 on log/wood"),]$Moss<-5
grcov2[which(grcov2$Moss=="<1"),]$Moss<-0
grcov2$Moss<-as.numeric(grcov2$Moss)
grcov2[which(grcov2$Wood=="<1"),]$Wood<-0
grcov2[which(grcov2$Wood=="1 (root)"),]$Root<-1
grcov2[which(grcov2$Wood=="1 (root)"),]$Wood<-0
grcov2[which(grcov2$Wood=="1 root, 1 trunk"),]$Root<-2
grcov2[which(grcov2$Wood=="1 root, 1 trunk"),]$Wood<-0
grcov2[which(grcov2$Wood=="10 (rotten)"),]$Wood<-10
grcov2[which(grcov2$Wood=="2 (rotten)"),]$Wood<-2
grcov2[which(grcov2$Wood=="50 (rotting)"),]$Wood<-50
grcov2[which(grcov2$Wood=="2 (8 = moss covered log)"),]$Wood<-2
grcov2[which(grcov2$Wood=="2 (nurse log)"),]$Wood<-2
grcov2[which(grcov2$Wood=="5 (covered in moss)"),]$Wood<-5
grcov2$Wood<-as.numeric(grcov2$Wood)
grcov2[which(grcov2$Plot=="B (formerly A)"),]$Plot<-"B"
grcov2[which(grcov2$Plot=="C (formerly B)"),]$Plot<-"C"
grcov2[which(grcov2$Plot=="D (formerly C)"),]$Plot<-"D"
grcov2[which(grcov2$Plot=="G (formerly D)"),]$Plot<-"G"
grcov2[which(grcov2$Plot=="H (formerly E)"),]$Plot<-"H"
grcov2[which(grcov2$Plot=="I (formerly F)"),]$Plot<-"I"
grcov2[which(grcov2$Plot=="J (formerly F)"),]$Plot<-"J"
grcov2[which(grcov2$Plot=="G (formerly E)"),]$Plot<-"G"
grcov2[which(grcov2$Plot=="H (formerly F)"),]$Plot<-"H"
grcov2[which(grcov2$Plot=="C (formerly A)"),]$Plot<-"C"
grcov2[which(grcov2$Plot=="D (formerly B)"),]$Plot<-"D"
grcov2[which(grcov2$Plot=="E (formerly C)"),]$Plot<-"E"
grcov2[which(grcov2$Plot=="H (formerly D)"),]$Plot<-"H"
grcov2[which(grcov2$Plot=="I (formerly E)"),]$Plot<-"C"
#are there trends in moss by canopy status, and/or by elevation?

tapply(grcov2$Moss,list(grcov2$Stand,grcov2$Canopy), mean, na.rm=T)
###Problem: moss cover is higher in gaps!
germdat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
#head(germdat)
#Get columns of grcov for Canopy and Block to match format of germdat so that i can merge the two files
colnames(grcov2)[3]<-"BlockNum"
grcov2$Block<-paste(grcov2$Stan,grcov2$BlockNum,sep="")
grcov2$Canopy<-NA
grcov2[grcov2$Gap.=="Y",]$Canopy<-"CompAbsent"
grcov2[grcov2$Gap.=="N",]$Canopy<-"CompPresent"
germdat2<-join(germdat, grcov2, by=c("Block","Plot","Canopy"), match="first")
colnames(germdat2)[5]<-"Elev"
dat<-subset(germdat2,select=c("Stand","Elev","Block","Name","Plot","Canopy","Understory","Comp","SpPlant","SpObs","Origin","Year","TotalGerms","SeedsAdded","TotalFails","Litter","Moss","Wood","Feces","Rush.1","Sand.Soil","Rock.1","Root","Bare.1"))
dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms=50
##First, select only rows in which seeds were added
addat<-dat[dat$SeedsAdded>0,]
#select out 2011 data for analysis
dat2011<-addat[addat$Year=="2011",]
#divide data by species
abamgermdat<-dat2011[dat2011$SpPlant=="ABAM",]
tsmegermdat<-dat2011[dat2011$SpPlant=="TSME",]
tshegermdat<-dat2011[dat2011$SpPlant=="TSHE",]
#remove NAs
#tsmegermdat<-tsmegermdat[1:min(which(is.na(tsmegermdat$Stand)))-1,]
#tshegermdat<-tshegermdat[1:min(which(is.na(tshegermdat$Stand)))-1,]
#abamgermdat<-abamgermdat[1:min(which(is.na(abamgermdat$Stand)))-1,]

abamgermdat$Stand=factor(abamgermdat$Stand)
tsmegermdat$Stand=factor(tsmegermdat$Stand)
tshegermdat$Stand=factor(tshegermdat$Stand)
abamgermdat$Origin=factor(abamgermdat$Origin)
tsmegermdat$Origin=factor(tsmegermdat$Origin)
tshegermdat$Origin=factor(tshegermdat$Origin)
tsmegermdat$TotalGerms=factor(tsmegermdat$TotalGerms)
tshegermdat$TotalGerms=factor(tshegermdat$TotalGerms)
abamgermdat$TotalGerms=factor(abamgermdat$TotalGerms)
tsmegermdat$Block=factor(tsmegermdat$Block)
tshegermdat$Block=factor(tshegermdat$Block)
abamgermdat$Block=factor(abamgermdat$Block)
tsmegermdat$Moss=as.numeric(tsmegermdat$Moss)
tshegermdat$Moss=as.numeric(tshegermdat$Moss)
abamgermdat$Moss=as.numeric(abamgermdat$Moss)
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)
#ABAM
#look at effect of moss on germination, just by itself
abam.moss<-glm(abamy~Stand*Canopy*as.numeric(Moss),data=abamgermdat, family=binomial)
summary(abam.moss)
Anova(abam.moss)#moss has negative effecton germination for ABAM
tshe.moss<-glm(tshey~as.numeric(Moss),data=tshegermdat, family=binomial)
summary(tshe.moss)
Anova(tshe.moss)#moss has positive effect on germination for TSHE
tsme.moss<-glm(tsmey~as.numeric(Moss),data=tsmegermdat, family=binomial)
summary(tsme.moss)
Anova(tsme.moss)#moss has no effect on germination for TSME
#so, moss seems to matter for some species but not others. since moss is higher in gaps than nongaps, I'll try 
#I'll try adding moss in to the model selection process
#ABAM
#make sure that datasets are the same, regardless of explanatory variables
abamgermdat<- abamgermdat[apply(abamgermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)

abammod<- glm(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod<- glmer(abamy~-1+Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory + (1|Block), data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
Anova(abam.mmod, test.statistic="Chisq", type="III")#model failed to converge
##Model selection, by hand
abammod1.moss<- glm(abamy~Stand+Origin+Canopy+Understory+Moss+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory, data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod1.moss<- glmer(abamy~Stand+Origin+Canopy+Understory+Moss+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+ (1|Block), data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod2.moss<- glmer(abamy~Stand*Origin +Moss+ (1|Block), data=abamgermdat, family=binomial)# still failed to converge
abammod2.moss<- glm(abamy~Stand*Origin+Moss, data=abamgermdat, family=binomial)#
abam.mmod2a.moss<- glmer(abamy~Canopy+Stand*Origin+Moss+ (1|Block), data=abamgermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
abammod2a.moss<- glm(abamy~Canopy+Stand*Origin+Moss, data=abamgermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
abammod3.moss<- glm(abamy~Canopy*Understory*Origin+Moss, data=abamgermdat, family=binomial)# 
abam.mmod3.moss<- glmer(abamy~Canopy*Understory*Origin+Moss+(1|Block), data=abamgermdat, family=binomial)# 
abammod4.moss<- glm(abamy~Canopy*Stand+Moss, data=abamgermdat, family=binomial)# 
abammod4a.moss<- glm(abamy~Origin+Canopy*Stand+Moss, data=abamgermdat, family=binomial)# 
abam.mmod4.moss<- glmer(abamy~Canopy*Stand+Moss+(1|Block), data=abamgermdat, family=binomial)# 
abam.mmod4a.moss<- glmer(abamy~Origin+Canopy*Stand+Moss+(1|Block), data=abamgermdat, family=binomial)# 
abammod5.moss<- glm(abamy~Understory*Stand+Moss, data=abamgermdat, family=binomial)# 
abam.mmod5.moss<- glmer(abamy~Understory*Stand+Moss+(1|Block), data=abamgermdat, family=binomial)# 
abammod6.moss<- glm(abamy~Comp*Stand+Moss, data=abamgermdat, family=binomial)# 
abam.mmod6.moss<- glmer(abamy~Comp*Stand+Moss + (1|Block), data=abamgermdat, family=binomial)# 
abammod1<- glm(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory, data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod1<- glmer(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+ (1|Block), data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod2<- glmer(abamy~Stand*Origin + (1|Block), data=abamgermdat, family=binomial)# still failed to converge
abammod2<- glm(abamy~Stand*Origin, data=abamgermdat, family=binomial)#
abam.mmod2a<- glmer(abamy~Canopy+Stand*Origin+ (1|Block), data=abamgermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
abammod2a<- glm(abamy~Canopy+Stand*Origin, data=abamgermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
abammod3<- glm(abamy~Canopy*Understory*Origin, data=abamgermdat, family=binomial)# 
abam.mmod3<- glmer(abamy~Canopy*Understory*Origin+(1|Block), data=abamgermdat, family=binomial)# 
abammod4<- glm(abamy~Canopy*Stand, data=abamgermdat, family=binomial)# 
abammod4a<- glm(abamy~Origin+Canopy*Stand, data=abamgermdat, family=binomial)# 
abam.mmod4<- glmer(abamy~Canopy*Stand+(1|Block), data=abamgermdat, family=binomial)# 
abam.mmod4a<- glmer(abamy~Origin+Canopy*Stand+(1|Block), data=abamgermdat, family=binomial)# 
abammod5<- glm(abamy~Understory*Stand, data=abamgermdat, family=binomial)# 
abam.mmod5<- glmer(abamy~Understory*Stand+(1|Block), data=abamgermdat, family=binomial)# 
abammod6<- glm(abamy~Comp*Stand, data=abamgermdat, family=binomial)# 
abam.mmod6<- glmer(abamy~Comp*Stand + (1|Block), data=abamgermdat, family=binomial)# 

AICctab(abammod1,abammod2,abammod3,abammod4,abammod5,abammod6,abammod2a,abammod4a,abammod1.moss,abammod2.moss,abammod3.moss,abammod4.moss,abammod5.moss,abammod6.moss,abammod2a.moss,abammod4a.moss,logLik=T)
AICctab(abam.mmod1,abam.mmod2,abam.mmod3,abam.mmod4,abam.mmod5,abam.mmod6,abam.mmod2a,abam.mmod4a,abam.mmod1.moss,abam.mmod2.moss,abam.mmod3.moss,abam.mmod4.moss,abam.mmod5.moss,abam.mmod6.moss,abam.mmod2a.moss,abam.mmod4a.moss,logLik=T)
#abammod2a.moss "wins"
Anova(abam.mmod2a.moss, type="III")
Anova(abammod2a.moss, type="III")
#TSME
#make sure that datasets are the same, regardless of explanatory variables
tsmegermdat<- tsmegermdat[apply(tsmegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)

tsmemod<- glm(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmod<- glmer(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsmemoda<- glm(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmoda<- glmer(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy + (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsmemod2<- glm(tsmey~Stand*Origin, data=tsmegermdat, family=binomial)# also checked for a 3way int
tsme.mmod2<- glmer(tsmey~Stand*Origin+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way int
tsmemod3<- glm(tsmey~Origin*Canopy, data=tsmegermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
tsme.mmod3<- glmer(tsmey~Canopy*Origin + (1|Block), data=tsmegermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
tsmemod4<- glm(tsmey~Origin*Understory, data=tsmegermdat, family=binomial)# 
tsme.mmod4<- glmer(tsmey~Origin*Understory+ (1|Block), data=tsmegermdat, family=binomial)# 
tsmemod5<- glm(tsmey~Stand*Canopy, data=tsmegermdat, family=binomial)# 
tsme.mmod5<- glmer(tsmey~Stand*Canopy+ (1|Block), data=tsmegermdat, family=binomial)# 
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)
#add moss
tsmemod.moss<- glm(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmod.moss<- glmer(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsmemoda.moss<- glm(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy+Moss, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmoda.moss<- glmer(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy+Moss + (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsmemod2.moss<- glm(tsmey~Stand*Origin+Moss, data=tsmegermdat, family=binomial)# also checked for a 3way int
tsme.mmod2.moss<- glmer(tsmey~Stand*Origin+Moss+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way int
tsmemod3.moss<- glm(tsmey~Origin*Canopy+Moss, data=tsmegermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
tsme.mmod3.moss<- glmer(tsmey~Canopy*Origin+Moss + (1|Block), data=tsmegermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
tsmemod4.moss<- glm(tsmey~Origin*Understory+Moss, data=tsmegermdat, family=binomial)# 
tsme.mmod4.moss<- glmer(tsmey~Origin*Understory+Moss+ (1|Block), data=tsmegermdat, family=binomial)# 
tsmemod5.moss<- glm(tsmey~Stand*Canopy+Moss, data=tsmegermdat, family=binomial)# 
tsme.mmod5.moss<- glmer(tsmey~Stand*Canopy+Moss+ (1|Block), data=tsmegermdat, family=binomial)# 

AICctab(tsmemod,tsmemoda,tsmemod2,tsmemod3,tsmemod4,tsmemod5,tsmemod.moss,tsmemoda.moss,tsmemod2.moss,tsmemod3.moss,tsmemod4.moss,tsmemod5.moss,logLik=T)
AICctab(tsme.mmod,tsme.mmoda,tsme.mmod2,tsme.mmod3,tsme.mmod4,tsme.mmod5,tsme.mmod.moss,tsme.mmoda.moss,tsme.mmod2.moss,tsme.mmod3.moss,tsme.mmod4.moss,tsme.mmod5.moss,logLik=T)

anova(tsme.mmod3,type="III")
#tsmemod3 has lowest AIC, so I'll use that model for the table
Anova(tsmemod3,test="LR", type="II")#without moss

#Tshe
tshegermdat<- tshegermdat[apply(tshegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

tshemod<- glm(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=tshegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig
tshe.mmod<- glmer(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory +(1|Block), data=tshegermdat, family=binomial, REML=FALSE)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig
tshemod2<- glm(tshey~Stand*Origin, data=tshegermdat,contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
tshe.mmod2<- glmer(tshey~Stand*Origin +(1|Block), data=tshegermdat, contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
tshemod2a<- glm(tshey~Stand*Origin+Canopy, data=tshegermdat, family=binomial)# 
tshe.mmod2a<- glmer(tshey~Stand*Origin+Canopy +(1|Block), data=tshegermdat, family=binomial, REML=FALSE)# 
tshemod3<- glm(tshey~Canopy*Understory*Origin, data=tshegermdat, family=binomial)# 
tshe.mmod3<- glmer(tshey~Canopy*Understory*Origin +(1|Block), data=tshegermdat, family=binomial)# 
tshemod4<- glm(tshey~Canopy*Stand, data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshe.mmod4<- glmer(tshey~Canopy*Stand +(1|Block), data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshemod5<- glm(tshey~Understory*Stand, data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshe.mmod5<- glmer(tshey~Understory*Stand +(1|Block), data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
##Add moss
tshemod.moss<- glm(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss, data=tshegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig
tshe.mmod.moss<- glmer(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss +(1|Block), data=tshegermdat, family=binomial, REML=FALSE)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig
tshemod2.moss<- glm(tshey~Stand*Origin+Moss, data=tshegermdat,contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
tshe.mmod2.moss<- glmer(tshey~Stand*Origin+Moss +(1|Block), data=tshegermdat, contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
tshemod2a.moss<- glm(tshey~Stand*Origin+Canopy+Moss, data=tshegermdat, family=binomial)# 
tshe.mmod2a.moss<- glmer(tshey~Stand*Origin+Canopy+Moss +(1|Block), data=tshegermdat, family=binomial, REML=FALSE)# 
tshemod3.moss<- glm(tshey~Canopy*Understory*Origin+Moss, data=tshegermdat, family=binomial)# 
tshe.mmod3.moss<- glmer(tshey~Canopy*Understory*Origin+Moss +(1|Block), data=tshegermdat, family=binomial)# 
tshemod4.moss<- glm(tshey~Canopy*Stand+Moss, data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshe.mmod4.moss<- glmer(tshey~Canopy*Stand+Moss +(1|Block), data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshemod5.moss<- glm(tshey~Understory*Stand+Moss, data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-
tshe.mmod5.moss<- glmer(tshey~Understory*Stand+Moss +(1|Block), data=tshegermdat, family=binomial)# couldn't fit model with both understory and canopy-

AICctab(tshemod,tshemod2,tshemod3,tshemod4,tshemod5,tshemod2a,tshemod.moss,tshemod2.moss,tshemod3.moss,tshemod4.moss,tshemod5.moss,tshemod2a.moss,logLik=T)
AICctab(tshe.mmod,tshe.mmod2,tshe.mmod3,tshe.mmod4,tshe.mmod5,tshe.mmod2a,tshe.mmod.moss,tshe.mmod2.moss,tshe.mmod3.moss,tshe.mmod4.moss,tshe.mmod5.moss,tshe.mmod2a.moss,logLik=T)
#mod2.moss is best fit for tshe
Anova(tshe.mmod2.moss, test="Chisq",type="III")
Anova(tshemod2.moss,test="LR", type="III")
