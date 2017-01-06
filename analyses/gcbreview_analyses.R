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
dim(tsmedat)#1500 rows
dim(tshedat)#960 rows
dim(abamdat)#1499 rows
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
##Figure 1 (Variance explained): use the above models for figure adn table
propdev.tsmegerm=read.csv("analyses/tsmegerm.propdev.csv",header=TRUE)
propdev.abamgerm=read.csv("analyses/abamgerm.propdev.csv",header=TRUE)
propdev.tshegerm=read.csv("analyses/tshegerm.propdev.csv",header=TRUE)
propdev.tsmesurv=cbind(rownames(anova(constmod.tsme)),round(anova(constmod.tsme)$Dev/(anova(constmod.tsme)$"-2*LL"[1]-anova(constmod.tsme)$"-2*LL"[12]),digits=3))
propdev.tshesurv=cbind(rownames(anova(constmod.tshe)),round(anova(constmod.tshe)$Dev/(anova(constmod.tshe)$"-2*LL"[1]-anova(constmod.tshe)$"-2*LL"[12]),digits=3))
propdev.abamsurv=cbind(rownames(anova(constmod.abam)),round(anova(constmod.abam)$Dev/(anova(constmod.abam)$"-2*LL"[1]-anova(constmod.abam)$"-2*LL"[12]),digits=3))
propdev.tsmehi=cbind(rownames(anova(consthimod.tsme)),round(anova(consthimod.tsme)$"Sum Sq"/sum(anova(consthimod.tsme)$"Sum Sq"),digits=3))
propdev.tshehi=cbind(rownames(anova(consthimod.tshe)),round(anova(consthimod.tshe)$"Sum Sq"/sum(anova(consthimod.tshe)$"Sum Sq"),digits=3))
propdev.abamhi=cbind(rownames(anova(consthimod.abam)),round(anova(consthimod.abam)$"Sum Sq"/sum(anova(consthimod.abam)$"Sum Sq"),digits=3))

propdev.tsmeall=cbind(propdev.tsmegerm,propdev.tsmesurv[2:12,2],propdev.tsmehi[,2])
propdev.tsmeall<-propdev.tsmeall[,-2]
colnames(propdev.tsmeall)=c("Predictor","Germination","Survival","Growth")
propdev.tsmeall[,2]=as.numeric(propdev.tsmeall[,2])
propdev.tsmeall[,3]=as.numeric(propdev.tsmeall[,3])
propdev.tsmeall[,4]=as.numeric(propdev.tsmeall[,4])
#ABAM
#sort germ to be in same order as surv and growth
ord=c(2,3,1,4,5,6,7,8,9,10,11)
propdev.abamgerm <- propdev.abamgerm[order(ord),] 
propdev.abamall=cbind(propdev.abamgerm,propdev.abamsurv[2:12,2],propdev.abamhi[,2])
propdev.abamall<-propdev.abamall[,-2]
colnames(propdev.abamall)=c("Predictor","Germination","Survival","Growth")
propdev.abamall[,2]=as.numeric(propdev.abamall[,2])
propdev.abamall[,3]=as.numeric(propdev.abamall[,3])
propdev.abamall[,4]=as.numeric(propdev.abamall[,4])
#now tshe
#sort germ to be in same order as surv and growth
ord=c(2,3,1,4,5,6,7,8,9,10,11)
propdev.tshegerm <- propdev.tshegerm[order(ord),] 
propdev.tsheall=cbind(propdev.tshegerm,propdev.tshesurv[2:12,2],propdev.tshehi[,2])
propdev.tsheall<-propdev.tsheall[,-2]
colnames(propdev.tsheall)=c("Predictor","Germination","Survival","Growth")
propdev.tsheall[,2]=as.numeric(propdev.tsheall[,2])
propdev.tsheall[,3]=as.numeric(propdev.tsheall[,3])
propdev.tsheall[,4]=as.numeric(propdev.tsheall[,4])

ord=c(1,5,2,3,9,6,7,10,11,4,8)
propdev.tsmeall2 <- propdev.tsmeall[order(ord),] 
propdev.tsheall2 <- propdev.tsheall[order(ord),] 
propdev.abamall2 <- propdev.abamall[order(ord),] 
#choose colors
brewer.pal(9,"Paired")
#quartz(height=6.5,width=4)
quartz(height=6.5,width=5)
par(mfrow=c(3,1),mai=c(.6,.6,.2,.6), omi=c(.7,.1,.2,.2))
#par(mfrow=c(3,1),mai=c(.6,.6,.1,.1))
#x=barplot(cbind(propdev.tsmeall[1:9,],c(rep(NA,times=9))), ylab=" ", col=c("darkblue", "darkred","greenyellow","greenyellow","purple4","turquoise4","turquoise4","darorange3","darkorange3"),cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="Vital rate", space=c(.5,.5,.5,1), ylim=c(0,1))
#x=barplot(cbind(propdev.tsmeall[1:9,],c(rep(NA,times=9))), ylab="", col=c("darkblue", "darkred","gold1","gold1","purple4","darkgreen","darkgreen","darkorange2","darkorange2"),cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="", space=c(.5,.5,.5,1), ylim=c(0,1))
colors=c("darkblue","greenyellow","greenyellow", "darkred","darkgreen","darkgreen","purple4","orange4","orange4")
colors2=c("#1F78B4","#33A02C","#33A02C","#33A02C","#E31A1C","#B2DF8A","#B2DF8A","#B2DF8A","#CAB2D6","#FDBF6F","#FDBF6F")
x=barplot(cbind(as.matrix(propdev.tsmeall2[1:11,2:4]),c(1/6,1/6,0,0,1/6,1/6,0,0,1/6,1/6,0)), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="", space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15), ylim=c(0,1), las=1)
#labeloc=c(c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+(as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+(as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])+(+as.numeric(propdev.tsmeall[8,2])+as.numeric(propdev.tsmeall[9,2])/2)))
#text(x[1],c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])/2),labels=propdev.tsmeall[1:9,1],pos=1)
#text(x[1],c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+(as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+(as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])+(+as.numeric(propdev.tsmeall[8,2])+as.numeric(propdev.tsmeall[9,2])/2)),labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),pos=1)
labelsloc=c(1/12,3/12,5/12,7/12,9/12,11/12)
text(2.5,labelsloc,labels=c("Elevation","Competition","Origin","Elevation*Competition","Elevation*Origin","Origin*Competition"),cex=.9, adj=0)
mtext("Tsuga mertensiana",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("a",side=3,line=.5, adj=-.1, cex=.9)
#ABAM
x=barplot(cbind(as.matrix(propdev.abamall2[1:11,2:4]),c(rep(NA,times=11))), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="", space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15),  ylim=c(0,1), las=1)
#labelsloc=c(.2,.351,.48,.63,.83,.97)
#text(2.3,labelsloc,labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),cex=.9, adj=1)
mtext("Abies amabilis",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("b",side=3,line=.5, adj=-.1, cex=.9)
mtext("Proportion variance explained", side=2,line=3,cex=.9)
#TSHE
x=barplot(cbind(as.matrix(propdev.tsheall2[1:11,2:4]),c(rep(NA,times=11))), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="Vital rate",space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15), ylim=c(0,1), las=1)
#labelsloc=c(.2,.351,.48,.63,.83,.97)
#text(2.3,labelsloc,labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),cex=.9, adj=1)
mtext("Tsuga heterophylla",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("c",side=3,line=.5, adj=-.1, cex=.9)

#Now add microclimate data, including GDD, snowduration, and mean light for each block
climdat<-read.csv("analyses/clim_gf.csv", header=TRUE)
#Combine the climate data with the transplant data columns that we want to use for survival and growth analysis
transdatsub<-subset(transdat,select=c("UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi"))
climdat$Canopy<-"CompAbsent"
climdat[which(climdat$canopy=="N"),]$Canopy<-"CompPresent"
climdat$Understory<-"CompAbsent"
climdat[which(climdat$understory=="C"),]$Understory<-"CompPresent"
#take average across both years?
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

#Next step is to modify the models so that the explanatory variables are GDD and light instead of range position/treatment and 
#Previous models:
#in previous preliminary analyses (available upon request), I used model selection to identify that the lognormal distirbution is best-fit for these data
#make sure that datasets are the same, regardless of explanatory variables
#"UniqueID","Species","Block","PlantedStand","OriginStand","Canopy","Understory","time1","time2","annhi","Canopy","Understory","GDD_total","GDD_totaln","meanGST","Light_Mean","Light_GDD","MAT","snow_dur","snow_dur_cont"
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
Anova(constmod.tsme, test.statistic="LR", type="III")
Anova(constmod.abam, test.statistic="LR", type="III")
Anova(constmod.tshe,test.statistic="LR", type="III")
#For new models, replace categorical competition and elevation variables with contuous microcliamte variables that we measured
#climate models
matmod.abam<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=abamsurvdat)
matmod.tsme<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
matmod.tshe<-survreg(Surv(time1,time2, type="interval2")~MAT_cent+MAT_cent:OriginStand, dist="lognormal", data=tshesurvdat)

gddmod.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gddmod.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gddmod.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_total_cent+GDD_total_cent:OriginStand, dist="lognormal", data=tshesurvdat)
gddmod2.abam<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gddmod2.tsme<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gddmod2.tshe<-survreg(Surv(time1,time2, type="interval2")~GDD_totaln_cent+GDD_totaln_cent:OriginStand, dist="lognormal", data=tshesurvdat)

gstmod.abam<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=abamsurvdat)
gstmod.tsme<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=tsmesurvdat)
gstmod.tshe<-survreg(Surv(time1,time2, type="interval2")~meanGST_cent+meanGST_cent:OriginStand, dist="lognormal", data=tshesurvdat)
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
nullmod.abam<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=abamsurvdat)
nullmod.tsme<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=tsmesurvdat)
nullmod.tshe<-survreg(Surv(time1,time2, type="interval2")~1, dist="lognormal", data=tshesurvdat)

AICtab(nullmod.abam,constmod.abam, matmod.abam,gddmod.abam,gddmod2.abam,sdmod.abam,sdmod2.abam,gstmod.abam,lightmod.abam,lightmod2.abam,constgddlightmod.abam,constgddlight2mod.abam,constgdd2lightmod.abam,constgdd2light2mod.abam,constgstlightmod.abam,constgstlight2mod.abam,constsdlightmod.abam,constsd2lightmod.abam,constsdlight2mod.abam,constsd2light2mod.abam,constmatlightmod.abam,constmatlight2mod.abam)#lowest aic for constmod.abam- 
AICtab(nullmod.tsme,constmod.tsme, matmod.tsme,gddmod.tsme,gddmod2.tsme,sdmod.tsme,sdmod2.tsme,gstmod.tsme,lightmod.tsme,lightmod2.tsme,constgddlightmod.tsme,constgddlight2mod.tsme,constgdd2lightmod.tsme,constgdd2light2mod.tsme,constgstlightmod.tsme,constgstlight2mod.tsme,constsdlightmod.tsme,constsd2lightmod.tsme,constsdlight2mod.tsme,constsd2light2mod.tsme,constmatlightmod.tsme,constmatlight2mod.tsme)#lowest aic for constmod.tsme- #
AICtab(nullmod.tshe,constmod.tshe, matmod.tshe,gddmod.tshe,gddmod2.tshe,sdmod.tshe,sdmod2.tshe,gstmod.tshe,lightmod.tshe,lightmod2.tshe,constgddlightmod.tshe,constgddlight2mod.tshe,constgdd2lightmod.tshe,constgdd2light2mod.tshe,constgstlightmod.tshe,constgstlight2mod.tshe,constsdlightmod.tshe,constsd2lightmod.tshe,constsdlight2mod.tshe,constsd2light2mod.tshe,constmatlightmod.tshe,constmatlight2mod.tshe)#lowest aic for constmod.tshe- 

#AICc(constmod.tshe)-AICc(sdmod.tshe);AICc(gddmod.tshe);AICc(sdmod.tshe);AICc(gstmod.tshe);AICc(lightmod.tshe);AICc(constgddlightmod.tshe);AICc(constgstlightmod.tshe);AICc(constsdlightmod.tshe)

#Now growth data
#make sure that datasets are the same, regardless of explanatory variables
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

#GDD
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


AICtab(nullhimod.abam,consthimod.abam, mathimod.abam,gddhimod.abam,gddhimod2.abam,sdhimod.abam,sdhimod2.abam,gsthimod.abam,lighthimod.abam,lighthimod2.abam,constgddlighthimod.abam,constgddlight2himod.abam,constgdd2lighthimod.abam,constgdd2light2himod.abam,constgstlighthimod.abam,constgstlight2himod.abam,constsdlighthimod.abam,constsd2lighthimod.abam,constsdlight2himod.abam,constsd2light2himod.abam,constmatlighthimod.abam,constmatlight2himod.abam)#lowest aic for consthihimod.abam- 
AICtab(nullhimod.tsme,consthimod.tsme, mathimod.tsme,gddhimod.tsme,gddhimod2.tsme,sdhimod.tsme,sdhimod2.tsme,gsthimod.tsme,lighthimod.tsme,lighthimod2.tsme,constgddlighthimod.tsme,constgddlight2himod.tsme,constgdd2lighthimod.tsme,constgdd2light2himod.tsme,constgstlighthimod.tsme,constgstlight2himod.tsme,constsdlighthimod.tsme,constsd2lighthimod.tsme,constsdlight2himod.tsme,constsd2light2himod.tsme,constmatlighthimod.tsme,constmatlight2himod.tsme)#lowest aic for consthihimod.tsme- #
AICtab(nullhimod.tshe,consthimod.tshe, mathimod.tshe,gddhimod.tshe,gddhimod2.tshe,sdhimod.tshe,sdhimod2.tshe,gsthimod.tshe,lighthimod.tshe,lighthimod2.tshe,constgddlighthimod.tshe,constgddlight2himod.tshe,constgdd2lighthimod.tshe,constgdd2light2himod.tshe,constgstlighthimod.tshe,constgstlight2himod.tshe,constsdlighthimod.tshe,constsd2lighthimod.tshe,constsdlight2himod.tshe,constsd2light2himod.tshe,constmatlighthimod.tshe,constmatlight2himod.tshe)#lowest aic for consthihimod.tshe- 
#Get coefficients from best-fit models

summary(constgdd2light2himod.tshe)
summary(constgdd2lighthimod.tsme)
summary(constgdd2light2himod.abam)
#Try to figure out why the coefficient is negative
quartz()
plot(tsmegrdat$annhi~tsmegrdat$meanGST_cent)
abline(coef(constgstlighthimod.tsme))
tapply(tsmegrdat$annhi,list(tsmegrdat$PlantedStand,tsmegrdat$Canopy,tsmegrdat$Understory), mean)
tapply(tsmegrdat$meanGST,list(tsmegrdat$PlantedStand,tsmegrdat$Canopy,tsmegrdat$Understory), mean)

#remake Figure 4 adding germination
#Figure 4.plot, with change in effects of comp vs no comp
tshegermdat$Elev<-as.factor(tshegermdat$Elev)
tsmegermdat$Elev<-as.factor(tsmegermdat$Elev)
abamgermdat$Elev<-as.factor(abamgermdat$Elev)

compsurvmod.tsme=survreg(Surv(time1,time2, type="interval2")~-1+PlantedStand*Canopy*Understory, dist="lognormal", data=tsmedat)
summary(compsurvmod.tsme)#LL: -2813.5 (Null LL:-2864.7)
compsurvmod.abam=survreg(Surv(time1,time2, type="interval2")~-1+PlantedStand*Canopy*Understory, dist="lognormal", data=abamdat)
summary(compsurvmod.abam)#LL: -2616.3 (Null LL:-2765.6)
compsurvmod.tshe=survreg(Surv(time1,time2, type="interval2")~-1+PlantedStand*Canopy*Understory, dist="lognormal", data=tshedat)
summary(compsurvmod.tshe)#LL: -1422.8 (Null LL:-1546.4)
himod.abam<-lmer(annhi ~ -1+PlantedStand*Canopy*Understory+(1|Block),REML=FALSE, data=abamdat)#
summary(himod.abam)
himod.tsme<-lmer(annhi ~ -1+PlantedStand*Canopy*Understory+(1|Block),REML=FALSE, data=tsmedat)
summary(himod.tsme)
himod.tshe<-lmer(annhi ~ -1+PlantedStand*Canopy*Understory+(1|Block),REML=FALSE, data=tshedat)#
summary(himod.tshe)
abamgermdat$Elev<-as.factor(abamgermdat$Elev)
germod.abam<-glmmadmb(abamy~-1+Elev*Canopy*Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
germod.tshe<-glmmadmb(tshey~-1+Elev*Canopy*Understory+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
germod.tsme<-glmmadmb(tsmey~-1+Elev*Canopy*Understory+(1|Block), data=tsmegermdat, family="binomial")

#get coefs and their ses
#germination
abamgermcan=c(coef(germod.abam)[6],coef(germod.abam)[6]+coef(germod.abam)[9:12])
abamgermund=c(coef(germod.abam)[7],coef(germod.abam)[7]+coef(germod.abam)[13:16])
abamgermcan.se=c(summary(germod.abam)$coef[6,2],summary(germod.abam)$coef[9:12,2])
abamgermund.se=c(summary(germod.abam)$coef[7,2],summary(germod.abam)$coef[13:16,2])
tshegermcan=c(coef(germod.tshe)[4],coef(germod.tshe)[4]+coef(germod.tshe)[7:8])
tshegermund=c(coef(germod.tshe)[5],coef(germod.tshe)[5]+coef(germod.tshe)[9:10])
tshegermcan.se=c(summary(germod.tshe)$coef[6,2],summary(germod.tshe)$coef[7:8,2])
tshegermund.se=c(summary(germod.tshe)$coef[7,2],summary(germod.tshe)$coef[9:10,2])
tsmegermcan=c(coef(germod.tsme)[6],coef(germod.tsme)[6]+coef(germod.tsme)[8:11])
tsmegermund=c(coef(germod.tsme)[7],coef(germod.tsme)[7]+coef(germod.tsme)[12:15])
tsmegermcan.se=c(summary(germod.tsme)$coef[6,2],summary(germod.tsme)$coef[8:11,2])
tsmegermund.se=c(summary(germod.tsme)$coef[7,2],summary(germod.tsme)$coef[12:15,2])

#survival
tsmesurvcan=c(coef(compsurvmod.tsme)[6],coef(compsurvmod.tsme)[6]+coef(compsurvmod.tsme)[8:11])
tsmesurvund=c(coef(compsurvmod.tsme)[7],coef(compsurvmod.tsme)[7]+coef(compsurvmod.tsme)[12:15])
tsmesurvcan.se=c(summary(compsurvmod.tsme)$table[6,2],summary(compsurvmod.tsme)$table[8:11,2])
tsmesurvund.se=c(summary(compsurvmod.tsme)$table[7,2],summary(compsurvmod.tsme)$table[12:15,2])
abamsurvcan=c(coef(compsurvmod.abam)[6],coef(compsurvmod.abam)[6]+coef(compsurvmod.abam)[8:11])
abamsurvund=c(coef(compsurvmod.abam)[7],coef(compsurvmod.abam)[7]+coef(compsurvmod.abam)[12:15])
abamsurvcan.se=c(summary(compsurvmod.abam)$table[6,2],summary(compsurvmod.abam)$table[8:11,2])
abamsurvund.se=c(summary(compsurvmod.abam)$table[7,2],summary(compsurvmod.abam)$table[12:15,2])
tshesurvcan=c(coef(compsurvmod.tshe)[4],coef(compsurvmod.tshe)[4]+coef(compsurvmod.tshe)[6:7])
tshesurvund=c(coef(compsurvmod.tshe)[5],coef(compsurvmod.tshe)[5]+coef(compsurvmod.tshe)[8:9])
tshesurvcan.se=c(summary(compsurvmod.tshe)$table[4,2],summary(compsurvmod.tshe)$table[6:7,2])
tshesurvund.se=c(summary(compsurvmod.tshe)$table[4,2],summary(compsurvmod.tshe)$table[8:9,2])
#HI
tsmehican=c(fixef(himod.tsme)[6],fixef(himod.tsme)[6]+fixef(himod.tsme)[8:11])
tsmehiund=c(fixef(himod.tsme)[7],fixef(himod.tsme)[7]+fixef(himod.tsme)[12:15])
tsmehican.se=c(summary(himod.tsme)$coef[6,2],summary(himod.tsme)$coef[8:11,2])
tsmehiund.se=c(summary(himod.tsme)$coef[7,2],summary(himod.tsme)$coef[12:15,2])
abamhican=c(fixef(himod.abam)[6],fixef(himod.abam)[6]+fixef(himod.abam)[8:11])
abamhiund=c(fixef(himod.abam)[7],fixef(himod.abam)[7]+fixef(himod.abam)[12:15])
abamhican.se=c(summary(himod.abam)$coef[6,2],summary(himod.abam)$coef[8:11,2])
abamhiund.se=c(summary(himod.abam)$coef[7,2],summary(himod.abam)$coef[12:15,2])
tshehican=c(fixef(himod.tshe)[4],fixef(himod.tshe)[4]+fixef(himod.tshe)[6:7])
tshehiund=c(fixef(himod.tshe)[5],fixef(himod.tshe)[5]+fixef(himod.tshe)[8:9])
tshehican.se=c(summary(himod.tshe)$coef[4,2],summary(himod.tshe)$coef[6:7,2])
tshehiund.se=c(summary(himod.tshe)$coef[4,2],summary(himod.tshe)$coef[8:9,2])

###the barplot
quartz(height=7,width=12)
par(mfcol=c(3,3),mai=c(.6,.8,.2,.1), omi=c(.7,.01,.2,.2))
#Germination
plottsme<-barplot(as.matrix(rbind(tsmegermund,tsmegermcan)),width=.9,ylab="",xlab="",names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("a",side=3,line=1, adj=-.1)
polygon(c(3.5,3.5,11,11),c(-1,1,1,-1),col="light gray", border="light gray")
par(new=T)
plottsme<-barplot(as.matrix(rbind(tsmegermund,tsmegermcan)),width=.9,ylab="",xlab="",names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottsme)
alltsme.lcl=rbind(tsmegermund-tsmegermund.se,tsmegermcan-tsmegermcan.se)
alltsme.ucl=rbind(tsmegermund+tsmegermund.se,tsmegermcan+tsmegermcan.se)
#add error bars
for (i in 1:length(alltsme.lcl)){
  arrows(x[i],alltsme.lcl[i],x[i],alltsme.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
legend(.3,1.1,legend=c("Understory", "Canopy"),bty="n",fill=c("palegreen1","darkgreen"), cex=1.2)
#ABAM
plotabam<-barplot(as.matrix(rbind(abamgermund,abamgermcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),xaxt='n', yaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.5,1.5,1.5,-1.5),col="light gray", border="light gray")
par(new=T)
plotabam<-barplot(as.matrix(rbind(abamgermund,abamgermcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("b",side=3,line=1, adj=-.1)
abline(h=0)
mtext("Effect of neighbors on germination", side=2, line=4.3)
mtext("(Difference in germination, relative to no neighbors, on a logit scale)", side=2, line=3, cex=.9)
x<-c(plotabam)
allabam.lcl=rbind(abamgermund-abamgermund.se,abamgermcan-abamgermcan.se)
allabam.ucl=rbind(abamgermund+abamgermund.se,abamgermcan+abamgermcan.se)
#add error bars
for (i in 1:length(allabam.lcl)){
  arrows(x[i],allabam.lcl[i],x[i],allabam.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(668)", "(704)","(1064)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##TSHE
tshegermund2=c(NA,NA,tshegermund)
tshegermcan2=c(NA,NA,tshegermcan)
plottshe<-barplot(as.matrix(rbind(tshegermund2,tshegermcan2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),xaxt='n', yaxt='n',ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.5,.5,.5,-1.5),col="light gray", border="light gray")
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("c",side=3,line=1, adj=-.1)
abline(h=0)
par(new=T)
plottshe<-barplot(as.matrix(rbind(tshegermund2,tshegermcan2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottshe)
tshegermund.se2=c(0,0,tshegermund.se)
tshegermcan.se2=c(0,0,tshegermcan.se)
alltshe.lcl2=rbind(tshegermund2-tshegermund.se2,tshegermcan2-tshegermcan.se2)
alltshe.ucl2=rbind(tshegermund2+tshegermund.se2,tshegermcan2+tshegermcan.se2)
#add error bars
for (i in 1:length(alltshe.lcl2)){
  arrows(x[i],alltshe.lcl2[i],x[i],alltshe.ucl2[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext("Location in range (m)",line=-13.5, adj=.6)
#Survival
#TSME
plottsme<-barplot(as.matrix(rbind(tsmesurvund,tsmesurvcan)),width=.9,ylab="",xlab="",names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("d",side=3,line=1, adj=-.1)
polygon(c(3.5,3.5,11,11),c(-1,1,1,-1),col="light gray", border="light gray")
par(new=T)
plottsme<-barplot(as.matrix(rbind(tsmesurvund,tsmesurvcan)),width=.9,ylab="",xlab="",names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottsme)
alltsme.lcl=rbind(tsmesurvund-tsmesurvund.se,tsmesurvcan-tsmesurvcan.se)
alltsme.ucl=rbind(tsmesurvund+tsmesurvund.se,tsmesurvcan+tsmesurvcan.se)
#add error bars
for (i in 1:length(alltsme.lcl)){
  arrows(x[i],alltsme.lcl[i],x[i],alltsme.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

#ABAM
plotabam<-barplot(as.matrix(rbind(abamsurvund,abamsurvcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.5,1.5),las=1,col=c("palegreen1","darkgreen"),xaxt='n', yaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.5,1.5,1.5,-1.5),col="light gray", border="light gray")
par(new=T)
plotabam<-barplot(as.matrix(rbind(abamsurvund,abamsurvcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.5,1.5),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("e",side=3,line=1, adj=-.1)
abline(h=0)
mtext("Effect of neighbors on survival time", side=2, line=4.3)
mtext("(Proportion, relative to no neighbors)", side=2, line=3, cex=.9)
x<-c(plotabam)
allabam.lcl=rbind(abamsurvund-abamsurvund.se,abamsurvcan-abamsurvcan.se)
allabam.ucl=rbind(abamsurvund+abamsurvund.se,abamsurvcan+abamsurvcan.se)
#add error bars
for (i in 1:length(allabam.lcl)){
  arrows(x[i],allabam.lcl[i],x[i],allabam.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(668)", "(704)","(1064)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##TSHE
tshesurvund2=c(NA,NA,tshesurvund)
tshesurvcan2=c(NA,NA,tshesurvcan)
plottshe<-barplot(as.matrix(rbind(tshesurvund2,tshesurvcan2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),xaxt='n', yaxt='n',ylim=c(-1.5,.5),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.5,.5,.5,-1.5),col="light gray", border="light gray")
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("f",side=3,line=1, adj=-.1)
abline(h=0)
par(new=T)
plottshe<-barplot(as.matrix(rbind(tshesurvund2,tshesurvcan2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.5,.5),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottshe)
tshesurvund.se2=c(0,0,tshesurvund.se)
tshesurvcan.se2=c(0,0,tshesurvcan.se)
alltshe.lcl2=rbind(tshesurvund2-tshesurvund.se2,tshesurvcan2-tshesurvcan.se2)
alltshe.ucl2=rbind(tshesurvund2+tshesurvund.se2,tshesurvcan2+tshesurvcan.se2)
#add error bars
for (i in 1:length(alltshe.lcl2)){
  arrows(x[i],alltshe.lcl2[i],x[i],alltshe.ucl2[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext("Location in range (m)",line=-13.5, adj=.6)

#Now ANNUAL HI
plottsme<-barplot(as.matrix(rbind(tsmehiund,tsmehican)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("g",side=3,line=1, adj=-.1)
polygon(c(3.5,3.5,11,11),c(-1,1,1,-1),col="light gray", border="light gray")
par(new=T)
plottsme<-barplot(as.matrix(rbind(tsmehiund,tsmehican)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottsme)
allhitsme.lcl=rbind(tsmehiund-tsmehiund.se,tsmehican-tsmehican.se)
allhitsme.ucl=rbind(tsmehiund+tsmehiund.se,tsmehican+tsmehican.se)
#add error bars
for (i in 1:length(allhitsme.lcl)){
  arrows(x[i],allhitsme.lcl[i],x[i],allhitsme.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##ABAM
plotabam<-barplot(as.matrix(rbind(abamhiund,abamhican)),ylab="",xlab="",names.arg=c("","","","",""),width=.9,ylim=c(-1.2,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
polygon(c(3.5,3.5,11,11),c(-1.2,1,1,-1.2),col="light gray", border="light gray")
par(new=T)
plotabam<-barplot(as.matrix(rbind(abamhiund,abamhican)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.2,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plotabam)
allhiabam.lcl=rbind(abamhiund-abamhiund.se,abamhican-abamhican.se)
allhiabam.ucl=rbind(abamhiund+abamhiund.se,abamhican+abamhican.se)
#668,704,1064,1197,1460,1605,1676
#add error bars
for (i in 1:length(allhiabam.lcl)){
  arrows(x[i],allhiabam.lcl[i],x[i],allhiabam.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(668)","(704)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("h",side=3,line=1, adj=-.1)
mtext("Effect of neighbors on growth", side=2, line=4.3)
mtext("(Difference in annual increment, relative to no neighbors)", side=2, line=3, cex=.9)
##TSHE
tshehiund2=c(NA,NA,tshehiund)
tshehican2=c(NA,NA,tshehican)
plottshe<-barplot(as.matrix(rbind(tshehiund2,tshehican2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),xaxt='n', yaxt='n',ylim=c(-1.7,.5),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.7,.5,.5,-1.7),col="light gray", border="light gray")
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("i",side=3,line=1, adj=-.1)
abline(h=0)
par(new=T)
plottshe<-barplot(as.matrix(rbind(tshehiund2,tshehican2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.7,.5),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
abline(h=0)
x<-c(plottshe)
tshehiund.se2=c(0,0,tshehiund.se)
tshehican.se2=c(0,0,tshehican.se)
alltshe.lcl2=rbind(tshehiund2-tshehiund.se2,tshehican2-tshehican.se2)
alltshe.ucl2=rbind(tshehiund2+tshehiund.se2,tshehican2+tshehican.se2)
#add error bars
for (i in 1:length(alltshe.lcl2)){
  arrows(x[i],alltshe.lcl2[i],x[i],alltshe.ucl2[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext("Location in range (m)",line=-13.5, adj=.6)

