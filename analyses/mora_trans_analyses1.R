######Code for Ettinger & HilleRisLambers "Assymmetric range shift dynamics" 
#####Submitted to GCB June 2016
#####Code by Ailene Ettinger, Ailene.Ettinger@tufts.edu; started 
######Germination, Survival and growth Analyses of 2013 Mt Rainier Transplants (survival/growth)  
#data from full experiment at the end of the 2013 growing season, at which point all seedlings were removed
setwd("~/GitHub/mora_transplant")
library(car)
library(lme4)
library(boot)
library(survival)
library(RColorBrewer)
library(KMsurv)
library(interval)
update.packages()
citation(package = "lme4", lib.loc = NULL, auto = NULL)
citation(package = "survival", lib.loc = NULL, auto = NULL)

###Germination 
dat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms=50
##First, select only rows in which seeds were added
addat=dat[dat$SeedsAdded>0,]
#select out 2011 data for analysis
dat2011<-addat[addat$Year=="2011",]
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
tsmegermdat$TotalGerms=factor(tsmegermdat$TotalGerms)
tshegermdat$TotalGerms=factor(tshegermdat$TotalGerms)
abamgermdat$TotalGerms=factor(abamgermdat$TotalGerms)
tsmegermdat$Block=factor(tsmegermdat$Block)
tshegermdat$Block=factor(tshegermdat$Block)
abamgermdat$Block=factor(abamgermdat$Block)
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)
##Using best-fit model (code for model selection for germination available upon request)
#TSME
tsme.mmod3<- glmer(tsmey~Origin*Canopy+(1|Block), data=tsmegermdat, family=binomial,contrasts=c(unordered="contr.sum", ordered="contr.poly"))# 
Anova(tsme.mmod3, type="III")
propdev.tsmegerm=cbind(rownames(anova(tsme.mmod3)),round(anova(tsme.mmod3)$"Sum Sq"/sum(anova(tsme.mmod3)$"Sum Sq"),digits=3))
#ABAM
abam.mmod2a<- glmer(abamy~Canopy+Stand*Origin+(1|Block), data=abamgermdat, contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
Anova(abam.mmod2a,type="III")
propdev.abamgerm=cbind(rownames(anova(abam.mmod2a)),round(anova(abam.mmod2a)$"Sum Sq"/sum(anova(abam.mmod2a)$"Sum Sq"),digits=3))
#TSHE
tshe.mmod2<- glmer(tshey~Stand*Origin +(1|Block), data=tshegermdat, contrasts=c(unordered="contr.sum", ordered="contr.poly"),family=binomial)# 
Anova(tshe.mmod2, test="Chisq",type="III")
propdev.tshegerm=cbind(rownames(anova(tshe.mmod2)),round(anova(tshe.mmod2)$"Sum Sq"/sum(anova(tshe.mmod2)$"Sum Sq"),digits=3))

##Survival and Growth data from transplants
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
#separate out species
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
####Look at root to crown height (RtCrnHeightDate5, measured only on the last census to avoid soil disturbance) and height from soil surface to apical bud tip on the same date (HeightDate5; this is the way that height was measured on alll other censuses)
cor(transdat$HeightDate5,transdat$RtCrnHeightDate5, use="pairwise.complete.obs")#0.867
summary(lm(transdat$HeightDate5~transdat$RtCrnHeightDate5))#p<0.001
#Check for elevational patterns in difference between these two heights
transdat$heightdif=transdat$RtCrnHeightDate5-transdat$HeightDate5
summary(lmer(heightdif~-1+PlantedStand + (1|Block), data=transdat))
mean(transdat$heightdif, na.rm=T)
sd(transdat$heightdif, na.rm=T)
###############################Species-specific survival regression models########## 
#in previous preliminary analyses (available upon request), I used model selection to identify that the lognormal distirbution is best-fit for these data
#ABAM
constmod.abam<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=abamdat)
#TSME
constmod.tsme<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tsmedat)
#TSHE
constmod.tshe<-survreg(Surv(time1,time2, type="interval2")~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory, dist="lognormal", data=tshedat)
Anova(constmod.tsme, test.statistic="LR", type="III")
Anova(constmod.abam, test.statistic="LR", type="III")
Anova(constmod.tshe,test.statistic="LR", type="III")

####Species-specific mixed models for annual height increment
#ABAM
consthimod.abam<-lmer(annhi ~ PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block),REML=FALSE, data=abamdat,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
consthimod.tsme<-lmer(annhi~ PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block),REML=FALSE, data=tsmedat,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
consthimod.tshe<-lmer(annhi~PlantedStand+OriginStand+Canopy+Understory+PlantedStand:OriginStand+PlantedStand:Canopy+PlantedStand:Understory+OriginStand:Canopy+OriginStand:Understory+Canopy:Understory+PlantedStand:Canopy:Understory+(1|Block),REML=FALSE, data=tshedat,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(consthimod.tsme, type="III")
Anova(consthimod.abam, type="III")
Anova(consthimod.tshe,type="III")
#Proportion variance explained for each vital rate
propdev.tsmesurv=cbind(rownames(anova(constmod.tsme)),round(anova(constmod.tsme)$Dev/(anova(constmod.tsme)$"-2*LL"[1]-anova(constmod.tsme)$"-2*LL"[12]),digits=3))
propdev.tshesurv=cbind(rownames(anova(constmod.tshe)),round(anova(constmod.tshe)$Dev/(anova(constmod.tshe)$"-2*LL"[1]-anova(constmod.tshe)$"-2*LL"[12]),digits=3))
propdev.abamsurv=cbind(rownames(anova(constmod.abam)),round(anova(constmod.abam)$Dev/(anova(constmod.abam)$"-2*LL"[1]-anova(constmod.abam)$"-2*LL"[12]),digits=3))
propdev.tsmehi=cbind(rownames(anova(consthimod.tsme)),round(anova(consthimod.tsme)$"Sum Sq"/sum(anova(consthimod.tsme)$"Sum Sq"),digits=3))
propdev.tshehi=cbind(rownames(anova(consthimod.tshe)),round(anova(consthimod.tshe)$"Sum Sq"/sum(anova(consthimod.tshe)$"Sum Sq"),digits=3))
propdev.abamhi=cbind(rownames(anova(consthimod.abam)),round(anova(consthimod.abam)$"Sum Sq"/sum(anova(consthimod.abam)$"Sum Sq"),digits=3))

#Figure 2 (deviance/variance explained by each factor in the models)
propdev.tsmeall=cbind(propdev.tsmesurv[2:12,1],c(0,propdev.tsmegerm[2,2],propdev.tsmegerm[1,2],0,0,0,0,propdev.tsmegerm[3,2],0,0,0),propdev.tsmesurv[2:12,2],propdev.tsmehi[,2])
colnames(propdev.tsmeall)=c("Predictor","Germination","Survival","Growth")
propdev.tsmeall[,2:4]=as.numeric(propdev.tsmeall[,2:4])
propdev.abamall=cbind(propdev.abamsurv[2:12,1],c(propdev.abamgerm[2,2],propdev.abamgerm[3,2],propdev.abamgerm[1,2],0,propdev.abamgerm[4,2],0,0,0,0,0,0),propdev.abamsurv[2:12,2],propdev.abamhi[,2])
colnames(propdev.abamall)=c("Predictor","Germination","Survival","Growth")
propdev.abamall[,2:4]=as.numeric(propdev.abamall[,2:4])
propdev.tsheall=cbind(propdev.tshesurv[2:12,1],c(propdev.tshegerm[1:2,2],0,0,propdev.tshegerm[3,2],0,0,0,0,0,0),propdev.tshesurv[2:12,2],propdev.tshehi[,2])
colnames(propdev.tsheall)=c("Predictor","Germination","Survival","Growth")
propdev.tsheall[,2:4]=as.numeric(propdev.tsheall[,2:4])
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
x=barplot(cbind(propdev.tsmeall2[1:11,2:4],c(1/6,1/6,0,0,1/6,1/6,0,0,1/6,1/6,0)), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="", space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15), ylim=c(0,1), las=1)
#labeloc=c(c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+(as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+(as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])+(+as.numeric(propdev.tsmeall[8,2])+as.numeric(propdev.tsmeall[9,2])/2)))
#text(x[1],c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])/2),labels=propdev.tsmeall[1:9,1],pos=1)
#text(x[1],c(as.numeric(propdev.tsmeall[1,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])/2,as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+(as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+(as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])/2),as.numeric(propdev.tsmeall[1,2])+as.numeric(propdev.tsmeall[2,2])+as.numeric(propdev.tsmeall[3,2])+as.numeric(propdev.tsmeall[4,2])+as.numeric(propdev.tsmeall[5,2])+as.numeric(propdev.tsmeall[6,2])+as.numeric(propdev.tsmeall[7,2])+(+as.numeric(propdev.tsmeall[8,2])+as.numeric(propdev.tsmeall[9,2])/2)),labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),pos=1)
labelsloc=c(1/12,3/12,5/12,7/12,9/12,11/12)
text(2.5,labelsloc,labels=c("Elevation","Competition","Origin","Elevation*Competition","Elevation*Origin","Origin*Competition"),cex=.9, adj=0)
mtext("Tsuga mertensiana",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("a",side=3,line=.5, adj=-.1, cex=.9)
#ABAM
x=barplot(cbind(propdev.abamall2[1:11,2:4],c(rep(NA,times=11))), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="", space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15),  ylim=c(0,1), las=1)
#labelsloc=c(.2,.351,.48,.63,.83,.97)
#text(2.3,labelsloc,labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),cex=.9, adj=1)
mtext("Abies amabilis",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("b",side=3,line=.5, adj=-.1, cex=.9)
mtext("Proportion variance explained", side=2,line=3,cex=.9)
#TSHE
x=barplot(cbind(propdev.tsheall2[1:11,2:4],c(rep(NA,times=11))), ylab="", col=colors2,cex.names=1.1,cex.lab=1.1,cex.axis=1.1,xlab="Vital rate",space=c(.5,.5,.5,.5),width=c(.5,.5,.5,.15), ylim=c(0,1), las=1)
#labelsloc=c(.2,.351,.48,.63,.83,.97)
#text(2.3,labelsloc,labels=c("Elevation","Origin","Competition","Elevation*Origin", "Elevation*Competition","Origin*Competition"),cex=.9, adj=1)
mtext("Tsuga heterophylla",side=3,line=.5, adj=0, font=3, cex=.9)
mtext("c",side=3,line=.5, adj=-.1, cex=.9)

#Fig. 1 hypothesis graph
quartz(height=7,width=5)
par(mfrow=c(3,1),mai=c(.6,.7,.2,.1), omi=c(.7,.1,.2,.2))
x<-c(0,50,100,150,200)
#EFfect of origin
#par(new=T)
perf<-c(.4,1,1.5,1,.4)
plot(perf~x,ylab="",type="p",bty="l",xlab="",pch=21,cex=1.5,las=1, cex.axis=1.2,bg="dark gray",col.axis="white", ylim=c(0,2),xlim=c(-35,200))
polygon(c(50,50,150,150),c(0,4, 4,0),col="light gray", border="light gray")
lines(perf~x,lty=1)
points(perf~x,pch=21,cex=1.5,bg="black")
axis(2, at = c(0,2), labels = c("0","+"), tick = FALSE, cex.axis=1,las=1)
mtext("Performance", side=2, adj=.8, line=2)
mtext("a", side=3, adj=-.15, line=.8)
mtext("Observed range",adj=0.6,cex=.9, line=-.2)
origin<-as.matrix(cbind(c(.6,.4,.2),c(1.5,1,.5),c(1.08,1.5,.98),c(.5,1,1.5),c(.2,.4,.6)))
plot(origin[1,]~x,ylab="",type="p",bty="l",xlab="",pch=21,cex=1.5,las=1, cex.axis=1.2,bg="darkred",col.axis="white", ylim=c(0,2),xlim=c(-35,200))
polygon(c(50,50,150,150),c(0,4, 4,0),col="light gray", border="light gray")
lines(origin[1,]~x,lty=1)
points(origin[1,]~x,pch=21,cex=1.5,bg="black")
lines(origin[2,]~x,lty=2)
points(origin[2,]~x,pch=21,cex=1.5,bg="dark gray")
lines(origin[3,]~x,lty=3)
points(origin[3,]~x,pch=21,cex=1.5,bg="white")
axis(2, at = c(0,2), labels = c("0","+"), tick = FALSE, cex.axis=1,las=1)
mtext("Performance", side=2, adj=.8, line=2)
mtext("b", side=3, adj=-.15, line=.8)
legend(-42,2.05,legend=c("Lower limit","Mid-range", "Upper limit"),bty="n",pch=21,pt.bg=c("black","dark gray","white"),lty=c(1,2,3),angle=45,cex=.8)
mtext("Seed origin",adj=0.04,cex=.9, line=-.6)
#effect of comptition
effect<-c(-2,-1.5,-1,-.5,-.1)
q=barplot(effect, col="dark gray", col.axis="white",ylim=c(-4,4), xlim=c(0,200), yaxt="n",width=30, space=.45)
polygon(c(72,72,158,158),c(4,-4, -4,4),col="light gray", border="light gray")
par(new=T)
barplot(effect, col="dark gray", col.axis="white",ylim=c(-4,4), xlim=c(0,200), width=30, space=.45)
mtext("Effect of neighbors", adj=1,side=2, line=2)
axis(2, at = c(-4,0,4), labels = c( "-","0","+"), tick = FALSE, cex.axis=1,las=1)
abline(h=0, lwd=1, lty=2)
mtext("c", side=3, adj=-.15, line=1.5)
mtext("Limit determined by climate", adj=1.1,side=3, line=.8,cex=.9)
arrows(158,2.5,158,4,length=.1,angle=45,code=1,lwd=2)
mtext("Limit determined by competition", adj=-.15,side=3, line=.8, cex=.9)
arrows(72,2.5,72,4,length=.1,angle=45,code=1, lwd=2)
axis(1, at = q, tick = T,labels = c( "Below","Lower","Mid", "Upper", "Above"), cex.axis=.9)
mtext("Range position", line=-14)

####Figure 3
###Survival & Height increment across all treatments, by source elevation, using estimates from models.
###Survival
tsmedat$PlantedStand=factor(tsmedat$PlantedStand)
tshedat$PlantedStand=factor(tshedat$PlantedStand)
abamdat$PlantedStand=factor(abamdat$PlantedStand)
tsmedat$OriginStand=factor(tsmedat$OriginStand)
tshedat$OriginStand=factor(tshedat$OriginStand)
abamdat$OriginStand=factor(abamdat$OriginStand)
tsmedat$Block=factor(tsmedat$Block)
tshedat$Block=factor(tshedat$Block)
abamdat$Block=factor(abamdat$Block)
abamdat$OriginRangeLoc<-NA
tsmedat$OriginRangeLoc<-NA
tshedat$OriginRangeLoc<-NA
abamdat[abamdat$OriginStand=="704",]$OriginRangeLoc<-"Lower"
abamdat[abamdat$OriginStand=="1603",]$OriginRangeLoc<-"Upper"
abamdat[abamdat$OriginStand=="1064"|abamdat$OriginStand=="1197",]$OriginRangeLoc<-"Mid"
tshedat[tshedat$OriginStand=="704",]$OriginRangeLoc<-"Mid"
tshedat[tshedat$OriginStand=="1197",]$OriginRangeLoc<-"Upper"
tsmedat[tsmedat$OriginStand=="1197",]$OriginRangeLoc<-"Lower"
tsmedat[tsmedat$OriginStand=="1460",]$OriginRangeLoc<-"Mid"
tsmedat[tsmedat$OriginStand=="1603",]$OriginRangeLoc<-"Upper"
#TSME
tsmesurv<-tapply(as.numeric(tsmedat$StatusDate5),list(tsmedat$Block,tsmedat$PlantedStand,tsmedat$OriginRangeLoc),sum)/tapply(as.numeric(tsmedat$StatusDate5),list(tsmedat$Block,tsmedat$PlantedStand,tsmedat$OriginRangeLoc),length)
tsmesurv.org=c(tsmesurv)
tsmeblock=c(rep(rownames(tsmesurv), times=15))
tsmeorg=c(rep("Low", times=125),rep("Mid", times=125),rep("Up", times=125)) 
tsmesurvorg=cbind(tsmeblock,tsmeorg,tsmesurv.org)
colnames(tsmesurvorg)=c("Block","OrgStand","Surv")
tsmesurvorg<-tsmesurvorg[-which(is.na(tsmesurvorg[,3])),]
PS=substr(tsmesurvorg[,1], 1, 4)
tsmemnsorg=tapply(as.numeric(tsmesurvorg[,3]),list(PS,tsmesurvorg[,2]),mean)
tsmesesorg=tapply(as.numeric(tsmesurvorg[,3]),list(PS,tsmesurvorg[,2]),sd)/sqrt(5)
###ABAM
abamsurv<-tapply(as.numeric(abamdat$StatusDate5),list(abamdat$Block,abamdat$PlantedStand,abamdat$OriginRangeLoc),sum)/tapply(as.numeric(abamdat$StatusDate5),list(abamdat$Block,abamdat$PlantedStand,abamdat$OriginRangeLoc),length)
abamsurv.org=c(abamsurv)
abamblock=c(rep(rownames(abamsurv), times=15))
abamorg=c(rep("Low", times=125),rep("Mid", times=125),rep("Up", times=125)) 
abamsurvorg=cbind(abamblock,abamorg,abamsurv.org)

colnames(abamsurvorg)=c("Block","OrgStand","Surv")
abamsurvorg<-abamsurvorg[-which(is.na(abamsurvorg[,3])),]
PS=substr(abamsurvorg[,1], 1, 4)
abammnsorg=tapply(as.numeric(abamsurvorg[,3]),list(PS,abamsurvorg[,2]),mean)
abamsesorg=tapply(as.numeric(abamsurvorg[,3]),list(PS,abamsurvorg[,2]),sd)/sqrt(5)
###TSHE
tshesurv<-tapply(as.numeric(tshedat$StatusDate5),list(tshedat$Block,tshedat$PlantedStand,tshedat$OriginRangeLoc),sum)/tapply(as.numeric(tshedat$StatusDate5),list(tshedat$Block,tshedat$PlantedStand,tshedat$OriginRangeLoc),length)
tshesurv.org=c(tshesurv)
tsheblock=c(rep(rownames(tshesurv), times=6))
tsheorg=c(rep("Mid", times=45),rep("Up", times=45)) 
tshesurvorg=cbind(tsheblock,tsheorg,tshesurv.org)
colnames(tshesurvorg)=c("Block","OrgStand","Surv")
tshesurvorg<-tshesurvorg[-which(is.na(tshesurvorg[,3])),]
PS=substr(tshesurvorg[,1], 1, 4)
tshemnsorg=tapply(as.numeric(tshesurvorg[,3]),list(PS,tshesurvorg[,2]),mean)
tshesesorg=tapply(as.numeric(tshesurvorg[,3]),list(PS,tshesurvorg[,2]),sd)/sqrt(5)
####HeightIncrement
tsmehi<-tapply(as.numeric(tsmedat$annhi),list(tsmedat$PlantedStand,tsmedat$OriginRangeLoc),mean,na.rm=T)
tsmehi.org=c(tsmehi)
tsmesehiorg=tapply(as.numeric(tsmedat$annhi),list(tsmedat$PlantedStand,tsmedat$OriginRangeLoc),sd,na.rm=T)/(sqrt(tapply(as.numeric(tsmedat$annhi),list(tsmedat$PlantedStand,tsmedat$OriginRangeLoc),length)))
#ABAM
abamhi<-tapply(as.numeric(abamdat$annhi),list(abamdat$PlantedStand,abamdat$OriginRangeLoc),mean,na.rm=T)
abamhi.org=c(abamhi)
abamsehiorg=tapply(as.numeric(abamdat$annhi),list(abamdat$PlantedStand,abamdat$OriginRangeLoc),sd,na.rm=T)/(sqrt(tapply(as.numeric(abamdat$annhi),list(abamdat$PlantedStand,abamdat$OriginRangeLoc),length)))
#TSHE
tshehi<-tapply(as.numeric(tshedat$annhi),list(tshedat$PlantedStand,tshedat$OriginRangeLoc),mean,na.rm=T)
tshehi.org=c(tshehi)
tshesehiorg=tapply(as.numeric(tshedat$annhi),list(tshedat$PlantedStand,tshedat$OriginRangeLoc),sd,na.rm=T)/(sqrt(tapply(as.numeric(tshedat$annhi),list(tshedat$PlantedStand,tshedat$OriginRangeLoc),length)))

#####Figure 3: plot of germination, survival and growth, by species
rownames(tsmemnsorg)=c(1460,1197,1064,1676,1603)
tsmemnsorg2=tsmemnsorg[order(rownames(tsmemnsorg)),] 
rownames(tsmesesorg)=c(1460,1197,1064,1676,1603)
tsmesesorg2=tsmesesorg[order(rownames(tsmesesorg)),] 
rownames(tshemnsorg)=c(3,2,1)
tshemnsorg3=tshemnsorg[order(rownames(tshemnsorg)),] 
tshemnsorg2=cbind(tshemnsorg3[,2],tshemnsorg3[,1])
rownames(tshesesorg)=c(3,2,1)
tshesesorg3=tshesesorg[order(rownames(tshesesorg)),] 
tshesesorg2=cbind(tshesesorg3[,2],tshesesorg3[,1])
rownames(abammnsorg)=c(3,5,1,4,2)
abammnsorg2=abammnsorg[order(rownames(abammnsorg)),] 
rownames(abamsesorg)=c(3,5,1,4,2)
abamsesorg3=abamsesorg[order(rownames(abamsesorg)),] 

quartz(height=7,width=12)
par(mfcol=c(3,3),mai=c(.6,.7,.2,.1), omi=c(.7,.01,.2,.1))
#Germination, by species
#TSME
tsmegerm<-tapply(as.numeric(tsmegermdat$TotalGerms),list(tsmegermdat$Block,tsmegermdat$Stand,tsmegermdat$Origin),sum)/tapply(as.numeric(tsmegermdat$SeedsAdded),list(tsmegermdat$Block,tsmegermdat$Stand,tsmegermdat$Origin),sum)
tsmegerm.org=c(tsmegerm)
length(tsmegerm.org)
tsmeblock=c(rep(rownames(tsmegerm), times=15))
tsmeorg=c(rep("1197", times=125),rep("1460", times=125),rep("1603", times=125)) 
tsmegermorg=cbind(tsmeblock,tsmeorg,tsmegerm.org)
colnames(tsmegermorg)=c("Block","OrgStand","PropGerm")
tsmegermorg<-tsmegermorg[-which(is.na(tsmegermorg[,3])),]
PS=substr(tsmegermorg[,1], 1, 4)
tsmemngermorg=tapply(as.numeric(tsmegermorg[,3]),list(PS,tsmegermorg[,2]),mean)
tsmegermseorg=tapply(as.numeric(tsmegermorg[,3]),list(PS,tsmegermorg[,2]),sd)/sqrt(5)
###ABAM
abamgerm<-tapply(as.numeric(abamgermdat$TotalGerms),list(abamgermdat$Block,abamgermdat$Stand,abamgermdat$Origin),sum)/tapply(as.numeric(abamgermdat$SeedsAdded),list(abamgermdat$Block,abamgermdat$Stand,abamgermdat$Origin),sum)
abamgerm.org=c(abamgerm)
abamblock=c(rep(rownames(abamgerm), times=10))
abamorg=c(rep("1197", times=125),rep("1603", times=125)) 
abamgermorg=cbind(abamblock,abamorg,abamgerm.org)
colnames(abamgermorg)=c("Block","OrgStand","PropGerm")
abamgermorg<-abamgermorg[-which(is.na(abamgermorg[,3])),]
PS=substr(abamgermorg[,1], 1, 4)
abammngermorg=tapply(as.numeric(abamgermorg[,3]),list(PS,abamgermorg[,2]),mean)
abamgermseorg=tapply(as.numeric(abamgermorg[,3]),list(PS,abamgermorg[,2]),sd)/sqrt(5)
###TSHE
tshegerm<-tapply(as.numeric(tshegermdat$TotalGerms),list(tshegermdat$Block,tshegermdat$Stand,tshegermdat$Origin),sum)/tapply(as.numeric(tshegermdat$SeedsAdded),list(tshegermdat$Block,tshegermdat$Stand,tshegermdat$Origin),sum)
tshegerm.org=c(tshegerm)
tsheblock=c(rep(rownames(tshegerm), times=6))
tsheorg=c(rep("704", times=45),rep("1197", times=45)) 
tshegermorg=cbind(tsheblock,tsheorg,tshegerm.org)
colnames(tshegermorg)=c("Block","OrgStand","Germ")
tshegermorg<-tshegermorg[-which(is.na(tshegermorg[,3])),]
PS=substr(tshegermorg[,1], 1, 4)
tshemngermorg=tapply(as.numeric(tshegermorg[,3]),list(PS,tshegermorg[,2]),mean)
tshegermseorg=tapply(as.numeric(tshegermorg[,3]),list(PS,tshegermorg[,2]),sd)/sqrt(5)
X<-c(1,2,3,4,5)
quartz(height=7,width=10)
par(mfcol=c(3,3), mai=c(.5,.7,.3,.1), omi=c(.8,.01,.3,.2))
#Germination
#TSME
plot(tsmemngermorg[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.05),type="p",bty="l",pch=21,cex=1.8,las=1,bg="darkred", cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
#polygon(c(2,2,4,4),c(0,.5,.5,0),col="light gray", border="light gray")
polygon(c(2,2,4,4),c(0,.05,.05,0),col="light gray", border="light gray")#for inconsistent y axes
#axis(2,at=c(0,.1,.2,.3,.4),las=1, cex.axis=1.3)
axis(2,at=c(0,.01,.02,.03,.04,.05),las=1, cex.axis=1.3)
mtext("Tsuga mertensiana",side=3,line=2, adj=0, font=3)
mtext("a",side=3,line=2, adj=-.1)
alltsmese<-c(tsmegermseorg)
alltsme<-c(tsmemngermorg)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(alltsme)){
  arrows(x[i],alltsme[i]-alltsmese[i],x[i],alltsme[i]+alltsmese[i],length=.05,angle=90,code=0)}
lines(tsmemngermorg[,1]~X,lty=3)
points(tsmemngermorg[,1]~X,pch=21,cex=1.8,bg="black")
lines(tsmemngermorg[,2]~X,lty=2)
points(tsmemngermorg[,2]~X,pch=21,cex=1.8,bg="dark gray")
lines(tsmemngermorg[,3]~X,lty=1)
points(tsmemngermorg[,3]~X,pch=21,cex=1.8,bg="white")
legend(1,.055,legend=c("Lower limit","Mid-range", "Upper limit"),bty="n",pch=21,pt.bg=c("black","dark gray","white"),angle=45,cex=1.1,lty=c(3,2,1), pt.cex=1.5)
axis(1, at = c(1,2,3,4,5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##ABAM
plot(abammngermorg[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.4),type="p",bty="l",pch=21,cex=1.8,las=1,bg="dark gray", cex.main=1.3)
polygon(c(2,2,4,4),c(0,.5,.5,0),col="light gray", border="light gray")
axis(2,at=c(0,.1,.2,.3,.4),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("b",side=3,line=1, adj=-.1)
allabamse<-c(abamgermseorg)
allabam<-c(abammngermorg)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(allabam)){
  arrows(x[i],allabam[i]-allabamse[i],x[i],allabam[i]+allabamse[i],length=.05,angle=90,code=0)}
lines(abammngermorg[,1]~X,lty=2)
points(abammngermorg[,1]~X,pch=21,cex=1.8,bg="dark gray")
lines(abammngermorg[,2]~X,lty=1)
points(abammngermorg[,2]~X,pch=21,cex=1.8,bg="white")

axis(1, at = c(1,2,3,4,5), labels = c( "(668)","(704)","(1064)","(1603)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Proportion Seeds Germinating", side=2, line=3)

##TSHE
tsheplot=plot(tshemngermorg[,1]~X[3:5],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.05),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)#ylim=c(0,.05) for inconsistent y axes
polygon(c(2,2,4,4),c(0,.06,.06,0),col="light gray", border="light gray")#for inconsistent y axes
axis(2,at=c(0,.01,.02,.03,.04,.05),las=1, cex.axis=1.3)#for inconsistent y axes
#polygon(c(2,2,4,4),c(0,.5,.5,0),col="light gray", border="light gray")
#axis(2,at=c(0,.1,.2,.3,.4),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("c",side=3,line=1, adj=-.1)
alltshese<-c(tshegermseorg)
alltshe<-c(tshemngermorg)
x<-c(3,4,5,3,4,5)
for (i in 1:length(alltshe)){
  arrows(x[i],alltshe[i]-alltshese[i],x[i],alltshe[i]+alltshese[i],length=.05,angle=90,code=0)}

lines(tshemngermorg[,1]~X[3:5],lty=2)
points(tshemngermorg[,1]~X[3:5],pch=21,cex=1.8,bg="dark gray")
lines(tshemngermorg[,2]~X[3:5],lty=1)
points(tshemngermorg[,2]~X[3:5],pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1,2,3,4,5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext ("Range position (m)", line=-13, cex=.9)

#Survival
#TSME
plot(tsmemnsorg2[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.5),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(0,.6,.6,0),col="light gray", border="light gray")
axis(2,at=c(0,.1,.2,.3,.4,.5),las=1, cex.axis=1.3)
mtext("Tsuga mertensiana",side=3,line=2, adj=0, font=3)
mtext("d",side=3,line=2, adj=-.1)
alltsmese<-c(tsmesesorg2)
alltsme<-c(tsmemnsorg2)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(alltsme)){
  arrows(x[i],alltsme[i]-alltsmese[i],x[i],alltsme[i]+alltsmese[i],length=.05,angle=90,code=0)}
lines(tsmemnsorg2[,1]~X,lty=3)
points(tsmemnsorg2[,1]~X,pch=21,cex=1.8,bg="black")
lines(tsmemnsorg2[,2]~X,lty=2)
points(tsmemnsorg2[,2]~X,pch=21,cex=1.8,bg="dark gray")
lines(tsmemnsorg2[,3]~X,lty=1)
points(tsmemnsorg2[,3]~X,pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##ABAM
plot(abammnsorg2[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.5),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(0,.6,.6,0),col="light gray", border="light gray")
axis(2,at=c(0,.1,.2,.3,.4,.5),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("e",side=3,line=1, adj=-.1)
allabamse<-c(abamsesorg2)
allabam<-c(abammnsorg2)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(allabam)){
  arrows(x[i],allabam[i]-allabamse[i],x[i],allabam[i]+allabamse[i],length=.05,angle=90,code=0)}
lines(abammnsorg2[,1]~X,lty=3)
points(abammnsorg2[,1]~X,pch=21,cex=1.8,bg="black")
lines(abammnsorg2[,2]~X,lty=2)
points(abammnsorg2[,2]~X,pch=21,cex=1.8,bg="dark gray")
lines(abammnsorg2[,3]~X,lty=1)
points(abammnsorg2[,3]~X,pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "(668)","(704)","(1064)","(1603)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Proportion Surviving Seedlings", side=2, line=3)

##TSHE
tsheplot=plot(tshemnsorg2[,1]~X[3:5],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,.5),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(0,.6,.6,0),col="light gray", border="light gray")
axis(2,at=c(0,.1,.2,.3,.4,.5),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("f",side=3,line=1, adj=-.1)
alltshese<-c(tshesesorg2)
alltshe<-c(tshemnsorg2)
x<-c(3,4,5,3,4,5)
for (i in 1:length(alltshe)){
  arrows(x[i],alltshe[i]-alltshese[i],x[i],alltshe[i]+alltshese[i],length=.05,angle=90,code=0)}

lines(tshemnsorg2[,1]~X[3:5],lty=2)
points(tshemnsorg2[,1]~X[3:5],pch=21,cex=1.8,bg="dark gray")
lines(tshemnsorg2[,2]~X[3:5],lty=1)
points(tshemnsorg2[,2]~X[3:5],pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1,2,3,4,5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext("Range position (m)", line=-13, cex=.9)

####Height Increment
plot(tsmehi[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,1),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(0,1,1,0),col="light gray", border="light gray")
axis(2,at=c(0,.2,.4,.6,.8,1),las=1, cex.axis=1.3)
mtext("Tsuga mertensiana",side=3,line=2, adj=0, font=3)
mtext("g",side=3,line=2, adj=-.1)
alltsmese<-c(tsmesehiorg)
alltsme<-c(tsmehi)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(alltsme)){
  arrows(x[i],alltsme[i]-alltsmese[i],x[i],alltsme[i]+alltsmese[i],length=.05,angle=90,code=0)}
lines(tsmehi[,1]~X,lty=3)
points(tsmehi[,1]~X,pch=21,cex=1.8,bg="black")
lines(tsmehi[,2]~X,lty=2)
points(tsmehi[,2]~X,pch=21,cex=1.8,bg="dark gray")
lines(tsmehi[,3]~X,lty=1)
points(tsmehi[,3]~X,pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "(1064)", "(1197)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)

##ABAM
plot(abamhi[,1]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(-1,1),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(-1,1,1,-1),col="light gray", border="light gray")
axis(2,las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("h",side=3,line=1, adj=-.1)
allabamse<-c(abamsehiorg)
allabam<-c(abamhi)
x<-c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
for (i in 1:length(allabam)){
  arrows(x[i],allabam[i]-allabamse[i],x[i],allabam[i]+allabamse[i],length=.05,angle=90,code=0)}
lines(abamhi[,1]~X,lty=3)
points(abamhi[,1]~X,pch=21,cex=1.8,bg="black")
lines(abamhi[,2]~X,lty=2)
points(abamhi[,2]~X,pch=21,cex=1.8,bg="dark gray")
lines(abamhi[,3]~X,lty=1)
points(abamhi[,3]~X,pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "(668)","(704)","(1064)","(1603)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Annual height increment (cm/year)", side=2, line=3)
##TSHE
plot(tshehi[,1]~X[3:5],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,1),type="p",bty="l",pch=21,cex=1.8,las=1,bg="black", cex.main=1.3)
polygon(c(2,2,4,4),c(0,1,1,0),col="light gray", border="light gray")
axis(2,at=c(0,.2,.4,.6,.8,1),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("i",side=3,line=1, adj=-.1)
alltshese<-c(tshesehiorg)
alltshe<-c(tshehi)
x<-c(3,4,5,3,4,5)
for (i in 1:length(alltshe)){
  arrows(x[i],alltshe[i]-alltshese[i],x[i],alltshe[i]+alltshese[i],length=.05,angle=90,code=0)}
lines(tshehi[,1]~X[3:5],lty=2)
points(tshehi[,1]~X[3:5],pch=21,cex=1.8,bg="dark gray")
lines(tshehi[,2]~X[3:5],lty=1)
points(tshehi[,2]~X[3:5],pch=21,cex=1.8,bg="white")
axis(1, at = c(1,2,3,4,5), labels = c( "", "","(668)", "(704)","(1064)"), tick = FALSE, cex.axis=1.1, line=-.5)
axis(1, at = c(1,2,3,4,5), labels = c( "Below","Lower","Mid", "Upper", "Above"), tick = FALSE, cex.axis=1.1, line=0.5)
mtext ("Range position (m)", line=-13, cex=.9)

#plot, with change in effects of comp vs no comp
compsurvmod.tsme=survreg(Surv(time1,time2, type="interval2")~-1+as.factor(PlantedStand)*as.factor(Canopy)*as.factor(Understory), dist="lognormal", data=tsmedat)
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
#get coefs and their ses
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

###Barplot
quartz(height=7,width=7)
par(mfcol=c(3,2),mai=c(.6,.8,.2,.1), omi=c(.7,.01,.2,.2))
plottsme<-barplot(as.matrix(rbind(tsmesurvund,tsmesurvcan)),width=.9,ylab="",xlab="",names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("a",side=3,line=1, adj=-.1)
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

##ABAM
plotabam<-barplot(as.matrix(rbind(abamsurvund,abamsurvcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.5,1.5),las=1,col=c("palegreen1","darkgreen"),xaxt='n', yaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.5,1.5,1.5,-1.5),col="light gray", border="light gray")
par(new=T)
plotabam<-barplot(as.matrix(rbind(abamsurvund,abamsurvcan)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1.5,1.5),las=1,col=c("palegreen1","darkgreen"),,beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("b",side=3,line=1, adj=-.1)
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
mtext("c",side=3,line=1, adj=-.1)
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
mtext("Range position (m)",line=-13.5, adj=.6)

#Now ANNUAL HI
plottsme<-barplot(as.matrix(rbind(tsmehiund,tsmehican)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),ylim=c(-1,1),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
mtext("Tsuga mertensiana",side=3,line=1, adj=0, font=3)
mtext("d",side=3,line=1, adj=-.1)
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
668,704,1064,1197,1460,1605,1676
#add error bars
for (i in 1:length(allhiabam.lcl)){
  arrows(x[i],allhiabam.lcl[i],x[i],allhiabam.ucl[i],length=.03,angle=90,code=0)}
axis(1, at = c(1.8,4.5,7.2,9.9,12.5), labels = c( "(668)","(704)","(1460)","(1605)","(1676)"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Abies amabilis",side=3,line=1, adj=0, font=3)
mtext("e",side=3,line=1, adj=-.1)
mtext("Effect of neighbors on growth", side=2, line=4.3)
mtext("(Difference in annual increment, relative to no neighbors)", side=2, line=3, cex=.9)
##TSHE
tshehiund2=c(NA,NA,tshehiund)
tshehican2=c(NA,NA,tshehican)
plottshe<-barplot(as.matrix(rbind(tshehiund2,tshehican2)),ylab="",xlab="",width=.9,names.arg=c("","","","",""),xaxt='n', yaxt='n',ylim=c(-1.7,.5),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
polygon(c(3.5,3.5,11,11),c(-1.7,.5,.5,-1.7),col="light gray", border="light gray")
mtext("Tsuga heterophylla",side=3,line=1, adj=0, font=3)
mtext("f",side=3,line=1, adj=-.1)
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
mtext("Range position (m)",line=-13.5, adj=.6)
####Microclimate
data<-read.csv("data/AllStands_clim1.csv", header=TRUE)
data$Elevation_mfact<-as.factor(data$Elevation_m)
data$Elevation_m<-as.numeric(data$Elevation_m)
data$Block<-as.factor(data$Block)
constmod.snow<-lmer(snow_cover_duration~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.snow, type="III")
constmod.gdd<-lmer(GDD_total~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.gdd, type="III")#
summary(constmod.gdd)
constmod.light<-lmer(Light_Mean~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.light, type="III")#
#########Soil analyses
soilsdat<-read.csv("data/SoilsData.csv", header=TRUE)
soilsdat$Elev_m<-as.factor(soilsdat$Elev_m)
soilmoist.lm<-lm(soilsdat$SoilMoist_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilmoist.lm)
soilC.lm<-lm(soilsdat$C_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilC.lm)
soilH.lm<-lm(soilsdat$H_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilH.lm)
soilN.lm<-lm(soilsdat$N_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilN.lm)
soilC.aov<-aov(soilsdat$C_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilC.aov)
soilH.aov<-aov(soilsdat$H_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilH.aov)
soilN.aov<-aov(soilsdat$N_pce~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(soilN.aov)
CanCov.lm<-lm(soilsdat$CanopyCover~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(CanCov.lm)
###summary of GDD by elevation & gap/nongap (means)
data<-read.csv("data/AllStands_clim1.csv", header=TRUE)
(meangdd<-tapply(data$GDD_total,list(data$Canopy,data$Elevation_m),mean))
data<-data[-48,]
data<-data[-88,]
(meanlight<-tapply(data$Light_Mean,list(data$Canopy,data$Elevation_m),mean))

#######summary of slope and aspect by elevation & gap/nongap (means)
spatdat<-read.csv("data/GapSizeData.csv", header=TRUE)
(meanslope<-tapply(spatdat$Slope...,list(spatdat$Canopy,spatdat$Stand),mean))
slope.lm<-lm(spatdat$Slope...~-1+spatdat$Stand*spatdat$Canopy)
summary(slope.lm)
slope.aov<-aov(spatdat$Slope...~-1+spatdat$Stand*spatdat$Canopy)
summary(slope.aov)
(meanaspect<-tapply(spatdat$Aspect..degrees,list(spatdat$Canopy,spatdat$Stand),mean, na.rm=T))
aspect.lm<-lm(spatdat$Aspect..degrees~-1+spatdat$Stand*spatdat$Canopy)
summary(aspect.lm)
aspect.aov<-aov(spatdat$Aspect..degrees~-1+spatdat$Stand*spatdat$Canopy)
summary(aspect.aov)
#PAR Data
pardat<-read.csv("PARdata.csv", header=TRUE)
justpar<-cbind((pardat[,9]),(pardat[,10]),(pardat[,11]))
parmn<-rowMeans(justpar, na.rm=TRUE)
(meanpar<-tapply(parmn,list(pardat$Canopy,pardat$Stand),mean, na.rm=TRUE))
parmn<-rowMeans(justpar)
par.lm<-lm(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.lm)
par.aov<-aov(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.aov)