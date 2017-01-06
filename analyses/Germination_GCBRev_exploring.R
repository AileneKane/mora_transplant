#New germination models
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library(dplyr)
library(plyr)
library(glmmADMB)
library(boot)
library(car)
library(bbmle)
setwd("~/git/mora_transplant")
dat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
head(dat);tail(dat)
dim(dat)#2939  13
dat$TotalGerms<-as.numeric(dat$TotalGerms)
dat$SeedsAdded<-as.numeric(dat$SeedsAdded)

#dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms<-50#replace any success that were 
##First, select only rows in which seeds were added
addat<-dat[dat$SeedsAdded>0,]
head(addat);tail(addat)
head(dat);tail(dat)
#dim(addat)#1080 plots
dat2011<-addat[addat$Year=="2011",]#
dat2012<-addat[addat$Year=="2012",]#
alldat2011<-dat[dat$Year=="2011",]#
alldat2012<-dat[dat$Year=="2012",]#
dim(dat2011)#620 plots
dim(dat2012)#460
dim(alldat2011)#1500 plots
dim(alldat2012)#1439

#select data by species
#for TSME, model selection show that moss does not affect germination, so interpret models without moss (to keep more data)
tsmegermdat<-dat[dat$SpPlant=="TSME",]
tsmegermdat$TotalGerms<-as.numeric(tsmegermdat$TotalGerms)
tsmegermdat$TotalFails<-as.numeric(tsmegermdat$TotalFails)
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)
tsmegermdat$Stand=factor(tsmegermdat$Stand)
tsmegermdat$Origin=factor(tsmegermdat$Origin)
tsmegermdat$Block=factor(tsmegermdat$Block)
tsmegermdat$Year=factor(tsmegermdat$Year)
tsmegermdat$Canopy=factor(tsmegermdat$Canopy)
tsmegermdat$Understory=factor(tsmegermdat$Understory)
tsmegermdat2<-subset(tsmegermdat,select=c(TotalGerms,TotalFails,Block,Year,Stand,Origin,Canopy,Understory))
tsmegermdat3<-na.omit(tsmegermdat2)
dim(tsmegermdat3)
tsmey3=cbind(tsmegermdat3$TotalGerms,tsmegermdat3$TotalFails)
#Try using ALL data with random effect of year
#consgermmod.tsme<-glmmadmb(tsmey3 ~ Stand+Origin+Canopy+Understory +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block)+(1|Year),family="binomial", data=tsmegermdat3)#now merge in microclimate variables
#above won't converge
#Fit ZIP model
consgermmod.tsme.zip<-glmmadmb(TotalGerms ~ Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="poisson", data=tsmegermdat, zeroInflation =TRUE)#
#consgermmod.tsme2<-glmmadmb(TotalGerms ~ Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="poisson", data=tsmegermdat, zeroInflation =FALSE)#doesn't fit
#consgermmod.tsme.nb<-glmmadmb(TotalGerms ~ Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="nbinom", data=tsmegermdat, zeroInflation =TRUE)#doesn't fit
#consgermmod.tsme.nb1<-glmmadmb(TotalGerms ~ Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="nbinom1", data=tsmegermdat, zeroInflation =TRUE)#doesn't fit
summary(consgermmod.tsme.zip)
anova(consgermmod.tsme.zip)
#Now proportion of deviance
nullgermod.tsme.zip<-glmmadmb(TotalGerms~1+(1|Block), data=tsmegermdat, family="poisson",zeroInflation =TRUE)
stand.tsme<-glmmadmb(tsmey ~ Stand+(1|Block),family="binomial", data=tsmegermdat)
standorg.tsme<-glmmadmb(tsmey ~ Stand  +Origin +(1|Block),family="binomial", data=tsmegermdat)
standorgcan<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +(1|Block),family="binomial", data=tsmegermdat)
standorgcanund<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint1<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint2<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint3<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+(1|Block),family="binomial", data=tsmegermdat, zeroInflation =TRUE)
standorgcanundint4<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint5<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint6<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint7<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
Stand<-anova(nullgermod.tsme,stand.tsme)$Deviance[2]
Org<-anova(stand.tsme,standorg.tsme)$Deviance[2]
Can<-anova(standorg.tsme,standorgcan)$Deviance[2]
Und<-anova(standorgcan,standorgcanund)$Deviance[2]
storgint<-anova(standorgcanund,standorgcanundint1)$Deviance[2]
stcanint<-anova(standorgcanundint1,standorgcanundint2)$Deviance[2]
stundint<-anova(standorgcanundint2,standorgcanundint3)$Deviance[2]
canorg<-anova(standorgcanundint3,standorgcanundint4)$Deviance[2]
undorg<-anova(standorgcanundint4,standorgcanundint5)$Deviance[2]
canund<-anova(standorgcanundint5,standorgcanundint6)$Deviance[2]
stcanund<-anova(standorgcanundint6,standorgcanundint7)$Deviance[2]
#threewayint<-anova(standorgcanundint6,consgermmod.tsme)$Deviance[2]#should be the same as the above, and it is
anovatab.tsme<-cbind(row.names(Anova(consgermmod.tsme)),c(Stand,Org,Can,Und,storgint,stcanint,stundint,canorg,undorg,canund,stcanund))
colnames(anovatab.tsme)<-c("Source","Deviance")
PropDev<-as.numeric(anovatab.tsme[,2])/sum(as.numeric(anovatab.tsme[,2]))
anovatab.tsme<-cbind(anovatab.tsme,PropDev)
write.csv(anovatab.tsme,"analyses/tsmegerm.propdev.csv", row.names = FALSE)
#For other two species, we need to add moss in
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
colnames(grcov2)[3]<-"BlockNum"
grcov2$Block<-paste(grcov2$Stand,grcov2$BlockNum,sep="")
grcov2$Canopy<-NA
grcov2[grcov2$Gap.=="Y",]$Canopy<-"CompAbsent"
grcov2[grcov2$Gap.=="N",]$Canopy<-"CompPresent"
colnames(dat2011)[2]<-"Elev"
germdat2<-join(dat2011, grcov2, by=c("Block","Plot","Canopy"), match="all")
dat<-subset(germdat2,select=c("Stand","Block","Plot","Canopy","Understory","SpPlant","Origin","Year","TotalGerms","SeedsAdded","TotalFails","Litter","Moss","Wood"))
abamgermdat<-dat[dat$SpPlant=="ABAM",]
#abamgermdat<-na.omit(abamgermdat)
abamgermdat$Stand=factor(abamgermdat$Stand)
abamgermdat$Origin=factor(abamgermdat$Origin)
abamgermdat$TotalGerms=as.numeric(abamgermdat$TotalGerms)
abamgermdat$Block=factor(abamgermdat$Block)
abamgermdat$Moss=as.numeric(abamgermdat$Moss)
abamgermdat$Moss_cent<-scale(abamgermdat$Moss)
abamgermdat$Year=as.factor(abamgermdat$Year)
abamgermdat$SeedsAdded=as.factor(abamgermdat$SeedsAdded)

abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)
germod.abam.moss<-glmmadmb(abamy~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
germod.abam.moss.zip<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+SeedsAdded+(1|Block)+(1|Year), data=abamgermdat, family="poisson", zeroInflation = TRUE)
#zip fits better than binomial, based on likelihood ratio test
nullgermod.abam<-glmmadmb(abamy~Moss_cent+(1|Block), data=abamgermdat, family="binomial")
stand.abam<-glmmadmb(abamy ~ Stand+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorg.abam<-glmmadmb(abamy ~ Stand  +Origin+Moss_cent +(1|Block),family="binomial", data=abamgermdat)
standorgcan<-glmmadmb(abamy ~ Stand  +Origin + Canopy+Moss_cent +(1|Block),family="binomial", data=abamgermdat)
standorgcanund<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint1<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint2<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint3<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint4<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint5<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint6<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
standorgcanundint7<-glmmadmb(abamy ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block),family="binomial", data=abamgermdat)
#germod.abam.moss
Stand<-anova(nullgermod.abam,stand.abam)$Deviance[2]
Org<-anova(stand.abam,standorg.abam)$Deviance[2]
Can<-anova(standorg.abam,standorgcan)$Deviance[2]
Und<-anova(standorgcan,standorgcanund)$Deviance[2]
storgint<-anova(standorgcanund,standorgcanundint1)$Deviance[2]
stcanint<-anova(standorgcanundint1,standorgcanundint2)$Deviance[2]
stundint<-anova(standorgcanundint2,standorgcanundint3)$Deviance[2]
canorg<-anova(standorgcanundint3,standorgcanundint4)$Deviance[2]
undorg<-anova(standorgcanundint4,standorgcanundint5)$Deviance[2]
canund<-anova(standorgcanundint5,standorgcanundint6)$Deviance[2]
stcanund<-anova(standorgcanundint6,standorgcanundint7)$Deviance[2]
#threewayint<-anova(standorgcanundint6,germod.abam.moss)$Deviance[2]#should be the same as the above and it is
abamanova<-Anova(germod.abam.moss)
abamanova<-abamanova[-5,]
anovatab.abam<-cbind(row.names(abamanova),c(Org,Can,Stand,Und,storgint,stcanint,stundint,canorg,undorg,canund,stcanund))
colnames(anovatab.abam)<-c("Source","Deviance")
PropDev<-as.numeric(anovatab.abam[,2])/sum(as.numeric(anovatab.abam[,2]))
anovatab.abam<-cbind(anovatab.abam,PropDev)
write.csv(anovatab.abam,"analyses/abamgerm.propdev.csv",row.names = FALSE)
#TSHE
tshegermdat<-dat[dat$SpPlant=="TSHE",]
tshegermdat<-na.omit(tshegermdat)
tshegermdat$Stand=factor(tshegermdat$Stand)
tshegermdat$Origin=factor(tshegermdat$Origin)
tshegermdat$TotalGerms=factor(tshegermdat$TotalGerms)
tshegermdat$Block=factor(tshegermdat$Block)
tshegermdat$Moss=as.numeric(tshegermdat$Moss)
tshegermdat$Moss_cent<-scale(tshegermdat$Moss)
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

germod.tshe.moss<-glmmadmb(tshey~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
nullgermod.tshe<-glmmadmb(tshey~Moss_cent+(1|Block), data=tshegermdat, family="binomial")
stand.tshe<-glmmadmb(tshey ~ Stand+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorg.tshe<-glmmadmb(tshey ~ Stand  +Origin+Moss_cent +(1|Block),family="binomial", data=tshegermdat)
standorgcan<-glmmadmb(tshey ~ Stand  +Origin + Canopy+Moss_cent +(1|Block),family="binomial", data=tshegermdat)
standorgcanund<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint1<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint2<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint3<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint4<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint5<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint6<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
standorgcanundint7<-glmmadmb(tshey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block),family="binomial", data=tshegermdat)
Stand<-anova(nullgermod.tshe,stand.tshe)$Deviance[2]
Org<-anova(stand.tshe,standorg.tshe)$Deviance[2]
Can<-anova(standorg.tshe,standorgcan)$Deviance[2]
Und<-anova(standorgcan,standorgcanund)$Deviance[2]
storgint<-anova(standorgcanund,standorgcanundint1)$Deviance[2]
stcanint<-anova(standorgcanundint1,standorgcanundint2)$Deviance[2]
stundint<-anova(standorgcanundint2,standorgcanundint3)$Deviance[2]
canorg<-anova(standorgcanundint3,standorgcanundint4)$Deviance[2]
undorg<-anova(standorgcanundint4,standorgcanundint5)$Deviance[2]
canund<-anova(standorgcanundint5,standorgcanundint6)$Deviance[2]
stcanund<-anova(standorgcanundint6,standorgcanundint7)$Deviance[2]
#threewayint<-anova(standorgcanundint6,germod.tshe.moss)$Deviance[2]#should be the same as the above and it is
tsheanova<-Anova(germod.tshe.moss)
tsheanova<-tsheanova[-5,]
anovatab.tshe<-cbind(row.names(tsheanova),c(Org,Can,Stand,Und,storgint,stcanint,stundint,canorg,undorg,canund,stcanund))
colnames(anovatab.tshe)<-c("Source","Deviance")
PropDev<-as.numeric(anovatab.tshe[,2])/sum(as.numeric(anovatab.tshe[,2]))
anovatab.tshe<-cbind(anovatab.tshe,PropDev)
write.csv(anovatab.tshe,"analyses/tshegerm.propdev.csv",row.names = FALSE)

###Model comparison with continuous microclimate variables
climdat<-read.csv("data/AllStands_clim1.csv", header=TRUE)
climdat$Block<-as.factor(climdat$Block)
#Combine the climate data with the transplant data columns that we want to use for survival and growth analysis
germdatsub<-subset(dat2011,select=c("SpPlant","Stand","Block","Origin","Canopy","Understory","SeedsAdded","TotalFails"))
climdatsub<-subset(climdat,select=c("Stand","Block","Elevation_m","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
colnames(climdatsub)[1]<-"Standname"
colnames(climdatsub)[3]<-"Stand"
alldat <- join(germdatsub, climdatsub, by=c("Stand","Block","Canopy","Understory"), match="first")
alldat<-subset(alldat,select=c("SpPlant","Stand","Block","Origin","Canopy","Understory","SeedsAdded","TotalFails","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
#Ok, so now I have climate data and germination in the same place. 
#separate out species
alldat$snow_cover_duration=as.numeric(alldat$snow_cover_duration)
alldat$GDD_total=as.numeric(alldat$GDD_total)
alldat$meanGST=as.numeric(alldat$meanGST)
alldat$meanGST=as.numeric(alldat$meanGST)
alldat$Light_Mean=as.numeric(alldat$Light_Mean)
#abamgermdat<-alldat[alldat$SpPlant=="ABAM",]
tsmegermdat<-alldat[alldat$SpPlant=="TSME",]
#tshegermdat<-alldat[alldat$SpPlant=="TSHE",]
tsmegermdat$Stand<-as.factor(tsmegermdat$Stand)
tsmegermdat$Block<-as.factor(tsmegermdat$Block)
tsmegermdat$snow_cover_duration_cent=scale(as.numeric(tsmegermdat$snow_cover_duration))
#tshegermdat$snow_cover_duration_cent=scale(as.numeric(tshegermdat$snow_cover_duration))
#abamgermdat$snow_cover_duration_cent=scale(as.numeric(abamgermdat$snow_cover_duration))
tsmegermdat$GDD_total_cent=scale(as.numeric(tsmegermdat$GDD_total))
#tshegermdat$GDD_total_cent=scale(as.numeric(tshegermdat$GDD_total))
#abamgermdat$GDD_total_cent=scale(as.numeric(abamgermdat$GDD_total))
tsmegermdat$meanGST_cent=scale(as.numeric(tsmegermdat$meanGST))
#tshegermdat$meanGST_cent=scale(as.numeric(tshegermdat$meanGST))
#abamgermdat$meanGST_cent=scale(as.numeric(abamgermdat$meanGST))
tsmegermdat$Light_Mean_cent=scale(as.numeric(tsmegermdat$Light_Mean))
#tshegermdat$Light_Mean_cent=scale(as.numeric(tshegermdat$Light_Mean))
#abamgermdat$Light_Mean_cent=scale(as.numeric(abamgermdat$Light_Mean))
#Now climate models
tsmegermdat<-na.omit(tsmegermdat)
tsmegermdat$tsmey<-cbind(tsmegermdat$SeedsAdded,tsmegermdat$TotalFails)
consgermod.tsme<-glmmadmb(tsmey ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
#gddgermod.abam<-glmmadmb(abamy~Origin+GDD_total_cent+Origin:GDD_total_cent+Moss_cent+(1|Block),data=abamgermdat, family="binomial")
gddgermod.tsme<-glmmadmb(tsmey~Origin+GDD_total_cent+Origin:GDD_total_cent+(1|Block),data=tsmegermdat, family="binomial")
#gddgermod.tshe<-glmmadmb(tshey~Origin+GDD_total_cent+Origin:GDD_total_cent+Moss_cent+(1|Block),data=tshegermdat, family="binomial")
#gstgermod.abam<-glmmadmb(abamy~Origin+meanGST_cent+Origin:meanGST_cent+Moss_cent+(1|Block),data=abamgermdat, family="binomial")
gstgermod.tsme<-glmmadmb(tsmey~Origin+meanGST_cent+Origin:meanGST_cent+(1|Block), data=tsmegermdat, family="binomial")
#gstgermod.tshe<-glmmadmb(tshey~Origin+meanGST_cent+Origin:meanGST_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
#sdgermod.abam<-glmmadmb(abamy~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
sdgermod.tsme<-glmmadmb(tsmey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+(1|Block), data=tsmegermdat, family="binomial")
#sdgermod.tshe<-glmmadmb(tshey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
#summary(sdgermod.tshe)
#light models
#lightgermod.abam<-glmmadmb(abamy~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
lightgermod.tsme<-glmmadmb(tsmey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
#lightgermod.tshe<-glmmadmb(tshey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

#Try replacing canopy and understory with "light" AND planted stand with GDD- this might tell us how much competition is due to light
#constgddlightgermod.abam<-glmmadmb(abamy~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constgddlightgermod.tsme<-glmmadmb(tsmey~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
#constgddlightgermod.tshe<-glmmadmb(tshey~GDD_total_cent+Origin+Light_Mean_cent+GDD_total_cent:Origin+GDD_total_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")

#constsdlightgermod.abam<-glmmadmb(abamy~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constsdlightgermod.tsme<-glmmadmb(tsmey~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
#constsdlightgermod.tshe<-glmmadmb(tshey~snow_cover_duration_cent+Origin+Light_Mean_cent+snow_cover_duration_cent:Origin+snow_cover_duration_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
#Try replacing canopy and understory with "light" AND planted stand with meanGST- this might tell us how much competition is due to light
#constgstlightgermod.abam<-glmmadmb(abamy~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
constgstlightgermod.tsme<-glmmadmb(tsmey~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
#constgstlightgermod.tshe<-glmmadmb(tshey~meanGST_cent+Origin+Light_Mean_cent+meanGST_cent:Origin+meanGST_cent:Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
nullgermod.tsme<-glmmadmb(tsmey~1+(1|Block), data=tsmegermdat, family="binomial")

#AICctab(constgermod.abam,sdgermod.abam,lightgermod.abam,gstgermod.abam,gddgermod.abam,constgstlightgermod.abam,constgddlightgermod.abam,constsdlightgermod.abam)#lowest aic for constgddhimod.abam
AICctab(nullgermod.tsme,consgermod.tsme,sdgermod.tsme,lightgermod.tsme,gstgermod.tsme,gddgermod.tsme,constgstlightgermod.tsme,constgddlightgermod.tsme,constsdlightgermod.tsme)#lowest aic for constsdhimod.tsme 
#AICctab(constgermod.tshe,sdgermod.tshe,lightgermod.tshe,gstgermod.tshe,gddgermod.tshe,constgstlightgermod.tshe,constgddlightgermod.tshe,constsdlightgermod.tshe)#lowest aic for constgddhimod.tshe 
summary(sdgermod.tsme)
anova(consgermod.tsme)
#model fitting for prop dev
consno3mod.tsme<-glmmadmb(tsmey ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
#TSME
setwd("~/Dropbox/Documents/Work/UW/Research/Mount Rainier/2014Germination(2012data)")
dat<-read.csv("MORAGermData20112012.csv", header=TRUE)
head(dat)
dim(dat)
dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms=50#replace any success that were 
##First, select only rows in which seeds were added
addat=dat[dat$SeedsAdded>0,]
dim(addat)#1281 plots...strange number...
#Next, let's get some summary stats on the data and look at i
#select out 2011 dat- this is most complete.
dat2011<-addat[addat$Year=="2011",]#
dim(dat2011)#820 rows

tsmegermdat2<-subset(tsmegermdat,selec=c(TotalGerms,TotalFails,Block,Stand,Origin,Canopy,Understory))
tsmegermdat3<-na.omit(tsmegermdat2)
consgermmod.tsme<-glmmadmb(tsmey ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat3)#now merge in microclimate variables
climdat<-read.csv("data/AllStands_clim1.csv", header=TRUE)
climdat$Block<-as.factor(climdat$Block)
#Combine the climate data with the transplant data columns that we want to use for survival and growth analysis
germdatsub<-subset(dat2011,select=c("SpPlant","Stand","Block","Origin","Canopy","Understory","SeedsAdded","TotalFails"))
climdatsub<-subset(climdat,select=c("Stand","Block","Elevation_m","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
colnames(climdatsub)[1]<-"Standname"
colnames(climdatsub)[3]<-"Stand"
alldat <- join(germdatsub, climdatsub, by=c("Stand","Block","Canopy","Understory"), match="first")
alldat<-subset(alldat,select=c("SpPlant","Stand","Block","Origin","Canopy","Understory","SeedsAdded","TotalFails","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
#Ok, so now I have climate data and germination in the same place. 
#separate out species
alldat$snow_cover_duration=as.numeric(alldat$snow_cover_duration)
alldat$GDD_total=as.numeric(alldat$GDD_total)
alldat$meanGST=as.numeric(alldat$meanGST)
alldat$meanGST=as.numeric(alldat$meanGST)
alldat$Light_Mean=as.numeric(alldat$Light_Mean)
tsmegermdat<-alldat[alldat$SpPlant=="TSME",]
#tshegermdat<-alldat[alldat$SpPlant=="TSHE",]
####Germination and "Seed bed" analysis for Reviewer 2. 
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
#tapply(grcov2$Moss,list(grcov2$Stand,grcov2$Canopy), mean, na.rm=T)
###Problem: moss cover is higher in gaps!
#clearly substrate has been shown to matter- 
#does it change the affect sign of the 
germdat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
germdat[which(germdat$TotalGerms>50 &germdat$SeedsAdded==50),]$TotalGerms=50
##First, select only rows in which seeds were added
addat<-germdat[germdat$SeedsAdded>0,]
#select out 2011 data for analysis
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
#make sure that datasets are the same, regardless of explanatory variables
abamgermdat<- abamgermdat[apply(abamgermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)

#TSME
#make sure that datasets are the same, regardless of explanatory variables
tsmegermdat<-tsmegermdat[apply(tsmegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)

#Tshe
tshegermdat<- tshegermdat[apply(tshegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

germod.abam<-glmmadmb(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=abamgermdat, family="binomial")
#germod.tsme<-glmmadmb(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="binomial")
germod.tshe<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="binomial")
anova(germod.tsme)
#Get columns of grcov for Canopy and Block to match format of germdat so that i can merge in the moss
colnames(grcov2)[3]<-"BlockNum"
grcov2$Block<-paste(grcov2$Stand,grcov2$BlockNum,sep="")
grcov2$Canopy<-NA
grcov2[grcov2$Gap.=="Y",]$Canopy<-"CompAbsent"
grcov2[grcov2$Gap.=="N",]$Canopy<-"CompPresent"
colnames(germdat)[2]<-"Elev"
germdat2<-join(germdat, grcov2, by=c("Block","Plot","Canopy"), match="all")
dat<-subset(germdat2,select=c("Stand","Block","Plot","Canopy","Understory","SpPlant","Origin","Year","TotalGerms","SeedsAdded","TotalFails","Litter","Moss","Wood"))
dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms=50
##First, select only rows in which seeds were added
addat<-dat[dat$SeedsAdded>0,]
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

#Tshe
tshegermdat<- tshegermdat[apply(tshegermdat, 1, function(x) all(!is.na(x))),] # only keep rows of all not na
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

#germod.abam<-glmmadmb(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=abamgermdat, family="binomial")
#summary(germod.abam)
abamgermdat$Moss_cent<-scale(abamgermdat$Moss)
germod.abam.moss<-glmmadmb(abamy~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
summary(germod.abam.moss)
#Anova(germod.abam.moss,type="III")
#AIC(germod.abam,germod.abam.moss)
#so, include moss for ABAM
#germod.tsme<-glmmadmb(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="binomial")
#summary(germod.tsme)
#Anova(germod.tsme,type="III")
#tsmegermdat$Moss_cent<-scale(tsmegermdat$Moss)
#germod.tsme.moss<-glmmadmb(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss+(1|Block), data=tsmegermdat, family="binomial")
#summary(germod.tsme.moss)
#Anova(germod.tsme.moss,type="III")
#AIC(germod.tsme,germod.tsme.moss)
#anova(germod.tsme,germod.tsme.moss)
#do not include moss for TSME- higher AIC
germod.tshe<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="binomial")
#summary(germod.tshe)
Anova(germod.tshe,type="III")
tshegermdat$Moss_cent<-scale(tshegermdat$Moss)
germod.tshe.moss<-glmmadmb(tshey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss+(1|Block), data=tshegermdat, family="binomial")
#summary(germod.tshe.moss)
Anova(germod.tshe.moss,type="III")
AIC(germod.tshe,germod.tshe.moss)
anova(germod.tshe,germod.tshe.moss)
#include moss for tshe model
###now compare these categorical models with continuous models
colnames(climdatsub)[3]<-"Stand"
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

summary(sdgermod.tshe)
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

AICctab(constgermod.abam,sdgermod.abam,lightgermod.abam,gstgermod.abam,gddgermod.abam,constgstlightgermod.abam,constgddlightgermod.abam,constsdlightgermod.abam)#lowest aic for constgddhimod.abam
AICctab(constgermod.tsme,sdgermod.tsme,lightgermod.tsme,gstgermod.tsme,gddgermod.tsme,constgstlightgermod.tsme,constgddlightgermod.tsme,constsdlightgermod.tsme)#lowest aic for constsdhimod.tsme 
AICctab(constgermod.tshe,sdgermod.tshe,lightgermod.tshe,gstgermod.tshe,gddgermod.tshe,constgstlightgermod.tshe,constgddlightgermod.tshe,constsdlightgermod.tshe)#lowest aic for constgddhimod.tshe 


##Greenhouse germination
ggermdat<-read.csv("data/AllGreenhouseGerm.csv", head=T)
ggermdat$grate<-ggermdat$SeedsGerm/ggermdat$SeedsSown
ggermdat$SeedsFail<-ggermdat$SeedsSown-ggermdat$SeedsGerm
ggermdat$yrsstore<-ggermdat$YearSown-ggermdat$YearCollected# number of years stored before sowing
ggermdat$Flat2<-paste(ggermdat$YearSown,ggermdat$Flat,sep=".")
ggermdat$YearCollected<-as.factor(ggermdat$YearCollected)
tsmeggerm<-ggermdat[ggermdat$Species=="TSME",]
tsheggerm<-ggermdat[ggermdat$Species=="TSHE",]
abamggerm<-ggermdat[ggermdat$Species=="ABAM",]
#models for supplemental
gmod.abam<-glm(cbind(abamggerm$SeedsGerm,abamggerm$SeedsFail)~SourceStand+YearCollected+yrsstore, data=abamggerm, family=binomial)
summary(gmod.abam)
Anova(gmod.abam)
gmod.tsme<-glm(cbind(tsmeggerm$SeedsGerm,tsmeggerm$SeedsFail)~SourceStand+YearCollected+yrsstore, data=tsmeggerm, family=binomial)
summary(gmod.tsme)
anova(gmod.tsme)
gmod.tshe<-glm(cbind(tsheggerm$SeedsGerm,tsheggerm$SeedsFail)~SourceStand+YearCollected+yrsstore, data=tsheggerm, family=binomial)
summary(gmod.tshe)
anova(gmod.tshe)
#years stored has negative effect on germination for TSHE but no other species
ggermdat09<-ggermdat[which(ggermdat$YearCollected==2009),]
ggermdat10<-ggermdat[which(ggermdat$YearCollected==2010),]
ggermdat09[which(ggermdat09$Elevation=="1191" & ggermdat09$Species=="ABAM"),]$Elevation<-1069
gmean09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), mean,na.rm=T)
gmean10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), mean,na.rm=T)
gsd09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), sd,na.rm=T)
gsd10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), sd,na.rm=T)
gn09<-tapply(ggermdat09$grate,list(ggermdat09$Species,ggermdat09$Elevation,ggermdat09$YearSown), length)
gn10<-tapply(ggermdat10$grate,list(ggermdat10$Species,ggermdat10$Elevation,ggermdat10$YearSown), length)
gse09<-gsd09/sqrt(gn09)
gse10<-gsd10/sqrt(gn10)

#now look at field germination data
addat$grate<-addat$TotalGerms/addat$SeedsAdded
fieldgmean<-tapply(addat$grate,list(addat$SpPlant,addat$Origin,addat$Year), mean,na.rm=T)
fieldsd<-tapply(addat$grate,list(addat$SpPlant,addat$Origin,addat$Year), sd,na.rm=T)
fieldn<-tapply(addat$grate,list(addat$SpPlant,addat$Origin,addat$Year), length)
fieldse<-fieldsd/sqrt(fieldn)


#Figure of germination in the greenhouse, for supplemental
quartz(height=7,width=7)
par(mfcol=c(3,2),mai=c(.6,.7,.2,.1), omi=c(.7,.01,.2,.1))
X<-c(1,2,3,4,5)
#TSME
plot(gmean09[3,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c(NA,NA,"black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Greenhouse",side=3,line=1.5, adj=0)
mtext("Tsuga mertensiana",side=3,line=.2, adj=0, font=3)
mtext("a",side=3,line=.2, adj=-.1)
for (i in 1:5){
  arrows(X[i],gmean09[3,i,]-gse09[3,i,],X[i],gmean09[3,i,]+gse09[3,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[3,i,1]-gse10[3,i,1],X[i],gmean10[3,i,1]+gse10[3,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[3,i,2]-gse10[3,i,2],X[i],gmean10[3,i,2]+gse10[3,i,2],length=.05,angle=90,code=0)}
points(gmean09[3,,]~X,pch=21,cex=1.5,bg=c(NA,NA,"black","gray","white"))
points(gmean10[3,,1]~X,pch=22,cex=1.5,bg=c(NA,NA,"black","gray","white"))
points(gmean10[3,,2]~X,pch=23,cex=1.5,bg=c(NA,NA,"black","gray","white"))
legend(1,.6,legend=c("2009","2010-1", "2010-2"),bty="n",pch=c(21,22,23),pt.bg=c("black"),angle=45,cex=1.1, pt.cex=1.5)

##ABAM
plot(gmean09[1,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c(NA,NA,"black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=.2, adj=0, font=3)
mtext("b",side=3,line=.2, adj=-.1)
for (i in 1:5){
  arrows(X[i],gmean09[1,i,]-gse09[1,i,],X[i],gmean09[1,i,]+gse09[1,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[1,i,1]-gse10[1,i,1],X[i],gmean10[1,i,1]+gse10[1,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[1,i,2]-gse10[1,i,2],X[i],gmean10[1,i,2]+gse10[1,i,2],length=.05,angle=90,code=0)}
points(gmean09[1,,]~X,pch=21,cex=1.5,bg=c("black","gray",NA,NA,"white"))
points(gmean10[1,,1]~X,pch=22,cex=1.5,bg=c("black","gray",NA,NA,"white"))
points(gmean10[1,,2]~X,pch=23,cex=1.5,bg=c("black","gray",NA,NA,"white"))
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
points(gmean09[2,,]~X,pch=21,cex=1.5,bg=c("gray",NA,"white",NA,NA))
points(gmean10[2,,1]~X,pch=22,cex=1.5,bg=c("gray",NA,"white",NA,NA))
points(gmean10[2,,2]~X,pch=23,cex=1.5,bg=c("gray",NA,"white",NA,NA))
mtext("Tsuga heterophylla",side=3,line=.2, adj=0, font=3)
mtext("c",side=3,line=.2, adj=-.1)
mtext ("Elevation of seed origin population (m)", line=-13, cex=.9)
axis(1, at = c(1,2,3,4,5), labels = colnames(gmean09), tick = FALSE, cex.axis=1.1, line=-.5)

#TSME
plot(fieldgmean[3,,1]~X[3:5],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Tsuga mertensiana",side=3,line=.2, adj=0, font=3)
mtext("d",side=3,line=.2, adj=-.1)
mtext("Mt. Rainier",side=3,line=1.5, adj=0)
arrows(X[3:5],fieldgmean[3,,1]-fieldse[3,,1],X[3:5],fieldgmean[3,,1]+fieldse[3,,1],length=.05,angle=90,code=0)
arrows(X[3:5],fieldgmean[3,,2]-fieldse[3,,2],X[3:5],fieldgmean[3,,2]+fieldse[3,,2],length=.05,angle=90,code=0)
points(fieldgmean[3,,1]~X[3:5],pch=21,cex=1.5,bg=c("black","gray","white"))
points(fieldgmean[3,,2]~X[3:5],pch=22,cex=1.5,bg=c("black","gray","white"))

##ABAM
elevs=c(1,2,5)

plot(fieldgmean[1,,1]~X[elevs],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=.2, adj=0, font=3)
mtext("e",side=3,line=.2, adj=-.1)
arrows(X[elevs],fieldgmean[1,,1]-fieldse[1,,1],X[elevs],fieldgmean[1,,1]+fieldse[1,,1],length=.05,angle=90,code=0)
arrows(X[elevs],fieldgmean[1,,2]-fieldse[1,,2],X[elevs],fieldgmean[1,,2]+fieldse[1,,2],length=.05,angle=90,code=0)

points(fieldgmean[1,,1]~X[elevs],pch=21,cex=1.5,bg=c("black","gray","white"))
points(fieldgmean[1,,2]~X[elevs],pch=22,cex=1.5,bg=c("black","gray","white"))

##TSHE
elevs=c(NA,1,3)
plot(fieldgmean[2,,1]~X[elevs],ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c(NA,"gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
  arrows(X[elevs],fieldgmean[2,,1]-fieldse[2,,1],X[elevs],fieldgmean[2,,1]+fieldse[2,,1],length=.05,angle=90,code=0)
  arrows(X[elevs],fieldgmean[2,,2]-fieldse[2,,2],X[elevs],fieldgmean[2,,2]+fieldse[2,,2],length=.05,angle=90,code=0)
points(fieldgmean[2,,1]~X[elevs],pch=21,cex=1.5,bg=c(NA,"gray","white"))
points(fieldgmean[2,,2]~X[elevs],pch=22,cex=1.5,bg=c(NA,"gray","white"))
mtext("Tsuga heterophylla",side=3,line=.2, adj=0, font=3)
mtext("f",side=3,line=.2, adj=-.1)
mtext ("Elevation of seed origin population (m)", line=-13, cex=.9)
axis(1, at = c(1,2,3,4,5), labels = colnames(gmean09), tick = FALSE, cex.axis=1.1, line=-.5)
