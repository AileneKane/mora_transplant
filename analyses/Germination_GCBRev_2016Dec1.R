#New germination models
rm(list=ls()) 
options(stringsAsFactors=FALSE)
library(dplyr)
library(plyr)
library(glmmADMB)
library(car)
setwd("~/git/mora_transplant")
dat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
head(dat);tail(dat)
dim(dat)#2939  13
dat$TotalGerms<-as.factor(dat$TotalGerms)
#dat$SeedsAdded<-as.numeric(dat$SeedsAdded)
#dat[which(dat$TotalGerms>50&dat$SeedsAdded==50),]$TotalGerms<-50#replace any success that were 
##First, select only rows in which seeds were added
addat<-dat[dat$SeedsAdded>0,]
#head(addat);tail(addat)
#dim(addat)#1080 plots
dat2011<-addat[addat$Year=="2011",]#
dat2012<-addat[addat$Year=="2012",]#
dim(dat2011)#620 plots
dim(dat2012)#460 plots
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
grcov2[which(grcov2$Plot=="I (formerly E)"),]$Plot<-"I"
colnames(grcov2)[3]<-"BlockNum"
grcov2$Block<-paste(grcov2$Stand,grcov2$BlockNum,sep="")
grcov2$Canopy<-NA
grcov2[grcov2$Gap.=="Y",]$Canopy<-"CompAbsent"
grcov2[grcov2$Gap.=="N",]$Canopy<-"CompPresent"
grcov3<-grcov2[-which(is.na(grcov2$Litter)),]
colnames(dat2011)[2]<-"Elev"
colnames(dat)[2]<-"Elev"
germdat2<-join(dat2011, grcov3, by=c("Block","Plot","Canopy"), match="all")
dim(germdat2[which(is.na(germdat2$Moss)),])#6 rows missing because of missing moss data for 2011 data; 354 missing for both years 
#the following sites are missing moss data:"AE101" "AE102" "AE103" "AE104" "AV064" "AV065" "HIGH1" "HIGH2"
#[9] "HIGH3" "HIGH5"
#rows 1,2,6,7,21:29,33:34,,39,41,312:313 missing moss, 335:336,382:383
dat2<-subset(germdat2,select=c("Stand","Elev","Block","Plot","Canopy","Understory","SpPlant","Origin","Year","TotalGerms","SeedsAdded","TotalFails","Litter","Moss","Wood"))
abamgermdat<-dat2[dat2$SpPlant=="ABAM",]#200 rows
abamgermdat<-na.omit(abamgermdat)#417 rows
abamgermdat$Stand<-as.factor(abamgermdat$Stand)
abamgermdat$Origin<-as.factor(abamgermdat$Origin)
abamgermdat$TotalGerms<-as.numeric(abamgermdat$TotalGerms)
abamgermdat$Block<-as.factor(abamgermdat$Block)
abamgermdat$Moss<-as.numeric(abamgermdat$Moss)
abamgermdat$Moss_cent<-scale(abamgermdat$Moss)
abamy=cbind(abamgermdat$TotalGerms,abamgermdat$TotalFails)
germod.abam.moss<-glmmadmb(abamy~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="binomial")
germod.abam<-glmmadmb(abamy~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=abamgermdat, family="binomial")
abamgermdat$TotalGerms<-as.numeric(abamgermdat$TotalGerms)
germod.abam.zip.moss<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=abamgermdat, family="poisson",zeroInflation = TRUE)
germod.abam.zip<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=abamgermdat, family="poisson",zeroInflation = TRUE)
AIC(germod.abam.moss,germod.abam,germod.abam.zip.moss,germod.abam.zip)
#zero-inflated poisson with could not be fit, but binomial model with moss wins
summary(germod.abam.moss)
Anova(germod.abam.moss, type="III")
Anova(germod.abam.zip.moss, type="III")

anova(germod.abam.zip,germod.abam.zip.moss)
anova(germod.abam,germod.abam.moss)
anova(germod.abam.moss,germod.abam.zip.moss)
###Add in microclimate data

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
tshegermdat<-dat2[dat2$SpPlant=="TSHE",]
tshegermdat<-na.omit(tshegermdat)
tshegermdat$Stand=factor(tshegermdat$Stand)
tshegermdat$Origin=factor(tshegermdat$Origin)
tshegermdat$TotalGerms=as.numeric(tshegermdat$TotalGerms)
tshegermdat$Block=as.factor(tshegermdat$Block)
tshegermdat$Moss=as.numeric(tshegermdat$Moss)
tshegermdat$Moss_cent<-scale(tshegermdat$Moss)
tshey=cbind(tshegermdat$TotalGerms,tshegermdat$TotalFails)

germod.tshe.moss<-glmmadmb(tshey~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
germod.tshe<-glmmadmb(tshey~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="binomial")
germod.tshe.zip.moss<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tshegermdat, family="poisson",zeroInflation = TRUE)
germod.tshe.zip<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tshegermdat, family="poisson",zeroInflation = TRUE)
AIC(germod.tshe.zip.moss,germod.tshe.zip,germod.tshe.moss,germod.tshe)
#binomial with moss wins! use this one for further analysis
summary(germod.tshe.moss)
Anova(germod.tshe.moss)
anova(germod.tshe.zip,germod.tshe.zip.moss)
anova(germod.tshe,germod.tshe.moss)
anova(germod.tshe.moss,germod.tshe.zip.moss)
#Add microclimate
tshegermdat <- join(tshegermdat, climdatsub3, by=c("Block","Canopy","Understory"), match="first")
#separate out species
tshegermdat$snow_cover_duration=as.numeric(tshegermdat$snow_dur)
tshegermdat$GDD_total=as.numeric(tshegermdat$GDD_total)
tshegermdat$meanGST=as.numeric(tshegermdat$meanGST)
tshegermdat$MAT=as.numeric(tshegermdat$MAT)
tshegermdat$Light_GDD=as.numeric(tshegermdat$Light_GDD)
tshegermdat$Light_Mean=as.numeric(tshegermdat$Light_Mean)
tshegermdat$Stand<-as.factor(tshegermdat$Stand)
tshegermdat$Block<-as.factor(tshegermdat$Block)
tshegermdat$snow_cover_duration_cent=scale(as.numeric(tshegermdat$snow_cover_duration))
tshegermdat$GDD_total_cent=scale(as.numeric(tshegermdat$GDD_total))
tshegermdat$meanGST_cent=scale(as.numeric(tshegermdat$meanGST))
tshegermdat$Light_Mean_cent=scale(as.numeric(tshegermdat$Light_Mean))
tshegermdat$Light_GDD_cent=scale(as.numeric(tshegermdat$Light_GDD))
tshegermdat$MAT=scale(as.numeric(tshegermdat$MAT))

#Now climate models
tshegermdat<-na.omit(tshegermdat)
tshegermdat$SeedsAdded<-as.factor(tshegermdat$SeedsAdded)
tshegermdat$tshey<-cbind(tshegermdat$SeedsAdded,tshegermdat$TotalFails)
consgermod.tshe<-glmmadmb(tshey ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tshegermdat)
gddgermod.tshe<-glmmadmb(tshey~Origin+GDD_total_cent+Origin:GDD_total_cent+(1|Block),data=tshegermdat, family="binomial")
gstgermod.tshe<-glmmadmb(tshey~Origin+meanGST_cent+Origin:meanGST_cent+(1|Block), data=tshegermdat, family="binomial")
lightgermod.tshe<-glmmadmb(tshey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tshegermdat, family="binomial")
light2germod.tshe<-glmmadmb(tshey~Origin+Light_GDD_cent+Origin:Light_GDD_cent+(1|Block), data=tshegermdat, family="binomial")
matgermod.tshe<-glmmadmb(tshey~Origin+MAT+Origin:MAT+(1|Block), data=tshegermdat, family="binomial")
sdgermod.tshe<-glmmadmb(tshey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+(1|Block), data=tshegermdat, family="binomial")
#interaction models
gddlightgermod.tshe<-glmmadmb(tshey~Origin+GDD_total_cent+Light_GDD_cent+Origin:GDD_total_cent+Light_GDD_cent:Origin+Light_GDD_cent:GDD_total_cent+(1|Block),data=tshegermdat, family="binomial")
gstlightgermod.tshe<-glmmadmb(tshey~Origin+meanGST_cent+Light_GDD_cent+Origin:meanGST_cent+Light_GDD_cent:Origin+Light_GDD_cent:meanGST_cent+(1|Block), data=tshegermdat, family="binomial")
matlightgermod.tshe<-glmmadmb(tshey~Origin+MAT+Light_GDD_cent+Origin:Light_GDD_cent+Origin:MAT+Light_GDD_cent:MAT+(1|Block), data=tshegermdat, family="binomial")
sdlightgermod.tshe<-glmmadmb(tshey~Origin+snow_cover_duration_cent+Light_GDD_cent+Light_GDD_cent:Origin+Origin:snow_cover_duration_cent+Light_GDD_cent:snow_cover_duration+(1|Block), data=tshegermdat, family="binomial")
nullgermod.tshe<-glmmadmb(tshey~1+Moss_cent+(1|Block), data=tshegermdat, family="binomial")
summary(nullgermod.tshe)
AICctab(nullgermod.tshe,consgermod.tshe,gddgermod.tshe,gstgermod.tshe,lightgermod.tshe,light2germod.tshe,matgermod.tshe,sdgermod.tshe,gddlightgermod.tshe,gstlightgermod.tshe,matlightgermod.tshe,sdlightgermod.tshe)
#



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

#TSME
tsmegermdat<-dat2[dat2$SpPlant=="TSME",]
tsmegermdat<-na.omit(tsmegermdat)
tsmegermdat$Stand=factor(tsmegermdat$Stand)
tsmegermdat$Origin=factor(tsmegermdat$Origin)
tsmegermdat$TotalGerms=as.numeric(tsmegermdat$TotalGerms)
tsmegermdat$Block=as.factor(tsmegermdat$Block)
tsmegermdat$Moss=as.numeric(tsmegermdat$Moss)
tsmegermdat$Moss_cent<-scale(tsmegermdat$Moss)
tsmey=cbind(tsmegermdat$TotalGerms,tsmegermdat$TotalFails)

germod.tsme.moss<-glmmadmb(tsmey~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tsmegermdat, family="binomial")
germod.tsme<-glmmadmb(tsmey~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="binomial")
germod.tsme.zip.moss<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent+(1|Block), data=tsmegermdat, family="poisson",zeroInflation = TRUE)
germod.tsme.zip<-glmmadmb(TotalGerms~Origin+Canopy+Stand+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block), data=tsmegermdat, family="poisson",zeroInflation = TRUE)
AIC(germod.tsme.zip.moss,germod.tsme.zip,germod.tsme.moss,germod.tsme)
#binomial with moss wins! use this one for further analysis
summary(germod.tsme)
Anova(germod.tsme, type="III")
anova(germod.tsme.zip,germod.tsme.zip.moss)
anova(germod.tsme,germod.tsme.moss)
anova(germod.tsme.moss,germod.tsme.zip.moss)
nullgermod.tsme<-glmmadmb(tsmey~(1|Block), data=tsmegermdat, family="binomial")
stand.tsme<-glmmadmb(tsmey ~ Stand+(1|Block),family="binomial", data=tsmegermdat)
standorg.tsme<-glmmadmb(tsmey ~ Stand  +Origin+(1|Block),family="binomial", data=tsmegermdat)
standorgcan<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +(1|Block),family="binomial", data=tsmegermdat)
standorgcanund<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint1<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint2<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint3<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint4<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint5<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint6<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
standorgcanundint7<-glmmadmb(tsmey ~ Stand  +Origin + Canopy +Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
Stand<-anova(nullgermod.tsme,stand.tsme)$Deviance[2]
Org<-anova(stand.tsme,standorg.tsme)$Deviance[2]
Can<-anova(standorg.tsme,standorgcan)$Deviance[2]
Und<-anova(standorgcan,standorgcanund)$Deviance[2]
storgint<-anova(standorgcanund,standorgcanundint1, type="III")$Deviance[2]
stcanint<-anova(standorgcanundint1,standorgcanundint2)$Deviance[2]
stundint<-anova(standorgcanundint2,standorgcanundint3)$Deviance[2]
canorg<-anova(standorgcanundint3,standorgcanundint4)$Deviance[2]
undorg<-anova(standorgcanundint4,standorgcanundint5)$Deviance[2]
canund<-anova(standorgcanundint5,standorgcanundint6)$Deviance[2]
stcanund<-anova(standorgcanundint6,standorgcanundint7)$Deviance[2]
#threewayint<-anova(standorgcanundint6,germod.tsme.moss)$Deviance[2]#should be the same as the above and it is
tsmeanova<-Anova(germod.tsme)
anovatab.tsme<-cbind(row.names(tsmeanova),c(Org,Can,Stand,Und,storgint,stcanint,stundint,canorg,undorg,canund,stcanund))
colnames(anovatab.tsme)<-c("Source","Deviance")
PropDev<-as.numeric(anovatab.tsme[,2])/sum(as.numeric(anovatab.tsme[,2]))
anovatab.tsme<-cbind(anovatab.tsme,PropDev)
write.csv(anovatab.tsme,"analyses/tsmegerm.propdev.csv",row.names = FALSE)
#
###Model comparison with continuous microclimate variables

#add microclimate data
tsmegermdat <- join(tsmegermdat, climdatsub3, by=c("Block","Canopy","Understory"), match="first")
#alldat<-subset(alldat,select=c("SpPlant","Stand","Block","Origin","Canopy","Understory","SeedsAdded","TotalFails","Canopy","Understory","snow_appearance_date","snow_disappearance_date","snow_cover_duration","meanGST","GDD_total","Light_Mean"))
#Ok, so now I have climate data and germination in the same place. 
#separate out species
tsmegermdat$snow_cover_duration=as.numeric(tsmegermdat$snow_dur)
tsmegermdat$GDD_total=as.numeric(tsmegermdat$GDD_total)
tsmegermdat$meanGST=as.numeric(tsmegermdat$meanGST)
tsmegermdat$MAT=as.numeric(tsmegermdat$MAT)
tsmegermdat$Light_GDD=as.numeric(tsmegermdat$Light_GDD)
tsmegermdat$Light_Mean=as.numeric(tsmegermdat$Light_Mean)
tsmegermdat$Stand<-as.factor(tsmegermdat$Stand)
tsmegermdat$Block<-as.factor(tsmegermdat$Block)
tsmegermdat$snow_cover_duration_cent=scale(as.numeric(tsmegermdat$snow_cover_duration))
tsmegermdat$GDD_total_cent=scale(as.numeric(tsmegermdat$GDD_total))
tsmegermdat$meanGST_cent=scale(as.numeric(tsmegermdat$meanGST))
tsmegermdat$Light_Mean_cent=scale(as.numeric(tsmegermdat$Light_Mean))
tsmegermdat$Light_GDD_cent=scale(as.numeric(tsmegermdat$Light_GDD))
tsmegermdat$MAT=scale(as.numeric(tsmegermdat$MAT))

#Now climate models
tsmegermdat<-na.omit(tsmegermdat)
tsmegermdat$SeedsAdded<-as.numeric(tsmegermdat$SeedsAdded)
tsmegermdat$tsmey<-cbind(tsmegermdat$SeedsAdded,tsmegermdat$TotalFails)
consgermod.tsme<-glmmadmb(tsmey ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+(1|Block),family="binomial", data=tsmegermdat)
gddgermod.tsme<-glmmadmb(tsmey~Origin+GDD_total_cent+Origin:GDD_total_cent+(1|Block),data=tsmegermdat, family="binomial")
gstgermod.tsme<-glmmadmb(tsmey~Origin+meanGST_cent+Origin:meanGST_cent+(1|Block), data=tsmegermdat, family="binomial")
lightgermod.tsme<-glmmadmb(tsmey~Origin+Light_Mean_cent+Origin:Light_Mean_cent+(1|Block), data=tsmegermdat, family="binomial")
light2germod.tsme<-glmmadmb(tsmey~Origin+Light_GDD_cent+Origin:Light_GDD_cent+(1|Block), data=tsmegermdat, family="binomial")
matgermod.tsme<-glmmadmb(tsmey~Origin+MAT+Origin:MAT+(1|Block), data=tsmegermdat, family="binomial")
sdgermod.tsme<-glmmadmb(tsmey~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+(1|Block), data=tsmegermdat, family="binomial")
#interaction models
gddlightgermod.tsme<-glmmadmb(tsmey~Origin+GDD_total_cent+Light_GDD_cent+Origin:GDD_total_cent+Light_GDD_cent:Origin+Light_GDD_cent:GDD_total_cent+(1|Block),data=tsmegermdat, family="binomial")
gstlightgermod.tsme<-glmmadmb(tsmey~Origin+meanGST_cent+Light_GDD_cent+Origin:meanGST_cent+Light_GDD_cent:Origin+Light_GDD_cent:meanGST_cent+(1|Block), data=tsmegermdat, family="binomial")
matlightgermod.tsme<-glmmadmb(tsmey~Origin+MAT+Light_GDD_cent+Origin:Light_GDD_cent+Origin:MAT+Light_GDD_cent:MAT+(1|Block), data=tsmegermdat, family="binomial")
sdlightgermod.tsme<-glmmadmb(tsmey~Origin+snow_cover_duration_cent+Light_GDD_cent+Light_GDD_cent:Origin+Origin:snow_cover_duration_cent+Light_GDD_cent:snow_cover_duration+(1|Block), data=tsmegermdat, family="binomial")
nullgermod.tsme<-glmmadmb(tsmey~1+(1|Block), data=tsmegermdat, family="binomial")

AICctab(nullgermod.tsme,consgermod.tsme,gddgermod.tsme,gstgermod.tsme,lightgermod.tsme,light2germod.tsme,matgermod.tsme,sdgermod.tsme,gddlightgermod.tsme,gstlightgermod.tsme,matlightgermod.tsme,sdlightgermod.tsme)
##ABAM
#add microclimate data
abamgermdat <- join(abamgermdat, climdatsub3, by=c("Block","Canopy","Understory"), match="first")

#separate out species
abamgermdat$snow_cover_duration=as.numeric(abamgermdat$snow_dur)
abamgermdat$GDD_total=as.numeric(abamgermdat$GDD_total)
abamgermdat$meanGST=as.numeric(abamgermdat$meanGST)
abamgermdat$MAT=as.numeric(abamgermdat$MAT)
abamgermdat$Light_GDD=as.numeric(abamgermdat$Light_GDD)
abamgermdat$Light_Mean=as.numeric(abamgermdat$Light_Mean)
abamgermdat$Stand<-as.factor(abamgermdat$Stand)
abamgermdat$Block<-as.factor(abamgermdat$Block)
abamgermdat$snow_cover_duration_cent=scale(as.numeric(abamgermdat$snow_cover_duration))
abamgermdat$GDD_total_cent=scale(as.numeric(abamgermdat$GDD_total))
abamgermdat$meanGST_cent=scale(as.numeric(abamgermdat$meanGST))
abamgermdat$Light_Mean_cent=scale(as.numeric(abamgermdat$Light_Mean))
abamgermdat$Light_GDD_cent=scale(as.numeric(abamgermdat$Light_GDD))
abamgermdat$MAT=scale(as.numeric(abamgermdat$MAT))

#Now climate models
abamgermdat<-na.omit(abamgermdat)
abamgermdat$SeedsAdded<-as.factor(abamgermdat$SeedsAdded)
abamgermdat$abamy<-cbind(abamgermdat$SeedsAdded,abamgermdat$TotalFails)
consgermod.abam<-glmmadmb(abamy ~ Stand  +Origin + Canopy + Understory  +Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+Moss_cent + (1|Block),family="binomial", data=abamgermdat)
gddgermod.abam<-glmmadmb(abamy~Origin+GDD_total_cent+Origin:GDD_total_cent+Moss_cent + (1|Block),data=abamgermdat, family="binomial")
gstgermod.abam<-glmmadmb(abamy~Origin+meanGST_cent+Origin:meanGST_cent+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
lightgermod.abam<-glmmadmb(abamy~Origin+Light_Mean_cent+Origin:Light_Mean_cent+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
light2germod.abam<-glmmadmb(abamy~Origin+Light_GDD_cent+Origin:Light_GDD_cent+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
matgermod.abam<-glmmadmb(abamy~Origin+MAT+Origin:MAT+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
sdgermod.abam<-glmmadmb(abamy~Origin+snow_cover_duration_cent+Origin:snow_cover_duration_cent+Moss_cent + (1|Block), data=abamgermdat, family="binomial")

#interaction models
gddlightgermod.abam<-glmmadmb(abamy~Origin+GDD_total_cent+Light_GDD_cent+Origin:GDD_total_cent+Light_GDD_cent:Origin+Light_GDD_cent:GDD_total_cent+Moss_cent + (1|Block),data=abamgermdat, family="binomial")
gstlightgermod.abam<-glmmadmb(abamy~Origin+meanGST_cent+Light_GDD_cent+Origin:meanGST_cent+Light_GDD_cent:Origin+Light_GDD_cent:meanGST_cent+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
matlightgermod.abam<-glmmadmb(abamy~Origin+MAT+Light_GDD_cent+Origin:Light_GDD_cent+Origin:MAT+Light_GDD_cent:MAT+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
sdlightgermod.abam<-glmmadmb(abamy~Origin+snow_cover_duration_cent+Light_GDD_cent+Light_GDD_cent:Origin+Origin:snow_cover_duration_cent+Light_GDD_cent:snow_cover_duration+Moss_cent + (1|Block), data=abamgermdat, family="binomial")
nullgermod.abam<-glmmadmb(abamy~Moss_cent+(1|Block), data=abamgermdat, family="binomial")

AICctab(nullgermod.abam,consgermod.abam,gddgermod.abam,gstgermod.abam,lightgermod.abam,light2germod.abam,matgermod.abam,sdgermod.abam,gddlightgermod.abam,gstlightgermod.abam,matlightgermod.abam,sdlightgermod.abam)


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
addat$TotalGerms<-as.numeric(addat$TotalGerms)
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
#for (i in 1:5){
#  arrows(X[i],gmean10[3,i,1]-gse10[3,i,1],X[i],gmean10[3,i,1]+gse10[3,i,1],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[3,i,2]-gse10[3,i,2],X[i],gmean10[3,i,2]+gse10[3,i,2],length=.05,angle=90,code=0)}
points(gmean09[3,,]~X,pch=21,cex=1.5,bg=c(NA,NA,"black","gray","white"))
#points(gmean10[3,,1]~X,pch=22,cex=1.5,bg=c(NA,NA,"black","gray","white"))
points(gmean10[3,,2]~X,pch=22,cex=1.5,bg=c(NA,NA,"black","gray","white"))
legend(1,.6,legend=c("2009","2010"),pch=c(21,22),pt.bg=c("black"),angle=45,cex=1.1, pt.cex=1.5)

##ABAM
plot(gmean09[1,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c(NA,NA,"black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
mtext("Abies amabilis",side=3,line=.2, adj=0, font=3)
mtext("b",side=3,line=.2, adj=-.1)
for (i in 1:5){
  arrows(X[i],gmean09[1,i,]-gse09[1,i,],X[i],gmean09[1,i,]+gse09[1,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[1,i,1]-gse10[1,i,1],X[i],gmean10[1,i,1]+gse10[1,i,1],length=.05,angle=90,code=0)}
#for (i in 1:5){
#  arrows(X[i],gmean10[1,i,2]-gse10[1,i,2],X[i],gmean10[1,i,2]+gse10[1,i,2],length=.05,angle=90,code=0)}
points(gmean09[1,,]~X,pch=21,cex=1.5,bg=c("black","gray",NA,NA,"white"))
points(gmean10[1,,1]~X,pch=22,cex=1.5,bg=c("black","gray",NA,NA,"white"))
#points(gmean10[1,,2]~X,pch=23,cex=1.5,bg=c("black","gray",NA,NA,"white"))
mtext("Proportion Seeds Germinating", side=2, line=3)

##TSHE
plot(gmean09[2,,]~X,ylab="",xlab="",xlim=c(1,5),col.axis="white",ylim=c(0,0.6),type="p",bty="l",pch=21,cex=1.5,las=1,bg=c("black","gray","white"), cex.main=1.3)#ylim=c(0,.03) for inconsistent yaxes
axis(2,at=c(0,.1,.2,.3,.4,.5,.6),las=1, cex.axis=1.3)
for (i in 1:5){
  arrows(X[i],gmean09[2,i,]-gse09[2,i,],X[i],gmean09[2,i,]+gse09[2,i,],length=.05,angle=90,code=0)}
for (i in 1:5){
  arrows(X[i],gmean10[2,i,1]-gse10[2,i,1],X[i],gmean10[2,i,1]+gse10[2,i,1],length=.05,angle=90,code=0)}
#for (i in 1:5){
#  arrows(X[i],gmean10[2,i,2]-gse10[2,i,2],X[i],gmean10[2,i,2]+gse10[2,i,2],length=.05,angle=90,code=0)}
points(gmean09[2,,]~X,pch=21,cex=1.5,bg=c("gray",NA,"white",NA,NA))
points(gmean10[2,,1]~X,pch=22,cex=1.5,bg=c("gray",NA,"white",NA,NA))
#points(gmean10[2,,2]~X,pch=23,cex=1.5,bg=c("gray",NA,"white",NA,NA))
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

#Microclimate analyses
data<-read.csv("clim_gf.csv", header=TRUE)
data$Elevation_m<-NA
#Add elevations of stands:
data[which(data$stand=="NISQ"),]$Elevation_m<-668
data[which(data$stand=="TO04"),]$Elevation_m<-704
data[which(data$stand=="AV06"),]$Elevation_m<-1064
data[which(data$stand=="AM16"),]$Elevation_m<-1197
data[which(data$stand=="AE10"),]$Elevation_m<-1460
data[which(data$stand=="PARA"),]$Elevation_m<-1603
data[which(data$stand=="HIGH"),]$Elevation_m<-1650
colnames(data)[21]<-"Block"
colnames(data)[16]<-"snow_cover_duration"
data$Elevation_mfact<-as.factor(data$Elevation_m)
data$Elevation_m<-as.numeric(data$Elevation_m)
data$Block<-as.factor(data$Block)
data$Canopy<-"CompAbsent"
data[which(data$canopy=="N"),]$Canopy<-"CompPresent"
data$Understory<-"CompAbsent"
data[which(data$understory=="C"),]$Understory<-"CompPresent"
data$Year<-as.factor(data$year)

constmod.snow<-lmer(snow_cover_duration~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.snow, type="III")
#constmod.snow2<-lmer(snow_dur_cont~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
#Anova(constmod.snow2, type="III")
#summary(constmod.snow)
#constmod.gdd<-lmer(GDD_total~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
#Anova(constmod.gdd, type="III")#
#summary(constmod.gdd)
constmod.gdd2<-lmer(GDD_totaln~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.gdd2, type="III")#
#summary(constmod.gdd2)
#constmod.light<-lmer(Light_Mean~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
#Anova(constmod.light, type="III")#
constmod.light2<-lmer(Light_GDD~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block)+(1|Year),REML=FALSE,  data=data,contrasts=c(unordered="contr.sum", ordered="contr.poly"))
Anova(constmod.light2, type="III")#
summary(constmod.light)
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
pardat<-read.csv("data/PARdata.csv", header=TRUE)
justpar<-cbind((pardat[,9]),(pardat[,10]),(pardat[,11]))
parmn<-rowMeans(justpar, na.rm=TRUE)
(meanpar<-tapply(parmn,list(pardat$Canopy,pardat$Stand),mean, na.rm=TRUE))
parmn<-rowMeans(justpar)
par.lm<-lm(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.lm)
par.aov<-aov(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.aov)

#Figure of microclimate
#Effect of canopy on light
light2can<-c(fixef(constmod.light2)[8],fixef(constmod.light2)[10:15])
light2und<-c(fixef(constmod.light2)[9],fixef(constmod.light2)[16:21])
light2all<-rbind(light2und,light2can)
light2canse<-c(summary(constmod.light2)$coefficients[8,2],summary(constmod.light2)$coefficients[10:15,2])
light2undse<-c(summary(constmod.light2)$coefficients[9,2],summary(constmod.light2)$coefficients[16:21,2])
light2allse<-rbind(light2undse,light2canse)
#Effect of canopy on snow duration
snowcan<-c(fixef(constmod.snow)[8],fixef(constmod.snow)[10:15])
snowund<-c(fixef(constmod.snow)[9],fixef(constmod.snow)[16:21])
snowall<-rbind(snowund,snowcan)
snowcanse<-c(summary(constmod.snow)$coefficients[8,2],summary(constmod.snow)$coefficients[10:15,2])
snowundse<-c(summary(constmod.snow)$coefficients[9,2],summary(constmod.snow)$coefficients[16:21,2])
snowallse<-rbind(snowundse,snowcanse)
#Effect of canopy on gdd
gddcan<-c(fixef(constmod.gdd2)[8],fixef(constmod.gdd2)[10:15])
gddund<-c(fixef(constmod.gdd2)[9],fixef(constmod.gdd2)[16:21])
gddall<-rbind(gddund,gddcan)
gddcanse<-c(summary(constmod.gdd2)$coefficients[8,2],summary(constmod.gdd2)$coefficients[10:15,2])
gddundse<-c(summary(constmod.gdd2)$coefficients[9,2],summary(constmod.gdd2)$coefficients[16:21,2])
gddallse<-rbind(gddundse,gddcanse)

quartz(height=7,width=5)
par(mfcol=c(3,1),mai=c(.6,.8,.2,.1), omi=c(.7,.01,.2,.2))
####Light
plotlt<-barplot(as.matrix(rbind(light2und,light2can)),ylab="",xlab="",width=.9,names.arg=c("","","","","","",""),ylim=c(-1500,600),las=1,col=c("palegreen1","darkgreen"),xaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
x<-c(plotlt)
for (i in 1:length(light2all)){
  arrows(x[i],light2all[i]-light2allse[i],x[i],light2all[i]+light2allse[i],length=0,angle=90,code=3)}
mtext("Light",side=3,line=1, adj=0)
mtext("a",side=3,line=1, adj=-.1)
abline(h=0)
mtext("on lux", side=2, line=3, cex=.9)
legend(1,800,legend=c("Understory", "Canopy"),bty="n",fill=c("palegreen1","darkgreen"), cex=1.2)

#Snow
plotsd<-barplot(as.matrix(rbind(snowund,snowcan)),ylab="",xlab="",width=.9,names.arg=c("","","","","","",""),ylim=c(-15,10),las=1,col=c("palegreen1","darkgreen"),xaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
x<-c(plotsd)
for (i in 1:length(snowall)){
  arrows(x[i],snowall[i]-snowallse[i],x[i],snowall[i]+snowallse[i],length=0,angle=90,code=3)}
mtext("Snow",side=3,line=2, adj=0)
mtext("b",side=3,line=2, adj=-.1)
abline(h=0)
mtext("on days", side=2, line=3, cex=.9)
x<-c(plotabam)
mtext("Effect of neighbors", side=2, line=4.3)
#GDD
plotgd<-barplot(as.matrix(rbind(gddund,gddcan)),ylab="",xlab="",width=.9,names.arg=c("","","","","","",""),ylim=c(-10,15),las=1,col=c("palegreen1","darkgreen"),xaxt='n',beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3)
x<-c(plotgd)
for (i in 1:length(gddall)){
  arrows(x[i],gddall[i]-gddallse[i],x[i],gddall[i]+gddallse[i],length=0,angle=90,code=3)}
mtext("Growing degree days",side=3,line=2, adj=0)
mtext("c",side=3,line=2, adj=-.1)
abline(h=0)
mtext("on days", side=2, line=3, cex=.9)
axis(1, at = c(1.8,4.5,7.2,9.9,12.6,15.3,18),, labels = c( "668", "704","1064","1194","1460","1603","1650"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Elevation (m)",line=-13.5, adj=.6)
