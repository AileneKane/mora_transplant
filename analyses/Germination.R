##Analysis of 2011&2012 Mt Rainier Germination DATA

#Located in Documents/UW/Research/Mount Rainier/2011TransplantExperiment
#2011 data from tranplant experiment	
library(boot)
library(car)
library(bbmle)
setwd("~/GitHub/mora_transplant")
dat<-read.csv("data/MORAGermData20112012.csv", header=TRUE)
#Include only plots where seeds were added
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
tsmegermdat<-tsmegermdat[1:min(which(is.na(tsmegermdat$Stand)))-1,]
tshegermdat<-tshegermdat[1:min(which(is.na(tshegermdat$Stand)))-1,]
abamgermdat<-abamgermdat[1:min(which(is.na(abamgermdat$Stand)))-1,]

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
#ABAM
#abammod<- glm(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod<- glmer(abamy~-1+Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory + (1|Block), data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
Anova(abam.mmod, test.statistic="Chisq", type="III")#model failed to converge
##Model selection, by hand
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
AICctab(abammod,abammod2,abammod3,abammod4,abammod5,abammod6,abammod2a,abammod4a,logLik=T)
AICctab(abam.mmod,abam.mmod2,abam.mmod3,abam.mmod4,abam.mmod5,abam.mmod6,abam.mmod2a,abam.mmod4a,logLik=T)
Anova(abam.mmod2a, type="III")
Anova(abammod2a, type="III")

#TSME
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
AICctab(tsmemod,tsmemoda,tsmemod2,tsmemod3,tsmemod4,tsmemod5,logLik=T)
AICctab(tsme.mmod,tsme.mmoda,tsme.mmod2,tsme.mmod3,tsme.mmod4,tsme.mmod5,logLik=T)
anova(tsme.mmod3,type="III")
#tsmemod3 has lowest AIC, so I'll use that model for the table
Anova(tsmemod3,test="LR", type="II")

#Tshe
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

AICctab(tshemod,tshemod2,tshemod3,tshemod4,tshemod5,tshemod2a,logLik=T)
AICctab(tshe.mmod,tshe.mmod2,tshe.mmod3,tshe.mmod4,tshe.mmod5,tshe.mmod2a,logLik=T)
#mod2 is best fit for tshe
Anova(tshe.mmod2, test="Chisq",type="III")
Anova(tshemod2,test="LR", type="III")
