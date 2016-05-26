##Analysis of 2011&2012 Mt Rainier Germination DATA

#Located in Documents/UW/Research/Mount Rainier/2011TransplantExperiment
#2011 data from tranplant experiment	
library(boot)
library(car)
library(bbmle)
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
abammod<- glm(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
abam.mmod<- glmer(abamy~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory + (1|Block), data=abamgermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig; warning: "glm.fit: fitted probabilities numerically 0 or 1 occurred"
Anova(abammod, test.statistic="LR", type="III")#
Anova(abam.mmod, test.statistic="Chisq", type="III")#
summary(abam.mmod)
anova(abammod)
round(coef(abammod),digits=3)
summary(abammod)
propdev.abamgerm=cbind(rownames(anova(abammod)),round(anova(abammod)$Dev/(anova(abammod)$"Resid. Dev"[1]-deviance(abammod)),digits=3))
propdev_fullmod.abamgerm=round((anova(abammod)$"Resid. Dev"[1]-deviance(abammod))/anova(abammod)$"Resid. Dev"[1], digits=3)

#germination in mid-origin seeds with comp absent
abammidorg=round(inv.logit(coef(abammod)[1:5]), digits=3)
#germination in upper-origin seeds with comp absent
abamuporg=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[6]+c(0,coef(abammod)[9:12])), digits=3)
#germ with canopy present, mid-origin seeds
abammidcan=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[7]+c(0,coef(abammod)[13:16])), digits=3)
#germ with canopy present,upper-origin seeds
abamupcan=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[7]+coef(abammod)[6]+c(0,coef(abammod)[13:16])+coef(abammod)[21]+c(0,coef(abammod)[9:12])), digits=3)
#germ with understory present, mid-origin seeds
abammidund=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[8]+c(0,coef(abammod)[17:20])), digits=3)
#germ with understory present,upper-origin seeds
abamupund=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[8]+coef(abammod)[6]+c(0,coef(abammod)[17:20])+coef(abammod)[22]+c(0,coef(abammod)[9:12])), digits=3)
#germ with both canopy and understory present, mid-origin seeds
abammidboth=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[7]+coef(abammod)[8]+c(0,coef(abammod)[17:20])+c(0,coef(abammod)[13:16])+coef(abammod)[23]+c(0,coef(abammod)[24:27])), digits=3)
#germ with understory present,upper-origin seeds
abamupboth=round(inv.logit(coef(abammod)[1:5]+coef(abammod)[7]+coef(abammod)[8]+c(0,coef(abammod)[17:20])+c(0,coef(abammod)[13:16])+coef(abammod)[23]+c(0,coef(abammod)[24:27])+coef(abammod)[6]+c(0,coef(abammod)[9:12])), digits=3)
#TSME
tsmemod<- glm(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmod<- glmer(tsmey~Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
summary(tsme.mmod)
Anova(tsmemod, test.statistic="LR", type="III")
Anova(tsme.mmod, test.statistic="Chisq", type="III")

propdev.tsmegerm=cbind(rownames(anova(tsmemod)),round(anova(tsmemod)$Dev/(anova(tsmemod)$"Resid. Dev"[1]-deviance(tsmemod)),digits=3))
propdev_fullmod.tsmegerm=round((anova(tsmemod)$"Resid. Dev"[1]-deviance(tsmemod))/anova(tsmemod)$"Resid. Dev"[1], digits=3)
##Can't use this model- too complex and doesnt' fit correctly- try less complex model
tsmemoda<- glm(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy, data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge
tsme.mmoda<- glmer(tsmey~Stand+Origin+Canopy+Stand:Origin+Stand:Canopy+Origin:Canopy + (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig, #this model has perfect linear separation problem so can't be used...standard errors=huge

summary(tsmemoda)
anova(tsmemoda)
propdev.tsmegerm=cbind(rownames(anova(tsmemoda)),round(anova(tsmemoda)$Dev/(anova(tsmemoda)$"Resid. Dev"[1]-deviance(tsmemoda)),digits=3))
propdev_fullmod.tsmegerm=round((anova(tsmemoda)$"Resid. Dev"[1]-deviance(tsmemoda))/anova(tsmemoda)$"Resid. Dev"[1], digits=3)

tsmemod2<- glm(tsmey~Stand*Origin, data=tsmegermdat, family=binomial)# also checked for a 3way int
tsme.mmod2<- glmer(tsmey~Stand*Origin+ (1|Block), data=tsmegermdat, family=binomial)# also checked for a 3way int

summary(tsmemod2)
Anova(tsmemod2, test.statistic="LR", type="III")#germination differed marginally by stand (p=0.08)
tsmemod3<- glm(tsmey~Origin*Canopy, data=tsmegermdat, family=binomial)# 
tsme.mmod3<- glmer(tsmey~Origin*Canopy + (1|Block), data=tsmegermdat, family=binomial)# 

summary(tsmemod3)
Anova(tsmemod3, test.statistic="LR", type="III")#differed significantly by origin
#tmse.orggerm=round(inv.logit(coef(tsmemod3.noint)), digits=4)#highest germ from lower limit origins
#round(inv.logit(confint(tsmemod3.noint)), digits=4)
tsmemod4<- glm(tsmey~Origin*Understory, data=tsmegermdat, family=binomial)# 
tsme.mmod4<- glmer(tsmey~Origin*Understory+ (1|Block), data=tsmegermdat, family=binomial)# 

summary(tsmemod4)
Anova(tsmemod4, test.statistic="LR", type="III")#Canopy had sig effect (p=0.004), no effect of understroy or interaction
tsmemod5<- glm(tsmey~Stand*Canopy, data=tsmegermdat, family=binomial)# 
tsme.mmod5<- glmer(tsmey~Stand*Canopy+ (1|Block), data=tsmegermdat, family=binomial)# 

summary(tsmemod5)#
Anova(tsmemod5, test.statistic="LR", type="III")#
library(bbmle)
AICctab(tsmemod,tsmemoda,tsmemod2,tsmemod3,tsmemod4,tsmemod5,logLik=T)
AICctab(tsme.mmod,tsme.mmoda,tsme.mmod2,tsme.mmod3,tsme.mmod4,tsme.mmod5,logLik=T)

#tsmemod3 has lowest AIC, so I'll use that model for the table

#Tshe
tshey=cbind(tshedat$TotalGerms,tshedat$TotalFails)
tshemod<- glm(tshey~-1+Stand+Origin+Canopy+Understory+Stand:Origin+Stand:Canopy+Stand:Understory+Origin:Canopy+Origin:Understory+Canopy:Understory+Stand:Canopy:Understory, data=tshedat, family=binomial)# also checked for a 3way interaction between +Origin:Canopy:Understory but not sig
summary(tshemod)##this model perfect linear separation problem so can't be used...standard errors=huge for 3way interaction
Anova(tshemod, test.statistics="LR")
propdev.tshegerm=cbind(rownames(anova(tshemod)),round(anova(tshemod)$Dev/(anova(tshemod)$"Resid. Dev"[1]-deviance(tshemod)),digits=3))
propdev.tshegerm
propdev_fullmod.tshegerm=round((anova(tshemod)$"Resid. Dev"[1]-deviance(tshemod))/anova(tshemod)$"Resid. Dev"[1], digits=3)
round(coef(tshemod),digits=3)
summary(tshemod)
#germination in low-origin seeds with comp absent
tsheloworg=round(inv.logit(coef(tshemod)[1:5]), digits=3)
#germination in mid-origin seeds with comp absent
tshemidorg=round(inv.logit(coef(tshemod)[1:5]+coef(tshemod)[6]+c(0,coef(tshemod)[10:13])), digits=3)
#germination in upper-origin seeds with comp absent
tsheuporg=round(inv.logit(coef(tshemod)[1:5]+coef(tshemod)[7]+c(0,coef(tshemod)[14:17])), digits=3)
#germ with canopy present, low-origin seeds

###perfect separation (no germinants) in many combinations of treatments, so i analyzed things separately. origin and elevation were important for all species, as was either canopy or understory, but the patterns were not consistent across species. For example, Abies had highest germinantion from upper range limit and when planted at upper and beyond upper range limits. Tsme, on the other hand, had highest germination in seeds frmo lower range limit population but the highest germination was still at the upper tange limit. canopy presence reduced germination in both species. tsuga heterophylla was not affected by canopy, but understory had a negative effect. this species had the highest germination when planted at upper range limits and from its middle origin
##want to see overall effect of origin, stand, canopy, and stand&und/stand*can interactions

tshemod2<- glm(tshey~Stand*Origin, data=tshedat, family=binomial)# 
summary(tshemod2)
Anova(tshemod2)#germination differed significantly by stand and origin, but no interaction
#tshemod2.noint<- glm(tshey~-1+Stand, data=tshedat, family=binomial)# to get mean germ for each stand
#tshe.standgerm=round(inv.logit(coef(tshemod2.noint)), digits=4)#highest germ at upper range limit
#round(inv.logit(confint(tshemod2.noint)), digits=4)#conf limits
tshemod2a<- glm(tshey~Stand*Origin+Canopy, data=tshedat, family=binomial)# 
summary(tshemod2a)
Anova(tshemod2a)

tshemod3<- glm(tshey~Canopy*Understory*Origin, data=tshedat, family=binomial)# 
summary(tshemod3)
Anova(tshemod3)#differed significantly by origin and understory, not by canopy
#tshemod3.noint<- glm(tshey~-1+Origin, data=tshedat, family=binomial)# 
#tshe.orggerm=round(inv.logit(coef(tshemod3.noint)), digits=4)#highest germ from lower limit origins
#round(inv.logit(confint(tshemod3.noint)), digits=4)

tshemod4<- glm(tshey~Canopy*Stand, data=tshedat, family=binomial)# couldn't fit model with both understory and canopy-
summary(tshemod4)
Anova(tshemod4)#stand sig, no effect of canopy, and interaction was marginal
#tshemod4.noint<- glm(tshey~-1+Canopy*Understory, data=tshedat, family=binomial)#
#tshe.compgerm=round(inv.logit(coef(tshemod4.noint)), digits=4)#highest germination withunderstory present 
#round(inv.logit(confint(tshemod4.noint)), digits=4)
tshemod5<- glm(tshey~Understory*Stand, data=tshedat, family=binomial)# couldn't fit model with both understory and canopy-
summary(tshemod5)
Anova(tshemod5)#understory had sig effect , no effect of canopy, and interaction was marginal
AICctab(tshemod,tshemod2,tshemod3,tshemod4,tshemod5,tshemod2a,logLik=T)


##Abam
abammod2<- glm(abamy~Stand*Origin, data=abamdat, family=binomial)# 
summary(abammod2)
Anova(abammod2)#germination differed significantly by stand
abammod2a<- glm(abamy~Canopy+Stand*Origin, data=abamdat, family=binomial)# 
summary(abammod2a)
Anova(abammod2a)#germi

#abammod2.noint<- glm(abamy~-1+Stand, data=abamdat, family=binomial)# to get mean germ for each stand
#abam.standgerm=round(inv.logit(coef(abammod2.noint)), digits=4)#highest germ at above upper range limit
#round(inv.logit(confint(abammod2.noint)), digits=4)#conf limits

abammod3<- glm(abamy~Canopy*Understory*Origin, data=abamdat, family=binomial)# 
summary(abammod3)
Anova(abammod3)#differed significantly by origin
#abammod3.noint<- glm(abamy~-1+Origin, data=abamdat, family=binomial)# 
#abam.orggerm=round(inv.logit(coef(abammod3.noint)), digits=4)#highest germ from upper limit origins
#round(inv.logit(confint(abammod3.noint)), digits=4)
abammod4<- glm(abamy~Canopy*Stand, data=abamdat, family=binomial)# 
summary(abammod4)
Anova(abammod4)#canopy had marginal effect (0.08), no effect of understory or interaction
#abammod4.noint<- glm(abamy~-1+Canopy*Understory, data=abamdat, family=binomial)#
#abam.compgerm=round(inv.logit(coef(abammod4.noint)), digits=4)#highest germination withunderstory present 
#round(inv.logit(confint(abammod4.noint)), digits=4)
abammod5<- glm(abamy~Understory*Stand, data=abamdat, family=binomial)# 
summary(abammod5)#this model has perfect linear separation problem so can't be used...standard errors=huge
abammod6<- glm(abamy~Comp*Stand, data=abamdat, family=binomial)# 
summary(abammod6)#this model has perfect linear separation problem so can't be used...the ses are huge as in some other models, though...but some pvalues of 1 and estimates are 1
AICctab(abammod,abammod2,abammod3,abammod4,abammod5,abammod6,abammod2a,logLik=T)


