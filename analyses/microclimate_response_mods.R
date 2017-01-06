###How do competition treatments affect microclimate?
library(lme4)
library(car)
setwd("~/Dropbox/Documents/Work/UW/Research/Mount Rainier/2014ManuscriptTransplant/microclimate")
data<-read.csv("AllStands_clim1.csv", header=TRUE)
data$Elevation_mfact<-as.factor(data$Elevation_m)
data$Elevation_m<-as.numeric(data$Elevation_m)
data$Block<-as.factor(data$Block)
head(data)
hist(data$snow_cover_duration)
plot(data$snow_cover_duration~data$Elevation_m)
plot(data$GDD_total~data$Elevation_m)
plot(data$Light_Mean~data$Elevation_m)
#Fit model for light, snow cover, and GDD with same explanatory variables used for seedlings (wihout origin)
constmod.snow<-lmer(snow_cover_duration~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
Anova(constmod.snow)#3 way interaction not significant, but understory is significant
summary(constmod.snow)
constmod.gdd<-lmer(GDD_total~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
Anova(constmod.gdd)#3 way interaction not significant, understory not significant
summary(constmod.gdd)
constmod.gdd_noint<-lmer(GDD_total~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.gdd_noint)
constmod.light<-lmer(Light_Mean~Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
Anova(constmod.light)#understory is not significant; 3 way interaction not significant
constmod.light_noint<-lmer(Light_Mean~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.light_noint)
constmod.light2<-lmer(Light_Mean~Elevation_m+Canopy+Understory+Elevation_m:Canopy+Elevation_m:Understory+Canopy:Understory+Elevation_m:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
Anova(constmod.light2)#understory is not significant; 3 way interaction not significant
##Relative abundance data for TSME, TSHE, ABAM, using file that Janneke sent
setwd("~/Dropbox/Documents/Work/UW/Research/Mount Rainier/2014ManuscriptTransplant")
badat=read.csv("MtRainierData_2.csv", header=T)
head(badat)
dim(badat)
#only want to include sites on the southwest side of mt rainier, because treeline is higher on the east side. so exclude utms with northings <5182900.89 and eastings <598104.17
badat2=badat[which(badat$UTME<599000),]
badat3=badat2[which(badat2$UTMN<5190000),]
head()
head(badat3)
max(badat$UTMN)
max(badat3$elev)
min(badat3$ABAM.ba)
max(badat3$ABAM.ba)
is.na(badat3$ABAM.ba)
min(badat3$TSME.ba)
max(badat3$TSME.ba)
min(badat3$TSHE.ba)
max(badat3$TSHE.ba)
unique(badat3$ABAM.ba)


#badat3=badat2
badat3$ABAM.dens=badat3$ABAM.ba/badat3$plotsize#basal area per m2
badat3$TSHE.dens=badat3$TSHE.ba/badat3$plotsize#basal area per m2
badat3$TSME.dens=badat3$TSME.ba/badat3$plotsize#basal area per m2


spline_tshe = smooth.spline(badat3$TSHE.dens~ badat3$elevation..m., nknots=6)
spline_tsme = smooth.spline(badat3$TSME.dens~ badat3$elevation..m., nknots=6)
spline_abam = smooth.spline(badat3$ABAM.dens~ badat3$elevation..m., nknots=6)

####Figure with focal species abundance, GDD, and light vs. elevation and by comp treatment
quartz(height=7,width=6)
par(mfcol=c(3,1), mai=c(.6,.7,.2,.5), omi=c(.7,.01,.2,.5))
#abundance of focal species
plot(badat3$elevation..m.,badat3$TSHE.dens,ylab="",xlab="Elevation (m)",xlim=c(600,1700),col.axis="white",ylim=c(0,55),type="p",bty="l",pch=21,cex=.8,cex.lab=1.3,las=1,col="white", cex.main=1.3)
axis(2,las=1, cex.axis=1.3)
axis(1,las=1, cex.axis=1.3)
mtext("Abundance",side=2,line=3.5,cex=.9)
mtext("(basal area/m2)",side=2,line=2.5,cex=.8)
mtext("a",side=3,line=1, adj=-.1)
points(badat3$elevation..m.,badat3$TSME.dens,pch=21,cex=.8,col="darkblue")
points(badat3$elevation..m.,badat3$ABAM.dens,pch=21,cex=.8,co="goldenrod")
points(badat3$elevation..m.,badat3$TSHE.dens,pch=21,cex=.8,co="darkred")
#legend(520,53,legend=c("Abies amabilis","Tsuga heterophylla", "Tsuga mertensiana"), lty=1, lwd=1.5,col=c("goldenrod", "darkred","darkblue"),cex=1.2, bty='n',text.font=3)
lines(spline_tshe,col="darkred",lwd=2)
lines(spline_tsme,col="darkblue",lwd=2)
lines(spline_abam,col="goldenrod",lwd=2)
mtext("A. amabilis",side=1,at=c(1620),line=-8,cex=.8, font=3,adj=0)
mtext("T. mertensiana",side=1,at=c(1620),line=-3,cex=.8, font=3, adj=0)
mtext("T. heterophylla",side=1,at=c(1600),line=-1,cex=.8, font=3,adj=0)
points(x=c(668,704,1064,1197,1460,1605,1676), y=c(rep(55,times=7)), pch=15, cex=1.5)

x=c(668,704,1064,1197,1460,1605,1676)
#now light
lightnocan=fixef(constmod.light_noint)[1:7]#modeled mean gdd by elevation, without canopy
lightnocan.CI= confint(constmod.gdd_noint)[3:9,]#95% confint for fixed effects
plot(lightnocan~x,ylab="",xlab="",xlim=c(500,2000),col.axis="white",ylim=c(200,2000),log="y",type="p",bty="l",pch=21,cex=1.8,cex.lab=1.3,las=1,bg="white", cex.main=1.3, cex.axis=1.3)
#axis(2,at=c(500,1000,1500,2000),las=1, cex.axis=1.3)
mtext("b",side=3,line=1, adj=-.1)
mtext("Visible light (lux)",side=2,line=3.5, cex=.9)
#for (i in 1:length(lightnocan)){
#arrows(x[i],lightnocan.CI[i,1],x[i],lightnocan.CI[i,2],length=0.03,angle=90,code=3)}
lines(lightnocan~x,lty=3)
points(lightnocan~x,pch=21,cex=1.8,bg="white")
#add points with canopy present
data$Canopy <-relevel(data$Canopy, ref = "CompPresent")
data$Understory <-relevel(data$Canopy, ref = "CompAbsent")
constmod.light_noint2<-lmer(Light_Mean~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
lightcan=fixef(constmod.light_noint2)[1:7]#modeled mean gdd by elevation, with canopy
#lightcan.CI= confint(constmod.light_noint2)[3:9,]#95% confint for fixed effects
#for (i in 1:length(lightcan)){
# arrows(x[i]+2,lightcan.CI[i,1],x[i],lightcan.CI[i,2],length=0.03,angle=90,code=3)}
lines(lightcan~c(x+5),lty=1)
points(lightcan~c(x+5),pch=21,cex=1.8,bg="darkgreen")
#add points with understory present
data$Canopy <-relevel(data$Canopy, ref = "CompAbsent")
data$Understory <-relevel(data$Understory, ref = "CompPresent")
constmod.light_noint3<-lmer(Light_Mean~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
lightund=fixef(constmod.light_noint3)[1:7]#modeled mean gdd by elevation, with canopy
#lightund.CI= confint(constmod.gdd_noint3)[3:9,]#95% confint for fixed effects
#for (i in 1:length(lightund)){
arrows(x[i],lightund.CI[i,1],x[i],lightund.CI[i,2],length=0.03,angle=90,code=3)}
lines(lightund~c(x+10),lty=2)
points(lightund~c(x+10),pch=21,cex=1.8,bg="palegreen1")
#legend(500,130,legend=c("Competition Removed", "Understory Present","Canopy Present"),bty="n",pch=21,pt.bg=c("white","palegreen1","darkgreen"), cex=1.3, pt.cex=1.3)
mtext("Competition Removed",side=1,at=c(1710),line=-9.4,cex=.8,adj=0)
mtext("Understory Present",side=1,at=c(1710),line=-8.1,cex=.8, adj=0)
mtext("Canopy Present",side=1,at=c(1710),line=-5,cex=.8,adj=0)
axis(2,las=1, cex.axis=1.3)

#growing degree days
gddnocan=fixef(constmod.gdd_noint)[1:7]#modeled mean gdd by elevation, without canopy
gddnocan.CI= confint(constmod.gdd_noint)[3:9,]#95% confint for fixed effects
plot(gddnocan~x,ylab="",xlab="Elevation (m)",xlim=c(500,2000),ylim=c(30,200),log="y",type="p",bty="l",pch=21,cex=1.8,cex.lab=1.3,las=1,bg="white", cex.main=1.3, cex.axis=1.3,)
mtext("Growing degree days",side=2,line=3.5,cex=.9)
mtext("c",side=3,line=1, adj=-.1)
#for (i in 1:length(gddnocan)){
 # arrows(x[i],gddnocan.CI[i,1],x[i],gddnocan.CI[i,2],length=0.03,angle=90,code=3)}
lines(gddnocan~x,lty=3)
points(gddnocan~x,pch=21,cex=1.8,bg="white")
#add points with canopy present
data$Canopy <-relevel(data$Canopy, ref = "CompPresent")
constmod.gdd_noint2<-lmer(GDD_total~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
gddcan=fixef(constmod.gdd_noint2)[1:7]#modeled mean gdd by elevation, with canopy
gddcan.CI= confint(constmod.gdd_noint2)[3:9,]#95% confint for fixed effects
#for (i in 1:length(gddcan)){
 # arrows(x[i]+2,gddcan.CI[i,1],x[i],gddcan.CI[i,2],length=0.03,angle=90,code=3)}
lines(gddcan~c(x+10),lty=1)
points(gddcan~c(x+10),pch=21,cex=1.8,bg="darkgreen")
#add points with understory present
data$Canopy <-relevel(data$Canopy, ref = "CompAbsent")
data$Understory <-relevel(data$Understory, ref = "CompPresent")
constmod.gdd_noint3<-lmer(GDD_total~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
gddund=fixef(constmod.gdd_noint3)[1:7]#modeled mean gdd by elevation, with canopy
gddund.CI= confint(constmod.gdd_noint3)[3:9,]#95% confint for fixed effects
#for (i in 1:length(gddund)){
 # arrows(x[i],gddund.CI[i,1],x[i],gddund.CI[i,2],length=0.03,angle=90,code=3)}
lines(gddund~c(x+20),lty=2)
points(gddund~c(x+20),pch=21,cex=1.8,bg="palegreen1")
#legend(500,130,legend=c("Competition Removed", "Understory Present","Canopy Present"),bty="n",pch=21,pt.bg=c("white","palegreen1","darkgreen"), cex=1.3, pt.cex=1.3)

##Try same figure, but with difference between no competition and understory/canopy (similar to figure for survival/growth figure)
quartz(height=7,width=6)
par(mfcol=c(3,1), mai=c(.6,.7,.2,.5), omi=c(.7,.01,.2,.5))
#abundance of focal species
plot(badat3$elevation..m.,badat3$TSHE.dens,ylab="",xlab="Elevation (m)",xlim=c(600,1700),col.axis="white",ylim=c(0,55),type="p",bty="l",pch=21,cex=.8,cex.lab=1.3,las=1,col="white", cex.main=1.3)
axis(2,las=1, cex.axis=1.3)
axis(1,las=1, cex.axis=1.3)
mtext("Abundance",side=2,line=3.5,cex=.9)
mtext("(basal area/m2)",side=2,line=2.5,cex=.8)
mtext("a",side=3,line=1, adj=-.1)
#points(badat3$elevation..m.,badat3$TSME.dens,pch=21,cex=.8,col="darkblue")
#points(badat3$elevation..m.,badat3$ABAM.dens,pch=21,cex=.8,co="goldenrod")
#points(badat3$elevation..m.,badat3$TSHE.dens,pch=21,cex=.8,co="darkred")
#legend(520,53,legend=c("Abies amabilis","Tsuga heterophylla", "Tsuga mertensiana"), lty=1, lwd=1.5,col=c("goldenrod", "darkred","darkblue"),cex=1.2, bty='n',text.font=3)
lines(spline_tshe,col="darkred",lwd=2)
lines(spline_tsme,col="darkblue",lwd=2)
lines(spline_abam,col="goldenrod",lwd=2)
mtext("A. amabilis",side=1,at=c(1620),line=-8,cex=.8, font=3,adj=0)
mtext("T. mertensiana",side=1,at=c(1620),line=-3,cex=.8, font=3, adj=0)
mtext("T. heterophylla",side=1,at=c(1600),line=-1,cex=.8, font=3,adj=0)
points(x=c(668,704,1064,1197,1460,1605,1676), y=c(rep(55,times=7)), pch=15, cex=1.5)

x=c(668,704,1064,1197,1460,1605,1676)
#now light: plot, with change in effects of comp vs no comp, as a barplot
#Fit models to get effect of elevation on light and snow duration in canopy vs gap
constmod.gdd_noint<-lmer(GDD_total~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.gdd_noint)
constmod.light_noint<-lmer(Light_Mean~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.light_noint)

#get coefs and their ses
gddcan=c(fixef(constmod.gdd_noint)[8],fixef(constmod.gdd_noint)[8]+fixef(constmod.gdd_noint)[10:15])
gddund=c(fixef(constmod.gdd_noint)[9],fixef(constmod.gdd_noint)[9]+fixef(constmod.gdd_noint)[16:21])
gddcan.se=c(summary(constmod.gdd_noint)$coef[8,2],summary(constmod.gdd_noint)$coef[10:15,2])
gddund.se=c(summary(constmod.gdd_noint)$coef[9,2],summary(constmod.gdd_noint)$coef[16:21,2])

ltcan=c(fixef(constmod.light_noint)[8],fixef(constmod.light_noint)[8]+fixef(constmod.light_noint)[10:15])
ltund=c(fixef(constmod.light_noint)[9],fixef(constmod.light_noint)[9]+fixef(constmod.light_noint)[16:21])
ltcan.se=c(summary(constmod.light_noint)$coef[8,2],summary(constmod.light_noint)$coef[10:15,2])
ltund.se=c(summary(constmod.light_noint)$coef[9,2],summary(constmod.light_noint)$coef[16:21,2])

#figure b: light
plotlt<-barplot(as.matrix(rbind(ltund,ltcan)),width=.9,ylab="",xlab="",names.arg=c("","","","","","",""),ylim=c(-1600,0),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
axis(2,las=1, cex.axis=1.1)
mtext("b",side=3,line=3, adj=-.1)
abline(h=0)
x2<-c(plotlt)
alllt.lcl=rbind(ltund-ltund.se,ltcan-ltcan.se)
alllt.ucl=rbind(ltund+ltund.se,ltcan+ltcan.se)
#add error bars
for (i in 1:length(alllt.lcl)){
  arrows(x2[i],alllt.lcl[i],x2[i],alllt.ucl[i],length=.03,angle=90,code=0)}
#legend(.3,1.1,legend=c("Understory", "Canopy"),bty="n",fill=c("palegreen1","darkgreen"), cex=1.2)
mtext("Effect of neighbors on light (lux)", side=2, line=4, cex=.9)
#mtext("(Relative to no neighbors)", side=2, line=3.1, cex=.8)
axis(1, at = c(1.8,4.5,7.2,9.9,12.5,15.3,18), labels = c( "668","704","1064", "1197","1460","1605","1676"), tick = FALSE, cex.axis=1.1, line=-.5)


#figure c: gdd
plotgdd<-barplot(as.matrix(rbind(gddund,gddcan)),width=.9,ylab="",xlab="",names.arg=c("","","","","","",""),ylim=c(-5,26),las=1,col=c("palegreen1","darkgreen"),beside=TRUE,cex.names=1.3,cex.lab=1,cex.main=1.5,cex.axis=1.3,xaxt='n', yaxt='n')
axis(2,las=1, cex.axis=1.1)
mtext("c",side=3,line=2 , adj=-.1)
abline(h=0)
x2<-c(plotgdd)
allgdd.lcl=rbind(gddund-gddund.se,gddcan-gddcan.se)
allgdd.ucl=rbind(gddund+gddund.se,gddcan+gddcan.se)
#add error bars
for (i in 1:length(allgdd.lcl)){
  arrows(x2[i],allgdd.lcl[i],x2[i],allgdd.ucl[i],length=.03,angle=90,code=0)}
#legend(.3,1.1,legend=c("Understory", "Canopy"),bty="n",fill=c("palegreen1","darkgreen"), cex=1.2)
mtext("Effect of neighbors on GDD", side=2, line=4, cex=.9)
#mtext("(Relative to no neighbors)", side=2, line=3.1, cex=.8)

axis(1, at = c(1.8,4.5,7.2,9.9,12.5,15.3,18), labels = c( "668","704","1064", "1197","1460","1605","1676"), tick = FALSE, cex.axis=1.1, line=-.5)
mtext("Elevation (m)",line=-13.5, adj=.5, cex=.9)


###Points instead of bars- i don't like this one.
##Try same figure, but with difference between no competition and understory/canopy (similar to figure for survival/growth figure)
quartz(height=7,width=6)
par(mfcol=c(3,1), mai=c(.6,.7,.2,.5), omi=c(.7,.01,.2,.5))
#abundance of focal species
plot(badat3$elevation..m.,badat3$TSHE.dens,ylab="",xlab="",xlim=c(500,2000),col.axis="white",ylim=c(0,40),type="p",bty="l",pch=21,cex=.8,cex.lab=1.3,las=1,col="white", cex.main=1.3)
axis(2,las=1, cex.axis=1.3)
mtext("Abundance",side=2,line=3.5,cex=.9)
mtext("(basal area/m2)",side=2,line=2.5,cex=.8)
mtext("a",side=3,line=1, adj=-.1)
#points(badat$elevation..m.,badat$TSME.dens,pch=21,cex=.8,col="darkblue")
#points(badat$elevation..m.,badat$ABAM.dens,pch=21,cex=.8,co="goldenrod")
#points(badat$elevation..m.,badat$TSHE.dens,pch=21,cex=.8,co="darkred")
#legend(520,53,legend=c("Abies amabilis","Tsuga heterophylla", "Tsuga mertensiana"), lty=1, lwd=1.5,col=c("goldenrod", "darkred","darkblue"),cex=1.2, bty='n',text.font=3)
lines(spline_tshe,col="darkred",lwd=2)
lines(spline_tsme,col="darkblue",lwd=2)
lines(spline_abam,col="goldenrod",lwd=2)
mtext("A. amabilis",side=1,at=c(1700),line=-8,cex=.8, font=3,adj=0)
mtext("T. mertensiana",side=1,at=c(1750),line=-4.4,cex=.8, font=3, adj=0)
mtext("T. heterophylla",side=1,at=c(1800),line=-1,cex=.8, font=3,adj=0)
#points(x=c(668,704,1064,1197,1460,1605,1676), y=c(rep(40,times=7)), pch=15, cex=1.5)

x=c(668,704,1064,1197,1460,1605,1676)
#now light
#plot, with change in effects of comp vs no comp
#Fit models to get effect of elevation on light and snow duration in canopy vs gap
constmod.gdd_noint<-lmer(GDD_total~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.gdd_noint)
constmod.light_noint<-lmer(Light_Mean~-1+Elevation_mfact+Canopy+Understory+Elevation_mfact:Canopy+Elevation_mfact:Understory+Canopy:Understory+Elevation_mfact:Canopy:Understory+(1|Block),REML=FALSE,  data=data)
summary(constmod.light_noint)

#get coefs and their ses
gddcan=c(fixef(constmod.gdd_noint)[8],fixef(constmod.gdd_noint)[8]+fixef(constmod.gdd_noint)[10:15])
gddund=c(fixef(constmod.gdd_noint)[9],fixef(constmod.gdd_noint)[9]+fixef(constmod.gdd_noint)[16:21])
gddcan.se=c(summary(constmod.gdd_noint)$coef[8,2],summary(constmod.gdd_noint)$coef[10:15,2])
gddund.se=c(summary(constmod.gdd_noint)$coef[9,2],summary(constmod.gdd_noint)$coef[16:21,2])

ltcan=c(fixef(constmod.light_noint)[8],fixef(constmod.light_noint)[8]+fixef(constmod.light_noint)[10:15])
ltund=c(fixef(constmod.light_noint)[9],fixef(constmod.light_noint)[9]+fixef(constmod.light_noint)[16:21])
ltcan.se=c(summary(constmod.light_noint)$coef[8,2],summary(constmod.light_noint)$coef[10:15,2])
ltund.se=c(summary(constmod.light_noint)$coef[9,2],summary(constmod.light_noint)$coef[16:21,2])


#figure b: light
plot(ltund~x,ylab="",xlab="",xlim=c(500,2000),col.axis="white",ylim=c(-1500,0),type="p",bty="l",pch=21,cex=1.8,las=1,bg="palegreen1", cex.main=1.3)
axis(2,las=1, cex.axis=1.3)
mtext("b",side=3,line=1, adj=-.1)
mtext("Effect of neighbors on visible light", side=2, line=3.5, cex=.9)
mtext("(Difference in lux, relative to no neighbors)", side=2, line=2.5, cex=.8)
abline(h=0)
for (i in 1:length(ltund)){
  arrows(x[i],ltund[i]+ltund.se[i],x[i],ltund[i]-ltund.se[i],length=0.05,angle=90,code=0)}
for (i in 1:length(ltcan)){
  arrows(x[i]+20,ltcan[i]+ltcan.se[i],x[i]+20,ltcan[i]-ltcan.se[i],length=0.05,angle=90,code=0)}
points(ltund~x,pch=21,cex=1.8,bg="palegreen1")
#add points with canopy present
points(ltcan~c(x+20),pch=21,cex=1.8,bg="darkgreen")
#legend(680,.7,legend=c("Understory", "Canopy"),bty="n",pch=21,pt.bg=c("palegreen1","darkgreen"), cex=1.3, pt.cex=1.3)

#figure c: gdd
plot(gddund~x,ylab="",xlab="",xlim=c(500,2000),col.axis="white",ylim=c(-5,25),type="p",bty="l",pch=21,cex=1.8,las=1,bg="palegreen1", cex.main=1.3)
axis(2,las=1, cex.axis=1.3)
mtext("c",side=3,line=1, adj=-.1)
mtext("Effect of neighbors on growing degree days", side=2, line=3.5, cex=.9)
mtext("(Difference, relative to no neighbors)", side=2, line=2.5, cex=.8)
abline(h=0)
for (i in 1:length(gddund)){
  arrows(x[i],gddund[i]+gddund.se[i],x[i],gddund[i]-gddund.se[i],length=0.05,angle=90,code=0)}
for (i in 1:length(gddcan)){
  arrows(x[i]+20,gddcan[i]+gddcan.se[i],x[i]+20,gddcan[i]-gddcan.se[i],length=0.05,angle=90,code=0)}
points(gddund~x,pch=21,cex=1.8,bg="palegreen1")
#add points with canopy present
points(gddcan~c(x+20),pch=21,cex=1.8,bg="darkgreen")
#legend(680,.7,legend=c("Understory", "Canopy"),bty="n",pch=21,pt.bg=c("palegreen1","darkgreen"), cex=1.3, pt.cex=1.3)










#Fit models to get effect of elevation on light and snow duration in canopy vs gap
mod2light<-lmer(Light_Mean~-1+Elevation_m+Canopy+(1|Block), data=data)#compared it to model with interaction, but this one had better fit
summary(mod2light)
mod1snow<-lmer(snow_cover_duration~-1+Elevation_m*Canopy+(1|Block), data=data,REML=FALSE)#compared it to model with interaction, but this one had better fit
summary(mod1snow)
mod2snow<-lmer(snow_cover_duration~-1+Elevation_m+Canopy+(1|Block), data=data, REML=FALSE)#compared it to model with interaction, but this one had better fit
summary(mod2snow)
mod1snow_fact<-lmer(snow_cover_duration~-1+Elevation_mfact*Canopy+(1|Block), data=data,REML=FALSE)#compared it to model with interaction, but this one had better fit
summary(mod1snow_fact)
mod2snow_fact<-lmer(snow_cover_duration~-1+Elevation_mfact+Canopy+(1|Block), data=data, REML=FALSE)#compared it to model with interaction, but this one had better fit
summary(mod2snow_fact)
mod3snow_fact<-lmer(snow_cover_duration~-1+Elevation_mfact+Canopy+Understory+(1|Block), data=data, REML=FALSE)#compared it to model with interaction, but this one had better fit
summary(mod3snow_fact)
anova(mod1snow,mod2snow,mod1snow_fact,mod2snow_fact,mod3snow_fact, test="Chi")
#lowest AIC= with factorial elevation, interaction between elevation and canopy (mod1snow_fact)

#Make graph showing light and snow duration in gaps versus under canopy, across elevation
quartz(width=8, height=6)
par(mai=c(1, 1, .5, .5), omi=c(.1,.1,.01,.01), mfrow=c(2,1)) #margin space (bottom, left, top, right)data=cit
plot(data$Elevation_m,data$snow_cover_duration,pch=21, bg=c("white","gray")[data$Canopy],xlab="Elevation (m)", ylab="Snow duration (days/year)", bty="l")
x<-seq(min(data$Elevation_m),max(data$Elevation_m), 1)
lines(x,fixef(mod2snow)[2]+fixef(mod2snow)[1]*x,lty=3,lwd=2)#slope and intercept from mod2snow model for gap, across all sites)
lines(x,fixef(mod2snow)[3]+fixef(mod2snow)[1]*x,lty=1,lwd=2)#slope and intercept from mod2snow model for nongap, across all sites)

plot(data$Elevation_m,data$Light_Mean,pch=21, bg=c("white","gray")[data$Canopy],xlab="Elevation (m)", ylab="Light (lux)", bty="l")
lines(x,fixef(mod2light)[2]+fixef(mod2light)[1]*x,lty=3,lwd=2)#slope and intercept from mod2light model for gap, across all sites)
lines(x,fixef(mod2light)[3]+fixef(mod2light)[1]*x,lty=1,lwd=2)#nongap
unique(data$Canopy)
interaction.plot(data$Elevation_m,data$Canopy,data$snow_cover_duration)
interaction.plot(data$Elevation_m,data$Understory,data$snow_cover_duration)
mod1<-glmer(snow_cover_duration~Elevation_m*Canopy*Understory+(1|Block), family=poisson, data=data)
summary(mod1)#AIC=1041 understory & all 3way interactions nonsignificant so remove, 
mod2<-glmer(snow_cover_duration~-1+Elevation_m+Canopy+(1|Block), family=poisson, data=data)
summary(mod2)#AIC=1018, Canopy presence has negative effect on snow duration (fewer days with snow when canopy is present), understory treatment does not affect microclimate
mod3<-glmer(snow_cover_duration~-1+Elevation_m*Canopy+(1|Block), family=poisson, data=data)
summary(mod3)#AIC=1020, Interaction does not provide much better fit so seems that it is not significant, so mod2 seems best...

boxplot(data$snow_cover_duration~data$Understory)#no difference
boxplot(data$snow_cover_duration~data$Canopy)
boxplot(data$snow_cover_duration~data$Elevation_m)
boxplot(data$GDD_total~data$Understory)#no difference
boxplot(data$GDD_total~data$Canopy)
boxplot(data$GDD_total~data$Elevation_m)
mod1gdd<-glmer(GDD_total~-1+Elevation_m*Canopy*Understory+ (1|Block), family=poisson, data=data)
summary(mod1gdd)#AIC=1343 understory & most 3way interactions nonsignificant so remove, 
mod2gdd<-glmer(GDD_total~-1+Elevation_m+Canopy+Understory+Elevation_m:Canopy+Canopy:Understory+Elevation_m:Understory+ (1|Block), family=poisson, data=data)
summary(mod2gdd)#AIC=1356, worse fit, & LRT sig different so interactions are important!
anova(mod1gdd,mod2gdd, test="Chi")
mod3gdd<-glmer(GDD_total~-1+Elevation_m*Canopy+ (1|Block), family=poisson, data=data)
summary(mod3gdd)#AIC=1387, Interaction provides better fit. seems to be best-fit model
anova(mod1gdd,mod3gdd, test="Chi")
mod4gdd<-glmer(GDD_total~-1+Elevation_m+Canopy+Understory+ (1|Block), family=poisson, data=data)
summary(mod4gdd)#AIC=1446, worse fit, tried also with understory, present- worse fit
anova(mod1gdd,mod4gdd, test="Chi")
nongapdat<-data[data$Canopy=="Present",]
gapdat<-data[data$Canopy=="Absent",]
boxplot(nongapdat$snow_cover_duration~nongapdat$Elevation_m, col="darkgreen")
par(new=T)
boxplot(gapdat$snow_cover_duration~gapdat$Elevation_m, col="white")
interaction.plot(data$Elevation_m,data$Canopy,data$GDD_total)

mod1light<-lmer(Light_Mean~-1+Elevation_m*Canopy*Understory+(1|Block), data=data)
summary(mod1light)#AIC=1343 understory & most 3way interactions nonsignificant so remove, 
mod2light<-lmer(Light_Mean~-1+Elevation_m+Canopy+Understory+Elevation_m:Canopy+Canopy:Understory+Elevation_m:Understory+(1|Block), data=data)
summary(mod2light)
anova(mod1light,mod2light,test="Chi")#no difference so go with mod2
mod3light<-lmer(Light_Mean~Elevation_m+Canopy+Understory+Elevation_m:Canopy+Canopy:Understory+(1|Block), data=data)
anova(mod2light,mod3light,test="Chi")#no difference so go with mod3
summary(mod3light)
mod4light<-lmer(Light_Mean~-1+Elevation_m+Canopy+Understory+Canopy:Understory+(1|Block), data=data)
anova(mod4light,mod3light,test="Chi")#there is a difference so go with mod3
mod5light<-lmer(Light_Mean~-1+Elevation_m+Canopy+Understory+Elevation_m:Canopy+(1|Block), data=data)
anova(mod5light,mod3light,test="Chi")#there is a difference, so stick with mod3
#try removing elevation
mod6light<-lmer(Light_Mean~Canopy+Understory+(1|Block), data=data)
anova(mod6light,mod3light,test="Chi")#mod 3 is better fit

meansnodurnongaps<-tapply(as.numeric(nongapdat$snow_cover_duration),list(nongapdat$Elevation_m),mean)
meansnodurgaps<-tapply(as.numeric(gapdat$snow_cover_duration),list(gapdat$Elevation_m),mean)
sdsnodurnongaps<-tapply(as.numeric(nongapdat$snow_cover_duration),list(nongapdat$Elevation_m),sd)
sdsnodurgaps<-tapply(as.numeric(gapdat$snow_cover_duration),list(gapdat$Elevation_m),sd)
sesnodurnongaps<-sdsnodurnongaps/sqrt(10)
sesnodurgaps<-sdsnodurgaps/sqrt(10)

quartz(width=8, height=6)
par(mai=c(1, 1, .5, .5), omi=c(.1,.1,.01,.01)) #margin space (bottom, left, top, right)data=cit
plot(names(meansnodurnongaps),meansnodurnongaps, pch=21,bg="darkgreen", xlab="Elevation (m)", ylab="Snow duration (days)", cex=1.4, cex.lab=1.2, cex.axis=1.2, ylim=c(0,270))
lines(names(meansnodurgaps),meansnodurgaps,lty=1)
lines(names(meansnodurnongaps),meansnodurnongaps,col="darkgreen", lty=1)
points(names(meansnodurgaps),meansnodurgaps, pch=21,bg="white", cex=1.2)
legend(1450,50,legend=c("Absent", "Present"), pch=c(21, 21), pt.bg=c("white", "darkgreen"), title="Canopy",cex=1.2, bty='n')

#ERror bars- so narrow, not worth including
#x<-as.numeric(names(sesnodurgaps))
#for (i in 1:length(sesnodurgaps)){
#	arrows(x[i],meansnodurgaps[i]-sesnodurgaps[i],x[i],meansnodurgaps[i]+sesnodurgaps[i],length=0.1,angle=90,code=3)}

#for (i in 1:length(sesnodurnongaps)){
#	arrows(x[i],meansnodurnongaps[i]-sesnodurnongaps[i],x[i],meansnodurnongaps[i]+sesnodurnongaps[i],length=0.1,angle=90,code=3)}

nongapnound<-nongapdat[nongapdat$Understory=="Absent",]
nongapund<-nongapdat[nongapdat$Understory=="Present",]
gapnound<-gapdat[gapdat$Understory=="Absent",]
gapund<-gapdat[gapdat$Understory=="Present",]
meansnodurnongaps<-tapply(as.numeric(nongapdat$snow_cover_duration),list(nongapdat$Elevation_m),mean)
meansnodurgaps<-tapply(as.numeric(gapdat$snow_cover_duration),list(gapdat$Elevation_m),mean)
sdsnodurnongaps<-tapply(as.numeric(nongapdat$snow_cover_duration),list(nongapdat$Elevation_m),sd)
sdsnodurgaps<-tapply(as.numeric(gapdat$snow_cover_duration),list(gapdat$Elevation_m),sd)
sesnodurnongaps<-sdsnodurnongaps/sqrt(10)
sesnodurgaps<-sdsnodurgaps/sqrt(10)



meangddnongaps<-tapply(as.numeric(nongapdat$GDD_total),list(nongapdat$Elevation_m),mean)
meangddgaps<-tapply(as.numeric(gapdat$GDD_total),list(gapdat$Elevation_m),mean)
sdgddnongaps<-tapply(as.numeric(nongapdat$GDD_total),list(nongapdat$Elevation_m),sd)
sdgddgaps<-tapply(as.numeric(gapdat$GDD_total),list(gapdat$Elevation_m),sd)
segddnongaps<-sdgddnongaps/sqrt(10)
segddgaps<-sdgddgaps/sqrt(10)

quartz(width=8, height=6)
par(mai=c(1, 1, .5, .5), omi=c(.1,.1,.01,.01)) #margin space (bottom, left, top, right)data=cit
plot(names(meangddnongaps),meangddnongaps, pch=21,bg="darkgreen", xlab="Elevation (m)", ylab="# Growing degree days", cex=1.4, cex.lab=1.2, cex.axis=1.2, ylim=c(0,200))
lines(names(meangddgaps),meangddgaps,lty=1)
lines(names(meangddnongaps),meangddnongaps,col="darkgreen", lty=1)
points(names(meangddgaps),meangddgaps, pch=21,bg="white", cex=1.2)
legend(1450,200,legend=c("Absent", "Present"), pch=c(21, 21), pt.bg=c("white", "darkgreen"), title="Canopy",cex=1.2, bty='n')

####Plot with difference between gaps and nongpas for snow duration and gdd. 
##EFfect of canopy competition on Snow duration, GDD, and light
#gdd:
gddmean<-tapply(as.numeric(data$GDD_total),list(data$Elevation_m),mean)
gddnongaps<-tapply(as.numeric(nongapdat$GDD_total),list(nongapdat$Block),mean)
gddgaps<-tapply(as.numeric(gapdat$GDD_total),list(gapdat$Block),mean)
gdddiff<-gddnongaps-gddgaps
sdnongaps<-tapply(as.numeric(nongapdat$snow_cover_duration),list(nongapdat$Block),mean)
sdgaps<-tapply(as.numeric(gapdat$snow_cover_duration),list(gapdat$Block),mean)
sddiff<-sdnongaps-sdgaps
blockelevs<-c(rep(1460,times=5),rep(1197,times=5),rep(1064,times=5),rep(1650,times=5),rep(668,times=5),rep(1603,times=5),rep(704 ,times=5) )  
#visible light:
lightmean<-tapply(as.numeric(data$Light_Mean),list(data$Elevation_m),mean)
lightnongaps<-tapply(as.numeric(nongapdat$Light_Mean),list(nongapdat$Block),mean)
lightgaps<-tapply(as.numeric(gapdat$Light_Mean),list(gapdat$Block),mean)
lightgaps<-lightgaps[-13]#because it =NA
lightnongaps<-lightnongaps[-22]#because it =NA
lightgaps<-lightgaps[-21]#because it =NA
lightnongaps<-lightnongaps[-13]#because it =NA for gaps

lightdiff<-lightnongaps-lightgaps

blockclim<-as.data.frame(cbind(blockelevs,gdddiff,sddiff,lightdiff))
blockclim$blockelevs<-as.factor(blockclim$blockelevs)
meangdddiff<-tapply(as.numeric(blockclim$gdddiff),list(as.factor(blockclim$blockelevs)),mean)
meansddiff<-tapply(as.numeric(blockclim$sddiff),list(blockclim$blockelevs),mean)
meanlightdiff<-tapply(as.numeric(blockclim$lightdiff),list(blockclim$blockelevs),mean, na.rm = TRUE)
sdgdd<-tapply(as.numeric(blockclim$gdddiff),list(blockclim$blockelevs),sd)
sdsd<-tapply(as.numeric(blockclim$sddiff),list(blockclim$blockelevs),sd)
sdlight<-tapply(as.numeric(blockclim$lightdiff),list(blockclim$blockelevs),sd,na.rm = TRUE)
segdd<-sdgdd/sqrt(5)
sesd<-sdsd/sqrt(5)
selight<-sdlight/sqrt(5)
sdmean<-tapply(as.numeric(data$snow_cover_duration),list(data$Elevation_m),mean)

gddperc<-meangdddiff/gddmean*100#difference in gdd as a percentage of average growing degree days across all treatments.
sdperc<-meansddiff/sdmean*100#difference in sd as a percentage of average growing degree days across all treatments.
lightperc<-meanlightdiff/lightmean*100#difference in light as a percentage of average growing degree days across all treatments.

quartz(height=7,width=6)
par(mfcol=c(3,1),mai=c(.8,.8,.5,.1), omi=c(.7,.01,.2,.2))
plotsd<-barplot(meansddiff,ylab="Effect of canopy on snow",xlab="",ylim=c(-28,4),las=1,col="green3",beside=TRUE,names.arg=c(names(meansddiff)),cex.names=1.2,cex.lab=1.2,cex.axis=1.2)
abline(h=0)
mtext("A.",side=3,line=3.5, adj=-.1)
#error bars
x<-c(plotsd)
for (i in 1:length(meansddiff)){
	arrows(x[i],meansddiff[i]-sesd[i],x[i],meansddiff[i]+sesd[i],length=0.05,angle=90,code=3)}

plotgdd<-barplot(meangdddiff,ylab="Effect of canopy on GDD",xlab="",ylim=c(-5,60),las=1,col="green3",beside=TRUE,names.arg=c(names(meansddiff)),cex.names=1.2,cex.lab=1.2,cex.axis=1.2)
abline(h=0)
mtext("B.",side=3,line=3.5, adj=-.1)
#error bars
x<-c(plotgdd)
for (i in 1:length(meangdddiff)){
	arrows(x[i],meangdddiff[i]-segdd[i],x[i],meangdddiff[i]+segdd[i],length=0.05,angle=90,code=3)}

plotlight<-barplot(meanlightdiff,ylab="Effect of canopy on light",xlab="Elevation (m)",ylim=c(-1200,50),las=1,col="green3",beside=TRUE,names.arg=c(names(meanlightdiff)),cex.names=1.3,cex.lab=1.2,cex.main=1.5,cex.axis=1.2)
abline(h=0)
mtext("C.",side=3,line=3.5, adj=-.1)

#error bars
x<-c(plotlight)
for (i in 1:length(meanlightdiff)){
	arrows(x[i],meanlightdiff[i]-selight[i],x[i],meanlightdiff[i]+selight[i],length=0.05,angle=90,code=3)}




#########Soil analyses
soilsdat<-read.csv("SoilsData.csv", header=TRUE)
soilsdat$Elev_m<-as.factor(soilsdat$Elev_m)
head(soilsdat)
boxplot(soilsdat$SoilMoist_pce~soilsdat$Elev_m)
boxplot(soilsdat$SoilMoist_pce~soilsdat$Canopy)
interaction.plot(soilsdat$Elev_m,soilsdat$Canopy,soilsdat$SoilMoist)
interaction.plot(soilsdat$Elev_m,soilsdat$Canopy,soilsdat$C_pce)
interaction.plot(soilsdat$Elev_m,soilsdat$Canopy,soilsdat$H_pce)
interaction.plot(soilsdat$Elev_m,soilsdat$Canopy,soilsdat$N_pce)
soilmoist.aov<-aov(soilsdat$SoilMoist_pce~soilsdat$Elev_m*soilsdat$Canopy)
summary(soilmoist.aov)
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
interaction.plot(soilsdat$Elev_m,soilsdat$Canopy,soilsdat$CanopyCover)
tapply(soilsdat$CanopyCover,list(soilsdat$Canopy,soilsdat$Elev_m),mean)
CanCov.lm<-lm(soilsdat$CanopyCover~-1+soilsdat$Elev_m*soilsdat$Canopy)
summary(CanCov.lm)


###summary of GDD by elevation & gap/nongap (means)
data<-read.csv("AllStands_clim1.csv", header=TRUE)

head(data)
(meangdd<-tapply(data$GDD_total,list(data$Canopy,data$Elevation_m),mean))

data<-data[-48,]
data<-data[-88,]

(meanlight<-tapply(data$Light_Mean,list(data$Canopy,data$Elevation_m),mean))

#######summary of slope and aspect by elevation & gap/nongap (means)
setwd("~/Dropbox/Documents/Work/UW/Research/Mount Rainier/2014ManuscriptTransplant/microclimate")
spatdat<-read.csv("GapSizeData.csv", header=TRUE)
head(spatdat)
#spatdat<-spatdat[-1,]
#spatdat<-spatdat[-5,]
#spatdat<-spatdat[-67,]
#spatdat<-spatdat[-67,]
(meanslope<-tapply(spatdat$Slope...,list(spatdat$Canopy,spatdat$Stand),mean))
slope.lm<-lm(spatdat$Slope...~-1+spatdat$Stand*spatdat$Canopy)
summary(slope.lm)
slope.aov<-aov(spatdat$Slope...~-1+spatdat$Stand*spatdat$Canopy)
summary(slope.aov)
#spatdat<-spatdat[-43,]
#spatdat<-spatdat[-60,]

(meanaspect<-tapply(spatdat$Aspect..degrees,list(spatdat$Canopy,spatdat$Stand),mean, na.rm=T))
aspect.lm<-lm(spatdat$Aspect..degrees~-1+spatdat$Stand*spatdat$Canopy)
summary(aspect.lm)
aspect.aov<-aov(spatdat$Aspect..degrees~-1+spatdat$Stand*spatdat$Canopy)
summary(aspect.aov)


#PAr Data
pardat<-read.csv("PARdata.csv", header=TRUE)
justpar<-cbind((pardat[,9]),(pardat[,10]),(pardat[,11]))
is.numeric(justpar)
parmn<-rowMeans(justpar, na.rm=TRUE)
(meanpar<-tapply(parmn,list(pardat$Canopy,pardat$Stand),mean, na.rm=TRUE))
parmn<-rowMeans(justpar)
par.lm<-lm(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.lm)
par.aov<-aov(parmn~-1+pardat$Stand*pardat$Canopy)
summary(par.aov)