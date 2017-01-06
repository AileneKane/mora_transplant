#Modeling climate data
#I realized that we are missing some microclimate data.
#The blocks with missing data are:
#HIGH-GX-1
#HIGH-GX-5
#HIGH-NC-1
#NISQ-NC-4
#Light for NISQ-NC-2_2012Summer.csv
#Light for AV06-GX-3_2012Summer.csv
#So, I need to use existing data to estimate missing data by fitting linear models
#Ian said to "use alinear model with:
##categorical year effect
##site effect,
##plot effect (possibly with a year:site interaction if the data justifies this). 
##We then used that model to predict missing values. 
##As long as you’ve got at least some measurements for 
##each site and year, this works like a charm. 
##Cross-validated prediction errors were typically on the order of 2-5 days.
setwd("~/Dropbox/Documents/Work/uw/Research/Mount Rainier/2014ManuscriptTransplant/microclimate")
library(lme4)
library(bbmle)
library(dplyr)
library(plyr)
clim2<-read.csv("ALLmicroclim20102013.csv", header=TRUE)
snow<-read.csv("snowmicroclim20102013.csv", header=TRUE)
clim<-join(clim2, snow, by=c("filesname","year","yrstartdate","yrenddate","cumna.temp","cumnumrows"), match="first")
head(clim)
clim$temppropna<-clim$cumna.temp/clim$cumnumrows
clim$lightpropna<-clim$cumna.light/clim$cumnumrows
#If more than 25% missing, should replace all derived values with NA
clim[which(clim$temppropna>0.25),]#after looking at these, i decided to replace all temp AND LIGHT values with NA for these, because light values also have high proportion of missing data or they look off.
clim[which(clim$temppropna>0.25),]$GDD_total<-NA
clim[which(clim$temppropna>0.25),]$GDD_totaln<-NA
clim[which(clim$temppropna>0.25),]$Light_Mean<-NA
clim[which(clim$temppropna>0.25),]$Light_GDD<-NA
clim[which(clim$temppropna>0.25),]$GST_Mean<-NA
clim[which(clim$temppropna>0.25),]$MAT<-NA
clim[which(clim$temppropna>0.25),]$snow_dur<-NA
clim[which(clim$temppropna>0.25),]$snow_dur_cont<-NA

#Add column for elevation, canopy status, understory status, and block
clim$stand<-substr(clim$filesname,1,4)
clim$block<-paste(substr(clim$filesname,1,4),substr(clim$filesname,9,9),sep="")
clim$canopy<-substr(clim$filesname,6,6)
clim$understory<-substr(clim$filesname,7,7)
hist(clim$GDD_total)#look at extreme values of GDD_total
hist(clim$Light_Mean)#look at extreme values of light
clim[which(clim$GDD_total>2000),]$GDD_total<-NA#replace extreme values with NA (these are likely a measurement error)
clim[which(clim$Light_Mean>3500),]$Light_Mean<-NA#replace extreme values with NA (these are likely a measurement error)
clim[which(clim$Light_GDD>3500),]$Light_GDD<-NA#replace extreme values with NA (these are likely a measurement error)
clim[which(clim$GST_Mean>20),]$GST_Mean<-NA#replace extreme values with NA (these are likely a measurement error)
clim[which(clim$MAT>25),]$MAT<-NA#replace extreme values with NA (these are likely a measurement error)
#save this climate data, prior to adding in predicted values
write.csv(clim,"clim.csv", row.names=FALSE)

#now fit a model for each climate variable of interest that is missing
clim$stand<-as.factor(clim$stand)
clim$block<-as.factor(clim$block)
clim$year<-as.factor(clim$year)
clim$canopy<-as.factor(clim$canopy)
clim$understory<-as.factor(clim$understory)
#GDD 
climsub<-subset(clim,select=c("stand","block","year","canopy","understory","GDD_total","GDD_totaln","Light_Mean","Light_GDD","GST_Mean","MAT"))
clim.nona<-na.omit(climsub)

gddmod<-lm(GDD_total~stand+block+year+canopy+understory,data=clim.nona)
gddmod2<-lm(GDD_total~stand+block+year+canopy*understory,data=clim.nona)
gddmod3<-lm(GDD_total~stand*year+block+canopy*understory,data=clim.nona)
gddmod4<-lm(GDD_total~stand*year+block*canopy*understory,data=clim.nona)
AICtab(gddmod,gddmod2,gddmod3,gddmod4)#gddmod4
#now try to predict missing values
predict.gdd<-predict(gddmod3,newdata=clim,type="response")
cbind(clim$GDD_total,predict.gdd)
plot(clim$GDD_total,predict.gdd)
abline(a=0,b=1,lty=1)
#GDDn
gddnmod<-lm(GDD_totaln~stand+block+year+canopy+understory,data=clim.nona)
gddnmod2<-lm(GDD_totaln~stand+block+year+canopy*understory,data=clim.nona)
gddnmod3<-lm(GDD_totaln~stand*year+block+canopy*understory,data=clim.nona)
gddnmod4<-lm(GDD_totaln~stand*year+block*canopy*understory,data=clim.nona)
AICtab(gddnmod,gddnmod2,gddnmod3,gddmod4)#gddnmod3
#now try to predict missing values
predict.gddn<-predict(gddnmod3,newdata=clim,type="response")
cbind(clim$GDD_totaln,predict.gddn)
plot(clim$GDD_totaln,predict.gddn)
abline(a=0,b=1,lty=1)
#Light
ltmod<-lm(Light_Mean~stand+block+year+canopy+understory,data=clim.nona)
ltmod2<-lm(Light_Mean~stand+block+year+canopy*understory,data=clim.nona)
ltmod3<-lm(Light_Mean~stand*year+block+canopy*understory,data=clim.nona)
ltmod4<-lm(Light_Mean~stand*year+block*canopy*understory,data=clim.nona)
AICtab(ltmod,ltmod2,ltmod3,ltmod4)#ltmod4
#now estimate missing values
predict.lt<-predict(ltmod4,newdata=clim,type="response")
cbind(clim$Light_Mean,predict.lt)
plot(clim$Light_Mean,predict.lt)
abline(a=0,b=1,lty=1)
#light during days when temps above 5C
ltgddmod<-lm(Light_GDD~stand+block+year+canopy+understory,data=clim.nona)
ltgddmod2<-lm(Light_GDD~stand+block+year+canopy*understory,data=clim.nona)
ltgddmod3<-lm(Light_GDD~stand*year+block+canopy*understory,data=clim.nona)
ltgddmod4<-lm(Light_GDD~stand*year+block*canopy*understory,data=clim.nona)
AICtab(ltgddmod,ltgddmod2,ltgddmod3,ltgddmod4)#ltgddmod4
#now estimate missing values
predict.ltgdd<-predict(ltgddmod4,newdata=clim,type="response")
cbind(clim$Light_GDD,predict.ltgdd)
plot(clim$Light_GDD,predict.ltgdd)
abline(a=0,b=1,lty=1)
#GST
gstmod<-lm(GST_Mean~stand+block+year+canopy+understory,data=clim.nona)
gstmod2<-lm(GST_Mean~stand+block+year+canopy*understory,data=clim.nona)
gstmod3<-lm(GST_Mean~stand*year+block+canopy*understory,data=clim.nona)
gstmod4<-lm(GST_Mean~stand*year+block*canopy*understory,data=clim.nona)
AICtab(gstmod,gstmod2,gstmod3,gstmod4)#gstmod4
#now estimate missing values
predict.gst<-predict(gstmod4,newdata=clim,type="response")
cbind(clim$GST_Mean,predict.gst)
plot(clim$GST_Mean,predict.gst)
abline(a=0,b=1,lty=1)
##MAT
matmod<-lm(MAT~stand+block+year+canopy+understory,data=clim.nona)
matmod2<-lm(MAT~stand+block+year+canopy*understory,data=clim.nona)
matmod3<-lm(MAT~stand*year+block+canopy*understory,data=clim.nona)
matmod4<-lm(MAT~stand*year+block*canopy*understory,data=clim.nona)
AICtab(matmod,matmod2,matmod3,matmod4)#matmod4
#now estimate missing values
predict.mat<-predict(matmod4,newdata=clim,type="response")
cbind(clim$MAT,predict.mat)
plot(clim$MAT,predict.mat)
abline(a=0,b=1,lty=1)
#Snow duration (snow_dur)
snowmod<-lm(snow_dur~stand+block+year+canopy+understory,data=clim)
snowmod2<-lm(snow_dur~stand+block+year+canopy*understory,data=clim)
snowmod3<-lm(snow_dur~stand*year+block+canopy*understory,data=clim)
snowmod4<-lm(snow_dur~stand*year+block*canopy*understory,data=clim)
AICtab(snowmod,snowmod2,snowmod3,snowmod4)#snowmod4
#now estimate missing values
predict.snow<-predict(snowmod4,newdata=clim,type="response")
cbind(clim$snow_dur,predict.snow)
plot(clim$snow_dur,predict.snow)
abline(a=0,b=1,lty=1)
#Snow duration assuming continual cover from first to alst date (snow_dur_cont)
snow2mod<-lm(snow_dur_cont~stand+block+year+canopy+understory,data=clim)
snow2mod2<-lm(snow_dur_cont~stand+block+year+canopy*understory,data=clim)
snow2mod3<-lm(snow_dur_cont~stand*year+block+canopy*understory,data=clim)
snow2mod4<-lm(snow_dur_cont~stand*year+block*canopy*understory,data=clim)
AICtab(snow2mod,snow2mod2,snow2mod3,snow2mod4)#snow2mod4
#now estimate missing values
predict.snow2<-predict(snow2mod4,newdata=clim,type="response")
cbind(clim$snow_dur_cont,predict.snow2)
plot(clim$snow_dur_cont,predict.snow2)
abline(a=0,b=1,lty=1)



#Add predicted climate (gap-filled climate)to dataframe of observed climate
clim.gf<-cbind(clim,predict.mat,predict.gst,predict.ltgdd,predict.lt,predict.gdd,predict.gddn,predict.snow,predict.snow2)
clim.gf[which(is.na(clim.gf$GDD_total)),]$GDD_total<-clim.gf[which(is.na(clim.gf$GDD_total)),]$predict.gdd
clim.gf[which(is.na(clim.gf$GDD_totaln)),]$GDD_totaln<-round(clim.gf[which(is.na(clim.gf$GDD_totaln)),]$predict.gddn,digits=0)
clim.gf[which(is.na(clim.gf$Light_Mean)),]$Light_Mean<-clim.gf[which(is.na(clim.gf$Light_Mean)),]$predict.lt
clim.gf[which(is.na(clim.gf$Light_GDD)),]$Light_GDD<-clim.gf[which(is.na(clim.gf$Light_GDD)),]$predict.ltgdd
clim.gf[which(is.na(clim.gf$GST_Mean)),]$GST_Mean<-clim.gf[which(is.na(clim.gf$GST_Mean)),]$predict.gst
clim.gf[which(is.na(clim.gf$MAT)),]$MAT<-clim.gf[which(is.na(clim.gf$MAT)),]$predict.mat
clim.gf[which(is.na(clim.gf$snow_dur)),]$snow_dur<-clim.gf[which(is.na(clim.gf$snow_dur)),]$predict.snow
clim.gf[which(is.na(clim.gf$snow_dur_cont)),]$snow_dur_cont<-clim.gf[which(is.na(clim.gf$snow_dur_cont)),]$predict.snow2


dim(clim.gf)
clim.gf
head(clim.gf)

write.csv(clim.gf,"clim_gf.csv", row.names=FALSE)
