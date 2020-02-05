#setwd("/Users/carlykenkel/Dropbox/NSFexpt")

#if library command below does not work
install.packages("piecewiseSEM")
source("http://bioconductor.org/biocLite.R")
biocLite("ggplot2")


library(ggplot2)
library(car)
library(multcomp)
library(gtools)
library(plyr)
library(quantreg)
library(calibrate)
library(MASS)
library(AICcmodavg)
library(e1071)
library(nlme)
library(lmmfit)
library(MCMCglmm)
library(labdsv)
library(vegan)
library(plotrix)
library(pgirmess)
library(gridExtra)
library(lme4)
source("RsquaredGLMM.R")

library(survival) # for survival analysis 
library(vegan) # for plotting ellipses in principal coordinates analysis plots
library(corrr) # for correlations
library(MCMCglmm) # for mcmcglmm stats
library(car) # for data transformations 
library(nlme) #for lme
library(MASS) # for stepAIC
library(coxme) # for mixed effects cox model
library(tidyverse) # for data wrangling and visualization
library(ggridges) # for ridge plots
library(reshape2) # for melt
library(corrplot) # for correlations
library(summarytools) # for dfSummary
library(ggplot2)


#############################

#Summary Function - execute this entire block of code to store

#############################

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###############################PAM data
setwd("/Users/carlykenkel/Dropbox/CarlsLab/ResearchProjects/NSF_RAPID_Irma")
datAll=read.csv("IrmaSurveys.csv")
summary(datAll)
head(datAll)
datNoSid<-subset(datAll,Species!="Ssid")
dat<-subset(datNoSid,InTransect=="Y")
summary(dat)


Buoy10=subset(dat,Site=="LK_Buoy10")
Buoy12=subset(dat,Site=="LK_Buoy12")
Buoy6=subset(dat,Site=="LK_Buoy6")


summary(Buoy10) #27 frags
summary(Buoy12) #11 frags
summary(Buoy6) #8 frags

par(mfrow=c(1,2))
plot(Attached~Site,dat)
plot(Attached~Zone,dat)

insh<-subset(dat,dat$Zone=="Inshore")
offsh<-subset(dat,dat$Zone=="Offshore")

#Do g-test for association of two categorical variables (dislodged by zone)
frag.df = matrix(c(133,106,102,47), nrow = 2)

frag.df
fisher.test(frag.df)
G.test(frag.df)

102/(102+133)
47/(47+106)


par(mfrow=c(2,1))
plot(Attached~Species, insh,ylab="inshore")
plot(Attached~Species, offsh,ylab="offshore")


frags<-subset(datNoSid,Attached=="Dislodged")

fragsCC<-frags[complete.cases(frags$tod),]

deadFrags=Surv(fragsCC$tod, fragsCC$dead)


#colorR <- c("black", "orangered", "orange", "red", "blue")
colorS <- c("blue", "lightblue", "purple")
site <- c("LK_Buoy10", "LK_Buoy12", "LK_Buoy6")

sur1site = survfit(deadFrags ~ Site, data=fragsCC)
sur1site
# 1-(0/25) # 100% survival of geno 1 (0% mortality)
# 1-(1/26) # 96% survival of geno 3 (4% mortality)


plot(sur1site, col=colorS, lwd=5, 
     xlab="Time in Months", ylab="Fraction Surviving", main="Fragment Survival by Site", ylim=c(0.7,1))
legend("bottomleft", col=colorS, legend=site, lwd=5, bty="n")


