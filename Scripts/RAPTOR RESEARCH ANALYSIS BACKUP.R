
##loading data 
rm(list = ls())

setwd("F:/MASTERS YEAR 1/Raptor Research")

deliveriesday<-read.csv("N DELIVERIES.csv")


##packages installed 
library(ggplot2)
library(Hmisc)
library(dplyr)
library(lme4)


#SUMMARY OF WORK FROM 28.10.19: PLOTTED MEANS CONFIDENCE INTERVALS AND SD FOR DELIVERIES PER DAY AND IVIs THEN DABBLED IN HOW TO CALCULATE MEANS AND SDs USING CODE AND MAKING THEM INTO THEIR OWN TABLES

##Deliveries Per Day 
#Plot of SITE/DELIVERIES PER DAY with means, confidence intervals (red) and SD (black)

fig = ggplot(deliveriesday, aes(x=SITE, y=N.DELIVERIES)) +
  stat_summary(fun.y = mean, geom = "point", show.legend = T) +
  stat_summary(fun.data = mean_cl_boot, geom = "crossbar", width=0.2,colour="red")+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", width=0.2)+  
  labs(x="SITE",y="DELIVERIES PER DAY")
fig

#mean and SD tables 
data <- read.csv("N DELIVERIES.csv")
deliveriesdaySD <- data %>% group_by(DATE) %>% summarize(sd = sd(N.DELIVERIES), mean = mean(N.DELIVERIES))


##IVIs
#Plot of SITE/IVIs with means, confidence intervals (red) and SD (black)
IVIs<-read.csv("IVIs.csv")

fig = ggplot(IVIs, aes(x=SITE, y=IVInum)) +
  stat_summary(fun.y = mean, geom = "point", show.legend = T) +
  stat_summary(fun.data = mean_cl_boot, geom = "crossbar", width=0.2,colour="red")+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", width=0.2)+  
  labs(x="SITE",y="IVIs")
fig

#mean and SD tables --- couldnt calculate means/SDs because of missing values
IVIs <- read.csv("IVIs.csv")
IVISD <- IVIs %>% group_by(DATE) %>% summarize(sd = sd(IVInum), mean = mean(IVInum))


#mean and SD tables --- saved as 'edited' version with all blanks changed to NA - now works 
IVIs <- read.csv("IVIsedited.csv")
IVISDedited <- IVIs %>% group_by(DATE) %>% summarize(sd = sd(IVInum), mean = mean(IVInum))


---------------------------------------------------------------------------------------------


##SUMMARY OF WORK 29.10.2019: RAPTOR RESEARCH CONFERENCE GRAPHS AND BASIC MODELS FOR DELIVERIES PER DAY AND INTERVISIT INTERVALS AGROSS NESTLING AGE AND NUMBER OF NESTLINGS (INCLUDING INTERACTIONS)


##running models and making graphs for deliveries per day 

 ##histogram deliveries per day 

hist(data$N.DELIVERIES)

data$logdeliveries <- log(data$N.DELIVERIES +0.01)
hist(data$logdeliveries)

  #m1 non-interaction
m1 <- lmer (logdeliveries~N.CHICKS + AGE.CHICKS + (1|SITE), data=data)
names(data)
summary(m1)

# m2 interaction 

m2 <- lmer (logdeliveries~N.CHICKS * AGE.CHICKS + (1|SITE), data=data)
summary(m2)

#when more chicks the effect of age of chicks is stronger - effect of age most pronounced when many offspring 

plot(as.factor(data$N.CHICKS), data$N.DELIVERIES, xlab="NUMBER OF NESTLINGS", ylab ="NUMBER OF VISITS PER DAY")
##use for poster 

plot(as.factor(data$AGE.CHICKS), data$N.DELIVERIES, xlab = "AGE OF CHICKS", ylab = "NUMBER OF VISITS PER DAY",xlim=c(1,29.5))
##LIMIT TO BETWEEN 1-30
##pull out nests that dont survive??

data$deliveriesperchick <- data$N.DELIVERIES/data$N.CHICKS
hist(data$deliveriesperchick)
plot(data$AGE.CHICKS, data$deliveriesperchick)

##tell ggplot to fit the best line (not necessarily linear) - also fix axis titles 
bestfit <- ggplot(data, aes(x=AGE.CHICKS, y=deliveriesperchick, colour = SITE)) + 
  geom_point() + 
  geom_smooth(method = "auto") +
  ##limit deliveries per day to above 0 beacuse will not be less than 0 deliveries 
  coord_cartesian(ylim= c(0, 11))

bestfit

## plotting log to see if i like it better than non-log
bestfit2 <- ggplot(data, aes(x=AGE.CHICKS, y=logdpc)) + 
  geom_point() + 
  geom_smooth(method = "auto") +
  theme_classic() +
  labs(x="\n AGE OF CHICKS", y="LOG (DELIVERIES PER CHICK)\n")


bestfit2

##correlation analysis for deliveries per chick 

cor.test(data$AGE.CHICKS, data$deliveriesperchick)
##very high correlation for ecological data - as chicks get older per chick feeding rate increases



##m3
data$logdpc <- log(data$deliveriesperchick +0.01)
hist (data$logdpc)

m3 <- lmer (logdpc ~ AGE.CHICKS + (1|SITE), data=data)
summary(m3)
##repeatable, significant
##overall linear effect - increasing at 0.01 per dat, sig as t value = 2.9

##m4
data$age2 <- (data$AGE.CHICKS)^2
m4 <- lmer (logdpc ~ age2 + AGE.CHICKS + (1|SITE), data=data)
##squared term for age of chicks
summary(m4)

##overall delivery rate per chick increasing as chicks get older (age chicks line), quadratic term also significant and negative meaning the lowest and highest and highest values of age have lower feeding rates than the middle (peak not super strong) and linear effect 



    ##running same stuff but for IVIs
##histogram IVIs
names(IVIs)

hist(IVIs$IVInum)

IVIs$logIVInum <- log(IVIs$IVInum +0.01)

hist(IVIs$logIVInum)

#IVI model 1 (IVIm1)-IVI non-interaction
IVIm1 <- lmer (logIVInum~N.CHICKS + AGE.CHICKS + (1|SITE), data=IVIs)
summary(IVIm1)

# IVI m2 interaction (IVIm2)- IVI interaction effect 

IVIm2 <- lmer (logIVInum~N.CHICKS * AGE.CHICKS + (1|SITE), data=IVIs)
summary(IVIm2)
##reapeatibility (variance of site/(variance of site+residual variance) is just 0.040014752 so low repeatbility of IVIs within each site
##when you increase number of chicks, effect of age is slightly weaker????

plot(as.factor(IVIs$N.CHICKS), IVIs$IVInum, xlab="NUMBER OF NESTLINGS", ylab ="INTER-VISIT INTERVALS")
##use for poster 

plot(as.factor(IVIs$AGE.CHICKS), IVIs$IVInum, xlab = "AGE OF NESTLINGS", ylab = "INTER-VISIT INTERVALS",xlim=c(1,28.5)) 

##next steps: calcuate numbers of prey in each guild and plot against nestlinga age and number of nestlings as above
##Plot prey type as function of provisioning stage (early, middle, late, 10 day intervals of age) - so that i can then explain with, e.g., "proportion of birds in diet increases with age"

  