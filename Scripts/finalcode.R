data <- read.csv("Data/Clean IVI years 2013-2019.csv")

library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(vegan)
library(forcats)
library(stringr)
library(arm)
library(dplyr)
#log transform IVI
hist(data$Real.IVI)
data$logIVI <- log(data$Real.IVI)
hist(data$logIVI)

#create julian date which starts at 1 for the first hatch date for the full dataset (2015-2019)
summary(data$julian2)
data$jhatchdate <- data$julian2-186

#julian for date in season for full dataset (2015-2019)
summary(data$jdate2)
data$dayinseason <- data$jdate2-186

#calculate chick age 
data$chickage <- data$dayinseason - data$jhatchdate
hist(data$chickage)
summary(data$chickage)

plot(as.factor(data$chickage), data$Real.IVI)

#filter to include only chickages 1-12 
younger12 <-subset(data, data$chickage<=12)
only <-subset(younger12, younger12$chickage>0)

#filter again to include only unsupplemented nests 
only_unsupp <- only %>% 
  filter(supplimented != "y")#supplemented is spelled wrong intentionally, it's also wrong in the .csv and I haven't corrected it 

#count guilds per chickage
guilds_chickage <- only_unsupp %>% 
  group_by(year, chickage, guild) %>% 
  summarize(count=n())

#####this calculates evenness for each nestling age each year, I just need to run for each year. Later I will learn mapping 
##2013
#div13 <- guilds_chickage %>% 
 # filter(year == 2013) %>% 
#  pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div13[is.na(div13)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d <- diversity(div13)
#e <- d/log(specnumber(div13))
#print(e)

##2014
#div14 <- guilds_chickage %>% 
#  filter(year == 2014) %>% 
 # pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div14[is.na(div14)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d1 <- diversity(div14)
#e1 <- d1/log(specnumber(div14))
#print(e1)

##2015
#div15 <- guilds_chickage %>% 
#  filter(year == 2015) %>% 
 # pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div15[is.na(div15)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d2 <- diversity(div15)
#e2 <- d2/log(specnumber(div15))
#print(e2)

##2016
#div16 <- guilds_chickage %>% 
 # filter(year == 2016) %>% 
#  pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div16[is.na(div16)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d3 <- diversity(div16)
#e3 <- d3/log(specnumber(div16))
#print(e3)

##2017
#div17 <- guilds_chickage %>% 
 # filter(year == 2017) %>% 
  #pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div17[is.na(div17)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d4 <- diversity(div17)
#e4 <- d4/log(specnumber(div17))
#print(e4)

##2018
#div18 <- guilds_chickage %>% 
#  filter(year == 2018) %>% 
 # pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div18[is.na(div18)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d5 <- diversity(div18)
#e5 <- d5/log(specnumber(div18))
#print(e5)

##2019
#div19 <- guilds_chickage %>% 
#  filter(year == 2019) %>% 
#  pivot_wider(id_cols = "chickage", names_from = "guild", values_from = "count")
#div19[is.na(div19)] <- 0 ####replace na with 0 because the NA exists because there are no waterfowl in 2013 apparently 
#d6 <- diversity(div19)
#e6 <- d2/log(specnumber(div19))
#print(e6)

##print all of them as a fake table
#e
#e1
#e2
#e3
#e4
#e5
#e6




#calculate evenness (using vegan package)
#H <- diversity(guilds_chickage$count)
#e <- H/log(specnumber(guilds_chickage$count))
#this just gives 1 value for evenness but I want evenness for each year 

#
#evenness_year <- function(year){
 # e1 <- (guilds_chickage[guilds_chickage$year==year,]) 
 # H <- diversity(e1$count)
 # e <- H/log(specnumber(e1$count))
#  return(e)
#}

   # evenness_year("2013")
   # evenness_year("2014") ##call 2 to check they're different 
   
    
   # for(i in unique(guilds_chickage$year)){
    #  evenness_per_year <- evenness_year(i)
    #  print(evenness_per_year)
  #  } #make it a for loop that prints for each year 


    
##same thing for number of chicks     
   # guilds_chickn <- only_unsupp %>% 
   #   group_by(year, chicks, guild) %>% 
    #  mutate(count=n())
    
    
    #calculate evenness (using vegan package)
   # H2 <- diversity(guilds_chickn$count)
   # e2 <- H/log(specnumber(guilds_chickn$count))
    #this just gives 1 value for evenness but I want evenness for each year 
    
    #
   # evenness_year_n <- function(year){
    #  e1 <- (guilds_chickn[guilds_chickn$year==year,]) 
    #  H2 <- diversity(e1$count)
#e2 <- H2/log(specnumber(e1$count))
   #   return(e2)
  #  }
    
   # evenness_year_n("2013")
   # evenness_year_n("2014") ##call 2 to check they're different 
    
    
  #  for(i in unique(guilds_chickn$year)){
   #   evenness_per_year_n <- evenness_year_n(i)
#print(evenness_per_year_n)
  #  } #make it a for loop that prints for each year 
    
    
   # head(guilds_chickage, 15)
    
    ##viewed these things to realise my data was in the wrong format for what I was doing 
#data(BCI)
#View(BCI)

#parental level decision - logIVI
#modelled as a function of the year (random effect) as factor 
#other random effects: chicks, chickage
#fixed effects: yearsite 

no_unknown_settings <- only %>% 
  filter(motion != "unknown") %>% 
  filter(trigger != "unknown")


all_cam_settings <- lmer(logIVI ~ -1 + as.factor(year) + chicks + chickage + trigger + motion + (1|yearsite) + (1|site), data=only)
summary(all_cam_settings)

no_cam_settings <- lmer(logIVI ~ -1 + as.factor(year) + chicks + chickage + (1|yearsite)+ (1|site), data=only)
summary(no_cam_settings)

only_trigger_settings <- lmer(logIVI ~ -1 + as.factor(year) + chicks + chickage + trigger + (1|yearsite) + (1|site), data=only)
summary(only_trigger_settings)

only_motion_settings <- lmer(logIVI ~ -1 + as.factor(year) + chicks + chickage + motion + (1|yearsite) + (1|site), data=only)


anova(all_cam_settings, only_trigger_settings, only_motion_settings, no_cam_settings)


##hi Erik - here's what I'm trying to do. 

###If they are able to adjust based on their provisioning rates, they shouldn’t show variance prone behaviour. If they are not able to adjust their provisioning rates, they are more likely to show variance prone behaviour. 
#If I think that year specific conditions shape the options that are available then the prediction is about how the reaction norm of rate and variance correlate across years. 
#In a year where they don’t have scope to respond to increased brood demand, there should be higher variance (more variance prone behaviour). 
#As the slope gets higher i.e. if they are able to respond more strongly to variation in brood demand then their variance should decrease. 
#I have a study where I can estimate these values across 7 different years to get slopes and variances for each year. 
#My explicit prediction is that if this year-specific phenomenon is true there should be a negative covariance between the reaction norm slope and the residuals. 
#The model will include interactions. When you are talking about reaction norm slopes – so how do peregrines respond to variation in brood demand. We are explicitly interested in that slope in different years then that is an interaction between brood size and year or chick age and year. 

#The model would look like:
  

only$year <- as.factor(only$year)#code year as factor 
str(only$year)##check it worked 

IVI_log <- lmer(logIVI ~ -1 + chicks:year + chickage:year + (1|yearsite) + (1|site), data=only)# this works
summary(IVI_log)$sigma^2#this gives me one value, need 7

hist(resid(IVI_log)) #residuals are fairly normal, good
summary(IVI_log)
plot(IVI_log) #heteroskedastic. idk if thats bad I vaguely know thats whats happening based on a 5 min youtube video i just watched 
#var(ui|xi)=f(xi) -> variance of errors given xi is some function of xi (depends on xi) -- direction not clear 
###is that even right? not sure tbf. its blue. someone said that was good. 


##stolen code below 
require(MCMCglmm)

smod<-sim(IVI_log,1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))




#####actual model attempts
library(nlme)
only_unsupp$site <- as.factor(only_unsupp$site) 
only_unsupp$yearsite <- as.factor(only_unsupp$yearsite)
only_unsupp$year <- as.factor(only_unsupp$year)

ctrl <- lmeControl(opt="optim");

a1<- lme(logIVI ~ chickage:year + chicks:year,   
          random = ~ 1 | yearsite / site,           
          weights = varIdent(form = ~ 1|year),
          control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
          data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object

a2<- lme(logIVI ~ chickage:year + chicks:year,   
          random = list (site=~1, yearsite=~1),           
          weights = varIdent(form = ~ 1|year),
          control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
          data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object

a3<- lme(logIVI ~ chickage:year + chicks:year,   
          random = list (site=~1, yearsite=~1),           
          control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
          data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object

a4 <- lme(logIVI ~ chickage + chicks + year,   
          random = list (site=~1, yearsite=~1),   
          weights = varIdent(form = ~ 1|year),
          data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object
summary(a4)
anova(a1, a2, a3, a4)


ctrl <- lmeControl(opt="optim");
a2<- lme(logIVI ~ chickage + chicks+ year,    
         random = list (site=~1, yearsite=~1),           
         weights = varIdent(form = ~ 1|year),
         control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
         data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object
summary(a2)


ctrl <- lmeControl(opt="optim");
a3<- lme(logIVI ~ chickage:year + chicks:year,   
         random = list (site=~1, yearsite=~1),           
         control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
         data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object
summary(a3)


anova(a1, a3)


#with cam settings 
a4<- lme(logIVI ~ chickage:year + chicks:year +trigger,   
         random = list (site=~1, yearsite=~1),           
         weights = varIdent(form = ~ 1|year),
         control=ctrl, #this was to fix the error message -> nlminb problem, convergence error code = 1 message = iteration limit reached without convergence (10)
         data = only_unsupp, 
         na.action=na.exclude) #this was to fix the error message -> Error in na.fail.default(list(year = c(2015L, 2015L, 2015L, 2015L, 2015L,  :missing values in object
summary(a4)
anova(a4)



##USE THIS ONE :)

a5<- lme(logIVI ~ -1+chickage:year + chicks:year,   
         random = ~1|site/yearsite,          
         weights = varIdent(form = ~ 1|year),
         data = only_unsupp, 
         na.action=na.exclude) 
summary(a5)
anova(a5)

# plot residuals
plot(a5, resid(., type = "n") ~ chickage | year)
plot(a5, resid(., type = "n") ~ chicks | year)




##indep years
y2019 <- subset(only_unsupp, only_unsupp$year=="2019")
y2018 <- subset(only_unsupp, only_unsupp$year=="2018")
y2017 <- subset(only_unsupp, only_unsupp$year=="2017")
y2016 <- subset(only_unsupp, only_unsupp$year=="2016")
y2015 <- subset(d, d$year=="2015")
y2014 <- subset(only_unsupp, only_unsupp$year=="2014")
y2013 <- subset(d, d$year=="2013")
  

y2013_m <- lme(logIVI ~ chickage_s + chicks_s,   
         random = ~1|site,          
         data = y2015, 
         na.action=na.exclude) 

nd_2013 <- expand.grid(chicks_s = seq(min(y2013$chicks_s, na.rm = TRUE), max(y2013$chicks_s, na.rm = TRUE), length = 4),
                  chickage_s = seq(min(y2013$chickage_s), max(y2013$chickage_s), length = 12), 
                  site = unique(y2013$site))

preds_2013 <- predict(y2013_m, nd_2013, level = 0:1)
y2013_pred <- cbind(nd_2013, preds_2013)
preds_nlme <- predict(m, nd, level = 0:2)


summary(y2013)
summary(m)




anova(y2013)

y2014<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2014, 
            na.action=na.exclude) 
summary(y2014)
anova(y2014)

y2015<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2015, 
            na.action=na.exclude) 
summary(y2015)
anova(y2015)

y2016<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2016, 
            na.action=na.exclude) 
summary(y2016)
anova(y2016)

y2017<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2017, 
            na.action=na.exclude) 
summary(y2017)
anova(y2017)

y2018<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2018, 
            na.action=na.exclude) 
summary(y2018)
anova(y2018)

y2019<- lme(logIVI ~ -1+chickage + chicks,   
            random = ~1|site,          
            data = y2019, 
            na.action=na.exclude) 
summary(y2019)
anova(y2019)




###graphs - playing 

IVIA <- ggplot(data=only_unsupp, aes(x=chickage, y=logIVI)) +
  geom_point(aes(colour=year)) +
  stat_smooth(method="lm", se= TRUE,formula=y~x) +
  facet_wrap(~year)
print(IVIA)

IVIB <- ggplot(data=only_unsupp, aes(x=chicks, y=logIVI)) +
  geom_point(aes(colour=year)) +
  stat_smooth(method="lm", se= TRUE,formula=y~x) +
  facet_wrap(~year)
print(IVIB)
