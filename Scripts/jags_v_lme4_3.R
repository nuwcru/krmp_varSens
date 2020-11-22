library(tidyverse)
library(lme4)
library(arm)

# function to calculate mode
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


# load data ---------------------------------------------------------------

d <- read_csv("data/Clean IVI years 2013-2019.csv") 
#glimpse(d)
d <- d %>% filter(supplimented == "n")
d$logIVI <- log(d$ivi)
d$hatchdate <- 
d$dayinseason <- d$jdate2-186
d$chickage <- d$dayinseason - d$jhatchdate
d <- d %>% filter(chickage < 13)
d$year <- as.factor(d$year)
d$site <- as.factor(d$site)
d$chickage <- d$chickage-1


# Jags summary ------------------------------------------------------------


# I've compressed the jags output into a parquet to save space
jags <- arrow::read_parquet("data/jags_out.parquet")



# these are the estimated parameters available to look at
  # sigma is the yearly residual variance
  # a is the random intercept for nest
  # g is the random intercept for yearsite
unique(jags$parameter)

jags <- all

# and you can even pull out specific levels of the parameter, so a nest site, or year
# these are all the combinations

# betas
jags %>% 
  filter(parameter == "beta" ) %>% 
  distinct(level)

# sigma - yearly residual variance
jags %>% 
  filter(parameter == "sigma" ) %>% 
  distinct(level)

# a - nest site
jags %>% 
  filter(parameter == "g" ) %>% 
  distinct(level)

# all values in the posterior distribution for each of these is in the "value" column
head(jags)

# we can summarize the posterior distributions, mode, 95, 80, and 50% credibles
jags_summary <- jags %>%
                  group_by(parameter, level) %>%
                  summarize(mode = getmode(value),
                         sd = sd(value)) %>%
                  mutate(upper95 = mode + (sd*1.96),
                         upper80 = mode + (sd*1.282),
                         upper50 = mode + (sd*0.674),
                         lower95 = mode - (sd*1.96),
                         lower80 = mode - (sd*1.282),
                         lower50 = mode - (sd*0.674))

# Beta chicks for each year
jags_summary %>% 
  filter(parameter == "beta") %>% filter(str_detect(level, "chicks")) 

# Beta chickage
jags_summary %>% 
  filter(parameter == "beta") %>% filter(str_detect(level, "chickage")) 

library(nuwcru)
# and plot chickage as an example
jags_summary %>% 
  filter(parameter == "beta") %>% filter(str_detect(level, "chickage")) %>%
ggplot() +
  geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
  geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
  geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
  geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
  geom_point(aes(y = 2013:2019, x = mode), shape = 21, fill = "white", colour = "#4A84BD", size = 3) +
  #scale_x_continuous(limits = c(-0.100,0.100)) +
  scale_y_continuous(breaks = 2013:2019) +
  ylab("") + xlab("Chickage : Year") + 
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())







# lme4 yearly models ----------------------------------------------

library(lme4)
# Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# model each year in lme4, and store in a list
year <- d %>% group_split(year)

model_outs <- lapply(year, function(x){ 
  lmer(logIVI ~ chicks + chickage + (1|site), x)})

# year is a list with 7 slots, each storing the model output of a given year
model_outs[[1]] # 2013 for example
model_outs[[2]] # 2014 as another example




library(arm)
# Simulate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# simulate all years in the list
sim <- lapply(model_outs, function(x){
  sim(x, n.sims=1000)})



# pull out yearly residual variance ~~~~~~~~~~~~~~~~~
yearly_resid <- lapply(sim, function(x){as.mcmc(x@sigma^2)})  # raw posteriors


nd <- data.frame(year     = 2013:2019,
                 post_mode= rep(NA, length(2013:2019)),
                 sd       = rep(NA, length(2013:2019)),
                 lower95  = rep(NA, length(2013:2019)),
                 lower80  = rep(NA, length(2013:2019)),
                 lower50  = rep(NA, length(2013:2019)),
                 upper95  = rep(NA, length(2013:2019)),
                 upper80  = rep(NA, length(2013:2019)),
                 upper50  = rep(NA, length(2013:2019)))

for (i in 1:length(yearly_resid)){
  nd[i,]$post_mode <- getmode(as.numeric(yearly_resid[[i]]))
  nd[i,]$sd <- sd(as.numeric(yearly_resid[[i]]))
  nd[i,]$lower95 <- nd[i,]$post_mode - (1.96*nd[i,]$sd)
  nd[i,]$lower80 <- nd[i,]$post_mode - (1.282*nd[i,]$sd)
  nd[i,]$lower50 <- nd[i,]$post_mode - (0.674*nd[i,]$sd)
  nd[i,]$upper95 <- nd[i,]$post_mode + (1.96*nd[i,]$sd)
  nd[i,]$upper80 <- nd[i,]$post_mode + (1.282*nd[i,]$sd)
  nd[i,]$upper50 <- nd[i,]$post_mode + (0.674*nd[i,]$sd)
}






# Pull out and summarize betas ~~~~~~~~~~~~~~~~~~~~~
betas <- lapply(sim, function(x){as.mcmc(x@fixef)}) # raw posteriors


temp <- list()
for (i in 1:length(betas)){
  temp[[i]] <- betas[[i]] %>% 
    tidybayes::gather_draws(`(Intercept)`, chicks, chickage) %>%
    group_by(.variable) %>%
    summarize(mode = getmode(.value),
              sd = sd(.value)) %>%
    mutate(upper95 = mode + (sd*1.96),
           upper80 = mode + (sd*1.282),
           upper50 = mode + (sd*0.674),
           lower95 = mode - (sd*1.96),
           lower80 = mode - (sd*1.282),
           lower50 = mode - (sd*0.674)) %>%
    mutate(year = c(2013:2019)[i])
}

betas_summary <- bind_rows(temp) # summary of beta posteriors
betas_summary %>% filter(.variable == "(Intercept)") %>%
  ggplot() +
  geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
  # lme4 models
  geom_segment(aes(x = lower95, xend = upper95, y = year, yend = year), colour = "#e2b6b6", size = 2, alpha = 0.65) +
  geom_segment(aes(x = lower80, xend = upper80, y = year, yend = year), colour = "#c76f6f", size = 2, alpha = 0.65) +
  geom_segment(aes(x = lower50, xend = upper50, y = year, yend = year), colour = "#7f1111", size = 2, alpha = 0.65) +
  geom_point(aes(x = mode, y = year), shape = 21, fill = "white", colour = "#7f1111", size = 2) +
  ylab("") + xlab("Intercept") +
  #scale_x_continuous(limits = c(0.57,1.43)) +
  scale_y_continuous(breaks = 2013:2019) +
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

  

# Didn't summarize random - we can do it later ~~~~~~~~~~~~~~
random <- lapply(sim, function(x){
  as.mcmc(as.vector(apply(x@ranef$site[,,1], 1, var)))})





# Comparisons -------------------------------------------------------------


# yearly residual variance comparison - jags in blue, lme4 in red
jags_summary %>% 
  filter(parameter == "sigma") %>%
  ggplot() +
    geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
    geom_vline(xintercept = 1, colour = grey6, linetype = "dashed") +
    geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
    geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
    geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
    geom_point(aes(y = 2013:2019, x = mode), shape = 21, fill = "white", colour = "#4A84BD", size = 2)+
    
    # lme4 models
    geom_segment(data = nd, aes(x = lower95, xend = upper95, y = year+0.2, yend = year+0.2), colour = "#e2b6b6", size = 2, alpha = 0.65) +
    geom_segment(data = nd, aes(x = lower80, xend = upper80, y = year+0.2, yend = year+0.2), colour = "#c76f6f", size = 2, alpha = 0.65) +
    geom_segment(data = nd, aes(x = lower50, xend = upper50, y = year+0.2, yend = year+0.2), colour = "#7f1111", size = 2, alpha = 0.65) +
    geom_point(data = nd, aes(x = post_mode, y = year+0.2), shape = 21, fill = "white", colour = "#7f1111", size = 2) +
    ylab("") + xlab("Sigma yearly estimate (95% CI)") +
    scale_x_continuous(limits = c(0.57,1.43)) +
    scale_y_continuous(breaks = 2013:2019) +
    theme_nuwcru() + 
    theme(panel.border = element_blank(),
          #axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())


# chickage betas
jags_summary %>% 
  filter(parameter == "beta") %>% filter(as.integer(level) > 1 & as.integer(level) < 9) %>%
  ggplot() +
  geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
  geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
  geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 1) +
  geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 1.5) +
  geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
  #geom_point(aes(y = 2013:2019, x = mode), shape = 21, fill = "white", colour = "#4A84BD", size = 2)+
  
  # lme4 models
  geom_segment(data = filter(betas_summary, .variable == "chickage"), aes(x = lower95, xend = upper95, y = year+0.2, yend = year+0.2), colour = "#e2b6b6", size = 1, alpha = 0.65) +
  geom_segment(data = filter(betas_summary, .variable == "chickage"), aes(x = lower80, xend = upper80, y = year+0.2, yend = year+0.2), colour = "#c76f6f", size = 1.5, alpha = 0.65) +
  geom_segment(data = filter(betas_summary, .variable == "chickage"), aes(x = lower50, xend = upper50, y = year+0.2, yend = year+0.2), colour = "#7f1111", size = 2, alpha = 0.65) +
  #geom_point(data = filter(betas_summary, .variable == "chickage"), aes(x = mode, y = year+0.2), shape = 21, fill = "white", colour = "#7f1111", size = 2) +
  ylab("") + xlab("Chickage:year") +
  #scale_x_continuous(limits = c(-0.15,0.15)) +
  scale_y_continuous(breaks = 2013:2019) +
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


# chicks - some differences here - jags is blue, lme4 is red
jags_summary %>% 
  filter(parameter == "beta") %>% filter(as.integer(level)>8) %>%
  ggplot() +
  geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
  geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
  geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 1) +
  geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 1.5) +
  geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
  geom_point(aes(y = 2013:2019, x = mode), shape = "|",  colour = "white", size = 5)+
  
  # lme4 models
  geom_segment(data = filter(betas_summary, .variable == "chicks"), aes(x = lower95, xend = upper95, y = year+0.2, yend = year+0.2), colour = "#e2b6b6", size = 1, alpha = 0.65) +
  geom_segment(data = filter(betas_summary, .variable == "chicks"), aes(x = lower80, xend = upper80, y = year+0.2, yend = year+0.2), colour = "#c76f6f", size = 1.5, alpha = 0.65) +
  geom_segment(data = filter(betas_summary, .variable == "chicks"), aes(x = lower50, xend = upper50, y = year+0.2, yend = year+0.2), colour = "#7f1111", size = 2, alpha = 0.65) +
  geom_point(data = filter(betas_summary, .variable == "chicks"), aes(x = mode, y = year+0.2), shape = 21,  colour = "white", size = 2) +
  ylab("") + xlab("Broodsize:Year") +
  #scale_x_continuous(limits = c(-0.15,0.15)) +
  scale_y_continuous(breaks = 2013:2019) +
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


# next step is to do posterior check comparisons (compare fitted with real values)
# but This is good for now. We can take this further after your meeting
