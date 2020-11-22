library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(R2jags)
d <- read_csv("Data/Clean IVI years 2013-2019.csv")
# Data Load/prep ----------------------------------------------------------

#glimpse(d)
d <- d %>% filter(supplimented == "n")
d$logIVI <- log(d$ivi)
d$hatch_date <- date(d$hatch_date)
d$chickage <- as.Date(d$date) - as.Date(d$hatch_date)
d <- d %>% filter(chickage < 13)
d$year <- as.factor(year(d$date))
d$site <- as.factor(d$site)
d$chickage <- d$chickage-1










# change to factors, and calculate unique levels of randoms
nestID <- as.factor(d$site)
n_nests <- length(levels(nestID))

yearsite_f <- as.factor(d$yearsite)
n_yearsites <- length(levels(yearsite_f))

n_years <- length(levels(d$year))

d$chickage <- d$chickage - 1

# model matrix used to store betas
X <- model.matrix(~ 1 + chickage:year + chicks:year, data = d)



# data for Jags model
jags_data <- list(y = d$logIVI,     # ivi 
                 years = d$year, # year identifier for variance
                 nest = nestID,              # random intercept for nest
                  year = d$year,
                  yearsite = yearsite_f,      # random intercept for yearsite
                  n_years = n_years,          # number of years
                  n_nests = n_nests,          # number of nests
                  n_yearsites = n_yearsites,  # number of unique yearsites
                  X = X,                      # intercept + covariates (model matrix)
                  N = nrow(d),      # sample size
                  K = ncol(X))                # Number of betas


#~~~~~~~~ Jags Model ~~~~~~~~~~~~#
sink("het_var.txt")
cat("
    model{
    
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,])  + g[yearsite[i]] + a[nest[i]] #+ year_int[year[i]]
        
        # store predicted values at each iteration in y_pred
        y_pred[i] ~ dnorm(mu[i], tau[years[i]])
        }
    
    # Priors ~~~~~~~~~~~~~~~
        
     # priors for betas
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
        
     # prior for residual variance weighting matrix
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        tau[i]   <- pow(sigma[i], -2)
        }

    
     # prior for random intercepts, a = site, g = yearsite
       # for (i in 1:n_years) {year_int[i] ~ dnorm(year_bar, sigma_year)}
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
        for (i in 1:n_yearsites) {g[i] ~ dnorm(g_bar, sigma_yearsite)}
    
     # prior for mean/variance of random intercepts
       # year_bar ~ dnorm(0, 1.5)
       # sigma_year ~ dexp(1)
        a_bar ~ dnorm(0, 1.5)
        sigma_nest ~ dexp(1)
        g_bar ~ dnorm(0, 1.5)
        sigma_yearsite ~ dexp(1)
        
    }
", fill = TRUE)
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Run Model ---------------------------------------------------------------

# Store draw information from the folowing parms
params <- c("beta", "sigma","g", "g_bar", "sigma_yearsite", "year_int")

het_m1   <- jags(data      = jags_data,
                 inits      = NULL,     # runs ok w/o starting values, change if things become more complex
                 parameters = params,
                 model      = "het_var.txt",
                 n.thin     = 10, 
                 n.chains   = 3,
                 n.burnin   = 4000,
                 n.iter     = 5000)

het_m2 <- update(het_m1, n.iter = 10000, n.thin = 10) 
het_m3 <- update(het_m2, n.iter = 20000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 20000, n.thin = 10) 


# convert our jags output into mcmc object
het_mcmc <- as.mcmc(het_m3)




# Model Diagnostics -------------------------------------------------------

# mixing
het_mcmc %>%
    window(thin=10) %>% 
    
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], sigma[i], g[i]) %>%
    
    # change this filter to look at mixing for the parameter of interest. See above line for options
    filter(.variable == "beta") %>% 
    ungroup() %>%
    mutate(term = ifelse(is.na(i), .variable, paste0(.variable,"[",i,"]"))) %>%
    ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
    scale_color_manual(values=c("#461220", "#b23a48", "#fcb9b2")) +
    geom_line(alpha=0.5) +
    facet_grid(term~., scale="free_y") +
    labs(color="chain", x="iteration") +
    nuwcru::theme_nuwcru()

# ACF

# more diagnostics, to do
x <- het_mcmc %>%
    window(thin=10) %>% 
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], sigma[i], g[i], a[i])

unique(x$.variable)
unique(test$i)

beta <- x %>% filter(.variable == "beta") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(beta$i))){levels(beta$i)[i] <- colnames(X)[i]}
levels(beta$i)[1] <- "intercept"

sigma <- x %>% filter(.variable == "sigma") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(sigma$i))){levels(sigma$i)[i] <- paste0("sigma", "_", levels(d$year)[i])}

a <- x %>% filter(.variable == "a") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(a$i))){levels(a$i)[i] <- unique(as.character(d$site))[i]}

g <- x %>% filter(.variable == "g") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(g$i))){levels(g$i)[i] <- unique(as.character(d$yearsite))[i]}

year <- x %>% filter(.variable == "year") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(year$i))){levels(year$i)[i] <- paste0("int", "_", levels(d$year)[i])}
unique(year$i)

all <- rbind(beta, sigma, a, g)
names(all) <- c("level", "chain", "iteration", "draw", "parameter", "value")



arrow::write_parquet(all, "data/jags_out.parquet")

mutate(i = case_when(
        i == 1 ~ "intercept",
        i == 2 ~ tolower(colnames(X)[2]),
        i == 3 ~ colnames(X)[3],
        i == 4 ~ colnames(X)[4],
        i == 5 ~ colnames(X)[5],
        i == 6 ~ colnames(X)[6],
        i == 7 ~ colnames(X)[7],
        i == 8 ~ colnames(X)[8],
        i == 9 ~ colnames(X)[9],
        i == 10 ~ colnames(X)[10],
        i == 11 ~ colnames(X)[11],
        i == 12 ~ colnames(X)[12],
        i == 13 ~ colnames(X)[13],
        i == 14 ~ colnames(X)[14],
        i == 15 ~ colnames(X)[15]
        ))

sigma <- x %>% filter(.variable == "sigma") %>% group_by(i) %>% tally()
    mutate(i = case_when(
        i == 1 ~ "intercept",
        i == 2 ~ tolower(colnames(X)[2]),
        i == 3 ~ colnames(X)[3],
        i == 4 ~ colnames(X)[4],
        i == 5 ~ colnames(X)[5],
        i == 6 ~ colnames(X)[6],
        i == 7 ~ colnames(X)[7],
        i == 8 ~ colnames(X)[8],
        i == 9 ~ colnames(X)[9],
        i == 10 ~ colnames(X)[10],
        i == 11 ~ colnames(X)[11],
        i == 12 ~ colnames(X)[12],
        i == 13 ~ colnames(X)[13],
        i == 14 ~ colnames(X)[14],
        i == 15 ~ colnames(X)[15]
    ))


# Effects Plotting --------------------------------------------------------

# join data set with model predictions
het_d <- het_mcmc %>%
    tidybayes::spread_draws(mu[i], y_pred[i]) %>%
    ungroup() %>%
    left_join(
        mutate(d, i = 1:n())) %>%
    mutate(resid = logIVI - mu)



# Predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log ivi ~ chickage
# dimension exported: 1200 x 804
ggplot(d, aes(x=chickage, y=logIVI)) +
    tidybayes::stat_lineribbon(data=het_d, aes(y=y_pred), colour = "#08519C") +
    scale_fill_brewer() + 
    geom_point(alpha = 0.5) +
    facet_wrap(~year) +
    ylab("log IVI") + xlab("Chickage") +
    theme_nuwcru() + facet_nuwcru() +
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_blank(), 
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank())
display.brewer.pal(n = 3)

# log ivi ~ broodsize
# dimension exported: 1200 x 804
ggplot(d, aes(x=chicks, y=logIVI)) +
    tidybayes::stat_lineribbon(data=het_d, aes(y=y_pred), colour = "#08519C") +
    scale_fill_brewer() + 
    geom_point(alpha = 0.5) +
    facet_wrap(~year) +
    ylab("log IVI") + xlab("brood size") +
    theme_nuwcru() + facet_nuwcru() +
    theme(panel.border = element_blank(), axis.line.y = element_blank(), axis.line.x = element_blank(), 
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank())



# Sigmas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sigma as a ratio to 2013 (reference)
# dimension exported: 1200 x 804
het_mcmc %>%
    tidybayes::spread_draws(sigma[i]) %>%
    group_by(i) %>%
    summarize(mean = getmode(sigma), sd = sd(sigma)) %>%
    mutate(upper95 = mean + (1.96 * sd), lower95 = mean - (1.96 * sd),
           upper80 = mean + (1.282 * sd), lower80 = mean - (1.282 * sd),
           upper50 = mean + (0.674 * sd), lower50 = mean - (0.674 * sd)) %>%
    ggplot() +
    geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
    geom_vline(xintercept = 1, colour = grey6, linetype = "dashed") +
    geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
    geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
    geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
    geom_point(aes(y = 2013:2019, x = mean), shape = "|",  colour = "white", size = 5)+
    
    # lme4 models
    geom_segment(data = nd, aes(x = lower95, xend = upper95, y = year+0.2, yend = year+0.2), colour = "#e2b6b6", size = 2, alpha = 0.65) +
    geom_segment(data = nd, aes(x = lower80, xend = upper80, y = year+0.2, yend = year+0.2), colour = "#c76f6f", size = 2, alpha = 0.65) +
    geom_segment(data = nd, aes(x = lower50, xend = upper50, y = year+0.2, yend = year+0.2), colour = "#7f1111", size = 2, alpha = 0.65) +
    geom_point(data = nd, aes(x = post_mode, y = year+0.2), shape = "|",  colour = "white", size = 5) +
    ylab("") + xlab("Sigma yearly estimate (95% CI)") +
    scale_x_continuous(limits = c(0.57,1.43)) +
    scale_y_continuous(breaks = 2013:2019) +
    theme_nuwcru() + 
    theme(panel.border = element_blank(),
          #axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
    

# Betas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gather betas, calculate mean and credibles
t <- het_mcmc[1] %>%
    tidybayes::spread_draws(beta[i]) %>%
    group_by(i) %>%
    summarize(mean = mean(beta), sd = sd(beta)) %>%
    mutate(upper95 = mean + (1.96 * sd), lower95 = mean - (1.96 * sd),
           upper80 = mean + (1.282 * sd), lower80 = mean - (1.282 * sd),
           upper50 = mean + (0.674 * sd), lower50 = mean - (0.674 * sd)) %>%
    mutate(param = colnames(X))

# betas for chickage
t[2:8,] %>% 
    ggplot() +
    geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
    geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
    geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
    geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
    geom_point(aes(y = 2013:2019, x = mean), shape = 21, fill = "white", colour = "black", size = 4) +
    scale_x_continuous(limits = c(-0.100,0.100)) +
    scale_y_continuous(breaks = 2013:2019) +
    ylab("") + xlab("Chickage : Year") + 
    theme_nuwcru() + 
    theme(panel.border = element_blank(),
          #axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())


# betas for broodsize
t[9:15,] %>% 
    ggplot() +
    geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
    geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
    geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
    geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
    geom_point(aes(y = 2013:2019, x = mean), shape = 21, fill = "white", colour = "black", size = 4) +
    scale_x_continuous(limits = c(-0.3,0.300)) +
    scale_y_continuous(breaks = 2013:2019) +
    ylab("") + xlab("Broodsize : Year") + 
    theme_nuwcru() + 
    theme(panel.border = element_blank(),
          #axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
