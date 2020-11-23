library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(R2jags)
library(tidybayes)

d <- read_csv("Data/Clean IVI years 2013-2019.csv")

# Data Load/prep ----------------------------------------------------------

# run ivi_scratch first
d <- d %>% filter(ivi < 2000)

d$year <- as.factor(year(d$date))
d$site <- as.factor(d$site)


ivi <- d %>% filter(ivi < 1440)

hist(ivi$ivi)

# change to factors, and calculate unique levels of randoms
nestID <- as.factor(d$site)
n_nests <- length(levels(nestID))

yearsite_f <- as.factor(d$yearsite)
n_yearsites <- length(levels(yearsite_f))

n_years <- length(levels(d$year))

d$chickage <- d$chickage - 1

# model matrix used to store betas

X <- model.matrix(~ 1 + chickage + chicks, data = d)

# data for Jags model
jags_data <- list(y = d$ivi,     # ivi 
                  years = d$year, # year identifier for variance
                  nest = nestID,              # random intercept for nest
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
    
    # Likelihood ###################################
    
        for (i in 1:N) {
          y[i]  ~ dgamma(r[years[i]], mu.eff[i])
          mu.eff[i]  <- r[years[i]] / mu[i] 
          log(mu[i]) <- inprod(beta[], X[i,]) + g[yearsite[i]] + a[nest[i]]
        }
    
    # Priors ######################################
        
     # priors for betas ###
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
        
     # prior for residual variance weighting matrix ###
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        r[i]     <- pow(sigma[i], -2)
        }

    
     # prior for random intercepts, a = site, g = yearsite ###
       # for (i in 1:n_years) {year_int[i] ~ dnorm(year_bar, sigma_year)}
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
        for (i in 1:n_yearsites) {g[i] ~ dnorm(g_bar, sigma_yearsite)}
    
     # prior for mean/variance of random intercepts ###
        a_bar ~ dnorm(0, 1.5)
        sigma_nest ~ dexp(1)
        g_bar ~ dnorm(0, 1.5)
        sigma_yearsite ~ dexp(1)
        
    }
", fill = TRUE)
sink()



# Store draw information from the folowing parms
params <- c("beta", "r","g", "g_bar", "sigma_yearsite")

het_m1   <- jags(data      = jags_data,
                 inits      = NULL,     # runs ok w/o starting values, change if things become more complex
                 parameters = params,
                 model      = "het_var.txt",
                 n.thin     = 10, 
                 n.chains   = 3,
                 n.burnin   = 4000,
                 n.iter     = 5000)



