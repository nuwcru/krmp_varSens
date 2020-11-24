library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(R2jags)
library(tidybayes)


d <- read_csv("Data/Clean IVI years 2013-2019.csv")


# Data Load/prep ----------------------------------------------------------

# make sure these transformations are what you want
d <- read_csv("Data/ivi_eh.csv") %>% 
    filter(supplimented == "n") %>%
    mutate(logIVI   = log(ivi),
           year     = as.factor(year),
           site     = as.factor(site),
           chickage = chickage -1,  
           site     = as.factor(site),
           yearsite_f = as.factor(yearsite))

d <- d %>% filter(!is.na(chicks) & !is.na(chickage))

# change to factors, and calculate unique levels of randoms

# lengths to use for jags
n_nests <- length(levels(d$site))
n_yearsites <- length(levels(d$yearsite_f))
n_years <- length(levels(d$year))



# model matrix used to store betas

X <- model.matrix(~ 1 + chickage:year + chicks:year, data = d)


jags_data <- list(y           = d$ivi,    # ivi 
                  years       = d$year,      # year identifier for variance
                  nest        = d$site,      # random intercept for nest
                  yearsite    = d$yearsite_f,  # random intercept for yearsite
                  n_years     = n_years,     # number of years
                  n_nests     = n_nests,     # number of nests
                  n_yearsites = n_yearsites, # number of unique yearsites
                  X           = X,           # intercept + covariates (model matrix)
                  N           = nrow(d),     # sample size
                  K           = ncol(X))     # Number of betas

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
        
     # prior for residual variance weighting matrix ###
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        r[i]     <- pow(sigma[i], -2)
        }

     # priors for betas ###
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
    
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
params <- c("beta", "r", "mu", "a", "g")

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


het_mcmc <- as.mcmc(het_m4)#change back to 4 later 


# for some reason ggplot isn't working for me... 
# mixing
het_mcmc %>%
    window(thin=10) %>% 
    
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], r[i]) %>%
    
    # change this filter to look at mixing for the parameter of interest. See above line for options
    filter(.variable == "beta") %>% 
    ungroup() %>%
    mutate(term = ifelse(is.na(i), .variable, paste0(.variable,"[",i,"]"))) %>%
    ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
    scale_color_manual(values=c("#461220", "#b23a48", "#fcb9b2")) +
    geom_line(alpha=0.5) +
    facet_grid(term~., scale="free_y") +
    labs(color="chain", x="iteration") +
    theme_nuwcru()


x <- het_mcmc %>%
    window(thin=10) %>% 
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], sigma[i], g[i])


beta <- x %>% filter(.variable == "beta") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(beta$i))){levels(beta$i)[i] <- colnames(X)[i]}
levels(beta$i)[1] <- "intercept"

sigma <- x %>% filter(.variable == "sigma") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(sigma$i))){levels(sigma$i)[i] <- paste0("sigma", "_", levels(d$year)[i])}

# a <- x %>% filter(.variable == "a") %>% mutate(i = as.factor(i))
# for (i in 1:length(unique(a$i))){levels(a$i)[i] <- unique(as.character(d$site))[i]}

g <- x %>% filter(.variable == "g") %>% mutate(i = as.factor(i))
for (i in 1:length(unique(g$i))){levels(g$i)[i] <- unique(as.character(d$yearsite))[i]}

# year <- x %>% filter(.variable == "year") %>% mutate(i = as.factor(i))
# for (i in 1:length(unique(year$i))){levels(year$i)[i] <- paste0("int", "_", levels(d$year)[i])}
# unique(year$i)

all <- rbind(beta, sigma, g)
names(all) <- c("level", "chain", "iteration", "draw", "parameter", "value")



arrow::write_parquet(all, "data/jags_out_gamma.parquet")

