# here we log transform IVI so that we can easily model heterogenous residuals.
# We include random intercepts for yearsite and site, and include interactions 
# between year and chickage/chicks


library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)


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
d$chicks


# model matrix used to store betas

X <- model.matrix(~ 1 + chickage + chicks, data = d)


jags_data <- list(y           = d$logIVI,    # ivi 
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
    
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,])  + g[yearsite[i]] + a[nest[i]] 
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
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
        for (i in 1:n_yearsites) {g[i] ~ dnorm(g_bar, sigma_yearsite)}
    
     # prior for mean/variance of random intercepts
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
params <- c("beta", "sigma","g", "g_bar", "sigma_yearsite")

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
het_mcmc <- as.mcmc(het_m4)#change back to 4 later 




# Model Diagnostics -------------------------------------------------------


# for some reason ggplot isn't working for me... 
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
    theme_nuwcru()

# ACF


# Convert to data frame ---------------------------------------------------

# this section converts the data into a data frame
# Use the resulting data from in the lme4_v_jags script to visualize model.

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



arrow::write_parquet(all, "data/jags_out_newivi.parquet")



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
