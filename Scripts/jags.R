library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)



# Data Load/prep ----------------------------------------------------------

# filter NA's out of covariates
only_unsupp <- only_unsupp %>% filter(!is.na(chicks) & !is.na(chickage))

# change to factors, and calculate unique levels of randoms
nestID <- as.factor(only_unsupp$site)
n_nests <- length(levels(nestID))

yearsite_f <- as.factor(only_unsupp$yearsite)
n_yearsites <- length(levels(yearsite_f))

n_years <- length(levels(only_unsupp$year))

# model matrix used to store betas
X <- model.matrix(~ 1 + chickage:year + chicks:year, data = only_unsupp)

head(X)

# data for Jags model
jags_data <- list(y = only_unsupp$logIVI,     # ivi 
                  years = only_unsupp$year_f, # year identifier for variance
                  nest = nestID,              # random intercept for nest
                  yearsite = yearsite_f,      # random intercept for yearsite
                  n_years = n_years,          # number of years
                  n_nests = n_nests,          # number of nests
                  n_yearsites = n_yearsites,  # number of unique yearsites
                  X = X,                      # intercept + covariates (model matrix)
                  N = nrow(only_unsupp),      # sample size
                  K = ncol(X))                # Number of betas



#~~~~~~~~ Jags Model ~~~~~~~~~~~~#
sink("het_var.txt")
cat("
    model{
    
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,]) + a[nest[i]] + g[yearsite[i]]
        
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
params <- c("beta", "sigma", "a", "g", "sigma_nest", "y_pred", "mu")

het_m1   <- jags(data      = jags_data,
                 inits      = NULL,     # runs ok w/o starting values, change if things become more complex
                 parameters = params,
                 model      = "het_var.txt",
                 n.thin     = 10, 
                 n.chains   = 3,
                 n.burnin   = 4000,
                 n.iter     = 5000)

het_m2 <- update(het_m1, n.iter = 10000, n.thin = 10) 
het_m3 <- update(het_m2, n.iter = 50000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 50000, n.thin = 10) 


# convert our jags output into mcmc object
het_mcmc <- as.mcmc(het_m3)

# Save and load model
save(het_mcmc, file = "Models/het_mcmc_m3.rda")
het_mcmc <- load("Models/het_mcmc_m3.rda")


# model matrix as a reminder
head(X)



# Model Diagnostics -------------------------------------------------------

# mixing
het_mcmc %>%
    window(thin=10) %>% 
    
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], sigma[i], a[i], g[i]) %>%
    
    # change this filter to look at mixing for the parameter of interest. See above line for options
    filter(.variable == "g") %>% 
    ungroup() %>%
    mutate(term = ifelse(is.na(i), .variable, paste0(.variable,"[",i,"]"))) %>%
    ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
    scale_color_manual(values=c("#461220", "#b23a48", "#fcb9b2")) +
    geom_line(alpha=0.5) +
    facet_grid(term~., scale="free_y", ncol = 2) +
    labs(color="chain", x="iteration") +
    theme_nuwcru()

# ACF

# more diagnostics, to do




# Effects Plotting --------------------------------------------------------

# join data set with model predictions
het_d <- het_mcmc %>%
    tidybayes::spread_draws(mu[i], y_pred[i]) %>%
    ungroup() %>%
    left_join(
        mutate(only_unsupp, i = 1:n())) %>%
    mutate(resid = logIVI - mu)



# Predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log ivi ~ chickage
# dimension exported: 1200 x 804
ggplot(only_unsupp, aes(x=chickage, y=logIVI)) +
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
ggplot(only_unsupp, aes(x=chicks, y=logIVI)) +
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
het_mcmc[1] %>%
    tidybayes::spread_draws(sigma[i]) %>%
    group_by(i) %>%
    summarize(mean = mean(sigma), sd = sd(sigma)) %>%
    mutate(upper95 = mean + (1.96 * sd), lower95 = mean - (1.96 * sd),
           upper80 = mean + (1.282 * sd), lower80 = mean - (1.282 * sd),
           upper50 = mean + (0.674 * sd), lower50 = mean - (0.674 * sd)) %>%
    ggplot() +
    geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
    geom_vline(xintercept = 1, colour = grey6, linetype = "dashed") +
    geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
    geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
    geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
    geom_point(aes(y = 2013:2019, x = mean), shape = 21, fill = "white", colour = "black", size = 4)+
    ylab("") + xlab("Sigma yearly estimate (95% CI)") +
    scale_x_continuous(limits = c(0.7,1.3)) +
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
