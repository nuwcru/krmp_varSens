library(R2jags)
library(nuwcru)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ concept ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# maybe a simplification, but my high level understanding of what you want is the following:
#     - model IVI as a function of chick age, and brood size, incorporate grouping variables as needed (nest ID etc.)
#     - Variation in IVI may depend on the year, causal mechanisms will be investigated post-hoc.
#           - investigate this yearly variation by allowing variance to be estimated seperately for f(year) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Simulate data  -----------------------------------------------------------

n_years  <- 7     # number of years
n_nests  <- 15    # number of nests per year
n_sample <- 12    # days of data per nest
N <- n_sample * n_nests * n_years


chickage   <- rep(1:12, n_nests*n_years)
broodsize  <- sample(c(1,2,3,4), replace = TRUE, size = n_nests * n_years) 
broodsize  <- rep(broodsize, each = n_sample)

# random nest intercepts
nest_effects <- c(-0.51735632,
                  2.09047711,
                  0.03626926,
                  -0.94967298,
                  1.07797951,
                  2.66709311,
                  -0.51708258,
                  -2.11048553,
                  0.44702011,
                  -1.21980210,
                  0.69325885,
                  -2.42639331,
                  2.74158690,
                  0.01023289,
                  -2.02312492)


# combine data for simplicity
d <- data.frame(nestID    = rep(as.character(1:n_nests), each = n_sample * n_years),
                b0_nest   = rep(nest_effects, each = n_sample * n_years),
                chickage  = chickage,
                broodsize = broodsize,
                year      = rep(as.character(2013:2019), each = n_nests * n_sample),
                eps_mean  = rep(0, length = N),
                eps_sigma = rep(c(1,3,2,5,3,6,4), each = n_nests * n_sample))

head(d)

# Parameters to be estimated
b1       <- -1    # beta for chickage
b2       <- -2    # beta for broodsize
intercept <- 15   # population mean

# other ways to do this, but I find a loop intuitve since it's generating each observation using the mechanisms we specified
# similar to ivi ~ 1 + chickage + broodsize + (1|nest) + epsilon_year

for (i in 1:N){ 
      d$ivi[i] <- intercept + d$b0_nest[i] +                  # overall population mean IVI + random intercept for nest
                  b1 * d$chickage[i] + b2 * d$broodsize[i] +  # fixed effects
                  rnorm(1, d$eps_mean[i], d$eps_sigma[i])     # draw residuals from year specific distributions
    }

  

# visualize simulated data ------------------------------------------------

### yearly variance ~~~~~~~~~~~~~
boxplot(ivi ~ year, d)

# raw points
d %>%
  ggplot() +
  geom_jitter(aes(x = year, y = ivi)) +
  xlab("") +
  theme_nuwcru()




### nest specific variance 
boxplot(ivi ~ nest, d)

d %>%
  ggplot() +
  geom_jitter(aes(x = nestID, y = ivi)) +
  xlab("") +
  theme_nuwcru()



### covariates ~~~~~~~~~~~~~~~~~~
# chickage ~ ivi
d %>% 
  ggplot() +
  geom_jitter(aes(x = chickage, y = ivi)) +
  theme_nuwcru()

# broodsize ~ ivi 
d %>% 
  ggplot() +
  geom_jitter(aes(x = broodsize, y = ivi)) +
  theme_nuwcru()


# Model simulated data ----------------------------------------------------
library(nlme)

# regular model with no covariance structure
m1 <- lme(ivi ~ chickage + broodsize,   # you can change these, but data was created with this model
          random = ~ 1|nestID,                # random intercepts for nest
          data = d)   

# normalized residuals for the covariates within each year 
plot(m1, resid(., type = "n") ~ broodsize | year)
plot(m1, resid(., type = "n") ~ chickage | year)

# residuals clearly differ among years which isn't good
#     - let's model it again, but allow variance to differ among years



# regular model _with_ variance dependant on year
m2 <- lme(ivi ~ chickage:year + broodsize:year,   # you can change these, but data was created with this model
         random = ~ 1 | nestID,           # random intercepts for nest
         weights = varIdent(form = ~ 1|year),
         data = d)   




# same plots as before
plot(m2, resid(., type = "n") ~ broodsize | year)
plot(m2, resid(., type = "n") ~ chickage | year)

d$resid1 <- resid(m1, type="normalized") 
d$fitted1 <- as.vector(m1$fitted[,1])
d$resid2 <- resid(m2, type="normalized") 
d$fitted2 <- as.vector(m2$fitted[,2])



d %>%
  ggplot() +
  geom_jitter(aes(x = chickage, y = fitted2)) +            # model fit
  geom_jitter(aes(x = chickage, y = ivi), colour = red4) + # true data
  facet_grid(.~year) +
  theme_nuwcru()


# examine fixed effects
# estimated fixed effects
m$coefficients$fixed

# true values (what we used to create the data)
b1       #chickage
b2       #broodsize     
intercept 

# examine random intercepts
# true values
nest_effects
# estimated
m$coefficients$random

d$resid <- as.vector(m$residuals[,1])
d$fitted <- as.vector(m$fitted[,1])
glimpse(d)

d %>%
  ggplot() +
  geom_jitter(aes(x = year, y = fitted)) +
  theme_nuwcru()

### don't go past here, this is for my own benefit and not ready yet.
# we may want to use a bayesian framework for the paper, but the above model is sufficient 
# to get you going and understanding the process













# JAGS --------------------------------------------------------------------

# I want to better understand this and hard code the covariance so i get a better feel for it
 # use variance matrix from:
# http://www.flutterbys.com.au/stats/tut/tut8.2b.html

X = model.matrix(~ 1 + chickage + broodsize + year, data = d)


jags_data <- list(y = d$ivi,                 # ivi
                  years = as.factor(d$year), # year identifier for variance
                  nest = as.factor(d$nestID),# random intercept for nest
                  n_nests = n_nests,         # number of nests
                  X = X,                     # intercept + covariates
                  N = N,                     # sample size
                  K = ncol(X))               # Number of betas



#~~~~~~~~ Jags Model ~~~~~~~~~~~~#
sink("het_var.txt")
cat("
    model{
    
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,]) + a[nest[i]]
    
        y_pred[i] ~ dnorm(mu[i], tau[years[i]])
        }
    
    # Priors ~~~~~~~~~~~~~~~
    
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)
        
        sigma[i] <- z[i] / sqrt(chSq[i])
        z[i]     ~ dnorm(0, 0.04)I(0,)
        chSq[i]  ~ dgamma(0.5, 0.5)
        tau[i]   <- pow(sigma[i], -2)
        }

    
    # prior for random intercept
        for (i in 1:n_nests) {a[i] ~ dnorm(0, tau_nest)}
    
    # prior for variance of random intercept
        sigma_nest ~ dunif(0,5)
        tau_nest <- 1 / (sigma_nest * sigma_nest)

    }
", fill = TRUE)
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Store draw information from the folowing parms
params <- c("beta", "sigma", "a", "sigma_nest", "y_pred", "mu")

het_m1   <- jags(data      = jags_data,
                inits      = NULL,
                parameters = params,
                model      = "het_var.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

het_m2 <- update(het_m1, n.iter = 10000, n.thin = 10) 
het_m3 <- update(het_m2, n.iter = 50000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 50000, n.thin = 10) 

het_mcmc <- as.mcmc(het_m4)

# good mixing for betas, not sure why sigma[8] and [9] are terrible...
het_mcmc %>%
  window(thin=10) %>% 
  tidybayes::gather_draws(beta[i], sigma[i], a[i]) %>%
  filter(.variable == "beta") %>%
  ungroup() %>%
  mutate(term = ifelse(is.na(i), .variable, paste0(.variable,"[",i,"]"))) %>%
  ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
  scale_color_manual(values=c("#461220", "#b23a48", "#fcb9b2")) +
  geom_line(alpha=0.5) +
  facet_grid(term~., scale="free_y") +
  labs(color="chain", x="iteration") +
  theme_nuwcru()



sleep_lm_pred = 
  het_mcmc[1] %>%
  tidybayes::spread_draws(mu[i], y_pred[i]) %>%
  ungroup() %>%
  left_join(
    mutate(sleep, i = 1:n())
  ) %>%
  mutate(resid = Reaction - mu)
