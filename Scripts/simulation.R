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

m <- lme(ivi ~ chickage + broodsize + year,   # you can change these, but data was created with this model
          random = ~ 1|nestID,                # random intercepts for nest
          weights = varIdent(form= ~ 1|year), # variance dependant on year
          data = d)   


# examine fixed effects
# compare these with the beta values we used to create the model
m$coefficients$fixed

# true values (what we used to create the data)
b1       #chickage
b2       #broodsize     
intercept 

# examine random intercepts
# true values
nest_effects
m$coefficients$random







### don't go past here, this is for my own benefit and not ready yet.
# we may want to use a bayesian framework for the paper, but the above model is sufficient 
# to get you going and understanding the process


# JAGS --------------------------------------------------------------------

# I want to better understand this and hard code the covariance so i understand it better
 # use variance matrix from:
# http://www.flutterbys.com.au/stats/tut/tut8.2b.html

X = model.matrix(~ 1 + chickage + broodsize + year, data = d)
year = as.factor(2013:2019)
year_mat = model.matrix(~ year - 1)

jags_data <- list(y = d$ivi,                 # ivi
                  years = 2012:2019, # year identifier for variance
                  n_years = n_years,         # number of years
                  nest = d$nestID,           # random intercept for nest
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
        }
    
    # Priors ~~~~~~~~~~~~~~~
    
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
        
        for (i in 1:n_years){
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


#Step 6: Run JAGS
params <- c("beta", "sigma", "sigma_nest")

het_m1   <- jags(data      = jags_data,
                inits      = NULL,
                parameters = params,
                model      = "het_var.txt",
                n.thin     = 10, 
                n.chains   = 3,
                n.burnin   = 4000,
                n.iter     = 5000)

