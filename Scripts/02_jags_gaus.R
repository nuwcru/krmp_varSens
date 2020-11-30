# We take jags_gaus.R another step. To investigate correlation between
# slope and intercepts in year, we want to nest random slopes
# for chicks + chickage in random intercept for year


library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(bayesm)
library(lme4)


# make sure these transformations are what you want
d <- read_csv("Data/ivi_eh.csv") %>% 
  filter(supplimented == "n") %>%
  mutate(logIVI   = log(ivi),
         year     = as.factor(year),
         site     = as.factor(site),
         # chickage = chickage -1,  
         site     = as.factor(site),
         yearsite_f = as.factor(yearsite))

d <- d %>% filter(!is.na(chicks) & !is.na(chickage) & !is.na(logIVI))
dim(d)
# change to factors, and calculate unique levels of randoms

# lengths to use for jags
n_nests <- length(levels(d$site))
n_yearsites <- length(levels(d$yearsite_f))
n_years <- length(levels(d$year))


jags_data <- list(y           = d$logIVI,    # ivi 
                  chickage    = d$chickage,   # vector of chickages
                  chicks      = d$chicks,     # vector of brood sizes
                  years       = d$year,      # year identifier for variance
                  n_years     = n_years,      # number of years
                  N           = nrow(d),      # sample size
                  W           = diag(3))      # wishart matrix


#~~~~~~~~ Jags Model ~~~~~~~~~~~~#
sink("Models/het_var_cor.txt")
cat("
    model{
    
    # Covariance between ranefs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Set up the means for the multivariate ranef distribution
    
          for (i in 1:3) {
            xi[i]     ~ dunif(0, 100)
            mu_raw[i] ~ dnorm(0, .0001)
            mu[i]     <- xi[i] * mu_raw[i]
            }
        
          mu_int      <- mu[1] # mu intercept 
          mu_chickage <- mu[2] # mu slope chickage
          mu_chicks   <- mu[3] # mu chicks
          
          Tau_B_raw[1:3, 1:3] ~ dwish(W[,], 4)
          Sigma_B_raw[1:3, 1:3] <- inverse(Tau_B_raw[,])
          
          for (i in 1:3) {
              sigma[i] <- xi[i] * sqrt(Sigma_B_raw[i, i])
          }
          
          sigma_int      <- sigma[1]
          sigma_chickage <- sigma[2]
          sigma_chicks   <- sigma[3]
          
          for (i in 1:3) { for (j in 1:3) {
            rho[i, j] <- Sigma_B_raw[i, j] / sqrt(Sigma_B_raw[i, i] * Sigma_B_raw[j, j])
          }}
          
          rho_int_chickage    <- rho[1, 2] # correlation year / chickage
          rho_int_chicks      <- rho[1, 3] # correlation year / chicks
          rho_chickage_chicks <- rho[2, 3] # correlation chickage / chicks
          
          for (j in 1:n_years) {
            B_raw_hat[j, 1] <- mu_raw[1]
            B_raw_hat[j, 2] <- mu_raw[2]
            B_raw_hat[j, 3] <- mu_raw[3]
            B_raw[j, 1:3] ~ dmnorm(B_raw_hat[j, ], Tau_B_raw[, ])
            a_year[j]     <- xi[1] * B_raw[j, 1]
            b_chickage[j] <- xi[2] * B_raw[j, 2]
            b_chicks[j]   <- xi[3] * B_raw[j, 3]
        }  
    
      # homogenous residuals
        sigma_res ~ dunif(0, 100) # Residual standard deviation
        tau_res   <- 1 / (sigma_res * sigma_res)
    
    
      # Likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         for (i in 1:N) {
         mu_obs[i] <- a_year[years[i]] + b_chickage[years[i]] * chickage + b_chicks[years[i]] * chicks
         y[i] ~ dnorm(mu_obs[i], tau_res)
         }

    
    # Priors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        for (i in 1:3) {
            xi_prior[i] ~ dunif(0, 100)
            mu_raw_prior[i] ~ dnorm(0, 0.0001)
            mu_prior[i] <- xi_prior[i] * mu_raw_prior[i]
        }
        
        mu_int_prior    <- mu_prior[1]
        mu_slope1_prior <- mu_prior[2]
        mu_slope2_prior <- mu_prior[3]
        Tau_B_raw_prior[1:3, 1:3] ~ dwish(W[,], 4)
        
        Sigma_B_raw_prior[1:3, 1:3] <- inverse(Tau_B_raw_prior[,])
        for (i in 1:3) {
            sigma_prior[i] <- xi_prior[i] * sqrt(Sigma_B_raw_prior[i, i])
        }
        
        sigma_int_prior      <- sigma_prior[1]
        sigma_chickage_prior <- sigma_prior[2]
        sigma_chicks_prior   <- sigma_prior[3]
        
        for (i in 1:3) { for (j in 1:3) {
            rho_prior[i, j] <- Sigma_B_raw_prior[i, j] / sqrt(Sigma_B_raw_prior[i, i] * Sigma_B_raw_prior[j, j])
        }}
        
        rho_int_chickage_prior    <- rho_prior[1, 2]
        rho_int_chicks_prior      <- rho_prior[1, 3]
        rho_chickage_chicks_prior <- rho_prior[2, 3]

     # het residuals - prior for residual variance weighting matrix
        # for (i in 1:n_years){
        # chSq[i]  ~ dgamma(0.5, 0.5)
        # z[i]     ~ dnorm(0, 0.04)I(0,)
        # sigma_resid[i] <- z[i] / sqrt(chSq[i])
        # tau[i]   <- pow(sigma_resid[i], -2)
        # }
    }
", fill = TRUE)
sink()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lme_fit <- lme4::lmer(logIVI ~ 1 + chickage + chicks + (1 + chickage + chicks | year), data = d)
var_vec <- apply(coef(lme_fit)$year, 2, var)

# the model has become sufficiently complicated - we likely need to specify starting values
# Inits function
inits <- function() {
  list(xi = rlnorm(3), 
       mu_raw = rnorm(3), 
       Tau_B_raw = rwishart(4, diag(3)*var_vec)$W, 
       mu_raw_prior = rnorm(3), 
       Tau_B_raw_prior = rwishart(4, diag(3)*var_vec)$W)
}

# Parameters to estimate
params <- c("mu", 
            "mu_int", 
            "mu_chickage", 
            "mu_chicks", 
            "sigma", 
            "sigma_int", 
            "sigma_chickage", 
            "sigma_chicks", 
            "sigma_resid",
            "rho", 
            "rho_int_chickage", 
            "rho_int_chicks", 
            "rho_chickage_chicks", 
            "alpha", 
            "b_chickage", 
            "b_chicks", 
            "sigma_res", 
            "mu_int_prior", 
            "mu_chickage_prior", 
            "mu_chicks_prior", 
            "sigma_int_prior", 
            "sigma_chickage_prior", 
            "sigma_chicks_prior", 
            "rho_prior", 
            "rho_int_chickage_prior", 
            "rho_int_chicks_prior", 
            "rho_chickage_chicks_prior")  


het_m1   <- jags(data      = jags_data,
                 inits      = inits,
                 parameters = params,
                 model      = "Models/het_var_cor.txt",
                 n.thin     = 10, 
                 n.chains   = 3,
                 n.burnin   = 4000,
                 n.iter     = 7000)
