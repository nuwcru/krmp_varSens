# We take jags_gaus.R another step. To investigate correlation between
# slope and intercepts in year, we want to nest random slopes
# for chicks + chickage in random intercept for year
# https://people.ucsc.edu/~abrsvn/general-random-effects-jags.html


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


# Jags model / run --------------------------------------------------------


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
    
    
      # Likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         for (i in 1:N) {
         mu_obs[i] <- a_year[years[i]] + b_chickage[years[i]] * chickage[i] + b_chicks[years[i]] * chicks[i]
         y[i] ~ dnorm(mu_obs[i], tau[years[i]])
         }

      # homogenous residuals - het residuals defined on line 132
       # sigma_res ~ dunif(0, 100) # Residual standard deviation
       # tau_res   <- 1 / (sigma_res * sigma_res)
    
     # het residuals - prior for residual variance weighting matrix
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma_resid[i] <- z[i] / sqrt(chSq[i])
        tau[i]   <- pow(sigma_resid[i], -2)
        }
    
    # Priors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        for (i in 1:3) {
            xi_prior[i] ~ dunif(0, 100)
            mu_raw_prior[i] ~ dnorm(0, 0.0001)
            mu_prior[i] <- xi_prior[i] * mu_raw_prior[i]
        }
        
        mu_int_prior    <- mu_prior[1]
        mu_chickage_prior <- mu_prior[2]
        mu_chicks_prior <- mu_prior[3]
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
            "b_chickage", 
            "b_chicks", 
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

het_m2 <- update(het_m1, n.iter = 10000, n.thin = 10) 
het_m3 <- update(het_m2, n.iter = 20000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 20000, n.thin = 10) 






# Trace Plots -------------------------------------------------------------

## Parameters ##
# sigma_resid - yearly residual variance

# off diagonal covariances between ranefs
# rho[1,3]
# rho[1,2]
# rho[2,1]
# rho[2,3]
# rho[3,1]
# rho[3,2]

# rho_int_chickage
# rho_int_chicks
# rho_chickage_chicks

het_mcmc <- as.mcmc(het_m3) #change back to 4 later 


het_mcmc %>%
  window(thin=10) %>% 
  
  # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
  tidybayes::gather_draws(b_chickage[i], b_chicks[i], rho_int_chickage, rho_int_chicks) %>%
  
  # change this filter to look at mixing for the parameter of interest. See above line for options
  filter(.variable == "rho_int_chicks") %>% 
  ungroup() %>%
  mutate(term = ifelse(is.na(i), .variable, paste0(.variable,"[",i,"]"))) %>%
  ggplot(aes(x=.iteration, y=.value, color=as.factor(.chain))) +
  scale_color_manual(values=c("#461220", "#b23a48", "#fcb9b2")) +
  geom_line(alpha=0.5) +
  facet_grid(term~., scale="free_y") +
  labs(color="chain", x="iteration") +
  theme_nuwcru()



# Rhat --------------------------------------------------------------------


rhat(het_mcmc) %>% 
  data.frame() %>% 
  rownames_to_column() %>% 
  set_names("parameter", "rhat") %>% 
  filter(parameter != "lp__") %>% 
  
  ggplot(aes(x = rhat, y = reorder(parameter, rhat))) + 
  geom_segment(aes(xend = 1, yend = parameter),
               color = "#EEDA9D") +
  geom_point(aes(color = rhat > 1), 
             size = 2) +
  scale_color_manual(values = c("#80A0C7", "#A65141")) +
  labs(x = NULL, y = NULL) +
  theme_nuwcru() +
  theme(legend.position = "none",
        axis.ticks.y    = element_blank(),
        axis.text.y     = element_text(hjust = 0))



# Ranef Covariance --------------------------------------------------------


posterior_samples(het_mcmc) %>%
  ggplot() +
  geom_density(aes(x = cor_year__Intercept__chickage),
               color = "transparent", fill = blue2, alpha = 5/10) +
  geom_density(aes(x = cor_year__Intercept__chicks),
               color = "transparent", fill = red2, alpha = 5/10) +
  annotate(geom = "text", x = 0.6, y = 1.75, 
           label = "a_year | B_chickage", color = blue2, family = "Courier") +
  annotate(geom = "text", x = 0.6, y = 2, 
           label = "a_year | B_chicks", color = red2, family = "Courier") +
  scale_y_continuous(limits = c(0, 2.5)) +
  labs(subtitle = "Posterior distributions for correlations\nbetween intercepts and slopes",
       x = "correlation") +
  theme_nuwcru()


# Parm Estimates ----------------------------------------------------------




# Convert to data to be exported if needed ~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








# Posterior Prediction ----------------------------------------------------

getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}



m1 <-
  coef(fit1)$year[ , 1, 1:3] %>%
  as_tibble() %>%               
  mutate(year = 2013:2019) %>%  
  select(year, everything())    


d %>%
  tidyr::expand(chickage = 1:12, chicks = 1:4, year = 2013:2019) %>%
  tidybayes::add_predicted_draws(fit1, n = 100) %>%
  group_by(chickage, year) %>%
  summarize(mode = mean(.prediction),
            sd = sd(.prediction)) %>%
  mutate(upper95 = mode + (sd*1.96),
         upper80 = mode + (sd*1.282),
         upper50 = mode + (sd*0.674),
         lower95 = mode - (sd*1.96),
         lower80 = mode - (sd*1.282),
         lower50 = mode - (sd*0.674)) %>%
  ggplot() +
  geom_ribbon(aes(x = chickage, ymin = lower95, ymax = upper95), fill = red5, alpha = 0.4) +
  geom_ribbon(aes(x = chickage, ymin = lower80, ymax = upper80), fill = red4, alpha = 0.5) +
  geom_ribbon(aes(x = chickage, ymin = lower50, ymax = upper50), fill = red3, alpha = 0.8) +
  geom_line(aes(x = chickage, y = mode), colour = red1) +
  facet_grid(~year) +
  theme_nuwcru() 





