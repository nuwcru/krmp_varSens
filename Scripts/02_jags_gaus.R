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
library(tidybayes)

getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


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

d <- read_csv("/Volumes/GoogleDrive/My Drive/NuWCRU/Analysis/emhedlin/pefa.surv/data/sites.csv") %>% 
  select(site = SiteID, lat = SiteLatitudeDD , long = SiteLongitudeDD) %>%
  mutate(site = as.factor(site)) %>% 
  right_join(d, by = "site")


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
                  nest        = d$site,
                  n_nests     = n_nests,
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
    
    # Prior for random intercept at the nest level, not inducing correlation
        for (i in 1:n_nests) {a_nest[i] ~ dnorm(a_bar_nest, a_sigma_nest)}
        a_bar_nest ~ dnorm(0, 1.5)
        a_sigma_nest ~ dexp(1)
    
      # Likelihood ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         for (i in 1:N) {
         mu_obs[i] <- a_nest[nest[i]] + a_year[years[i]] + b_chickage[years[i]] * chickage[i] + b_chicks[years[i]] * chicks[i]
         y[i] ~ dnorm(mu_obs[i], tau_res)
         }
      # Residuals    
        for (i in 1:N) {res[i] <- y[i] - mu_obs[i]}

      # homogenous residuals - het residuals defined on line 132
       sigma_res ~ dunif(0, 100) # Residual standard deviation
       tau_res   <- 1 / (sigma_res * sigma_res)
    
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

# estimate starting values with lme4
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
            "res",
            #"sigma_resid",
            "rho", 
            "rho_int_chickage", 
            "rho_int_chicks", 
            "rho_chickage_chicks", 
            "a_bar_nest",
            "a_sigma_nest",
            "a_nest",
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
het_m3 <- update(het_m2, n.iter = 10000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 20000, n.thin = 10) 
het_m5 <- update(het_m4, n.iter = 20000, n.thin = 10) 






# Model Fitting Diagnostics -----------------------------------------------


# * Trace Plots -------------------------------------------------------------

## Parameters ##
# sigma_resid - yearly residual variance


# rho_int_chickage      - correlation between yearly intercept and beta_chickage
# rho_int_chicks        - correlation between yearly intercept and chicks
# rho_chickage_chicks   - correlation between beta_chicks and beta_chickage

het_mcmc <- as.mcmc(het_m2) #change back to 4 later 


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



# * Rhat --------------------------------------------------------------------

?rhat
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


# N Eff Size --------------------------------------------------------------

# to fill in


# Model Inference ---------------------------------------------------------


# * Residuals -------------------------------------------------------------

resids <- het_mcmc %>%
  spread_draws(res[i]) %>%
  group_by(i) %>%
  summarize(mode = getmode(res), sd = sd(res)) %>%
  mutate(upper95 = mode + (1.96 * sd), lower95 = mode - (1.96 * sd),
         upper80 = mode + (1.282 * sd), lower80 = mode - (1.282 * sd),
         upper50 = mode + (0.674 * sd), lower50 = mode - (0.674 * sd))


nd <- cbind(d, resids)

nd %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = blue5) +
  geom_segment(aes(x = logIVI, xend = logIVI, y = lower95, yend = upper95), colour = "#DFEBF7", size = 1.3) +
  geom_segment(aes(x = logIVI, xend = logIVI, y = lower80, yend = upper80), colour = "#A5CADF", size = 1.3) +
  geom_segment(aes(x = logIVI, xend = logIVI, y = lower50, yend = upper50), colour = "#4A84BD", size = 1.3) +
  geom_point(aes(x = logIVI, y = mode), shape = 21, colour = blue1, size = 1) +
  xlab("y_i") + ylab("Residual") +
  facet_grid(year~.) +
  theme_nuwcru() + facet_nuwcru()

# * Ranef Covariance --------------------------------------------------------
het_mcmc %>%
  window(thin=10) %>% 
  
  # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
  tidybayes::spread_draws(rho_int_chickage, rho_int_chicks) %>%

  ggplot() +
  geom_density(aes(x = rho_int_chickage),
               color = "transparent", fill = blue2, alpha = 5/10) +
  geom_density(aes(x = rho_int_chicks),
               color = "transparent", fill = red2, alpha = 5/10) +
  annotate(geom = "text", x = 0.6, y = 1.75, 
           label = "a_year | B_chickage", color = blue2, family = "Courier") +
  annotate(geom = "text", x = 0.6, y = 2, 
           label = "a_year | B_chicks", color = red2, family = "Courier") +
  scale_y_continuous(limits = c(0, 2.5)) +
  labs(subtitle = "Posterior distributions for correlations\nbetween intercepts and slopes",
       x = "correlation") +
  theme_nuwcru()


# * Parm Estimates ----------------------------------------------------------




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





# Sigmas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sigma as a ratio to 2013 (reference)
# dimension exported: 1200 x 804
het_mcmc %>%
  tidybayes::spread_draws(sigma_resid[i]) %>%
  group_by(i) %>%
  summarize(mode = getmode(sigma_resid), sd = sd(sigma_resid)) %>%
  mutate(upper95 = mode + (1.96 * sd), lower95 = mode - (1.96 * sd),
         upper80 = mode + (1.282 * sd), lower80 = mode - (1.282 * sd),
         upper50 = mode + (0.674 * sd), lower50 = mode - (0.674 * sd)) %>%
  ggplot() +
  geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
  geom_vline(xintercept = 1, colour = grey6, linetype = "dashed") +
  geom_segment(aes(x = lower95, xend = upper95, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
  geom_segment(aes(x = lower80, xend = upper80, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
  geom_segment(aes(x = lower50, xend = upper50, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
  geom_point(aes(y = 2013:2019, x = mode), shape = "|",  colour = "white", size = 5) +
  ylab("") + xlab("Yearly Sigma (residual variance)") +
  scale_x_continuous(limits = c(0.75,1.43)) +
  scale_y_continuous(breaks = 2013:2019) +
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


# Betas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sigma as a ratio to 2013 (reference)
# dimension exported: 1200 x 804
het_mcmc %>%
  tidybayes::spread_draws(b_chickage[i], b_chicks[i]) %>%
  group_by(i) %>%
  summarize(mode_chickage = getmode(b_chickage), sd_chickage = sd(b_chickage),
            mode_chicks = getmode(b_chicks), sd_chicks = sd(b_chicks)) %>%
  mutate(upper95_chickage = mode_chickage + (1.96 * sd_chickage),  lower95_chickage = mode_chickage - (1.96 * sd_chickage),
         upper80_chickage = mode_chickage + (1.282 * sd_chickage), lower80_chickage = mode_chickage - (1.282 * sd_chickage),
         upper50_chickage = mode_chickage + (0.674 * sd_chickage), lower50_chickage = mode_chickage - (0.674 * sd_chickage)) %>%
  mutate(upper95_chicks = mode_chicks + (1.96 * sd_chicks),  lower95_chicks = mode_chicks - (1.96 * sd_chicks),
         upper80_chicks = mode_chicks + (1.282 * sd_chicks), lower80_chicks = mode_chicks - (1.282 * sd_chicks),
         upper50_chicks = mode_chicks + (0.674 * sd_chicks), lower50_chicks = mode_chicks - (0.674 * sd_chicks)) %>%

  ggplot() +
  geom_hline(yintercept = 2013:2019, colour = grey8, alpha = 0.5) +
  geom_vline(xintercept = 0, colour = grey6, linetype = "dashed") +
  
  # Plot chickage betas
  geom_segment(aes(x = lower95_chickage, xend = upper95_chickage, y = 2013:2019, yend = 2013:2019), colour = "#DFEBF7", size = 2) +
  geom_segment(aes(x = lower80_chickage, xend = upper80_chickage, y = 2013:2019, yend = 2013:2019), colour = "#A5CADF", size = 2) +
  geom_segment(aes(x = lower50_chickage, xend = upper50_chickage, y = 2013:2019, yend = 2013:2019), colour = "#4A84BD", size = 2) +
  geom_point(aes(y = 2013:2019, x = mode_chickage), shape = "|",  colour = "white", size = 5) +
  
  # Plot Chicks betas
  #geom_segment(aes(x = lower95_chicks, xend = upper95_chicks, y = 2013.2:2019.2, yend = 2013.2:2019.2), colour = "#e2b6b6", size = 2) +
  #geom_segment(aes(x = lower80_chicks, xend = upper80_chicks, y = 2013.2:2019.2, yend = 2013.2:2019.2), colour = "#c76f6f", size = 2) +
  #geom_segment(aes(x = lower50_chicks, xend = upper50_chicks, y = 2013.2:2019.2, yend = 2013.2:2019.2), colour = "#7f1111", size = 2) +
  #geom_point(aes(y = 2013.2:2019.2, x = mode_chicks), shape = "|",  colour = "white", size = 5) +
  
  
  ylab("") + xlab("
                  chickage") +
  #scale_x_continuous(limits = c(0.75,1.43)) +
  scale_y_continuous(breaks = 2013:2019) +
  theme_nuwcru() + 
  theme(panel.border = element_blank(),
        #axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())





