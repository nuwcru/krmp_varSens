
# This script changes the model to brms, a package that converts lme4-like syntax to stan.
# fully bayesian, but unforunately we can't model heterogenous residuals, so if that's needed, we'll have to 
# go back to jags (or Stan). But it's quick and painless to compare models using brms, and maybe easier for you 
# to understand.


install.packages("tictoc")
install.packages("commandstanr")
library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(brms)
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

# add lat long locations
d <- read_csv("/Volumes/GoogleDrive/My Drive/NuWCRU/Analysis/emhedlin/pefa.surv/data/sites.csv") %>% 
  select(site = SiteID, lat = SiteLatitudeDD , long = SiteLongitudeDD) %>%
  mutate(site = as.factor(site)) %>% 
  right_join(d, by = "site")


# Models ------------------------------------------------------------------


# chicks + chickage as fixed, but we also will allow their slopes to vary for each year. 
# We'll also allow correlation between the intercepts and slopes as per kim.
# using the log transformed to keep things simple

tictoc::tic()
fit1 <- brm(logIVI ~ 1 +  (1 + chicks + chickage | year ),
            data = d, 
            # you can delete this prior section, brms will automatically pick some uninformative priors for you
            prior = c(#set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.98),
            backend = "cmdstanr",
            cores = 8)
tictoc::toc()

fit2 <- brm(logIVI ~ 1 +  chickage + (1 + chicks + chickage | year ),
            data = d, 
            # you can delete this prior section, brms will automatically pick some uninformative priors for you
            prior = c(#set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 2000, 
            chains = 4,
            control = list(adapt_delta = 0.98),
            cores = 8)

fit2 <- brm(logIVI ~ 1 + (1 + chicks + chickage | year ),
            data = d, 
            # you can delete this prior section, brms will automatically pick some uninformative priors for you
            prior = c(#set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.98),
            cores = 8)

waic(fit2)

# examine residuals ~~~~~~~~~~~~~~~~~~~~~~~~~
# work in progress

x <- d %>%
  add_residual_draws(fit1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  facet_grid(~year) +
  stat_pointinterval() +
  theme_nuwcru() + 
  facet_nuwcru()
    







loc <- read_csv("/Volumes/GoogleDrive/My Drive/NuWCRU/Analysis/emhedlin/pefa.surv/data/sites.csv") %>% 
  select(site = SiteID, lat = SiteLatitudeDD , long = SiteLongitudeDD) %>%
  mutate(site = as.factor(site))

nd <- nd %>% left_join(loc, by = "site")

# checking for spatial patterns in resids
nd %>% group_by(site, lat, long) %>% summarize(mean_resid = mean(abs(Estimate))) %>%
  ggplot() +
  geom_point(aes(x = long, y = lat, size = mean_resid), shape = 21) +
  # facet_grid(~year) +
  theme_nuwcru() + facet_nuwcru()





# trace plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

posterior_samples(fit1, add_chain = T) %>% 
  select(-lp__) %>% 
  gather(key, value, -chain, -iter) %>% 
  mutate(chain = as.character(chain)) %>% 
  
  ggplot(aes(x = iter, y = value, group = chain, color = chain)) +
  geom_line(size = 1/15) +
  scale_color_manual(values = c(red1, red2, red3, red4)) +
  scale_x_continuous(NULL, breaks = c(1001, 5000)) +
  ylab(NULL) +
  theme_nuwcru() +
  theme(legend.position  = c(.825, .06),
        legend.direction = "horizontal") +
  facet_wrap(~key, ncol = 3, scales = "free_y")



# Rhat

rhat(fit1) %>% 
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





# Covariance between intercept and slopes ---------------------------------


posterior_samples(fit1) %>%
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



ranef(fit1)


# Posterior prediction ----------------------------------------------------
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

