
# This script changes the model to brms, a package that converts lme4-like syntax to stan.
# fully bayesian, but unforunately we can't model heterogenous residuals, so if that's needed, we'll have to 
# go back to jags (or Stan). But it's quick and painless to compare models using brms, and maybe easier for you 
# to understand.


# install.packages("tictoc")
# install_cmdstan()
library(cmdstanr)
# install_cmdstan()
set_cmdstan_path()

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
d <- read_csv("Data/sites.csv") %>% 
  select(site = SiteID, lat = SiteLatitudeDD , long = SiteLongitudeDD) %>%
  mutate(site = as.factor(site)) %>% 
  right_join(d, by = "site")


# Models ------------------------------------------------------------------


# chicks + chickage as fixed, but we also will allow their slopes to vary for each year. 
# We'll also allow correlation between the intercepts and slopes as per kim.
# using the log transformed to keep things simple

# 2 models, both identical except m2 will allow heterogenous residual variance by year



# * Priors ----------------------------------------------------------------
# Sample from the prior only do get a sense of potential likelihood estimations




# Homogenous variance by year
m1_cor_prior <- brm(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ),   # linear model with random slopes for chicks and chickage nested in year
            data = d, 
            prior = c(set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(5, 1)", class = "Intercept"),
                      set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                      set_prior("cauchy(0,1)", class = "sd")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.99),
            backend = "cmdstanr", 
            cores = 8)



# Heterogenous variance by year
 m2_cor_prior <- brm(bf(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ), sigma ~ year),   # linear model with random slopes for chicks and chickage nested in year
            data = d, 
            prior = c(set_prior("normal(0, 1)", class = "b"),
                      set_prior("normal(5, 1)", class = "Intercept"),
                      set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                      set_prior("cauchy(0,1)", class = "sd")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.99),
            backend = "cmdstanr", 
            cores = 8)

 
 


# * Prior predictive check ------------------------------------------------

## Chicks variation
# possible values given the prior specifications - 
d %>%
  tidyr::expand(chickage = 5, chicks = 0:5, year = factor(2013:2019)) %>%
  tidybayes::add_fitted_draws(m1_cor_prior, n = 100) %>% 
  arrange(year, .draw) %>%
  group_by(.draw, year) %>%
  ggplot() +
  geom_line(aes(x = chicks, y = .value, group = .draw), alpha = .2) +
  facet_grid(~year) +
  scale_y_continuous(limits = c(-20, 20)) +
  ggtitle("Plausible curves before seeing data") +
  theme_nuwcru()

## chickage variation
# possible values given the prior specifications - 
d %>%
  tidyr::expand(chickage = 0:12, chicks = 2, year = factor(2013:2019)) %>%
  tidybayes::add_fitted_draws(m1_cor_prior, n = 100) %>% 
  arrange(year, .draw) %>%
  group_by(.draw, year) %>%
  ggplot() +
  geom_line(aes(x = chickage, y = .value, group = .draw), alpha = .2) +
  facet_grid(~year) +
  scale_y_continuous(limits = c(-500, 500)) +
  ggtitle("Plausible curves before seeing data") +
  theme_nuwcru()





# * Fit Model -------------------------------------------------------------


m1_cor <- brm(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ),   # linear model with random slopes for chicks and chickage nested in year
            data = d, 
            prior = c(set_prior("normal(-0.1, 1)", class = "b"),
                      set_prior("normal(5, 2)", class = "Intercept"),
                      set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                      set_prior("cauchy(0,2)", class = "sd")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.99),
            backend = "cmdstanr", 
            cores = 8)


m2_cor <- brm(bf(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ), sigma ~ year),   # linear model with random slopes for chicks and chickage nested in year
                    data = d, 
                    prior = c(set_prior("normal(0, 1)", class = "b"),
                              set_prior("normal(5, 1)", class = "Intercept"),
                              set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                              set_prior("cauchy(0,1)", class = "sd")),
                    warmup = 1000, 
                    iter = 5000, 
                    chains = 4,
                    control = list(adapt_delta = 0.99),
                    # backend = "cmdstanr", # R crashing with cmdstanr
                    cores = 8)





# examine residuals ~~~~~~~~~~~~~~~~~~~~~~~~~
# work in progress

resids <- residuals(m2_cor, probs = c(0.05, 0.95))
x <- cbind(d, resids)

fitted <- fitted(m2_cor)
x$fitted <- fitted[,1]

# something very wrong about 2019. That diagnoal line indicates an issue in the data. Error perfectly scales with logIVI in some case
x %>%
  ggplot() +
  geom_point(aes(x = logIVI, y = Estimate), alpha = 0.2, colour = blue2) +
  facet_grid(~year) +
  theme_nuwcru() +
  facet_nuwcru()

# list the sites in 2019
x %>% filter(year == "2019") %>% distinct(site)

# we can see there's an issue with site 8, this is not a normal pattern
x %>% filter(year == "2019") %>%
  ggplot() +
  geom_point(aes(x = logIVI, y = Estimate), colour = grey6) +
  geom_point(data = filter(x, site == 8 & year == "2019"), aes(x = logIVI, y = Estimate), colour = blue2) +
  facet_grid(~year) +
  theme_nuwcru() +
  facet_nuwcru()

# let's fit a linear model to every site. If we get large absurdly large p values, we'll
# know there's a strange pattern associated with that site.
x %>% filter(year == "2019") %>%
  group_by(site) %>%
  do(fit_resid = broom::tidy(lm(logIVI ~ Estimate, data = .))) %>% 
  unnest(fit_resid) %>%
  filter(term == "Estimate")

# from the above, we can see sites 8 and 151 have incredible tight linear patterns
# between the logIVI and the model residuals. This isn't right. let's plot both sites
x %>% filter(year == "2019") %>%
  ggplot() +
  geom_point(aes(x = logIVI, y = Estimate), colour = grey6) +
  geom_point(data = filter(x, site == 151 & year == "2019"), aes(x = logIVI, y = Estimate), alpha = 0.5, colour = red2) +
  geom_point(data = filter(x, site == 8 & year == "2019"), aes(x = logIVI, y = Estimate), alpha = 0.5, colour = blue2) +
  facet_grid(~year) +
  theme_nuwcru() +
  facet_nuwcru()

# logIVI by date, highlight the problem sites
x %>% filter(year == "2019") %>%
  ggplot() +
  geom_point(aes(y = logIVI, x = date), colour = grey6) +
  geom_point(data = filter(x, site == 151 & year == "2019"), aes(y = logIVI, x = date), colour = red2) +
  geom_point(data = filter(x, site == 8 & year == "2019"), aes(y = logIVI, x = date), colour = blue2) +
  theme_nuwcru()

d %>%
  add_fitted_draws(m1_cor) 


# trace plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

posterior_samples(m1_cor, add_chain = T) %>% 
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

rhat(m1_cor) %>% 
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


posterior_samples(m1_cor) %>%
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






# Posterior prediction ----------------------------------------------------
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

# betas
coef(m1_cor)$year[ , 1, 1:3] %>%
  as_tibble() %>%               
  mutate(year = 2013:2019) %>%  
  select(year, everything())    

## Chickage
d %>%
  tidyr::expand(chickage = 1:12, chicks = 2, year = 2013:2019) %>%
  tidybayes::add_predicted_draws(m1_cor, n = 100) %>%
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


# Brood Size - Chicks
# Look at the slope for brood size while holding chickage constant (specified on line 288)
chicks_fit <- d %>%
  tidyr::expand(chickage = 5, chicks = 1:4, year = factor(2013:2019)) %>%
  tidybayes::add_fitted_draws(m1_cor, n = 100) 

chicks_fit %>%
ggplot() +
  geom_line(aes(x = chicks, y = .value, group = .draw), alpha = .2) +
  geom_point(data = filter(d, chickage == 5), aes(x = chicks, y = logIVI), colour = blue2) +
  facet_grid(~year) +
  scale_y_continuous(limits = c(2.5, 7.5)) +
  ggtitle("Plausible curves before seeing data") +
  theme_nuwcru()