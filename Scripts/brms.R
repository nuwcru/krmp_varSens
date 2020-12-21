
# This script changes the model to brms, a package that converts lme4-like syntax to stan.
# fully bayesian, but unforunately we can't model heterogenous residuals, so if that's needed, we'll have to 
# go back to jags (or Stan). But it's quick and painless to compare models using brms, and maybe easier for you 
# to understand.


# install.packages("tictoc")
# install_cmdstan()
library(cmdstanr)
install_cmdstan()
set_cmdstan_path()

library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)
library(lubridate)


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

weather <- read_csv("Data/weather data 2013-2019.csv")
weather$date <- dmy(weather$date)


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
            prior = c(set_prior("normal(0, 0.5)", class = "b", coef = "chickage"),
                      set_prior("normal(-0.1, 0.5)", class = "b", coef = "chicks"),
                      set_prior("normal(5, 1)", class = "Intercept"),
                      set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                      set_prior("cauchy(0,1)", class = "sd")),
            warmup = 1000, 
            iter = 5000, 
            chains = 4,
            control = list(adapt_delta = 0.99),
            sample_prior = "only",
            # backend = "cmdstanr", 
            cores = 8)






# * * Prior predictive check ------------------------------------------------

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
  ggtitle("Possible curves from prior information only") +
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
  scale_y_continuous(limits = c(-10, 10)) +
  ggtitle("Possible curves from prior information only") +
  theme_nuwcru()





# * Fit Model -------------------------------------------------------------

# linear model with random slopes for chicks and chickage nested in year
m1_cor <- brm(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ),   
            data = d, 
            prior = c(set_prior("normal(-0.1, 1)", class = "b"),
                      set_prior("normal(5, 2)", class = "Intercept"),
                      set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                      set_prior("cauchy(0,2)", class = "sd")),
            warmup = 1000, 
            iter = 5000, 
            chains = 3,
            cores = 3,
            control = list(adapt_delta = 0.99),
            backend = "cmdstanr", 
            threads = threading(4)) # within chain parallelization, number of threads + number of cores should match total cores you intend to use


# linear model with random slopes for chicks and chickage nested in year, heterogenous residual variance
# log(sigma_i) = beta_0 + beta_1*year 
m2_cor <- brm(bf(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year ), 
                 sigma ~ year),   
                    data = d, 
                    prior = c(set_prior("normal(0, 0.5)", class = "b", coef = "chickage"),
                              set_prior("normal(-0.1, 0.5)", class = "b", coef = "chicks"),
                              set_prior("normal(5, 1)", class = "Intercept"),
                              set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                              set_prior("cauchy(0,1)", class = "sd")),
                    warmup = 1000, 
                    iter = 5000, 
                    chains = 3,
                    cores  = 3,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr", # R crashing with cmdstanr
                    threads = threading(4))


# linear model with random slopes for chicks and chickage nested in year, random intercept for site, and resid varying by year.
m3_cor <- brm(bf(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year) + (1 | site), 
                 sigma ~ year),   
                    data = d, 
                    prior = c(set_prior("normal(0, 0.5)", class = "b", coef = "chickage"),
                              set_prior("normal(-0.1, 0.5)", class = "b", coef = "chicks"),
                              set_prior("normal(5, 1)", class = "Intercept"),
                              set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
                              set_prior("cauchy(0,1)", class = "sd")),
                    warmup = 1000, 
                    iter = 5000, 
                    chains = 3,
                    cores  = 3,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr", # R crashing with cmdstanr
                    threads = threading(4))

# m4_cor <- brm(bf(logIVI ~ 1 + chicks + chickage + (1 + chicks + chickage | year) + (1 | site), 
#                  sigma ~ site:year),   
#                     data = d, 
#                     prior = c(set_prior("normal(0, 0.5)", class = "b", coef = "chickage"),
#                               set_prior("normal(-0.1, 0.5)", class = "b", coef = "chicks"),
#                               set_prior("normal(5, 1)", class = "Intercept"),
#                               set_prior("lkj(2)", class = "cor"),    # delete this line if you don't want to estimate correlation among ranef
#                               set_prior("cauchy(0,1)", class = "sd")),
#                     warmup = 1000, 
#                     iter = 5000, 
#                     chains = 3,
#                     cores  = 3,
#                     control = list(adapt_delta = 0.99),
#                     backend = "cmdstanr", # R crashing with cmdstanr
#                     threads = threading(4))

pp_check(m1_cor)

y <- posterior_predict(m2_cor)
y_true <- d$logIVI
loo1 <- loo(m2_cor, save_psis = TRUE, cores = 4)
psis1 <- loo1$psis_object
lw <- weights(psis1)
bayesplot::ppc_loo_pit_overlay(y_true, y, lw = lw)

 # examine residuals ~~~~~~~~~~~~~~~~~~~~~~~~~
# work in progress

resids <- residuals(m1_cor, probs = c(0.05, 0.95))
res <- cbind(d, resids)

# residuals by chickage
res %>%
  ggplot() +
  geom_jitter(aes(x = chickage, y = Estimate),  alpha = 0.6, shape = 21, fill = blue2, colour = blue4) +
  facet_grid(~year) +
  xlab("") + ylab("") +
  theme_nuwcru() + facet_nuwcru()

# residuals by chicks
res %>%
  ggplot() +
  geom_jitter(aes(x = chicks, y = Estimate),  alpha = 0.6, shape = 21, fill = blue2, colour = blue4) +
  facet_grid(~year) +
  xlab("") + ylab("") +
  theme_nuwcru() + facet_nuwcru()

# look for Spatial pattern in residuals
res %>% 
  group_by(lat, long, site, year) %>% 
  summarize(scale = mean(abs(Estimate)), # the absolute error for a specific site
            dir = sum(Estimate)) %>%    # the direction of the error (neg or pos)
  ggplot() +
  geom_point(aes(x = long, y = lat, size = (1+scale)^4, fill = dir), shape = 21, colour = grey4) +
  scale_fill_distiller(palette = "Spectral") + 
  xlab("") + ylab("") +
  facet_grid(~year) +
  theme_nuwcru() + facet_nuwcru() + theme(legend.position = "bottom")


# Nothing really obvious stands out


fitted <- fitted(m1_cor)
res$fitted <- fitted[,1]

# something very wrong about 2019. That diagnoal line indicates an issue in the data. Error perfectly scales with logIVI in some case
res %>% filter(year == 2019) %>%
  ggplot() +
  geom_point(aes(x = Estimate, y = logIVI), alpha = 0.5, colour = grey6) +
  geom_point(data = filter(res, year == 2019 & site == "151"), aes(x = Estimate, y = logIVI), colour = blue3) +
  # facet_grid(~year) +
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

posterior_samples(m3_cor, add_chain = T) %>% 
  select(-lp__) %>% 
  gather(key, value, -chain, -iter) %>% 
  mutate(chain = as.character(chain)) %>%
  filter(key %in% c("b_sigma_year2014","b_sigma_year2015","b_sigma_year2016","b_sigma_year2017", "b_sigma_year2018", "b_sigma_year2019")) %>%

  ggplot(aes(x = iter, y = value, group = chain, color = chain)) +
  geom_line(size = 1/15) +
  scale_color_manual(values = c(red1, red3, red4)) +
  scale_x_continuous(NULL, breaks = c(1001, 5000)) +
  ylab(NULL) +
  theme_nuwcru() + facet_nuwcru() +
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



pp_check(m1_cor)

# Covariance between intercept and slopes ---------------------------------


posterior_samples(m2_cor) %>%
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
  tidybayes::add_predicted_draws(m3_cor, n = 100) %>%
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
  tidybayes::add_fitted_draws(m2_cor, n = 100) 

chicks_fit %>%
ggplot() +
  geom_line(aes(x = chicks, y = .value, group = .draw), alpha = .2) +
  geom_point(data = filter(d, chickage == 5), aes(x = chicks, y = logIVI), colour = blue2) +
  facet_grid(~year) +
  scale_y_continuous(limits = c(2.5, 7.5)) +
  ggtitle("Plausible curves after seeing data") +
  theme_nuwcru()



# visualize sigmas --------------------------------------------------------
install.packages("latex2exp")
library(latex2exp)
sigmas <- exp(posterior_samples(m2_cor, "^b_sigma_"))

sigmas %>% pivot_longer(cols = everything()) %>%
ggplot(aes(value)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = blue2)+
  geom_histogram(binwidth = 0.01, color = grey6, fill = grey2, alpha = 0.1) +
  facet_grid(~name) +
  xlab(TeX("$ \\sigma_{year} $")) + ylab("") +
  ggtitle(TeX("Posterior distributions of yearly $ \\sigma $")) +
  theme_nuwcru() + facet_nuwcru()

x %>% group_by(year) %>% summarize(sd = sd(Estimate))

x %>%
  ggplot() +
  geom_point(aes(x=chicks, y = logIVI), colour = grey6) +
  geom_point(aes(x=chicks, y = fitted), alpha = 0.2) +
  facet_grid(~year)





pp_check(m1_cor, type = "error_scatter_avg_vs_x", size = 1.1, alpha = 0.5,
         x = "chicks")
