library(R2jags)
library(nuwcru)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)


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

d <- d %>% filter(!is.na(chicks) & !is.na(chickage) & !is.na(ivi))

# change to factors, and calculate unique levels of randoms

# lengths to use for jags
n_nests <- length(levels(d$site))
n_yearsites <- length(levels(d$yearsite_f))
n_years <- length(levels(d$year))



# model matrix used to store betas

X <- model.matrix(~ 1 + chickage:year + chicks:year, data = d)


jags_data <- list(y           = d$ivi,    # ivi 
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
sink("het_var_gamma.txt")
cat("
    model{
    
    # Likelihood ###################################
    
        for (i in 1:N) {
          y[i]  ~ dgamma(r[years[i]], mu.eff[i])
          mu.eff[i]  <- r[years[i]] / mu[i] 
          log(mu[i]) <- inprod(beta[], X[i,]) + g[yearsite[i]] + a[nest[i]]
        }
    
    # Priors ######################################
        
     # prior for residual variance weighting matrix ###
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        r[i]     <- pow(sigma[i], -2)
        }

     # priors for betas ###
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
    
     # prior for random intercepts, a = site, g = yearsite ###
       # for (i in 1:n_years) {year_int[i] ~ dnorm(year_bar, sigma_year)}
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
        for (i in 1:n_yearsites) {g[i] ~ dnorm(g_bar, sigma_yearsite)}
    
     # prior for mean/variance of random intercepts ###
        a_bar ~ dnorm(0, 1.5)
        sigma_nest ~ dexp(1)
        g_bar ~ dnorm(0, 1.5)
        sigma_yearsite ~ dexp(1)
        
    }
", fill = TRUE)
sink()



# Store draw information from the folowing parms
params <- c("beta", "r", "mu", "a", "g")

het_m1   <- jags(data      = jags_data,
                 inits      = NULL,     # runs ok w/o starting values, change if things become more complex
                 parameters = params,
                 model      = "het_var_gamma.txt",
                 n.thin     = 10, 
                 n.chains   = 3,
                 n.burnin   = 4000,
                 n.iter     = 5000)

het_m2 <- update(het_m1, n.iter = 10000, n.thin = 10) 
het_m3 <- update(het_m2, n.iter = 20000, n.thin = 10) 
het_m4 <- update(het_m3, n.iter = 20000, n.thin = 10) 


het_mcmc <- as.mcmc(het_m4)#change back to 4 later 


# for some reason ggplot isn't working for me... 
# mixing
het_mcmc %>%
    window(thin=10) %>% 
    
    # sigma = the yearly residual variance, a = random intercept for site, g = random intercept for yearsite
    tidybayes::gather_draws(beta[i], r[i]) %>%
    
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


devtools::install_github("brooke-watson/BRRR")
library(BRRR)
skrrrahh(15)
library(brms)



# brms --------------------------------------------------------------------



# sample from the prior ---------------------------------------------------


m1_prior<- brm(ivi ~ 1 + (1 + chicks + chickage | year ),
            data = d, 
            family = lognormal(),
            prior = c(set_prior("normal(0,.001)", class = "b"),
                      set_prior("cauchy(0,.001)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 2000, 
            chains = 4,
            sample_prior = "only", 
            control = list(adapt_delta = 0.95))

get_prior(ivi ~ 1 + chicks + chickage + (1|year),
          data = d, family = lognormal())
# Plot prior
draws_prior <- d %>%
    tidyr::expand(chickage = 1:12, chicks = 1:4, year = 2013:2019) %>%
    tidybayes::add_fitted_draws(m1_prior, n = 100)


plot(density(draws_prior$.value))

draws_prior %>% group_by(chicks) %>% summarize(mean = mean(.value, na.rm = TRUE))

p1 <- ggplot(draws_prior) +
    aes(x = chickage, y = .value) +
    geom_line(aes(group = .draw), alpha = .2) +
    theme(
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()
    ) + 
    expand_limits(y = 0:1) +
    ggtitle("Plausible curves before seeing data")
p1


library(rnaturalearth)
install.packages("rnaturalearth")
library("rnaturalearthdata")

# fit the model -----------------------------------------------------------

d <-
    d %>% 
    mutate(lat_adj  = lat  * 0.11132,
           lon2_adj = lon2 * 0.11132)







fit1 <- brm(ivi ~ 1 + (1 + chicks + chickage | year ),
            sigma ~ year,
            data = d, 
            family = lognormal(),
            prior = c(#set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      prior("cauchy(0, 2)", class = "sigma"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 2000, 
            chains = 4,
            control = list(adapt_delta = 0.98))

fit2 <- brm(ivi ~ 1 + chickage + (1 + chicks + chickage | year ),
            data = d, 
            family = lognormal(),
            prior = c(set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      #prior("cauchy(0, 2)", class = "sigma"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, 
            iter = 2000, 
            chains = 4,
            control = list(adapt_delta = 0.98))
waic(fit1)
# examine residuals ~~~~~~~~~~~~~~~~~~~~~~~~~

resid <- residuals(fit1,
             method = "posterior_predict",
             type = "pearson")


nd <- d %>% filter(!is.na(chicks)) %>%filter(!is.na(ivi))
nd <- cbind(nd, resid)

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





# trace plots
post <- posterior_samples(fit1, add_chain = T)

post %>% 
    select(-lp__) %>% 
    gather(key, value, -chain, -iter) %>% 
    mutate(chain = as.character(chain)) %>% 
    
    ggplot(aes(x = iter, y = value, group = chain, color = chain)) +
    geom_line(size = 1/15) +
    scale_color_manual(values = c("#80A0C7", "#B1934A", "#A65141", "#EEDA9D")) +
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




# correlations among randoms
post <- posterior_samples(fit1)

post %>%
    ggplot() +
    geom_density(aes(x = cor_year__chicks__chickage),
                 color = "transparent", fill = "#A65141", alpha = 3/10) +
   geom_density(aes(x = cor_year__Intercept__chickage),
             color = "transparent", fill = "#A65141", alpha = 3/10) +
    geom_density(aes(x = cor_year__Intercept__chicks),
                 color = "transparent", fill = "#A65141", alpha = 3/10) +
    annotate(geom = "text", x = -0.2, y = 1.1, 
             label = "a_year | B_chickage", color = "#A65141", family = "Courier") +
    annotate(geom = "text", x = 0.75, y = 1.5, 
             label = "B_chicks | B_age", color = "#A65141", family = "Courier") +
    annotate(geom = "text", x = -.7, y = 1.5, 
             label = "a_year | B_chicks", color = "#A65141", family = "Courier") +
    scale_y_continuous(limits = c(0, 2.5)) +
    labs(subtitle = "Posterior distributions for correlations\nbetween intercepts and slopes",
         x = "correlation") +
    theme_nuwcru()






m1 <-
    coef(fit1)$year[ , 1, 1:3] %>%
    as_tibble() %>%               
    mutate(year = 2013:2019) %>%  
    select(year, everything())    

draws_posterior <- d %>%
    tidyr::expand(chickage = 1:12, chicks = 1:4, year = 2013:2019) %>%
    tidybayes::add_fitted_draws(fit1, n = 100)

library(nuwcru)
ggplot(draws_posterior) +
    aes(x = chickage, y = .value) +
    geom_point(
        aes(y = ivi), 
        color = red2, size = 2, 
        data = d, alpha = 0.2
    ) +
    geom_point(aes(group = .draw), alpha = .2) +
    facet_grid(~year) +
    theme_nuwcru() +
    scale_y_continuous(limits = c(0, 1440)) +
    expand_limits(y = 0:1) +
    ggtitle("Plausible curves after seeing data")


unique(d$site)

d %>%
    tidyr::expand(chickage = 1:12, chicks = 1:4, year = 2013:2019) %>%
    tidybayes::add_fitted_draws(fit1, n = 100) %>% 
    group_by(year, chickage) %>% 
    summarize(mean = mean(.value), sd = sd(.value)) %>%
    mutate(UL95 = mean + 1.96*sd,
           LL95 = mean - 1.96*sd,
           UL50 = mean + 0.674*sd,
           LL50 = mean - 0.674*sd) %>%
    ggplot() +
    #geom_line(aes(x = chickage, y = ivi),  colour = grey[7]) +
    geom_point(data = d, aes(x = chickage, y = ivi), shape = 21, colour = grey7 ) +
    geom_point(data = filter(d, site == "19"), aes(x = chickage, y = ivi), colour = grey3 ) +
    geom_ribbon(aes(x = chickage, ymin = LL95, ymax = UL95), fill = red3, alpha = 0.3) +
    geom_ribbon(aes(x = chickage, ymin = LL50, ymax = UL50), fill = red3, alpha = 0.8) +
    geom_line(aes(x = chickage, y = mean, group = year), colour = red2, linetype = "dashed") +
    geom_line(aes(x = chickage, y = mean, group = year), colour = red2) +
    xlab("") + ylab("first snowfree day (julian)") +
    scale_y_continuous(limits = c(0,1440)) +
    facet_wrap(~year, nrow = 1) +
    theme_nuwcru()       
    
residuals(fit1)

p3