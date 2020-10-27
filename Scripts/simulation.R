library(R2jags)
library(nuwcru)

# concept ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# maybe a simplification, but my high level understanding of what you want is the following:
#     - model IVI as a function of chick age, and brood size, incorporate grouping variables as needed (nest ID etc.)
#     - Variation in IVI may depend on the year, causal mechanisms will be investigated post-hoc.
#           - investigate this yearly variation by allowing variance to be estimated seperately for f(year) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Model

# ivi_i = Normal(mu_i, sigma^2*w) - IVI are normally distributed with mean = mu and variation = sigma^2 multiplied by a 
#                                    weighting matrix corresponding to year
# mu_i = a_nest + b_1 * chickage_i + b_2 * broodsize        - expected values are determined by linear function described by population level mean, 
#                                                                  random intercept for each nest, and chich age / brood size
# a = Normal(a-bar, sigma_a)
# a-bar = Normal(0, 1.5)
# sigma_a = Exponential(1)
# B_chickage = Normal(0, 1) 
# B_broodsize = Normal(0, 1) 




# Simulate data  -----------------------------------------------------------

n_years <- 7     # number of years
n_nests <- 15    # number of nests per year
n_sample <- 12   # days of data per nest
N <- n_sample * n_nests * n_years


chickage         <- rep(1:12, n_nests*n_years)
broodsize        <- sample(c(1,2,3,4), replace = TRUE, size = n_nests * n_years) 
broodsize        <- rep(broodsize, each = n_sample)
eps              <- data.frame(mean = c(1, 3, -2, -3, 5, 1),
                               var  = c(1, 2, -.5, 4, 6, 3))


d <- data.frame(nestID    = rep(as.character(1:n_nests), each = n_sample * n_years),
                b0_nest   = rep(rnorm(n_nests, 0, 2), each = n_sample * n_years),
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

  

