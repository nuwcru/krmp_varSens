## Variance Sensitivity

#### Contents
[Modeling Real Data](#Heterogenous-residual-variance-by-year) |
[Modeling simulated data](#Simulations) |


### Intro
Investigating variance sensitivity among breeding Peregrine Falcons in response to environmental conditions and brood demand.

### File Organization

```
├── Data
│     ├── Clean IVI years 2013-2019.csv  <- data used for analysis - cleaned by Rebekah
│     ├── 
│     ├── 
│     ├── 
│     ├── 
│     └── 
├── Scripts 
│     ├──  simulation.R         <- simulate data to ensure variance is being captured properly by whatever model we use
│     │                            - nlme model for simplicity / jags for thorough understanding and benefits from bayesian framework
│     ├──  jags.R               <- Real Data using similar model structure to that tested in simulation.R
│     └──  ivi_scratch.R        <- Programmatically calculate IVI
└── Figures
      ├── beta_chickage.jpg     <- from jags.R - het resid variance by year
      ├── betas_broodsize.jpg   <- from jags.R - het resid variance by year
      ├── betas_broodsize.jpg   <- from jags.R - het resid variance by year
      ├── sigma_CIs.jpg         <- from jags.R - het resid variance by year
      ├── chickage-year.jpg     <- from jags.R - het resid variance by year
      └── broodsize-year.jpg    <- from jags.R - het resid variance by year

```
<br />

<br />

# Heterogenous residual variance by year

``` 
How do you know that your model is right? 
When the residuals contain no information.  

How do you know that your model is good enough?  
When the residuals contain no information that you can resolve.

-box and tiao
```

<br />

<br />

The code for this model is in ```scripts/jags.R```

```r
# Fixed effects
X <- model.matrix(~ 1 + chickage:year_f + chicks:year_f, data = only_unsupp)

# all data going into model below
jags_data <- list(y = only_unsupp$logIVI,     # ivi 
                  years = only_unsupp$year_f, # year identifier for variance
                  nest = nestID,              # random intercept for nest
                  yearsite = yearsite_f,      # random intercept for yearsite
                  n_years = n_years,          # number of years
                  n_nests = n_nests,          # number of nests
                  n_yearsites = n_yearsites,  # number of unique yearsites
                  X = X,                      # intercept + covariates (model matrix)
                  N = nrow(only_unsupp),      # sample size
                  K = ncol(X))                # Number of betas
```
<br />

## Model
```
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,]) + a[nest[i]] + g[yearsite[i]]
        
        # store predicted values at each iteration in y_pred
        y_pred[i] ~ dnorm(mu[i], tau[years[i]])
        }
    
    # Priors ~~~~~~~~~~~~~~~
        
     # priors for betas
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
        
     # prior for residual variance weighting matrix
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        tau[i]   <- pow(sigma[i], -2)
        }

    
     # prior for random intercepts, a = site, g = yearsite
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
        for (i in 1:n_yearsites) {g[i] ~ dnorm(g_bar, sigma_yearsite)}
    
     # prior for mean/variance of random intercepts
        a_bar ~ dnorm(0, 1.5)
        sigma_nest ~ dexp(1)
        g_bar ~ dnorm(0, 1.5)
        sigma_yearsite ~ dexp(1)
        
```

For the sake of brevity, I won't show diagnostic plots, but mixing went relatively well. Random intercepts struggled a little, so I'll have to look at that, and the intercept wasn't great. Everything else was good though. I'm using a very similar model structure to the one posted in the simulation section at the bottom of this page. That model estimated my simulated parameters values accurately.

Kim, the model output is converted to an mcmc object in the script, which I think is similar to what MCMCglmm produces. So if you have a particular workflow with model diagnostics/plotting you like, it should transfer over.

<br />

## Betas

<p align="center">
  <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/Figures/beta_chickage.jpeg">
   <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/Figures/betas_broodsize.jpeg">
</p>

<br />

## Sigma - yearly residual variance

logIVI_i ~ N(mu_i, sigma^2_year)

Below we're plotting the estimated sigma (+/- credibles) as a ratio to the reference group which is 2013 in this case. Estimates to the right of 1 (dashed line) had greater residual variance than 2013. Estimates to the left had less residual variance. Some interesting results, maybe we can explore different ways of modeling these residuals. Still model them at the yearly level, but instead of using year, generate a metric that summarizes the yearly conditions. This may be more informative and match residual variance patterns better than a binned category like year: logIVI_i ~ N(mu_i, sigma^2_yearly-conditions)

<p align="center">
  <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/Figures/sigma_CIs.jpeg">
</p>

<br />

## Predicted vs. real IVI

Real logIVI values are plotted as points, and the relationship with the respective covariate is plotted as a trend with credible intervals.

<p align="center">
  <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/Figures/chickage-year.jpeg">
   <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/Figures/broodsize-year.jpeg">
</p>

<br />
<br />

## Next steps
* directions from Kim resulting from meeting 
* convert IVI to original scale so we can more intuitively interpret results. Change distribution to lognormal or gamma to handle original scale
* Model residuals with a metric of yearly environmental conditions, rather than flat category of year
* investigate options for estimating energetics of provisioning events

<br />
<br />

# Simulations

#### Gaussian distributed IVI - heterogenous residual variance among years

Simulate data and specify covariate effect sizes. We want residual variance to depend on year, so we'll draw residuals from a normal distribution N(0, sigma_year). The specific sigma associated with each year is defined in the vector, ```eps_sigma```.

##### simulate data
```r
b1       <- -1                # beta for chickage
b2       <- -2                # beta for broodsize
intercept <- 15               # intercept
eps_sigma = c(1,5,2,9,3,8,3)  # yearly variation in IVI

for (i in 1:nrow(d)){ 
      d$ivi[i] <- intercept + d$b0_nest[i] + b1 * d$chickage[i] + b2 * d$broodsize[i] + rnorm(1, 0, d$eps_sigma[i] + chickage*0.3)     
    }
```

##### nlme
```r
m <- lme(ivi ~ chickage:year + broodsize:year,   
         random = ~ 1 | nestID,           
         weights = varIdent(form = ~ 1|year),
         data = d)
 ```

##### JAGS

The variance structure is a bit hard to follow with all the priors, but what we're looking for is:

<p align="left">
  <img width="200" src="https://github.com/nuwcru/krmp_varSens/blob/master/documents/small.png">
</p>

```w``` is a weighting matrix, where the weights associated with each year are determined by the inverse of the ratio of each year to the first (reference) group.

All paramters and associated priors (all diffuse):

<p align="left">
  <img width="400" src="https://github.com/nuwcru/krmp_varSens/blob/master/documents/het_var-model.png">
</p>


And the JAGS model (for simplicity, all betas are packaged in a matrix together, and so are the covariate vectors):

```model{
    
    # Likelihood ~~~~~~~~~~~~~
    
        for (i in 1:N) {
          y[i]  ~ dnorm(mu[i], tau[years[i]])
          mu[i] <- inprod(beta[], X[i,]) + a[nest[i]]
    
        y_pred[i] ~ dnorm(mu[i], tau[years[i]])
        }
    
    # Priors ~~~~~~~~~~~~~~~
    
        for (i in 1:K) {beta[i] ~ dnorm(0,0.001)}
        
        for (i in 1:n_years){
        chSq[i]  ~ dgamma(0.5, 0.5)
        z[i]     ~ dnorm(0, 0.04)I(0,)
        sigma[i] <- z[i] / sqrt(chSq[i])
        tau[i]   <- pow(sigma[i], -2)
        }

    
    # prior for random intercept
        for (i in 1:n_nests) {a[i] ~ dnorm(a_bar, sigma_nest)}
    
    # prior for variance of random intercept
        a_bar ~ dnorm(0, 1.5)
        sigma_nest ~ dexp(1)
        
    }
```

Real IVI values plotted as points, mean predicted value from jags model, and credible intervals surrounding the mean. X axis is chickage.
<p align="left">
  <img width="800" src="https://github.com/nuwcru/krmp_varSens/blob/master/documents/yearly_variance.png">
</p>



### Contributors
* Rebekah McKinnon
* Kim Mathot
* Kevin Hawkshaw
* Erik Hedlin
* Alastair Franke


