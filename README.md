# Variance Sensitivity

### Intro
Investigating variance sensitivity among breeding Peregrine Falcons in response to environmental conditions and brood demand.

### File Organization

```
├── Data
│     ├── 
│     ├── 
│     ├── 
│     ├── 
│     ├── 
│     └── 
├── Scripts 
│     ├──  simulation.R    <- simulate data to ensure variance is being captured properly by whatever model we use
│     │                            - nlme model for simplicity / jags for thorough understanding and benefits from bayesian framework
│     ├── 
│     └── 
└── Figures
      ├── 
      └── 

```

### Simulations

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

Where ```w``` is a weighting matrix, where the weights associated with each year are determined by the inverse of the ratio of each year to the first (reference) group.

All paramters and associated priors (all diffuse):

<p align="left">
  <img width="400" src="https://github.com/nuwcru/krmp_varSens/blob/master/documents/het_var-model.png">
</p>


And the JAGS model:

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




### Contributors
* Rebekah McKinnon
* Kim Mathot
* Kevin Hawkshaw
* Erik Hedlin
* Alastair Franke


