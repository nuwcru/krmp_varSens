
## Variance Sensitivity

### To Do
- [x] clean repo
- [x] re-calculate IVI
- [x] translate jags to gamma distribution
- [x] incorporate random slope / intercept for year with correlation
- [ ] examine brood mass options



### Intro
Investigating variance sensitivity among breeding Peregrine Falcons in response to environmental conditions and brood demand.

### File Organization

```
├── Data
│     ├── Clean IVI years 2013-2019.csv  <- data used for analysis - cleaned by Rebekah
│     ├── ivi_eh.csv                     <- new ivi data with fixed errors
│     ├── ____ .parquet                  <- all parquet files are jags outputs converted to dataframes for easy manipulation/viz
│     ├── 
│     ├── 
│     └── 
├── Scripts 
│     ├──  simulation.R         <- simulate data to ensure variance is being captured properly by whatever model we use
│     │                            - nlme model for simplicity / jags for thorough understanding and benefits from bayesian framework
│     ├──  jags.R               <- Real Data using similar model structure to that tested in simulation.R
│     ├──  jags_v_lme4.R        <- Comparisons between lme4 and jags models
│     ├──  jags_gamma.R         <- Gamma distribution, heterogenous resids, random intercepts 
│     ├──  jags_gaus.R          <- log transformed ivi, heterogenous resids, random intercepts
│     ├──  02_jags_gaus.R       <- log transformed ivi, het resids, random intercepts, random slopes and correlation between ranefs
│     ├──  brms.R               <- 02_jags_gaus.R model minus het resids. Quick modeling of ranef covariance using lme4 syntax, with a stan backend
│     └──  ivi_scratch.R        <- Programmatically calculate IVI, produces ivi_eh.csv
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

### Current status

**02_jags_gaus.R** runs the most up-to-date model based on comments from Kim / Rebekah.

The model includes 
* **random intercept** for year and nest site, and **random slopes** for chickage and chicks nested in years. 
* The model estimates ranef from a **covariance matrix of slopes/intercepts** which allows us to estimate the level of covariance between yearly intercepts and the effects of chickage/chicks. 
* The model also estimates **independent residual variance by year**, but this can be commented out if need be. It's important to note that estimating residual variance independantly by year eats into the variation explained by other parameters such as chicks and chickage. This is demonstrated below.


# Heterogenous Residuals

This is the 02_jags_gaus model with heterogenous residuals turned on. 

*Residual variance by year. The estimate visualized here is sigma, and associated uncertainty in the estimate of sigma stemming from the shape of the posterior:*

<p float="center">
  <img src="Figures/het_resid.png" width="600" />
</p>

*The effect of chickage within the same model where yearly variation in the effect of chickage is captured by the random slope.*

<p float="center">
  <img src="Figures/het_resid_chickage.png" width="600" />
</p>


*Covariance between yearly intercept and random slopes for chicks and chickage*


<p float="center">
  <img src="Figures/ranef_cov.png" width="600" />
</p>




# Homogenous Residuals

If we assume residual variance is homogenous across years, varying slopes for chickage then steps up to capture that yearly variation.

<p float="center">
  <img src="Figures/hom_resid_chickage.png" width="600" />
</p>

We plot the residuals by year below to see if there are yearly patterns not explained in the model. It looks pretty uniform

<p float="center">
  <img src="Figures/hom_residuals.png" width="800" />
</p>

Still moderate covariance between intercepts and slope, and the covariance is in the direction you would predict, lower intercept leads to more pronounced decrease in IVI as brood size increases. 
<p float="center">
  <img src="Figures/hom_resid_covariance.png" width="600" />
</p>





### Conclusion

Residual variance among years is captured by letting our slopes and intercepts vary. Moderate covariance between yearly intercepts and slopes.

### Contributors
* Rebekah McKinnon
* Kim Mathot
* Kevin Hawkshaw
* Erik Hedlin
* Alastair Franke


