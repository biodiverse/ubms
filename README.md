ubms: Unmarked Bayesian Models with Stan
========================================

Experimental `R` package for fitting hierarchical models of animal
occurrence and abundance. The package has a formula-based interface
compatible with
[unmarked](https://cran.r-project.org/web/packages/unmarked/index.html),
but the model is fit using MCMC with [Stan](https://mc-stan.org/)
instead of using maximum likelihood.

Advantages over `unmarked`:

1.  Easily obtain posterior distributions of parameters and derived
    parameters
2.  Include random effects in parameter formulas (same syntax as `lme4`)

Disadvantages over `unmarked`:

1.  MCMC will always be slower than maximum likelihood
2.  Potential convergence issues
3.  Only one model (single-season occupancy) is available for now

### Example

``` r
#Simulate occupancy data including random effect on occupancy
set.seed(123)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

re_fac <- factor(sample(letters[1:10], 500, replace=T))
dat_occ$group <- re_fac
re <- rnorm(10, 0, 0.25)
re_idx <- as.numeric(re_fac)

idx <- 1
for (i in 1:500){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$x1[i] + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1, 
                    plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}

library(ubms)

#Create unmarked frame
umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

#Fit model with random intercept
(fm <- stan_occu(~x2 ~x1 + (1|group), umf, refresh=0))
```

    ## 
    ## Call:
    ## stan_occu(formula = ~x2 ~ x1 + (1 | group), data = umf, refresh = 0)
    ## 
    ## Occupancy:
    ##               Estimate    SD    2.5%  97.5% n_eff Rhat
    ## (Intercept)      0.486 0.148  0.1905  0.785  2298    1
    ## x1              -0.942 0.121 -1.1916 -0.714  3760    1
    ## sigma [group]    0.318 0.184  0.0598  0.744   673    1
    ## 
    ## Detection:
    ##             Estimate     SD  2.5% 97.5% n_eff  Rhat
    ## (Intercept)    0.247 0.0562 0.133 0.357  3848 1.000
    ## x2             0.492 0.0598 0.376 0.612  3756 0.999
    ## 
    ## WAIC: 2522.421
