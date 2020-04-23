# ubms: Unmarked Bayesian Models with Stan

In-development `R` package for fitting Bayesian hierarchical models of
animal occurrence and abundance. The package has a formula-based
interface compatible with
[unmarked](https://cran.r-project.org/web/packages/unmarked/index.html),
but the model is fit using MCMC with [Stan](https://mc-stan.org/)
instead of using maximum likelihood.

Advantages compared to `unmarked`:

1.  Obtain posterior distributions of parameters and derived parameters
2.  Include random effects in parameter formulas (same syntax as `lme4`)

Disadvantages compared to `unmarked`:

1.  MCMC is slower than maximum likelihood
2.  Only a few models available for now, and some classes of models
    might never be practical
3.  Potential convergence issues
4.  Does not handle missing values

### Example

``` r
#Simulate occupancy data including random effect on occupancy
set.seed(123)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

re_fac <- factor(sample(letters[1:26], 500, replace=T))
dat_occ$group <- re_fac
re <- rnorm(26, 0, 1.2)
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
    ##                 Estimate    SD   2.5%  97.5% n_eff  Rhat
    ## (Intercept)        0.314 0.288 -0.298  0.870  1118 1.001
    ## x1                -0.463 0.117 -0.690 -0.242  8090 0.999
    ## sigma [1|group]    1.392 0.287  0.937  2.063  2723 1.001
    ## 
    ## Detection:
    ##             Estimate     SD  2.5% 97.5% n_eff  Rhat
    ## (Intercept)    0.382 0.0608 0.263 0.503  6634 0.999
    ## x2             0.587 0.0615 0.469 0.708  9528 1.000
    ## 
    ## WAIC: 2267.501
