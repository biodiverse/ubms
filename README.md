# ubms: Unmarked Bayesian Models with Stan

[![Build
Status](https://travis-ci.org/kenkellner/ubms.svg?branch=master)](https://travis-ci.org/kenkellner/ubms)
[![codecov](https://codecov.io/gh/kenkellner/ubms/branch/master/graph/badge.svg)](https://codecov.io/gh/kenkellner/ubms)

In-development `R` package for fitting Bayesian hierarchical models of
animal occurrence and abundance. The package has a formula-based
interface compatible with
[unmarked](https://cran.r-project.org/web/packages/unmarked/index.html),
but the model is fit using MCMC with [Stan](https://mc-stan.org/)
instead of using maximum likelihood. Right now there are Stan versions
of unmarked functions `occu`, `occuRN`, `colext`, and `pcount`
(`stan_occu`, `stan_occuRN`, `stan_colext`, `stan_pcount`).

Advantages compared to `unmarked`:

1.  Obtain posterior distributions of parameters and derived parameters
2.  Include random effects in parameter formulas (same syntax as `lme4`)

Disadvantages compared to `unmarked`:

1.  MCMC is slower than maximum likelihood
2.  Limited selection of model types
3.  Potential convergence issues

### Example

Simulate occupancy data including a random effect on occupancy:

``` r
library(ubms)

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
```

Create `unmarked` frame:

``` r
umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)
```

Fit a model with a random intercept:

``` r
options(mc.cores=3) #number of parallel cores to use
(fm <- stan_occu(~x2 ~x1 + (1|group), umf, refresh=0))
```

    ## 
    ## Call:
    ## stan_occu(formula = ~x2 ~ x1 + (1 | group), data = umf, refresh = 0)
    ## 
    ## Occupancy:
    ##                 Estimate    SD   2.5%  97.5% n_eff  Rhat
    ## (Intercept)        0.323 0.304 -0.293  0.905  1074 1.001
    ## x1                -0.466 0.117 -0.695 -0.244  5627 0.999
    ## sigma [1|group]    1.404 0.278  0.946  2.022  3014 1.002
    ## 
    ## Detection:
    ##             Estimate     SD  2.5% 97.5% n_eff Rhat
    ## (Intercept)    0.381 0.0591 0.266 0.496  6364    1
    ## x2             0.587 0.0603 0.471 0.708  6832    1
    ## 
    ## LOOIC: 2268.389

Examine residuals for occupancy and detection submodels (following
[Wright et al. 2019](https://doi.org/10.1002/ecy.2703)). Each panel
represents a draw from the posterior distribution.

``` r
plot(fm)
```

![](README_figs/README-resids-1.png)<!-- -->

Assess model goodness-of-fit with a posterior predictive check, using
the MacKenzie-Bailey chi-square test:

``` r
(fm_fit <- gof(fm, draws=500, quiet=TRUE))
```

    ## MacKenzie-Bailey Chi-square 
    ## Point estimate = 30.054
    ## Posterior predictive p = 0.482

``` r
plot(fm_fit)
```

![](README_figs/README-gof-1.png)<!-- -->
