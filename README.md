# ubms: Unmarked Bayesian Models with Stan

In-development `R` package for fitting Bayesian hierarchical models of
animal occurrence and abundance. The package has a formula-based
interface compatible with
[unmarked](https://cran.r-project.org/web/packages/unmarked/index.html),
but the model is fit using MCMC with [Stan](https://mc-stan.org/)
instead of using maximum likelihood. Right now there are Stan versions
of unmarked functions `occu`, `occuRN`, and `pcount` (`stan_occu`,
`stan_occuRN`, `stan_pcount`).

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
    ##                 Estimate    SD   2.5%  97.5% n_eff Rhat
    ## (Intercept)        0.329 0.310 -0.278  0.957  1184    1
    ## x1                -0.464 0.117 -0.692 -0.238  5339    1
    ## sigma [1|group]    1.406 0.293  0.942  2.094  2432    1
    ## 
    ## Detection:
    ##             Estimate     SD  2.5% 97.5% n_eff  Rhat
    ## (Intercept)    0.383 0.0592 0.267 0.500  5807 1.000
    ## x2             0.586 0.0629 0.463 0.714  6628 0.999
    ## 
    ## WAIC: 2266.865

Assess model goodness-of-fit with a posterior predictive check, using
the MacKenzie-Bailey chi-square test:

``` r
(fm_fit <- gof(fm, draws=500, quiet=TRUE))
```

    ## MacKenzie-Bailey Chi-square 
    ## Point estimate = 30.225
    ## Posterior predictive p = 0.474

``` r
plot(fm_fit)
```

![](README_figs/README-gof-1.png)<!-- -->
