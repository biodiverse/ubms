# ubms: Unmarked Bayesian Models with Stan

<!-- badges: start -->

[![R build
status](https://github.com/hmecology/ubms/workflows/R-CMD-check/badge.svg)](https://github.com/hmecology/ubms/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/ubms)](https://cran.r-project.org/web/packages/ubms/index.html)
<!-- badges: end -->

`ubms` is an `R` package for fitting Bayesian hierarchical models of
animal occurrence and abundance. The package has a formula-based
interface compatible with
[unmarked](https://cran.r-project.org/web/packages/unmarked/index.html),
but the model is fit using MCMC with [Stan](https://mc-stan.org/)
instead of using maximum likelihood. Currently there are Stan versions
of unmarked functions `occu`, `occuRN`, `colext`, `occuTTD`, `pcount`,
`distsamp`, and `multinomPois`. These functions follow the `stan_`
prefix naming format established by
[rstanarm](https://cran.r-project.org/web/packages/rstanarm/index.html).
For example, the Stan version of the `unmarked` function `occu` is
`stan_occu`.

Advantages compared to `unmarked`:

1.  Obtain posterior distributions of parameters and derived parameters
2.  Include random effects in parameter formulas (same syntax as `lme4`)
3.  Assess model fit using WAIC and LOO via the
    [loo](https://cran.r-project.org/web/packages/loo/index.html)
    package

Disadvantages compared to `unmarked`:

1.  MCMC is slower than maximum likelihood
2.  Not all model types are supported
3.  Potential for convergence issues

## Installation

`ubms` is on
[CRAN](https://cran.r-project.org/web/packages/ubms/index.html):

``` r
install.packages("ubms")
```

Alternatively, the latest development version can be installed from
Github:

``` r
# install.packages("devtools")
devtools::install_github("hmecology/ubms")
```

## Package Overview

A detailed vignette for the package is available [here](https://cran.r-project.org/web/packages/ubms/vignettes/ubms.html).
