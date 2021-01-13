## Test environments

* Local: Ubuntu 18.04 (R-release and R-devel), Windows (R-release)
* Github actions CI: Ubuntu 20.04 R-release, Ubuntu 20.04 R-devel, Windows R-release, MacOS R-release
* Win-builder (R-devel)

## R CMD check results (R-devel on Ubuntu 18.04)

There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Ken Kellner <contact@kenkellner.com>’

New submission

This is the first time I have submitted this package.

* checking installed package size ... NOTE
  installed size is 63.0Mb
  sub-directories of 1Mb or more:
    R      1.4Mb
    libs  60.6Mb

As per the guidelines for developing R packages that interface with Stan (https://mc-stan.org/rstantools/articles/developer-guidelines.html), all Stan models I have developed for this package are pre-compiled when the package is built. This allows the models to be run on Mac/Windows without a C++ compiler available and speeds up model runtime at the cost of increasing installed package size. I have followed suggested best practices to combine my models into as few Stan model files as possible to speed up compilation and minimize size.

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

I believe GNU make is required to use Rcpp and compiled Stan code.