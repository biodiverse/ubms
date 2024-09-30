# ubms 1.2.7

* Package links now go to biodiverse Github org
* Use reformulas instead of lme4
* Return eta from ranef when submodel has spatial component
* Support terra format rasters in predict
* Minor bugfixes

# ubms 1.2.6

* Fix deprecated Stan code

# ubms 1.2.5

* Fix issue with version number not being a string
* Adding missing roxygen tag

# ubms 1.2.4

* Adjust Makevars to fix problems with clang 16

# ubms 1.2.3

* Fix bug in Stan code revealed by new StanHeaders
* Drop C++ specification per CRAN instructions
* Remove deprecated ggplot2 functions
* Small updates to overview vignette

# ubms 1.2.2

* Fix broken URLs in vignettes
* New kfold function
* Add ability to not save log likelihoods in posterior
* Add new Laplace prior (thanks Justin Cally)
* Better handling of factors with levels that don't appear in data
* Fix CRAN clang warnings
* Misc. bugfixes

# ubms 1.1.0

* Support for restricted spatial regression (RSR)
* Add vignette for RSR
* Update default priors and allow setting custom priors
* Support for offsets
* New utility functions and bugfixes

# ubms 1.0.2

* Make default priors less informative
* Fix issues in `sim_z` C++ code that triggered clang-UBSAN errors
* Fix various minor issues with checks/tests on R-oldrel 

# ubms 1.0.1

* Added vignette comparing ubms and JAGS output
* Added configuration for pkgdown site
* Small adjustments to wording in DESCRIPTION and docs for CRAN submission

# ubms 0.1.9

* Add time-to-detection occupancy model (`stan_occuTTD`)
* Wrote examples for all fitting functions
* Many new tests added
* Wrote a vignette demonstrating random effects

# ubms 0.1.8

* Added distance-sampling model (`stan_distsamp`)
* Added multinomial-Poisson mixture model (`stan_multinomPois`)
* Wrote overview vignette

# ubms 0.1.4

* Added dynamic occupancy model (`stan_colext`)
* Consolidated single-season models into one Stan file to speed up compilation
* Using rasters as `predict` input is now supported
