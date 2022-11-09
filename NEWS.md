# ubms v1.2.0

* New kfold function
* Add ability to not save log likelihoods in posterior
* Add new Laplace prior (thanks Justin Cally)
* Fix CRAN clang warnings
* Misc. bugfixes

# ubms v1.1.0

* Support for restricted spatial regression (RSR)
* Add vignette for RSR
* Update default priors and allow setting custom priors
* Support for offsets
* New utility functions and bugfixes

# ubms v1.0.2

* Make default priors less informative
* Fix issues in `sim_z` C++ code that triggered clang-UBSAN errors
* Fix various minor issues with checks/tests on R-oldrel 

# ubms v1.0.1

* Added vignette comparing ubms and JAGS output
* Added configuration for pkgdown site
* Small adjustments to wording in DESCRIPTION and docs for CRAN submission

# ubms v0.1.9

* Add time-to-detection occupancy model (`stan_occuTTD`)
* Wrote examples for all fitting functions
* Many new tests added
* Wrote a vignette demonstrating random effects

# ubms v0.1.8

* Added distance-sampling model (`stan_distsamp`)
* Added multinomial-Poisson mixture model (`stan_multinomPois`)
* Wrote overview vignette

# ubms v0.1.4

* Added dynamic occupancy model (`stan_colext`)
* Consolidated single-season models into one Stan file to speed up compilation
* Using rasters as `predict` input is now supported
