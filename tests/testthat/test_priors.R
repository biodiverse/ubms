context("Prior specification and processing")

test_that("normal prior can be specified",{
  n <- normal()
  expect_equal(n, list(dist=1, par1=0, par2=2.5, par3=0, autoscale=TRUE))
  n <- normal(2, 0.5, autoscale=FALSE)
  expect_equal(n, list(dist=1, par1=2, par2=0.5, par3=0, autoscale=FALSE))
  n <- normal(c(0.5,1), c(0.1,0.2))
  expect_equal(n, list(dist=1, par1=c(0.5,1), par2=c(0.1,0.2), par3=0, autoscale=TRUE))
  n <- normal(c(0.5,1), 0.1) # no error
  expect_error(normal(c(0.5,1), c(0.1,0.2,0.3)))
  expect_error(normal(0.5, -0.2))
})

test_that("uniform prior can be specified",{
  n <- uniform()
  expect_equal(n, list(dist=2, par1=-5, par2=5, par3=0, autoscale=FALSE))
  n <- uniform(0.1, 0.5)
  expect_equal(n, list(dist=2, par1=0.1, par2=0.5, par3=0, autoscale=FALSE))
  n <- uniform(c(0.5,1), c(1,1.5))
  expect_equal(n, list(dist=2, par1=c(0.5,1), par2=c(1,1.5), par3=0, autoscale=FALSE))
  expect_error(uniform(c(0.5,1), c(1)))
  expect_error(uniform(2, -0.2))
})

test_that("null prior can be specified",{
  expect_equal(null_prior(), list(dist=0, par1=0, par2=0, par3=0, autoscale=FALSE))
})

test_that("expand_prior replicates parameters as needed",{
  n <- normal()
  expect_equal(expand_prior(n, 3),
               list(dist=1, par1=rep(0,3), par2=rep(2.5,3), par3=rep(0,3),
                    autoscale=TRUE))
  n <- normal(0, c(1,2))
  expect_equal(expand_prior(n, 2),
               list(dist=1, par1=c(0,0), par2=c(1,2), par3=c(0,0), autoscale=TRUE))
  expect_error(expand_prior(n, 3))
})

test_that("autoscale_prior scales prior by sd of covariate",{
  Xmat <- matrix(c(1,3,1,5), nrow=2)
  n <- expand_prior(normal(), 2)
  expect_equal(autoscale_prior(n, Xmat),
               list(dist=1, par1=c(0,0), par2=c(1.7678, 0.8839), par3=c(0,0),
                    autoscale=TRUE), tol=1e-4)
  n <- expand_prior(normal(autoscale=FALSE), 2)
  expect_equal(n, autoscale_prior(n, Xmat))
  expect_error(autoscale_prior(expand_prior(normal(), 3), Xmat))
})

test_that("autoscale_prior handles NAs",{
  Xmat <- matrix(c(1,3,NA,0,1,NA), nrow=2)
  n <- expand_prior(normal(), 3)
  expect_equal(autoscale_prior(n, Xmat),
               list(dist=1, par1=c(0,0,0), par2=c(1.76777,2.5,2.5), par3=c(0,0,0),
                    autoscale=TRUE), tol=1e-4)
  n <- expand_prior(normal(autoscale=FALSE), 3)
  expect_equal(n, autoscale_prior(n, Xmat))
})

test_that("process_coef_prior can expand and scale priors",{
  Xmat <- matrix(c(1,3,1,5), nrow=2)
  p <- normal()
  p2 <- process_coef_prior(p, Xmat)
  expect_equal(p2, autoscale_prior(expand_prior(p, nrow(Xmat)), Xmat))

  p <- normal(autoscale=FALSE)
  Xmat <- matrix(c(1,1),ncol=1)
  expect_equal(p, process_coef_prior(p, Xmat))

  Xmat <- matrix(c(1,3,1,5), nrow=2)
  colnames(Xmat) <- c("(Intercept)", "x1")
  expect_equal(process_coef_prior(normal(), Xmat),
               list(dist=1, par1=0, par2=0.88388, par3=0, autoscale=TRUE),
               tol=1e-4)

  # When there are no covariates, only an intercept
  Xmat <- Xmat[,-2,drop=FALSE]
  expect_equal(process_coef_prior(normal(), Xmat),
               list(dist=0, par1=NA, par2=NA, par3=NA, autoscale=TRUE))
})

test_that("process_int_prior identifies when there is no intercept",{
  Xmat <- matrix(c(1,3,1,5), nrow=2)
  colnames(Xmat) <- c("(Intercept)", "x1")
  expect_equal(process_int_prior(normal(), Xmat),
               list(dist=1,par1=0,par2=2.5,par3=0,autoscale=FALSE))

  Xmat <- matrix(c(1,3,1,5), nrow=2)
  colnames(Xmat) <- c("x2", "x1")
  expect_equal(process_int_prior(normal(), Xmat),
               list(dist=0, par1=NA, par2=NA, par3=NA, autoscale=FALSE))

  expect_error(process_int_prior(normal(c(0,0), 2.5), Xmat))
})

test_that("process_priors ubmsSubmodel method works",{
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,3,5)), ~x1, "plogis",
                     prior_intercept=normal(0,10), prior_coef=normal(0,2.5))
  expect_equal(process_priors(sm),
               list(prior_dist=c(1,1), prior_pars=matrix(c(0,10,0,0,1.25,0),nrow=3)))

  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,3,5)), ~x1, "plogis",
                     prior_intercept=normal(0,10), prior_coef=uniform(0,2.5))
  expect_equal(process_priors(sm),
               list(prior_dist=c(1,2), prior_pars=matrix(c(0,10,0,0,2.5,0),nrow=3)))

  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,3,5)), ~x1, "plogis",
                     prior_intercept=normal(0,10), prior_coef=student_t(1,0,2.5))
  expect_equal(process_priors(sm),
               list(prior_dist=c(1,3), prior_pars=matrix(c(0,10,0,0,1.25,1),nrow=3)))

  # when no covariates
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,3,5)), ~1, "plogis",
                     prior_intercept=normal(0,10), prior_coef=student_t(1,0,2.5))
  expect_equal(process_priors(sm),
               list(prior_dist=c(1,0), prior_pars=matrix(c(0,10,0),nrow=3)))

  # when no intercept
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,3,5)), ~x1-1, "plogis",
                     prior_intercept=normal(0,10), prior_coef=student_t(1,0,2.5))
  expect_equal(process_priors(sm),
               list(prior_dist=c(0,3), prior_pars=matrix(c(0,1.25,1),nrow=3)))

})

test_that("process_prior ubmsSubmodelScalar method works",{
  sm <- ubmsSubmodelScalar("Det", "det", "exp", prior_intercept=normal(0,10))
  expect_equal(process_priors(sm),
               list(prior_dist=c(1,0), prior_pars=matrix(c(0,10,0), nrow=3)))
})
