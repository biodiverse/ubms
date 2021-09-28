context("stan_occuRN function and methods")

skip_on_cran()

#Simulate dataset
set.seed(123)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
N <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

idx <- 1
for (i in 1:500){
  N[i] <- rpois(1, exp(b[1]+b[2]*dat_occ$x1[i]))
  for (j in 1:5){
    r <- plogis(b[3] + b[4]*dat_p$x2[idx])
    p <- 1 - (1-r)^N[i]
    y[i,j] <- rbinom(1, 1, p)
    idx <- idx + 1
  }
}

umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occuRN(~x2~x1, umf[1:10,], K=15,
                                    chains=2, iter=100, refresh=0))

fit_na <- suppressWarnings(stan_occuRN(~x2~x1, umf2[1:10,], K=15,
                                       chains=2, iter=100, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})

skip_if(!good_fit, "Test setup failed")

test_that("stan_occuRN output structure is correct",{
  expect_is(fit, "ubmsFitOccuRN")
  expect_is(fit, "ubmsFitOccu")
  expect_equal(nsamples(fit), 100)
})

test_that("stan_occuRN produces accurate results",{
  skip_on_cran()
  skip_on_ci()
  skip_on_covr()
  set.seed(123)
  fit_long <- suppressWarnings(stan_occuRN(~x2~x1, umf[1:200,], K=15, chains=2,
                                           iter=200, refresh=0))
  fit_unm <- occuRN(~x2~x1, umf[1:200,], K=15)
  #similar to truth
  expect_RMSE(coef(fit_long), b, 0.1)
  #similar to unmarked
  expect_RMSE(coef(fit_long), coef(fit_unm), 0.02)
  #similar to previous known values
  expect_RMSE(coef(fit_long), c(0.4838,-0.6449,0.2749,0.5012), 0.05)
})

test_that("stan_occuRN handles NA values",{
  expect_is(coef(fit_na), "numeric")
})

test_that("ubmsFitOccuRN gof method works",{
  set.seed(123)
  g <- gof(fit, draws=5, quiet=TRUE)
  expect_between(g@estimate, 30, 50)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitOccuRN gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_na, draws=5, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("ubmsFitOccuRN predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(10, 4))
  expect_between(pr[1,1], 0.5, 3.5)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(10*obsNum(umf2),4))
  expect_between(pr[1,1], 0, 1)
  #with newdata
  nd <- data.frame(x1=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_between(pr[1,1], 0.5, 3.5)
})

test_that("ubmsFitOccuRN sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 10))
  expect_between(mean(zz), 1, 3)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=5)
  expect_equivalent(zz, pz)
})

test_that("ubmsFitOccuRN sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), 10*obsNum(umf)))
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=5)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitOccuRN work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,10))
  expect_true(all(is.na(zna[,1])))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3, 10*obsNum(umf2)))
  expect_equal(sum(is.na(yna[1,])), 6)
  expect_equal(sum(is.na(yna[2,])), 6)
})

test_that("Posterior linear pred methods work for ubmsFitOccuRN",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 10))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccuRN",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=5)
  ft2 <- ubms_fitted(fit, "det", draws=5)
  expect_equal(dim(ft), c(5, 10))
  expect_equal(dim(ft2), c(5, 10*obsNum(umf)))

  res <- ubms_residuals(fit, "state", draws=5)
  res2 <- ubms_residuals(fit, "det", draws=5)
  expect_equal(dim(res), c(5, 10))
  expect_equal(dim(res2), c(5, 10*obsNum(umf)))

  pdf(NULL)
  rp <- plot_residuals(fit, "state")
  rp2 <- plot_residuals(fit, "det")
  rp3 <- ubms_plot(fit)
  mp <- plot_marginal(fit, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(rp2, "gg")
  expect_is(rp3, "gtable")
  expect_is(mp, "gg")
})

test_that("occuRN spatial works", {
  skip_on_cran()
  umf2 <- umf
  umf2@siteCovs$x <- runif(numSites(umf2), 0, 10)
  umf2@siteCovs$y <- runif(numSites(umf2), 0, 10)
  fit_spat <- suppressMessages(suppressWarnings(stan_occuRN(~1~x1+RSR(x,y,1),
                umf2[1:20,], K=15, chains=2, iter=50, refresh=0)))
  expect_is(fit_spat@submodels@submodels$state, "ubmsSubmodelSpatial")
  expect_equal(names(coef(fit_spat))[3], "state[RSR [tau]]")

  ps <- plot_spatial(fit_spat)
  expect_is(ps, "gg")
})
