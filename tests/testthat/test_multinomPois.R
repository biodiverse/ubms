context("stan_multinomPois function and methods")

skip_on_cran()

#Simulate dataset
set.seed(567)
nSites <- 50
lambda <- 10
p1 <- 0.5
p2 <- 0.3
cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
set.seed(9023)
N <- rpois(nSites, lambda)
y <- matrix(NA, nSites, 3)
for(i in 1:nSites) {
  y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
}

# Fit model
observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
umf_double <- suppressWarnings(unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
                          siteCovs=data.frame(x=rnorm(nSites)),
                          type="double"))
umf_double_na <- umf_double
umf_double_na@y[1,] <- NA
umf_double_na@y[2,1] <- NA


good_fit <- TRUE
tryCatch({
fit_double <- suppressWarnings(stan_multinomPois(~observer-1 ~x, umf_double[1:10,],
                               chains=2, iter=100, refresh=0))

fit_double_na <- suppressWarnings(stan_multinomPois(~observer-1 ~x, umf_double_na[1:10,],
                               chains=2, iter=100, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})

## Removal
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
              siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
              type = "removal")

ovenFrame_na <- ovenFrame
ovenFrame_na@y[1,] <- NA
ovenFrame_na@y[2,1] <- NA

tryCatch({
fit_rem <- suppressWarnings(stan_multinomPois(~1~ufc, ovenFrame[1:20,],
                            chains=2, iter=200, refresh=0))

fit_rem_na <- suppressWarnings(stan_multinomPois(~1~ufc, ovenFrame_na[1:10,],
                            chains=2, iter=200, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})

skip_if(!good_fit, "Test setup failed")

test_that("stan_multinomPois output structure is correct",{
  expect_is(fit_double, "ubmsFitMultinomPois")
  expect_is(fit_double@response, "ubmsResponseMultinomPois")
  expect_equal(nsamples(fit_double), 100)
  expect_is(fit_rem, "ubmsFitMultinomPois")
})

test_that("stan_multinomPois produces accurate results",{
  skip_on_cran()
  skip_on_ci()
  skip_on_covr()
  set.seed(123)
  fit_double_long <- suppressWarnings(stan_multinomPois(~observer-1 ~x, umf_double,
                               chains=3, iter=300, refresh=0))
  fit_rem_long <- suppressWarnings(stan_multinomPois(~1~ufc, ovenFrame,
                               chains=3, iter=300, refresh=0))
  um_rem <- multinomPois(~1~ufc, ovenFrame)
  um_double <- multinomPois(~observer-1~x, umf_double)

  b <- c(log(lambda), 0, log(0.5/(1-0.5)), log(0.3/(1-0.3)))
  #similar to trutih
  expect_RMSE(coef(fit_double_long), b, 0.2)
  #similar to unmarked
  expect_RMSE(coef(fit_double_long), coef(um_double), 0.05)
  #similar to previous known values
  expect_RMSE(coef(fit_double_long), c(2.2268,0.1133,0.1853,-0.5620), 0.05)

  #Removal
  expect_RMSE(coef(fit_rem_long), coef(um_rem),  0.03)
  expect_RMSE(coef(fit_rem_long), c(0.11398,0.17389,0.26453), 0.05)
})

test_that("stan_multinomPois handles NA values",{
  expect_is(coef(fit_double_na), "numeric")
  expect_is(coef(fit_rem_na), "numeric")
})

test_that("ubmsFitMultinomPois gof method works",{ ##here
  set.seed(123)
  g <- gof(fit_double, draws=5, quiet=TRUE)
  expect_between(g@estimate, 10, 35)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitMultinomPois gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_double_na, draws=5, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("ubmsFitMultinomPois predict method works",{
  pr <- predict(fit_double_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(10, 4))
  expect_between(pr[1,1], 5, 15)
  pr <- predict(fit_double_na, "det")
  expect_equal(dim(pr), c(10*obsNum(umf_double),4))
  expect_between(pr[1,1], 0, 1)
  #with newdata
  nd <- data.frame(x=c(0,1))
  pr <- predict(fit_double_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_between(pr[1,1], 5, 15)
})

test_that("ubmsFitMultinomPois sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit_double, 5)
  zz <- sim_z(fit_double, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 10))
  expect_between(mean(zz), 5, 15)

  set.seed(123)
  pz <- posterior_predict(fit_double, "z", draws=5)
  expect_equivalent(zz, pz)
})

test_that("stan_multinomPois sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit_double, 5)
  yy <- sim_y(fit_double, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), 10*3))
  expect_between(mean(yy), 1, 10)
  set.seed(123)
  py <- posterior_predict(fit_double, "y", draws=5)
  expect_equivalent(yy, py)

  #set.seed(123)
  #samples <- get_samples(fit_rem, 5)
  #yy <- sim_y(fit_rem, samples, re.form=NULL)
  #expect_is(yy, "matrix")
  #expect_equal(dim(yy), c(length(samples), 20*obsNum(ovenFrame)))
  #expect_true(between(mean(yy), 0, 5))
})

test_that("Posterior sim methods for ubmsFitMultinomPois work with NAs",{
  zna <- posterior_predict(fit_double_na, "z", draws=3)
  expect_equal(dim(zna), c(3,10))
  expect_true(any(is.na(zna)))
  yna <- posterior_predict(fit_double_na, "y", draws=3)
  expect_equal(dim(yna), c(3,10*3))
  expect_equal(sum(is.na(yna[1,])), 6)
  expect_equal(sum(is.na(yna[2,])), 6)
})

test_that("Posterior linear pred methods work for ubmsFitMultinomPois",{
  set.seed(123)
  samples <- get_samples(fit_double, 3)
  lp1 <- sim_lp(fit_double, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 10))
  set.seed(123)
  pl <- posterior_linpred(fit_double, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccu",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit_double, "state", draws=5)
  ft2 <- ubms_fitted(fit_double, "det", draws=5)
  expect_equal(dim(ft), c(5,10))
  expect_equal(dim(ft2), c(5,30))

  res <- ubms_residuals(fit_double, "state", draws=5)
  res2 <- ubms_residuals(fit_double, "det", draws=5)
  expect_equal(dim(res), c(5,10))
  expect_equal(dim(res2), c(5,30))

  pdf(NULL)
  rp <- plot_residuals(fit_double, "state")
  rp2 <- plot_residuals(fit_double, "det")
  rp3 <- ubms_plot(fit_double)
  mp <- plot_marginal(fit_double, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(rp2, "gg")
  expect_is(rp3, "gtable")
  expect_is(mp, "gg")
})

test_that("get_pifun_type returns correct value",{
  expect_equal("double", get_pifun_type(umf_double))
  expect_equal("removal", get_pifun_type(ovenFrame))
  umf_broken <- umf_double
  umf_broken@piFun <- "fake"
  expect_error(get_pifun_type(umf_broken))
})

test_that("ubmsResponseMultinomPois find_missing method works",{
  resp <- ubmsResponseMultinomPois(getY(umf_double)[1:10,], "double", "P")
  subs <- fit_double@submodels
  expect_true(!all(is.na(find_missing(resp, subs))))

  resp <- ubmsResponseMultinomPois(getY(umf_double_na)[1:10,], "double", "P")
  subs <- fit_double_na@submodels
  expect_equal(find_missing(resp, subs), c(rep(TRUE,6),rep(FALSE,24)))

  resp <- ubmsResponseMultinomPois(getY(ovenFrame_na)[1:10,], "removal", "P")
  subs <- fit_rem_na@submodels
  expect_equal(find_missing(resp, subs), c(rep(TRUE,8),rep(FALSE,10*4-8)))
})

test_that("ubmsResponseMultinomPois update_missing method works",{
  resp <- ubmsResponseMultinomPois(getY(umf_double_na)[1:10,], "double", "P")
  subs <- fit_double_na@submodels
  subs@submodels$state@missing <- rep(FALSE, 10)
  subs@submodels$det@missing <- rep(FALSE, length(subs@submodels$det@missing))
  out <- update_missing(subs, resp)
  expect_equal(out@submodels$state@missing, c(rep(TRUE, 2), rep(FALSE, 8)))
  expect_equal(out@submodels$det@missing, c(rep(TRUE, 4), rep(FALSE, 16)))

  resp <- ubmsResponseMultinomPois(getY(ovenFrame_na)[1:10,], "removal", "P")
  subs <- fit_rem_na@submodels
  subs@submodels$state@missing <- rep(FALSE, numSites(ovenFrame_na))
  subs@submodels$det@missing <- rep(FALSE, length(subs@submodels$det@missing))
  out <- update_missing(subs, resp)
  expect_equal(out@submodels$state@missing, c(rep(TRUE, 2), rep(FALSE, 8)))
  expect_equal(out@submodels$det@missing, c(rep(TRUE, 8), rep(FALSE, 32)))
})

test_that("get_pi_for_multinom function works",{
  pi_out <- get_pi_for_multinom(fit_double, 1:3)
  expect_is(pi_out, "array")
  expect_equal(dim(pi_out), c(10,3+1,3))
  expect_equal(rowSums(pi_out[,,1]), rep(1,10))

  pi_out <- get_pi_for_multinom(fit_rem, 1:3)
  expect_is(pi_out, "array")
  expect_equal(dim(pi_out), c(20,4+1,3))
  expect_equal(rowSums(pi_out[,,1]), rep(1,20))
})

test_that("getP and sim_p for ubmsFitMultinomPois work",{
  p <- sim_p(fit_double, 1:3)
  expect_equal(dim(p), c(3,10*3))
  p <- sim_p(fit_double_na, 1:3)
  expect_equal(dim(p), c(3,10*3))

  gp <- getP(fit_double, 3)
  expect_equal(dim(gp), c(10,3,3))
  gp <- getP(fit_double_na, 3)
  expect_equal(dim(gp), c(10,3,3))

  p <- sim_p(fit_rem, 1:3)
  expect_equal(dim(p), c(3,20*4))
  p <- sim_p(fit_rem_na, 1:3)
  expect_equal(dim(p), c(3,10*4))

  gp <- getP(fit_rem, 3)
  expect_equal(dim(gp), c(20,4,3))
  gp <- getP(fit_rem_na, 3)
  expect_equal(dim(gp), c(10,4,3))
})

test_that("multinomPois spatial works", {
  skip_on_cran()
  umf2 <- umf_double
  umf2@siteCovs$x1 <- umf2@siteCovs$x
  umf2@siteCovs$x <- runif(numSites(umf2), 0, 10)
  umf2@siteCovs$y <- runif(numSites(umf2), 0, 10)
  fit_spat <- suppressMessages(suppressWarnings(stan_multinomPois(~1~x1+RSR(x,y,1),
                umf2[1:20,], chains=2, iter=100, refresh=0)))
  expect_is(fit_spat@submodels@submodels$state, "ubmsSubmodelSpatial")
  expect_equal(names(coef(fit_spat))[3], "state[RSR [tau]]")

  ps <- plot_spatial(fit_spat)
  expect_is(ps, "gg")
})
