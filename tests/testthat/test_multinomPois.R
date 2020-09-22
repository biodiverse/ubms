context("stan_multinomPois function and methods")

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

fit_double <- suppressWarnings(stan_multinomPois(~observer-1 ~x, umf_double,
                               chains=3, iter=300, refresh=0))

um_double <- multinomPois(~observer-1~x, umf_double)

umf_double_na <- umf_double
umf_double_na@y[1,] <- NA
umf_double_na@y[2,1] <- NA

set.seed(123)
fit_double_na <- suppressWarnings(stan_multinomPois(~observer-1 ~x, umf_double_na,
                               chains=3, iter=300, refresh=0))

## Removal
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
              siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
              type = "removal")

set.seed(123)
fit_rem <- suppressWarnings(stan_multinomPois(~1~ufc, ovenFrame,
                            chains=3, iter=300, refresh=0))

um_rem <- multinomPois(~1~ufc, ovenFrame)

ovenFrame_na <- ovenFrame
ovenFrame_na@y[1,] <- NA
ovenFrame_na@y[2,1] <- NA

set.seed(123)
fit_rem_na <- suppressWarnings(stan_multinomPois(~1~ufc, ovenFrame_na,
                            chains=3, iter=300, refresh=0))


test_that("stan_multinomPois output structure is correct",{
  expect_is(fit_double, "ubmsFitMultinomPois")
  expect_is(fit_double@response, "ubmsResponseMultinomPois")
  expect_equal(nsamples(fit_double), 450)
  expect_is(fit_rem, "ubmsFitMultinomPois")
})

test_that("stan_multinomPois produces accurate results",{
  b <- c(log(lambda), 0, log(0.5/(1-0.5)), log(0.3/(1-0.3)))
  #similar to truth
  expect_equal(as.vector(coef(fit_double)), b, tol=0.3)
  #similar to unmarked
  expect_equivalent(as.vector(coef(fit_double)), coef(um_double), tol=0.05)
  #similar to previous known values
  expect_equal(as.vector(coef(fit_double)), c(2.22254,0.11881,0.17918,-0.56747), tol=1e-4)

  #Removal
  expect_equivalent(as.vector(coef(fit_rem)), coef(um_rem), tol=0.05)
  expect_equal(as.vector(coef(fit_rem)), c(0.11004,0.17587,0.24698), tol=1e-4)
})

test_that("stan_multinomPois handles NA values",{
  expect_equal(as.vector(coef(fit_double)), as.vector(coef(fit_double_na)), tol=0.1)
  expect_equal(as.vector(coef(fit_rem)), as.vector(coef(fit_rem_na)), tol=0.1)
})

test_that("ubmsFitMultinomPois gof method works",{ ##here
  set.seed(123)
  g <- gof(fit_double, draws=10, quiet=TRUE)
  expect_equal(g@estimate, 155.3, tol=0.1)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitMultinomPois gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_double_na, draws=10, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("ubmsFitMultinomPois predict method works",{
  pr <- predict(fit_double_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(numSites(umf_double), 4))
  expect_equivalent(pr[1,1], 9.50708, tol=0.01)
  pr <- predict(fit_double_na, "det")
  expect_equal(dim(pr), c(numSites(umf_double)*obsNum(umf_double),4))
  expect_equivalent(pr[1,1], 0.5487, tol=0.01)
  #with newdata
  nd <- data.frame(x=c(0,1))
  pr <- predict(fit_double_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equivalent(pr[1,1], 9.28689, tol=0.01)
})

test_that("ubmsFitMultinomPois sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit_double, 10)
  zz <- sim_z(fit_double, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), numSites(umf_double)))
  expect_equal(mean(zz), 9.834)

  set.seed(123)
  pz <- posterior_predict(fit_double, "z", draws=10)
  expect_equivalent(zz, pz)
})

test_that("stan_occu sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit_double, 10)
  yy <- sim_y(fit_double, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), numSites(umf_double)*3))
  expect_equal(mean(yy), mean(umf_double@y), tol=0.1)
  set.seed(123)
  py <- posterior_predict(fit_double, "y", draws=10)
  expect_equivalent(yy, py)

  set.seed(123)
  samples <- get_samples(fit_rem, 10)
  yy <- sim_y(fit_rem, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), numSites(ovenFrame)*obsNum(ovenFrame)))
  expect_equal(mean(yy), mean(ovenFrame@y), 0.1)
})

test_that("Posterior sim methods for ubmsFitMultinomPois work with NAs",{
  zna <- posterior_predict(fit_double_na, "z", draws=3)
  expect_equal(dim(zna), c(3,50))
  expect_true(any(is.na(zna)))
  yna <- posterior_predict(fit_double_na, "y", draws=3)
  expect_equal(dim(yna), c(3,50*3))
  expect_equal(sum(is.na(yna[1,])), 6)
  expect_equal(sum(is.na(yna[2,])), 6)
})

test_that("Posterior linear pred methods work for ubmsFitMultinomPois",{
  set.seed(123)
  samples <- get_samples(fit_double, 3)
  lp1 <- sim_lp(fit_double, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), numSites(umf_double)))
  set.seed(123)
  pl <- posterior_linpred(fit_double, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccu",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit_double, "state", draws=10)
  ft2 <- ubms_fitted(fit_double, "det", draws=10)
  expect_equal(dim(ft), c(10,50))
  expect_equal(dim(ft2), c(10,150))

  res <- ubms_residuals(fit_double, "state", draws=10)
  res2 <- ubms_residuals(fit_double, "det", draws=10)
  expect_equal(dim(res), c(10,50))
  expect_equal(dim(res2), c(10,150))

  pdf(NULL)
  rp <- plot_residuals(fit_double, "state")
  rp2 <- plot_residuals(fit_double, "det")
  rp3 <- ubms_plot(fit_double)
  mp <- plot_marginal(fit_double, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(rp2, "gg")
  expect_is(rp3, "gtable")
  expect_is(mp, "gtable")
})

test_that("get_pifun_type returns correct value",{
  expect_equal("double", get_pifun_type(umf_double))
  expect_equal("removal", get_pifun_type(ovenFrame))
  umf_broken <- umf_double
  umf_broken@piFun <- "fake"
  expect_error(get_pifun_type(umf_broken))
})

test_that("ubmsResponseMultinomPois find_missing method works",{
  resp <- ubmsResponseMultinomPois(getY(umf_double), "double", "P")
  subs <- fit_double@submodels
  expect_true(!all(is.na(find_missing(resp, subs))))

  resp <- ubmsResponseMultinomPois(getY(umf_double_na), "double", "P")
  subs <- fit_double_na@submodels
  expect_equal(find_missing(resp, subs), c(rep(TRUE,6),rep(FALSE,144)))

  resp <- ubmsResponseMultinomPois(getY(ovenFrame_na), "removal", "P")
  subs <- fit_rem_na@submodels
  expect_equal(find_missing(resp, subs), c(rep(TRUE,8),rep(FALSE,70*4-8)))
})

test_that("ubmsResponseMultinomPois update_missing method works",{
  resp <- ubmsResponseMultinomPois(getY(umf_double_na), "double", "P")
  subs <- fit_double_na@submodels
  subs@submodels$state@missing <- rep(FALSE, 50)
  subs@submodels$det@missing <- rep(FALSE, length(subs@submodels$det@missing))
  out <- update_missing(subs, resp)
  expect_equal(out@submodels$state@missing, c(rep(TRUE, 2), rep(FALSE, 48)))
  expect_equal(out@submodels$det@missing, c(rep(TRUE, 4), rep(FALSE, 96)))

  resp <- ubmsResponseMultinomPois(getY(ovenFrame_na), "removal", "P")
  subs <- fit_rem_na@submodels
  subs@submodels$state@missing <- rep(FALSE, numSites(ovenFrame_na))
  subs@submodels$det@missing <- rep(FALSE, length(subs@submodels$det@missing))
  out <- update_missing(subs, resp)
  expect_equal(out@submodels$state@missing, c(rep(TRUE, 2), rep(FALSE, 68)))
  expect_equal(out@submodels$det@missing, c(rep(TRUE, 8), rep(FALSE, 272)))
})

test_that("get_pi_for_multinom function works",{
  pi_out <- get_pi_for_multinom(fit_double, 1:3)
  expect_is(pi_out, "array")
  expect_equal(dim(pi_out), c(50,3+1,3))
  expect_equal(rowSums(pi_out[,,1]), rep(1,50))

  pi_out <- get_pi_for_multinom(fit_rem, 1:3)
  expect_is(pi_out, "array")
  expect_equal(dim(pi_out), c(70,4+1,3))
  expect_equal(rowSums(pi_out[,,1]), rep(1,70))
})

test_that("getP and sim_p for ubmsFitMultinomPois work",{
  p <- sim_p(fit_double, 1:3)
  expect_equal(dim(p), c(3,50*3))
  p <- sim_p(fit_double_na, 1:3)
  expect_equal(dim(p), c(3,50*3))

  gp <- getP(fit_double, 3)
  expect_equal(dim(gp), c(50,3,3))
  gp <- getP(fit_double_na, 3)
  expect_equal(dim(gp), c(50,3,3))

  p <- sim_p(fit_rem, 1:3)
  expect_equal(dim(p), c(3,70*4))
  p <- sim_p(fit_rem_na, 1:3)
  expect_equal(dim(p), c(3,70*4))

  gp <- getP(fit_rem, 3)
  expect_equal(dim(gp), c(70,4,3))
  gp <- getP(fit_rem_na, 3)
  expect_equal(dim(gp), c(70,4,3))
})
