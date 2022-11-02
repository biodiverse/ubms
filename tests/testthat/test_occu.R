context("stan_occu function and methods")

skip_on_cran()

#Simulate dataset
set.seed(567)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)
b <- c(0.4, -0.5, 0, 0.5)

#re_fac <- factor(sample(letters[1:26], 500, replace=T))
#dat_occ$group <- re_fac
#re <- rnorm(26, 0, 1.2)
#re_idx <- as.numeric(re_fac)

idx <- 1
for (i in 1:500){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$x1[i]))# + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1, plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}

umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occu(~x2~x1, umf[1:10,], chains=2,
                                  iter=100, refresh=0))

fit_na <- suppressWarnings(stan_occu(~x2~x1, umf2[1:10,], chains=2,
                                     iter=100, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("stan_occu output structure is correct",{
  expect_is(fit, "ubmsFitOccu")
  expect_equal(nsamples(fit), 100)
})

test_that("stan_occu produces accurate results",{
  skip_on_ci()
  skip_on_cran()
  skip_on_covr()
  set.seed(123)
  fit_long <- suppressWarnings(stan_occu(~x2~x1, umf[1:100,], chains=3,
                                  iter=300, refresh=0))
  fit_unm <- occu(~x2~x1, umf[1:100,])
  #similar to truth
  expect_RMSE(coef(fit_long), b, 0.25)
  #similar to unmarked
  expect_RMSE(coef(fit_long), coef(fit_unm), 0.05)
  #similar to previous known values
  expect_RMSE(coef(fit_long), c(0.66842,-0.71230,0.04183,0.45068), 0.1)
})

test_that("stan_occu handles NA values",{
  expect_RMSE(coef(fit), coef(fit_na), 2)
})

test_that("extract_log_lik method works",{
  ll <- extract_log_lik(fit)
  expect_is(ll, "matrix")
  expect_equal(dim(ll), c(100/2 * 2, numSites(fit@data)))
  expect_between(sum(ll), -3500, -3200)
})

test_that("extract_log_lik works when there are missing values and random effects",{
  skip_on_cran()
  skip_on_ci()
  umf3 <- umf2
  umf3@siteCovs$group <- sample(letters[1:5], nrow(umf2@siteCovs), replace=TRUE)
  fit_na <- suppressWarnings(stan_occu(~x2+(1|group)~1, umf3[1:10,], chains=2,
                                     iter=50, refresh=0))
  expect_is(fit_na, "ubmsFitOccu")
  ll <- extract_log_lik(fit_na)
  expect_is(ll, "matrix")
  expect_equal(dim(ll), c(50,9))
})

test_that("log_lik argument controls saving log_lik parameter", {
  skip_on_cran()
  set.seed(123)
  fit <- suppressWarnings(stan_occu(~x2~x1, umf[1:10,], chains=2,
                                  iter=100, refresh=0))
  set.seed(123)
  fit2 <- suppressWarnings(stan_occu(~x2~x1, umf[1:10,], chains=2,
                                  iter=100, refresh=0, log_lik=FALSE))

  expect_equal(fit@loo$estimates, fit2@loo$estimates)
  expect_true("log_lik" %in% fit@stanfit@sim$pars_oi)
  expect_false("log_lik" %in% fit2@stanfit@sim$pars_oi)

})

test_that("ubmsFitOccu gof method works",{
  set.seed(123)
  g <- gof(fit, draws=5, quiet=TRUE)
  expect_between(g@estimate, 15, 40)
  out <- capture.output(g)
  expect_equal(out[1], "MacKenzie-Bailey Chi-square ")
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
  #Check progress bar works
  # this test doesn't work non-interactively ??
  #out_pb <- capture.output(g <- gof(fit, draws=5))
  #expect_true(grepl("elapsed", out_pb))
  out_pb <- capture.output(g <- gof(fit, draws=5, quiet=TRUE))
  expect_equal(out_pb, character(0))
})

test_that("ubmsFitOccu gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_na, draws=5, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("stan_occu predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(10, 4))
  expect_between(pr[1,1], 0, 1)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(50,4))
  expect_between(pr[1,1], 0, 1)
  #with newdata
  nd <- data.frame(x1=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_between(pr[1,1], 0, 1)
})

test_that("stan_occu sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 10))
  expect_true(all(zz %in% c(0,1)))
  expect_between(mean(zz), 0, 1)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=5)
  expect_equivalent(zz, pz)
})

test_that("stan_occu sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), 50))
  expect_equal(max(yy), 1)
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=5)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitOccu work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,10))
  expect_true(all(is.na(zna[,1])))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3,50))
  expect_equal(sum(is.na(yna[1,])), 6)
  expect_equal(sum(is.na(yna[2,])), 6)
})

test_that("Posterior linear pred methods work for ubmsFitOccu",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 10))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccu",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=5)
  ft2 <- ubms_fitted(fit, "det", draws=5)
  expect_equal(dim(ft), c(5,10))
  expect_equal(dim(ft2), c(5,50))

  res <- ubms_residuals(fit, "state", draws=5)
  res2 <- ubms_residuals(fit, "det", draws=5)
  expect_equal(dim(res), c(5,10))
  expect_equal(dim(res2), c(5,50))

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

test_that("stan_occu kfold method works",{

set.seed(123)
kf <- kfold(fit, K=2, quiet=TRUE)
expect_is(kf, "elpd_generic")
expect_equal(kf$estimates[1,1], -37.3802, tol=1e-4)

expect_error(kfold(fit, K=2, folds=rep(3,10)))
expect_error(kfold(fit, K=2, folds=rep(1,5)))

})
