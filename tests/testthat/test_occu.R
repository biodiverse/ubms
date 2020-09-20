context("stan_occu function and methods")

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

fit <- suppressWarnings(stan_occu(~x2~x1, umf[1:100,], chains=3,
                                  iter=300, refresh=0))
fit_unm <- occu(~x2~x1, umf[1:100,])

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

fit_na <- suppressWarnings(stan_occu(~x2~x1, umf2[1:100,], chains=3,
                                     iter=300, refresh=0))

test_that("stan_occu output structure is correct",{
  expect_is(fit, "ubmsFitOccu")
  expect_equal(nsamples(fit), 450)
})

test_that("stan_occu produces accurate results",{
  #similar to truth
  expect_equal(as.vector(coef(fit)), b, tol=0.3)
  #similar to unmarked
  expect_equivalent(as.vector(coef(fit)), coef(fit_unm), tol=0.05)
  #similar to previous known values
  expect_equal(as.vector(coef(fit)), c(0.6564,-0.7333,0.0347,0.4315), tol=1e-4)
})

test_that("stan_occu handles NA values",{
  expect_equal(as.vector(coef(fit)), as.vector(coef(fit_na)), tol=0.1)
})

test_that("ubmsFitOccu gof method works",{
  set.seed(123)
  g <- gof(fit, draws=10, quiet=TRUE)
  expect_equal(g@estimate, 25.8909, tol=1e-4)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitOccu gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_na, draws=10, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("stan_occu predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(100, 4))
  expect_equivalent(pr[1,1], 0.5373, tol=0.01)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(500,4))
  expect_equivalent(pr[1,1], 0.6124, tol=0.01)
  #with newdata
  nd <- data.frame(x1=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equivalent(pr[1,1], 0.65298, tol=0.01)
})

test_that("stan_occu sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 100))
  expect_equal(unique(as.vector(zz)), c(1,0))
  expect_equal(mean(zz), 0.647)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=10)
  expect_equivalent(zz, pz)
})

test_that("stan_occu sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), 500))
  expect_equal(unique(as.vector(yy)), c(1,0))
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=10)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitOccu work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,100))
  expect_false(any(is.na(zna)))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3,500))
  expect_equal(sum(is.na(yna[1,])), 6)
  expect_equal(sum(is.na(yna[2,])), 6)
})

test_that("Posterior linear pred methods work for ubmsFitOccu",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 100))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccu",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=10)
  ft2 <- ubms_fitted(fit, "det", draws=10)
  expect_equal(dim(ft), c(10,100))
  expect_equal(dim(ft2), c(10,500))

  res <- ubms_residuals(fit, "state", draws=10)
  res2 <- ubms_residuals(fit, "det", draws=10)
  expect_equal(dim(res), c(10,100))
  expect_equal(dim(res2), c(10,500))

  pdf(NULL)
  rp <- plot_residuals(fit, "state")
  rp2 <- plot_residuals(fit, "det")
  rp3 <- ubms_plot(fit)
  mp <- plot_marginal(fit, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(rp2, "gg")
  expect_is(rp3, "gtable")
  expect_is(mp, "gtable")
})
