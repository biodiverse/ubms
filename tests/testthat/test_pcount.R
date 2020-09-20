context("stan_pcount function and methods")

#Simulate dataset
set.seed(123)
M <- 50
J <- 5
beta <- c(1, 0.5, 0.2, -0.4)

sc <- data.frame(x1=rnorm(M), x2=sample(letters[1:26],M,replace=T),
                stringsAsFactors=TRUE)
oc <- data.frame(x3=rnorm(M*J))

#lambda <- exp(beta[1] + beta[2]*sc$x1)

#sig_x2 <- 0.3
#rx2 <- rnorm(26, 0, sig_x2)
#rx_idx <- as.numeric(sc$x2)

lambda <- exp(beta[1] + beta[2]*sc$x1) #+ rx2[rx_idx])
N <- rpois(M, lambda)
p <- plogis(beta[3] + beta[4]*oc$x3)

y <- matrix(NA, M, J)
idx <- 1
for (i in 1:M){
  y[i,] <- rbinom(J, N[i], p[idx:(idx+J-1)])
  idx <- idx + J
}

umf <- unmarkedFramePCount(y=y,siteCovs=sc, obsCovs=oc)

fit <- suppressWarnings(stan_pcount(~x3~x1, umf, K=15, chains=2, iter=200, refresh=0))

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

fit_na <- suppressWarnings(stan_pcount(~x3~x1, umf2, K=15, chains=2, iter=200, refresh=0))

fit_unm <- pcount(~x3~x1, umf, K=15)

test_that("stan_pcount output structure is correct",{
  expect_is(fit, "ubmsFitPcount")
  expect_is(fit, "ubmsFitAbun")
  expect_equal(nsamples(fit), 200)
})

test_that("stan_pcount produces accurate results",{
  #similar to truth
  expect_equal(as.vector(coef(fit)), beta, tol=0.2)
  #similar to unmarked
  expect_equivalent(as.vector(coef(fit)), coef(fit_unm), tol=0.05)
  #similar to previous known values
  expect_equal(as.vector(coef(fit)), c(0.94395,0.54603,0.05662,-0.36615), tol=1e-4)
})

test_that("stan_pcount handles NA values",{
  expect_equal(as.vector(coef(fit)), as.vector(coef(fit_na)), tol=0.1)
})

test_that("ubmsFitPcount gof method works",{
  set.seed(123)
  g <- gof(fit, draws=10, quiet=TRUE)
  expect_equal(g@estimate, 224.699, tol=1e-4)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitPcount gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_na, draws=10, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("ubmsFitPcount predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(numSites(umf2), 4))
  expect_equivalent(pr[1,1], 1.97867, tol=0.01)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(numSites(umf2)*obsNum(umf2),4))
  expect_equivalent(pr[1,1], 0.52625, tol=0.01)
  #with newdata
  nd <- data.frame(x1=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equivalent(pr[1,1], 2.68866, tol=0.01)
})

test_that("ubmsFitPcount sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), numSites(umf)))
  expect_equal(mean(zz), 3.032)
  expect_equal(colMeans(zz), N, tol=0.2)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=10)
  expect_equivalent(zz, pz)
})

test_that("ubmsFitPcount sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), numSites(umf)*obsNum(umf)))
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=10)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitPcount work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,numSites(umf2)))
  expect_true(all(is.na(zna[,1])))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3, numSites(umf2)*obsNum(umf2)))
  expect_equal(sum(is.na(yna[1,])), 5)
  expect_equal(sum(is.na(yna[2,])), 5)
})

test_that("Posterior linear pred methods work for ubmsFitPcount",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), numSites(umf)))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitPcount",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=10)
  ft2 <- ubms_fitted(fit, "det", draws=10)
  expect_equal(dim(ft), c(10,numSites(umf)))
  expect_equal(dim(ft2), c(10,numSites(umf)*obsNum(umf)))

  res <- ubms_residuals(fit, "state", draws=10)
  res2 <- ubms_residuals(fit, "det", draws=10)
  expect_equal(dim(res), c(10,numSites(umf)))
  expect_equal(dim(res2), c(10,numSites(umf)*obsNum(umf)))

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
