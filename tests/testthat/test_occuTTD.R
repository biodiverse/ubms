context("stan_occuTTD function and methods")

skip_on_cran()

set.seed(123)

N <- 500; J <- 1
scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                         forest=runif(N,0,1),
                         wind=runif(N,0,1))

beta_psi <- c(-0.69, 0.71)
psi <- plogis(cbind(1, scovs$elev) %*% beta_psi)
z <- rbinom(N, 1, psi)

#Simulate detection
Tmax <- 10 #Same survey length for all observations
beta_lam <- c(-2, 0.7)
rate <- exp(cbind(1, scovs$wind) %*% beta_lam)
ttd <- rexp(N, rate)
ttd[z==0] <- Tmax #Censor at unoccupied sites
ttd[ttd>Tmax] <- Tmax #Censor when ttd was greater than survey length

#Build unmarkedFrame
umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

set.seed(123)
ocovs <- data.frame(obs=rep(c('A','B'),N))
Tmax <- 10
rateB <- exp(cbind(1, scovs$wind) %*% beta_lam + 0.2)
rate2 <- as.numeric(t(cbind(rate, rateB)))
ttd <- rexp(N*2, rate2)
ttd[ttd>Tmax] <- Tmax
ttd <- matrix(ttd, nrow=N, byrow=T)
ttd[z==0,] <- Tmax

umf_2obs <- suppressWarnings(unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax,
                                 siteCovs=scovs, obsCovs=ocovs))

good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occuTTD(~elev, detformula=~wind,
                        data=umf[1:10,], chains=2, iter=100, refresh=0))

fit_na <- suppressWarnings(stan_occuTTD(~elev, detformula=~wind,
                            data=umf2[1:10,], chains=2, iter=100, refresh=0))

fit_2obs <- suppressWarnings(stan_occuTTD(~elev, detformula=~wind,
                             data=umf_2obs[1:10,], chains=2, iter=100, refresh=0))

fit_weib <- suppressWarnings(stan_occuTTD(~elev, detformula=~wind,
                             data=umf[1:10,], ttdDist="weibull",
                             chains=2, iter=100, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("stan_occuTTD output structure is correct",{
  expect_is(fit, "ubmsFitOccuTTD")
  expect_equal(nsamples(fit), 100)
})

test_that("stan_occuTTD produces accurate results",{
  skip_on_ci()
  skip_on_cran()
  skip_on_covr()
  set.seed(123)
  fit_long <- suppressWarnings(stan_occuTTD(~elev, detformula=~wind,
                               data=umf, chains=3, iter=300, refresh=0))
  fit_unm <- occuTTD(~elev, detformula=~wind, data=umf)
  #similar to truth
  expect_RMSE(coef(fit_long), c(beta_psi, beta_lam), 0.4)
  #similar to unmarked
  expect_RMSE(coef(fit_long), coef(fit_unm), 0.1)
  #similar to previous known values
  expect_RMSE(coef(fit_long), c(-0.805,0.714,-1.759,0.694), 0.05)
})

test_that("stan_occuTTD handles NA values",{
  expect_is(coef(fit_na), "numeric")
})

test_that("ubmsFitOccuTTD gof method gives error",{
  expect_error(gof(fit, draws=5, quiet=TRUE))
})

test_that("stan_occuTTD predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(10, 4))
  expect_between(pr[1,1], 0, 1)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(10,4))
  expect_between(pr[1,1], 0, 10)
  #with newdata
  nd <- data.frame(elev=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_between(pr[1,1], 0, 1)
})

test_that("stan_occuTTD getP method works",{
  p <- getP(fit, draws=3)
  expect_equal(dim(p), c(10,1,3))
  p2 <- getP(fit_weib, draws=3)
  expect_equal(dim(p2), c(10,1,3))
  pna <- getP(fit_na, draws=3)
  expect_equal(dim(pna), c(10,1,3))
  expect_true(all(is.na(pna[1:2,1,1])))
})

test_that("stan_occuTTD sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 10))
  expect_equal(unique(as.vector(zz)), c(0,1))
  expect_between(mean(zz), 0, 0.5)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=5)
  expect_equivalent(zz, pz)
})

test_that("stan_occuTTD sim_z method warns when >1 obs per site",{
  expect_warning(posterior_predict(fit_2obs, "z", draws=3))
})

test_that("stan_occuTTD sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 5)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), 10))
  expect_equal(max(yy), max(umf@surveyLength))
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=5)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitOccu work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,10))
  expect_true(all(is.na(zna[,1])))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3,10))
  expect_equal(sum(is.na(yna[1,])), 2)
  expect_equal(sum(is.na(yna[2,])), 2)
})

test_that("Posterior linear pred methods work for ubmsFitOccuTTD",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 10))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitOccuTTD",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=5)
  ft2 <- ubms_fitted(fit, "det", draws=5)
  expect_equal(dim(ft), c(5,10))
  expect_equal(dim(ft2), c(5,10))

  res <- ubms_residuals(fit, "state", draws=5)
  res2 <- ubms_residuals(fit, "det", draws=5)
  expect_equal(dim(res), c(5,10))
  expect_equal(dim(res2), c(5,10))

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

test_that("occuTTD spatial works", {
  skip_on_cran()
  umf2 <- umf
  umf2@siteCovs$x <- runif(numSites(umf2), 0, 10)
  umf2@siteCovs$y <- runif(numSites(umf2), 0, 10)

  fit_spat <- suppressMessages(suppressWarnings(stan_occuTTD(~elev+RSR(x,y,1),
                                                             detformula=~wind,
                        data=umf2[1:20,], chains=2, iter=100, refresh=0)))

  expect_is(fit_spat@submodels@submodels$state, "ubmsSubmodelSpatial")
  expect_equal(names(coef(fit_spat))[3], "state[RSR [tau]]")

  ps <- plot_spatial(fit_spat)
  expect_is(ps, "gg")
})
