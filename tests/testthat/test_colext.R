context("stan_colext function and methods")

#Simulate dataset
set.seed(123)
M <- 150; T <- 3; J <- 3
z <- matrix(NA, M, T)

z[,1] <- rbinom(M, 1, 0.3)

for (i in 1:M){
  for (t in 1:(T-1)){
    if(z[i, t] == 1) z[i,(t+1)] <- rbinom(1, 1, 0.75)
    if(z[i, t] == 0) z[i,(t+1)] <- rbinom(1, 1, 0.6)
  }
}

x3 <- matrix(rnorm(M*T*J), M, T*J)
pest <- plogis(0 + x3*0.7)
zrep <- z[,rep(1:T, each=J)]
y <- matrix(NA, M, T*J)
for (i in 1:M){
  idx <- 1
  for (t in 1:T){
    for (j in 1:J){
      y[i,idx] <- rbinom(1, 1, pest[i,idx]*zrep[i,idx])
      idx <- idx + 1
    }
  }
}

sc <- data.frame(x1=factor(sample(letters[1:10], M, replace=TRUE)),
                 x2=rnorm(M))
oc <- list(x3=x3)
ysc <- data.frame(x4=rnorm(M*T))
umf <- unmarkedMultFrame(y, numPrimary=T, siteCovs=sc, obsCovs=oc,
                         yearlySiteCovs=ysc)

fit <- suppressWarnings(stan_colext(~x2,~x4,~1,~1, umf, chains=2, iter=300,
                                    refresh=0))

umf2 <- umf
umf2@y[1,] <- NA
umf2@y[2,1] <- NA

fit_na <- suppressWarnings(stan_colext(~x2,~x4,~1,~1, umf2, chains=2, iter=300,
                                       refresh=0))

fit_unm <- colext(~x2,~x4,~1,~1, umf)

test_that("stan_pcount output structure is correct",{
  expect_is(fit, "ubmsFitColext")
  expect_is(fit, "ubmsFitOccu")
  expect_equal(nsamples(fit), 300)
})

test_that("stan_colext produces accurate results",{
  #similar to truth
  beta <- c(log(0.3/0.7), 0, log(0.6/0.4), 0, log(0.25/0.75), 0)
  expect_equal(as.vector(coef(fit)), beta, tol=0.2)
  #similar to unmarked
  expect_equivalent(as.vector(coef(fit)), coef(fit_unm), tol=0.07)
  #similar to previous known values
  expect_equal(as.vector(coef(fit)), c(-0.84321,-0.038997,0.587029,
                                       0.1401878,-1.061324,0.0281296), tol=1e-4)
})

test_that("stan_colext handles NA values",{
  expect_equal(as.vector(coef(fit)), as.vector(coef(fit_na)), tol=0.1)
})

test_that("ubmsFitColext gof method works",{
  set.seed(123)
  g <- gof(fit, draws=10, quiet=TRUE)
  expect_equal(g@estimate, 16.8974, tol=1e-4)
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g)
  dev.off()
  expect_is(pg, "gg")
})

test_that("ubmsFitColext gof method works with missing values",{
  set.seed(123)
  g <- gof(fit_na, draws=10, quiet=TRUE)
  expect_is(g, "ubmsGOF")
})

test_that("ubmsFitColext predict method works",{
  pr <- predict(fit_na, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(numSites(umf2), 4))
  expect_equivalent(pr[1,1], 0.335, tol=0.01)
  pr <- predict(fit_na, "det")
  expect_equal(dim(pr), c(numSites(umf)*obsNum(umf),4))
  expect_equivalent(pr[1,1], 0.5106, tol=0.01)
  #with newdata
  nd <- data.frame(x2=c(0,1))
  pr <- predict(fit_na, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equivalent(pr[1,1], 0.3052, tol=0.01)
})

test_that("ubmsFitColext sim_z method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  zz <- sim_z(fit, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), numSites(umf)*umf@numPrimary))
  expect_equal(mean(zz), 0.5549, tol=1e-4)
  expect_equal(max(zz), 1)

  set.seed(123)
  pz <- posterior_predict(fit, "z", draws=10)
  expect_equivalent(zz, pz)
})

test_that("ubmsFitColext sim_y method works",{
  set.seed(123)
  samples <- get_samples(fit, 10)
  yy <- sim_y(fit, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(length(samples), numSites(umf)*obsNum(umf)))
  expect_equal(max(yy), 1)
  set.seed(123)
  py <- posterior_predict(fit, "y", draws=10)
  expect_equivalent(yy, py)
})

test_that("Posterior sim methods for ubmsFitColext work with NAs",{
  zna <- posterior_predict(fit_na, "z", draws=3)
  expect_equal(dim(zna), c(3,numSites(umf2)*umf2@numPrimary))
  expect_true(all(!is.na(zna[,1])))
  yna <- posterior_predict(fit_na, "y", draws=3)
  expect_equal(dim(yna), c(3, numSites(umf2)*obsNum(umf2)))
  expect_equal(sum(is.na(yna[1,])), 10)
  expect_equal(sum(is.na(yna[2,])), 10)
})

test_that("Posterior linear pred methods work for ubmsFitColext",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  lp1 <- sim_lp(fit, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), numSites(umf)))
  set.seed(123)
  pl <- posterior_linpred(fit, draws=3, submodel="state")
})

test_that("Fitted/residual methods work with ubmsFitColext",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  ft <- ubms_fitted(fit, "state", draws=10)
  ft2 <- ubms_fitted(fit, "det", draws=10)
  expect_equal(dim(ft), c(10,numSites(umf)*umf@numPrimary))
  expect_equal(dim(ft2), c(10,numSites(umf)*obsNum(umf)))

  res <- ubms_residuals(fit, "state", draws=10)
  res2 <- ubms_residuals(fit, "det", draws=10)
  expect_equal(dim(res), c(10,numSites(umf)*umf@numPrimary))
  expect_equal(dim(res2), c(10,numSites(umf)*obsNum(umf)))

  pdf(NULL)
  rp <- plot_residuals(fit, "state")
  mp <- plot_marginal(fit, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(mp, "gtable")
})

test_that("projected function and sim_state works with ubmsFitColext",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  pro <- sim_projected(fit, samples, NULL)
  expect_equal(dim(pro), c(3, numSites(umf)*umf@numPrimary))
  expect_equal(mean(pro), 0.557, tol=0.01)
  set.seed(123)
  expect_equivalent(projected(fit, 3), pro)

  pro_na <- sim_projected(fit_na, samples, NULL)
  expect_equal(dim(pro), dim(pro_na))

  set.seed(123)
  expect_equivalent(sim_state(fit, samples), pro)

  set.seed(123)
  expect_equal(dim(sim_state(fit_na, samples)), dim(pro))
})

test_that("turnover function works with ubmsFitColext",{
  set.seed(123)
  samples <- get_samples(fit, 3)
  turn <- sim_turnover(fit, samples, NULL)
  expect_equal(dim(turn), c(3, numSites(umf)*(umf@numPrimary-1)))
  expect_equal(mean(turn), 0.4817, tol=0.01)
  set.seed(123)
  expect_equivalent(turnover(fit, 3), turn)

  turn_na <- sim_turnover(fit_na, samples, NULL)
  expect_equal(dim(turn), dim(turn_na))
})
