context("ubmsFit methods")

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
fit <- suppressWarnings(stan_occu(~x3~x1+(1|x2), umf,
                                  chains=2, iter=40, refresh=0))

test_that("fit is a ubmsFit object",{
  expect_is(fit, "ubmsFit")
})

test_that("[ method works for ubmsFit",{
  expect_is(fit["state"], "ubmsSubmodel")
  expect_equal(fit["state"]@type, "state")
  expect_error(fit["test"])
})

test_that("nsamples method works for ubmsFit",{
  expect_equal(nsamples(fit), (40/2)*2)
})

test_that("show method works for ubmsFit",{
  printed <- utils::capture.output(fit)
  expect_equal(printed[2], "Call:")
  expect_equal(printed[6], "Occupancy:")
  expect_equal(printed[12], "Detection:")
  expect_equal(printed[17], "LOOIC: 4.752")
})

test_that("summary method works for ubmsFit",{
  expect_equal(round(summary(fit, "state"), 3),
    structure(list(mean = c(-0.03, -0.674, 23.072),
    se_mean= c(0.503,1.153, 8.869), sd = c(2.557, 3.458, 21.484),
    `2.5%` = c(-6.19,-6.198, 2.999), `25%` = c(-1.315, -2.923, 6.526),
    `50%`= c(0.411,-1.177, 16.558), `75%` = c(1.65, 0.995, 35.19),
    `97.5%`= c(4.076,6.211, 73.71), n_eff = c(25.853, 8.998, 5.868),
    Rhat = c(1.104,1.143, 1.26)),
    row.names = c("(Intercept)", "x1", "sigma [1|x2]"),
    class = "data.frame")
  )
  expect_equal(round(summary(fit, "det"),3),
  structure(list(mean = c(1.318, -2.819), se_mean = c(0.157, 0.315),
                 sd = c(0.964, 1.626), `2.5%` = c(-0.208, -5.66),
                 `25%` = c(0.616,-3.823), `50%` = c(1.21, -2.603),
                 `75%` = c(2.031, -1.534), `97.5%` = c(3.163,-0.561),
                 n_eff = c(37.871, 26.728), Rhat = c(0.973, 0.981)),
            row.names = c("(Intercept)","x3"), class = "data.frame")

  )
})

test_that("loo method works for ubmsFit",{
  lout <- suppressWarnings(loo(fit))
  expect_is(lout, "psis_loo")
  expect_equivalent(lout$estimates[,1],
                    c(-2.3760,0.5779,4.7520), tol=1e-4)
})

test_that("waic method works for ubmsFit",{
  wout <- suppressWarnings(waic(fit))
  expect_is(wout, "waic")
  expect_equivalent(wout$estimates[,1],
                    c(-2.6566, 0.8585, 5.3132), tol=1e-3)
})

test_that("extract method works for ubmsFit",{
  ex1 <- extract(fit, "beta_state")
  expect_is(ex1, "list")
  expect_equal(names(ex1), "beta_state")
  expect_equal(dim(ex1[[1]]), c(40,2))
  ex_all <- extract(fit)
  expect_is(ex_all, "list")
  expect_equal(names(ex_all), c("beta_state","beta_det","b_state",
                                "sigma_state","log_lik","lp__"))
})

test_that("traceplot method works for ubmsFit",{
  #need a vdiff test for this eventually
  tr <- traceplot(fit, pars="beta_state")
  expect_is(tr, "gg")
})
