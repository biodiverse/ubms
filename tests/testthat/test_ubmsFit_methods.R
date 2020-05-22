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
  expect_equal(printed[17], "LOOIC: 5.664")
})

test_that("summary method works for ubmsFit",{
  expect_equal(summary(fit, "state"),
    structure(list(mean = c(0.190718324171947, -1.02725030475444,
    314.838236809268), se_mean = c(0.364024063773018, 0.811952045654714,
    312.21313804252), sd = c(2.91406662318696, 2.72374523311318,
    444.172946007566), `2.5%` = c(-5.73254793281281, -5.79076638574038,
    2.99860731644845), `25%` = c(-1.45888024739656, -2.68275463657314,
    6.52616270794227), `50%` = c(0.411161792691752, -0.974109349810946,
    75.2311281787547), `75%` = c(1.650285839359, 0.690270318180331,
    507.037490482349), `97.5%` = c(5.71826675399397, 4.44994016205417,
    1472.40640510548), n_eff = c(64.0823996531185, 11.2531007128209,
    2.02395968062711), Rhat = c(0.970272787156061, 1.09195150515215,
    1.43312409652205)), row.names = c("(Intercept)", "x1", "sigma [1|x2]"
    ), class = "data.frame"))
  expect_equal(summary(fit, "det"),
    structure(list(mean = c(0.978648704264532, -2.43493485114631),
    se_mean = c(0.341197662951039, 0.154285303291694), sd = c(0.949404137581837,
    1.23507673671525), `2.5%` = c(-0.446671916854944, -5.02782121198986
    ), `25%` = c(0.206897014832347, -3.2569772018939), `50%` = c(0.961396415468061,
    -2.37656296376859), `75%` = c(1.80352107292484, -1.59412308092421
    ), `97.5%` = c(2.68389851003062, -0.140303241987431), n_eff = c(7.74265921347573,
    64.0823996531185), Rhat = c(1.10152012550616, 0.970375238670973
    )), class = "data.frame", row.names = c("(Intercept)", "x3"
    ))
  )
})

test_that("loo method works for ubmsFit",{
  lout <- suppressWarnings(loo(fit))
  expect_is(lout, "psis_loo")
  expect_equal(lout$estimates, 
    structure(c(-2.83219181573505, 0.856174775795061, 5.6643836314701,
    1.01621662398154, 0.174951669900272, 2.03243324796308), .Dim = 3:2, 
    .Dimnames = list(c("elpd_loo", "p_loo", "looic"), c("Estimate", "SE")))
  )
})

test_that("waic method works for ubmsFit",{
  wout <- suppressWarnings(waic(fit))
  expect_is(wout, "waic")
  expect_equivalent(wout$estimates[,1], 
                    c(-3.091303, 1.115286, 6.182607), tol=1e-6)
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
