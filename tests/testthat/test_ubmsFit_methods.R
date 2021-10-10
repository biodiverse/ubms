context("ubmsFit methods")

skip_on_cran()

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model

good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occu(~x3~x1+(1|x2), umf,
                                  chains=2, iter=40, refresh=0))

fit_coef <- suppressWarnings(stan_occu(~x2~x1, umf, chains=2,
                                  iter=40, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

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
  expect_equal(printed[6], "Occupancy (logit-scale):")
  expect_equal(printed[12], "Detection (logit-scale):")
  expect_true(grepl("LOOIC:", printed[17]))
})

test_that("coef method works for ubmsFit",{
  co <- coef(fit)
  expect_is(co, "numeric")
  expect_equal(as.numeric(co), c(summary(fit, "state")$mean,
                                 summary(fit, "det")$mean))

  #When equal # of params in each submodel
  co <- coef(fit_coef)
  expect_is(co, "numeric")
  expect_equal(as.numeric(co), c(summary(fit_coef, "state")$mean,
                                 summary(fit_coef, "det")$mean))
})

test_that("summary method works for ubmsFit",{
  sum_fit <- summary(fit, "state")
  expect_is(sum_fit, "data.frame")
  expect_equal(rownames(sum_fit), c("(Intercept)", "x1", "sigma [1|x2]"))
  expect_equal(colnames(sum_fit), c("mean","se_mean","sd","2.5%","25%",
                                    "50%","75%","97.5%","n_eff","Rhat"))
  expect_true(all(sapply(sum_fit, inherits, "numeric")))
})

test_that("getY method works for ubmsFit",{
  expect_equal(fit@data@y, getY(fit))
})

test_that("loo method works for ubmsFit",{
  lout <- suppressWarnings(loo(fit))
  expect_is(lout, "psis_loo")
})

test_that("waic method works for ubmsFit",{
  wout <- suppressWarnings(waic(fit))
  expect_is(wout, "waic")
})

test_that("extract method works for ubmsFit",{
  ex1 <- extract(fit, "beta_state")
  expect_is(ex1, "list")
  expect_equal(names(ex1), "beta_state")
  expect_equal(dim(ex1[[1]]), c(40,2))
  ex_all <- extract(fit)
  expect_is(ex_all, "list")
  expect_equal(names(ex_all), c("beta_state","b_state","sigma_state",
                                "beta_det","log_lik","lp__"))
})

test_that("traceplot method works for ubmsFit",{
  #need a vdiff test for this eventually
  tr <- traceplot(fit, pars="beta_state")
  expect_is(tr, "gg")
})

test_that("predicting map works for ubmsFit",{
  r <- raster::raster(matrix(rnorm(30), ncol=5, nrow=6))
  names(r) <- "x1"
  r2 <- r
  names(r2) <- "x2"
  rs <- raster::stack(r, r2)
  pr_rast <- predict(fit, "state", newdata=r, re.form=NA)
  expect_is(pr_rast, "RasterBrick")
  expect_equal(length(pr_rast), 30*4)
  pr_rast2 <- predict(fit, "state", newdata=rs, re.form=NA)
  expect_is(pr_rast2, "RasterBrick")
  expect_equal(length(pr_rast2), 30*4)
})

test_that("getP method works for ubmsFit",{
  gp <- getP(fit)
  expect_equal(dim(gp), c(3,3,40))
})

test_that("get_elapsed_time method works for ubmsFit",{
  et <- get_elapsed_time(fit)
  expect_is(et, "matrix")
  expect_equal(colnames(et), c("warmup","sample"))
  expect_equal(rownames(et), c("chain:1","chain:2"))
})

test_that("get_runtime calculates runtime for display",{
  fit2 <- fit
  options(mc.cores=2)
  attr(fit2@stanfit@sim$samples[[1]],"elapsed_time") <- c(warmup=49,sample=50)
  rt <- get_runtime(fit2)
  expect_equal(rt, "99.000 sec")
  attr(fit2@stanfit@sim$samples[[1]],"elapsed_time") <- c(warmup=100,sample=50)
  expect_equal(get_runtime(fit2), "2.500 min")
  attr(fit2@stanfit@sim$samples[[1]],"elapsed_time") <- c(warmup=10000,sample=50)
  expect_equal(get_runtime(fit2), "2.792 hr")

  # non-parallel
  options(mc.cores=1)
  attr(fit2@stanfit@sim$samples[[1]],"elapsed_time") <- c(warmup=10,sample=10)
  attr(fit2@stanfit@sim$samples[[2]],"elapsed_time") <- c(warmup=10,sample=10)
  expect_equal(get_runtime(fit2), "40.000 sec")
})

test_that("plot_effects creates gg or grid object",{
  #Multiple covariates
  pdf(NULL)
  mp <- plot_posteriors(fit)
  expect_is(mp, "gg")
  mp2 <- plot_posteriors(fit, density=TRUE)
  expect_is(mp2, "gg")
  mp3 <- plot_posteriors(fit, "lp__")
  expect_is(mp3, "gg")
  expect_error(plot_posteriors(fit, pars="fake"))
  dev.off()
})

test_that("get_stancode method works",{
  out <- get_stancode(fit)
  expect_is(out, "character")
})
