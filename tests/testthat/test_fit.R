context("Build ubmsFit object")

skip_on_cran()

#Set up a submodel list
covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")),
                   x3=factor(c("a","b","c")))
pint <- uniform(-5,5)
pcof <- normal(0,2.5)
psig <- gamma(1,1)
sm1 <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis", pint, pcof, psig)
sm2 <- ubmsSubmodel("Occ", "state", covs, ~x1, "plogis", pint, pcof, psig)
sm3 <- ubmsSubmodel("Occ", "state", covs, ~x1+(1|x2), "plogis", pint, pcof, psig)
sm4 <- ubmsSubmodel("Occ", "state", covs, ~(1|x2)+(1|x3), "plogis", pint, pcof, psig)

#Set up response and stan inputs
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3))
umf <- process_umf(umf)
dm <- ubmsSubmodel("Det", "det", obsCovs(umf), ~1, "plogis", pint, pcof, psig)
sm <- ubmsSubmodel("Occ", "state", siteCovs(umf), ~1, "plogis", pint, pcof, psig)
sl <- ubmsSubmodelList(sm,dm)
resp <- ubmsResponse(umf@y,"binomial","binomial",max_primary=1)

#Build a stanfit object
set.seed(123)
inp <- build_stan_inputs("occu", resp, sl, log_lik=FALSE)

good_fit <- TRUE
tryCatch({
sf <- suppressWarnings(rstan::sampling(stanmodels[["single_season"]], inp$stan_data,
                        inp$pars,chains=2, iter=40, refresh=0))
test <- process_stanfit(sf, sl)
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("ubmsFit object is constructed correctly",{
  ufit <- suppressWarnings(ubmsFit("occu",
                                   as.call(str2lang("stan_occu(formula = ~1 ~ 1, data = umf, chains = 2, iter = 40)")),
                                   umf, resp, sl,
                                   chains=2, iter=40, refresh=0))
  expect_true(inherits(ufit, "ubmsFit"))
  expect_true(inherits(ufit@stanfit, "stanfit"))
  expect_equal(ufit@data, umf)
  expect_equivalent(ufit@submodels, sl)
  expect_equivalent(ufit@response, resp)
  expect_true(inherits(ufit@loo, "psis_loo"))
})

test_that("fit_class generates class from model name",{
  expect_equal(fit_class("occu"), "ubmsFitOccu")
})

test_that("remove_placeholders removes placeholder submodels from list",{
  ps <- placeholderSubmodel("fake")
  list_ps <- ubmsSubmodelList(sm, dm, ps)
  list_remove <- remove_placeholders(list_ps)
  expect_equal(list_remove, sl)
})

#test_that("get_loo generates loo object from stanfit",{
#  loo_obj <- suppressWarnings(get_loo(sf))
#  expect_true(inherits(loo_obj, "psis_loo"))
#})

test_that("fit_model builds model correctly",{
  ufit <- suppressWarnings(
    fit_model("occu", resp, sl, log_lik=FALSE, chains=2, iter=20, refresh=0))
  expect_true(inherits(ufit, "stanfit"))
  nms <- stanfit_names(sl)
  expect_equal(ufit@sim$fnames_oi[1:length(nms)], nms)
  #Check MCMC options are passed
  arg <- ufit@stan_args
  expect_equal(length(arg), 2)
  expect_equal(arg[[1]]$warmup, 10)
  expect_equal(arg[[1]]$iter, 20)
  expect_equal(arg[[1]]$refresh, 0)
})

test_that("specific model name is shown in console output",{
  # e.g. 'occu' instead of 'single_season'
  out <- capture.output(ufit <- suppressWarnings(
    fit_model("occu", resp, sl, log_lik=FALSE, chains=2, iter=20)))
  expect_true(any(grepl("occu", out)))
  expect_false(any(grepl("single_season", out)))
})

test_that("process_stanfit cleans up fitted stan model",{
  expect_true(inherits(sf, "stanfit"))
  psf <- process_stanfit(sf, sl)
  nms <- stanfit_names(sl)
  expect_equal(psf@sim$fnames_oi[1:length(nms)], nms)

  #test model failure
  sf2 <- sf
  sf2@mode = 1L
  expect_error(process_stanfit(sf2, sl))
})

test_that("stanfit_names for params in stan object are correct",{
  sl <- ubmsSubmodelList(sm1, sm2)
  expect_equal(stanfit_names(sl),
    c("beta_det[(Intercept)]", "beta_det[x1]", "b_det[(Intercept) x2:a]",
      "b_det[(Intercept) x2:b]", "b_det[(Intercept) x2:c]",
      "sigma_det[sigma [1|x2]]", "beta_state[(Intercept)]",
      "beta_state[x1]"))
})

test_that("beta names for ubmsSubmodel are generated correctly",{
  expect_equal(beta_names(sm1), c("(Intercept)", "x1"))
  expect_equal(beta_names(sm4), c("(Intercept)"))
})

test_that("b names for ubmsSubmodel are generated correctly",{
  # No random effect
  expect_true(is.na(b_names(sm2)))
  # One random effect
  expect_equal(b_names(sm1), paste0("(Intercept) x2:",levels(covs$x2)))
  # Two random effects
  expect_equal(b_names(sm4), c(paste0("(Intercept) x2:",levels(covs$x2)),
                              paste0("(Intercept) x3:", levels(covs$x3))))
})

test_that("sigma names for stanfit are generated correctly",{
  # No random effect
  expect_true(is.na(sigma_names(sm2)))
  # One random effect
  expect_equal(sigma_names(sm1), "sigma [1|x2]")
  # Two random effects
  expect_equal(sigma_names(sm4), c("sigma [1|x2]", "sigma [1|x3]"))
})
