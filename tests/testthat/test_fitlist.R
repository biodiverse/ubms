context("fitList creation and methods")

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

fit_null <- suppressWarnings(stan_occu(~1~1, umf,
                                  chains=2, iter=40, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("fitList of ubmsFit objects is created correctly",{
  #no names
  fl1 <- fitList(fit, fit_null)
  expect_is(fl1, "ubmsFitList")
  expect_equal(names(fl1@models), c("fit","fit_null"))
  expect_equal(length(fl1@models), 2)
  expect_equal(fl1@models[[1]], fit)
  expect_equal(fl1@models[[2]], fit_null)

  #Named
  fl2 <- fitList(mod1=fit, mod2=fit_null)
  expect_is(fl2, "ubmsFitList")
  expect_equal(names(fl2@models), c("mod1", "mod2"))
})

test_that("fitlist can be created from a list",{
  mod_list <- list(fit=fit, fit_null=fit_null)
  fl <- fitList(fit, fit_null)
  fl2 <- fitList(mod_list)
  expect_is(fl2, "ubmsFitList")
  expect_equal(fl, fl2)
  mod_list2 <- list(fit, fit_null)
  fl3 <- fitList(mod_list2)
  expect_equal(names(fl3@models), c("mod1","mod2"))
  fl4 <- fitList(fits=mod_list)
  expect_equal(fl4, fl)
})

test_that("fitList for unmarkedFits passes through to unmarked",{
  unm_mod <- occu(~1~1, umf)
  unm_list <- list(mod1=unm_mod, mod2=unm_mod)
  ufl <- fitList(fits=unm_list)
  expect_is(ufl, "unmarkedFitList")
  ufl2 <- fitList(unm_mod, unm_mod)
  expect_is(ufl2, "unmarkedFitList")
})

test_that("modSel creates selection table for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  ms <- suppressWarnings(modSel(fl1))
  expect_is(ms, "data.frame")
  expect_equal(names(ms), c("elpd","nparam","elpd_diff","se_diff"))
  expect_equal(ms$elpd_diff[1], 0)
})

test_that("[ method works for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  expect_equal(fl1[1], list(fit=fit))
})

test_that("[[ method works for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  expect_equal(fl1[[1]], fit)
})

test_that("$ method works for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  expect_equal(fl1$fit, fit)
})

test_that("names method works for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  expect_equal(names(fl1), c("fit","fit_null"))
  fl2 <- fitList(mod1=fit, mod2=fit_null)
  expect_equal(names(fl2), c("mod1", "mod2"))
})
