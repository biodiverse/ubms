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

test_that("modSel creates selection table for ubmsFitList",{
  fl1 <- fitList(fit, fit_null)
  ms <- suppressWarnings(modSel(fl1))
  expect_is(ms, "data.frame")
  expect_equal(names(ms), c("elpd","nparam","elpd_diff","se_diff","weight"))
  expect_equal(ms$elpd_diff[1], 0)
  expect_equal(ms$weight[1], 1, tol=1e-7)
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
