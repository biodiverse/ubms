context("fitList creation and methods")

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
fit <- suppressWarnings(stan_occu(~x3~x1+(1|x2), umf, 
                                  chains=2, iter=40, refresh=0))

fit_null <- suppressWarnings(stan_occu(~1~1, umf, 
                                  chains=2, iter=40, refresh=0))

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
  expect_equal(ms,
    structure(list(elpd = c(-2.83219181573505, -6.99295560810481),
    nparam = c(0.856174775795061, 1.05728725245336), 
    elpd_diff = c(0,-4.16076379236976), se_diff = c(0, 0.19600789058965), 
    weight = c(0.999999973484196,2.6515803885907e-08)), 
              row.names = c("fit", "fit_null"), class ="data.frame"))
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
