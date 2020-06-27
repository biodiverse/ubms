context("Utility functions")

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
fit <- suppressWarnings(stan_occu(~x3~x1+(1|x2), umf,
                                  chains=2, iter=40, refresh=0))

test_that("get_samples generates sample indices",{
  set.seed(123)
  gs1 <- get_samples(fit, 10)
  expect_equal(gs1, c(31,15,14,3,38,25,26,27,32,5))
  #Make sure seed gets same samples
  set.seed(123)
  expect_equal(get_samples(fit, 10), gs1)
  #Test when no draws provided
  expect_equal(get_samples(fit, NULL), 1:nsamples(fit))
  #Test when too many draws provided
  expect_equal(get_samples(fit, 50), 1:nsamples(fit))
})

test_that("submodel_types gets names of submodels",{
  expect_equal(submodel_types(fit), c("state","det"))
})

test_that("check_type throws error when bad submodel type is given",{
  expect_error(check_type("state", submodel_types(fit)), NA)
  expect_error(check_type("fake", submodel_types(fit)))
})

test_that("Theme function produces ggplot theme",{
  theme_object <- plot_theme()
  expect_is(theme_object, "theme")
  expect_is(theme_object, "gg")
})
