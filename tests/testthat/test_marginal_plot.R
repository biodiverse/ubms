context("Marginal effect plots")

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","b")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
fit <- suppressWarnings(stan_occu(~x3~x1+x2, umf, 
                                  chains=2, iter=40, refresh=0))
fit2 <- suppressWarnings(stan_occu(~1~1, umf, 
                                  chains=2, iter=40, refresh=0))

test_that("marginal_plot creates grid object",{
  #Multiple covariates
  pdf(NULL)
  mp <- marginal_plot(fit, "state")
  expect_is(mp, "gtable")
  #Submodel with single covariate
  mp2 <- marginal_plot(fit, "det")
  expect_is(mp2, "gtable")
  #Specific covariate
  mp3 <- marginal_plot(fit, "state", "x1")
  expect_is(mp3, "gtable")
  #No covariates
  expect_error(marginal_plot(fit2, "state"))
  dev.off()
})

test_that("marginal_covariate_plot builds plot for either cov type",{
  pdf(NULL)
  mcp <- marginal_covariate_plot(fit,"state","x1")
  expect_is(mcp, "gg")
  mcp2 <- marginal_covariate_plot(fit,"state","x2")
  expect_is(mcp, "gg")
  #Covariate not in submodel
  expect_error(marginal_covariate_plot(fit,"state","fake"))
  dev.off()
})

test_that("marg_numeric_plot builds plot for numeric covariate",{
  pdf(NULL)
  mnp <- marg_numeric_plot(fit, "state", "x1", c(0.025,0.975))
  expect_is(mnp, "gg")
  expect_error(marg_numeric_plot(fit,"state","x2",c(0.025,0.975)))
  dev.off()
})

test_that("marg_factor_plot builds plot for factor covariate",{
  pdf(NULL)
  mfp <- marg_factor_plot(fit, "state", "x2", c(0.025,0.975))
  expect_is(mfp, "gg")
  expect_error(marg_factor_plot(fit,"state","x1",c(0.025,0.975)))
  dev.off()
})

test_that("get_mean_df gets baseline data frame for plot data",{
  mean_df <- get_mean_df(fit["state"])
  expect_equal(mean_df, data.frame(x1=0.2560184, x2="a"), tol=1e-7)
})

test_that("col_is_factor identifies factor column",{
  expect_true(col_is_factor("x2", sc))
  expect_false(col_is_factor("x1", sc))
})

test_that("Marginal plot data is generated correctly",{
  nd <- data.frame(x1=c(-1,1), x2=c("a","b"))

  #Continuous
  mdata <- get_margplot_data(fit, "state", "x1", c(0.025,0.975),
                             samples=1:3, nd)
  expect_equal(mdata,
    structure(list(covariate = c(-1, 1), mn = c(0.980940262011558,
    0.39301938256072), lower = c(0.969065755954778, 0.194268798455224
    ), upper = c(0.987901780495942, 0.705277052076411)), class = "data.frame", 
              row.names = c(NA,-2L)))
  #Factor
  mdata2 <- get_margplot_data(fit, "state", "x2", c(0.025,0.975),
                              samples=1:3, nd)
  expect_equal(mdata2,
    structure(list(covariate = c("a", "b"), mn = c(0.980940262011558,
    0.39301938256072), lower = c(0.969065755954778, 0.194268798455224), 
                   upper = c(0.987901780495942, 0.705277052076411)), 
              class = "data.frame", row.names = c(NA,-2L)))
})

test_that("Theme function produces ggplot theme",{
  theme_object <- marg_theme()
  expect_is(theme_object, "theme")
  expect_is(theme_object, "gg")
})
