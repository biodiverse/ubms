context("Residuals and residual plots")

on_mac <- tolower(Sys.info()[["sysname"]]) == "darwin"
on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
skip_if(on_mac & on_cran, "On CRAN mac")

#sim_res methods tested in fitting function test scripts

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(9), x2=factor(sample(c("a","b","c"),9,replace=T)))
oc <- data.frame(x3=rnorm(27))
umf <- unmarkedFrameOccu(y=matrix(rep(c(1,0,0,1,1,0,0,1,0), 3), nrow=9),
        siteCovs=sc, obsCovs=oc)
#Fit model
good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occuRN(~x3~x1, umf,
                                  chains=2, iter=40, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("residuals generates matrix of correct structure",{
  residuals <- getMethod("residuals", "ubmsFit") #why? only an issue in tests
  r <- residuals(fit, "state", draws=3)
  expect_is(r, "matrix")
  expect_equal(dim(r), c(3,9))
  r <- residuals(fit, "state", draws=1)
  expect_is(r, "matrix")
  expect_equal(dim(r), c(1,9))
  r <- residuals(fit, "state")
  expect_equal(dim(r), c(40, 9))
})

test_that("plot_residuals generates correct plot",{
  pdf(NULL)
  pl1 <- plot_residuals(fit, "state")
  pl2 <- plot_residuals(fit, "det")
  pl3 <- plot_residuals(fit, "state", covariate="x1")
  dev.off()

  #State should be Pearson
  expect_is(pl1, "gg")
  pl1_build <- ggplot2::ggplot_build(pl1)
  expect_equal(pl1_build$plot$labels$y, "Pearson residual")

  #Det should be binned
  expect_is(pl2, "gg")
  pl2_build <- ggplot2::ggplot_build(pl2)
  expect_equal(pl2_build$plot$labels$y, "Mean binned residual")

  #Covariate plot
  expect_is(pl3, "gg")
  pl3_build <- ggplot2::ggplot_build(pl3)
  expect_equal(pl3_build$plot$labels$x, "x1 value")
})

test_that("plot_pearson_residuals generates ggplot",{
  x <- matrix(rnorm(100*10), nrow=10)
  res <- matrix(rnorm(100*10), nrow=10)
  pdf(NULL)
  out <- plot_pearson_residuals(x, res, "dummy xlab", "test")
  dev.off()
  expect_is(out, "gg")
})

test_that("plot_binned_residuals generates ggplot",{
  x <- matrix(rnorm(100*10), nrow=10)
  res <- matrix(rnorm(100*10), nrow=10)
  pdf(NULL)
  out <- plot_binned_residuals(x, res, "dummy xlab", "test", NULL)
  dev.off()
  expect_is(out, "gg")
})

test_that("get_binned_residuals separates residuals into appropriate bins",{
  set.seed(123)
  x <- rnorm(10)
  y <- rnorm(10)
  res <- get_binned_residuals(x, y, ind=1)
  expect_is(res, "data.frame")
  expect_equal(names(res), c("x_bar","y_bar","y_lo","y_hi","ind"))
  expect_true(all(res$ind==1))
  expect_equal(nrow(res), 4)
  expect_equivalent(res[1,1:2], c(-0.97595,-0.63263), tol=1e-4)
})

test_that("get_binned_residuals handles NAs",{
  set.seed(123)
  x <- rnorm(10)
  x[1] <- NA
  y <- rnorm(10)
  y[2] <- NA
  expect_error(get_binned_residuals(x, y, ind=1)) #Not enough data points

  set.seed(123)
  x <- rnorm(20)
  x[1] <- NA
  y <- rnorm(20)
  y[2] <- NA
  res <- get_binned_residuals(x, y, ind=1)
  expect_is(res, "data.frame")
  expect_equal(nrow(res), 5)
})

test_that("get_breaks finds correct cut points",{
  set.seed(123)
  x <- rnorm(100)
  br <- get_breaks(x, 10)
  expect_is(br, "list")
  expect_equal(br$nbins, 10)
  expect_equal(length(br$x_binned), length(x))
  expect_equal(br$x_binned[1:5], c(3,4,10,6,6))
  expect_equal(length(unique(br$x_binned)), 10)

  #Not enough bins
  expect_error(get_breaks(x, 2))

  #Fewer bins required
  set.seed(123)
  x <- c(rnorm(10,0,1), rnorm(10,1,1))
  br2 <- get_breaks(x, 25)
  expect_equal(br2$nbins, 20)
  expect_equal(length(br2$x_binned), length(x))
  expect_equal(br2$x_binned[1:5], c(4,6,16,7,8))
})

