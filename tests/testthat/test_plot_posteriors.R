context("Posterior plots")

skip_on_cran()

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","b")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
good_fit <- TRUE
tryCatch({
fit <- suppressWarnings(stan_occu(~x3~x1+x2, umf,
                                  chains=2, iter=40, refresh=0))
fit2 <- suppressWarnings(stan_occu(~1~x1+(1|x2), umf,
                                  chains=2, iter=40, refresh=0))
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

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
  mp4 <- plot_posteriors(fit2)
  expect_is(mp4, "gg")
  dev.off()
})

