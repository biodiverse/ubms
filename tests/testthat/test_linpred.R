context("Get posterior of linear predictor")

on_mac <- tolower(Sys.info()[["sysname"]]) == "darwin"
on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
skip_if(on_mac & on_cran, "On CRAN mac")

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
}, error=function(e){
  good_fit <<- FALSE
})
skip_if(!good_fit, "Test setup failed")

test_that("predict correctly wraps sim_lp",{
  expect_error(predict(fit, "fake"))
  pr <- predict(fit, "state")
  expect_is(pr, "data.frame")
  expect_equal(names(pr), c("Predicted","SD","2.5%","97.5%"))
  expect_equal(dim(pr), c(3,4))

  #Newdata
  nd <- data.frame(x1=0)
  #Missing random effect
  expect_error(predict(fit, "state", newdata=nd))
  #Should work now
  pr2 <- predict(fit, "state", newdata=nd, re.form=NA)
  expect_is(pr2, "data.frame")
  expect_equal(dim(pr2), c(1,4))
  #Change level
  pr3 <- predict(fit, "state", newdata=nd, re.form=NA, level=0.8)
  expect_equal(names(pr3)[3:4], c("10%","90%"))
})

test_that("posterior_linpred correctly wraps sim_lp",{
  expect_error(posterior_linpred(fit, "fake"))
  expect_equal(posterior_linpred(fit, FALSE,"state",NULL,
                                 NULL,NULL),
               sim_lp(fit,"state",FALSE,NULL,1:40,NULL))
  #Check that smaller number of draws works
  set.seed(123)
  pl <- posterior_linpred(fit, FALSE,"state",NULL,draws=3,NULL)
  set.seed(123)
  lp <- sim_lp(fit,"state",FALSE,NULL,get_samples(fit,3),NULL)
  expect_equal(pl, lp)
})

test_that("sim_lp generates posterior correctly",{
  #For all iterations
  samp <- get_samples(fit, NULL)
  lp <- sim_lp(fit, "state", transform=FALSE, newdata=NULL,
               samp, re.form=NULL)
  expect_is(lp, "matrix")
  expect_equal(dim(lp), c(40,3))

  #For a few iterations
  samp <- c(1,2,3)
  lp <- sim_lp(fit, "state", transform=FALSE, newdata=NULL,
               samp, re.form=NULL)
  expect_is(lp, "matrix")
  expect_equal(dim(lp), c(3,3))
})

test_that("sim_lp transforms parameters correctly",{

  samp <- c(1,2,3)
  lp <- sim_lp(fit, "state", transform=FALSE, newdata=NULL,
               samp, re.form=NULL)

  lpt_compare <- do.call(fit["state"]@link, list(lp))

  lpt  <- sim_lp(fit, "state", transform=TRUE, newdata=NULL,
               samp, re.form=NULL)
  expect_equal(lpt, lpt_compare)
})

test_that("sim_lp can remove random effects from linear predictor",{
  samp <- c(1,2,3)
  lp_rand <- sim_lp(fit, "state", transform=FALSE, newdata=NULL,
               samp, re.form=NULL)
  lp_norand <- sim_lp(fit, "state", transform=FALSE, newdata=NULL,
               samp, re.form=NA)
  expect_equal(dim(lp_norand), dim(lp_rand))
  expect_true(all(lp_rand != lp_norand))

  #Unchanged when submodel has no random effects
  lp_rand <- sim_lp(fit, "det", transform=FALSE, newdata=NULL,
               samp, re.form=NULL)
  lp_norand <- sim_lp(fit, "det", transform=FALSE, newdata=NULL,
               samp, re.form=NA)
  expect_equal(dim(lp_norand), dim(lp_rand))
  expect_equal(lp_rand, lp_norand)
})

test_that("sim_lp handles newdata",{
  nd <- data.frame(x1=c(0,1))
  samp <- c(1,2,3)
  #Missing covariate
  expect_error(sim_lp(fit, "state", transform=FALSE, newdata=nd,
               samp, re.form=NULL))

  nd$x2 <- factor(c("a","b"), levels=c("a","b","c"))
  lp <- sim_lp(fit,"state",transform=FALSE, newdata=nd, samp,
               re.form=NULL)
  expect_equal(dim(lp), c(3,2))
})
