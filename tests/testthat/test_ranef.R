context("ubmsFit ranef method")

#Setup umf
set.seed(123)
sc <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")),
                 x4=factor(c("d","e","f")))
oc <- data.frame(x3=rnorm(9))
umf <- unmarkedFrameOccu(y=matrix(c(1,0,0,1,1,0,0,1,0), nrow=3),
        siteCovs=sc, obsCovs=oc)
#Fit model
fit <- suppressWarnings(stan_occu(~x3~x1+(1|x2), umf,
                                  chains=2, iter=40, refresh=0))
fit2 <- suppressWarnings(stan_occu(~x3~x1+(1+x1||x2), umf,
                                  chains=2, iter=40, refresh=0))
fit3 <- suppressWarnings(stan_occu(~x3~x1+(1|x2)+(1|x4), umf,
                                  chains=2, iter=40, refresh=0))
fit4 <- suppressWarnings(stan_occu(~x3~x1+(1|x2)-1, umf,
                                  chains=2, iter=40, refresh=0))

test_that("ranef on submodel without random effect errors",{
  expect_error(ranef(fit, "det"))
})

test_that("ranef on submodel with one random effect works",{
  r <- ranef(fit, "state")
  expect_is(r, "list")
  expect_equal(names(r), "x2")
  expect_is(r$x2, "list")
  expect_equal(names(r$x2), "(Intercept)")
  expect_equal(length(r$x2$`(Intercept)`), 3)
})

test_that("ranef works with means parameterization",{
  #Might want to do a longer run here in future to match to fit1 results
  r_mn <- ranef(fit4, "state")
  expect_is(r_mn, "list")
})

test_that("ranef summary works with one random effect",{
  r <- ranef(fit, "state", summary=TRUE)
  df <- r$x2$`(Intercept)`
  expect_is(df, "data.frame")
  expect_equal(names(df), c("Estimate","SD","2.5%","97.5%"))
  expect_equal(dim(df), c(3,4))
})

test_that("ranef works with both random slope/intercept",{
  r <- ranef(fit2, "state")
  expect_equal(names(r$x2), c("(Intercept)", "x1"))
  expect_true(all(sapply(r, function(x) inherits(x, "list"))))
  r2 <- ranef(fit2, "state", summary=TRUE)
  expect_equal(names(r2$x2), c("(Intercept)", "x1"))
  expect_true(all(sapply(r2$x2, inherits, "data.frame")))
})

test_that("ranef works with multiple random effects",{
  r <- ranef(fit3, "state")
  expect_equal(names(r), c("x2","x4"))
  expect_equivalent(lapply(r, names), c("(Intercept)","(Intercept)"))
  expect_true(all(sapply(r, inherits, "list")))
  r2 <- ranef(fit3, "state", summary=TRUE)
})

test_that("combine_same_name combines lists properly",{
  test1 <- list(a=list(a1=c(1,1)), a=list(a2=c(2,2)), b=list(b1=c(1,1)))
  expect_equal(length(test1$a), 1)
  comb1 <- combine_same_name(test1)
  expect_equal(length(comb1), 2)
  expect_equal(names(comb1$a), c("a1","a2"))
  expect_equal(names(comb1$b), "b1")
})
