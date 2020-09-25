context("Submodel construction and methods")

test_that("ubmsSubmodel can be built",{
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,2,3)), ~x1, "plogis")
  expect_true(inherits(sm, "ubmsSubmodel"))
  expect_equal(sm@data, data.frame(x1=c(1,2,3)))
  expect_equal(sm@type, "det")
  expect_equal(sm@link, "plogis")
  expect_equal(do.call(sm@link, list(0)), 0.5)
  expect_equivalent(sm@missing, rep(FALSE, 3))
})

test_that("ubmsSubmodelTransition can be built",{
  sm <- ubmsSubmodelTransition("Col", "col", data.frame(x1=c(1,2,3)), ~x1,
                               "plogis", 3)
  expect_true(inherits(sm, "ubmsSubmodelTransition"))
  expect_equal(nrow(sm@data), 2)
})

test_that("ubmsSubmodelTransition errors if NAs in yearlySiteCovs",{
  sm <- ubmsSubmodelTransition("Col", "col", data.frame(x1=c(1,2,3)), ~x1,
                               "plogis", 3)
  sm@data$x1[1] <- NA
  expect_error(ubmsSubmodelTransition("Col", "col", data.frame(x1=c(1,NA,3)), ~x1,
                            "plogis", 3))
  expect_error(ubmsSubmodel("Col", "col", data.frame(x1=c(1,NA,3)), ~x1,
                            "plogis"), NA)
})

test_that("ubmsSubmodelScalar built correctly",{
  ss <- ubmsSubmodelScalar("Fake", "fake", "plogis")
  expect_is(ss, "ubmsSubmodelScalar")
  expect_equal(ss@data, data.frame(X1=1))
  expect_equal(ss@formula, ~1)
  expect_equivalent(ss@missing, c(FALSE))
})

test_that("placeholderSubmodel creates blank submodel",{
  ps <- placeholderSubmodel("fake")
  expect_is(ps, "ubmsSubmodel")
  expect_equal(ps@data, data.frame())
  expect_equal(ps@formula, ~1)
  expect_equal(ps@link, "identity")
})

test_that("drop_final_year removes final year of yearly site covs",{
  M <- 5; T <- 3
  test_df <- data.frame(x1=rnorm(M*T),
                        x2=factor(rep(c("1","2","3"), M)))
  expect_equal(levels(test_df$x2), c("1","2","3"))

  dr <- drop_final_year(test_df, T)
  expect_is(dr, "data.frame")
  expect_equal(dim(dr), c(M*(T-1), 2))
  expect_equal(levels(dr$x2), c("1","2"))
})

test_that("Missing values are detected when submodel is built",{
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(NA,2,3)), ~x1, "plogis")
  expect_equivalent(sm@missing, c(TRUE,FALSE,FALSE))
})

test_that("model.matrix method works",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  mm <- model.matrix(sm)
  ref <- structure(c(1, 1, 1, 1, 2, 3), .Dim = 3:2,
                   .Dimnames =list(c("1","2", "3"), c("(Intercept)", "x1")),
                   assign = 0:1)
  expect_equal(mm, ref)
})

test_that("model.matrix works with newdata",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  nd <- data.frame(x1=4, x2=5)
  mm <- model.matrix(sm, nd)
  ref <- structure(c(1, 4), .Dim = 1:2, .Dimnames = list("1", c("(Intercept)",
        "x1")), assign = 0:1)
  expect_equal(mm, ref)

  #Missing covariate
  nd <- data.frame(x2=5)
  expect_error(model.matrix(sm, nd))
})

test_that("model.matrix errors if variables in newdata are missing",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  nd <- data.frame(x1=4, x2=5)
  expect_error(model.matrix(sm, nd), NA)
  nd <- data.frame(x3=4)
  expect_error(model.matrix(sm, nd))
})

test_that("check_newdata identifies missing variables in newdata",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+x2, "plogis")
  nd <- data.frame(x3=4)
  expect_error(check_newdata(nd, sm@formula), "newdata: x1, x2")
})

test_that("model.matrix handles NAs",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+x2, "plogis")
  mm <- model.matrix(sm)
  ref <- structure(c(1, 1, 1, 1, 2, 3, NA, 2, 4), .Dim = c(3L, 3L),
        .Dimnames = list(c("1", "2", "3"), c("(Intercept)", "x1", "x2")),
        assign = 0:2)
  expect_equal(mm, ref)
  sm@missing[1] <- TRUE
  mm <- model.matrix(sm, na.rm=TRUE)
  ref <- structure(c(1, 1, 2, 3, 2, 4), .Dim = 2:3, .Dimnames =
                   list(c("2","3"), c("(Intercept)", "x1", "x2")))
  expect_equal(mm, ref)
})

test_that("model.matrix catches invalid factor levels", {
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+x2, "plogis")
  nd <- data.frame(x1=1, x2="a")
  mm <- ubms:::model.matrix(sm, nd)
  ref <- structure(c(1, 1, 0, 0), .Dim = c(1L, 4L),
        .Dimnames =list("1",c("(Intercept)", "x1", "x2b", "x2c")),
        assign = c(0L, 1L,2L, 2L), contrasts = list(x2 = "contr.treatment"))
  expect_equal(mm, ref)
  nd <- data.frame(x1=1, x2="d")
  expect_error(model.matrix(sm, nd))
})

test_that("model.matrix handles functions in formulas", {
  covs <- data.frame(x1=c(2,1))
  sm <- ubmsSubmodel("Det", "det", covs, ~scale(x1), "plogis")
  mm <- model.matrix(sm)
  ref <- structure(c(1, 1, 0.707106781186547, -0.707106781186547),
                   .Dim = c(2L,2L), .Dimnames = list(c("1", "2"),
                  c("(Intercept)", "scale(x1)")), assign = 0:1)
  expect_equal(mm,ref)
  expect_equal(model.matrix(sm, covs), ref)
})

test_that("get_xlev gets levels from factor columns",{
  ref_df <- data.frame(x1=factor(c("a","a","b")),x2=c(1,2,3),
                      x3=factor(c("c","d","e")))
  mf <- model.frame(~x1+x2, ref_df)
  lvls <- get_xlev(ref_df, mf)
  expect_equal(lvls, list(x1=c("a","b")))
})

test_that("get_reTrms works",{
  dat <- data.frame(x1=rnorm(3), x2=c("a","b","c"))
  rt <- get_reTrms(~(1|x2), dat)
  expect_is(rt, "list")
  expect_equal(names(rt)[c(1,4,8)], c("Zt","Gp","cnms"))
})

test_that("get_reTrms works with NAs",{
  dat <- data.frame(x1=c(NA,1,2), x2=c("a","b",NA))
  rt <- get_reTrms(~(1|x2), dat)
  expect_is(rt, "list")
  expect_equal(names(rt)[c(1,4,8)], c("Zt","Gp","cnms"))
  rt <- get_reTrms(~(1+x1|x2), dat)
  expect_is(rt, "list")
  expect_equal(names(rt)[c(1,4,8)], c("Zt","Gp","cnms"))
})

test_that("correct Z matrix is built",{
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4),x3=c('a','a','b'))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")

  #No random effect
  expect_equal(Z_matrix(sm), matrix(0, nrow=0,ncol=0))

  #Random intercept
  sm <- ubmsSubmodel("Det", "det", covs, ~(1|x3), "plogis")
  ref <- structure(c(1, 1, 0, 0, 0, 1), .Dim = 3:2,
                   .Dimnames =list(c("1","2", "3"), c("a", "b")))
  expect_equal(Z_matrix(sm), ref)

  #Random slope + intercept
  sm <- ubmsSubmodel("Det", "det", covs, ~(1+x1||x3), "plogis")
  ref <- structure(c(1, 1, 0, 0, 0, 1, 1, 2, 0, 0, 0, 3), .Dim =3:4,
                   .Dimnames = list(c("1", "2", "3"), c("a", "b", "a", "b")))
  expect_equal(Z_matrix(sm), ref)
})

test_that("Missing values are handled when Z matrix is built", {
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4),x3=c('a',NA,'b'))
  sm <- ubmsSubmodel("Det", "det", covs, ~(1|x3), "plogis")

  #Random intercept
  ref <- structure(c(1, 0, 0, 0, 0, 1), .Dim = 3:2,
                   .Dimnames =list(c("1","2", "3"), c("a", "b")))

  #Random slope + intercept
  sm <- ubmsSubmodel("Det", "det", covs, ~(1+x1||x3), "plogis")
  ref <- structure(c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 3), .Dim =3:4,
                   .Dimnames = list(c("1", "2", "3"), c("a", "b", "a", "b")))
  expect_equal(Z_matrix(sm), ref)

  #With missing removed
  sm@missing[1] <- TRUE
  expect_equal(Z_matrix(sm), ref)
  ref <- structure(c(0, 0, 0, 1, 0, 0, 0, 3), .Dim = c(2L, 4L),
        .Dimnames = list(c("2", "3"), c("a", "b", "a", "b")))
  expect_equal(Z_matrix(sm, na.rm=T), ref)
})

test_that("Generating Z matrix from newdata works", {
  covs <- data.frame(x1=c(1,2,3),x2=c(NA,2,4),x3=c('a','a','b'))
  sm <- ubmsSubmodel("Det", "det", covs, ~(1|x3), "plogis")
  nd <- data.frame(x1=4, x3='b')
  ref <- structure(1, .Dim = c(1L, 1L), .Dimnames = list("1", "b"))
  expect_equal(Z_matrix(sm, nd), ref)
})

test_that("beta names are correct",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+x2, "plogis")
  expect_equal(beta_names(sm), c("(Intercept)", "x1", "x2b", "x2c"))
})

test_that("b names are correct",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")),
                     x3=factor(c("d","e","d")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  expect_equal(b_names(sm), NA_character_)

  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  expect_equal(b_names(sm),
                   paste("(Intercept)", c("x2:a", "x2:b", "x2:c")))

  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1+x1||x2), "plogis")
  expect_equal(b_names(sm),
                   c(paste("(Intercept)", c("x2:a", "x2:b", "x2:c")),
                   c(paste("x1", c("x2:a", "x2:b", "x2:c")))))

  sm <- ubmsSubmodel("Det", "det", covs, ~(1|x2) + (1|x3), "plogis")
  expect_equal(b_names(sm),
               c(paste("(Intercept)", c("x2:a","x2:b","x2:c","x3:d","x3:e"))))

  #Test when there is a random slope for a factor
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(x1+x3||x2), "plogis")
  expect_equal(b_names(sm),
                   c(paste("(Intercept)", c("x2:a", "x2:b", "x2:c")),
                   paste("x1", c("x2:a", "x2:b", "x2:c")),
                   c("x3d x2:a", "x3e x2:a", "x3d x2:b", "x3e x2:b",
                     "x3d x2:c", "x3e x2:c")))
})

test_that("sigma names are correct", {
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  expect_equal(sigma_names(sm), NA_character_)

  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  expect_equal(sigma_names(sm),"sigma [1|x2]")

  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1+x1||x2), "plogis")
  expect_equal(sigma_names(sm), paste0("sigma [", c("1|x2]", "x1|x2]")))
})

test_that("submodels with random effects are identified",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  expect_false(has_random(sm))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  expect_true(has_random(sm))
})

test_that("submodels with intercepts are identified",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1-1, "plogis")
  expect_false(has_intercept(sm))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1, "plogis")
  expect_true(has_intercept(sm))
})

test_that("generating parameter names works",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  expect_equal(b_par(sm), "b_det")
  expect_equal(sig_par(sm), "sigma_det")
  expect_equal(beta_par(sm), "beta_det")
})

test_that("check_formula identifies unsupported formulas",{
  dat <- data.frame(y=rnorm(10), x1=rnorm(10), x2=rnorm(10),
                    x3=sample(letters[1:3], 10, replace=T),
                    x4=sample(letters[4:6], 10, replace=T))
  expect_error(check_formula(~x1+(1|x2), dat), NA)
  expect_error(check_formula(~x1+(x1|x2),dat))
  expect_error(check_formula(~x1+(x1||x2),dat), NA)
  expect_error(check_formula(~(1|x1 / x2), dat))
  expect_error(check_formula(~(1|x1 : x2), dat))
})

test_that("submodel list creation works", {
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm1 <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm2 <- ubmsSubmodel("Det", "state", covs, ~x1+(1|x2), "plogis")
  sl <- ubmsSubmodelList(sm1, sm2)
  expect_is(sl, "ubmsSubmodelList")
  expect_equal(names(sl@submodels), c("det","state"))
})

test_that("submodel [ works",{
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm1 <- ubmsSubmodel("Det", "det", covs, ~x1+(1|x2), "plogis")
  covs <- data.frame(x1=rnorm(3), x2=factor(c("a","b","c")))
  sm2 <- ubmsSubmodel("Det", "state", covs, ~x1+(1|x2), "plogis")
  sl <- ubmsSubmodelList(sm1, sm2)
  expect_equal(sl["state"], sl@submodels[[2]])
  expect_error(sl["fake"])
})
