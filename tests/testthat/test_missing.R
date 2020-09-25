context("Missing value handling")

sc <- data.frame(x1=rnorm(3), group=factor(c("a","b","a")))
oc <- data.frame(x2=rnorm(9))
state <- ubmsSubmodel("Occ", "state", sc, ~x1+(1|group), "plogis")
det <- ubmsSubmodel("Det", "det", oc, ~x2, "plogis")
sl <- ubmsSubmodelList(state, det)
y <- matrix(c(1,0,0,1,1,1,0,0,1), nrow=3, byrow=T)
resp <- ubmsResponse(y, "binomial", "binomial")

sc2 <- sc
sc2$x1[1] <- NA
oc2 <- oc
oc2$x2[4] <- NA
state2 <- ubmsSubmodel("Occ", "state", sc2, ~x1+(1|group), "plogis")
det2 <- ubmsSubmodel("Det", "det", oc2, ~x2, "plogis")
sl2 <- ubmsSubmodelList(state2, det2)
y2 <- matrix(c(1,0,0,1,1,1,NA,0,1), nrow=3, byrow=T)
resp2 <- ubmsResponse(y2, "binomial", "binomial")

test_that("update_missing works for submodel lists",{
  #No changes
  expect_equivalent(update_missing(sl, resp), sl)

  #With missing
  sl_new <- update_missing(sl2, resp2)
  expect_equal(sl_new["state"]@missing, c(T,F,F))
  expect_equal(sl_new["det"]@missing, c(T,T,T,T,F,F,T,F,F))
})

test_that("setting submodel missing attribute works",{
  state_new <- submodel_set_missing(sl['state'],
                c(T,T,T,F,F,F,T,F,F), resp)
  expect_equal(state_new@missing, c(T,F,F))

  det_new <- submodel_set_missing(sl['det'],
                c(T,T,T,F,F,F,T,F,F), resp)
  expect_equal(det_new@missing, c(T,T,T,F,F,F,T,F,F))
})

test_that("submodel_set_missing works with transition-type parameters",{
  ysc <- data.frame(x3 = rnorm(9))
  col <- ubmsSubmodelTransition("Col", "col", ysc, ~x3, "plogis", 3)
  y3 <- y2
  y3[1,] <- NA
  resp3 <- ubmsResponse(y3, "binomial", "binomial", 3)

  new_miss <- submodel_set_missing(col, rep(FALSE,9), resp3)@missing
  expect_equal(new_miss, c(TRUE,TRUE,rep(FALSE,4)))
})

test_that("error thrown if dimensions of missing slot changes",{
  ysc <- data.frame(x3 = rnorm(12))
  col <- ubmsSubmodelTransition("Col", "col", ysc, ~x3, "plogis", 3)
  y3 <- y2
  y3[1,] <- NA
  resp3 <- ubmsResponse(y3, "binomial", "binomial", 3)
  expect_error(submodel_set_missing(col, rep(FALSE,9), resp3))
})

test_that("update_missing works for response object",{
  #No changes
  expect_equivalent(update_missing(resp, sl), resp)

  #With missing
  resp_new <- update_missing(resp2, sl2)
  expect_equivalent(resp_new@missing, c(T,T,T,T,F,F,T,F,F))
})

test_that("find_missing identifies missing values",{
  miss1 <- find_missing(resp, sl)
  expect_equivalent(miss1, rep(FALSE, 9))

  miss2 <- find_missing(resp2, sl2)
  expect_equivalent(miss2, c(T,T,T,T,F,F,T,F,F))
})

test_that("find_missing ignores transition-type parameters",{
  ysc <- data.frame(x3 = rnorm(3))
  col <- ubmsSubmodel("Col", "col", ysc, ~x3, "plogis")
  col2 <- ubmsSubmodelTransition("Col", "col", ysc, ~x3, "plogis", 3)
  sl3 <- ubmsSubmodelList(state2, det2, col)
  sl4 <- ubmsSubmodelList(state2, det2, col2)
  expect_equal(find_missing(resp, sl3), find_missing(resp, sl4))
  expect_equal(find_missing(resp, sl2), find_missing(resp, sl4))
})

test_that("expand_model_matrix works correctly",{

  expand1 <- expand_model_matrix(sl["state"], resp)
  expected <- model.matrix(sl["state"])[rep(1:3, each=3),]
  expect_equivalent(expand1, expected)

  expand2 <- expand_model_matrix(sl["det"], resp)
  expected <- model.matrix(sl["det"])
  expect_equivalent(expand2, expected)
})

test_that("get_row_reps calculates correct replication factor",{
  expect_equal(get_row_reps(sl["state"], resp), 3)
  expect_equal(get_row_reps(sl["det"], resp), 1)

  y3 <- matrix(c(1,0,0,1), nrow=2)
  resp3 <- ubmsResponse(y3, "binomial", "binomial")
  expect_error(get_row_reps(sl["det"], resp3))
  expect_error(get_row_reps(sl["state"], resp3))
})

test_that("ubmsSubmodelScalar is returned unchanged by update_missing",{
  sub_scalar <- ubmsSubmodelScalar("Fake", "fake", "plogis")
  sl_scalar <- ubmsSubmodelList(state, det, sub_scalar)
  um <- update_missing(sl_scalar, resp2)
  expect_equivalent(um@submodels$fake, sub_scalar)
})

test_that("Placeholder submodel is returned unchanged by update_missing",{
  sub_place <- placeholderSubmodel("fake")
  expect_true(is_placeholder(sub_place))
  sl_place <- ubmsSubmodelList(state, det, sub_place)
  um <- update_missing(sl_place, resp2)
  expect_equivalent(um@submodels$fake, sub_place)
})
