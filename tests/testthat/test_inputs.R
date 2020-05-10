context("Stan input generation")

test_that("get_group_vars returns number of grouping variables",{
  expect_equal(get_group_vars(~x), 0) 
  expect_equal(get_group_vars(~(1|x)), 1)
  expect_equal(get_group_vars( ~(1|x) + (1|y)), 2)
})

test_that("get_nrandom returns number of levels of each grouping variable",{
  dat <- data.frame(x=factor(c("a","b","c")), y=factor("d","e"))
  expect_equal(get_nrandom(~x, dat), as.array(0))
  expect_equal(get_nrandom(~(1|x), dat), as.array(3))
  form <- ~(1|x) + (1|y)
  expect_equal(get_nrandom(form, dat), as.array(c(3,1)))
})

test_that("split_formula works",{
  inp <- ~1~1
  expect_equal(split_formula(inp), list(~1, ~1))
  inp <- ~x1~x2
  expect_equal(split_formula(inp), list(~x1, ~x2))
  inp <- ~x1+(1|x2)~x3
  expect_equal(split_formula(inp), list(~x1+(1|x2), ~x3))
  inp <- ~x1
  expect_error(split_formula(inp))
  inp <- y~x1
  expect_error(split_formula(inp))
})

