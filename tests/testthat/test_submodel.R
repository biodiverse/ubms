context("Test ubmsSubmodel class")

test_that("ubmsSubmodel can be built",{
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,2,3)), ~x1, "plogis")
  expect_true(inherits(sm, "ubmsSubmodel"))
  expect_equal(sm@data, data.frame(x1=c(1,2,3)))
  expect_equal(sm@type, "det")
  expect_equal(sm@link, "plogis")
  expect_equal(do.call(sm@link, list(0)), 0.5)
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
