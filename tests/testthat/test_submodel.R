context("Test ubmsSubmodel class")

test_that("ubmsSubmodel can be built",{
  sm <- ubmsSubmodel("Det", "det", data.frame(x1=c(1,2,3)), ~x1, "plogis")
  expect_true(inherits(sm, "ubmsSubmodel"))
  expect_equal(sm@data, data.frame(x1=c(1,2,3)))
  expect_equal(sm@type, "det")
  expect_equal(sm@link, "plogis")
  expect_equal(do.call(sm@link, list(0)), 0.5)
  expect_equal(sm@X, model.matrix(~x1, sm@data))
  expect_equal(sm@Z, matrix(0,0,0))
  expect_equal(sm@beta_names, colnames(model.matrix(~x1, sm@data)))
  expect_true(all(is.na(c(sm@b_names, sm@sigma_names))))
})

context("Test building design matrices")

test_that("get_X builds correct design matrices",{
  ref <- model.matrix(~x1+x2, data.frame(x1=1:3,x2=factor(letters[1:3])))
  expect_equal(ref,get_X(~x1+x2, data.frame(x1=1:3,x2=factor(letters[1:3]))))
  expect_equal(ref,get_X(~x1+x2+(1|x3), 
                         data.frame(x1=1:3,x2=factor(letters[1:3]))))
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

test_that("get_Z builds correct random effects matrices",{
  dat <- data.frame(y=rnorm(10), x1=rnorm(10), x2=rnorm(10),
                    x3=sample(letters[1:3], 10, replace=T),
                    x4=sample(letters[4:6], 10, replace=T))
  ref <- t(as.matrix(lme4::glFormula(y~x1 + (1|x3), data=dat)$reTrms$Zt))

  z1 <- get_Z(~x1+(1|x3), dat)
  expect_equal(ref, z1)
  expect_equal(ncol(z1), length(unique(dat$x3)))

})

