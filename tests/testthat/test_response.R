context("Response object construction and methods")

test_that("Response object is created properly",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_is(resp, "ubmsResponse")
  expect_equal(resp@y, y)
  expect_equal(resp@max_primary, 1)
  expect_equal(resp@max_obs, 3)
  expect_equal(resp@K, max(y, na.rm=T)+20)
  expect_equal(resp@missing, is.na(as.vector(t(y))))
})

test_that("get_max_obs calculates max # of observations",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  y2 <- matrix(c(1,0,0,0,1,1,0,0), nrow=2, byrow=T)

  expect_equal(get_max_obs(y, 1), 3)
  expect_error(get_max_obs(y, 2))
  expect_equal(get_max_obs(y2, 2), 2)
})

test_that("get_K takes input and generates valid K",{
  y <- matrix(c(1,NA,0,2,4,1,0,0,3), nrow=3, byrow=T)
  
  expect_equal(get_K(y), max(y,na.rm=T)+20)
  expect_equal(get_K(y, 5), 5)
  expect_error(get_K(y, 3))
})

test_that("tranpose method works",{
  y <- matrix(c(1,0,0,1,1,1,0,0,1), nrow=3, byrow=T)
  y2 <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(t(resp), t(y))
  resp@missing <- c(FALSE,TRUE, rep(FALSE, 7))
  expect_equal(t(resp), t(y2))
})

test_that("get_n_sites gets count of sites with data",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_n_sites(resp), 3)

  y <- matrix(c(NA,NA,NA,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_n_sites(resp), 2)
})

test_that("get_n_obs gets correct # of obs for each site",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_n_obs(resp), matrix(c(2,3,3), nrow=1))

  y <- matrix(c(NA,NA,NA,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_n_obs(resp), matrix(c(3,3),nrow=1))
})

test_that("as_vector converts response object to y vec",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equivalent(as_vector(resp), as.vector(t(y)))
  expect_equivalent(as_vector(resp, na.rm=TRUE), 
                    na.omit(as.vector(t(y))))
})

test_that("get_Kmin finds minimum K value for each site",{
  y <- matrix(c(1,NA,0,2,1,1,2,2,3), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_Kmin(resp), c(1,2,3))

  y <- matrix(c(NA,NA,NA,2,1,1,2,2,3), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_Kmin(resp), c(2,3))
})
