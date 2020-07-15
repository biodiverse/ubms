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

test_that("Error is thrown if y is not a matrix",{
  y <- c(1,0,1)
  expect_error(ubmsResponse(y, "binomial", "binomial"))
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
  resp <- ubmsResponse(y, "P", "P")

  expect_equal(get_K(resp), max(y,na.rm=T)+20)
  expect_equal(get_K(resp, 5), 5)
  expect_error(get_K(resp, 3))
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
  expect_equal(get_n_obs(resp), matrix(c(2,3,3), ncol=1))

  y <- matrix(c(NA,NA,NA,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_n_obs(resp), matrix(c(3,3),ncol=1))
})

test_that("per_sampled generates logical matrix of sampled periods",{
  M <- 3; J <- 3; T <- 4
  y1 = matrix(rbinom(M*T*J,1,0.5),M,T*J)
  resp1 <- ubmsResponse(y1, "binomial", "binomial", T)
  expect_equal(per_sampled(resp1), matrix(TRUE, nrow=3, ncol=4))

  y2 <- y1
  y2[2,1:3] <- NA
  y2[3,10:12] <- NA
  resp2 <- ubmsResponse(y2, "binomial", "binomial", T)
  ps <- per_sampled(resp2)
  expect_equal(ps, matrix(c(TRUE,FALSE,TRUE,rep(TRUE,8),FALSE),nrow=3))

  y3 <- y1[1,,drop=FALSE]
  resp3 <- ubmsResponse(y3, "binomial", "binomial", T)
  expect_equal(per_sampled(resp3), matrix(TRUE, nrow=1, ncol=4))

  y4 <- y1
  y4[3,] <- NA
  resp4 <- ubmsResponse(y4, "binomial", "binomial", T)
  expect_equal(per_sampled(resp4), matrix(TRUE, nrow=2, ncol=4))

  resp5 <- ubmsResponse(y1, "binomial", "binomial", 1)
  expect_equal(per_sampled(resp5), matrix(TRUE, nrow=3, ncol=1))
})

test_that("which_per_sampled identifies indices of sampled periods",{
  M <- 3; J <- 3; T <- 4
  y = matrix(rbinom(M*T*J,1,0.5),M,T*J)

  y2 <- y
  y2[2,1:3] <- NA
  y2[3,10:12] <- NA
  resp2 <- ubmsResponse(y2, "binomial", "binomial", T)
  expect_equal(which_per_sampled(resp2),
              c(1,2,3,4, 2,3,4, 1,2,3))

  y3 <- y
  y3[1,] <- NA
  resp3 <- ubmsResponse(y3, "binomial", "binomial", T)
  expect_equal(which_per_sampled(resp3), c(1,2,3,4,1,2,3,4))

  resp5 <- ubmsResponse(y, "binomial", "binomial", 1)
  expect_equal(which_per_sampled(resp5), rep(1,M))

})

test_that("get_n_pers gets correct # of primary pers by site",{
  M <- 3; J <- 3; T <- 4
  y1 = matrix(rbinom(M*T*J,1,0.5),M,T*J)
  resp1 <- ubmsResponse(y1, "binomial", "binomial", T)
  np <- get_n_pers(resp1)
  expect_equal(np, rep(T, M))

  y2 <- y1
  y2[2,1:3] <- NA
  y2[3,10:12] <- NA
  resp2 <- ubmsResponse(y2, "binomial", "binomial", T)
  np2 <- get_n_pers(resp2)
  expect_equal(np2, c(4,3,3))

  y3 <- y1
  y3[1,] <- NA
  resp3 <- ubmsResponse(y3, "binomial", "binomial", T)
  expect_equal(get_n_pers(resp3), c(4,4))

  y4 <- y1[1,,drop=FALSE]
  resp4 <- ubmsResponse(y4, "binomial", "binomial", T)
  expect_equal(get_n_pers(resp4), 4)

  resp5 <- ubmsResponse(y1, "binomial", "binomial", 1)
  expect_equal(get_n_pers(resp5), rep(1,M))
})

test_that("generate_inds creates start-stop indices from count vector",{
  cv1 <- c(3,4,5)
  expect_equivalent(generate_inds(cv1),
                    matrix(c(1,3,4,7,8,12), nrow=3, byrow=T))
  cv2 <- c(3)
  expect_equivalent(generate_inds(cv2), matrix(c(1,3), nrow=1))
  #This shouldn't happen, but worth checking
  cv3 <- c(3,0,5)
  expect_error(generate_inds(cv3))
})

test_that("get_subset_inds works correctly",{
  M <- 3; J <- 3; T <- 4
  y1 = matrix(rbinom(M*T*J,1,0.5),M,T*J)
  resp1 <- ubmsResponse(y1, "binomial", "binomial", T)
  ind1 <- get_subset_inds(resp1)
  expect_equivalent(ind1,
                    matrix(c(1,12, 1,4, 1,3,
                             13,24, 5,8, 4,6,
                             25,36, 9,12, 7,9), nrow=3, byrow=T))

  y2 <- y1
  y2[2,1:3] <- NA
  y2[3,10:12] <- NA
  resp2 <- ubmsResponse(y2, "binomial", "binomial", T)
  ind2 <- get_subset_inds(resp2)
  expect_equivalent(ind2,
                    matrix(c(1,12, 1,4, 1,3,
                             13,21, 5,7, 4,6,
                             22,30, 8,10, 7,9), nrow=3, byrow=T))

  y3 <- y1[1,,drop=FALSE]
  resp3 <- ubmsResponse(y3, "binomial", "binomial", T)
  ind3 <- get_subset_inds(resp3)
  expect_equivalent(ind3, matrix(c(1,12,1,4,1,3), nrow=2))

  resp4 <- ubmsResponse(y1, "binomial", "binomial", 1)
  ind4 <- get_subset_inds(resp4)
  expect_equivalent(ind4,
                    matrix(c(1,12, 1,1, 1,1,
                             13,24, 2,2, 2,2,
                             25,36, 3,3, 3,3), nrow=3, byrow=T))

  y2[3,] <- NA
  resp5 <- ubmsResponse(y2, "binomial", "binomial", T)
  ind5 <- get_subset_inds(resp5)
  expect_equivalent(ind5,
                    matrix(c(1,12, 1,4, 1,3,
                             13,21, 5,7, 4,6), nrow=2, byrow=T))
})

test_that("as_vector converts response object to y vec",{
  y <- matrix(c(1,NA,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equivalent(as_vector(resp), as.vector(t(y)))
  expect_equivalent(as_vector(resp, na.rm=TRUE),
                    na.omit(as.vector(t(y))))
})

test_that("get_Kmin finds minimum K by site and primary period",{
  y <- matrix(c(1,NA,0,2,1,1,2,2,3), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_Kmin(resp), matrix(c(1,2,3), ncol=1))

  y <- matrix(c(NA,NA,NA,2,1,1,2,2,3), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  expect_equal(get_Kmin(resp), matrix(c(2,3),ncol=1))

  y <- matrix(c(1,0,0,0,1,1,NA,1), nrow=2, byrow=T)
  resp <- ubmsResponse(y, "binomial", "binomial", 2)
  expect_equal(get_Kmin(resp), matrix(c(1,1,0,1), nrow=2))

  y <- matrix(c(1,0,0,0), nrow=1, byrow=T)
  resp <- ubmsResponse(y, "binomial", "binomial", 2)
  expect_equal(get_Kmin(resp), matrix(c(1,0), nrow=1))
})
