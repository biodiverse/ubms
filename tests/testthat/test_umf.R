context("unmarkedFrame processing")

set.seed(567)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

re_fac <- factor(sample(letters[1:26], 500, replace=T))
dat_occ$group <- re_fac
re <- rnorm(26, 0, 1.2)
re_idx <- as.numeric(re_fac)

obs_fac <- factor(sample(letters[1:20], 500*5, replace=T))
obs <- rnorm(20, 0, 0.5)
obs_idx <- as.numeric(obs_fac)
dat_p$obs <- obs_fac

idx <- 1
for (i in 1:500){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$x1[i] + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1,
                    plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}
umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

test_that("cov data frames in unmarkedFrame are cleaned up",{

  umf2 <- process_umf(umf)
  expect_equal(siteCovs(umf), siteCovs(umf2))
  expect_equal(names(obsCovs(umf2)),
               c(names(obsCovs(umf)), names(siteCovs(umf))))
  expect_equal(nrow(obsCovs(umf)), nrow(obsCovs(umf2)))

  umf3 <- umf
  umf3@siteCovs <- NULL
  umf4 <- process_umf(umf3)
  expect_equal(siteCovs(umf4), data.frame(.dummy=rep(NA,500)))

  umf3 <- umf
  umf3@obsCovs <- NULL
  umf4 <- process_umf(umf3)
  expect_equal(names(obsCovs(umf4)), c(".dummy","x1", "group"))
  expect_equal(nrow(obsCovs(umf4)), 2500)
})


test_that("process_covs generates correct data frames",{

  expect_equal(process_covs(NULL, 5), data.frame(.dummy=rep(NA,5)))

  df1 <- data.frame(x1=rnorm(3), x2=rnorm(3))
  expect_equal(process_covs(df1, nrow(df1)), df1)

  df2 <- data.frame(x3=rnorm(1))
  expect_equal(process_covs(df1, nrow(df1), df2),
               cbind(df1, df2[c(1,1,1),,drop=FALSE]))

  expect_equal(process_covs(NULL, nrow(df1), df2),
               cbind(data.frame(.dummy=rep(NA, 3)),
                     df2[c(1,1,1),,drop=FALSE]))
})
