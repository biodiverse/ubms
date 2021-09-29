context("Spatial modeling functions")

skip_on_cran()

#Simulate dataset
set.seed(567)
dat_occ <- data.frame(cov1=rnorm(500), x=runif(500, 0,10), y=runif(500,0,10))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)
b <- c(0.4, -0.5, 0, 0.5)

#re_fac <- factor(sample(letters[1:26], 500, replace=T))
#dat_occ$group <- re_fac
#re <- rnorm(26, 0, 1.2)
#re_idx <- as.numeric(re_fac)

idx <- 1
for (i in sample(1:500, 300, replace=FALSE)){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$cov1[i]))# + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1, plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}

umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

fit <- suppressMessages(suppressWarnings(stan_occu(~1~cov1+RSR(x,y,1),
                              data=umf[1:20,], chains=2, iter=200, refresh=0)))
fit2 <- suppressWarnings(stan_occu(~1~1, data=umf[1:10,], chains=2, iter=200, refresh=0))

test_that("spatial model output structure", {
  expect_is(fit, "ubmsFitOccu")
  expect_true(has_spatial(fit@submodels@submodels$state))
  expect_equal(names(coef(fit))[3], "state[RSR [tau]]")
})

test_that("methods for spatial model work", {
  pr <- suppressMessages(predict(fit, "state"))
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(20,4))

  nd <- data.frame(cov1=c(0,1))
  expect_error(suppressMessages(predict(fit, "state", newdata=nd)))
  pr <- suppressMessages(predict(fit, "state", newdata=nd, re.form=NA))
  expect_equal(dim(pr), c(2,4))

  ss <- suppressMessages(sim_state(fit, samples=1:2))
  expect_equal(dim(ss), c(2,13))

  expect_warning(ppred <- suppressMessages(posterior_predict(fit, "z", draws=2)))
  expect_equal(dim(ppred), c(2, 13))

  expect_warning(ppred <- suppressMessages(posterior_predict(fit, "y", draws=2)))
  expect_equal(dim(ppred), c(2, 65))

  fitted <- getMethod("fitted", "ubmsFit") #why? only an issue in tests
  ft <- suppressMessages(fitted(fit, "state", draws=2))
  expect_is(ft, "matrix")
  expect_equal(dim(ft), c(2, 13))
})

test_that("RSR() generates spatial matrices", {

  rsr_out <- RSR(dat_occ$x, dat_occ$y, threshold=1)
  rsr_out2 <- RSR(dat_occ$x, dat_occ$y, threshold=5)
  expect_is(rsr_out, "list")
  expect_equal(names(rsr_out), c("A","Q","n_eig","coords"))
  expect_equal(as.matrix(rsr_out$coords),
               as.matrix(dat_occ[,c("x","y")]))
  expect_equal(rsr_out$Q[1], -sum(rsr_out$Q[1,2:500]))
  expect_true(sum(diag(rsr_out$Q)) < sum(diag(rsr_out2$Q)))

  expect_equal(rsr_out$n_eig, nrow(dat_occ)*0.1)
  rsr_out3 <- RSR(dat_occ$x, dat_occ$y, threshold=1, moran_cut=100)
  expect_equal(rsr_out3$n_eig, 100)
  expect_error(RSR(dat_occ$x, dat_occ$y, threshold=1, moran_cut=1000))

  expect_equal(dim(rsr_out$A), c(nrow(dat_occ), nrow(dat_occ)))
  expect_true(max(rsr_out$A)==1)

})

test_that("RSR() can generate a plot", {
  pdf(NULL)
  rsr_out <- RSR(dat_occ$x, dat_occ$y, threshold=1)
  gg <- RSR(dat_occ$x, dat_occ$y, threshold=1, plot_site=1)
  expect_is(gg, "gg")
  dev.off()
})

test_that("RSR info can be extracted from submodel", {

  sm <- fit@submodels@submodels$state
  inf <- get_rsr_info(sm)
  # will not match straight output from RSR() due to re-sorting  by
  # missing sites
  expect_is(inf, "list")
  expect_equal(names(inf), c("A","Q","n_eig","coords"))

})

test_that("remove_RSR removes spatial component of formula", {

  nf <- remove_RSR(fit@submodels@submodels$state@formula)
  expect_equal(as.formula(~cov1), nf)
  expect_equal(~1, remove_RSR(~RSR(x,y,1)))

  expect_error(remove_RSR(~cov1+RSR(x,y,1)+(1|fake)))
  expect_error(remove_RSR(~cov1+(1|fake)+RSR(x,y,1)))

})

test_that("has_spatial identifies spatial submodels", {
  expect_true(has_spatial(fit@submodels@submodels$state))
  expect_false(has_spatial(fit@submodels@submodels$det))
})

test_that("has_spatial works on lists of formulas", {
  expect_true(has_spatial(list(det=~1,state=~RSR(x,y,1))))
  expect_error(has_spatial(list(state=~1,det=~RSR(x,y,1))))
  expect_error(has_spatial(list(state=~RSR(x,y,1),det=~RSR(x,y,1))))
  expect_error(has_spatial(list(det=~1,state=~RSR(x,y,1)),support=FALSE))
})

test_that("construction of ubmsSubmodelSpatial objects", {
  ex <- extract_missing_sites(umf)
  sm <- ubmsSubmodelSpatial("Test","test", ex$umf@siteCovs, ~1+RSR(x,y,1), "plogis",
                            uniform(-5,5), normal(0,2.5), gamma(1,1),
                            ex$sites_augment, ex$data_aug)
  expect_is(sm, "ubmsSubmodelSpatial")
})

test_that("extract_missing_sites identifies augmented sites", {
  es <- extract_missing_sites(umf)
  expect_is(es, "list")
  expect_equivalent(es$sites_augment, apply(umf@y, 1, function(x) all(is.na(x))))
  expect_equal(nrow(es$data_aug), sum(es$sites_augment))
  expect_equal(nrow(siteCovs(es$umf)), numSites(umf) - sum(es$sites_augment))
  expect_true(!any(apply(es$umf@y, 1, function(x) all(is.na(x)))))

  # error on NAs
  umf2 <- umf
  umf2@siteCovs$cov1[1] <- NA
  expect_error(extract_missing_sites(umf2))
})

test_that("spatial_matrices builds correct RSR matrices", {

  sm <- fit@submodels@submodels$state

  mats <- suppressMessages(spatial_matrices(sm))
  expect_is(mats, "list")
  n_eig <- get_rsr_info(sm)$n_eig
  expect_equal(dim(mats$Qalpha), c(n_eig, n_eig))
  expect_equal(mats$Qalpha[1,1:2], c(0.08389,0.44826), tol=1e-4)
  expect_equal(dim(mats$Kmat), c(20, n_eig))
  expect_equal(mats$Kmat[1,1:2], c(-0.0197,0.0250), tol=1e-4)
  expect_equal(mats$n_eigen, n_eig)
})

test_that("get_pars method for ubmsSubmodelSpatial adds tau param", {
  sm <- fit@submodels@submodels$state
  expect_equal(get_pars(sm), c("beta_state","b_state","tau"))
})

test_that("get_stan_data for ubmsSubmodelSpatial includes spatial data", {
  sm <- fit@submodels@submodels$state
  dat <- suppressMessages(get_stan_data(sm))
  expect_is(dat, "list")
  expect_equal(dat$n_random_state[1], 2)
  expect_true(all(c("Kmat","Qalpha","n_eigen","n_aug_sites","X_aug","offset_aug") %in%
                  names(dat)))
  expect_equal(nrow(dat$X_aug), 7)
  expect_equal(length(dat$offset_aug), 7)
})

test_that("stanfit_names returns correct names for ubmsSubmodelSpatial", {
  sm <- fit@submodels@submodels$state
  sname <- stanfit_names(sm)
  expect_equal(length(sname), 5)
  expect_equal(sname[5], "tau")
})

test_that("plot_spatial returns ggplot", {
  pdf(NULL)
  gg1 <- suppressMessages(plot_spatial(fit, "state"))
  expect_is(gg1, "gg")
  gg2 <- suppressMessages(plot_spatial(fit, "eta"))
  expect_is(gg2, "gg")
  gg3 <- suppressMessages(plot_spatial(fit, "state", sites=TRUE))
  expect_is(gg3, "gg")
  dev.off()
  expect_error(plot_spatial(umf))
  expect_error(plot_spatial(fit2))
})
