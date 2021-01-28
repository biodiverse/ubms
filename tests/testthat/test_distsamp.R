context("stan_distsamp function and methods")

skip_on_cran() #for now
on_mac <- tolower(Sys.info()[["sysname"]]) == "darwin"
on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
skip_if(on_mac & on_cran, "On CRAN mac")

#Line transect fits
data(linetran)
ltUMF <- with(linetran, {
        unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
        siteCovs = data.frame(Length, area, habitat),
        dist.breaks = c(0, 5, 10, 15, 20),
        tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
        })

linetran_big <- do.call("rbind", lapply(1:12, function(x) linetran))
ltUMF_big <- with(linetran_big, {
        unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
        siteCovs = data.frame(Length, area, habitat),
        dist.breaks = c(0, 5, 10, 15, 20),
        tlength = linetran_big$Length * 1000, survey = "line", unitsIn = "m")
        })

good_fit <- TRUE
tryCatch({
fit_line_hn <- suppressWarnings(stan_distsamp(~1~habitat, ltUMF, chains=2,
                                              iter=200, refresh=0))
fit_line_exp <- suppressWarnings(stan_distsamp(~1~habitat, ltUMF, keyfun="exp",
                                              chains=2, iter=200, refresh=0))
fit_line_haz <- suppressWarnings(stan_distsamp(~1~habitat, ltUMF, keyfun="hazard",
                                               chains=2, iter=15, refresh=0))
fit_line_abun <- suppressWarnings(stan_distsamp(~1~habitat, ltUMF, output="abund",
                                                chains=2, iter=200, refresh=0))
line_mods <- list(fit_line_hn, fit_line_exp, fit_line_haz, fit_line_abun)
}, error=function(e){
  good_fit <<- FALSE
})

ltUMF_na <- ltUMF
ltUMF_na@y[1,] <- NA
ltUMF_na@y[2,1] <- NA

data(pointtran)
ptUMF <- with(pointtran, {
             unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5),
             siteCovs = data.frame(area, habitat),
             dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
             })
pointtran_big <- do.call("rbind", lapply(1:4, function(x) pointtran))
ptUMF_big <- with(pointtran_big, {
             unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5),
             siteCovs = data.frame(area, habitat),
             dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
             })

tryCatch({
fit_pt_hn <- suppressWarnings(stan_distsamp(~1~habitat, ptUMF, chains=2,
                                              iter=200, refresh=0))
fit_pt_exp <- suppressWarnings(stan_distsamp(~1~habitat, ptUMF, keyfun="exp",
                                              chains=2, iter=200, refresh=0))
fit_pt_haz <- suppressWarnings(stan_distsamp(~1~habitat, ptUMF, keyfun="hazard",
                                               chains=2, iter=15, refresh=0))
point_mods <- list(fit_pt_hn, fit_pt_exp, fit_pt_haz)
}, error=function(e){
  good_fit <<- FALSE
})

skip_if(!good_fit, "Test setup failed")

test_that("stan_distsamp output structure is correct",{
  expect_true(all(sapply(line_mods, function(x) class(x)[1])=="ubmsFitDistsamp"))
  expect_true(all(sapply(line_mods,
                         function(x) class(x@response)[1])=="ubmsResponseDistsamp"))
  expect_equal(sapply(line_mods, function(x) x@response@y_dist),
              c("halfnorm", "exp", "hazard", "halfnorm"))
  expect_equal(nsamples(fit_line_hn), 200)

  expect_true(all(sapply(point_mods, function(x) class(x)[1])=="ubmsFitDistsamp"))
  expect_equal(sapply(point_mods, function(x) x@response@y_dist),
              c("halfnorm", "exp", "hazard"))
})

test_that("stan_distsamp produces accurate results",{
  skip_on_cran()
  skip_on_ci()
  skip_on_covr()

  #Line transect
  set.seed(123)
  stan_mod <- suppressWarnings(stan_distsamp(~1~1, ltUMF_big,
                                             chains=2, iter=200, refresh=0))
  um_mod <- distsamp(~1~1, ltUMF_big)
  expect_equivalent(coef(stan_mod), coef(um_mod), tol=0.05)

  stan_mod <- suppressWarnings(stan_distsamp(~1~1, ltUMF_big, keyfun="exp",
                                             chains=2, iter=200, refresh=0))
  um_mod <- distsamp(~1~1, ltUMF_big, keyfun="exp")
  expect_equivalent(coef(stan_mod), coef(um_mod), tol=0.05)

  stan_mod <- suppressWarnings(stan_distsamp(~1~1, ltUMF_big, output="abund",
                            chains=2, iter=200, refresh=0))
  um_mod <- distsamp(~1~1, ltUMF_big, output="abund")
  expect_equivalent(coef(stan_mod), coef(um_mod), tol=0.05)

  #Point
  set.seed(123)
  stan_mod <- suppressWarnings(stan_distsamp(~1~1, ptUMF, chains=2,
                                             iter=200, refresh=0))
  um_mod <- distsamp(~1~1, ptUMF)
  expect_equivalent(coef(stan_mod), coef(um_mod), tol=0.05)

  stan_mod <- suppressWarnings(stan_distsamp(~1~1, ptUMF_big, keyfun="exp",
                                             chains=2, iter=200, refresh=0))
  expect_equivalent(coef(stan_mod), c(4.97957, 2.063640), tol=0.05)
  #unmarked model fails here

  #Need hazard tests here at some point, maybe when I can speed it up
})

test_that("stan_distsamp handles NA values",{
  expect_error(stan_distsamp(~1~1, ltUMF_na))
})

test_that("ubmsFitDistsamp gof method works",{
  set.seed(123)
  g <- lapply(line_mods, function(x) gof(x, draws=5, quiet=TRUE))
  expect_true(all(sapply(g, function(x) class(x))=="ubmsGOF"))
  gof_plot_method <- methods::getMethod("plot", "ubmsGOF")
  pdf(NULL)
  pg <- gof_plot_method(g[[1]])
  dev.off()
  expect_is(pg, "gg")

  set.seed(123)
  g <- lapply(point_mods, function(x) gof(x, draws=5, quiet=TRUE))
  expect_true(all(sapply(g, function(x) class(x))=="ubmsGOF"))
})

test_that("ubmsFitDistsamp predict method works",{
  pr <- predict(fit_line_hn, "state")
  expect_is(pr, "data.frame")
  expect_equal(dim(pr), c(12, 4))
  expect_true(between(pr[1,1], 0, 3))
  pr <- predict(fit_line_hn, "det")
  expect_equal(dim(pr), c(12,4))
  #expect_true(between(pr[1,1], 5, 20))
  #with newdata
  nd <- data.frame(habitat=c("A","B"))
  pr <- predict(fit_line_hn, "state", newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_true(between(pr[1,1], 0, 3))
})

test_that("ubmsFitDistsamp sim_z method works",{
  set.seed(123)
  samples <- 1:3
  zz <- sim_z(fit_line_hn, samples, re.form=NULL)
  expect_is(zz, "matrix")
  expect_equal(dim(zz), c(length(samples), 12))
  expect_true(between(mean(zz), 10, 20))

  set.seed(123)
  pz <- posterior_predict(fit_line_hn, "z", draws=3)
  expect_equivalent(dim(zz), dim(pz))

  zlist <- lapply(line_mods, function(x) sim_z(x, samples, re.form=NULL))
  expect_true(all(sapply(zlist, function(x) all(dim(x)==dim(pz)))))

  zlist2 <- lapply(point_mods, function(x) sim_z(x, samples, re.form=NULL))
  expect_true(all(sapply(zlist2, function(x) all(dim(x)==c(3,30)))))
})

test_that("ubmsFitDistsamp sim_y method works",{
  set.seed(123)
  samples <- 1:3
  yy <- sim_y(fit_line_hn, samples, re.form=NULL)
  expect_is(yy, "matrix")
  expect_equal(dim(yy), c(3, 48))
  expect_true(between(mean(yy), 1.5, 3.5))
  set.seed(123)
  py <- posterior_predict(fit_line_hn, "y", draws=3)
  expect_equivalent(dim(yy), dim(py))

  #Test abundance model
  yabun <- sim_y(fit_line_abun, samples, re.form=NULL)
  expect_true(between(mean(yy), 1.5, 3.5))

  ylist <- lapply(line_mods, function(x) sim_y(x, samples, re.form=NULL))
  expect_true(all(sapply(ylist, function(x) all(dim(x)==dim(py)))))

  ylist2 <- lapply(point_mods, function(x) sim_y(x, samples, re.form=NULL))
  expect_true(all(sapply(ylist2, function(x) all(dim(x)==c(3,150)))))
})

test_that("Posterior linear pred methods work for ubmsFitDistsamp",{
  set.seed(123)
  samples <- get_samples(fit_line_hn, 3)
  lp1 <- sim_lp(fit_line_abun, "state", transform=TRUE, samples=samples,
                newdata=NULL, re.form=NULL)
  expect_equal(dim(lp1), c(length(samples), 12))
  set.seed(123)
  pl <- posterior_linpred(fit_line_hn, draws=3, submodel="det")
})

test_that("Fitted/residual methods work with ubmsFitDistsamp",{
  ubms_fitted <- methods::getMethod("fitted", "ubmsFit")
  ubms_residuals <- methods::getMethod("residuals", "ubmsFit")
  ubms_plot <- methods::getMethod("plot", "ubmsFit")

  set.seed(123)
  ft <- ubms_fitted(fit_line_hn, "state", draws=5)
  ft2 <- ubms_fitted(fit_line_hn, "det", draws=5)
  ft3 <- ubms_fitted(fit_line_abun, "state", draws=5)
  expect_equal(dim(ft), c(5,12))
  expect_equal(dim(ft2), c(5,48))
  expect_equal(dim(ft3), c(5,12))
  expect_equal(mean(ft), mean(ft3), tol=5)

  res <- ubms_residuals(fit_line_hn, "state", draws=5)
  res2 <- ubms_residuals(fit_line_hn, "det", draws=5)
  expect_equal(dim(res), c(5,12))
  expect_equal(dim(res2), c(5,48))

  pdf(NULL)
  rp <- plot_residuals(fit_line_hn, "state")
  rp2 <- plot_residuals(fit_line_hn, "det")
  rp3 <- ubms_plot(fit_line_hn)
  mp <- plot_marginal(fit_line_hn, "state")
  dev.off()

  expect_is(rp, "gg")
  expect_is(rp2, "gg")
  expect_is(rp3, "gtable")
  expect_is(mp, "gtable")
})

test_that("sim_state works for ubmsFitDistsamp",{
  sdens <- sim_state(fit_line_hn, 1:3)
  adens <- sim_state(fit_line_abun, 1:3)
  expect_equal(mean(sdens), mean(adens), tol=4)
})

test_that("getP and sim_p for ubmsFitDistsamp work",{

  plist <- lapply(line_mods, function(x) sim_p(x, 1:3))
  expect_true(all(sapply(plist, function(x) all(dim(x)==c(3,48)))))

  plist2 <- lapply(point_mods, function(x) sim_p(x, 1:3))
  expect_true(all(sapply(plist2, function(x) all(dim(x)==c(3,150)))))

  gplist <- lapply(line_mods, function(x) getP(x, 3))
  expect_true(all(sapply(gplist, function(x) all(dim(x)==c(12,4,3)))))

  gplist2 <- lapply(point_mods, function(x) getP(x, 3))
  expect_true(all(sapply(gplist2, function(x) all(dim(x)==c(30,5,3)))))
})

test_that("hist function works for ubmsFitDistsamp",{

  ubms_hist <- methods::getMethod("hist", "ubmsFitDistsamp")

  line_hist <- lapply(line_mods[c(1,2,4)], ubms_hist, draws=3)
  expect_true(all(sapply(line_hist, inherits, "gg")))
  point_hist <- lapply(point_mods[1:2], ubms_hist, draws=3)
  expect_true(all(sapply(point_hist, inherits, "gg")))

  fit_line_pcov <- suppressWarnings(stan_distsamp(~habitat~1, ltUMF, chains=2,
                                    iter=200, refresh=0))
  expect_warning(hwarn <- ubms_hist(fit_line_pcov))
  expect_is(hwarn, "gg")


  pdf(NULL)
  line_hist[[1]]
  line_hist[[2]]
  line_hist[[3]]
  point_hist[[1]]
  point_hist[[2]]
  hwarn
  dev.off()
})

test_that("ubmsResponseDistsamp methods work",{

  resp <- fit_line_hn@response
  expect_equal(get_K(resp), max(rowSums(fit_line_hn@data@y)) + 50)
  expect_equal(get_K(resp, 40), 40)
  expect_error(get_K(resp, 10))

  resp2 <- resp
  resp2@y[1,] <- NA
  expect_equal(length(get_Kmin(resp2)), nrow(resp2@y)-1)

  resp3 <- resp
  resp3@units_in <- "km"
  resp3@units_out <- "kmsq"
  expect_equal(mean(get_area_adjust(resp)), 16, tol=1e-3)
  expect_equal(mean(get_area_adjust(resp3)), 16e4, tol=1)
})
