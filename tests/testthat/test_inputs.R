context("Stan input generation")

test_that("stan inputs are built correctly", {
  sc <- data.frame(x1=rnorm(5), group=factor(c("a","b","a","b","a")))
  state <- ubmsSubmodel("Occ", "state", sc, ~x1+(1|group), "plogis",
                        uniform(-5,5), normal(0,2.5), gamma(1,1))
  det <- ubmsSubmodel("Det", "det", sc, ~x1, "plogis",
                      uniform(-5,5), normal(0,2.5), gamma(1,1))
  sl <- ubmsSubmodelList(state, det)
  y <- matrix(c(1,0,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")

  inp <- build_stan_inputs("occu", resp, sl, log_lik=FALSE)
  expect_is(inp, "list")
  expect_equal(names(inp), c("stan_data", "pars"))

  gs1 <- name_to_modelcode("occu")
  gs2 <- get_stan_data(resp)
  gs3 <- get_stan_data(state)
  gs4 <- get_stan_data(det)
  gs5 <- list(prior_dist_shape=c(0,0,0), prior_pars_shape=matrix(rep(0,6),nrow=3),
              prior_dist_scale=c(0,0,0), prior_pars_scale=matrix(rep(0,6),nrow=3))
  gs_all <- c(gs1,gs2,gs3,gs4,gs5)
  expect_equal(inp$stan_data, gs_all)
})

test_that("parameter list for stan is generated correctly",{
  sc <- data.frame(x1=rnorm(5), group=factor(c("a","b","a","b","a")))
  state <- ubmsSubmodel("Occ", "state", sc, ~x1+(1|group), "plogis",
                       uniform(-5,5), normal(0,2.5), gamma(1,1))
  det <- ubmsSubmodel("Det", "det", sc, ~x1, "plogis",
                      uniform(-5,5), normal(0,2.5), gamma(1,1))
  sl <- ubmsSubmodelList(state, det)
  expect_equal(get_pars(det), "beta_det")
  expect_equal(get_pars(state), c("beta_state", "b_state", "sigma_state"))
  expect_equal(get_pars(sl, "occu", log_lik=FALSE),
               c("beta_state", "b_state", "sigma_state", "beta_det"))
  expect_equal(get_pars(sl, "occu", log_lik=TRUE),
               c("beta_state", "b_state", "sigma_state", "beta_det", "log_lik"))
})

test_that("get_stan_data pulls necessary info from response object",{
  y <- matrix(c(1,0,0,1,1,1,0,0,1), nrow=3, byrow=T)
  resp <- ubmsResponse(y, "binomial", "P")
  dat <- get_stan_data(resp)
  expect_is(dat, "list")
  expect_equal(names(dat), c("y", "y_dist", "z_dist", "M", "T", "Tsamp",
                             "Tsamp_size", "J", "R", "si", "K", "Kmin",
                             "aux1","aux2","aux3","n_aux1","n_aux2","n_aux3"))
  expect_equal(dat$y, as.vector(t(y)))
  expect_equal(dat$y_dist, 0)
  expect_equal(dat$z_dist, 1)
  expect_equal(dat$M, 3)
  expect_equal(dat$T, 1)
  expect_equal(dat$Tsamp, c(1,1,1))
  expect_equal(dat$Tsamp_size, 3)
  expect_equal(dat$J, matrix(c(3,3,3), ncol=1))
  expect_equal(dat$R, 9)
  expect_equal(dim(dat$si), c(3, 6))
  expect_equal(dat$K, max(y) + 20)
  expect_equal(dat$Kmin, matrix(c(1,1,1), ncol=1))
})

test_that("dist_code returns integer code for distribution", {
  expect_equal(dist_code("binomial"), 0)
  expect_equal(dist_code("P"),1)
})

test_that("get_stan_data pulls necessary info from submodel",{
  sc <- data.frame(x1=-1:3)
  submod <- ubmsSubmodel("Occ", "state", sc, ~x1, "plogis",
                         uniform(-5,5), normal(0,2.5), gamma(1,1))
  dat <- get_stan_data(submod)
  expect_is(dat, "list")
  expect_equal(names(dat),
               paste0(c("X","offset", "n_obs","n_fixed","n_group_vars",
                        "has_random","n_random", "Zdim", "Zw", "Zv", "Zu",
                        "prior_dist", "prior_pars"),
                      "_", submod@type))
  expect_equivalent(dat[[1]], model.matrix(submod))
  expect_equal(dat[[2]], model_offset(submod))
  expect_equal(dat[[3]], nrow(model.matrix(submod)))
  expect_equal(dat[[4]], ncol(model.matrix(submod)))
  expect_equal(dat[[5]], get_group_vars(submod@formula))
  expect_equal(dat[[6]], has_random(submod))
  expect_equivalent(dat[8:11], get_sparse_Z(Z_matrix(submod)))
  expect_equivalent(dat[[12]], c(2,1,5))
  expect_equivalent(dat[[13]], matrix(c(-5,5,0,0,1.581139,0,1,1,0),nrow=3), tol=1e-6)
})

test_that("get_stan_data pulls only priors from scalar submodel",{
  ss <- ubmsSubmodelScalar("Fake", "fake", "plogis", normal(0,2.5))
  pr <- process_priors(ss)
  names(pr) <- paste0(names(pr),"_fake")
  expect_equal(get_stan_data(ss), pr)
})

test_that("get_sparse_Z collapses Z into sparse parts",{

  Z1 <- matrix(0, nrow=0, ncol=0)
  Z2 <- matrix(c(1,0,0,0, 0,0,1,0, 0,0,0,0, 0,0,1,0), nrow=4)

  expect_equal(get_sparse_Z(Z1),
               list(Zdim=c(0,0,1,1,1), Zw=as.array(0),
                    Zv=as.array(0), Zu=as.array(0)))

  expect_equal(get_sparse_Z(Z2),
               list(Zdim=c(nrow(Z2), ncol(Z2), w=3, v=3, u=5),
                    Zw=c(1,1,1), Zv=c(1,2,4),
                    Zu=c(1,2,2,4,4)))
})

test_that("get_group_vars returns number of grouping variables",{
  expect_equal(get_group_vars(~x), 0)
  expect_equal(get_group_vars(~(1|x)), 1)
  expect_equal(get_group_vars( ~(1|x) + (1|y)), 2)
})

test_that("get_nrandom returns number of levels of each grouping variable",{
  dat <- data.frame(x=factor(c("a","b","c")), y=factor(c("d","e","e")))
  expect_equal(get_nrandom(~x, dat), as.array(0))
  expect_equal(get_nrandom(~(1|x), dat), as.array(3))
  form <- ~(1|x) + (1|y)
  expect_equal(get_nrandom(form, dat), as.array(c(3,2)))
  form <- ~(1|x/y)
  expect_error(get_nrandom(form, dat),
               "Nested random effects (using / and :) are not supported", fixed=TRUE)
  form <- ~(1|x:y)
  expect_error(get_nrandom(form, dat),
               "Nested random effects (using / and :) are not supported", fixed=TRUE)
})

test_that("get_nrandom handles situation where a factor level is missing in data",{
  dat <- data.frame(x=factor(c("a","b","c"),levels=letters[1:4]))
  expect_equivalent(get_nrandom(~(1|x), dat), 4)
})

test_that("split_formula works",{
  inp <- ~1~1
  expect_equal(split_formula(inp), list(det=~1, state=~1))
  inp <- ~x1~x2
  expect_equal(split_formula(inp), list(det=~x1, state=~x2))
  inp <- ~x1+(1|x2)~x3
  expect_equal(split_formula(inp), list(det=~x1+(1|x2), state=~x3))
  inp <- ~1~verylong1 + verylong2 + verylong3 + verylong4 + verylong5 + verylong6 +
    verylong7 + verylong8 + verylong9 + verylong10 +
    (1|fake)
  expect_silent(split_formula(inp))
  expect_equal(split_formula(inp), list(det=~1, state=~verylong1 + verylong2 + verylong3 + verylong4 + verylong5 + verylong6 +
    verylong7 + verylong8 + verylong9 + verylong10 +
    (1|fake)))
  inp <- ~x1
  expect_error(split_formula(inp))
  inp <- y~x1
  expect_error(split_formula(inp))
})

