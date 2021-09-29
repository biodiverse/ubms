#' Fit the Multinomial-Poisson Mixture Model
#'
#' This function fits the multinomial-Poisson mixture model, useful for
#' data collected via survey methods such as removal or double observer sampling.
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFrameMPois}} object
#' @param prior_intercept_state Prior distribution for the intercept of the
#'  state (abundance) model; see \code{?priors} for options
#' @param prior_coef_state Prior distribution for the regression coefficients of
#'  the state model
#' @param prior_intercept_det Prior distribution for the intercept of the
#'  detection probability model
#' @param prior_coef_det Prior distribution for the regression coefficients of
#'  the detection model
#' @param prior_sigma Prior distribution on random effect standard deviations
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitMultinomPois} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(ovendata)
#' ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
#'                                 siteCovs=ovendata.list$covariates,
#'                                 type="removal")
#'
#' oven_fit <- stan_multinomPois(~1~scale(ufc), ovenFrame, chains=3, iter=300)
#' }
#'
#' @seealso \code{\link{multinomPois}}, \code{\link{unmarkedFrameMPois}}
#' @export
stan_multinomPois <- function(formula,
                              data,
                              prior_intercept_state = normal(0, 5),
                              prior_coef_state = normal(0, 2.5),
                              prior_intercept_det = logistic(0, 1),
                              prior_coef_det = logistic(0, 1),
                              prior_sigma = gamma(1, 1),
                              ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)
  pifun_type <- get_pifun_type(umf)

  if(has_spatial(forms)){
    split_umf <- extract_missing_sites(umf)
    umf <- split_umf$umf
    state <- ubmsSubmodelSpatial("Abundance", "state", siteCovs(umf), forms[[2]],
                                 "exp", prior_intercept_state, prior_coef_state,
                                 prior_sigma,
                                 split_umf$sites_augment, split_umf$data_aug)
  } else {
    state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp",
                          prior_intercept_state, prior_coef_state, prior_sigma)
  }

  response <- ubmsResponseMultinomPois(getY(umf), pifun_type, "P")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis",
                      prior_intercept_det, prior_coef_det, prior_sigma)
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("multinomPois", match.call(), data, response, submodels, ...)
}

get_pifun_type <- function(umf){
  input <- umf@piFun
  if(! input %in% c("doublePiFun","removalPiFun")){
    stop("Custom pi functions not supported")
  }
  switch(input, doublePiFun = {"double"}, removalPiFun = {"removal"})
}

#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitMultinomPois", contains = "ubmsFitAbun")


#Response class----------------------------------------------------------------

#' @include response.R
setClass("ubmsResponseMultinomPois", contains="ubmsResponse")

ubmsResponseMultinomPois <- function(y, y_dist, z_dist, max_primary = 1, K=NULL){
  stopifnot(inherits(y, "matrix"))
  out <- new("ubmsResponseMultinomPois", y = y, y_dist= y_dist, z_dist = z_dist,
          max_primary = max_primary, max_obs = get_max_obs(y, max_primary))
  out@missing <- is.na(as_vector(out))
  out@K <- get_K(out, K)
  out
}

#Bespoke find_missing method for multinomPois models
#Special situation with double observer sampling
#Also remove site if any observation is missing

#' @include missing.R
setMethod("find_missing", "ubmsResponseMultinomPois",
  function(object, submodels, ...){

  type <- object@y_dist

  y <- as_vector(object)
  M <- ncol(t(object))
  sc <- expand_model_matrix(submodels@submodels$state, object)

  if(type == "removal"){
    oc <- expand_model_matrix(submodels@submodels$det, object)
  } else if (type == "double"){
    oc <- model.matrix(submodels@submodels$det)
    oc <- oc[rep(1:nrow(oc), rep(c(1,2), M)),]
  }

  comb <- cbind(y, sc, oc)
  miss <- apply(comb, 1, function(x) any(is.na(x)))
  miss <- matrix(miss, nrow=nrow(t(object)))
  remove_sites <- apply(miss, 2, function(x) any(x))
  rep(remove_sites, each=object@max_obs)
})

setMethod("update_missing", c("ubmsSubmodelList", "ubmsResponseMultinomPois"),
  function(object, object2){
  response <- object2
  is_na <- find_missing(response, object)
  na_mat <- matrix(is_na, nrow=response@max_obs)
  object@submodels$state@missing <- apply(na_mat, 2, any)

  if(response@y_dist == "double"){
    obs_miss <- as.vector(na_mat[-3,])
  } else if(response@y_dist == "removal"){
    obs_miss <- as.vector(na_mat)
  }
  stopifnot(length(obs_miss) == length(object@submodels$det@missing))
  object@submodels$det@missing <- obs_miss

  object
})


#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitMultinomPois", function(object, samples, re.form, K=NULL, ...){
  resp <- object@response
  y <- resp@y

  lam_post <- t(sim_lp(object, "state", transform=TRUE, newdata=NULL,
                       re.form=re.form, samples=samples))
  lam_post[object["state"]@missing] <- NA
  p_post <- get_pi_for_multinom(object, samples)

  K <- object@response@K
  Kmin <- apply(y, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE)))
  kvals <- 0:K

  t(simz_multinom(y, lam_post, p_post, K, Kmin, kvals))
})

setMethod("sim_y", "ubmsFitMultinomPois",
          function(object, samples, re.form, z=NULL, K=NULL, ...){
  nsamp <- length(samples)
  M <- nrow(object@response@y)
  J <- object@response@max_obs
  z <- process_z(object, samples, re.form, z)
  p <- get_pi_for_multinom(object, samples)

  out <- array(NA, c(J, M, nsamp))
  for (i in 1:nsamp){
    for (m in 1:M){
      if(is.na(z[m,i])) next
      out[,m,i] <- stats::rmultinom(n=1, size=z[m,i], prob=p[m,,i])[1:J]
    }
  }
  matrix(out, nrow=nsamp, byrow=TRUE)
})

get_pi_for_multinom <- function(object, samples){
  resp <- object@response
  nsamp <- length(samples)
  M <- nrow(object@response@y)
  J <- resp@max_obs
  pi_fun <- switch(resp@y_dist, double = unmarked::doublePiFun,
                   removal = unmarked::removalPiFun)

  p_raw <- t(sim_lp(object, "det", transform=TRUE, newdata=NULL,
                  re.form=NULL, samples=samples))
  p_raw[object["det"]@missing] <- NA
  p_raw <- array(p_raw, c(nrow(p_raw) / M, M, nsamp))
  p_raw <- aperm(p_raw, c(2,1,3))

  p_post <- array(NA, c(M, J+1, nsamp))
  for (i in 1:nsamp){
    p_post[,1:J,i] <- pi_fun(p_raw[,,i])
    p_post[,J+1,i] <- 1 - rowSums(p_post[,1:J,i])
  }
  p_post
}


#Get detection probability-----------------------------------------------------

#' @include posterior_linpred.R
setMethod("sim_p", "ubmsFitMultinomPois", function(object, samples, ...){

  J <- object@response@max_obs
  p_array <- get_pi_for_multinom(object, samples)[,1:J,]
  p_array <- aperm(p_array, c(2,1,3))
  p_mat <- matrix(p_array, ncol=length(samples))
  t(p_mat)
})


#Method for fitted values------------------------------------------------------

setMethod("sim_fitted", "ubmsFitMultinomPois", function(object, submodel, samples, ...){
  stopifnot(submodel %in% c("state", "det"))
  if(submodel == "state"){
    lp <- sim_lp(object, submodel, transform=TRUE, newdata=NULL,
                 samples=samples, re.form=NULL)
    return(lp)
  }

  p <- sim_p(object, samples=samples)
  J <- object@response@max_obs
  z <- sim_z(object, samples, re.form=NULL)
  z <- z[, rep(1:ncol(z), each=J)]
  out <- z * p
  out
})
