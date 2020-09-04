#' Multinomial-Poisson Mixture Model
#'
#' This function fits the multinomial-Poisson mixture model useful for
#' data collected via survey methods such as removal or double observer sampling
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFrameMultinomPois}} object
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitMultinomPois} object describing the model fit.
#'
#' @seealso \code{\link{multinomPois}}, \code{\link{unmarkedFrameMultinomPois}}
#' @export
stan_multinomPois <- function(formula, data, ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)
  pifun_type <- get_pifun_type(umf)

  response <- ubmsResponseMultinomPois(getY(umf), pifun_type, "P")
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
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
setClass("ubmsFitMultinomPois", contains = "ubmsFit")


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

#Goodness-of-fit---------------------------------------------------------------

#' @include gof.R


#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
