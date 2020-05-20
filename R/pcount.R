#' Fit the N-mixture model of Royle (2004)
#'
#' This function fits the single season N-mixture model of
#' Royle et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFramePCount}} object
#' @param K Integer upper index of integration for N-mixture. This should be
#'  set high enough so that it does not affect the parameter estimates. 
#'  Note that computation time will increase with K.
#' @param mixture Character specifying mixture: "P" is only option currently.
#' @param ... Arguments passed to the \code{\link{stan}} call, such as 
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitPcount} object describing the model fit.
#'
#' @seealso \code{\link{pcount}}, \code{\link{unmarkedFramePCount}}
#' @export
stan_pcount <- function(formula, data, K=NULL, mixture="P", ...){
  
  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist=mixture, K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("pcount", match.call(), data, response, submodels, ...)
}

#' @include fit.R 
setClass("ubmsFitPcount", contains = "ubmsFit")
