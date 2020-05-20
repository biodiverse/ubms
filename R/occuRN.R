#' Fit the Occupancy Model of Royle and Nichols (2003)
#'
#' Fit the occupancy model of Royle and Nichols (2003), which relates 
#' probability of detection of the species to the number of individuals 
#' available for detection at each site.
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFrameOccu}} object
#' @param K Integer upper index of integration for N-mixture. This should be
#'  set high enough so that it does not affect the parameter estimates. 
#'  Note that computation time will increase with K.
#' @param ... Arguments passed to the \code{\link{stan}} call, such as 
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccuRN} object describing the model fit.
#'
#' @seealso \code{\link{occuRN}}, \code{\link{unmarkedFrameOccu}}
#' @export
stan_occuRN <- function(formula, data, K=20, ...){
  
  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist="P", K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occuRN", match.call(), data, response, submodels, ...)
}

#' @include fit.R 
setClass("ubmsFitOccuRN", contains = "ubmsFitOccu")
