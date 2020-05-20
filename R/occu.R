#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of
#' MacKenzie et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link{unmarkedFrameOccu}} object
#' @param ... Arguments passed to the \code{\link{stan}} call, such as 
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccu} object describing the model fit.
#'
#' @seealso \code{\link{occu}}, \code{\link{unmarkedFrameOccu}}
#' @export
stan_occu <- function(formula, data, ...){
  
  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), forms[[2]], "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occu", match.call(), data, response, submodels, ...)
}

#' @include fit.R 
setClass("ubmsFitOccu", contains = "ubmsFit")
