#' @include fit.R 
setClass("ubmsFitOccu", contains = "ubmsFit")

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

