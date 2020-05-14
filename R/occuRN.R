#' @include fit.R 
setClass("ubmsFitOccuRN", contains = "ubmsFitOccu")

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

