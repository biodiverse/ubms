#' @include fit.R 
setClass("ubmsFitOccuRN", contains = "ubmsFitOccu")

#' @export
stan_occuRN <- function(formula, data, K=20, ...){
  
  forms <- split_formula(formula)

  #Need to process data first
  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist="P", K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(data), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)
 
  fit <- fit_model("occuRN", response, submodels, ...) 

  new("ubmsFitOccuRN", call=match.call(), data=data, stanfit=fit, 
      response=response, submodels=submodels, loo=get_loo(fit))
}

