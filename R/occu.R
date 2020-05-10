#' @include fit.R 
setClass("ubmsFitOccu", contains = "ubmsFit")

#' @export
stan_occu <- function(formula, data, ...){
  
  forms <- split_formula(formula)

  #Need to process data first
  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(data), forms[[2]], "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  fit <- fit_model("occu", response, submodels, ...) 

  new("ubmsFitOccu", call=match.call(), data=data, stanfit=fit, 
      response=response, submodels=submodels, loo=get_loo(fit))
}

