#' @include fit.R 
setClass("ubmsFitPcount", contains = "ubmsFit")

#' @export
stan_pcount <- function(formula, data, K=NULL, mixture="P", ...){
  
  forms <- split_formula(formula)

  #Need to process data first
  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist=mixture, K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(data), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)
  
  fit <- fit_model("pcount", response, submodels, ...) 

  new("ubmsFitPcount", call=match.call(), data=data, stanfit=fit, 
      response=response, submodels=submodels, loo=get_loo(fit))
}

