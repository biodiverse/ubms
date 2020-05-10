#' @include fit.R 
setClass("ubmsFitOccu", contains = "ubmsFit")

#' @export
stan_occu <- function(formula, data, ...){
  
  pformula <- split_formula(formula)[[1]]
  psiformula <- split_formula(formula)[[2]]

  #Need to process data first
  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(data), psiformula, "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, det)
  inp <- build_stan_inputs(submodels, response)

  fit <- sampling(stanmodels$occu, data=inp$stan_data, pars=inp$pars, ...)
  fit <- process_stanfit(fit, submodels)

  new("ubmsFitOccu", call=match.call(), data=data, stanfit=fit, 
      loo=get_loo(fit), submodels=submodels)
}

