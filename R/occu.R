#' @include fit.R 
setClass("ubmsFitOccu", contains = "ubmsFit")

#' @export
stan_occu <- function(formula, data, ...){
  
  pformula <- split_formula(formula)[[1]]
  psiformula <- split_formula(formula)[[2]]

  #Need to process data first
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(data), psiformula, "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, det)

  inp <- build_stan_inputs(data, submodels)

  fit <- sampling(stanmodels$occu, data=inp$stan_data, pars=inp$pars, ...)

  fit <- process_stanfit(fit, submodels)
  submodels <- add_estimates(submodels, fit)

  new("ubmsFitOccu", call=match.call(), data=data, stanfit=fit, 
      WAIC=get_waic(fit), submodels=submodels)
}

