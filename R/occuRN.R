#' @include fit.R 
setClass("ubmsFitOccuRN", contains = "ubmsFit")

#' @export
stan_occuRN <- function(formula, data, K=15, ...){
  
  pformula <- split_formula(formula)[[1]]
  lamformula <- split_formula(formula)[[2]]

  #Need to process data first
  state <- ubmsSubmodel("Abundance", "state", siteCovs(data), lamformula, "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, det)

  inp <- build_stan_inputs(data, submodels, K=K)

  fit <- sampling(stanmodels$occuRN, data=inp$stan_data, pars=inp$pars, ...)

  fit <- process_stanfit(fit, submodels)
  submodels <- add_estimates(submodels, fit)

  new("ubmsFitOccuRN", call=match.call(), data=data, stanfit=fit, 
      WAIC=get_waic(fit), submodels=submodels)
}

