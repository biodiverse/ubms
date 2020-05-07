#' @include fit.R 
setClass("ubmsFitOccuRN", contains = "ubmsFitOccu")

#' @export
stan_occuRN <- function(formula, data, K=20, ...){
  
  pformula <- split_formula(formula)[[1]]
  lamformula <- split_formula(formula)[[2]]

  #Need to process data first
  state <- ubmsSubmodel("Abundance", "state", siteCovs(data), lamformula, "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, det)
  submodels <- find_missing(submodels, data)
  inp <- build_stan_inputs(submodels, data, K=K)

  fit <- sampling(stanmodels$occuRN, data=inp$stan_data, pars=inp$pars, ...)
  fit <- process_stanfit(fit, submodels)

  new("ubmsFitOccuRN", call=match.call(), data=data, stanfit=fit, 
      loo=get_loo(fit), submodels=submodels)
}

