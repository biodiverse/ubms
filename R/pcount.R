#' @include fit.R 
setClass("ubmsFitPcount", contains = "ubmsFit")

#' @export
stan_pcount <- function(formula, data, K=NULL, mixture="P", ...){
  
  pformula <- split_formula(formula)[[1]]
  lambdaformula <- split_formula(formula)[[2]]

  #Need to process data first
  state <- ubmsSubmodel("Abundance", "state", siteCovs(data), lambdaformula, "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, det)

  inp <- build_stan_inputs(data, submodels, K=K, mixture=mixture)

  fit <- sampling(stanmodels$pcount, data=inp$stan_data, pars=inp$pars, ...)

  fit <- process_stanfit(fit, submodels)
  submodels <- add_estimates(submodels, fit)

  new("ubmsFitPcount", call=match.call(), data=data, stanfit=fit, 
      WAIC=get_waic(fit), submodels=submodels)
}

