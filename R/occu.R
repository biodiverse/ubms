#' @export
stan_occu <- function(formula, data, ...){
  
  pformula <- as.formula(formula[[2]])
  psiformula <- as.formula(paste0(formula[[1]], deparse(formula[[3]])))

  #Need to process data first
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(data), psiformula, "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(data), pformula, "plogis")

  inp <- build_stan_inputs(data, state, det)

  fit <- sampling(stanmodels$occupancy, data=inp$stan_data, pars=inp$pars, ...)
  check_stanfit(fit) 
  #Combine submodels and add estimate summary
  submodels <- ubmsSubmodelList(state, det)
  submodels <- add_estimates(submodels, fit)

  new("ubmsFit", call=match.call(), psiformula=psiformula, 
      pformula=pformula, data=data, stanfit=fit, 
      WAIC=loo::waic(loo::extract_log_lik(fit)),
      submodels=submodels)
}

