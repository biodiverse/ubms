#' @export
stan_occu <- function(formula, data, ...){
  
  pformula <- as.formula(formula[[2]])
  psiformula <- as.formula(paste0(formula[[1]],
                           deparse(formula[[3]])))

  #Need to process data first
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "logit")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), pformula, "logit")

  inp <- build_stan_inputs(umf, state, det)

  fit <- rstan::sampling(stanmodels$occupancy, data=inp$stan_data, 
                         pars=inp$pars, ...)
  
  #Combine submodels and add estimate summary
  submodels <- ubmsSubmodelList(state, det)
  submodels <- add_estimates(submodels, fit)

  new("ubmsFit", call=match.call(), psiformula=psiformula, 
      pformula=pformula, data=data, stanfit=fit, 
      WAIC=loo::waic(loo::extract_log_lik(fit)),
      submodels=submodels)
}

