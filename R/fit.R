setClass("umStanEstimate",
  slots=c(name="character", 
          fixed="data.frame",
          random="data.frame"),
  prototype=list(name=NA_character_,
                 fixed=data.frame(),
                 random=data.frame())
)

umStanEstimate <- function(object, mod_name, fixed_param, fixed_names,
                          random_param=NULL, random_names=NULL){

  fixed <- rstan::summary(object, fixed_param)$summary
  fixed <- as.data.frame(fixed)
  rownames(fixed) <- fixed_names

  random <- data.frame()
  if(!is.null(random_names)){
    random <- rstan::summary(object, random_param)$summary
    random <- as.data.frame(random)
    rownames(random) <- random_names
  }

  new("umStanEstimate", name=mod_name, fixed=fixed, random=random)

}

setClass("umStanEstimateList",
  slots=c(estimates="list")
)

umStanEstimateList <- function(...){
  new("umStanEstimateList", estimates=list(...))
}

#' @importClassesFrom rstan stanfit
#' @importFrom loo waic

setOldClass("waic")

setClass("umFitStan",
  slots=c(call="call",
          psiformula="formula",
          pformula="formula",
          data="unmarkedFrame",
          stanfit="stanfit",
          WAIC="waic",
          estimates="umStanEstimateList"
          )
)

setMethod("show", "umStanEstimate", function(object){
  
  out_df <- object@fixed[,c(1,3,4,8:10)]

  if(nrow(object@random)>0){
    out_df <- rbind(out_df, object@random[,c(1,3,4,8:10)])
  }
  names(out_df)[1:2] <- c("Estimate", "SD")

  cat(paste0(object@name,":\n"))
  print(out_df, digits=3)

})

setMethod("show", "umStanEstimateList", function(object){

  for (est in object@estimates){
    show(est)
    cat("\n")
  }

})

setMethod("show", "umFitStan", function(object){
  
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  show(object@estimates)
  cat(paste0("WAIC: ", round(object@WAIC$estimates[3,1], 3)))
  cat("\n")

})

#' @export
occuStan <- function(formula, data, ...){
  
  pformula <- as.formula(formula[[2]])
  psiformula <- as.formula(paste0(formula[[1]],
                           deparse(formula[[3]])))

  inp <- build_inputs(psiformula, pformula, umf)

  fit <- rstan::sampling(stanmodels$occupancy, data=inp$stan_data, 
                     pars=inp$params, ...)

  occ_mod <- umStanEstimate(fit, "Occupancy", "beta_occ", inp$fixed_names$occ,
                            "sigma_occ", inp$random_names$occ)

  det_mod <- umStanEstimate(fit, "Detection", "beta_det", inp$fixed_names$det,
                            "sigma_det", inp$random_names$det)

  new("umFitStan", call=match.call(), psiformula=psiformula, 
      pformula=pformula, data=data, stanfit=fit, 
      WAIC=loo::waic(loo::extract_log_lik(fit)),
      estimates=umStanEstimateList(occ_mod, det_mod))
}

