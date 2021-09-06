setOldClass("psis_loo")

#' @importClassesFrom rstan stanfit
#' @include submodel.R
#' @include response.R
setClass("ubmsFit",
  slots=c(call="call",
          data="unmarkedFrame",
          stanfit="stanfit",
          response="ubmsResponse",
          submodels="ubmsSubmodelList",
          loo="psis_loo"
          )
)

# Child class for occupancy-type models
setClass("ubmsFitOccu", contains = "ubmsFit")

# Child class for abundance/N-mixture type models
setClass("ubmsFitAbun", contains = "ubmsFit")

ubmsFit <- function(model, call, data, response, submodels, ...){
  #Find missing
  response <- update_missing(response, submodels)
  submodels <- update_missing(submodels, response)

  #Fit model
  fit <- fit_model(model, response, submodels, ...)

  #Remove placeholder submodels
  submodels <- remove_placeholders(submodels)

  #Construct output
  new(fit_class(model), call=call, data=data, stanfit=fit,
      response=response, submodels=submodels, loo=get_loo(fit))
}

remove_placeholders <- function(submodels){
  not_place <- !sapply(submodels@submodels, is_placeholder)
  submodels@submodels <- submodels@submodels[not_place]
  submodels
}

fit_class <- function(mod){
  cap <- paste0(toupper(substr(mod,1,1)), substr(mod,2,nchar(mod)))
  paste0("ubmsFit",cap)
}

#Fit stan model
#' @include inputs.R
fit_model <- function(name, response, submodels, ...){
  model <- name_to_stanmodel(name, submodels)
  inp <- build_stan_inputs(name, response, submodels)
  mod <- stanmodels[[model]]
  mod@model_name <- name
  fit <- sampling(mod, data=inp$stan_data, pars=inp$pars, ...)
  process_stanfit(fit, submodels)
}

name_to_stanmodel <- function(name, submodels){
  has_spatial <- any(sapply(submodels@submodels,
                            function(x) inherits(x, "ubmsSubmodelSpatial")))
  if(has_spatial) return("spatial")
  if(name == "colext") return("colext")
  return("single_season")
}

#Do some cleanup on stanfit object
process_stanfit <- function(object, submodels){
  if(object@mode == 2L || object@mode == 1L){
    stop("Fitting model failed", call.=FALSE)
  }
  new_names <- stanfit_names(submodels)
  object@sim$fnames_oi[1:length(new_names)] <- new_names
  object
}

setGeneric("stanfit_names", function(object, ...) standardGeneric("stanfit_names"))

setMethod("stanfit_names", "ubmsSubmodel", function(object, ...){
  out <- paste0("beta_",object@type,'[',beta_names(object),']')
  if(has_random(object)){
    bn <- paste0("b_",object@type,"[",b_names(object),']')
    sn <- paste0("sigma_",object@type,"[",sigma_names(object),"]")
    out <- c(out, bn, sn)
  }
  out
})

setMethod("stanfit_names", "ubmsSubmodelList", function(object, ...){
  out <- unlist(lapply(object@submodels, stanfit_names))
  out <- unname(out)
  out
})

get_loo <- function(stanfit, cores=getOption("mc.cores", 1)){
  loglik <- loo::extract_log_lik(stanfit, merge_chains=FALSE)
  r_eff <- loo::relative_eff(exp(loglik), cores=cores)
  loo::loo(loglik, r_eff=r_eff, cores=cores)
}
