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

ubmsFit <- function(model, call, data, response, submodels, ...){
  #Find missing
  response <- update_missing(response, submodels)
  submodels <- update_missing(submodels, response)
  
  #Fit model
  fit <- fit_model(model, response, submodels, ...)
  
  #Construct output
  new(fit_class(model), call=call, data=data, stanfit=fit,
      response=response, submodels=submodels, loo=get_loo(fit))
}

fit_class <- function(mod){
  cap <- paste0(toupper(substr(mod,1,1)), substr(mod,2,nchar(mod)))
  paste0("ubmsFit",cap)
}


#' @importFrom rstan extract
#' @export
setMethod("extract", "ubmsFit", 
  function(object, pars, permuted=TRUE, inc_warmup=FALSE, include=TRUE){
  rstan::extract(object@stanfit, pars, permuted, inc_warmup, include)
})

#Needs to be moved?
submodel_types <- function(object){
  names(object@submodels@submodels)
}


#' @importFrom rstan traceplot
#' @export
setMethod("traceplot", "ubmsFit", function(object, ...){
  rstan::traceplot(object@stanfit, ...)
})

#Fit stan model
#' @include inputs.R
fit_model <- function(name, response, submodels, ...){
  inp <- build_stan_inputs(response, submodels)  
  fit <- sampling(stanmodels[[name]], data=inp$stan_data, pars=inp$pars, ...)
  process_stanfit(fit, submodels)
}

#Do some cleanup on stanfit object
process_stanfit <- function(object, submodels){
  if(object@mode == 2L || object@mode == 1L){
    stop("Fitting model failed", call.=FALSE)
  }
  new_names <- stanfit_names(submodels@submodels)
  object@sim$fnames_oi[1:length(new_names)] <- new_names
  object
}

#Functions for renaming parameters in stanfit output object
stanfit_names <- function(submodels){
  c(stanfit_beta_names(submodels),
    stanfit_b_names(submodels),
    stanfit_sigma_names(submodels))
}

stanfit_beta_names <- function(submodels){
  types <- names(submodels)
  pars <- lapply(submodels, beta_names)
  names(pars) <- NULL
  for (i in 1:length(pars)){
    pars[[i]] <- paste0("beta_",types[i],"[",pars[[i]],"]")
  }
  unlist(pars)
}

stanfit_b_names <- function(submodels){
  types <- names(submodels)
  pars <- lapply(submodels, b_names)
  names(pars) <- NULL
  for (i in 1:length(pars)){
    if(all(is.na(pars[[i]]))){
      pars[[i]] <- character(0)
    } else {
      pars[[i]] <- paste0("b_",types[i],"[",pars[[i]],"]")
    }
  }
  unlist(pars)
}

stanfit_sigma_names <- function(submodels){
  nm <- unlist(lapply(submodels, sigma_names))
  if(all(is.na(nm))) return(character(0))
  nm <- nm[!is.na(nm)]
  for (i in 1:length(nm)){
    nm[i] <- gsub(" ", paste0("_",names(nm)[i]), nm, fixed=TRUE)
  }
  nm <- gsub("1|", "(Intercept) ", nm, fixed=TRUE)
  names(nm) <- NULL
  nm
}

get_loo <- function(stanfit, cores=getOption("mc.cores", 1)){
  loglik <- loo::extract_log_lik(stanfit, merge_chains=FALSE)
  r_eff <- loo::relative_eff(exp(loglik), cores=cores)
  loo::loo(loglik, r_eff=r_eff, cores=cores)
}
