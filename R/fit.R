setOldClass("waic")

#' @importClassesFrom rstan stanfit
#' @importFrom loo waic
#' @include submodel.R
#' @export
setClass("ubmsFit",
  slots=c(call="call",
          psiformula="formula",
          pformula="formula",
          data="unmarkedFrame",
          stanfit="stanfit",
          WAIC="waic",
          submodels="ubmsSubmodelList"
          )
)

setMethod("show", "ubmsFit", function(object){
  
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  show(object@submodels)
  cat(paste0("WAIC: ", round(object@WAIC$estimates[3,1], 3)))
  cat("\n")

})

#' @importFrom rstan extract
#' @export
setMethod("extract", "ubmsFit", 
  function(object, pars, permuted=TRUE, inc_warmup=FALSE,
          include=TRUE){
  rstan::extract(object@stanfit, pars, permuted, inc_warmup, include)
})

setMethod("[", c("ubmsFit", "character", "missing", "missing"),
  function(x, i){

  types <- names(x@submodels@submodels)
  if(! i %in% types){
    stop(paste("Possible types are:", paste(types, collapse=", ")),
         call. = FALSE)
  }
  x@submodels@submodels[[i]]
})

#' @importFrom rstan traceplot
#' @export
setMethod("traceplot", "ubmsFit", function(object, ...){
  rstan::traceplot(object@stanfit, ...)
})

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
  pars <- lapply(submodels, function(x) x@beta_names)
  names(pars) <- NULL
  for (i in 1:length(pars)){
    pars[[i]] <- paste0("beta_",types[i],"[",pars[[i]],"]")
  }
  unlist(pars)
}

stanfit_b_names <- function(submodels){
  types <- names(submodels)
  pars <- lapply(submodels, function(x) x@b_names)
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
  nm <- unlist(lapply(submodels, function(x) x@sigma_names))
  if(all(is.na(nm))) return(character(0))
  nm <- nm[!is.na(nm)]
  for (i in 1:length(nm)){
    nm[i] <- gsub(" ", paste0("_",names(nm)[i]), nm, fixed=TRUE)
  }
  nm <- gsub("1|", "(Intercept) ", nm, fixed=TRUE)
  names(nm) <- NULL
  nm
}
