setOldClass("waic")

#' @importClassesFrom rstan stanfit
#' @importFrom loo waic
#' @include submodel.R
setClass("ubmsFit",
  slots=c(call="call",
          data="unmarkedFrame",
          stanfit="stanfit",
          WAIC="waic",
          submodels="ubmsSubmodelList"
          )
)

#' @importFrom rstan extract
#' @export
setMethod("extract", "ubmsFit", 
  function(object, pars, permuted=TRUE, inc_warmup=FALSE, include=TRUE){
  rstan::extract(object@stanfit, pars, permuted, inc_warmup, include)
})

setMethod("[", c("ubmsFit", "character", "missing", "missing"),
  function(x, i){
  x@submodels[i]
})

submodel_types <- function(object){
  names(object@submodels@submodels)
}

get_summary <- function(object, type){
  sm <- object[type]
  out <- rstan::summary(object@stanfit, beta_par(sm))
  out <- as.data.frame(out$summary)
  rownames(out) <- beta_names(sm)

  if(has_random(sm)){
    random <- rstan::summary(object@stanfit, sig_par(sm))
    random <- as.data.frame(random$summary)
    rownames(random) <- sigma_names(sm)
    out <- rbind(out, random)
  }
  out
}

#' @export
setMethod("summary", "ubmsFit", function(object, type, ...){
  get_summary(object, type)
})

setMethod("show", "ubmsFit", function(object){
  
  cat("\nCall:\n")
  print(object@call)
  cat("\n")

  for (i in submodel_types(object)){
    to_print <- get_summary(object, i)[,c(1,3,4,8:10)]
    names(to_print)[1:2] <- c("Estimate", "SD")
    cat(paste0(object[i]@name,":\n"))
    print(to_print, digits=3)
    cat("\n")
  }

  cat(paste0("WAIC: ", round(object@WAIC$estimates[3,1], 3)))
  cat("\n")
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

get_waic <- function(stanfit){
  loo::waic(loo::extract_log_lik(stanfit))
}
