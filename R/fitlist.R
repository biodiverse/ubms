setClass("ubmsFitList", slots = c(models = "list"))

setGeneric("fitList", function(...){
  unmarked::fitList(...)
})

#' Create a List of ubmsFit Models
#'
#' Create a list of ubmsFit models
#'
#' @param ... \code{ubmsFit} model objects, preferably named, or a list
#'  of such models
#'
#' @return An object of class \code{ubmsFitList} containing the list of models
#'
#' @aliases fitList fitList,list-method
#' @export
setMethod("fitList", "ubmsFit", function(...){
  mods <- list(...)
  mod_names <- names(mods)
  if(is.null(mod_names)) mod_names <- rep("", length(mods))
  mod_names[mod_names==""] <- NA
  obj_names=sapply(substitute(list(...))[-1], deparse)
  mod_names[is.na(mod_names)] <- obj_names[is.na(mod_names)]
  names(mods) <- mod_names
  fitList(mods)
})

setMethod("fitList", "list", function(...){
  mods <- list(...)[[1]]
  if(!inherits(mods[[1]],"ubmsFit", )) return(unmarked::fitList(fits=mods))
  if(is.null(names(mods))) names(mods) <- paste0("mod", 1:length(mods))
  new("ubmsFitList", models=mods)
})

#' Model Selection For a List of ubmsFit Models
#'
#' Construct a model selection table from a \code{ubmsFitList}
#'
#' @param object An object of class \code{ubmsFitList}
#' @param ... Currently ignored
#'
#' @return A \code{data.frame} of model fit information with one row per
#'  model in the input \code{ubmsFitList}. Models are ranked in descending
#'  order by expected log pointwise predictive density (\code{elpd}).
#' @seealso \code{\link[loo]{loo}}, \code{\link[loo]{loo_compare}}
#'
#' @importFrom loo loo_compare loo_model_weights
#' @importFrom unmarked modSel
#' @export
setMethod("modSel", "ubmsFitList", function(object, ...){
  #loos <- lapply(object@models, loo, ...)
  loos <- lapply(object@models, function(x) x@loo)
  elpd <- sapply(loos, function(x) x$estimates[1])
  p_loo <- sapply(loos, function(x) x$estimates[2])
  compare <- loo::loo_compare(loos)[names(elpd),]
  wts <- as.vector(loo::loo_model_weights(loos))
  out <- data.frame(elpd=elpd, nparam=p_loo, elpd_diff=compare[,1],
                    se_diff=compare[,2], weight=wts)
  out[order(out$elpd_diff, decreasing=TRUE),]
})
