setClass("ubmsFitList", slots = c(models = "list"))

setGeneric("fitList", function(...){
  unmarked::fitList(...)
})

#' @export
setMethod("fitList", "ubmsFit", function(...){
  mods <- list(...)
	if (is.null(names(mods))) names(mods) <- ""
	names(mods)[names(mods) == ""] <- NA
	mc <- sys.calls()[[1]][-1]
  modnames <- ifelse(is.na(names(mods)), as.character(mc), names(mods))
  out <- new("ubmsFitList", models=mods)
  names(out@models) <- modnames
  out
})

setMethod("$", "ubmsFitList", function(x, name){
  x@models[[name]]
})

setMethod("names", "ubmsFitList", function(x){
  names(x@models)
})

setMethod("[[", c("ubmsFitList", "numeric","missing"), 
          function(x, i) x@models[[i]])

setMethod("[", c("ubmsFitList", "numeric","missing","missing"), 
          function(x, i) x@models[i])

#' @importFrom loo loo_compare loo_model_weights
#' @export
setMethod("modSel", "ubmsFitList", function(object, ...){
  loos <- lapply(object@models, loo, ...)
  elpd <- sapply(loos, function(x) x$estimates[1])
  p_loo <- sapply(loos, function(x) x$estimates[2]) 
  compare <- loo::loo_compare(loos)[names(elpd),]
  wts <- as.vector(loo::loo_model_weights(loos))
  out <- data.frame(elpd=elpd, nparam=p_loo, elpd_diff=compare[,1], 
                    se_diff=compare[,2], weight=wts)
  out[order(out$elpd_diff, decreasing=TRUE),]
})
