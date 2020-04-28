#' Methods for ubmsFit objects

#' These methods are similar to methods defined for objects of class \code{lm},
#' \code{merMod}, \code{stanreg}. The most important methods are documented
#' separately. Links to those pages are provided in the \strong{See Also}
#' section below.
#'
#' @name ubmsFit-methods
#'
#' @param object,x A fitted model of class \code{ubmsFit}
#' @param type,i A submodel type, e.g. \code{"det"} for the detection model
#' @param ... Currently ignored
#' 
#' @seealso 
#' \itemize{
#'  \item The \code{\link{posterior_predict}} for drawing from the posterior
#'    distribution of the observed outcome or unobserved latent state
#' }
#'
#' @include fit.R
NULL

#' @rdname ubmsFit-methods
#' @export
setMethod("[", c("ubmsFit", "character", "missing", "missing"),
  function(x, i){
  x@submodels[i]
})

#' @rdname ubmsFit-methods
#' @importFrom rstantools nsamples
#' @export
setMethod("nsamples", "ubmsFit", function(object, ...){
  sum(object@stanfit@sim$n_save - object@stanfit@sim$warmup)
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

#' @rdname ubmsFit-methods
#' @export
setMethod("summary", "ubmsFit", function(object, type, ...){
  get_summary(object, type)
})

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
