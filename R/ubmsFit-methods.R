#' Methods for ubmsFit objects

#' These methods are similar to methods defined for objects of class \code{lm},
#' \code{merMod}, \code{stanreg}. The most important methods are documented
#' separately. Links to those pages are provided in the \strong{See Also}
#' section below.
#'
#' @name ubmsFit-methods
#'
#' @param object,x A fitted model of class \code{ubmsFit}
#' @param submodel,i A submodel, e.g. \code{"det"} for the detection model
#' @param cores Number of parallel cores to use in calculation
#' @param y Currently ignored
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
#' @importFrom gridExtra grid.arrange
#' @export
setMethod("plot", "ubmsFit", function(x, ...){  
  submods <- submodel_types(x)
  pl <- lapply(submods, function(s) plot_residuals(x, s, draws=6))
  grid.arrange(grobs=pl) 
})


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
    to_print <- summary(object, i)[,c(1,3,4,8:10)]
    names(to_print)[1:2] <- c("Estimate", "SD")
    cat(paste0(object[i]@name,":\n"))
    print(to_print, digits=3)
    cat("\n")
  }

  cat(paste0("LOOIC: ", round(object@loo$estimates[3,1], 3)))
  cat("\n")
})

#' @rdname ubmsFit-methods
#' @importFrom unmarked summary
#' @export
setMethod("summary", "ubmsFit", function(object, submodel, ...){
  sm <- object[submodel]
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
})

#' @rdname ubmsFit-methods
#' @importFrom loo loo
#' @export
setMethod("loo", "ubmsFit", function(x, ..., cores=getOption("mc.cores", 1)){
  loglik <- loo::extract_log_lik(x@stanfit, merge_chains=FALSE)
  r_eff <- loo::relative_eff(exp(loglik), cores=cores)
  loo::loo(loglik, r_eff=r_eff, cores=cores)
})


#' @rdname ubmsFit-methods
#' @importFrom loo waic
#' @export
setMethod("waic", "ubmsFit", function(x, ...){
  loglik <- loo::extract_log_lik(x@stanfit)
  loo::waic(loglik)
})


#' Extract Samples From a ubmsFit Model
#' 
#' Extract samples from a \code{ubmsFit} model
#'
#' @param object A \code{ubmsFit} object
#' @param pars An optional character vector providing parameter
#'  names of interest. If not specified, all parameters are used
#' @param permuted Logical. If \code{TRUE}, draws are permuted
#'  and merged; if \code{FALSE}, the original order is kept
#' @param inc_warmup Logical. If \code{TRUE}, warmup iterations
#'  are included; if \code{FALSE} they are discarded.
#' @param include Logical. If \code{TRUE} provided parameter names
#'  in \code{pars} are kept; if \code{FALSE} they are excluded.
#'
#' @return If \code{permuted=TRUE}, a list; if \code{permuted=FALSE}, 
#'  an array.
#'
#' @importFrom rstan extract
#' @export
setMethod("extract", "ubmsFit", 
  function(object, pars, permuted=TRUE, inc_warmup=FALSE, include=TRUE){
  rstan::extract(object@stanfit, pars, permuted, inc_warmup, include)
})

#' Markov Chain Traceplots
#' 
#' Draws traceplots for chains from a \code{ubmsFit} object
#' 
#' @param object A \code{ubmsFit} object
#' @param ... Arguments passed to \code{rstan::traceplot}

#' @importFrom rstan traceplot
#' @export
setMethod("traceplot", "ubmsFit", function(object, ...){
  rstan::traceplot(object@stanfit, ...)
})

