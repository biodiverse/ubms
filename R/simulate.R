#' @include fit.R simz.R simy.R
#' @export
setMethod("simulate", "ubmsFit",
          function(object, nsim=NULL, param=c("y","z"), samples=NULL, ...){
  param <- match.arg(param, c("y", "z"))
  keep <- get_sample_inds(object@stanfit, nsim, samples)
  z <- sim_z(object, keep)
  if(identical(param, "z")) return(z) 
  sim_y(object, z, keep)
})

get_sample_inds <- function(stanfit, nsim, samples){
  nsamples <- nrow(as.matrix(stanfit))
  if(!is.null(nsim) & !is.null(samples)){
    stop("Can't provide both nsim and samples", call.=FALSE)
  }
  if(!is.null(nsim) && nsim < nsamples){
    return(sample(1:nsamples, nsim, replace=FALSE))
  } else if(!is.null(samples)){
    smax <- max(samples)
    if(smax > nsamples || any(samples < 1)){
      stop("Invalid vector of sample indices", call.=FALSE)
    }
    return(samples)
  }
  1:nsamples
}
