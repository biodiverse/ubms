#' Extract Model Residuals
#'
#' Extract residuals for a given submodel from a \code{ubmsFit} object.
#' Residuals are calculated separately for each submodel
#' using the posterior predictive distribution of the latent state z, 
#' following Wright et al. (2019). 
#' 
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to get residuals for, for example \code{"det"}
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param ... Currently ignored
#'
#' @return A matrix of residual values with dimension \code{draws} by 
#'   observations. Note that for occupancy models, calculation of residuals
#'   for the detection submodel is conditional on \eqn{z = 1}, so residuals
#'   for an observation in a posterior draw where \eqn{z = 0} are assigned 
#'   value \code{NA} (Wright et al. 2019).
#' 
#' @references Wright, W. J., Irvine, K. M., & Higgs, M. D. (2019). Identifying 
#'   occupancy model inadequacies: can residuals separately assess detection 
#'   and presence? Ecology, 100(6), e02703.
#' 
#' @include fit.R
#' @importFrom stats residuals
#' @export
setMethod("residuals", "ubmsFit", function(object, submodel, draws=NULL, ...){
  samples <- get_samples(object, draws)
  sim_res(object, submodel, samples)
})


#Internal function for calculating residuals
setGeneric("sim_res", function(object, ...) standardGeneric("sim_res"))

setMethod("sim_res", "ubmsFit", function(object, ...){
  stop("No available method for this fit type", call.=FALSE)
})

#' @include occu.R
setMethod("sim_res", "ubmsFitOccu", function(object, submodel, samples, ...){

  z <- sim_z(object, samples=samples, re.form=NULL)
  lp <- sim_lp(object, submodel, samples=samples, transform=TRUE,
               newdata=NULL, re.form=NULL)
  
  if(identical(submodel, "state")){
    res <- z - lp
  } else if(identical(submodel, "det")){
    y <- object@data@y
    J <- ncol(y)
    ylong <- as.vector(t(y))
    zrep <- z[, rep(1:ncol(z), each=J)]
    z1_mask <- zrep == 1
    res <- matrix(rep(ylong, each=nrow(lp)), nrow=nrow(lp)) - lp
    res[!z1_mask] <- NA #residuals conditional on z = 1
  }
  res
})
