#' Extract Fitted Values
#'
#' Extract fitted values for a given submodel from a \code{ubmsFit} object.
#' Fitted values are calculated separately for each submodel
#' using the posterior predictive distribution of the latent state z,
#' following Wright et al. (2019).
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to get fitted values for, options are \code{"state"}
#'  or \code{"det"}
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param ... Currently ignored
#'
#' @return A matrix of fitted values with dimension \code{draws} by
#'   observations. Note that calculation of fitted values
#'   for the detection submodel is conditional on \eqn{z > 0}, so fitted values
#'   for an observation in a posterior draw where \eqn{z = 0} are assigned
#'   value \code{NA} (Wright et al. 2019).
#'
#' @references Wright, W. J., Irvine, K. M., & Higgs, M. D. (2019). Identifying
#'   occupancy model inadequacies: can residuals separately assess detection
#'   and presence? Ecology, 100(6), e02703.
#'
#' @include fit.R
#' @importFrom stats fitted
#' @export
setMethod("fitted", "ubmsFit", function(object, submodel, draws=NULL, ...){
  samples <- get_samples(object, draws)
  sim_fitted(object, submodel, samples)
})

setGeneric("sim_fitted", function(object, ...) standardGeneric("sim_fitted"))

setMethod("sim_fitted", "ubmsFit", function(object, submodel, samples, ...){
  stopifnot(submodel %in% c("state", "det"))
  lp <- sim_lp(object, submodel, transform=TRUE, newdata=NULL, samples=samples,
               re.form=NULL)
  if(submodel == "state") return(lp)

  #Detection, fitted values conditional on z = 1
  J <- object@response@max_obs
  z <- sim_z(object, samples, re.form=NULL)
  z <- z[, rep(1:ncol(z), each=J)]
  lp[z == 0] <- NA
  lp
})
