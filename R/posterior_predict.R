#' Draw from the posterior predictive distribution
#'
#' Draw from the posterior predictive distribution after fitting a model.
#' You can draw from the posterior of the observed outcome \eqn{y} or
#' the latent unobserved state \eqn{z}.
#'
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param param Either \code{"y"} for the observed outcome or \code{"z"}
#'   for the unobserved latent state
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'   effects") are included. If \code{NA}, they are ignored
#' @param ... Currently ignored
#'
#' @return A matrix of simulations from the posterior predictive distribution.
#'   If \code{param = "z"}, the dimensions are \code{draws} by number of sites
#'   (or sites x primary periods in site-major order for dynamic models). If
#'   \code{param = "y"}, the dimensions are \code{draws} by sites x observations
#'   (or sites x primary periods x observations for dynamic models).
#'
#' @aliases posterior_predict
#' @method posterior_predict ubmsFit
#' @importFrom rstantools posterior_predict
#' @include fit.R
#' @export
setMethod("posterior_predict", "ubmsFit",
          function(object, param=c("y","z"), draws=NULL, re.form=NULL, ...){

  param <- match.arg(param, c("y", "z"))
  nsamp <- nsamples(object)
  samp_inds <- get_samples(object, draws)

  switch(param,
         "z" = sim_z(object, samples=samp_inds, re.form=re.form),
         "y" = sim_y(object, samples=samp_inds, re.form=re.form))
})


## simulation methods ---------------------------------------------------------

#Simulate the unobserved latent state
setGeneric("sim_z", function(object, ...) standardGeneric("sim_z"))

#Simulate the observed outcome
setGeneric("sim_y", function(object, ...) standardGeneric("sim_y"))


# utils -----------------------------------------------------------------------

# Generate z latent state matrix if it isn't provided
process_z <- function(object, samples, re.form, z){
  if(is.null(z)){
    z <- t(sim_z(object, samples=samples, re.form=re.form))
  } else {
    z <- t(z)
  }
  z
}
