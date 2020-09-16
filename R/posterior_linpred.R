#' Posterior Distribution of the Linear Predictor
#'
#' Extract posterior draws of the linear predictor for a \code{ubmsFit}
#' submodel, possibly transformed by the inverse-link function.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param transform Should the linear predictor be transformed using the
#'  inverse link function?
#' @param submodel The name of the submodel, as a character string, for
#'  which to calculate the linear predictor
#' @param newdata Optional data frame of newdata to use when calculating the linear
#'  predictor. If not provided, the model matrix is used.
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'   effects") are included. If \code{NA}, they are ignored
#' @param ... Currently ignored
#'
#' @return A matrix of simulations from the posterior predictive distribution
#'  of the linear predictor. The dimensions are \code{draws} by number of
#'  linear predictor values (e.g., number of sites or number of observations).
#'
#' @aliases posterior_linpred
#' @method posterior_linpred ubmsFit
#' @importFrom rstantools posterior_linpred
#' @include fit.R
#' @export
setMethod("posterior_linpred", "ubmsFit",
          function(object, transform=FALSE, submodel, newdata=NULL, draws=NULL,
                   re.form=NULL, ...){

  check_type(submodel, submodel_types(object))
  samp_inds <- get_samples(object, draws)

  sim_lp(object, submodel=submodel, transform=transform, newdata=newdata,
         samples=samp_inds, re.form=re.form)
})


setGeneric("sim_lp", function(object, ...) standardGeneric("sim_lp"))

setMethod("sim_lp", "ubmsFit", function(object, submodel, transform, newdata,
                                        samples, re.form, ...){
  sm <- object[submodel]
  beta <- extract(object, beta_par(sm))[[1]]
  lp <- model.matrix(sm, newdata) %*% t(beta[samples,,drop=FALSE])

  if(has_random(sm) & is.null(re.form)){
    b <- extract(object, b_par(sm))[[1]]
    lp <- lp + Z_matrix(sm, newdata) %*% t(b[samples,,drop=FALSE])
  }
  if(transform) lp <- do.call(sm@link, list(lp))
  t(unname(lp))
})

setGeneric("sim_p", function(object, samples, ...) standardGeneric("sim_p"))
