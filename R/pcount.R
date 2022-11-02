#' Fit the N-mixture model of Royle (2004)
#'
#' This function fits the single season N-mixture model of
#' Royle et al. (2004).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFramePCount}} object
#' @param K Integer upper index of integration for N-mixture. This should be
#'  set high enough so that it does not affect the parameter estimates.
#'  Note that computation time will increase with K.
#' @param mixture Character specifying mixture: "P" is only option currently.
#' @param prior_intercept_state Prior distribution for the intercept of the
#'  state (abundance) model; see \code{?priors} for options
#' @param prior_coef_state Prior distribution for the regression coefficients of
#'  the state model
#' @param prior_intercept_det Prior distribution for the intercept of the
#'  detection probability model
#' @param prior_coef_det Prior distribution for the regression coefficients of
#'  the detection model
#' @param prior_sigma Prior distribution on random effect standard deviations
#' @param log_lik If \code{TRUE}, Stan will save pointwise log-likelihood values
#'  in the output. This can greatly increase the size of the model. If
#'  \code{FALSE}, the values are calculated post-hoc from the posteriors
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitPcount} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(mallard)
#' mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs=mallard.site)
#'
#' (fm_mallard <- stan_pcount(~1~elev+forest, mallardUMF, K=30,
#'                            chains=3, iter=300))
#' }
#'
#' @references Royle JA. 2004. N-mixture models for estimating populaiton size
#'  from spatially replicated counts. Biometrics 60: 105-108.
#'
#' @seealso \code{\link{pcount}}, \code{\link{unmarkedFramePCount}}
#' @export
stan_pcount <- function(formula,
                        data,
                        K=NULL,
                        mixture="P",
                        prior_intercept_state = normal(0, 5),
                        prior_coef_state = normal(0, 2.5),
                        prior_intercept_det = logistic(0, 1),
                        prior_coef_det = logistic(0, 1),
                        prior_sigma = gamma(1, 1),
                        log_lik = TRUE,
                        ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)

  if(has_spatial(forms)){
    split_umf <- extract_missing_sites(umf)
    umf <- split_umf$umf
    state <- ubmsSubmodelSpatial("Abundance", "state", siteCovs(umf), forms[[2]],
                                 "exp", prior_intercept_state, prior_coef_state,
                                 prior_sigma,
                                 split_umf$sites_augment, split_umf$data_aug)
  } else {
    state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp",
                          prior_intercept_state, prior_coef_state, prior_sigma)
  }

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist=mixture, K=K)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis",
                      prior_intercept_det, prior_coef_det, prior_sigma)
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("pcount", match.call(), data, response, submodels, log_lik, ...)
}

#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitPcount", contains = "ubmsFitAbun")


#Method for fitted values------------------------------------------------------

#' @include fitted.R
setMethod("sim_fitted", "ubmsFitPcount",
          function(object, submodel, samples, ...){
  if(identical(submodel,"det")){
    p <- sim_lp(object, submodel, transform=TRUE, newdata=NULL,
                 samples=samples, re.form=NULL)
    z <- sim_z(object, samples, re.form=NULL)
    J <- object@response@max_obs
    z <- z[, rep(1:ncol(z), each=J)]
    return(z * p)
  }

  callNextMethod(object, submodel, samples, ...)
})


#Goodness-of-fit---------------------------------------------------------------

#' @describeIn gof
#' A goodness-of-fit test for N-mixture type models based on Pearson's chi-square.
#' @include gof.R
setMethod("gof", "ubmsFitAbun", function(object, draws=NULL, quiet=FALSE, ...){
  sim_gof(object, draws, Nmix_chisq, "N-mixture Chi-square", quiet)
})


#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitPcount", function(object, samples, re.form, ...){

  p_post <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                   samples=samples, re.form=re.form))
  p_post[object["det"]@missing] <- NA
  lam_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=re.form))
  lam_post[object["state"]@missing] <- NA

  M <- nrow(lam_post)
  J <- nrow(p_post) / M

  p_post <- array(p_post, c(J,M,length(samples)))
  p_post <- aperm(p_post, c(2,1,3))

  y <- getY(object@data)
  K <- object@response@K
  Kmin <- apply(y, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE)))

  t(simz_pcount(y, lam_post, p_post, K, Kmin, 0:K))
})


setMethod("sim_y", "ubmsFitPcount", function(object, samples, re.form, z=NULL, ...){
  nsamples <- length(samples)
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)

  z <- process_z(object, samples, re.form, z)
  p <- sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form)
  p <- as.vector(t(p))
  N <- as.vector(z[rep(1:nrow(z), each=J),])

  y_sim <- rep(NA, length(p))
  not_na <- !is.na(p) & !is.na(N)
  y_sim[not_na] <- rbinom(sum(not_na), N[not_na], p[not_na])
  matrix(y_sim, nrow=length(samples), ncol=M*J, byrow=TRUE)
})
