#' Fit the Occupancy Model of Royle and Nichols (2003)
#'
#' Fit the occupancy model of Royle and Nichols (2003), which relates
#' probability of detection of the species to the number of individuals
#' available for detection at each site.
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFrameOccu}} object
#' @param K Integer upper index of integration for N-mixture. This should be
#'  set high enough so that it does not affect the parameter estimates.
#'  Note that computation time will increase with K.
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
#' @return \code{ubmsFitOccuRN} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(birds)
#' woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)
#' #Add a site covariate
#' siteCovs(woodthrushUMF) <- data.frame(cov1=rnorm(numSites(woodthrushUMF)))
#'
#' (fm_wood <- stan_occuRN(~1~cov1, woodthrushUMF, chains=3, iter=300))
#' }
#'
#' @references Royle JA, Nichols JD. 2003. Estimating abundance from
#'  repeated presence-absence data or point counts. Ecology 84: 777-790.
#'
#' @seealso \code{\link{occuRN}}, \code{\link{unmarkedFrameOccu}}
#' @export
stan_occuRN <- function(formula,
                        data,
                        K=20,
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

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist="P", K=K)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis",
                      prior_intercept_det, prior_coef_det, prior_sigma)
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occuRN", match.call(), data, response, submodels, log_lik, ...)
}


#Fit object--------------------------------------------------------------------

#' @include occu.R
setClass("ubmsFitOccuRN", contains = "ubmsFitOccu")


#Method for fitted values------------------------------------------------------

#' @include fitted.R
#' @importFrom methods callNextMethod
setMethod("sim_fitted", "ubmsFitOccuRN",
          function(object, submodel, samples, ...){
  if(identical(submodel,"det")){
    lp <- sim_lp(object, submodel, transform=TRUE, newdata=NULL,
                 samples=samples, re.form=NULL)
    z <- sim_z(object, samples, re.form=NULL)
    J <- object@response@max_obs
    z <- z[, rep(1:ncol(z), each=J)]
    p <- 1 - (1 - lp)^z
    p[z == 0] <- NA
    return(p)
  }

  callNextMethod(object, submodel, samples, ...)
})


#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitOccuRN", function(object, samples, re.form, ...){

  r_post <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                   samples=samples, re.form=re.form))
  r_post[object["det"]@missing] <- NA
  lam_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=re.form))
  lam_post[object["state"]@missing] <- NA

  M <- nrow(lam_post)
  J <- nrow(r_post) / M

  r_post <- array(r_post, c(J,M,length(samples)))
  r_post <- aperm(r_post, c(2,1,3))

  y <- getY(object@data)
  K <- object@response@K
  Kmin <- apply(y, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE)))

  t(simz_occuRN(y, lam_post, r_post, K, Kmin, 0:K))
})


setMethod("sim_y", "ubmsFitOccuRN", function(object, samples, re.form, z=NULL, ...){
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)

  z <- process_z(object, samples, re.form, z)
  r <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form))
  r[object["det"]@missing] <- NA
  N <- z[rep(1:nrow(z), each=J),]
  p <- as.vector(1 - (1-r)^N)

  y_sim <- rep(NA, length(p))
  not_na <- !is.na(p)
  y_sim[not_na] <- rbinom(sum(not_na), 1, p[not_na])
  matrix(y_sim, nrow=length(samples), ncol=M*J, byrow=TRUE)
})
