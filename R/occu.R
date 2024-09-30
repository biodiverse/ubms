#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of
#' MacKenzie et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link[unmarked]{unmarkedFrameOccu}} object
#' @param prior_intercept_state Prior distribution for the intercept of the
#'  state (occupancy probability) model; see \code{?priors} for options
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
#' @param ... Arguments passed to the \code{\link[rstan]{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccu} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(frogs)
#' pferUMF <- unmarkedFrameOccu(pfer.bin)
#'
#' #Add some covariates
#' siteCovs(pferUMF) <- data.frame(cov1=rnorm(numSites(pferUMF)))
#'
#' #Fit model
#' (fm <- stan_occu(~1~cov1, pferUMF, chains=3, iter=300))
#' }
#'
#' @references MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA,
#'  Langtimm CA. 2002. Estimating site occupancy rates when detection
#'  probabilities are less than one. Ecology 83: 2248-2255.
#'
#' @seealso \code{\link[unmarked]{occu}}, \code{\link[unmarked]{unmarkedFrameOccu}}
#' @include fit.R
#' @export
stan_occu <- function(formula,
                      data,
                      prior_intercept_state = logistic(0, 1),
                      prior_coef_state = logistic(0, 1),
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
    state <- ubmsSubmodelSpatial("Occupancy", "state", siteCovs(umf), forms[[2]],
                                 "plogis", prior_intercept_state, prior_coef_state, prior_sigma,
                                 split_umf$sites_augment, split_umf$data_aug)

  } else {
    state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), forms[[2]],
                          "plogis", prior_intercept_state, prior_coef_state, prior_sigma)
  }

  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis",
                      prior_intercept_det, prior_coef_det, prior_sigma)
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occu", match.call(), data, response, submodels, log_lik, ...)
}


#Goodness-of-fit---------------------------------------------------------------

#' @describeIn gof
#' Applies the MacKenzie-Bailey chi-square goodness of fit test for
#' ocupancy models (MacKenzie and Bailey 2004).
#' @references MacKenzie, D. I., & Bailey, L. L. (2004). Assessing the
#'  fit of site-occupancy models. Journal of Agricultural, Biological,
#'  and Environmental Statistics, 9(3), 300-318.
#' @include gof.R
setMethod("gof", "ubmsFitOccu", function(object, draws=NULL, quiet=FALSE, ...){
  sim_gof(object, draws, mb_chisq, "MacKenzie-Bailey Chi-square", quiet)
})


#Methods to simulate posterior predictive distributions------------------------

setGeneric("knownZ", function(object, ...) standardGeneric("knownZ"))

setMethod("knownZ", "ubmsFitOccu", function(object, ...){
  apply(object@data@y, 1, function(x) sum(x, na.rm=T)>0)
})

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitOccu", function(object, samples, re.form, ...){

  p_post <- t(sim_p(object, samples))
  psi_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                       samples=samples, re.form=re.form))
  known_z <- knownZ(object)
  if(has_spatial(object["state"])){
    warning("Output only includes sites that were sampled at least once", call.=FALSE)
    sites_noaug <- !object["state"]@sites_aug
    psi_post <- psi_post[sites_noaug,,drop=FALSE]
    known_z <- known_z[sites_noaug]
  }
  psi_post[object["state"]@missing] <- NA

  M <- nrow(psi_post)
  J <- nrow(p_post) / M
  nsamp <- length(samples)

  p_post <- array(p_post, c(J,M,nsamp))
  p_post <- aperm(p_post, c(2,1,3))

  z_post <- matrix(NA, M, nsamp)

  z_post[known_z,] <- 1
  unkZ <- which(!known_z)

  q_post <- 1 - p_post

  for (i in 1:nsamp){
    psi <- psi_post[unkZ, i, drop=FALSE]
    qT <- apply(q_post[unkZ,,i,drop=FALSE], 1,
                function(x) ifelse(all(is.na(x)), NA, prod(x, na.rm=TRUE)))
    psi_con <- psi * qT / (psi * qT + (1-psi))
    to_sim <- unkZ[!is.na(psi_con)]
    z_post[to_sim, i] <- rbinom(length(to_sim),1,psi_con[!is.na(psi_con)])
  }

  t(z_post)
})


setMethod("sim_y", "ubmsFitOccu", function(object, samples, re.form, z=NULL, ...){
  nsamp <- length(samples)
  z <- process_z(object, samples, re.form, z)
  p <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form))
  p[object@response@missing] <- NA

  T <- object@response@max_primary
  M <- nrow(z) / T
  J <- object@response@max_obs

  zp <- z[rep(1:nrow(z), each=J),,drop=FALSE] * p

  y_sim <- suppressWarnings(rbinom(M*J*T*nsamp, 1, as.vector(zp)))
  matrix(y_sim, nrow=nsamp, ncol=M*J*T, byrow=TRUE)
})
