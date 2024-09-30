#' Fit Time-to-detection Occupancy Models
#'
#' Fit time-to-detection occupancy models of Garrard et al.
#' (2008, 2013). Time-to-detection can be modeled with either an exponential
#' or Weibull distribution.
#'
#' @param psiformula Right-hand sided formula for the initial probability of
#'  occupancy at each site.
#' @param gammaformula Right-hand sided formula for colonization probability.
#'  Currently ignored as dynamic models are not yet supported.
#' @param epsilonformula Right-hand sided formula for extinction probability.
#'  Currently ignored as dynamic models are not yet supported.
#' @param detformula Right-hand sided formula for mean time-to-detection.
#' @param data \code{unmarkedFrameOccuTTD} object that supplies the data
#'  (see \code{\link[unmarked]{unmarkedFrameOccuTTD}}).
#' @param ttdDist Distribution to use for time-to-detection; either
#'  \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
#'  which adds an additional shape parameter \eqn{k}.
#' @param linkPsi Link function for the occupancy model. Only option is
#'  \code{"logit"} for now, in the future \code{"cloglog"}
#'  will be supported for the complimentary log-log link.
#' @param prior_intercept_state Prior distribution for the intercept of the
#'  state (occupancy probability) model; see \code{?priors} for options
#' @param prior_coef_state Prior distribution for the regression coefficients of
#'  the state model
#' @param prior_intercept_det Prior distribution for the intercept of the
#'  time-to-detection model
#' @param prior_coef_det Prior distribution for the regression coefficients of
#'  the time-to-detection model
#' @param prior_intercept_shape Prior distribution for the intercept of the
#'  shape parameter (i.e., log(shape)) for Weibull TTD models
#' @param prior_sigma Prior distribution on random effect standard deviations
#' @param log_lik If \code{TRUE}, Stan will save pointwise log-likelihood values
#'  in the output. This can greatly increase the size of the model. If
#'  \code{FALSE}, the values are calculated post-hoc from the posteriors
#' @param ... Arguments passed to the \code{\link[rstan]{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccuTTD} object describing the model fit.
#'
#' @examples
#' \donttest{
#' #Simulate data
#' N <- 500; J <- 1
#' scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
#'                     forest=runif(N,0,1),
#'                     wind=runif(N,0,1))
#' beta_psi <- c(-0.69, 0.71, -0.5)
#' psi <- plogis(cbind(1, scovs$elev, scovs$forest) %*% beta_psi)
#' z <- rbinom(N, 1, psi)
#'
#' Tmax <- 10 #Same survey length for all observations
#' beta_lam <- c(-2, -0.2, 0.7)
#' rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
#' ttd <- rexp(N, rate)
#' ttd[z==0] <- Tmax #Censor at unoccupied sites
#' ttd[ttd>Tmax] <- Tmax #Censor when ttd was greater than survey length
#'
#' #Build unmarkedFrame
#' umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)
#'
#' #Fit model
#' (fit <- stan_occuTTD(psiformula=~elev+forest, detformula=~elev+wind,
#'                      data=umf, chains=3, iter=300))
#' }
#'
#' @references Garrard, G.E., Bekessy, S.A., McCarthy, M.A. and Wintle, B.A. 2008.
#' When have we looked hard enough? A novel method for setting minimum survey effort
#' protocols for flora surveys. Austral Ecology 33: 986-998.
#'
#' Garrard, G.E., McCarthy, M.A., Williams, N.S., Bekessy, S.A. and Wintle,
#' B.A. 2013. A general model of detectability using species traits. Methods in
#' Ecology and Evolution 4: 45-52.
#'
#' Kery, Marc, and J. Andrew Royle. 2016. Applied Hierarchical Modeling in
#' Ecology, Volume 1. Academic Press.
#'
#' @seealso \code{\link[unmarked]{occuTTD}}, \code{\link[unmarked]{unmarkedFrameOccuTTD}}
#' @include fit.R
#' @export
stan_occuTTD <- function(psiformula=~1,
                         gammaformula=~1,
                         epsilonformula=~1,
                         detformula=~1,
                         data,
                         ttdDist=c("exp", "weibull"),
                         linkPsi=c("logit"),
                         prior_intercept_state = logistic(0, 1),
                         prior_coef_state = logistic(0, 1),
                         prior_intercept_det = normal(0, 5),
                         prior_coef_det = normal(0, 2.5),
                         prior_intercept_shape = normal(0,2.5),
                         prior_sigma = gamma(1, 1),
                         log_lik = TRUE,
                         ...){

  if(data@numPrimary > 1) stop("Dynamic models not yet supported", call.=FALSE)
  umf <- process_umf(data)
  ttdDist <- match.arg(ttdDist)
  linkPsi <- match.arg(linkPsi)

  forms <- list(state=psiformula, det=detformula)
  if(has_spatial(forms)){
    split_umf <- extract_missing_sites(umf)
    umf <- split_umf$umf
    state <- ubmsSubmodelSpatial("Occupancy", "state", siteCovs(umf), psiformula,
                                 "plogis", prior_intercept_state, prior_coef_state,
                                 prior_sigma,
                                 split_umf$sites_augment, split_umf$data_aug)
  } else {
    state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "plogis",
                          prior_intercept_state, prior_coef_state, prior_sigma)
  }

  response <- ubmsResponseOccuTTD(umf, ttdDist)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), detformula, "exp",
                      prior_intercept_det, prior_coef_det, prior_sigma)

  shape <- placeholderSubmodel("shape")
  if(ttdDist=="weibull"){
    shape <- ubmsSubmodelScalar("Shape", "shape", "exp", prior_intercept_shape)
  }

  submodels <- ubmsSubmodelList(state, det, shape)

  ubmsFit("occuTTD", match.call(), data, response, submodels, log_lik, ...)
}


#Output object-----------------------------------------------------------------

#' @include occu.R
setClass("ubmsFitOccuTTD", contains = "ubmsFitOccu")


#Response class----------------------------------------------------------------

#' @include response.R
setClass("ubmsResponseOccuTTD", contains="ubmsResponse",
         slots=c(surveyLength="matrix"))

ubmsResponseOccuTTD <- function(umf, ttdDist){
  out <- ubmsResponse(getY(umf), ttdDist, "binomial", umf@numPrimary,
                      ceiling(max(getY(umf), na.rm=TRUE)))
  out <- as(out, "ubmsResponseOccuTTD")
  out@surveyLength <- umf@surveyLength
  out
}

#Bundle delta (indicator if observation occurred before max time) for Stan
#Also include y as aux2 (can't use normal y slot because y is not an integer)

#' @include inputs.R
setMethod("get_auxiliary_data", "ubmsResponseOccuTTD", function(object, ...){
  yt <- t(object)
  delta <- ifelse(yt < t(object@surveyLength), 1, 0)
  keep <- apply(delta, 2, function(x) !all(is.na(x)))
  delta <- delta[,keep,drop=FALSE]
  delta <- as.vector(delta)
  delta <- delta[!is.na(delta)]
  y <- as_vector(object, na.rm=TRUE) #Need to pass to aux because it's not integer
  stopifnot(length(delta) == length(y))
  list(aux1=delta,
       aux2=y,
       aux3=numeric(0),
       n_aux1=length(delta), n_aux2=length(y), n_aux3=0)
})

setMethod("get_stan_data", "ubmsResponseOccuTTD", function(object, ...){
  out <- methods::callNextMethod(object, ...)
  out$y <- out$aux1
  out$Kmin[] <- 0
  out
})


#Get detection probability-----------------------------------------------------

#' @include posterior_linpred.R
setMethod("sim_p", "ubmsFitOccuTTD", function(object, samples, ...){
  resp <- object@response
  tdist <- resp@y_dist
  tmax <- as.vector(t(resp@surveyLength))
  lam <- sim_lp(object, "det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=NULL)
  lam[,object["det"]@missing] <- NA

  if(tdist == "exp"){
    out <- sapply(1:length(samples), function(i){
      stats::pexp(tmax, lam[i,])
    })
  } else if(tdist == "weibull"){
    shape <- t(sim_lp(object, "shape", transform=TRUE, newdata=NULL,
                    samples=samples, re.form=NULL))
    out <- sapply(1:length(samples), function(i){
      stats::pweibull(tmax, shape[i], 1/lam[i,])
    })
  }

  t(out)
})


#Method for fitted values------------------------------------------------------

#' @include fitted.R
setMethod("sim_fitted", "ubmsFitOccuTTD",
          function(object, submodel, samples, ...){
  if(identical(submodel,"det")){
    lam <- sim_lp(object, submodel, transform=TRUE, newdata=NULL,
                 samples=samples, re.form=NULL)
    z <- sim_z(object, samples, re.form=NULL)
    J <- object@response@max_obs
    z <- z[, rep(1:ncol(z), each=J)]
    lam[z==0] <- NA
    return(1/lam)
  }

  callNextMethod(object, submodel, samples, ...)
})

#Methods to simulate posterior predictive distributions------------------------

setMethod("knownZ", "ubmsFitOccuTTD", function(object, ...){
  if(object@response@max_obs > 1){
    warning("Posteriors for y/z may be biased when there are multiple observations per site", call.=FALSE)
  }
  ybin <- ifelse(object@data@y < object@data@surveyLength,1,0)
  apply(ybin, 1, function(x) sum(x, na.rm=T)>0)
})

#' @include posterior_predict.R
setMethod("sim_y", "ubmsFitOccuTTD", function(object, samples, re.form, z=NULL, ...){
  nsamp <- length(samples)
  M <- nrow(object@response@y)
  J <- object@response@max_obs
  T <- object@response@max_primary
  tmax <- matrix(rep(object@response@surveyLength, nsamp), nrow=nsamp, byrow=T)

  z <- process_z(object, samples, re.form, z)
  z <- t(z[rep(1:nrow(z), each=J),,drop=FALSE])
  lam <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form))
  lam[object@response@missing] <- NA

  y_sim <- suppressWarnings(stats::rexp(M*J*T*nsamp, as.vector(lam)))
  y_sim[is.nan(y_sim)] <- NA
  y_sim <- matrix(y_sim, nrow=nsamp, ncol=M*J*T, byrow=TRUE)
  y_sim[z==0&!is.na(y_sim)] <- tmax[z==0&!is.na(y_sim)]
  y_sim[!is.na(y_sim) & y_sim>tmax] <- tmax[!is.na(y_sim) & y_sim>tmax]
  y_sim
})

#Goodness-of-fit---------------------------------------------------------------

#' @include gof.R
setMethod("sim_gof", "ubmsFitOccuTTD", function(object, ...){
  stop("No goodness-of-fit test for this model type available", call.=FALSE)
})
