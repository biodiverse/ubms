#' Fit the MacKenzie et al. (2003) Dynamic Occupancy Model
#'
#' This function fits the dynamic occupancy model of
#' MacKenzie et al. (2003).
#'
#' @param psiformula Right-hand sided formula for the initial probability of
#'                   occupancy at each site
#' @param gammaformula Right-hand sided formula for colonization probability
#' @param epsilonformula Right-hand sided formula for extinction probability
#' @param pformula Right-hand sided formula for detection probability
#' @param data A \code{\link{unmarkedMultFrame}} object
#' @param prior_intercept_psi Prior distribution for the intercept of the
#'  psi (initial occupancy probability) model; see \code{?priors} for options
#' @param prior_coef_psi Prior distribution for the regression coefficients of
#'  the psi model
#' @param prior_intercept_gamma Prior distribution on intercept for
#'  colonization probability
#' @param prior_coef_gamma Prior distribution on regression coefficients for
#'  colonization probability
#' @param prior_intercept_eps Prior distribution on intercept for
#'  extinction probability
#' @param prior_coef_eps Prior distribution on regression coefficients for
#'  extinction probability
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
#' @return \code{ubmsFitColext} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(frogs)
#' umf <- formatMult(masspcru)
#' umf@y[umf@y > 1] <- 1 #convert counts to presence/absence
#' umf <- umf[1:100,] #Use only 100 sites
#'
#' fit_frog <- stan_colext(~1, ~1, ~1, ~1, umf, chains=3, iter=300)
#' }
#'
#' @references MacKenzie DI, Nicholas JD, Hines JE, Knutson MG, Franklin AB.
#'             2003. Ecology 84: 2200-2207.
#'
#' @seealso \code{\link{colext}}, \code{\link{unmarkedMultFrame}}
#' @export
stan_colext <- function(psiformula = ~1,
                        gammaformula = ~1,
                        epsilonformula = ~1,
                        pformula = ~1,
                        data,
                        prior_intercept_psi = logistic(0, 1),
                        prior_coef_psi = logistic(0, 1),
                        prior_intercept_gamma = logistic(0, 1),
                        prior_coef_gamma = logistic(0, 1),
                        prior_intercept_eps = logistic(0, 1),
                        prior_coef_eps = logistic(0, 1),
                        prior_intercept_det = logistic(0, 1),
                        prior_coef_det = logistic(0, 1),
                        prior_sigma = gamma(1, 1),
                        log_lik = TRUE,
                        ...){

  umf <- process_umf(data)

  has_spatial(list(psi=psiformula, gamma=gammaformula, eps=epsilonformula,
                   p=pformula), support=FALSE)

  response <- ubmsResponse(getY(umf), "binomial", "binomial", umf@numPrimary)
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "plogis",
                       prior_intercept_psi, prior_coef_psi, prior_sigma)
  col <- ubmsSubmodelTransition("Colonization", "col", yearlySiteCovs(umf),
                      gammaformula, "plogis", T=umf@numPrimary,
                      prior_intercept_gamma, prior_coef_gamma, prior_sigma)
  ext <- ubmsSubmodelTransition("Extinction", "ext", yearlySiteCovs(umf),
                      epsilonformula, "plogis", T=umf@numPrimary,
                      prior_intercept_eps, prior_coef_eps, prior_sigma)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), pformula, "plogis",
                      prior_intercept_det, prior_coef_det, prior_sigma)
  submodels <- ubmsSubmodelList(state, col, ext, det)

  ubmsFit("colext", match.call(), data, response, submodels, log_lik, ...)
}


#Output object-----------------------------------------------------------------

#' @include occu.R
setClass("ubmsFitColext", contains = "ubmsFitOccu")


#Function for projected psi across years---------------------------------------

#' Projected Occupancy Trajectories
#'
#' Generate posterior draws of occupancy probability for all sites and primary
#' periods, i.e. the projected trajectory (Weir et al. 2009).
#'
#' @param object A fitted dynamic occupancy model of class inheriting \code{ubmsFit}
#' @param draws Number of draws from the posterior to use in the check
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'  effects") are included. If \code{NA}, they are ignored
#' @param ... Currently ignored
#'
#' @return  A matrix of occupancy values from the posterior predictive distribution.
#'   The dimensions are \code{draws} by number of sites x primary periods
#'   in site-major order.
#'
#' @references Weir LA, Fiske IJ, Royle J. 2009. Trends in Anuran
#'   Occupancy from Northeastern States of the North American Amphibian
#'   Monitoring Program. Herpetological Conservation and Biology.
#'   4: 389-402.
#'
#' @aliases projected
#' @seealso \code{\link{stan_colext}}
#' @export
setGeneric("projected", function(object, ...) standardGeneric("projected"))

#' @rdname projected
#' @method projected ubmsFitColext
setMethod("projected", "ubmsFitColext", function(object, draws=NULL,
                                                 re.form=NULL, ...){
  samples <- get_samples(object, draws)
  sim_projected(object, samples, re.form)
})

sim_projected <- function(object, samples, re.form){

  pp <- sapply(c("state","col","ext"), function(x){
        t(sim_lp(object, x, transform=TRUE, newdata=NULL, samples=samples,
                     re.form=re.form))
        }, simplify=FALSE)

  M <- nrow(pp$state)
  T <- object@response@max_primary
  nsamp <- length(samples)

  inv_psi <- 1 - pp$state

  out <- matrix(NA, M*T, nsamp)

  for (s in 1:nsamp){
    phi_raw <- cbind(1-pp$col[,s], pp$ext[,s], pp$col[,s], 1-pp$ext[,s])
    idx <- tidx <- 1
    for (i in 1:M){
      for (t in 1:T){
        if(t==1) zprob <- c(inv_psi[i,s], pp$state[i,s])
        else zprob <- zprob %*% matrix(phi_raw[tidx,], nrow=2)
        out[idx,s] <- zprob[2]
        idx <- idx + 1
        if(t>1) tidx <- tidx + 1
      }
    }
  }
  t(out)
}


#Function for turnover---------------------------------------------------------

#' Turnover Probability
#'
#' Generate posterior draws of turnover probability from dynamic occupancy models.
#' Turnover is calculated for each site and each primary period after the first.
#'
#' @param object A fitted dynamic occupancy model of class inheriting \code{ubmsFit}
#' @param draws Number of draws from the posterior to use in the check
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'  effects") are included. If \code{NA}, they are ignored
#' @param ... Currently ignored
#'
#' @return  A matrix of turnover values from the posterior predictive distribution.
#'   The dimensions are \code{draws} by number of sites x (primary periods - 1)
#'   in site-major order.
#'
#' @aliases turnover
#' @seealso \code{\link{stan_colext}}
#' @export
setGeneric("turnover", function(object, ...) standardGeneric("turnover"))

#' @rdname turnover
#' @method turnover ubmsFitColext
setMethod("turnover", "ubmsFitColext", function(object, draws, re.form=NULL, ...){
  samples <- get_samples(object, draws)
  sim_turnover(object, samples, re.form)
})

sim_turnover <- function(object, samples, re.form){
  M <- nrow(object@response@y)
  T <- object@response@max_primary
  nsamp <- length(samples)

  psi_hat <- sim_projected(object, samples, re.form=re.form)
  col <- sim_lp(object, "col", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form)

  out <- matrix(NA, M*(T-1), nsamp)
  for (s in 1:nsamp){
    idx <- 1
    psi_mat <- matrix(psi_hat[s,], nrow=T)
    inv_psi_mat <- 1 - psi_mat
    col_mat <- matrix(col[s,], nrow=(T-1))
    for (i in 1:M){
      for (t in 2:T){
        out[idx, s] <- col_mat[(t-1),i]*inv_psi_mat[(t-1),i] / psi_mat[t,i]
        idx <- idx + 1
      }
    }
  }
  t(out)
}

#Method for fitted values------------------------------------------------------

#' @include fitted.R
setMethod("sim_fitted", "ubmsFitColext",
          function(object, submodel, samples, ...){
  if(submodel == "state") return(sim_projected(object, samples, re.form=NULL))
  callNextMethod(object, submodel, samples, ...)
})


#Goodness-of-fit---------------------------------------------------------------

#' @include posterior_linpred.R
setMethod("sim_state", "ubmsFitColext", function(object, samples, ...){
  sim_projected(object, samples, re.form=NULL)
})

#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitColext", function(object, samples, re.form, ...){

  pp <- sapply(submodel_types(object), function(x){
        t(sim_lp(object, x, transform=TRUE, newdata=NULL, samples=samples,
                     re.form=re.form))
        }, simplify=FALSE)

  M <- nrow(pp$state)
  T <- object@response@max_primary
  J <- object@response@max_obs
  nsamp <- length(samples)

  yt <- matrix(t(object@response), ncol=M*T)
  knownZ <- apply(yt, 2, function(x) sum(x, na.rm=TRUE) > 0)

  q <- 1 - pp$det
  inv_psi <- 1 - pp$state
  z_post <- matrix(NA, M*T, nsamp)

  update_zprob <- function(zprob, psub){
    if(any(!is.na(psub))){
      qT <- prod(psub, na.rm=TRUE)
      zcon <- zprob[2] * qT / (zprob[2] * qT + zprob[1])
      zprob <- c(1-zcon, zcon)
    }
    zprob
  }

  for (s in 1:nsamp){
    phi_raw <- cbind(1-pp$col[,s], pp$ext[,s], pp$col[,s], 1-pp$ext[,s])
    idx <- tidx <- pidx <- 1
    for (i in 1:M){
      for (t in 1:T){
        if(knownZ[idx]){
          z_post[idx, s] <- 1
          zprob <- c(0,1)
        } else {
          if(t==1) zprob <- c(inv_psi[i,s], pp$state[i,s])
          else zprob <- zprob %*% matrix(phi_raw[tidx,], nrow=2)
          zprob <- update_zprob(zprob, q[pidx:(pidx+J-1), s])
          z_post[idx, s] <- rbinom(1, 1, zprob[2])
        }
        idx <- idx + 1
        pidx <- pidx + J
        if(t>1) tidx <- tidx + 1
      }
    }
  }
  t(z_post)
})
