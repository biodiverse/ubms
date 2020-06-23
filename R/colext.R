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
#' @param ... Arguments passed to the \code{\link{stan}} call, such as 
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitColext} object describing the model fit.
#'
#' @references MacKenzie DI, Nicholas JD, Hines JE, Knutson MG, Franklin AB.
#'             2003. Ecology 84: 2200-2207.
#'
#' @seealso \code{\link{colext}}, \code{\link{unmarkedMultFrame}}
#' @export
stan_colext <- function(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, 
                        pformula = ~1, data, ...){

  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), "binomial", "binomial", umf@numPrimary)
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "plogis")
  col <- ubmsSubmodel("Colonization", "col", yearlySiteCovs(umf), 
                      gammaformula, "plogis", transition=TRUE)
  ext <- ubmsSubmodel("Extinction", "ext", yearlySiteCovs(umf), 
                      epsilonformula, "plogis", transition=TRUE)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, col, ext, det)

  ubmsFit("colext", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include occu.R 
setClass("ubmsFitColext", contains = "ubmsFitOccu")


#Method for fitted values------------------------------------------------------

#' @include fitted.R
setMethod("sim_fitted", "ubmsFitColext", 
          function(object, submodel, samples, ...){
  if(submodel == "state") return(colext_occ_prob(object, samples))
  callNextMethod(object, submodel, samples, ...)
})

#This can probably be used with the z function below
colext_occ_prob <- function(object, samples){
  
  pp <- sapply(c("state","col","ext"), function(x){
        t(ubms:::sim_lp(object, x, transform=TRUE, newdata=NULL, samples=samples,
                     re.form=NULL))
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


#Goodness-of-fit---------------------------------------------------------------

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
