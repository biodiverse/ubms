#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of
#' MacKenzie et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link{unmarkedFrameOccu}} object
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccu} object describing the model fit.
#'
#' @seealso \code{\link{occu}}, \code{\link{unmarkedFrameOccu}}
#' @export
stan_occu <- function(formula, data, ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), forms[[2]], "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occu", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitOccu", contains = "ubmsFit")


#Goodness-of-fit---------------------------------------------------------------

#' @describeIn gof
#' Applies the MacKenzie-Bailey chi-square goodness of fit test for
#' ocupancy models (MacKenzie and Bailey 2004).
#' @references MacKenzie, D. I., & Bailey, L. L. (2004). Assessing the
#'  fit of site-occupancy models. Journal of Agricultural, Biological,
#'  and Environmental Statistics, 9(3), 300-318.
#' @include gof.R
setMethod("gof", "ubmsFitOccu", function(object, draws=NULL, quiet=FALSE, ...){

  samples <- get_samples(object, draws)
  draws <- length(samples)

  psi <- sim_state(object, samples)
  p <- sim_lp(object, transform=TRUE, submodel="det", newdata=NULL,
                     samples, re.form=NULL)

  M <- get_n_sites(object@response)
  T <- object@response@max_primary
  R <- T * object@response@max_obs
  ysim <- sim_y(object, samples, re.form=NULL)
  ysim <- array(ysim, c(draws,R,M))
  ysim <- aperm(ysim, c(3,2,1))

  mb_obs <- mb_sim <- rep(NA, draws)
  if(!quiet) pb <- utils::txtProgressBar(min = 0, max = draws, style = 3)
  object_star <- object
  for (i in 1:draws){
    mb_obs[i] <- mb_chisq(object, psi[i,], p[i,])
    object_star@response <- ubmsResponse(ysim[,,i], object@response@y_dist,
                                         object@response@z_dist, max_primary=T)
    #mb_chisq handles replicating NAs
    mb_sim[i] <- mb_chisq(object_star, psi[i,], p[i,])
    if(!quiet) utils::setTxtProgressBar(pb, i)
  }
  if(!quiet) close(pb)
  ubmsGOF("MacKenzie-Bailey Chi-square", data.frame(obs=mb_obs, sim=mb_sim))
})

setGeneric("sim_state", function(object, ...) standardGeneric("sim_state"))

setMethod("sim_state", "ubmsFitOccu", function(object, samples, ...){
  sim_lp(object, transform=TRUE, submodel="state", newdata=NULL,
         samples, re.form=NULL)
})

#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitOccu", function(object, samples, re.form, ...){

  p_post <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=re.form))
  psi_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                       samples=samples, re.form=re.form))

  M <- nrow(psi_post)
  J <- nrow(p_post) / M
  nsamp <- length(samples)

  p_post <- array(p_post, c(J,M,nsamp))
  p_post <- aperm(p_post, c(2,1,3))

  z_post <- matrix(NA, M, nsamp)

  knownZ <- apply(object@data@y, 1, function(x) sum(x, na.rm=T)>0)

  z_post[knownZ,] <- 1

  unkZ <- which(!knownZ)

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
  M <- get_n_sites(object@response)
  J <- object@response@max_obs
  T <- object@response@max_primary

  z <- process_z(object, samples, re.form, z)
  p <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form))

  zp <- z[rep(1:nrow(z), each=J),,drop=FALSE] * p

  y_sim <- suppressWarnings(rbinom(M*J*T*nsamp, 1, as.vector(zp)))
  matrix(y_sim, nrow=nsamp, ncol=M*J*T, byrow=TRUE)
})
