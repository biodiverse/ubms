#' Fit the N-mixture model of Royle (2004)
#'
#' This function fits the single season N-mixture model of
#' Royle et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and abundance in that order
#' @param data A \code{\link{unmarkedFramePCount}} object
#' @param K Integer upper index of integration for N-mixture. This should be
#'  set high enough so that it does not affect the parameter estimates.
#'  Note that computation time will increase with K.
#' @param mixture Character specifying mixture: "P" is only option currently.
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitPcount} object describing the model fit.
#'
#' @seealso \code{\link{pcount}}, \code{\link{unmarkedFramePCount}}
#' @export
stan_pcount <- function(formula, data, K=NULL, mixture="P", ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist=mixture, K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("pcount", match.call(), data, response, submodels, ...)
}

#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitPcount", contains = "ubmsFit")


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
#' A goodness-of-fit test for N-mixture models based on Pearson's chi-square.
#' @include gof.R
setMethod("gof", "ubmsFitPcount", function(object, draws=NULL, quiet=FALSE, ...){

  samples <- get_samples(object, draws)
  draws <- length(samples)

  lam <- sim_lp(object, "state", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=NULL)
  p <- sim_lp(object, "det", transform=TRUE, newdata=NULL, samples=samples,
                   re.form=NULL)

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
    mb_obs[i] <- Nmix_chisq(object, lam[i,], p[i,])
    object_star@response <- ubmsResponse(ysim[,,i], object@response@y_dist,
                                         object@response@z_dist, max_primary=T)
    mb_sim[i] <- Nmix_chisq(object_star, lam[i,], p[i,])
    if(!quiet) utils::setTxtProgressBar(pb, i)
  }
  if(!quiet) close(pb)
  ubmsGOF("N-mixture Chi-square", data.frame(obs=mb_obs, sim=mb_sim))
})

Nmix_chisq <- function(object, lambda, p){
  J <- object@response@max_obs
  lambda <- lambda[rep(1:length(lambda), each=J)]
  stopifnot(length(lambda) == length(p))
  fit <- lambda*p
  obs <- as_vector(object@response, na.rm=FALSE)
  stopifnot(length(obs) == length(fit))
  sum((obs - fit)^2/fit, na.rm=TRUE)
}

#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitPcount", function(object, samples, re.form, ...){

  p_post <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                   samples=samples, re.form=re.form))
  lam_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=re.form))

  M <- nrow(lam_post)
  J <- nrow(p_post) / M

  p_post <- array(p_post, c(J,M,length(samples)))
  p_post <- aperm(p_post, c(2,1,3))

  y <- getY(object@data)
  K <- object@response@K
  Kmin <- get_Kmin(object@response)[,1]

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
