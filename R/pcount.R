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


#Goodness-of-fit---------------------------------------------------------------


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
  Kinfo <- get_K(y, object@call[["K"]])

  t(simz_pcount(y, lam_post, p_post, Kinfo$K, Kinfo$Kmin, 0:Kinfo$K))
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
