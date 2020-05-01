#' Draw from the posterior predictive distribution
#'
#' Draw from the posterior predictive distribution after fitting a model.
#' You can draw from the posterior of the observed outcome \eqn{y} or
#' the latent unobserved state \eqn{z}.
#'
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param param Either \code{"y"} for the observed outcome or \code{"z"}
#'   for the unobserved latent state
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'   effects") are included. If \code{NA}, they are ignored
#' @param ... Currently ignored
#'
#' @return A matrix of simulations from the posterior predictive distribution.
#'   If \code{param = "z"}, the dimensions are \code{draws} by number of sites
#'   (or sites x primary periods in site-major order for dynamic models). If
#'   \code{param = "y"}, the dimensions are \code{draws} by sites x observations
#'   (or sites x primary periods x observations for dynamic models).
#'
#' @aliases posterior_predict
#' @method posterior_predict ubmsFit
#' @importFrom rstantools posterior_predict
#' @include fit.R
#' @export
setMethod("posterior_predict", "ubmsFit", 
          function(object, param=c("y","z"), draws=NULL, re.form=NULL, ...){
 
  param <- match.arg(param, c("y", "z"))
  nsamp <- nsamples(object)
  samp_inds <- get_samples(object, draws)

  switch(param, 
         "z" = sim_z(object, samples=samp_inds, re.form=re.form),
         "y" = sim_y(object, samples=samp_inds, re.form=re.form))
})


## simulation methods ---------------------------------------------------------

#Simulate the unobserved latent state
setGeneric("sim_z", function(object, ...) standardGeneric("sim_z"))

#Simulate the observed outcome
setGeneric("sim_y", function(object, ...) standardGeneric("sim_y"))


## occu -----------------------------------------------------------------------

#' @include occu.R
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
    psi <- psi_post[unkZ, i]
    qT <- apply(q_post[unkZ,,i],1, prod)
    psi_con <- psi * qT / (psi * qT + (1-psi))
    to_sim <- unkZ[!is.na(psi_con)]
    z_post[to_sim, i] <- rbinom(length(to_sim),1,psi_con[!is.na(psi_con)])
  }

  t(z_post)
})


setMethod("sim_y", "ubmsFitOccu", function(object, samples, re.form, z=NULL, ...){  
  nsamp <- length(samples)
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)

  z <- process_z(object, samples, re.form, z)
  p <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL, 
                samples=samples, re.form=re.form))
  
  zp <- z[rep(1:nrow(z), each=J),] * p

  y_sim <- suppressWarnings(rbinom(M*J*nsamp, 1, as.vector(zp)))  
  matrix(y_sim, nrow=nsamp, ncol=M*J, byrow=TRUE)
})


# occuRN ----------------------------------------------------------------------

#' @include occuRN.R
setMethod("sim_z", "ubmsFitOccuRN", function(object, samples, re.form, ...){

  r_post <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                   samples=samples, re.form=re.form))
  lam_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=re.form))

  M <- nrow(lam_post)
  J <- nrow(r_post) / M

  r_post <- array(r_post, c(J,M,length(samples)))
  r_post <- aperm(r_post, c(2,1,3)) 
  
  y <- getY(object@data)
  Kmin <- apply(y, 1, function(x) max(x, na.rm=TRUE))
  K <- object@call[["K"]]
  K <- ifelse(is.null(K), 20, K)

  t(simz_occuRN(y, lam_post, r_post, K, Kmin, 0:K))
})


setMethod("sim_y", "ubmsFitOccuRN", function(object, samples, re.form, z=NULL, ...){
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)
  
  z <- process_z(object, samples, re.form, z)
  r <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL, 
                samples=samples, re.form=re.form))
  N <- z[rep(1:nrow(z), each=J),]
  p <- as.vector(1 - (1-r)^N)
  
  y_sim <- rep(NA, length(p))
  not_na <- !is.na(p)
  y_sim[not_na] <- rbinom(sum(not_na), 1, p[not_na]) 
  matrix(y_sim, nrow=length(samples), ncol=M*J, byrow=TRUE)
})


# pcount ----------------------------------------------------------------------

#' @include pcount.R
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


# utils -----------------------------------------------------------------------

# Generate z latent state matrix if it isn't provided
process_z <- function(object, samples, re.form, z){
  if(is.null(z)){
    z <- t(sim_z(object, samples=samples, re.form=re.form))
  } else {
    z <- t(z)
  }
  z
}
