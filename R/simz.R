#' @include predict.R
setGeneric("sim_z", function(object, ...) standardGeneric("sim_z"))

#' @include occu.R
setMethod("sim_z", "ubmsFitOccu", function(object, samples, ...){  
  nsamples <- length(samples)
  p_post <- predict(object, 'det', summary=FALSE)[,samples,drop=FALSE]
  psi_post <- predict(object, 'state', summary=FALSE)[,samples,drop=FALSE]

  M <- nrow(psi_post)
  J <- nrow(p_post) / M

  p_post <- array(p_post, c(J,M,nsamples))
  p_post <- aperm(p_post, c(2,1,3)) 

  z_post <- matrix(NA, M, nsamples) 

  knownZ <- apply(object@data@y, 1, function(x) sum(x, na.rm=T)>0)
  
  z_post[knownZ,] <- 1

  unkZ <- which(!knownZ)

  q_post <- 1 - p_post

  for (i in 1:nsamples){
    psi <- psi_post[unkZ, i]
    qT <- apply(q_post[unkZ,,i],1, prod)
    psi_con <- psi * qT / (psi * qT + (1-psi))
    to_sim <- unkZ[!is.na(psi_con)]
    z_post[to_sim, i] <- rbinom(length(to_sim),1,psi_con[!is.na(psi_con)])
  }

  z_post
})


#' @include pcount.R
setMethod("sim_z", "ubmsFitPcount", function(object, samples, ...){
  nsamples <- length(samples)
  p_post <- predict(object, 'det', summary=FALSE)[,samples,drop=FALSE]
  lam_post <- predict(object, 'state', summary=FALSE)[,samples,drop=FALSE]
  
  M <- nrow(lam_post)
  J <- nrow(p_post) / M

  p_post <- array(p_post, c(J,M,nsamples))
  p_post <- aperm(p_post, c(2,1,3)) 
  
  y <- getY(object@data)
  Kinfo <- get_K(y, object@call[["K"]])

  simz_pcount(y, lam_post, p_post, Kinfo$K, Kinfo$Kmin, 0:Kinfo$K)
})


#' @include occuRN.R
setMethod("sim_z", "ubmsFitOccuRN", function(object, samples, ...){
  nsamples <- length(samples)
  r_post <- predict(object, 'det', summary=FALSE)[,samples,drop=FALSE]
  lam_post <- predict(object, 'state', summary=FALSE)[,samples,drop=FALSE]
  
  M <- nrow(lam_post)
  J <- nrow(r_post) / M

  r_post <- array(r_post, c(J,M,nsamples))
  r_post <- aperm(r_post, c(2,1,3)) 
  
  y <- getY(object@data)
  Kmin <- apply(y, 1, function(x) max(x, na.rm=TRUE))
  K <- object@call[["K"]]
  K <- ifelse(is.null(K), 20, K)

  simz_occuRN(y, lam_post, r_post, K, Kmin, 0:K)
})
