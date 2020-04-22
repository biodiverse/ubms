#' @export
setGeneric("getZ", function(object, ...) standardGeneric("getZ"))

#' @include occu.R
setMethod("getZ", "ubmsFitOccu", function(object, ...){

  p_post <- predict(object, 'det', summary=FALSE)
  psi_post <- predict(object, 'state', summary=FALSE)
  
  M <- nrow(psi_post)
  J <- nrow(p_post) / M

  nsamples <- dim(psi_post)[2]

  p_post <- array(p_post, c(J,M,nsamples))
  p_post <- aperm(p_post, c(2,1,3)) 


  Zpost <- matrix(NA, M, nsamples) 

  knownZ <- apply(object@data@y, 1, function(x) sum(x, na.rm=T)>0)
  
  Zpost[knownZ,] <- 1

  unkZ <- which(!knownZ)
  q_post <- 1 - p_post

  for (i in 1:nsamples){
    psi <- psi_post[unkZ, i]
    qT <- apply(q_post[unkZ,,i],1, prod)
    psi_con <- psi * qT / (psi * qT + (1-psi))
    Zpost[unkZ, i] <- rbinom(length(unkZ),1,psi_con)
  }

  Zpost


})

#' @include pcount.R
setMethod("getZ", "ubmsFitPcount", function(object, ...){

  p_post <- predict(object, 'det', summary=FALSE)
  lam_post <- predict(object, 'state', summary=FALSE)
  
  M <- nrow(lam_post)
  J <- nrow(p_post) / M
  nsamples <- dim(lam_post)[2]

  p_post <- array(p_post, c(J,M,nsamples))
  p_post <- aperm(p_post, c(2,1,3)) 
  
  y <- getY(object@data)
  Kinfo <- get_K(y, object@call[["K"]])

  getZ_pcount(y, lam_post, p_post, Kinfo$K, Kinfo$Kmin, 0:Kinfo$K)

})
