#' @include predict.R
setGeneric("sim_y", function(object, z, samples, ...) standardGeneric("sim_y"))

#' @include occu.R
setMethod("sim_y", "ubmsFitOccu", function(object, z, samples, ...){
  nsamples <- length(samples)
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)

  p <- predict(object, "det", summary=FALSE)[,samples,drop=FALSE]
  zp <- z[rep(1:nrow(z), each=J),] * p

  y_sim <- suppressWarnings(rbinom(M*J*nsamples, 1, as.vector(zp)))
  y_sim <- array(y_sim, c(J, M, nsamples))
  aperm(y_sim, c(2,1,3))
})


#' @include pcount.R
setMethod("sim_y", "ubmsFitPcount", function(object, z, samples, ...){  
  nsamples <- length(samples)
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)
  
  p <- as.vector(predict(object, "det", summary=FALSE)[,samples,drop=FALSE])  
  N <- as.vector(z[rep(1:nrow(z), each=J),])

  y_sim <- rep(NA, length(p))
  not_na <- !is.na(p) & !is.na(N)
  y_sim[not_na] <- rbinom(sum(not_na), N[not_na], p[not_na])
  y_sim <- array(y_sim, c(J, M, nsamples))
  aperm(y_sim, c(2,1,3))
})


#' @include occuRN.R
setMethod("sim_y", "ubmsFitOccuRN", function(object, z, samples, ...){
  nsamples <- length(samples)
  y <- getY(object@data)
  M <- nrow(y)
  J <- ncol(y)

  r <- predict(object, "det", summary=FALSE)[,samples,drop=FALSE]
  N <- z[rep(1:nrow(z), each=J),]
  p <- as.vector(1 - (1-r)^N)
  
  y_sim <- rep(NA, length(p))
  not_na <- !is.na(p)
  y_sim[not_na] <- rbinom(sum(not_na), 1, p[not_na])
  y_sim <- array(y_sim, c(J, M, nsamples))
  aperm(y_sim, c(2,1,3))
})
