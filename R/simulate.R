#' @include predict.R
#' @export
setMethod("simulate", "ubmsFit", function(object, nsim=1, seed = NULL, ...){
  
  det_post <- predict(object, 'det', summary=FALSE)
  state_post <- predict(object, 'state', summary=FALSE)
  M <- nrow(state_post)
  J <- nrow(det_post) / M
  nsamples <- dim(state_post)[2]
  
  if(missing(nsim) || nsim >= nsamples){
    samples <- 1:nsamples
  } else {
    samples <- sample(1:nsamples, nsim) 
  }

  z <- rep(0, M)
  zlong <- ylong <- rep(0, M*J)
  ysims <- array(0, c(M,J,length(samples)))
  for (i in 1:length(samples)){
    z <- rbinom(M, 1, state_post[,samples[i]])
    zlong <- rep(z, each=J)
    ylong <- rbinom(M*J, 1, det_post[,samples[i]] * zlong)
    ysims[,,i] <- t(matrix(ylong, J, M))
  }
  drop(ysims)
})
