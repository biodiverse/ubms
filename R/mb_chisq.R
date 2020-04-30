setGeneric("mb_chisq", function(object, ...) standardGeneric("mb_chisq"))

#' @include occu.R
setMethod("mb_chisq", "ubmsFitOccu", function(object, psi, p){
  
  cohorts <- get_cohorts(object)
  dat <- cohorts$data
  inds <- cohorts$inds
  site_idx <- rep(1:nrow(getY(object@data)), each=ncol(getY(object@data)))
  
  chisq <- 0
  for (i in 1:length(inds)){
    psi_sub <- psi[inds[[i]]]
    p_sub <- p[which(site_idx %in% inds[[i]])]
    N <- length(psi_sub)
    obs <- get_obs_counts(dat[[i]])
    expect <- get_exp_counts(object, obs, psi_sub, p_sub)
    chisq <- chisq + sum((obs - expect)^2 / expect, na.rm=TRUE)
    chisq <- chisq + max(0, N - sum(expect, na.rm=TRUE))
  }
  chisq
})

get_cohorts <- function(object){
  umf <- object@data
  yt <- as.vector(t(getY(umf)))
  yt[object["det"]@missing] <- NA
  y <- matrix(yt, nrow=nrow(umf@y), byrow=TRUE)
  na_string <- apply(is.na(y)*1, 1, paste, collapse="")  
  #Unique NA placement excluding sites with all NA
  unq <- unique(na_string)
  unq <- unq[grep("0", unq)]
  umf@y <- y
  inds <- lapply(unq, function(x) which(na_string == x))
  cohorts <- lapply(inds, function(x) umf[x,])
  list(inds=inds, data=cohorts)
}

get_obs_counts <- function(umf){
  eh_string <- apply(getY(umf), 1, paste, collapse=" ")
  table(eh_string)
}


setGeneric("get_exp_counts", function(object, obs, ...){
             standardGeneric("get_exp_counts")})

setMethod("get_exp_counts", "ubmsFitOccu", function(object, obs, psi, p){  
  eh_obs <- names(obs)
  obs_mat <- sapply(eh_obs, function(x) strsplit(x, " "))
  obs_mat <- suppressWarnings(sapply(obs_mat, as.numeric))
  nd <- apply(obs_mat, 2, function(x) 1 - max(x, na.rm=TRUE))  
  counts_expect <- exp_counts_occu(obs_mat, nd, psi, p)[,1]
  names(counts_expect) <- eh_obs
  counts_expect
})
