setGeneric("mb_chisq", function(object, ...) standardGeneric("mb_chisq"))

#' @include fit.R
setMethod("mb_chisq", "ubmsFitOccu", function(object, state, p){
  cohorts <- get_cohorts(object)
  dat <- cohorts$data
  inds <- cohorts$inds
  chisq <- 0
  J <- object@response@max_obs
  for (i in 1:length(inds)){
    state_sub <- state[inds[[i]]]
    idx <- c(sapply(inds[[i]], function(k) seq.int((k - 1) * J + 1, length.out=J)))
    p_sub <- p[idx]
    N <- length(state_sub)
    obs <- get_obs_counts(dat[[i]])
    expect <- get_exp_counts(object, obs, state_sub, p_sub)
    chisq <- chisq + sum((obs - expect)^2 / expect, na.rm=TRUE)
    chisq <- chisq + max(0, N - sum(expect, na.rm=TRUE))
  }
  chisq
})

setMethod("mb_chisq", "ubmsFitColext", function(object, state, p){

  T <- object@response@max_primary
  s_resp <- split_response_by_T(object@response)
  s_state <- split_psi_by_T(state, object@response)
  s_p <- split_p_by_T(p, object@response)

  mb_vals <- rep(NA, T)
  for (i in 1:T){
    object@response <- s_resp[[i]]
    mb_vals[i] <- callNextMethod(object, s_state[[i]], s_p[[i]])
  }
  sum(mb_vals)
})

get_cohorts <- function(object){
  yt <- t(object@response)
  na_string <- apply(is.na(yt)*1, 2, paste, collapse="")
  #Unique NA placement excluding sites with all NA
  unq <- unique(na_string)
  unq <- unq[grep("0", unq)]
  inds <- lapply(unq, function(x) which(na_string == x))
  cohorts <- lapply(inds, function(x) yt[,x,drop=FALSE])
  list(inds=inds, data=cohorts)
}

get_obs_counts <- function(yt){
  eh_string <- apply(yt, 2, paste, collapse=" ")
  table(eh_string)
}

setGeneric("get_exp_counts", function(object, obs, ...){
             standardGeneric("get_exp_counts")})

#' @include occu.R
setMethod("get_exp_counts", "ubmsFitOccu", function(object, obs, psi, p){
  eh_obs <- names(obs)
  obs_mat <- sapply(eh_obs, function(x) strsplit(x, " "))
  obs_mat <- suppressWarnings(sapply(obs_mat, as.numeric))
  nd <- apply(obs_mat, 2, function(x) 1 - max(x, na.rm=TRUE))
  counts_expect <- exp_counts_occu(obs_mat, nd, psi, p)[,1]
  names(counts_expect) <- eh_obs
  counts_expect
})

#' @include occuRN.R
setMethod("get_exp_counts", "ubmsFitOccuRN", function(object, obs, lam, p){
  eh_obs <- names(obs)
  obs_mat <- sapply(eh_obs, function(x) strsplit(x, " "))
  obs_mat <- suppressWarnings(sapply(obs_mat, as.numeric))
  Kmin <- apply(obs_mat, 2, function(x) max(x, na.rm=TRUE))
  counts_expect <- exp_counts_occuRN(obs_mat, Kmin, lam, p)[,1]
  names(counts_expect) <- eh_obs
  counts_expect
})

split_response_by_T <- function(resp){
  yt <- t(resp)
  J <- resp@max_obs
  lapply(1:resp@max_primary, function(i){
    cols <- (1:J) + (J*(i-1))
    y <- t(yt[cols,,drop=FALSE])
    ubmsResponse(y, resp@y_dist, resp@z_dist, max_primary=1)
  })
}

split_psi_by_T <- function(psi, resp){
  M <- nrow(resp@y)
  T <- resp@max_primary
  stopifnot(length(psi) == M*T)
  lapply(1:T, function(i){
    psi[seq(i, M*T, by=T)]
  })
}

split_p_by_T <- function(p, resp){
  M <- nrow(resp@y)
  T <- resp@max_primary
  J <- resp@max_obs
  R <- M*T*J
  stopifnot(length(p) == R)
  ind_mat <- matrix(1:R, nrow=J)
  lapply(1:T, function(i){
    cols <- seq(i, M*T, by=T)
    p[as.vector(ind_mat[,cols])]
  })
}
