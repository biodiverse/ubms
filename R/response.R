setClass("ubmsResponse",
  slots = c(
    y = "matrix",
    y_dist = "character",
    z_dist = "character",
    max_primary = "numeric",
    max_obs = "numeric",
    K = "numeric",
    missing = "logical"
  ),
  prototype = list(
    y = matrix(0),
    y_dist = NA_character_,
    z_dist = NA_character_,
    max_primary = 1,
    max_obs = 0,
    K = 1,
    missing = logical(0)
  )
)

ubmsResponse <- function(y, y_dist, z_dist, max_primary = 1, K=NULL){
  stopifnot(inherits(y, "matrix"))
  out <- new("ubmsResponse", y = y, y_dist= y_dist, z_dist = z_dist,
          max_primary = max_primary, max_obs = get_max_obs(y, max_primary))
  out@missing <- is.na(as_vector(out))
  out@K <- get_K(out, K)
  out
}

get_max_obs <- function(y, max_primary){
  if(ncol(y) %% max_primary != 0){
    stop("Primary periods must have same number of secondary periods",
         call.=FALSE)
  }
  ncol(y) / max_primary
}

setGeneric("get_K", function(object, K=NULL, ...) standardGeneric("get_K"))

setMethod("get_K", "ubmsResponse", function(object, K=NULL){
  ymax <- max(object@y, na.rm=TRUE)
  if(is.null(K)) return(ymax + 20)
  if(K < ymax) stop("K must be larger than max y value", call.=FALSE)
  K
})

setMethod("t", "ubmsResponse", function(x){
  yt <- t(x@y)
  yt[x@missing] <- NA
  yt
})

setGeneric("get_n_sites", function(object, ...) standardGeneric("get_n_sites"))

setMethod("get_n_sites", "ubmsResponse", function(object){
  keep_sites <- apply(t(object), 2, function(x) any(!is.na(x)))
  sum(keep_sites)
})

setGeneric("get_n_obs", function(object, ...) standardGeneric("get_n_obs"))

setMethod("get_n_obs", "ubmsResponse", function(object){
  yt <- matrix(t(object), nrow=object@max_obs)
  obs <- apply(yt, 2, function(x) sum(!is.na(x)))
  out <- matrix(obs, ncol=object@max_primary, byrow=TRUE)
  keep_sites <- rowSums(out) > 0
  out[keep_sites,, drop=FALSE]
})

setGeneric("per_sampled", function(object, ...){
             standardGeneric("per_sampled")})

setMethod("per_sampled", "ubmsResponse", function(object){
  yt <- matrix(t(object), nrow=object@max_obs)
  pers <- apply(yt, 2, function(x) any(!is.na(x)))
  pers <- matrix(pers, ncol=object@max_primary, byrow=TRUE)
  keep <- rowSums(pers) > 0
  pers[keep,,drop=FALSE]
})

which_per_sampled <- function(object){
  ps <- per_sampled(object)
  as.vector(unlist(apply(ps, 1, which)))
}

get_n_pers <- function(object){
  ps <- per_sampled(object)
  rowSums(ps)
}

setGeneric("get_subset_inds", function(object, ...) standardGeneric("get_subset_inds"))

setMethod("get_subset_inds", "ubmsResponse", function(object){
  #Indices for y/detection
  nJ <- rowSums(get_n_obs(object))
  #Indices for periods sampled
  nT <- get_n_pers(object)
  #Indices for primary-period level parameters
  #minimum of one phi to handle single-season models
  nP <- rep(max(object@max_primary - 1, 1),
                get_n_sites(object))

  do.call("cbind", lapply(list(nJ, nT, nP), generate_inds))
})

generate_inds <- function(count_vec){
  stopifnot(! 0 %in% count_vec)
  end <- Reduce(`+`, count_vec, accumulate=TRUE)
  cbind(end - count_vec + 1, end)
}

setGeneric("as_vector", function(x, ...) as.vector(x))
setMethod("as_vector", "ubmsResponse", function(x, na.rm=FALSE){
  out <- as.vector(t(x))
  if(na.rm) out <- out[!is.na(out)]
  out
})

setGeneric("get_Kmin", function(object) standardGeneric("get_Kmin"))

setMethod("get_Kmin", "ubmsResponse", function(object){
  yt <- t(object)
  keep <- apply(yt, 2, function(x) !all(is.na(x)))

  yt <- matrix(yt, nrow=object@max_obs)
  out <- apply(yt, 2, function(x){
            if(all(is.na(x))) return(0)
            max(x, na.rm=TRUE)
          })
  out <- matrix(out, ncol=object@max_primary, byrow=TRUE)
  out[keep,,drop=FALSE]
})
