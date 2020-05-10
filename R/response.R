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
  out <- new("ubmsResponse", y = y, y_dist= y_dist, z_dist = z_dist,
          max_primary = max_primary,
          max_obs = get_max_obs(y, max_primary),
          K = get_K(y, K)
      )
  out@missing <- is.na(as_vector(out))
  out
}

get_max_obs <- function(y, max_primary){
  if(ncol(y) %% max_primary != 0){
    stop("Primary periods must have same number of secondary periods", 
         call.=FALSE) 
  }
  ncol(y) / max_primary
}

get_K <- function(y, K=NULL){
  ymax <- max(y, na.rm=TRUE)
  if(is.null(K)) return(ymax + 20)
  if(K < ymax) stop("K must be larger than max y value", call.=FALSE)
  K
}

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
  out <- matrix(obs, nrow=object@max_primary)
  keep_sites <- colSums(out) > 0
  out[, keep_sites, drop=FALSE]
})

setGeneric("as_vector", function(x, ...) as.vector(x))
setMethod("as_vector", "ubmsResponse", function(x, na.rm=FALSE){  
  out <- as.vector(t(x))
  if(na.rm) out <- out[!is.na(out)]
  out
})

setGeneric("get_Kmin", function(object) standardGeneric("get_Kmin"))

setMethod("get_Kmin", "ubmsResponse", function(object){
  yt <- t(object)
  out <- apply(yt, 2, function(x){
           if(all(is.na(x))) return(NA)
           max(x, na.rm=TRUE)
          })
  out[!is.na(out)]
})
