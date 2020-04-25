build_stan_inputs <- function(umf, submodels, ...){
  
  clean <- remove_NA(umf, submodels)
  umf <- clean$umf
  submodels <- clean$submodels

  y_data <- get_y_data(umf, ...)

  submodels <- unname(submodels@submodels)
  types <- sapply(submodels, function(x) x@type)
  submodel_data <- lapply(submodels, get_stan_data)
  submodel_data <- do.call("c", submodel_data)
  
  stan_data <- c(y_data, submodel_data)

  pars <- get_pars(submodels)
  
  list(stan_data=stan_data, pars=pars)
}

get_pars <- function(submodels){

  types <- sapply(submodels, function(x) x@type)
  pars <- paste0("beta_", types)
  for (i in submodels){
    if(has_random(i)){
      pars <- c(pars, b_par(i), sig_par(i))
    }
  }
  c(pars, "log_lik")
}

setGeneric("get_y_data", function(object, ...){
             standardGeneric("get_y_data")})

setMethod("get_y_data", "unmarkedFrameOccu", function(object, K, ...){
  y <- getY(object)
  J <- apply(y, 1, function(x) sum(!is.na(x)))
  no_detects <- apply(y, 1, function(x) ifelse(sum(x, na.rm=TRUE)>0, 0, 1))
  ylong <- as.vector(t(y))
  ylong <- ylong[!is.na(ylong)]

  out <- list(y=ylong, M=nrow(y), J=J)  
  if(missing(K)){
    #Regular occupancy
    return(c(out, no_detects=list(no_detects)))
  } 
  #RN occupancy
  c(out, list(K=ifelse(is.null(K), 20, K), Kmin=1-no_detects))
})

setMethod("get_y_data", "unmarkedFramePCount", 
          function(object, K=NULL, mixture="P", ...){
  y <- getY(object)
  J <- apply(y, 1, function(x) sum(!is.na(x)))
  ylong <- as.vector(t(y))
  ylong <- ylong[!is.na(ylong)]
  Kinfo <- get_K(y, K)
  mixture <- switch(mixture, P={1})
  c(list(y=ylong, M=nrow(y), J=J, mixture=mixture), Kinfo)
})


setGeneric("remove_NA", 
           function(object, submodels, ...) standardGeneric("remove_NA"))

setMethod("remove_NA", c("unmarkedFrame", "ubmsSubmodelList"),
          function(object, submodels, ...){

  y <- getY(object)
  M <- nrow(y)
  J <- ncol(y)

  state <- submodels@submodels$state
  det <- submodels@submodels$det

  ylong <- as.vector(t(y))
  site_idx <- rep(1:M, each=J)
  state_long <- state@X[site_idx,]

  comb <- cbind(y=ylong, state_long, det@X)
  keep_obs <- apply(comb, 1, function(x) !any(is.na(x)))
  
  ylong[!keep_obs] <- NA
  det@X <- det@X[keep_obs,]
  if(has_random(det)) det@Z <- det@Z[keep_obs,]
  
  keep_sites <- unique(site_idx[keep_obs]) 
  state@X <- state@X[keep_sites,]
  if(has_random(state)) state@Z <- state@Z[keep_sites,] 
  ymat <- matrix(ylong, nrow=M, byrow=TRUE)[keep_sites,]
  umf@y <- ymat

  list(umf=umf, submodels=ubmsSubmodelList(state, det))
})


get_K <- function(y, K=NULL){
  ymax <- max(y, na.rm=TRUE)
  if(is.null(K)){
    K <- ymax + 20
  } else {
    if(K < ymax){
      stop("K must be larger than the maximum observed count", call.=FALSE)
    }
  }
  Kmin <- apply(y, 1, max, na.rm=TRUE)
  list(K=K, Kmin=Kmin)
}

setGeneric("get_stan_data", function(object, ...){
             standardGeneric("get_stan_data")})

#' @include submodel.R
setMethod("get_stan_data", "ubmsSubmodel", function(object, ...){
  n_group_vars <- get_group_vars(object@formula)
  has_rand <- has_random(object)
  n_random <- get_nrandom(object@formula, object@data)
  Zinfo <- get_sparse_Z(object@Z)
  out <- list(X=object@X, n_fixed=ncol(object@X), 
              n_group_vars=n_group_vars, has_random=has_rand,
              n_random=n_random)
  out <- c(out, Zinfo)
  names(out) <- paste0(names(out), "_", object@type)
  out
})

get_sparse_Z <- function(Z){
  if(all(dim(Z)==c(0,0))){
    return(list(Zdim=c(0,0,1,1,1), Zw=as.array(0), 
                Zv=as.array(0), Zu=as.array(0)))
  }
  wvu <- rstan::extract_sparse_parts(Matrix::Matrix(Z))
  #Zdim = Z rows, Z cols, length w, length v, length u
  Zdim <- c(dim(Z), sapply(wvu, length))
  list(Zdim=Zdim, Zw=wvu$w, Zv=wvu$v, Zu=wvu$u)
}

get_group_vars <- function(formula){
  rand <- lme4::findbars(formula)
  ifelse(is.null(rand), 0, length(rand))
}

get_nrandom <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    length(unique(data[[col_nm]]))
  })
  as.array(out)
}

split_formula <- function(formula){ 
  p1 <- as.formula(formula[[2]])
  p2 <- as.formula(paste0(formula[[1]], deparse(formula[[3]])))
  list(p1, p2)
}
