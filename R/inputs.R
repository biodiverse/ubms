build_stan_inputs <- function(umf, submodels, ...){

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

setMethod("get_y_data", "unmarkedFrame", function(object, ...){
  y <- getY(object)
  list(y=y, M=nrow(y), J=ncol(y))
})

setMethod("get_y_data", "unmarkedFrameOccu", function(object, K, ...){
  out <- callNextMethod(object, ...)
  
  no_detects <- apply(out$y, 1, function(x) ifelse(sum(x)>0, 0, 1)) 
  if(missing(K)){
    #Regular occupancy
    return(c(out, no_detects=list(no_detects)))
  } 
  #RN occupancy
  c(out, list(K=ifelse(is.null(K), 20, K), Kmin=1-no_detects))
})

setMethod("get_y_data", "unmarkedFramePCount", 
          function(object, K=NULL, mixture="P", ...){
  out <- callNextMethod(object, ...)
  Kinfo <- get_K(out$y, K)
  mixture <- switch(mixture, P={1})
  c(out, Kinfo, mixture=list(mixture))
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
