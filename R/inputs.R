build_stan_inputs <- function(submodels, umf, ...){
 
  y_data <- get_y_data(umf, submodels, ...)

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


setGeneric("get_y_data", function(object, ...) standardGeneric("get_y_data"))

setMethod("get_y_data", "unmarkedFrameOccu", function(object, submod, K=NULL, ...){
  yt <- get_yt(object, submod)
  J <- apply(yt, 2, function(x) sum(!is.na(x)))
  ylong <- stats::na.omit(as.vector(yt))
  out <- list(y=ylong, M=ncol(yt), J=J)
  c(out, get_K(yt, K))
})

setMethod("get_y_data", "unmarkedFramePCount", 
          function(object, submod, K=NULL, mixture="P", ...){
  yt <- get_yt(object, submod)
  J <- apply(yt, 2, function(x) sum(!is.na(x)))
  Kinfo <- get_K(yt, K)
  ylong <- stats::na.omit(as.vector(yt))
  mixture <- switch(mixture, P={1})
  c(list(y=ylong, M=ncol(yt), J=J, mixture=mixture), Kinfo)
})


setGeneric("get_yt", function(object, ...) standardGeneric("get_yt"))

setMethod("get_yt", "unmarkedFrame", function(object, submod, ...){
  yt <- t(getY(object))
  yt[submod["det"]@missing] <- NA
  if(has_missing(submod["state"])){
    yt <- yt[, -submod["state"]@missing, drop=FALSE]
  }
  yt
})


get_K <- function(yt, K=NULL){
  ymax <- max(yt, na.rm=TRUE)
  if(is.null(K)){
    K <- ymax + 20
  } else {
    if(K < ymax){
      stop("K must be larger than the maximum observed count", call.=FALSE)
    }
  }
  Kmin <- apply(yt, 2, max, na.rm=TRUE)
  list(K=K, Kmin=Kmin)
}

setGeneric("get_stan_data", function(object, ...){
             standardGeneric("get_stan_data")})

#' @include submodel.R
setMethod("get_stan_data", "ubmsSubmodel", function(object, ...){
  n_group_vars <- get_group_vars(object@formula)
  has_rand <- has_random(object)
  n_random <- get_nrandom(object@formula, object@data)
  Zinfo <- get_sparse_Z(Z_matrix(object, na.rm=TRUE))
  X <- model.matrix(object, na.rm=TRUE)
  out <- list(X=X, n_fixed=ncol(X), n_group_vars=n_group_vars, 
              has_random=has_rand, n_random=n_random)
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
