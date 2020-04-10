build_stan_inputs <- function(umf, ...){

  y_data <- get_y_data(umf)

  submodels <- list(...)
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

setMethod("get_y_data", "unmarkedFrameOccu", function(object, ...){
  
  y <- getY(object)
  no_detects <- apply(y, 1, function(x) ifelse(sum(x)>0, 0, 1))
  
  list(y=y, M=nrow(y), J=ncol(y), no_detects=no_detects)
})


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
