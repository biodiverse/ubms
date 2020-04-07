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

  pars <- c(paste0("beta_", types), "log_lik")
  for (i in submodels){
    if(get_group_vars(i@formula) > 0){
      pars <- c(pars, paste0(c("b_", "sigma_"), i@type)) 
    }
  }
  pars
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
  has_random <- ifelse(n_group_vars > 0, 1, 0)
  n_random <- get_nrandom(object@formula, object@data)

  out <- list(X=object@X, n_fixed=ncol(object@X), Z=object@Z, 
              n_group_vars=n_group_vars, has_random=has_random,
              n_random=n_random)
  names(out) <- paste0(names(out), "_", object@type)
  out
})

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
