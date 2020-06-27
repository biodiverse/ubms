 build_stan_inputs <- function(response, submodels, ...){

  y_data <- get_stan_data(response)

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

setGeneric("get_stan_data", function(object, ...){
             standardGeneric("get_stan_data")})


#' @include response.R
setMethod("get_stan_data", "ubmsResponse", function(object, ...){
  list(y = as_vector(object, na.rm=TRUE),
       y_dist = dist_code(object@y_dist),
       z_dist = dist_code(object@z_dist),
       M = get_n_sites(object),
       T = object@max_primary,
       Tsamp = which_per_sampled(object),
       Tsamp_size = length(which_per_sampled(object)),
       J = get_n_obs(object),
       R = sum(get_n_obs(object)),
       si = get_subset_inds(object),
       K = object@K,
       Kmin = get_Kmin(object))
})

dist_code <- function(dist){
  switch(dist, binomial = {0}, Poisson = {1}, P = {1})
}


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
  if(length(formula) != 3) stop("Double right-hand side formula required")
  p1 <- as.formula(formula[[2]])
  p2 <- as.formula(paste0(formula[[1]], deparse(formula[[3]])))
  list(p1, p2)
}
