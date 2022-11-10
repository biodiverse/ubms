 build_stan_inputs <- function(name, response, submodels, log_lik, ...){

  model_code <- name_to_modelcode(name)
  y_data <- get_stan_data(response)

  pars <- get_pars(submodels, name, log_lik)
  submodels <- unname(submodels@submodels)
  types <- sapply(submodels, function(x) x@type)
  submodel_data <- lapply(submodels, get_stan_data)
  submodel_data <- do.call("c", submodel_data)

  # Change this later?
  submodel_data <- add_placeholder_priors(submodel_data, types)

  stan_data <- c(model_code, y_data, submodel_data)

  list(stan_data=stan_data, pars=pars)
}

name_to_modelcode <- function(name){
  list(model_code=switch(name, occu={0}, occuRN={1}, pcount={2}, colext={3},
                         distsamp={4}, multinomPois={5}, occuTTD={6}))
}

# Add prior info for submodel not being used in a given model
# Placeholder info is still needed to satisfy the Stan data block
add_placeholder_priors <- function(submodel_data, types){
  # this is hacky
  if(! "shape" %in% types){
    submodel_data$prior_dist_shape <- c(0,0,0)
    submodel_data$prior_pars_shape <- matrix(rep(0,6), nrow=3)
  }
  if(! "scale" %in% types){
    submodel_data$prior_dist_scale <- c(0,0,0)
    submodel_data$prior_pars_scale <- matrix(rep(0,6), nrow=3)
  }
  submodel_data
}

setGeneric("get_pars", function(object, ...) standardGeneric("get_pars"))

setMethod("get_pars", "ubmsSubmodelList", function(object, model_name, log_lik, ...){
  #Remove placeholder submodels - we don't want to save those parameters
  submodels <- object@submodels
  submodels <- submodels[!sapply(submodels, is_placeholder)]
  submodels <- unname(submodels)
  out <- unlist(lapply(submodels, get_pars))

  if(model_name == "distsamp") log_lik <- TRUE
  if(any(sapply(submodels, has_spatial))) log_lik <- TRUE

  if(log_lik) out <- c(out, "log_lik")
  out
})

setMethod("get_pars", "ubmsSubmodel", function(object, ...){
  out <- paste0("beta_", object@type)
  if(has_random(object)){
    out <- c(out, b_par(object), sig_par(object))
  }
  out
})


setGeneric("get_stan_data", function(object, ...){
             standardGeneric("get_stan_data")})


#' @include response.R
setMethod("get_stan_data", "ubmsResponse", function(object, ...){
  out <- list(y = as_vector(object, na.rm=TRUE),
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
  c(out, get_auxiliary_data(object))
})

dist_code <- function(dist){
  switch(dist, binomial = {0}, P = {1}, NB = {2}, double = {0}, removal = {1},
         halfnorm={0}, exp={1}, hazard={2}, weibull={3}
  )
}

setGeneric("get_auxiliary_data", function(object, ...){
             standardGeneric("get_auxiliary_data")})

setMethod("get_auxiliary_data", "ubmsResponse", function(object, ...){
  list(aux1=numeric(0), aux2=numeric(0), aux3=numeric(0),
       n_aux1=0, n_aux2=0, n_aux3=0)
})

setMethod("get_stan_data", "ubmsSubmodelScalar", function(object, ...){
  out <- process_priors(object)
  names(out) <- paste0(names(out), "_", object@type)
  out
})

#' @include submodel.R
setMethod("get_stan_data", "ubmsSubmodel", function(object, ...){
  n_group_vars <- get_group_vars(object@formula)
  has_rand <- has_random(object)
  n_random <- get_nrandom(object@formula, object@data)
  Zinfo <- get_sparse_Z(Z_matrix(object, na.rm=TRUE))
  X <- model.matrix(object, na.rm=TRUE, warn=TRUE)
  priors <- process_priors(object)
  off <- model_offset(object, na.rm=TRUE)
  out <- list(X=X, offset=off, n_obs=nrow(X), n_fixed=ncol(X),
              n_group_vars=n_group_vars, has_random=has_rand, n_random=n_random)
  out <- c(out, Zinfo, priors)
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
  check_formula(formula, data)
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    length(levels(as.factor(data[[col_nm]])))
  })
  as.array(out)
}

split_formula <- function(formula){
  if(length(formula) != 3) stop("Double right-hand side formula required")
  char <- lapply(formula, function(x){
            paste(deparse(x), collapse="")
          })
  p1 <- as.formula(char[[2]])
  p2 <- as.formula(paste("~", char[[3]]))
  list(det=p1, state=p2)
}
