get_X <- function(formula, data, target_rows){
  if(is.null(data)){
    data <- data.frame(dummy=rep(0,target_rows))
  }
  fixed <- lme4::nobars(formula)
  model.matrix(fixed, data)
}

get_fixed_names <- function(dm_list){

  lapply(dm_list, colnames)

}

get_rand_names <- function(...){

  forms <- list(...)
  lapply(forms, function(x){
    rand <- lme4::findbars(x[[1]])
    if(length(rand)==0) return(NULL)

    sapply(rand, function(y){
      col_nm <- as.character(y[[3]])
      dat <- x[[2]][col_nm]
      paste0("sigma [",col_nm,"]")
    })

  })

}

convert_form <- function(x){
  as.formula(paste0("dummy", paste(as.character(x),collapse="")))
}

get_Z <- function(formula, data){

  rand <- lme4::findbars(formula)
  if(is.null(rand)) return( matrix(0,0,0))

  #For compatibility with rhs-only formulas
  new_data <- cbind(dummy=0, data)
  Z <- lme4::glFormula(convert_form(formula), new_data)$reTrms$Zt
  t(as.matrix(Z))

}

get_grpvars <- function(formula){
  rand <- lme4::findbars(formula)
  ifelse(is.null(rand), 0, length(rand))
}

get_levels <- function(formula, data){

  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    dat <- data[[col_nm]]
    length(unique(dat))
  })

  as.array(out)
}

build_inputs <- function(psiformula, pformula, umf){

  y <- getY(umf)
  M <- nrow(y)
  J <- ncol(y)
  no_detects <- apply(y, 1, function(x) ifelse(sum(x)>0, 0, 1))

  X_occ <- get_X(psiformula, siteCovs(umf), M)
  X_det <- get_X(pformula, obsCovs(umf), M*J)
  nFP_occ <- ncol(X_occ)
  nFP_det <- ncol(X_det)

  Z_occ <- get_Z(psiformula, siteCovs(umf))
  Z_det <- get_Z(pformula, obsCov(umf))

  n_grpvars_occ <- get_grpvars(psiformula)
  occ_has_random <- ifelse(n_grpvars_occ>0, 1, 0)

  n_grpvars_det <- get_grpvars(pformula)
  det_has_random <- ifelse(n_grpvars_det>0, 1, 0)

  nRE_occ <- get_levels(psiformula, siteCovs(umf))
  nRE_det <- get_levels(pformula, obsCovs(umf))

  stan_data <- mget(c("y", "M", "J", "no_detects", "X_occ", "X_det", 
                      "nFP_occ", "nFP_det", "Z_occ", "Z_det",
                      "n_grpvars_occ", "n_grpvars_det", "occ_has_random", 
                      "det_has_random", "nRE_occ", "nRE_det"))

  params <- c("beta_occ", "beta_det", "log_lik")

  if(occ_has_random){
    params <- c(params, "sigma_occ", "b_occ")
  }
  if(det_has_random){
    params <- c(params, "sigma_det", "b_det")
  }

  fixed_names <- get_fixed_names(list(occ=X_occ, det=X_det))
  random_names <- get_rand_names(occ=list(psiformula, siteCovs(umf)), 
                                 det=list(pformula, obsCovs(umf)))

  list(params=params, stan_data=stan_data, fixed_names=fixed_names, 
       random_names=random_names)

}
