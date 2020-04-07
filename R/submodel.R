setClass("ubmsSubmodel",
  slots = c(
    name = "character",
    type = "character",
    data = "data.frame",
    formula = "formula",
    link = "character",
    X = "matrix",
    Z = "matrix",
    beta_names = "character",
    b_names = "character",
    sigma_names = "character",
    fixed_estimates = "data.frame",
    random_estimates = "data.frame"
  ),
  prototype = list(
    name = NA_character_,
    type = NA_character_,
    data = data.frame(),
    formula = ~1,
    link = NA_character_,
    X = matrix(1),
    Z = matrix(0,0,0),
    beta_names = NA_character_,
    b_names = "NA_character_",
    sigma_names = "NA_character_",
    fixed_estimates = data.frame(),
    random_estimates = data.frame()
  )
)

ubmsSubmodel <- function(name, type, data, formula, link){

  X <- get_X(formula, data)
  beta_names <- colnames(X)
  
  Z <- get_Z(formula, data)
  b_names <- get_Z_names(Z)

  sigma_names <- get_sigma_names(formula, data)
  
  new("ubmsSubmodel", name=name, type=type, data=data, formula=formula,
      link=link, X=X, Z=Z, beta_names=beta_names, b_names=b_names,
      sigma_names=sigma_names)
}

get_X <- function(formula, data){
  fixed <- lme4::nobars(formula)
  model.matrix(fixed, data)
}

get_Z <- function(formula, data){

  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(matrix(0,0,0))

  #For compatibility with rhs-only formulas
  new_data <- cbind(dummy=0, data)
  new_formula <- as.formula(paste0("dummy", 
                            paste(as.character(formula), collapse="")))

  Z <- lme4::glFormula(new_formula, new_data)$reTrms$Zt
  t(as.matrix(Z))
}

get_Z_names <- function(Z){
  out <- colnames(Z)
  if(is.null(out)) out <- NA_character_
  out
}

get_sigma_names <- function(formula, data){

  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(NA_character_)

  sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    dat <- data[col_nm] #check it is in data
    paste0("sigma [", col_nm, "]")
  })

}

setMethod("model.matrix", "ubmsSubmodel", function(object, newdata, ...){

  if(missing(newdata)) return(object@X)

  data <- object@data
  formula <- lme4::nobars(object@formula)
  fac_col <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_col, levels)
  mf <- model.frame(formula, data)
  xlevs <- xlevs[names(xlevs) %in% names(mf)]
  model.matrix(formula, model.frame(stats::terms(mf), newdata, 
                                    na.action=stats::na.pass, xlev=xlevs))

})

setGeneric("add_estimates", function(object, stanfit, ...){
  standardGeneric("add_estimates")
})

setMethod("add_estimates", "ubmsSubmodel", function(object, stanfit, ...){
  fixed <- rstan::summary(stanfit, beta_par(object))
  fixed <- as.data.frame(fixed$summary)
  rownames(fixed) <- object@beta_names
  object@fixed_estimates <- fixed

  if(!is.na(object@sigma_names)){
    random <- rstan::summary(stanfit, sig_par(object))
    random <- as.data.frame(random$summary)
    rownames(random) <- object@sigma_names
    object@random_estimates <- random
  }
  return(object)
})

setMethod("show", "ubmsSubmodel", function(object){
  
  if(nrow(object@fixed_estimates) == 0){
    cat(paste0(object@name,":\n"))
    cat("No estimates yet\n")
    return(invisible)
  }
  out_df <- object@fixed_estimates[,c(1,3,4,8:10)]

  if(nrow(object@random_estimates)>0){
    out_df <- rbind(out_df, object@random_estimates[,c(1,3,4,8:10)])
  }
  names(out_df)[1:2] <- c("Estimate", "SD")

  cat(paste0(object@name,":\n"))
  print(out_df, digits=3)
})

setGeneric("has_random", function(object){
  standardGeneric("has_random")
})

setMethod("has_random", "ubmsSubmodel", function(object){
  !is.null(lme4::findbars(object@formula))
})

#Quickly generate parameter names from ubmsSubmodel
b_par <- function(object){
  paste0("b_", object@type)
}

beta_par <- function(object){
  paste0("beta_", object@type)
}

sig_par <- function(object){
  paste0("sigma_", object@type)
}

setClass("ubmsSubmodelList", slots=c(submodels="list"),
         prototype=list(submodels=list()))

ubmsSubmodelList <- function(...){
  submodels <- list(...)
  names(submodels) <- sapply(submodels, function(x) x@type)
  new("ubmsSubmodelList", submodels=submodels)
}

setMethod("add_estimates", "ubmsSubmodelList", function(object, stanfit, ...){
  for (i in 1:length(object@submodels)){
    object@submodels[[i]] <- add_estimates(object@submodels[[i]], stanfit)
  }
  object
})

setMethod("show", "ubmsSubmodelList", function(object){

  for (est in object@submodels){
    show(est)
    cat("\n")
  }

})
