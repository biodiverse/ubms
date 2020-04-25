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
  b_names <- get_b_names(formula, data)

  sigma_names <- get_sigma_names(formula, data)
  
  new("ubmsSubmodel", name=name, type=type, data=data, formula=formula,
      link=link, X=X, Z=Z, beta_names=beta_names, b_names=b_names,
      sigma_names=sigma_names)
}

get_X <- function(formula, data){
  fixed <- lme4::nobars(formula)
  mf <- model.frame(fixed, data, na.action=stats::na.pass)
  model.matrix(fixed, mf)
}

get_Z <- function(formula, data){
  check_formula(formula, data)
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(matrix(0,0,0))
  Zt <- get_reTrms(formula, data)$Zt
  t(as.matrix(Zt))
}

get_reTrms <- function(formula, data){
  #For compatibility with rhs-only formulas
  #new_data <- cbind(dummy=0, data)
  #new_formula <- as.formula(paste0("dummy", 
  #                          paste(as.character(formula), collapse="")))  
  fb <- lme4::findbars(formula)
  mf <- model.frame(lme4::subbars(formula), data, na.action=stats::na.pass)
  lme4::mkReTrms(fb, mf)
}

check_formula <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(invisible())
 
  char <- paste(deparse(formula))
  if(grepl(":|/", char)){
    stop("Nested random effects (using / and :) are not supported",
         call.=FALSE)
  }
  theta <- get_reTrms(formula, data)$theta
  if(0 %in% theta){
    stop("Correlated slopes and intercepts are not supported. Use || instead of |.",
         call.=FALSE)
  }
}

get_b_names <- function(formula, data){
  if(is.null(lme4::findbars(formula))) return(NA_character_)
  group <- get_reTrms(formula, data)
  group_nms <- names(group$cnms)
  z_nms <- character()
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      z_nms <- c(z_nms, paste0(nms_i, ":", levels(group$flist[[nm]])))
    } else {
      z_nms <- c(z_nms, c(t(sapply(paste0(nms_i), paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  z_nms  
}

get_sigma_names <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(NA_character_)  
  nms <- get_reTrms(formula, data)$cnms
  nms <- paste0(unlist(nms), "|", names(nms)) 
  nms <- gsub("(Intercept)", "1", nms, fixed=TRUE)
  paste0("sigma [", nms, "]")
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
