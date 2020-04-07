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
    sigma_names = "character"
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
    sigma_names = "NA_character_"
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
