setOldClass("waic")

#' @importClassesFrom rstan stanfit
#' @importFrom loo waic
#' @include submodel.R
setClass("ubmsFit",
  slots=c(call="call",
          psiformula="formula",
          pformula="formula",
          data="unmarkedFrame",
          stanfit="stanfit",
          WAIC="waic",
          submodels="ubmsSubmodelList"
          )
)

setMethod("show", "ubmsFit", function(object){
  
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  show(object@submodels)
  cat(paste0("WAIC: ", round(object@WAIC$estimates[3,1], 3)))
  cat("\n")

})
