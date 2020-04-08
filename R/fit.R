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

#' @importFrom rstan extract
#' @export
setMethod("extract", "ubmsFit", 
  function(object, pars, permuted=TRUE, inc_warmup=FALSE,
          include=TRUE){
  rstan::extract(object@stanfit, pars, permuted, inc_warmup, include)
})

setMethod("[", c("ubmsFit", "character", "missing", "missing"),
  function(x, i){

  types <- names(x@submodels@submodels)
  if(! i %in% types){
    stop(paste("Possible types are:", paste(types, collapse=", ")),
         call. = FALSE)
  }
  x@submodels@submodels[[i]]
})

check_stanfit <- function(object){
  if(object@mode == 2L || object@mode == 1L){
    stop("Fitting model failed", call.=FALSE)
  }
}
