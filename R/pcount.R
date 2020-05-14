#' @include fit.R 
setClass("ubmsFitPcount", contains = "ubmsFit")

#' @export
stan_pcount <- function(formula, data, K=NULL, mixture="P", ...){
  
  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), y_dist="binomial", z_dist=mixture, K=K)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")

  ubmsFit("pcount", match.call(), data, response, submodels, ...)
}

