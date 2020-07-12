#' @export
stan_distsamp <- function(formula, data, keyfun=c("halfnorm", "exp"), ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)
  keyfun <- match.arg(keyfun)
  det_param <- switch(keyfun, halfnorm={"Scale"}, exp={"Rate"})

  response <- ubmsResponseDistsamp(getY(umf), "P", survey = umf@survey,
                                   dist_breaks=umf@dist.breaks,
                                   tlength=umf@tlength, keyfun=keyfun)
  state <- ubmsSubmodel("Abundance", "state", siteCovs(umf), forms[[2]], "exp")
  det <- ubmsSubmodel(det_param, "det", siteCovs(umf), forms[[1]], "exp")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("distsamp", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitDistsamp", contains = "ubmsFit")


#Response class----------------------------------------------------------------

#' @include response.R
setClass("ubmsResponseDistsamp", contains="ubmsResponse",
         slots = c(survey="character", dist_breaks="numeric", tlength="numeric",
                   keyfun="character"))

ubmsResponseDistsamp <- function(y, z_dist, max_primary = 1, K=NULL,
                                 survey, dist_breaks, tlength, keyfun){
  out <- ubmsResponse(y, "P", z_dist, max_primary, K)
  out <- as(out, "ubmsResponseDistsamp")
  out@survey <- survey
  out@dist_breaks <- dist_breaks
  out@tlength <- tlength
  out@keyfun <- keyfun
  out
}

#' @include inputs.R
setMethod("get_stan_data", "ubmsResponseDistsamp", function(object, ...){
  out <- callNextMethod(object, ...)
  c(out,
    list(point=ifelse(object@survey=="point",1,0),
    dist_breaks=object@dist_breaks,
    keyfun = switch(object@keyfun, halfnorm={0}, exp={1}),
    conv_const = get_conv_const(object)))
})

#Builds the values in the denominator of the detection part of the likelihood
get_conv_const <- function(resp){
  db <- resp@dist_breaks
  w <- diff(db)
  if(resp@survey=="line"){
    a <- t(sapply(resp@tlength, function(x) x*w))
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 1 / t(w*t(u))
  } else {
    a <- pi*db^2
    a <- c(a[1], diff(a))
    M <- ubms:::get_n_sites(resp)
    a <- matrix(rep(a, M), nrow=M, byrow=TRUE)
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 2 * pi / (a * u)
  }
  as.vector(t(out))
}

#Goodness-of-fit---------------------------------------------------------------

#Methods to simulate posterior predictive distributions------------------------
