#' Fit the Royle et al. (2004) Distance Sampling Model
#'
#' This function fits the hierarchical distance sampling model of Royle
#' et al. (2004) to line or point transect data recorded in discerete
#' distance intervals.
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link{unmarkedFrameDS}} object
#' @param keyfun One of the following detection functions:
#'  \code{"halfnorm"} for half-normal or \code{"exp"} for negative exponential
#' @param output Model either density \code{"density"} or abundance \code{"abund"}
#' @param unitsOut Units of density. Either \code{"ha"} or \code{"kmsq"} for
#'  hectares and square kilometers, respectively.
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitDistsamp} object describing the model fit.
#'
#' @references Royle, J. A., D. K. Dawson, & Bates, S. (2004). Modeling
#'  abundance effects in distance sampling. Ecology 85: 1591-1597.
#'
#' @seealso \code{\link{distsamp}}, \code{\link{unmarkedFrameDS}}
#' @export
stan_distsamp <- function(formula, data, keyfun=c("halfnorm", "exp"),
                          output=c("density", "abund"), unitsOut=c("ha", "kmsq"),
                          ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)
  keyfun <- match.arg(keyfun)
  unitsOut <- match.arg(unitsOut)
  output <- match.arg(output)
  state_param <- switch(output, density={"Density"}, abund={"Abundance"})
  det_param <- switch(keyfun, halfnorm={"Scale"}, exp={"Rate"})

  response <- ubmsResponseDistsamp(data, "P", keyfun, output, unitsOut)
  state <- ubmsSubmodel(state_param, "state", siteCovs(umf), forms[[2]], "exp")
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
                   keyfun="character", output="character",
                   units_in="character", units_out="character"))

ubmsResponseDistsamp <- function(umf, z_dist, keyfun, output, units_out, K=NULL){
  max_primary <- ifelse(methods::.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  out <- ubmsResponse(umf@y, "P", z_dist, max_primary, K)
  out <- as(out, "ubmsResponseDistsamp")
  out@survey <- umf@survey; out@keyfun <- keyfun; out@output <- output
  out@dist_breaks <- umf@dist.breaks; out@tlength <- umf@tlength
  out@units_in <- umf@unitsIn; out@units_out <- units_out
  out
}


#' @include inputs.R
setMethod("get_stan_data", "ubmsResponseDistsamp", function(object, ...){
  out <- callNextMethod(object, ...)
  c(out,
    list(point=ifelse(object@survey=="point",1,0),
    db=object@dist_breaks,
    keyfun = switch(object@keyfun, halfnorm={0}, exp={1}),
    conv_const = get_conv_const(object)))
})

#Builds the values in the denominator of the detection part of the likelihood
get_conv_const <- function(resp){
  db <- resp@dist_breaks
  M <- get_n_sites(resp)

  if(resp@survey=="line"){
    w <- diff(db)
    a <- get_area(resp)
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 1 / matrix(rep(w, M), ncol=length(w), byrow=TRUE) * u
  } else {
    a <- get_area(resp)
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 2 * pi / a * u
  }
  out <- as.vector(t(out))
  area_adjust <- get_area_adjust(resp)
  out * rep(area_adjust, each=length(db)-1)
}

get_area <- function(resp){
  db <- resp@dist_breaks
  if(resp@survey=="line"){
    w <- diff(db)
    a <- t(sapply(resp@tlength, function(x) x*w))
  } else {
    a <- diff(pi*db^2)
    M <- get_n_sites(resp)
    a <- matrix(rep(a, M), nrow=M, byrow=TRUE)
  }
  a
}

get_area_adjust <- function(resp){
  if(resp@output == "abund") return(rep(1, get_n_sites(resp)))
  area <- get_area(resp)
  switch(resp@survey,
    line = A <- rowSums(area) * 2,
    point = A <- rowSums(area))
  switch(resp@units_in,
    m = A <- A / 1e6,
    km = A <- A)
  switch(resp@units_out,
    ha = A <- A * 100,
    kmsq = A <- A)
  A
}

#Goodness-of-fit---------------------------------------------------------------

#Methods to simulate posterior predictive distributions------------------------
