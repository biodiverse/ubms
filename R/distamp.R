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
#' @references Royle, J. A., Dawson, D. K., & Bates, S. (2004). Modeling
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

#Get detection probability-----------------------------------------------------

#' @importFrom unmarked getP
setMethod("getP", "ubmsFitDistsamp", function(object, draws=NULL, ...){
  samples <- get_samples(object, draws)
  resp <- object@response
  praw <- t(sim_p(object, samples))
  praw <- array(praw, c(resp@max_obs, get_n_sites(resp), length(samples)))
  aperm(praw, c(2,1,3))
})

setGeneric("sim_p", function(object, samples, ...) standardGeneric("sim_p"))

setMethod("sim_p", "ubmsFitDistsamp", function(object, samples, ...){
  resp <- object@response
  resp@output <- "abund" #Don't adjust for area
  db <- resp@dist_breaks
  conv_const <- get_conv_const(resp)
  inds <- get_subset_inds(resp)[,1:2]
  param1 <- sim_lp(object, "det", transform=TRUE, newdata=NULL,
                          samples=samples, re.form=NULL)

  stopifnot(ncol(param1) == get_n_sites(resp))
  param2 <- NULL
  if(resp@keyfun=="halfnorm"){
    pfun <- ifelse(resp@survey=="line", distprob_normal_line, distprob_normal_point)
  }
  out <- sapply(1:length(samples), function(i){
    pfun(param1[i,], param2[i,], db, conv_const, inds)
  })
  t(out)
})

distprob_normal_line <- function(sigma, param2, db, conv_const, inds){
  out <- sapply(1:length(sigma), function(i){
    int <-  pnorm(db[-1], 0, sd=sigma) - pnorm(db[-length(db)], 0, sd=sigma[i])
    int <- int / dnorm(0, 0, sd=sigma[i])
    if(!is.null(conv_const)) int <- int * conv_const[inds[i,1]:inds[i,2]]
    int
  })
  as.vector(out)
}

distprob_normal_point <- function(sigma, param2, db, conv_const, inds){
  out <- sapply(1:length(sigma), function(i){
    s2 <- sigma[i]^2
    a <- db[-1]; b <- db[-length(db)]
    int <- s2 * ((1 - exp(-a*a / (2*s2))) - (1-exp(-b*b / (2*s2))))
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

#Histogram---------------------------------------------------------------------

#' @importFrom graphics hist
#' @importFrom ggplot2 geom_histogram
setMethod("hist", "ubmsFitDistsamp", function(x, draws=30, ...){
  samples <- get_samples(x, draws)
  hist_data <- get_hist_data(x)
  mean_line <- get_mean_line(x)

  out <- ggplot(hist_data, aes_string(x="x")) +
    geom_histogram(aes_string(y="..density.."),fill='transparent',
                   col='black',breaks=x@response@dist_breaks)
  if(draws > 0){
    sample_lines <- get_sample_lines(x, samples)
    out <- out +
      geom_line(data=sample_lines, aes_string(x="x", y="val", group="ind"),
              alpha=0.3)
  }
  out +
    geom_line(data=mean_line, aes_string(x="x", y="val"), col='red') +
    labs(x=paste0("Distance (", x@response@units_in,")"), y="Density") +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text=element_text(size=12), axis.title=element_text(size=14))
})

get_hist_data <- function(fit){
  resp <- fit@response
  db <- fit@response@dist_breaks
  bin <- db[-length(db)] + diff(db)/2
  counts <- colSums(resp@y, na.rm=TRUE)
  data.frame(x=rep(bin, times=counts))
}

get_mean_line <- function(fit){
  db <- fit@response@dist_breaks
  xseq <- seq(db[1], db[length(db)], length.out=1000)
  par1 <- summary(fit, "det")
  if(nrow(par1) > 1) warning("Ignoring covariate effects", call.=FALSE)
  par1 <- exp(par1[1,1])
  detfun <- get_detfun(fit)
  data.frame(x=xseq, val=detfun(xseq, par1))
}

get_sample_lines <- function(fit, samples){
  db <- fit@response@dist_breaks
  xseq <- seq(db[1], db[length(db)], length.out=1000)
  par1 <- exp(extract(fit, "beta_det")[[1]][,1])
  detfun <- get_detfun(fit)
  sample_lines <- lapply(samples, function(i){
                         data.frame(x=xseq, val=detfun(xseq, par1[i]),ind=i)
                      })
  do.call("rbind", sample_lines)
}

get_detfun <- function(fit){
  resp <- fit@response
  if(resp@keyfun == "halfnorm"){
    out <- ifelse(resp@survey=="line", unmarked::dxhn, unmarked::drhn)
  } else if(resp@keyfun == "exp"){
    out <- ifelse(resp@survey=="line", unmarked::dxexp, unmarked::drexp)
  }
  out
}
