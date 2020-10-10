#' Fit Time-to-detection Occupancy Models
#'
#' Fit time-to-detection occupancy models of Garrard et al.
#' (2008, 2013). Time-to-detection can be modeled with either an exponential
#' or Weibull distribution.
#'
#' @param psiformula Right-hand sided formula for the initial probability of
#'  occupancy at each site.
#' @param gammaformula Right-hand sided formula for colonization probability.
#'  Currently ignored as dynamic models are not yet supported.
#' @param epsilonformula Right-hand sided formula for extinction probability.
#'  Currently ignored as dynamic models are not yet supported.
#' @param detformula Right-hand sided formula for mean time-to-detection.
#' @param data \code{unmarkedFrameOccuTTD} object that supplies the data
#'  (see \code{\link{unmarkedFrameOccuTTD}}).
#' @param ttdDist Distribution to use for time-to-detection; either
#'  \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
#'  which adds an additional shape parameter \eqn{k}.
#' @param linkPsi Link function for the occupancy model. Only option is
#'  \code{"logit"} for now, in the future \code{"cloglog"}
#'  will be supported for the complimentary log-log link.
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccuTTD} object describing the model fit.
#'
#' @examples
#' \dontrun{
#' }
#'
#' @references Garrard, G.E., Bekessy, S.A., McCarthy, M.A. and Wintle, B.A. 2008.
#' When have we looked hard enough? A novel method for setting minimum survey effort
#' protocols for flora surveys. Austral Ecology 33: 986-998.
#'
#' Garrard, G.E., McCarthy, M.A., Williams, N.S., Bekessy, S.A. and Wintle,
#' B.A. 2013. A general model of detectability using species traits. Methods in
#' Ecology and Evolution 4: 45-52.
#'
#' Kery, Marc, and J. Andrew Royle. 2016. Applied Hierarchical Modeling in
#' Ecology, Volume 1. Academic Press.
#'
#' @seealso \code{\link{occuTTD}}, \code{\link{unmarkedFrameOccuTTD}}
#' @include fit.R
#' @export
stan_occuTTD <- function(psiformula=~1, gammaformula=~1, epsilonformula=~1,
                         detformula=~1, data, ttdDist=c("exp", "weibull"),
                         linkPsi=c("logit"), ...){

  if(data@numPrimary > 1) stop("Dynamic models not yet supported", call.=FALSE)
  umf <- process_umf(data)
  ttdDist <- match.arg(ttdDist)
  linkPsi <- match.arg(linkPsi)

  response <- ubmsResponseOccuTTD(umf, ttdDist)
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), detformula, "exp")

  shape <- placeholderSubmodel("shape")
  if(ttdDist=="weibull"){
    shape <- ubmsSubmodelScalar("Shape", "shape", "exp")
  }

  submodels <- ubmsSubmodelList(state, det, shape)

  ubmsFit("occuTTD", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include occu.R
setClass("ubmsFitOccuTTD", contains = "ubmsFitOccu")


#Response class----------------------------------------------------------------

#' @include response.R
setClass("ubmsResponseOccuTTD", contains="ubmsResponse",
         slots=c(surveyLength="matrix"))

ubmsResponseOccuTTD <- function(umf, ttdDist){
  out <- ubmsResponse(getY(umf), ttdDist, "binomial", umf@numPrimary,
                      ceiling(max(getY(umf), na.rm=TRUE)))
  out <- as(out, "ubmsResponseOccuTTD")
  out@surveyLength <- umf@surveyLength
  out
}

#Bundle delta (indicator if observation occurred before max time) for Stan
#Also include y as aux2 (can't use normal y slot because y is not an integer)

#' @include inputs.R
setMethod("get_auxiliary_data", "ubmsResponseOccuTTD", function(object, ...){
  yt <- t(object)
  delta <- ifelse(yt < t(object@surveyLength), 1, 0)
  keep <- apply(delta, 2, function(x) !all(is.na(x)))
  delta <- delta[,keep,drop=FALSE]
  delta <- as.vector(delta)
  delta <- delta[!is.na(delta)]
  y <- as_vector(object, na.rm=TRUE) #Need to pass to aux because it's not integer
  stopifnot(length(delta) == length(y))
  list(aux1=delta,
       aux2=y,
       aux3=numeric(0),
       n_aux1=length(delta), n_aux2=length(y), n_aux3=0)
})

setMethod("get_stan_data", "ubmsResponseOccuTTD", function(object, ...){
  out <- methods::callNextMethod(object, ...)
  out$y <- out$aux1
  out$Kmin[] <- 0
  out
})


#Get detection probability-----------------------------------------------------

#' @importFrom unmarked getP
setMethod("getP", "ubmsFitOccuTTD", function(object, draws=NULL, ...){
  samples <- get_samples(object, draws)
  resp <- object@response
  praw <- t(sim_p(object, samples))
  praw <- array(praw, c(resp@max_obs, nrow(resp@y), length(samples)))
  aperm(praw, c(2,1,3))
})

setMethod("sim_p", "ubmsFitOccuTTD", function(object, samples, ...){
  resp <- object@response
  tdist <- resp@y_dist
  tmax <- as.vector(t(resp@surveyLength))
  lam <- sim_lp(object, "det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=NULL)
  lam[,object["det"]@missing] <- NA

  if(tdist == "exp"){
    out <- sapply(1:length(samples), function(i){
      stats::pexp(tmax, lam[i,])
    })
  } else if(tdist == "weibull"){
    shape <- t(sim_lp(object, "shape", transform=TRUE, newdata=NULL,
                    samples=samples, re.form=NULL))
    out <- sapply(1:length(samples), function(i){
      stats::pweibull(tmax, shape[i], 1/lam[i,])
    })
  }

  t(out)
})


#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R


#Goodness-of-fit---------------------------------------------------------------

#' @include gof.R

