#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of
#' MacKenzie et al. (2002).
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link{unmarkedFrameOccu}} object
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitOccu} object describing the model fit.
#'
#' @examples
#' \donttest{
#' data(frogs)
#' pferUMF <- unmarkedFrameOccu(pfer.bin)
#'
#' #Add some covariates
#' siteCovs(pferUMF) <- data.frame(cov1=rnorm(numSites(pferUMF)))
#'
#' #Fit model
#' (fm <- stan_occu(~1~cov1, pferUMF, chains=3, iter=300))
#' }
#'
#' @references MacKenzie DI, Nichols JD, Lachman GB, Droege S, Royle JA,
#'  Langtimm CA. 2002. Estimating site occupancy rates when detection
#'  probabilities are less than one. Ecology 83: 2248-2255.
#'
#' @seealso \code{\link{occu}}, \code{\link{unmarkedFrameOccu}}
#' @include fit.R
#' @export
stan_occu <- function(formula, data, ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), "binomial", "binomial")
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), forms[[2]], "plogis")
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), forms[[1]], "plogis")
  submodels <- ubmsSubmodelList(state, det)

  ubmsFit("occu", match.call(), data, response, submodels, ...)
}


#Goodness-of-fit---------------------------------------------------------------

#' @describeIn gof
#' Applies the MacKenzie-Bailey chi-square goodness of fit test for
#' ocupancy models (MacKenzie and Bailey 2004).
#' @references MacKenzie, D. I., & Bailey, L. L. (2004). Assessing the
#'  fit of site-occupancy models. Journal of Agricultural, Biological,
#'  and Environmental Statistics, 9(3), 300-318.
#' @include gof.R
setMethod("gof", "ubmsFitOccu", function(object, draws=NULL, quiet=FALSE, ...){
  sim_gof(object, draws, mb_chisq, "MacKenzie-Bailey Chi-square", quiet)
})


#Methods to simulate posterior predictive distributions------------------------

setGeneric("knownZ", function(object, ...) standardGeneric("knownZ"))

setMethod("knownZ", "ubmsFitOccu", function(object, ...){
  apply(object@data@y, 1, function(x) sum(x, na.rm=T)>0)
})

#' @include posterior_predict.R
setMethod("sim_z", "ubmsFitOccu", function(object, samples, re.form, ...){

  p_post <- t(sim_p(object, samples))
  psi_post <- t(sim_lp(object, submodel="state", transform=TRUE, newdata=NULL,
                       samples=samples, re.form=re.form))
  psi_post[object["state"]@missing] <- NA

  M <- nrow(psi_post)
  J <- nrow(p_post) / M
  nsamp <- length(samples)

  p_post <- array(p_post, c(J,M,nsamp))
  p_post <- aperm(p_post, c(2,1,3))

  z_post <- matrix(NA, M, nsamp)

  known_z <- knownZ(object)
  z_post[known_z,] <- 1
  unkZ <- which(!known_z)

  q_post <- 1 - p_post

  for (i in 1:nsamp){
    psi <- psi_post[unkZ, i, drop=FALSE]
    qT <- apply(q_post[unkZ,,i,drop=FALSE], 1,
                function(x) ifelse(all(is.na(x)), NA, prod(x, na.rm=TRUE)))
    psi_con <- psi * qT / (psi * qT + (1-psi))
    to_sim <- unkZ[!is.na(psi_con)]
    z_post[to_sim, i] <- rbinom(length(to_sim),1,psi_con[!is.na(psi_con)])
  }

  t(z_post)
})


setMethod("sim_y", "ubmsFitOccu", function(object, samples, re.form, z=NULL, ...){
  nsamp <- length(samples)
  M <- nrow(object@response@y)
  J <- object@response@max_obs
  T <- object@response@max_primary

  z <- process_z(object, samples, re.form, z)
  p <- t(sim_lp(object, submodel="det", transform=TRUE, newdata=NULL,
                samples=samples, re.form=re.form))
  p[object@response@missing] <- NA

  zp <- z[rep(1:nrow(z), each=J),,drop=FALSE] * p

  y_sim <- suppressWarnings(rbinom(M*J*T*nsamp, 1, as.vector(zp)))
  matrix(y_sim, nrow=nsamp, ncol=M*J*T, byrow=TRUE)
})
