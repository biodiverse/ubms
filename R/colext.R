#' Fit the MacKenzie et al. (2003) Dynamic Occupancy Model
#'
#' This function fits the dynamic occupancy model of
#' MacKenzie et al. (2003).
#'
#' @param psiformula Right-hand sided formula for the initial probability of
#'                   occupancy at each site
#' @param gammaformula Right-hand sided formula for colonization probability
#' @param epsilonformula Right-hand sided formula for extinction probability
#' @param pformula Right-hand sided formula for detection probability
#' @param data A \code{\link{unmarkedMultFrame}} object
#' @param ... Arguments passed to the \code{\link{stan}} call, such as 
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitColext} object describing the model fit.
#'
#' @references MacKenzie DI, Nicholas JD, Hines JE, Knutson MG, Franklin AB.
#'             2003. Ecology 84: 2200-2207.
#'
#' @seealso \code{\link{colext}}, \code{\link{unmarkedMultFrame}}
#' @export
stan_colext <- function(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, 
                        pformula = ~1, data, ...){

  umf <- process_umf(data)

  response <- ubmsResponse(getY(umf), "binomial", "binomial", umf@numPrimary)
  state <- ubmsSubmodel("Occupancy", "state", siteCovs(umf), psiformula, "plogis")
  col <- ubmsSubmodel("Colonization", "col", yearlySiteCovs(umf), 
                      gammaformula, "plogis", transition=TRUE)
  ext <- ubmsSubmodel("Extinction", "ext", yearlySiteCovs(umf), 
                      epsilonformula, "plogis", transition=TRUE)
  det <- ubmsSubmodel("Detection", "det", obsCovs(umf), pformula, "plogis")
  submodels <- ubmsSubmodelList(state, col, ext, det)

  ubmsFit("colext", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include fit.R 
setClass("ubmsFitColext", contains = "ubmsFit")


#Method to simulate residuals--------------------------------------------------

#' @include residuals.R


#Goodness-of-fit---------------------------------------------------------------

#' @describeIn gof 
#' Applies the MacKenzie-Bailey chi-square goodness of fit test for
#' ocupancy models (MacKenzie and Bailey 2004).
#' @references MacKenzie, D. I., & Bailey, L. L. (2004). Assessing the 
#'  fit of site-occupancy models. Journal of Agricultural, Biological, 
#'  and Environmental Statistics, 9(3), 300-318.
#' @include gof.R

#Methods to simulate posterior predictive distributions------------------------

#' @include posterior_predict.R
