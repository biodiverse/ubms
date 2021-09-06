check_prior_input <- function(input){
  #stopifnot(length(input) %in% c(1,2))
  stopifnot(length(input) == 2)
  #if(length(input) == 1){
  #  check_names <- "coef"
  #} else {
  check_names <- c("intercept","coef")
  #}
  if(!all(check_names %in% names(input))){
    stop(paste0("List of priors for each submodel must have named elements: ",
                paste(check_names, collapse=", ")),
         call.=FALSE)
  }
}

get_prior_info <- function(input, Xmat){

  check_prior_input(input)

  locations <- array(rep(input$coef$location, ncol(Xmat)))
  scales <- array(rep(input$coef$scale, ncol(Xmat)))
  prior_type <- c(1,1) # only normal supported at moment

  if(input$coef$autoscale){
    for (i in 1:ncol(Xmat)){

      if(all(unique(Xmat[,i]) %in% c(0,1))){
        next
      }
      scales[i] <- scales[i] * 1/sd(Xmat[,i])
    }
  }

  if(!is.null(input$intercept) && colnames(Xmat)[1] == "(Intercept)"){
    locations[1] <- input$intercept$location
    scales[1] <- input$intercept$scale
  }

  list(locations=locations, scales=scales, prior_type=prior_type)

}

check_missing_prior <- function(priors, type){
  if(is.null(priors)){
    stop(paste("Missing prior information for submodel",type), call.=FALSE)
  }
  priors
}

#' Specify normal prior distribution
#'
#' Controls settings for normal prior distributions on an intercept or
#' on coefficients for a submodel.
#'
#' @param location The mean of the distribution
#' @param scale The standard deviation of the distribution
#' @param autoscale If \code{TRUE}, ubms will automatically adjust priors
#'  for each regression coefficient relative to its corresponding covariate x.
#'  Specifically, the prior for a given coefficient will be divided by
#'  sd(x). This helps account for covariates with very different magnitudes
#'  in the same model. If your data are already standardized (e.g. with use of
#'  \code{scale()}), this will have minimal effect as sd(x) will be
#'  approximately 1. Standardizing your covariates is highly recommended.
#'
#' @return A \code{list} containing prior settings used internally by \code{ubms}.
#'
#' @examples
#' normal()
#'
#' @export
normal <- function(location=0, scale=10, autoscale=TRUE){
  stopifnot(scale > 0)
  list(dist="normal", location=location, scale=scale, autoscale=autoscale)
}

#' Default priors for single-season models
#'
#' Convenience function providing default priors for most \code{ubms} models.
#' Also useful for seeing correct list structure for custom priors.
#'
#' @return A list of lists of default prior specifications, one entry
#'  per submodel.
#'
#' @examples
#' default_priors()
#'
#' @export
default_priors <- function(){
  out <- lapply(1:2, function(x){
           list(intercept=normal(0, 10, autoscale=TRUE),
           coef=normal(0, 2.5, autoscale=TRUE))
         })
  names(out) <- c("state","det")
  out
}
