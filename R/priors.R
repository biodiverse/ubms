#' Prior distributions
#'
#' Specify prior distributions and associated parameters for use in
#' \code{ubms} models.
#'
#' @param location The mean of the distribution. If setting the priors for
#'  regression coefficients, this can be a single value, or multiple values,
#'  one per coefficient
#' @param scale The standard deviation of the distribution. If setting the priors
#'  for regression coefficients, this can be a single value, or multiple values,
#'  one per coefficient
#' @param lower The lower bound for the uniform distribution
#' @param upper The upper bound for the uniform distribution
#' @param df The number of degrees of freedom for the Student-t distribution
#' @param shape The gamma distribution shape parameter
#' @param rate The gamma distribution rate parameter (1/scale)
#' @param autoscale If \code{TRUE}, ubms will automatically adjust priors
#'  for each regression coefficient relative to its corresponding covariate x.
#'  Specifically, the prior for a given coefficient will be divided by
#'  sd(x). This helps account for covariates with very different magnitudes
#'  in the same model. If your data are already standardized (e.g. with use of
#'  \code{scale()}), this will have minimal effect as sd(x) will be
#'  approximately 1. Standardizing your covariates is highly recommended.
#'
#' @name priors
#'
#' @return A \code{list} containing prior settings used internally by \code{ubms}.
#'
#' @examples
#' normal()
#'
NULL

#' @rdname priors
#' @export
normal <- function(location=0, scale=2.5, autoscale=TRUE){
  stopifnot(all(scale > 0))
  if((length(location) > 1) & (length(scale) > 1)){
    stopifnot(length(location) == length(scale))
  }
  list(dist=1, par1=location, par2=scale, par3=0, autoscale=autoscale)
}

#' @rdname priors
#' @export
uniform <- function(lower=-5, upper=5){
  stopifnot(length(lower) == length(upper))
  stopifnot(all(lower < upper))
  list(dist=2, par1=lower, par2=upper, par3=0, autoscale=FALSE)
}

#' @rdname priors
#' @export
student_t <- function(df=1, location=0, scale=2.5, autoscale=TRUE){
  stopifnot(all(scale > 0))
  stopifnot(all(df > 0))
  if((length(location) > 1) & (length(scale) > 1)){
    stopifnot(length(location) == length(scale))
  }
  list(dist=3, par1=location, par2=scale, par3=df, autoscale=autoscale)
}

#' @rdname priors
#' @export
logistic <- function(location=0, scale=1){
    stopifnot(all(scale > 0))
  if((length(location) > 1) & (length(scale) > 1)){
    stopifnot(length(location) == length(scale))
  }
  list(dist=4, par1=location, par2=scale, par3=0, autoscale=FALSE)
}

#' @rdname priors
#' @export
cauchy <- function(location=0, scale=2.5, autoscale=TRUE){
  stopifnot(all(scale > 0))
  if((length(location) > 1) & (length(scale) > 1)){
    stopifnot(length(location) == length(scale))
  }
  list(dist=3, par1=location, par2=scale, par3=1, autoscale=autoscale)
}

#' @rdname priors
#' @export
gamma <- function(shape=1, rate=1){
  stopifnot(length(shape==1) & length(rate==1))
  stopifnot(shape > 0 & rate > 0)
  list(dist=5, par1=shape, par2=rate, par3=0, autoscale=FALSE)
}

null_prior <- function(){
  list(dist=0, par1=0, par2=0, par3=0, autoscale=FALSE)
}

expand_prior <- function(prior, np){
  rep_prior <- lapply(prior[c("par1","par2","par3")], function(x){
    if(length(x) > 1){
      stopifnot(length(x) == np)
    } else {
      x <- rep(x, np)
    }
    x
  })
  prior[c("par1","par2","par3")] <- rep_prior
  prior
}

autoscale_prior <- function(prior, Xmat){
  if(!prior$autoscale) return(prior)
  # par2 is always the 'scale' parameter
  stopifnot(ncol(Xmat) == length(prior$par2))

  for (i in 1:ncol(Xmat)){
    # skip if dummy variable
    if(all(unique(stats::na.omit(Xmat[,i])) %in% c(0,1))) next
    prior$par2[i] <- prior$par2[i] * 1/stats::sd(Xmat[,i], na.rm=TRUE)
  }
  prior
}

process_coef_prior <- function(prior, Xmat){
  if(prior$dist == 5){
    stop("gamma distribution can only be used with prior_sigma")
  }
  # remove intercept
  if("(Intercept)" %in% colnames(Xmat)) Xmat <- Xmat[,-1,drop=FALSE]
  # if no coefs
  if(ncol(Xmat)==0){
    prior$dist <- 0
    prior$par1 <- prior$par2 <- prior$par3 <- NA
    return(prior)
  }
  # expand to correct length
  prior <- expand_prior(prior, ncol(Xmat))
  # adjust scale if requested
  autoscale_prior(prior, Xmat)
}

process_int_prior <- function(prior, Xmat){
  stopifnot(length(unlist(prior[c("par1","par2","par3")])) == 3)
  if(prior$dist == 5){
    stop("gamma distribution can only be used with prior_sigma")
  }
  prior$autoscale <- FALSE
  # If there's an intercept, don't do anything
  if("(Intercept)" %in% colnames(Xmat)) return(prior)

  prior$dist <- 0
  prior$par1 <- prior$par2 <- prior$par3 <- NA
  prior
}

process_sigma_prior <- function(prior){
  prior
}

setGeneric("process_priors", function(submodel, ...) standardGeneric("process_priors"))

#' @include submodel.R
setMethod("process_priors", "ubmsSubmodel", function(submodel){
  Xmat <- model.matrix(submodel)
  coef_prior <- process_coef_prior(submodel@prior_coef, Xmat)
  int_prior <- process_int_prior(submodel@prior_intercept, Xmat)
  sigma_prior <- process_sigma_prior(submodel@prior_sigma)
  out <- mapply(function(x,y,z) c(x,y,z), int_prior, coef_prior, sigma_prior, SIMPLIFY=FALSE)
  prior_pars <- matrix(unlist(out[paste0("par",1:3)]), nrow=3, byrow=TRUE)
  drop_cols <- apply(prior_pars, 2, function(x) any(is.na(x)))
  prior_pars <- prior_pars[,!drop_cols,drop=FALSE]
  list(prior_dist = out$dist, prior_pars = prior_pars)
})

setMethod("process_priors", "ubmsSubmodelScalar", function(submodel){
  Xmat <- model.matrix(submodel)
  int_prior <- process_int_prior(submodel@prior_intercept, Xmat)
  prior_pars <- matrix(c(unlist(int_prior[paste0("par",1:3)]), rep(0,3)),
                       nrow=3)
  list(prior_dist = c(int_prior$dist,0,0), prior_pars = prior_pars)
})
