#' Predict parameter values from a fitted model
#'
#' This method generates predicted parameter values for the original dataset
#' or a new dataset using the posterior distribution. Standard deviation and 
#' a customizable uncertainty interval are also calculated.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param type Submodel type to predict from, for example \code{"det"}
#' @param newdata Optional data frame of covariates to generate
#'   predictions from. If not provided (the default), predictions are 
#'   generated from the original data
#' @param transform If \code{TRUE}, back-transform the predictions to their
#'   original scale
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'   effects") are included. If \code{NA}, they are ignored
#' @param level Probability mass to include in the uncertainty interval
#' @param ... Currently ignored
#'
#' @return A data frame with one row per prediction and four columns: 
#'   1) Predicted point estimates (posterior means),
#'   2) Standard deviation of the posterior,
#'   3-4) Lower and upper bounds of the specified uncertainty interval
#'
#'   For parameters with more than one dimension, the rows are in site-major
#'   order, or site-year-observation for dynamic models.
#'
#' @aliases predict
#' @method predict ubmsFit
#' @seealso posterior_linpred, posterior_interval
#' @include fit.R
#' @importFrom unmarked predict
#' @export
setMethod("predict", "ubmsFit",
          function(object, type, newdata=NULL, transform=TRUE, 
                   re.form=NULL, level=0.95, ...){
  
  samples <- 1:nsamples(object)
  q <- c((1-level)/2, level+(1-level)/2)
  
  lp <- sim_lp(object, type=type, transform=transform, newdata=newdata,
               samples=samples, re.form=re.form)
  
  stats <- apply(lp, 2, function(x){
        quant <- stats::quantile(x, q, na.rm=TRUE)
        c(Predicted = mean(x), SD = stats::sd(x), quant)
      })

  as.data.frame(t(stats)) 
})
