#' Predict parameter values from a fitted model
#'
#' This method generates predicted parameter values for the original dataset
#' or a new dataset using the posterior distribution. Standard deviation and
#' a customizable uncertainty interval are also calculated.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to predict from, for example \code{"det"}
#' @param newdata Optional data frame, SpatRaster, or RasterStack of covariates to generate
#'   predictions from. If not provided (the default), predictions are
#'   generated from the original data
#' @param transform If \code{TRUE}, back-transform the predictions to their
#'   original scale
#' @param re.form If \code{NULL}, any estimated group-level parameters ("random
#'   effects") are included. If \code{NA}, they are ignored
#' @param level Probability mass to include in the uncertainty interval
#' @param ... Currently ignored
#'
#' @return If \code{newdata} was a data frame: A data frame with one row per
#'   prediction and four columns:
#'   1) Predicted point estimates (posterior means),
#'   2) Standard deviation of the posterior,
#'   3-4) Lower and upper bounds of the specified uncertainty interval
#'
#'   For parameters with more than one dimension, the rows are in site-major
#'   order, or site-year-observation for dynamic models.
#'
#'   If \code{newdata} was a SpatRaster/RasterStack, returns a SpatRaster/RasterStack 
#'   with four layers corresponding to the four columns above with the same projection
#'   as the original SpatRaster/RasterStack.
#'
#' @aliases predict
#' @method predict ubmsFit
#' @seealso posterior_linpred, posterior_interval
#' @include fit.R
#' @importFrom unmarked predict
#' @export
setMethod("predict", "ubmsFit",
          function(object, submodel, newdata=NULL, transform=TRUE,
                   re.form=NULL, level=0.95, ...){

  if(inherits(newdata, c("RasterLayer", "RasterStack"))){
    return(predict_raster(object, submodel, newdata, transform, re.form, level))
  } else if(inherits(newdata, "SpatRaster")){
    return(predict_terra(object, submodel, newdata, transform, re.form, level))
  }

  samples <- 1:nsamples(object)
  q <- c((1-level)/2, level+(1-level)/2)

  lp <- sim_lp(object, submodel=submodel, transform=transform, newdata=newdata,
               samples=samples, re.form=re.form)

  stats <- apply(lp, 2, function(x){
        quant <- stats::quantile(x, q, na.rm=TRUE)
        c(Predicted = mean(x), SD = stats::sd(x), quant)
      })

  as.data.frame(t(stats))
})

predict_raster <- function(object, submodel, inp_rast, transform, re.form, level){
  if(!requireNamespace("raster", quietly=TRUE)){
    stop('Package "raster" is not installed', call.=FALSE)
  }
  df_dat <- raster::as.data.frame(inp_rast, xy=TRUE)
  out <- cbind(df_dat[,1:2,drop=FALSE],
               predict(object, submodel, newdata=df_dat, transform=transform,
                       re.form=re.form, level=level))
  out <- raster::rasterFromXYZ(out, crs=raster::crs(inp_rast))
  names(out)[3:4] <- c("Lower","Upper")
  out
}

predict_terra <- function(object, submodel, inp_rast, transform, re.form, level){
  if(!requireNamespace("terra", quietly=TRUE)){
    stop('Package "terra" is not installed', call.=FALSE)
  }
  df_dat <- terra::as.data.frame(inp_rast, xy=TRUE)
  out <- cbind(df_dat[,1:2,drop=FALSE],
               predict(object, submodel, newdata=df_dat, transform=transform,
                       re.form=re.form, level=level))
  out <- terra::rast(out, type="xyz", crs=terra::crs(inp_rast))
  names(out)[3:4] <- c("Lower","Upper")
  out
}
