#' Extract Model Residuals
#'
#' Extract residuals for a given submodel from a \code{ubmsFit} object.
#' Residuals are calculated separately for each submodel
#' using the posterior predictive distribution of the latent state z, 
#' following Wright et al. (2019). 
#' 
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to get residuals for, for example \code{"det"}
#' @param draws An integer indicating the number of draws to return. The
#'   default and maximum number of draws is the size of the posterior sample.
#' @param ... Currently ignored
#'
#' @return A matrix of residual values with dimension \code{draws} by 
#'   observations. Note that calculation of residuals
#'   for the detection submodel is conditional on \eqn{z > 0}, so residuals
#'   for an observation in a posterior draw where \eqn{z = 0} are assigned 
#'   value \code{NA} (Wright et al. 2019).
#' 
#' @references Wright, W. J., Irvine, K. M., & Higgs, M. D. (2019). Identifying 
#'   occupancy model inadequacies: can residuals separately assess detection 
#'   and presence? Ecology, 100(6), e02703.
#' 
#' @include fit.R
#' @importFrom stats residuals
#' @export
setMethod("residuals", "ubmsFit", function(object, submodel, draws=NULL, ...){
  samples <- get_samples(object, draws)
  sim_res(object, submodel, samples)
})


#Internal function for calculating residuals
setGeneric("sim_res", function(object, ...) standardGeneric("sim_res"))

setMethod("sim_res", "ubmsFit", function(object, ...){
  stop("No available method for this fit type", call.=FALSE)
})

#' @include occu.R
setMethod("sim_res", "ubmsFitOccu", function(object, submodel, samples, ...){

  lp <- sim_lp(object, submodel, samples=samples, transform=TRUE,
               newdata=NULL, re.form=NULL)
  z <- sim_z(object, samples=samples, re.form=NULL)

  if(identical(submodel, "state")){
    res <- z - lp
  } else if(identical(submodel, "det")){
    y <- object@data@y
    J <- ncol(y)
    ylong <- as.vector(t(y))
    zrep <- z[, rep(1:ncol(z), each=J)]
    z1_mask <- zrep == 1
    res <- matrix(rep(ylong, each=nrow(lp)), nrow=nrow(lp)) - lp
    res[!z1_mask] <- NA #residuals conditional on z = 1
  }
  res
})


setMethod("sim_res", "ubmsFitOccuRN", function(object, submodel, samples, ...){

  lp <- sim_lp(object, submodel, samples=samples, transform=TRUE,
               newdata=NULL, re.form=NULL)
  z <- sim_z(object, samples=samples, re.form=NULL)

  if(identical(submodel, "state")){
    res <- z - lp
  } else if(identical(submodel, "det")){
    y <- object@data@y
    J <- ncol(y)
    ylong <- as.vector(t(y))
    zrep <- z[, rep(1:ncol(z), each=J)]
    p <- 1 - (1 - lp)^zrep
    z1_mask <- zrep > 0
    res <- matrix(rep(ylong, each=nrow(p)), nrow=nrow(p)) - p
    res[!z1_mask] <- NA #residuals conditional on z > 0
  }
  res
})


# Plotting functions

setGeneric("plot_residuals", function(object, ...) standardGeneric("plot_residuals"))

#' Plot Model Residuals
#'
#' Plot residuals for a submodel from a \code{ubmsFit} object, for multiple 
#' posterior draws. By default, residuals are plotted against fitted values.
#' When the submodel has a binomial response (e.g., detection models), regular 
#' residual plots are not typically informative. Instead, the residuals and 
#' fitted values are divided into bins based on fitted value and the averages 
#' are plotted. For a count response (e.g., Poisson), Pearson residuals are 
#' calculated. To plot residuals against values of a particular covariate instead
#' of the fitted values, supply the name of the covariate (as a string) to the 
#' \code{covariate} argument. 
#' 
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to plot residuals for, for example \code{"det"}
#' @param covariate If specified, plot residuals against values of a covariate.
#'   Covariate name should be provided as a string. If \code{NULL}, 
#'   residuals are plotted against predicted values.
#' @param draws An integer indicating the number of posterior draws to use. 
#'   Separate plots are generated for each draw, so this number should be
#'   relatively small. The default and maximum number of draws is the size of 
#'   the posterior sample.
#' @param nbins For submodels with a binomial response, manually set the number 
#'   of bins to use
#' @param ... Currently ignored
#'
#' @return A \code{ggplot} of residuals vs. fitted values or covariate values, 
#'   with one panel per posterior draw. For binned residual plots, the shaded area
#'   represents plus/minus two standard deviations around the mean residual. 
#'   If the model is true, we would expect about 95% of the binned residuals to
#'   fall within this area.
#' 
#' @aliases plot_residuals
#' @seealso \code{\link{residuals}} 
#'
#' @export
setMethod("plot_residuals", "ubmsFit", function(object, submodel, covariate=NULL, 
                                                draws=9, nbins=NULL, ...){
  
  if(identical(object[submodel]@link, "plogis")){
    return(plot_binned_residuals(object, submodel, covariate, draws, nbins))
  }
  plot_pearson_residuals(object, submodel, covariate, draws)

})

#' @importFrom ggplot2 facet_wrap geom_hline
plot_pearson_residuals <- function(object, submodel, covariate=NULL, draws=9){

  samples <- get_samples(object, draws)
  res <- sim_res(object, submodel, samples)
  lp <- sim_lp(object, submodel, samples=samples, transform=TRUE,
                      newdata=NULL, re.form=NULL)
  res <- res / sqrt(lp) #Pearson residuals
  if(is.null(covariate)){
    x <- lp
    xlab <- "Predicted value"
  } else {
    x <- object[submodel]@data[[covariate]]
    x <- matrix(rep(x, each=nrow(res)), nrow=nrow(res))
    xlab <- paste(covariate, "value")
  }

  pl_dat <- lapply(1:draws, function(i){
              data.frame(x = x[i,], y= res[i,], ind=i)
            })
  pl_dat <- do.call("rbind", pl_dat)
  
  ggplot(data=pl_dat, aes_string(x="x", y="y")) +
    geom_hline(aes(yintercept=0), linetype=2) +
    geom_point() +
    facet_wrap("ind") +
    ggtitle(paste(object[submodel]@name, "submodel residuals plot")) +
    labs(x=xlab, y="Pearson residual") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.background=element_blank(),
          strip.text=element_blank())
}

#' @importFrom ggplot2 geom_ribbon geom_line
plot_binned_residuals <- function(object, submodel, covariate=NULL, draws=9, 
                                  nbins=NULL){
  samples <- get_samples(object, draws)
  res <- sim_res(object, submodel, samples)
  if(is.null(covariate)){
    x <- sim_lp(object, submodel, samples=samples, transform=TRUE,
                      newdata=NULL, re.form=NULL)
    xlab <- "Mean predicted value"
  } else {
    x <- object[submodel]@data[[covariate]]
    x <- matrix(rep(x, each=nrow(res)), nrow=nrow(res))
    xlab <- paste("Mean", covariate, "value")
  }

  pl_dat <- lapply(1:length(samples), function(i){
                  get_binned_residuals(x[i,], res[i,], i, nbins)})
  pl_dat <- do.call("rbind", pl_dat)

  ggplot(data=pl_dat, aes_string(x="xbar", y="ybar")) +
    geom_ribbon(aes_string(ymin="y_lo", ymax="y_hi"), alpha=0.1) +
    geom_hline(aes(yintercept=0), linetype=2) +
    geom_line(aes_string(y="y_hi"), col='gray', size=1.1) +
    geom_line(aes_string(y="y_lo"), col='gray', size=1.1)+
    geom_point() +
    facet_wrap("ind") +
    ggtitle(paste(object[submodel]@name, "submodel residuals plot")) +
    labs(x=xlab, y="Mean binned residual") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.background=element_blank(),
          strip.text=element_blank()) 
}

get_binned_residuals <- function(x, y, ind, nbins=NULL, ...){
  
  na <- is.na(x) | is.na(y)
  x <- x[!na]
  y <- y[!na]
  if(is.null(nbins)) nbins <- sqrt(length(x))

  bad_break <- TRUE
  while(bad_break){
    tryCatch({
      if(nbins < 4) stop("Couldn't find working breakpoints", call.=FALSE)      
      breaks.index <- floor(length(x)*(1:(nbins-1))/nbins)
      breaks <- c (-Inf, sort(x)[breaks.index], Inf)
      x.binned <- as.numeric(cut(x, breaks))
      bad_break <- FALSE
    }, error=function(e){
        nbins <<- nbins - 1
      }
    )
  }
 
  output <- NULL
  for (i in 1:nbins){
    items <- (1:length(x))[x.binned==i]
    x.range <- range(x[items], na.rm=T)
    xbar <- mean(x[items], na.rm=T)
    ybar <- mean(y[items], na.rm=T)
    n <- length(stats::na.omit(items))
    sdev <- stats::sd(y[items], na.rm=T)
    se <- sdev/sqrt(n)
    output <- rbind(output, c(xbar, ybar, n, x.range, -1.96*se, 1.96*se))
  }
  colnames(output) <- c("xbar", "ybar", "n", "x.lo", "x.hi", "y_lo","y_hi")
  output <- as.data.frame(output)
  output$ind <- ind
  output
}
