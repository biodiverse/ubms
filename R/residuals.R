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

setMethod("sim_res", "ubmsFit", function(object, submodel, samples,
                                         type=c("absolute", "pearson"), ...){
  type <- match.arg(type, c("absolute", "pearson"))

  fit <- sim_fitted(object, submodel, samples)

  if(identical(submodel, "state")){
    res <- sim_z(object, samples=samples, re.form=NULL) - fit
  } else if(identical(submodel, "det")){
    y <- as_vector(object@response)
    y <- matrix(rep(y, each=nrow(fit)), nrow=nrow(fit))
    res <- y - fit
  }

  if(identical(type, "pearson")) res <- res / sqrt(fit)
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

setMethod("plot_residuals", "ubmsFit", function(object, submodel, covariate=NULL,
                                                draws=9, nbins=NULL, ...){
  type <- ifelse(object[submodel]@link=="plogis", "absolute", "pearson")
  samples <- get_samples(object, draws)
  res <- sim_res(object, submodel, samples, type=type)
  x <- sim_fitted(object, submodel, samples=samples)

  name <- object[submodel]@name
  xlab <- "Predicted value"
  if(!is.null(covariate)){
    x <- object[submodel]@data[[covariate]]
    x <- matrix(rep(x, each=nrow(res)), nrow=nrow(res))
    xlab <- paste(covariate, "value")
  }

  if(identical(type, "absolute")){
    return(plot_binned_residuals(x, res, xlab, name, nbins))
  }
  plot_pearson_residuals(x, res, xlab, name)
})

#' @importFrom ggplot2 facet_wrap geom_hline
plot_pearson_residuals <- function(x, res, xlab, name){
  pl_dat <- lapply(1:nrow(res), function(i){
              data.frame(x = x[i,], y= res[i,], ind=i)
            })
  pl_dat <- do.call("rbind", pl_dat)
  pl_dat <- pl_dat[stats::complete.cases(pl_dat),]
  
  x <- sym("x"); y <- sym("y")
  ggplot(data=pl_dat, aes(x={{x}}, y={{y}})) +
    geom_hline(aes(yintercept=0), linetype=2) +
    geom_point() +
    facet_wrap("ind") +
    ggtitle(paste(name, "submodel residuals plot")) +
    labs(x=xlab, y="Pearson residual") +
    plot_theme()
}

#' @importFrom ggplot2 geom_ribbon geom_line
plot_binned_residuals <- function(x, res, xlab, name, nbins){
  pl_dat <- lapply(1:nrow(res), function(i){
                  get_binned_residuals(x[i,], res[i,], i, nbins)})
  pl_dat <- do.call("rbind", pl_dat)
  pl_dat <- pl_dat[stats::complete.cases(pl_dat),]
  
  x_bar <- sym("x_bar"); y_bar <- sym("y_bar")
  y_lo <- sym("y_lo"); y_hi <- sym("y_hi")
  ggplot(data=pl_dat, aes(x={{x_bar}}, y={{y_bar}})) +
    geom_ribbon(aes(ymin={{y_lo}}, ymax={{y_hi}}), alpha=0.1) +
    geom_hline(aes(yintercept=0), linetype=2) +
    geom_line(aes(y={{y_hi}}), col='gray', linewidth=1.1) +
    geom_line(aes(y={{y_lo}}), col='gray', linewidth=1.1)+
    geom_point() +
    facet_wrap("ind") +
    ggtitle(paste(name, "submodel residuals plot")) +
    labs(x=xlab, y="Mean binned residual") +
    plot_theme()
}

get_binned_residuals <- function(x, y, ind, nbins=NULL, ...){

  na <- is.na(x) | is.na(y)
  x <- x[!na]
  y <- y[!na]

  breaks <- get_breaks(x, nbins)

  output <- lapply(1:breaks$nbins, function(i){
    in_bin <- breaks$x_binned == i
    se <- stats::sd(y[in_bin]) / sqrt(sum(in_bin))
    data.frame(x_bar = mean(x[in_bin]), y_bar = mean(y[in_bin]),
               y_lo = -1.96*se, y_hi= 1.96*se)
  })
  output <- do.call("rbind", output)
  output$ind <- ind
  output
}

get_breaks <- function(x, nbins){

  if(is.null(nbins)) nbins <- ceiling(sqrt(length(x)))

  bad_break <- TRUE
  while(bad_break){
    if(nbins < 4) stop("Couldn't find working breakpoints", call.=FALSE)
    tryCatch({
      breaks_index <- floor(length(x)*(1:(nbins-1))/nbins)
      breaks <- c (-Inf, sort(x)[breaks_index], Inf)
      x_binned <- as.numeric(cut(x, breaks))
      bad_break <- FALSE
    }, error=function(e){
        nbins <<- nbins - 1
      }
    )
  }
  list(nbins=length(unique(x_binned)), x_binned=x_binned)
}
