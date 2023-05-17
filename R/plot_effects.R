setGeneric("plot_effects", function(object, ...) standardGeneric("plot_effects"))

#' Plot Marginal Effects of Covariates
#'
#' Generates marginal fixed effects plots of one or more covariates from a
#' \code{ubmsFit} submodel. For each plot, the focal covariate is allowed to
#' vary across its range (or possible discrete values, for a factor), while
#' the other covariates are held at their medians or reference levels.
#' Random effects are ignored.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to get plots for, for example \code{"det"}
#' @param covariate Plot a specific covariate; provide the name as a string
#' @param level Probability mass to include in the uncertainty interval
#' @param draws Number of draws from the posterior to use when generating the
#'  plots. If fewer than \code{draws} are available, they are all used
#' @param smooth Smoother span (f) value used for LOWESS smoothing of the
#'  upper and lower credible interval bounds for a continuous covariate.
#'  No smoothing is done if \code{NULL}. A reasonable value to try is 0.05.
#'  The larger the value, the smoother the curve. As with all smoothing, use
#'  with caution
#' @param ... Currently ignored
#'
#' @return A \code{ggplot} if a single covariate is plotted, or an object
#'  of class \code{grob} if there are multiple covariates/panels
#'
#' @aliases plot_effects plot_marginal
#' @include fit.R
#' @importFrom grid textGrob gpar
#' @importFrom ggplot2 geom_errorbar
#' @export
setMethod("plot_effects", "ubmsFit", function(object, submodel, covariate=NULL,
                                               level=0.95, draws=1000, smooth=NULL, ...){

  sm <- object[submodel]

  if(!is.null(covariate)){
    plots <- list(marginal_covariate_plot(object, submodel, covariate, level,
                                          draws, smooth))
  } else {
    vars <- all.vars(remove_offset(lme4::nobars(sm@formula)))
    if(length(vars) == 0) stop("No covariates in this submodel")
    plots <- lapply(vars, function(x){
                    marginal_covariate_plot(object, submodel, x, level,
                                            draws, smooth)
                    })
  }

  nplots <- length(plots)
  if(nplots == 1){
    out_plot <- plots[[1]] +
      theme(axis.title.y=ggplot2::element_text(angle=90)) +
      ggplot2::labs(y = sm@name)
    return(out_plot)
  }

  dims <- c(ceiling(sqrt(nplots)), round(sqrt(nplots)))
  gridExtra::grid.arrange(grobs=plots, nrow=dims[1], ncol=dims[2],
                          left=grid::textGrob(sm@name, rot=90, vjust=0.5,
                          gp=grid::gpar(fontsize=14)))
})

setGeneric("plot_marginal", function(object, ...) standardGeneric("plot_marginal"))

#' @rdname plot_effects-ubmsFit-method
#' @export
setMethod("plot_marginal", "ubmsFit", function(object, submodel, covariate=NULL,
                                               level=0.95, draws=1000, smooth=NULL, ...){
  plot_effects(object, submodel, covariate, level, draws, smooth, ...)
})

marginal_covariate_plot <- function(object, submodel, covariate, level=0.95,
                                    draws, smooth=NULL){
  sm <- object[submodel]
  quant <- c((1-level)/2, level+(1-level)/2)
  vars <- all.vars(lme4::nobars(sm@formula))
  stopifnot(covariate %in% vars)
  is_factor <- inherits(sm@data[[covariate]], "factor")
  if(is_factor){
    return(marg_factor_plot(object, submodel, covariate, quant, draws))
  }
  marg_numeric_plot(object, submodel, covariate, quant, draws, smooth)
}

marg_numeric_plot <- function(object, submodel, covariate, quant,
                              draws, smooth=NULL){
  sm <- object[submodel]
  samples <- get_samples(object, draws)
  newdata <- get_baseline_df(sm)[rep(1, 1000),,drop=FALSE]
  var_range <- range(sm@data[[covariate]], na.rm=TRUE)
  newdata[[covariate]] <- seq(var_range[1], var_range[2], length.out=1000)

  plot_df <- get_margplot_data(object, submodel, covariate, quant,
                               samples, newdata)

  # Smooth upper and lower CI bounds if requested
  if(!is.null(smooth)){
    plot_df$lower <- stats::lowess(plot_df$lower, f=smooth)[[2]]
    plot_df$upper <- stats::lowess(plot_df$upper, f=smooth)[[2]]
  }
  
  covariate <- sym("covariate"); mn <- sym("mn")
  lower <- sym("lower"); upper <- sym("upper")
  ggplot(data=plot_df, aes(x={{covariate}}, y={{mn}})) +
    geom_ribbon(aes(ymin={{lower}}, ymax={{upper}}), alpha=0.3) +
    geom_line() +
    labs(x = covariate, y = sm@name) +
    plot_theme() +
    theme(axis.title.y=element_blank())
}

marg_factor_plot <- function(object, submodel, covariate, quant, draws){
  sm <- object[submodel]
  samples <- get_samples(object, draws)
  nlev <- nlevels(sm@data[[covariate]])
  newdata <- get_baseline_df(sm)[rep(1, nlev),,drop=FALSE]
  newdata[[covariate]] <- factor(levels(sm@data[[covariate]]),
                                 levels=levels(sm@data[[covariate]]))

  plot_df <- get_margplot_data(object, submodel, covariate, quant,
                               samples, newdata)
  
  covariate <- sym("covariate"); mn <- sym("mn")
  lower <- sym("lower"); upper <- sym("upper")
  ggplot(data=plot_df, aes(x={{covariate}}, y={{mn}})) +
    geom_errorbar(aes(ymin={{lower}}, ymax={{upper}}), width=0.4) +
    geom_point(size=2) +
    labs(x = covariate, y = sm@name) +
    plot_theme() +
    theme(axis.title.y=element_blank())
}

get_baseline_df <- function(submodel){
  vars <- all.vars(lme4::nobars(submodel@formula))
  out <- lapply(vars, function(x, data){
                 if(col_is_factor(x, data)){
                   return(factor(levels(data[[x]])[1], levels=levels(data[[x]])))
                 }
                 return(stats::median(data[[x]], na.rm=TRUE))
                 }, data=submodel@data)
  out <- as.data.frame(out)
  names(out) <- vars
  out
}

col_is_factor <- function(x, data) inherits(data[[x]], "factor")

get_margplot_data <- function(object, submodel, covariate, quant,
                              samples, newdata){
  post <- sim_lp(object, submodel, transform=TRUE, newdata=newdata,
                        samples=samples, re.form=NA)
  mn <- apply(post, 2, mean, na.rm=TRUE)
  bounds <- t(apply(post, 2, stats::quantile, quant, na.rm=TRUE))
  colnames(bounds) <- c("lower", "upper")
  data.frame(covariate = newdata[[covariate]], mn = mn, bounds)
}

