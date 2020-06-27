setGeneric("plot_marginal", function(object, ...) standardGeneric("plot_marginal"))

#' Plot Marginal Effects of Covariates
#'
#' Generates marginal fixed effects plots of one or more covariates from a
#' \code{ubmsFit} submodel. For each plot, the focal covariate is allowed to
#' vary across its range (or possible discrete values, for a factor), while
#' the other covariates are held at their means or reference levels.
#' Random effects are ignored.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel Submodel to get plots for, for example \code{"det"}
#' @param covariate Plot a specific covariate; provide the name as a string
#' @param level Probability mass to include in the uncertainty interval
#' @param ... Currently ignored
#'
#' @return A \code{ggplot}, potentially with multiple panels, one per covariate
#'
#' @aliases plot_marginal
#' @include fit.R
#' @importFrom grid textGrob gpar
#' @importFrom ggplot2 geom_errorbar
#' @export
setMethod("plot_marginal", "ubmsFit", function(object, submodel, covariate=NULL,
                                               level=0.95, ...){

  sm <- object[submodel]

  if(!is.null(covariate)){
    plots <- list(marginal_covariate_plot(object, submodel, covariate, level))
  } else {
    vars <- all.vars(lme4::nobars(sm@formula))
    if(length(vars) == 0) stop("No covariates in this submodel")
    plots <- lapply(vars, function(x){
                    marginal_covariate_plot(object, submodel, x, level)
                    })
  }

  dims <- c(1,1)
  nplots <- length(plots)
  if(nplots > 1){
    dims <- c(ceiling(sqrt(nplots)), round(sqrt(nplots)))
  }

  gridExtra::grid.arrange(grobs=plots, nrow=dims[1], ncol=dims[2],
                          left=grid::textGrob(sm@name, rot=90, vjust=0.5,
                          gp=grid::gpar(fontsize=14)))
})

marginal_covariate_plot <- function(object, submodel, covariate, level=0.95){
  sm <- object[submodel]
  quant <- c((1-level)/2, level+(1-level)/2)
  vars <- all.vars(lme4::nobars(sm@formula))
  stopifnot(covariate %in% vars)
  is_factor <- inherits(sm@data[[covariate]], "factor")
  if(is_factor){
    return(marg_factor_plot(object, submodel, covariate, quant))
  }
  marg_numeric_plot(object, submodel, covariate, quant)
}

marg_numeric_plot <- function(object, submodel, covariate, quant){
  sm <- object[submodel]
  samples <- get_samples(object, 1000)
  newdata <- get_mean_df(sm)[rep(1, 1000),,drop=FALSE]
  var_range <- range(sm@data[[covariate]], na.rm=TRUE)
  newdata[[covariate]] <- seq(var_range[1], var_range[2], length.out=1000)

  plot_df <- get_margplot_data(object, submodel, covariate, quant,
                               samples, newdata)

  ggplot(data=plot_df, aes_string(x="covariate", y="mn")) +
    geom_ribbon(aes_string(ymin="lower", ymax="upper"), alpha=0.3) +
    geom_line() +
    labs(x = covariate, y = sm@name) +
    plot_theme() +
    theme(axis.title.y=element_blank())
}

marg_factor_plot <- function(object, submodel, covariate, quant){
  sm <- object[submodel]
  samples <- get_samples(object, 1000)
  nlev <- nlevels(sm@data[[covariate]])
  newdata <- get_mean_df(sm)[rep(1, nlev),,drop=FALSE]
  newdata[[covariate]] <- levels(sm@data[[covariate]])

  plot_df <- get_margplot_data(object, submodel, covariate, quant,
                               samples, newdata)

  ggplot(data=plot_df, aes_string(x="covariate", y="mn")) +
    geom_errorbar(aes_string(ymin="lower", ymax="upper"), width=0.4) +
    geom_point(size=2) +
    labs(x = covariate, y = sm@name) +
    plot_theme() +
    theme(axis.title.y=element_blank())
}

get_mean_df <- function(submodel){
  vars <- all.vars(lme4::nobars(submodel@formula))
  out <- lapply(vars, function(x, data){
                   ifelse(col_is_factor(x, data),
                          levels(data[[x]])[1], mean(data[[x]], na.rm=TRUE))},
                    data=submodel@data)
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

