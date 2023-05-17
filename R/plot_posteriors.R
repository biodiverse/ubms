setGeneric("plot_posteriors", function(object, ...) standardGeneric("plot_posteriors"))

#' Plot Posterior Distributions
#'
#' Plot posterior distributions for selected parameters. Posteriors can be
#' represented as density plots or histograms.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param pars A character vector of parameter names to include in the plot
#'  Look at \code{names(object@stanfit)} for the complete list of possible
#'  parameter names. If \code{NULL}, posteriors are shown for all parameters
#'  in the model summary output
#' @param density If \code{TRUE}, show posteriors as density plots (one per
#'  chain). If \code{FALSE}, show posteriors as histograms of samples from
#'  all chains combined
#' @param ... Arguments passed to \code{ggplot2::stat_density} for density
#'  plots, or \code{ggplot2::geom_histogram} for histograms. For example, you
#'  can supply argument \code{bins} to control the number of histogram bins
#'
#' @return A \code{ggplot}
#'
#' @aliases plot_posteriors
#' @include fit.R
#' @importFrom ggplot2 stat_density labs facet_wrap geom_histogram
#' @export
setMethod("plot_posteriors", "ubmsFit", function(object, pars=NULL, density=FALSE, ...){

  if(is.null(pars)) pars <- names(object)[!grepl("b_", names(object))]

  not_in_mod <- ! pars %in% names(object@stanfit)
  if(any(not_in_mod)){
    stop(paste0("Parameter(s) ", pars[not_in_mod], " were not found in model",
                collapse=", "), call.=FALSE)
  }

  # work around problem in rstan::check_pars
  object@stanfit@sim$fnames_oi <- gsub(" ", "", object@stanfit@sim$fnames_oi)
  pars <- gsub(" ", "", pars)

  samp <- extract(object, pars, permuted=FALSE)
  pars <- dimnames(samp)$parameters
  samp <- lapply(1:dim(samp)[3], function(i) samp[,,i])

  plot_dat <- lapply(1:length(samp), function(i){
    samp_i <- samp[[i]]
    nsamp <- nrow(samp_i)
    data.frame(parameter=pars[i],
               samples=as.vector(samp_i),
               chain=rep(as.character(1:ncol(samp_i)), each=nsamp))
  })
  plot_dat <- do.call(rbind, plot_dat)
  plot_dat$parameter <- factor(plot_dat$parameter, levels=pars)

  if(density){
    samples <- sym("samples"); chain <- sym("chain")
    out <- ggplot(data=plot_dat) +
      stat_density(aes(x={{samples}}, col={{chain}}), geom="line",
                   position="identity", ...) +
      labs(x="Value", y="Density") +
      facet_wrap("parameter", scales="free") +
      plot_theme() +
      theme(strip.text=element_text(size=12),
            strip.background=element_rect("white"))
  } else {
    samples <- sym("samples")
    out <- ggplot(data=plot_dat) +
      geom_histogram(aes(x={{samples}}),fill='gray90',col='black', ...) +
      labs(x="Value", y="Count") +
      facet_wrap("parameter", scales="free") +
      plot_theme() +
      theme(strip.text=element_text(size=12),
            strip.background=element_rect("white"))
  }
  out
})
