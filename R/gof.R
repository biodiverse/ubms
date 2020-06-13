setClass("ubmsGOF", slots=c(statistic="character", estimate="numeric",
                            post_pred_p="numeric", samples="data.frame"))

#Constructor
ubmsGOF <- function(stat, samples){
  new("ubmsGOF", statistic=stat, estimate = mean(samples$obs), 
      samples = samples, post_pred_p = mean(samples$sim > samples$obs))
}

setMethod("show", "ubmsGOF", function(object){
  cat(object@statistic, "\n")
  cat(paste0("Point estimate = ", round(object@estimate,3), "\n"))
  cat(paste0("Posterior predictive p = ", round(object@post_pred_p,3),"\n"))
})

#' @importFrom ggplot2 ggplot aes geom_abline geom_point theme_bw labs
#' @importFrom ggplot2 facet_wrap theme element_blank element_text element_rect
#' @importFrom ggplot2 geom_label unit aes_string ggtitle
setMethod("plot", "ubmsGOF", function(x, ...){  
  ppval <- data.frame(lab=paste("P =", round(x@post_pred_p, 3)))
  ggplot(x@samples, aes_string(x="obs", y="sim")) +
    geom_abline(aes(intercept=0, slope=1),size=1.2, col='red') +
    geom_point(alpha=0.4) +
    ggtitle(paste("Posterior predictive check:", x@statistic)) +
    theme_bw() +
    labs(y="Simulated data", x="Observed data") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=14),
          strip.background=element_rect(fill="transparent"),
          strip.text=element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    geom_label(data=ppval, aes_string(x=-Inf, y=Inf, label="lab"),
              hjust=-0.2, vjust=1.4, size=5,
              fill='white', label.size=0, 
              label.padding=unit(0.1, "lines"))
})

#' Check model goodness-of-fit
#' 
#' Goodness-of-fit tests for \code{ubmsFit} models using posterior predictive
#' checks
#' @param object A fitted model of class \code{ubmsFit}
#' @param draws Number of draws from the posterior to use in the check
#' @param quiet If \code{TRUE}, suppress progress bar
#' @param ... Currently ignored
#'
#' @aliases plot,ubmsGOF,ANY-method
#'
#' @return An object of class \code{ubmsGOF} containing statistics calculated
#' from the posterior predictive distribution.
#' @export
setGeneric("gof", function(object, draws=NULL, ...){
             standardGeneric("gof")})
