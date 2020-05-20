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

#' @describeIn gof 
#' Applies the MacKenzie-Bailey chi-square goodness of fit test for
#' ocupancy models (MacKenzie and Bailey 2004).
#' @references MacKenzie, D. I., & Bailey, L. L. (2004). Assessing the 
#'  fit of site-occupancy models. Journal of Agricultural, Biological, 
#'  and Environmental Statistics, 9(3), 300-318.
#' @include occu.R
setMethod("gof", "ubmsFitOccu", function(object, draws=NULL, quiet=FALSE, ...){

  samples <- get_samples(object, draws)
  draws <- length(samples)
  
  psi <- sim_lp(object, transform=TRUE, submodel="state", newdata=NULL, 
                       samples, re.form=NULL)
  p <- sim_lp(object, transform=TRUE, submodel="det", newdata=NULL, 
                     samples, re.form=NULL)

  yobs <- getY(object@data)
  M <- nrow(yobs)
  J <- ncol(yobs)
  ysim <- sim_y(object, samples, re.form=NULL)
  ysim <- array(ysim, c(draws,J,M))
  ysim <- aperm(ysim, c(3,2,1))

  mb_obs <- mb_sim <- rep(NA, draws)
  if(!quiet) pb <- utils::txtProgressBar(min = 0, max = draws, style = 3)
  object_star <- object
  for (i in 1:draws){
    mb_obs[i] <- mb_chisq(object, psi[i,], p[i,])
    object_star@data@y <- ysim[,,i]
    #mb_chisq handles replicating NAs
    mb_sim[i] <- mb_chisq(object_star, psi[i,], p[i,])
    if(!quiet) utils::setTxtProgressBar(pb, i)
  }
  if(!quiet) close(pb)
  ubmsGOF("MacKenzie-Bailey Chi-square", data.frame(obs=mb_obs, sim=mb_sim))
})
