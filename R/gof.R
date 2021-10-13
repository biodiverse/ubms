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
    geom_label(data=ppval, aes_string(x=-Inf, y=Inf, label="lab"),
               hjust=-0.2, vjust=1.4, size=5,
               fill='white', label.size=0,
               label.padding=unit(0.1, "lines")) +
    ggtitle(paste("Posterior predictive check:", x@statistic)) +
    labs(y="Simulated data", x="Observed data") +
    plot_theme() +
    theme(strip.background=element_rect(fill="transparent"),
          strip.text=element_text(size=14))
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
#' @importFrom pbapply pboptions pblapply
#' @export
setGeneric("gof", function(object, draws=NULL, ...){
             standardGeneric("gof")})


# Generic function for goodness-of-fit test------------------------------------
setGeneric("sim_gof", function(object, ...) standardGeneric("sim_gof"))

setMethod("sim_gof", "ubmsFit", function(object, draws, func, name, quiet=FALSE, ...){
  samples <- get_samples(object, draws)
  draws <- length(samples)

  state <- sim_state(object, samples)
  p <- sim_p(object, samples)

  #M <- ncol(t(object@response)) #Should fix get_n_sites method for this
  T <- object@response@max_primary
  M <- ncol(state) / T
  R <- T * object@response@max_obs
  ysim <- suppressMessages(sim_y(object, samples, re.form=NULL))
  stopifnot(ncol(ysim) == M*R)
  ysim <- array(ysim, c(draws,R,M))
  ysim <- aperm(ysim, c(3,2,1))

  object_star <- object
  op <- pbapply::pboptions()
  if(quiet) pbapply::pboptions(type = "none")
  out <- pbapply::pblapply(1:draws, function(i){
    stat_obs <- func(object, state[i,], p[i,])
    object_star@response <- ubmsResponse(ysim[,,i], object@response@y_dist,
                                         object@response@z_dist, max_primary=T)
    stat_sim <- func(object_star, state[i,], p[i,])
    c(stat_obs, stat_sim)
  })
  pbapply::pboptions(op)
  out <- do.call(rbind, out)

  ubmsGOF(name, data.frame(obs=out[,1], sim=out[,2]))
})

# N-mixture Chi-square test for abundance models-------------------------------

Nmix_chisq <- function(object, lambda, p){
  J <- object@response@max_obs
  lambda <- lambda[rep(1:length(lambda), each=J)]
  stopifnot(length(lambda) == length(p))
  fit <- lambda*p
  obs <- as_vector(object@response, na.rm=FALSE)
  stopifnot(length(obs) == length(fit))
  sum((obs - fit)^2/fit, na.rm=TRUE)
}
