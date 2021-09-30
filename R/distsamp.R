#' Fit the Royle et al. (2004) Distance Sampling Model
#'
#' This function fits the hierarchical distance sampling model of Royle
#' et al. (2004) to line or point transect data recorded in discerete
#' distance intervals.
#'
#' @param formula Double right-hand side formula describing covariates of
#'  detection and occupancy in that order
#' @param data A \code{\link{unmarkedFrameDS}} object
#' @param keyfun One of the following detection functions:
#'  \code{"halfnorm"} for half-normal, \code{"exp"} for negative exponential,
#'  or \code{"hazard"} for hazard-rate (see warning below)
#' @param output Model either density \code{"density"} or abundance \code{"abund"}
#' @param unitsOut Units of density. Either \code{"ha"} or \code{"kmsq"} for
#'  hectares and square kilometers, respectively
#' @param prior_intercept_state Prior distribution for the intercept of the
#'  state (abundance) model; see \code{?priors} for options
#' @param prior_coef_state Prior distribution for the regression coefficients of
#'  the state model
#' @param prior_intercept_det Prior distribution for the intercept of the
#'  detection probability model
#' @param prior_coef_det Prior distribution for the regression coefficients of
#'  the detection model
#' @param prior_intercept_scale Prior distribution for the intercept of the
#'  scale parameter (i.e., log(scale)) for Hazard-rate models
#' @param prior_sigma Prior distribution on random effect standard deviations
#' @param ... Arguments passed to the \code{\link{stan}} call, such as
#'  number of chains \code{chains} or iterations \code{iter}
#'
#' @return \code{ubmsFitDistsamp} object describing the model fit.
#'
#' @note Values of `dist.breaks` in the `unmarkedFrameDS` should be as small
#'  as possible (<10) to facilitate convergence. Consider converting `unitsIn` from
#'  meters to kilometers, for example. See example below.
#'
#' @section Warning: Use of the hazard-rate key function (\code{"hazard"})
#'  typically requires a large sample size in order to get good parameter
#'  estimates. If you have a relatively small number of points/transects (<100),
#'  you should be cautious with the resulting models. Check your results against
#'  estimates from \code{unmarked}, which doesn't require as much data to get
#'  good estimates of the hazard-rate shape and scale parameters.
#'
#' @examples
#' \donttest{
#' data(issj)
#' #Note use of km instead of m for distance breaks
#' jayUMF <- unmarkedFrameDS(y=as.matrix(issj[,1:3]),
#'                           siteCovs=issj[,c("elevation","forest")],
#'                           dist.breaks=c(0,0.1,0.2,0.3),
#'                           unitsIn="km", survey="point")
#'
#' fm_jay <- stan_distsamp(~1~scale(elevation), jayUMF, chains=3, iter=300)
#' }
#'
#' @references Royle, J. A., Dawson, D. K., & Bates, S. (2004). Modeling
#'  abundance effects in distance sampling. Ecology 85: 1591-1597.
#'
#' @seealso \code{\link{distsamp}}, \code{\link{unmarkedFrameDS}}
#' @export
stan_distsamp <- function(formula,
                          data,
                          keyfun=c("halfnorm", "exp", "hazard"),
                          output=c("density", "abund"),
                          unitsOut=c("ha", "kmsq"),
                          prior_intercept_state = normal(0, 5),
                          prior_coef_state = normal(0, 2.5),
                          prior_intercept_det = normal(0, 5),
                          prior_coef_det = normal(0, 2.5),
                          prior_intercept_scale = normal(0,2.5),
                          prior_sigma = gamma(1, 1),
                          ...){

  forms <- split_formula(formula)
  umf <- process_umf(data)
  keyfun <- match.arg(keyfun)
  unitsOut <- match.arg(unitsOut)
  output <- match.arg(output)
  state_param <- switch(output, density={"Density"}, abund={"Abundance"})
  det_param <- switch(keyfun, halfnorm={"Scale"}, exp={"Rate"},
                      hazard={"Shape"})

  if(has_spatial(forms)){
    split_umf <- extract_missing_sites(umf)
    umf <- split_umf$umf
    state <- ubmsSubmodelSpatial(state_param, "state", siteCovs(umf), forms[[2]],
                                 "exp", prior_intercept_state, prior_coef_state,
                                 prior_sigma,
                                 split_umf$sites_augment, split_umf$data_aug)
  } else {
    state <- ubmsSubmodel(state_param, "state", siteCovs(umf), forms[[2]],
                          "exp", prior_intercept_state, prior_coef_state, prior_sigma)
  }

  response <- ubmsResponseDistsamp(data, keyfun, "P", output, unitsOut)
  det <- ubmsSubmodel(det_param, "det", siteCovs(umf), forms[[1]],
                      "exp", prior_intercept_det, prior_coef_det, prior_sigma)

  scale <- placeholderSubmodel("scale")
  if(keyfun=="hazard"){
    warning("Hazard key function may perform poorly with small sample sizes",
            call.=FALSE)
    scale <- ubmsSubmodelScalar("Scale", "scale", "exp", prior_intercept_scale)
  }

  submodels <- ubmsSubmodelList(state, det, scale)

  ubmsFit("distsamp", match.call(), data, response, submodels, ...)
}


#Output object-----------------------------------------------------------------

#' @include fit.R
setClass("ubmsFitDistsamp", contains = "ubmsFitAbun")


#Response class----------------------------------------------------------------

#' @include response.R
setClass("ubmsResponseDistsamp", contains="ubmsResponse",
         slots = c(survey="character", dist_breaks="numeric", tlength="numeric",
                   output="character", units_in="character", units_out="character"))

ubmsResponseDistsamp <- function(umf, y_dist, z_dist, output, units_out, K=NULL){
  max_primary <- ifelse(methods::.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  out <- ubmsResponse(umf@y, y_dist, z_dist, max_primary, K)
  out <- as(out, "ubmsResponseDistsamp")
  out@K <- get_K(out, K)
  out@survey <- umf@survey; out@output <- output
  out@dist_breaks <- umf@dist.breaks; out@tlength <- umf@tlength
  out@units_in <- umf@unitsIn; out@units_out <- units_out
  out
}

#NA handling not currently fully implemented; just used to send error
#' @include missing.R
setMethod("find_missing", "ubmsResponseDistsamp", function(object, submodels, ...){
  sub_mm <- submodels@submodels
  #Don't include transition parameters when finding missing values
  sub_mm <- sub_mm[!sapply(sub_mm, inherits, "ubmsSubmodelTransition")]
  #Don't include placeholder parmeters either
  sub_mm <- sub_mm[!sapply(sub_mm, is_placeholder)]
  sub_mm <- lapply(sub_mm, expand_model_matrix, object)
  comb <- cbind(as_vector(object), do.call("cbind", sub_mm))
  miss <- apply(comb, 1, function(x) any(is.na(x)))
  if(any(miss)){
    stop("Missing values are not allowed in y or covariates", call.=FALSE)
  }
  return(miss)
  #miss <- matrix(miss, nrow=nrow(t(object)))
  #remove_sites <- apply(miss, 2, function(x) any(x))
  #rep(remove_sites, each=object@max_obs)
})

#' @include inputs.R
setMethod("get_auxiliary_data", "ubmsResponseDistsamp", function(object, ...){
  out <- list(aux1=c(ifelse(object@survey=="point", 1, 0), 0), #Hack b/c stan won't handle 1-element vector
              aux2=object@dist_breaks,
              aux3=get_conv_const(object))
  c(out, list(n_aux1=2, n_aux2=length(out$aux2), n_aux3=length(out$aux3)))
})

#Builds the values in the denominator of the detection part of the likelihood
get_conv_const <- function(resp){
  db <- resp@dist_breaks
  M <- nrow(resp@y)

  if(resp@survey=="line"){
    w <- diff(db)
    a <- get_area(resp)
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 1 / matrix(rep(w, M), ncol=length(w), byrow=TRUE) * u
  } else {
    a <- get_area(resp)
    u <- t(apply(a, 1, function(x) x/sum(x)))
    out <- 2 * pi / a * u
  }
  out <- as.vector(t(out))
  area_adjust <- get_area_adjust(resp)
  out * rep(area_adjust, each=length(db)-1)
}

get_area <- function(resp){
  db <- resp@dist_breaks
  if(resp@survey=="line"){
    w <- diff(db)
    a <- t(sapply(resp@tlength, function(x) x*w))
  } else {
    a <- diff(pi*db^2)
    M <- nrow(resp@y)
    a <- matrix(rep(a, M), nrow=M, byrow=TRUE)
  }
  a
}

get_area_adjust <- function(resp){
  if(resp@output == "abund") return(rep(1, nrow(resp@y)))
  area <- get_area(resp)
  switch(resp@survey,
    line = A <- rowSums(area) * 2,
    point = A <- rowSums(area))
  switch(resp@units_in,
    m = A <- A / 1e6,
    km = A <- A)
  switch(resp@units_out,
    ha = A <- A * 100,
    kmsq = A <- A)
  A
}

#There is some duplicate code with the base class
#that could be cleaned up here maybe
setMethod("get_K", "ubmsResponseDistsamp", function(object, K=NULL){
  ysum <- rowSums(object@y, na.rm=TRUE)
  ymax <- max(ysum, na.rm=TRUE)
  if(is.null(K)) return(ymax + 50)
  if(K < ymax) stop("K must be larger than max y value", call.=FALSE)
  K
})

setMethod("get_Kmin", "ubmsResponseDistsamp", function(object){
  yt <- t(object)
  keep <- apply(yt, 2, function(x) !all(is.na(x)))

  yt <- matrix(yt, nrow=object@max_obs)
  out <- apply(yt, 2, function(x){
            if(all(is.na(x))) return(0)
            sum(x, na.rm=TRUE)
          })
  out <- matrix(out, ncol=object@max_primary, byrow=TRUE)
  out[keep,,drop=FALSE]
})


#Methods to simulate posterior predictive distributions------------------------

setMethod("sim_z", "ubmsFitDistsamp", function(object, samples, re.form, K=NULL, ...){
  resp <- object@response
  y <- resp@y

  lam_post <- t(sim_lp(object, "state", transform=TRUE, newdata=NULL,
                       re.form=re.form, samples=samples))
  if(resp@output == "density"){
    lam_post <- lam_post * get_area_adjust(resp)
  }
  p_post <- get_p_for_multinom(object, samples)

  K <- get_K(resp, K)
  Kmin <- get_Kmin(resp)
  kvals <- 0:K

  t(simz_multinom(y, lam_post, p_post, K, Kmin, kvals))
})

setMethod("sim_y", "ubmsFitDistsamp",
          function(object, samples, re.form, z=NULL, K=NULL, ...){
  nsamp <- length(samples)
  M <- get_n_sites(object@response)
  J <- object@response@max_obs
  z <- process_z(object, samples, re.form, z)
  p <- get_p_for_multinom(object, samples)

  out <- array(NA, c(J,M, nsamp))
  for (i in 1:nsamp){
    for (m in 1:M){
      out[,m,i] <- stats::rmultinom(n=1, size=z[m,i], prob=p[m,,i])[1:J]
    }
  }
  matrix(out, nrow=nsamp, byrow=TRUE)
})

get_p_for_multinom <- function(object, samples){
  M <- get_n_sites(object@response)
  J <- object@response@max_obs
  nsamp <- length(samples)
  p_post_raw <- t(sim_p(object, samples=samples))
  p_post_raw <- array(p_post_raw, c(J, M, nsamp))
  p_post_raw <- aperm(p_post_raw, c(2,1,3))

  p_post <- array(NA, c(M, J+1, nsamp))
  for (i in 1:nsamp){
    p_post[,1:J,i] <- p_post_raw[,1:J,i]
    p_post[,J+1,i] <- 1 - rowSums(p_post_raw[,1:J,i], na.rm=TRUE)
  }
  p_post
}

#Get detection probability-----------------------------------------------------

#' @include posterior_linpred.R
setMethod("sim_p", "ubmsFitDistsamp", function(object, samples, ...){
  resp <- object@response
  resp@output <- "abund" #Don't adjust for area
  db <- resp@dist_breaks
  conv_const <- get_conv_const(resp)
  inds <- get_subset_inds(resp)[,1:2]
  param1 <- sim_lp(object, "det", transform=TRUE, newdata=NULL,
                          samples=samples, re.form=NULL)

  stopifnot(ncol(param1) == get_n_sites(resp))
  param2 <- NULL
  if(resp@y_dist=="halfnorm"){
    pfun <- ifelse(resp@survey=="line", distprob_normal_line, distprob_normal_point)
  } else if(resp@y_dist == "exp"){
    pfun <- ifelse(resp@survey=="line", distprob_exp_line, distprob_exp_point)
  } else if(resp@y_dist == "hazard"){
    param2 <- sim_lp(object, "scale", transform=TRUE, newdata=NULL,
                     samples=samples, re.form=NULL)
    pfun <- ifelse(resp@survey=="line", distprob_haz_line, distprob_haz_point)
  }
  out <- sapply(1:length(samples), function(i){
    pfun(param1[i,], param2[i,1], db, conv_const, inds)
  })
  t(out)
})

distprob_normal_line <- function(sigma, param2, db, conv_const, inds){
  out <- sapply(1:length(sigma), function(i){
    int <-  stats::pnorm(db[-1], 0, sd=sigma[i]) -
            stats::pnorm(db[-length(db)], 0, sd=sigma[i])
    int <- int / stats::dnorm(0, 0, sd=sigma[i])
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

distprob_normal_point <- function(sigma, param2, db, conv_const, inds){
  out <- sapply(1:length(sigma), function(i){
    s2 <- sigma[i]^2
    a <- db[-1]; b <- db[-length(db)]
    int <- s2 * ((1 - exp(-a*a / (2*s2))) - (1-exp(-b*b / (2*s2))))
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

distprob_exp_line <- function(rate, param2, db, conv_const, inds){
  out <- sapply(1:length(rate), function(i){
    a <- db[-length(db)]; b <- db[-1];
    int <- rate[i] * (exp(-a/rate[i]) - exp(-b/rate[i]))
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

distprob_exp_point <- function(rate, param2, db, conv_const, inds){
  out <- sapply(1:length(rate), function(i){
    a <- db[-length(db)]; b <- db[-1]
    int <- rate[i] * exp(-a/rate[i]) * (a+rate[i]) -
           rate[i] * exp(-b/rate[i]) * (b+rate[i])
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

distprob_haz_line <- function(shape, scale, db, conv_const, inds){
  out <- sapply(1:length(shape), function(i){
    a <- db[-length(db)]; b <- db[-1];
    int <- numeric(length(a))
    for (j in 1:length(int)){
      int[j] <- stats::integrate(unmarked::gxhaz, a[j], b[j], shape=shape[i], scale=scale[1])$value
    }
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}

distprob_haz_point <- function(shape, scale, db, conv_const, inds){
  out <- sapply(1:length(shape), function(i){
    a <- db[-length(db)]; b <- db[-1]
    int <- numeric(length(a))
    for (j in 1:length(int)){
      int[j] <- stats::integrate(unmarked::grhaz, a[j], b[j],
                                 shape=shape[i], scale=scale)$value
    }
    int * conv_const[inds[i,1]:inds[i,2]]
  })
  as.vector(out)
}


#Goodness-of-fit---------------------------------------------------------------

#' @include posterior_linpred.R
setMethod("sim_state", "ubmsFitDistsamp", function(object, samples, ...){
    lp <- methods::callNextMethod(object, samples, ...)
    if(object@response@output == "density"){
      lp <- t(t(lp) * get_area_adjust(object@response))
    }
    lp
})


#Method for fitted values------------------------------------------------------

setMethod("sim_fitted", "ubmsFitDistsamp", function(object, submodel, samples, ...){
  stopifnot(submodel %in% c("state", "det"))
  if(submodel == "state"){
    lp <- sim_lp(object, submodel, transform=TRUE, newdata=NULL, samples=samples,
                        re.form=NULL)
    if(object@response@output == "density"){
      lp <- t(t(lp) * get_area_adjust(object@response))
    }
    return(lp)
  }

  p <- sim_p(object, samples=samples)
  J <- object@response@max_obs
  z <- sim_z(object, samples, re.form=NULL)
  z <- z[, rep(1:ncol(z), each=J)]
  out <- z * p
  #out[z == 0] <- NA
  out
})


#Histogram---------------------------------------------------------------------

#' @importFrom graphics hist
#' @importFrom ggplot2 geom_histogram
setMethod("hist", "ubmsFitDistsamp", function(x, draws=30, ...){
  samples <- get_samples(x, draws)
  hist_data <- get_hist_data(x)
  mean_line <- get_mean_line(x)

  out <- ggplot(hist_data, aes_string(x="x")) +
    geom_histogram(aes_string(y="..density.."),fill='transparent',
                   col='black',breaks=x@response@dist_breaks)

  #Adjust the histogram height to match the density line
  bar_height <- ggplot2::ggplot_build(out)$data[[1]]$y[1]
  adj_factor <- max(mean_line$val, na.rm=TRUE) / bar_height

  out <- ggplot(hist_data, aes_string(x="x")) +
    geom_histogram(aes_string(y=paste0("..density..*",adj_factor)),fill='transparent',
                   col='black',breaks=x@response@dist_breaks)

  if(draws > 0){
    sample_lines <- get_sample_lines(x, samples)
    out <- out +
      geom_line(data=sample_lines, aes_string(x="x", y="val", group="ind"),
              alpha=0.3)
  }
  out +
    geom_line(data=mean_line, aes_string(x="x", y="val"), col='red') +
    labs(x=paste0("Distance (", x@response@units_in,")"), y="Density") +
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text=element_text(size=12), axis.title=element_text(size=14))
})

get_hist_data <- function(fit){
  resp <- fit@response
  db <- fit@response@dist_breaks
  bin <- db[-length(db)] + diff(db)/2
  counts <- colSums(resp@y, na.rm=TRUE)
  data.frame(x=rep(bin, times=counts))
}

get_mean_line <- function(fit){
  db <- fit@response@dist_breaks
  xseq <- seq(db[1], db[length(db)], length.out=1000)
  par1 <- summary(fit, "det")
  if(nrow(par1) > 1) warning("Ignoring covariate effects", call.=FALSE)
  par1 <- exp(par1[1,1])
  detfun <- get_detfun(fit)
  if(fit@response@y_dist == "hazard"){
    par2 <- exp(summary(fit, "scale")[1,1])
    val <- detfun(xseq, par1, par2)
  } else {
    val <- detfun(xseq, par1)
  }
  data.frame(x=xseq, val=val)
}

get_sample_lines <- function(fit, samples){
  db <- fit@response@dist_breaks
  xseq <- seq(db[1], db[length(db)], length.out=1000)
  par1 <- exp(extract(fit, "beta_det")[[1]][,1])
  detfun <- get_detfun(fit)
  df_try <- function(...){
    tryCatch(detfun(...), error=function(e) NA)
  }
  if(fit@response@y_dist=="hazard"){
    par2 <- exp(extract(fit, "beta_scale")[[1]])
    sample_lines <- lapply(samples, function(i){
                         data.frame(x=xseq, val=df_try(xseq, par1[i], par2[i]),ind=i)
                      })
  } else {
    sample_lines <- lapply(samples, function(i){
                         data.frame(x=xseq, val=df_try(xseq, par1[i]),ind=i)
                      })
  }
  out <- do.call("rbind", sample_lines)
  out[!is.na(out$val),]
}

get_detfun <- function(fit){
  resp <- fit@response
  if(resp@y_dist == "halfnorm"){
    out <- ifelse(resp@survey=="line", unmarked::dxhn, unmarked::drhn)
  } else if(resp@y_dist == "exp"){
    out <- ifelse(resp@survey=="line", unmarked::dxexp, unmarked::drexp)
  } else if(resp@y_dist == "hazard"){
    out <- ifelse(resp@survey=="line", unmarked::dxhaz, unmarked::drhaz)
  }
  out
}
