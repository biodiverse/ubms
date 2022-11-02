# Calculate parameter from X, Z, offset
# Can't use posterior_linpred for this because it doesn't drop NAs
# Inps are output from get_stan_inputs()
# Right now the permute argument isn't used
calculate_par <- function(object, inps, submodel, permute=FALSE){
  X <- inps$stan_data[[paste0("X_",submodel)]]
  beta <- extract_posterior(object, paste0("beta_", submodel))
  lp <- X %*% t(beta)
  if(inps$stan_data[[paste0("has_random_",submodel)]]){
    Z <- Z_matrix(inps$submodels[submodel], na.rm=TRUE)
    b <- extract_posterior(object, paste0("b_", submodel))
    lp <- lp + Z %*% t(b)
  }
  lp <- lp + inps$stan_data[[paste0("offset_",submodel)]]
  do.call(inps$submodels[submodel]@link, list(lp))
}

extract_posterior <- function(object, par, permute=FALSE){
  if(permute) return(extract(object@stanfit, par)[[1]])
  ex <- extract(object@stanfit, par, permuted=FALSE)
  smp <- lapply(1:ncol(ex), function(x){
                  sapply(1:dim(ex)[3], function(i) ex[,x,i,drop=FALSE])
                })
  do.call("rbind", smp)
}

setGeneric("get_loglik", function(object, ...) standardGeneric("get_loglik"))

#' @include fit.R
setMethod("get_loglik", "ubmsFit", function(object, inps, ...){
  stop("Method not available for this fit type", call.=FALSE)
})

#' @include occu.R
setMethod("get_loglik", "ubmsFitOccu", function(object, inps, ...){
  psi <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_occu(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                  psi, p, inps$stan_data$Kmin)
})

#' @include pcount.R
setMethod("get_loglik", "ubmsFitPcount", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_pcount(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                    lam, p, inps$stan_data$K, inps$stan_data$Kmin)
})

#' @include occuRN.R
setMethod("get_loglik", "ubmsFitOccuRN", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_occuRN(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                    lam, p, inps$stan_data$K, inps$stan_data$Kmin)
})

#' @include multinomPois.R
setMethod("get_loglik", "ubmsFitMultinomPois", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_multinomPois(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                          lam, p, inps$stan_data$y_dist)
})

#' @include colext.R
setMethod("get_loglik", "ubmsFitColext", function(object, inps, ...){
  psi <- calculate_par(object, inps, submodel="state")
  psicube <- array(NA, c(dim(psi), 2))
  psicube[,,1] <- 1 - psi
  psicube[,,2] <- psi
  psicube <- aperm(psicube, c(1,3,2))

  col <- calculate_par(object, inps, submodel="col")
  ext <- calculate_par(object, inps, submodel="ext")
  phicube <- array(NA, c(dim(col), 4))
  phicube[,,1] <- 1 - col
  phicube[,,2] <- col
  phicube[,,3] <- ext
  phicube[,,4] <- 1 - ext
  phicube <- aperm(phicube, c(1,3,2))
  pmat <- calculate_par(object, inps, submodel="det")

  get_loglik_colext(inps$stan_data$y, inps$stan_data$M, inps$stan_data$Tsamp-1,
                     inps$stan_data$J, inps$stan_data$si-1, psicube,
                     phicube, pmat, 1 - inps$stan_data$Kmin)
})

#' @include occuTTD.R
setMethod("get_loglik", "ubmsFitOccuTTD", function(object, inps, ...){
  psi <- calculate_par(object, inps, submodel="state")
  lam <- calculate_par(object, inps, submodel="det")
  k <- rep(0, ncol(psi))
  if(inps$stan_data$y_dist == 3){
    k <- extract(object, "beta_shape", permute=FALSE)
    k <- exp(as.vector(k))
  }
  get_loglik_occuTTD(inps$stan_data$aux2, inps$stan_data$M, inps$stan_data$si-1,
                     psi, lam, k, inps$stan_data$aux1, inps$stan_data$y_dist)
})

#' Extract Pointwise Log-likelihood From Model
#'
#' Extracts the pointwise log-likelihood matrix or array from a model.
#' This is useful as an input for functions in the \code{loo} package.
#' If called on a \code{ubmsFit} object, the log-likelihood matrix or array
#' is calculated using the posterior distributions of the parameters; the
#' log-likelihood values are not actually saved inside the model object.
#' If called on a \code{stanfit} object, \code{loo::extract_log_lik} is used.
#' In this case, the log-likelihood must be saved as one of the output
#' parameters in Stan.
#'
#' @param object A \code{ubmsFit} or \code{stanfit} object
#' @param parameter_name The name of the log-likelihood parameter in the
#'  Stan output; ignored when \code{object} is a \code{ubmsFit}
#' @param merge_chains If \code{TRUE} (the default), all Markov chains are
#'  merged together (i.e., stacked) and a matrix is returned. If ‘FALSE’
#'  they are kept separate and an array is returned.
#'
#' @return A matrix (samples x sites) or array (samples x chains x sites)
#'
#' @aliases extract_log_lik,ubmsFit-method extract_log_lik,ubmsFitDistsamp-method
#' @export
setGeneric("extract_log_lik",
           function(object, parameter_name="log_lik", merge_chains=TRUE){
             loo::extract_log_lik(object, parameter_name, merge_chains)
           })

setMethod("extract_log_lik", "ubmsFit",
          function(object, parameter_name = "log_lik", merge_chains=TRUE){
  #if log_lik was saved as an output parameter, just return that
  if(loglik_saved(object)){
    return(loo::extract_log_lik(object@stanfit, parameter_name, merge_chains))
  }
  inps <- rebuild_inputs(object)
  ll <- get_loglik(object, inps, merge_chains)
  if(!merge_chains){
    nchain <- length(rstan::get_inits(object@stanfit))
    niter <- nrow(ll) / nchain
    ll <- array(ll, c(niter, nchain, ncol(ll)))
  }
  ll
})

loglik_saved <- function(fit){
  "log_lik" %in% fit@stanfit@sim$pars_oi
}
