#' K-fold Cross-validation of a ubmsFit Model
#'
#' @param x A \code{ubmsFit} model
#' @param K Number of folds into which the data will be partitioned
#' @param folds An optional vector with length equal to the number of sites in the data and containing integers from 1 to K, to manually assign sites to folds. You should use this if you plan to compare multiple models, since the folds for each model should be identical. You can use \code{loo::kfold_split_random} to generate this vector
#' @param quiet If \code{TRUE}, suppress progress bar
#' @param ... Currently ignored
#'
#' @return An object of class \code{elpd_generic} that is compatible with \code{loo::loo_compare}
#'
#' @include fit.R
#' @importFrom loo kfold
#' @export
setMethod("kfold", "ubmsFit", function(x, K=10, folds=NULL, quiet=FALSE, ...){
  if(is.null(folds)){
    folds <- loo::kfold_split_random(K, numSites(x@data))
  } else {
    stopifnot(length(folds) == numSites(x@data))
    stopifnot(max(folds) == K)
  }

  op <- pbapply::pboptions()
  if(quiet) pbapply::pboptions(type = "none")
  ll_unsort <- pbapply::pblapply(1:K, ll_fold, x, folds)
  pbapply::pboptions(op)
  ll_unsort <- do.call("cbind", ll_unsort)

  # Sort sites back to original order
  sort_ind <- unlist(lapply(1:K, function(i) which(folds==i)))
  ll <- matrix(NA, nrow=nrow(ll_unsort), ncol=ncol(ll_unsort))
  ll[,sort_ind] <- ll_unsort
  ll <- ll[,apply(ll, 2, function(x) !any(is.na(x))),drop=FALSE]

  loo::elpd(ll)

})

ll_fold <- function(i, object, folds){
  cl <- object@call
  train_data <- object@data[which(!folds==i),]
  cl$data <- train_data
  cl$refresh <- 0
  refit <- suppressWarnings(eval(cl))

  test_data <- object@data[which(folds==i),]
  cl$data <- test_data
  cl$return_inputs <- TRUE
  inps <- eval(cl)
  refit@submodels <- inps$submodels
  ll <- get_loglik(refit, inps)

  # Handle missing sites
  out <- matrix(NA, nrow=nrow(ll), ncol=sum(folds==i))
  not_missing <- which(folds==i) %in% which(!removed_sites(object))
  out[,not_missing] <- ll
  out
}

# Calculate parameter from X, Z, offset
# Can't use posterior_linpred for this because it doesn't drop NAs
# Inps are output from get_stan_inputs()
calculate_par <- function(object, inps, submodel){
  X <- inps$stan_data[[paste0("X_",submodel)]]
  beta <- extract(object, paste0("beta_", submodel))[[1]]
  lp <- X %*% t(beta)
  if(inps$stan_data[[paste0("has_random_",submodel)]]){
    Z <- Z_matrix(inps$submodels[submodel])
    b <- extract(object, paste0("b_", submodel))[[1]]
    lp <- lp + Z %*% t(b)
  }
  lp <- lp + inps$stan_data[[paste0("offset_",submodel)]]
  do.call(inps$submodels[submodel]@link, list(lp))
}

setGeneric("get_loglik", function(object, ...) standardGeneric("get_loglik"))

setMethod("get_loglik", "ubmsFit", function(object, inps, ...){
  stop("Method not available for this fit type", call.=FALSE)
})

setMethod("get_loglik", "ubmsFitOccu", function(object, inps, ...){
  psi <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_occu(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                  psi, p, inps$stan_data$Kmin)
})

setMethod("get_loglik", "ubmsFitPcount", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_pcount(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                    lam, p, inps$stan_data$K, inps$stan_data$Kmin)
})

setMethod("get_loglik", "ubmsFitOccuRN", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_occuRN(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                    lam, p, inps$stan_data$K, inps$stan_data$Kmin)
})

setMethod("get_loglik", "ubmsFitMultinomPois", function(object, inps, ...){
  lam <- calculate_par(object, inps, submodel="state")
  p <- calculate_par(object, inps, submodel="det")
  get_loglik_multinomPois(inps$stan_data$y, inps$stan_data$M, inps$stan_data$si-1,
                          lam, p, inps$stan_data$y_dist)
})

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
