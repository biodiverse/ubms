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
  ll <- get_loglik(refit, inps$stan_data)

  # Handle missing sites
  out <- matrix(NA, nrow=nrow(ll), ncol=sum(folds==i))
  not_missing <- which(folds==i) %in% which(!removed_sites(object))
  out[,not_missing] <- ll
  out
}

setGeneric("get_loglik", function(object, ...) standardGeneric("get_loglik"))

setMethod("get_loglik", "ubmsFit", function(object, inps, ...){
  stop("Method not available for this fit type", call.=FALSE)
})

setMethod("get_loglik", "ubmsFitOccu", function(object, inps, ...){
  psi <- t(posterior_linpred(object, submodel="state", transform=TRUE))
  p <- t(posterior_linpred(object, submodel="det", transform=TRUE))
  get_loglik_occu(inps$y, inps$M, inps$si-1, psi, p, inps$Kmin)
})

setMethod("get_loglik", "ubmsFitPcount", function(object, inps, ...){
  lam <- t(posterior_linpred(object, submodel="state", transform=TRUE))
  p <- t(posterior_linpred(object, submodel="det", transform=TRUE))
  get_loglik_pcount(inps$y, inps$M, inps$si-1, lam, p, inps$K, inps$Kmin)
})

setMethod("get_loglik", "ubmsFitOccuRN", function(object, inps, ...){
  lam <- t(posterior_linpred(object, submodel="state", transform=TRUE))
  p <- t(posterior_linpred(object, submodel="det", transform=TRUE))
  get_loglik_occuRN(inps$y, inps$M, inps$si-1, lam, p, inps$K, inps$Kmin)
})
