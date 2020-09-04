setGeneric("process_umf", function(umf, ...) standardGeneric("process_umf"))

#' @importFrom unmarked yearlySiteCovs
#' @importFrom unmarked "siteCovs<-" "obsCovs<-" "yearlySiteCovs<-"
setMethod("process_umf", "unmarkedFrame", function(umf){

  nsites <- unmarked::numSites(umf)
  nprimary <- ifelse(methods::.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  nobs <- num_obs(umf)

  siteCovs(umf) <- process_covs(siteCovs(umf), nsites)

  ysc <- siteCovs(umf)
  if(nprimary > 1){
    ysc <- process_covs(yearlySiteCovs(umf), nsites*nprimary, siteCovs(umf))
    yearlySiteCovs(umf) <- ysc
  }

  if(!inherits(umf, "unmarkedFrameDS")){
    obsCovs(umf) <- process_covs(obsCovs(umf), nsites*nprimary*nobs, ysc)
  }

  umf
})

setGeneric("num_obs", function(umf, ...) standardGeneric("num_obs"))

setMethod("num_obs", "unmarkedFrame", function(umf){
  nprimary <- ifelse(methods::.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  ncol(umf@y) / nprimary
})

setMethod("num_obs", "unmarkedFrameMPois", function(umf){
  nrow(umf@obsToY)
})

process_covs <- function(covs, n_row, inherit=NULL){

  if(is.null(covs)) covs <- data.frame(.dummy=rep(NA, n_row))

  if(!is.null(inherit) && !identical(names(inherit), ".dummy")){
    rep_rows <- n_row / nrow(inherit)
    expand_df <- inherit[rep(1:nrow(inherit), each=rep_rows),,drop=FALSE]
    covs <- cbind(covs, expand_df)
  }
  covs
}
