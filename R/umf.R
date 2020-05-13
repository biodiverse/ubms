setGeneric("process_umf", function(umf, ...) standardGeneric("process_umf"))

setMethod("process_umf", "unmarkedFrame", function(umf){

  nsites <- unmarked::numSites(umf)
  nprimary <- ifelse(.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  nobs <- ncol(umf@y) / nprimary  
  
  siteCovs(umf) <- process_covs(siteCovs(umf), nsites)

  ysc <- siteCovs(umf)
  if(nprimary > 1){
    ysc <- process_covs(yearlySiteCovs(umf), nsites*nprimary, siteCovs(umf))
    yearlySiteCovs(umf) <- ysc
  }

  obsCovs(umf) <- process_covs(obsCovs(umf), nsites*nprimary*nobs, ysc)

  umf
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

