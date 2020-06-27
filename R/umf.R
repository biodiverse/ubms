setGeneric("process_umf", function(umf, ...) standardGeneric("process_umf"))

#' @importFrom unmarked yearlySiteCovs
#' @importFrom unmarked "siteCovs<-" "obsCovs<-" "yearlySiteCovs<-"
setMethod("process_umf", "unmarkedFrame", function(umf){

  nsites <- unmarked::numSites(umf)
  nprimary <- ifelse(methods::.hasSlot(umf, "numPrimary"), umf@numPrimary, 1)
  nobs <- ncol(umf@y) / nprimary

  siteCovs(umf) <- process_covs(siteCovs(umf), nsites)

  ysc <- siteCovs(umf)
  if(nprimary > 1){
    ysc <- process_covs(yearlySiteCovs(umf), nsites*nprimary, siteCovs(umf))
  }
  obsCovs(umf) <- process_covs(obsCovs(umf), nsites*nprimary*nobs, ysc)

  if(nprimary > 1){
    yearlySiteCovs(umf) <- drop_final_year(ysc, nsites, nprimary)
  }

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

drop_final_year <- function(yr_df, nsites, nprimary){
  to_drop <- seq(nprimary, nsites*nprimary, by=nprimary)
  yr_df <- yr_df[-to_drop,,drop=FALSE]
  droplevels(yr_df)
}
