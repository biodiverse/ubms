#' @importFrom rstantools posterior_linpred
#' @include fit.R
#' @export
setMethod("posterior_linpred", "ubmsFit", 
          function(object, transform=FALSE, type, newdata=NULL, draws=NULL, 
                   re.form=NULL, ...){
 
  check_type(type, submodel_types(object))
  samp_inds <- get_samples(object, draws)

  sim_lp(object, type=type, transform=transform, newdata=newdata, 
         samples=samp_inds, re.form=re.form) 
})


setGeneric("sim_lp", function(object, ...) standardGeneric("sim_lp"))

setMethod("sim_lp", "ubmsFit", function(object, type, transform, newdata, 
                                        samples, re.form, ...){  
  sm <- object[type]
  beta <- extract(object, beta_par(sm))[[1]]
  lp <- model.matrix(sm, newdata) %*% t(beta[samples,,drop=FALSE])
 
  if(has_random(sm) & is.null(re.form)){
    b <- extract(object, b_par(sm))[[1]]
    lp <- lp + Z_matrix(sm, newdata) %*% t(b[samples,,drop=FALSE])
  }
  if(transform) lp <- do.call(sm@link, list(lp))
  t(unname(lp))
})


check_type <- function(type, possible_types){
  if(! type %in% possible_types){
    stop(paste("Type must be one of:", paste(possible_types, collapse=", ")),
         call.=FALSE)
  }
}
