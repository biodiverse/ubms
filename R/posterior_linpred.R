#' @importFrom rstantools posterior_linpred
#' @include fit.R
#' @export
setMethod("posterior_linpred", "ubmsFit", 
          function(object, transform=FALSE, submodel, newdata=NULL, draws=NULL, 
                   re.form=NULL, ...){
 
  check_type(submodel, submodel_types(object))
  samp_inds <- get_samples(object, draws)

  sim_lp(object, submodel=submodel, transform=transform, newdata=newdata, 
         samples=samp_inds, re.form=re.form) 
})


setGeneric("sim_lp", function(object, ...) standardGeneric("sim_lp"))

setMethod("sim_lp", "ubmsFit", function(object, submodel, transform, newdata, 
                                        samples, re.form, ...){  
  sm <- object[submodel]
  beta <- extract(object, beta_par(sm))[[1]]
  lp <- model.matrix(sm, newdata) %*% t(beta[samples,,drop=FALSE])
 
  if(has_random(sm) & is.null(re.form)){
    b <- extract(object, b_par(sm))[[1]]
    lp <- lp + Z_matrix(sm, newdata) %*% t(b[samples,,drop=FALSE])
  }
  if(transform) lp <- do.call(sm@link, list(lp))
  t(unname(lp))
})
