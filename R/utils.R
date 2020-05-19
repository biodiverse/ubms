get_samples <- function(fit, draws){
  nsamp <- nsamples(fit)
  out <- 1:nsamp
  if(!is.null(draws) && draws < nsamp){
    out <- sample(out, draws)
  }
  out
}

submodel_types <- function(object){
  names(object@submodels@submodels)
}

check_type <- function(type, possible_types){
  if(! type %in% possible_types){
    stop(paste("Submodel must be one of:", paste(possible_types, collapse=", ")),
         call.=FALSE)
  }
}
