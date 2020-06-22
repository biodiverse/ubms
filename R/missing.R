setGeneric("update_missing", function(object, ...) 
           standardGeneric("update_missing"))

#' @include submodel.R
setMethod("update_missing", "ubmsSubmodelList", function(object, response){
  is_na <- find_missing(response, object)
  mods <- lapply(object@submodels, submodel_set_missing, is_na, response)
  object@submodels <- mods
  object
})

submodel_set_missing <- function(submodel, miss_vec, response){
  #For transition parameters, only drop values for sites with no detections
  if(submodel@transition){
    yt <- t(response)
    drop_sites <- apply(yt, 2, function(x) all(is.na(x)))
    T <- response@max_primary
    new_miss_vec <- rep(drop_sites, each=(T-1))
  } else{
    reps <- get_row_reps(submodel, response)
    mat <- matrix(miss_vec, nrow=reps)
    new_miss_vec <- apply(mat, 2, function(x) all(x))
  }
  stopifnot(length(submodel@missing) == length(new_miss_vec))
  submodel@missing <- new_miss_vec
  submodel
}

#' @include response.R
setMethod("update_missing", "ubmsResponse", function(object, submodels){
  object@missing <- find_missing(object, submodels)
  object
})

find_missing <- function(response, submodels){
  sub_mm <- submodels@submodels
  #Don't include transition parameters when finding missing values
  sub_mm <- sub_mm[sapply(sub_mm, function(x) !x@transition)]
  sub_mm <- lapply(sub_mm, expand_model_matrix, response)
  comb <- cbind(as_vector(response), do.call("cbind", sub_mm))
  apply(comb, 1, function(x) any(is.na(x)))
}

expand_model_matrix <- function(submodel, response){
  row_reps <- get_row_reps(submodel, response)
  mm <- model.matrix(submodel)
  mm[rep(1:nrow(mm), each=row_reps),,drop=FALSE]
}

get_row_reps <- function(submodel, response){
  yl <- length(as_vector(response))
  n <- nrow(model.matrix(submodel))
  if(yl %% n != 0) stop("Invalid model matrix dimensions", call.=FALSE)
  yl / n
}
