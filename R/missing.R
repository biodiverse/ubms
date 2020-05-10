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
  reps <- get_row_reps(submodel, response)
  mat <- matrix(miss_vec, nrow=reps)
  submodel@missing <- apply(mat, 2, function(x) all(x))
  submodel
}

#' @include response.R
setMethod("update_missing", "ubmsResponse", function(object, submodels){
  object@missing <- find_missing(object, submodels)
  object
})

find_missing <- function(response, submodels){
  sub_mm <- lapply(submodels@submodels, expand_model_matrix, response)
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
