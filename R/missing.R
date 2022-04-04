setGeneric("update_missing", function(object, object2, ...)
           standardGeneric("update_missing"))

#' @include submodel.R

#setMethod("update_missing", "ubmsSubmodelList", function(object, response){
#  is_na <- find_missing(response, object)
#  mods <- lapply(object@submodels, submodel_set_missing, is_na, response)
#  object@submodels <- mods
#  object
#})

setMethod("update_missing", c("ubmsSubmodelList", "ubmsResponse"),
  function(object, object2){
  response <- object2
  is_na <- find_missing(response, object)
  mods <- lapply(object@submodels, submodel_set_missing, is_na, response)
  object@submodels <- mods
  object
})

submodel_set_missing <- function(submodel, miss_vec, response){
  if(is_placeholder(submodel)) return(submodel)
  new_miss_vec <- new_missing_vec(submodel, miss_vec, response)
  stopifnot(length(submodel@missing) == length(new_miss_vec))
  submodel@missing <- new_miss_vec
  submodel
}

setGeneric("new_missing_vec", function(object, miss_vec, response, ...){
  standardGeneric("new_missing_vec")
})

setMethod("new_missing_vec", "ubmsSubmodel",
          function(object, miss_vec, response, ...){
  reps <- get_row_reps(object, response)
  mat <- matrix(miss_vec, nrow=reps)
  apply(mat, 2, function(x) all(x))
})

#For transition parameters, only drop values for sites with no detections
setMethod("new_missing_vec", "ubmsSubmodelTransition",
          function(object, miss_vec, response, ...){
  drop_sites <- apply(t(response), 2, function(x) all(is.na(x)))
  rep(drop_sites, each=(response@max_primary-1))
})

#' @include response.R
setMethod("update_missing", c("ubmsResponse", "ubmsSubmodelList"),
  function(object, object2){

  submodels <- object2
  object@missing <- find_missing(object, submodels)
  object
})

setGeneric("find_missing", function(object, ...) standardGeneric("find_missing"))

setMethod("find_missing", "ubmsResponse", function(object, submodels, ...){
  sub_mm <- submodels@submodels
  #Don't include transition parameters when finding missing values
  sub_mm <- sub_mm[!sapply(sub_mm, inherits, "ubmsSubmodelTransition")]
  #Don't include placeholder parmeters either
  sub_mm <- sub_mm[!sapply(sub_mm, is_placeholder)]
  sub_mm <- lapply(sub_mm, expand_model_matrix, object)
  comb <- cbind(as_vector(object), do.call("cbind", sub_mm))
  apply(comb, 1, function(x) any(is.na(x)))
})

#find_missing <- function(response, submodels){
  #sub_mm <- submodels@submodels
  #Don't include transition parameters when finding missing values
  #sub_mm <- sub_mm[!sapply(sub_mm, inherits, "ubmsSubmodelTransition")]
  #sub_mm <- lapply(sub_mm, expand_model_matrix, response)
  #comb <- cbind(as_vector(response), do.call("cbind", sub_mm))
  #apply(comb, 1, function(x) any(is.na(x)))
#}

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

#' @include fit.R
removed_sites <- function(object){
  object['state']@missing
}
