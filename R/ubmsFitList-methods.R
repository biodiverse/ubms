#' Extractors for ubmsFitList objects

#' Extract parts of ubmsFitList objects.
#'
#' @name ubmsFitList-extractors
#'
#' @param x A list of \code{ubmsFit} models of class \code{ubmsFitList}
#' @param name,i The names or indices of \code{ubmsFit} models in the \code{ubmsFitList}
#'
#' @return A \code{ubmsFit} object or list of such objects.
#'
NULL

#' @rdname ubmsFitList-extractors
#' @export
setMethod("$", "ubmsFitList", function(x, name){
  x@models[[name]]
})

#' @rdname ubmsFitList-extractors
#' @export
setMethod("[[", c("ubmsFitList", "numeric","missing"),
          function(x, i) x@models[[i]])

#' @rdname ubmsFitList-extractors
#' @export
setMethod("[", c("ubmsFitList", "numeric","missing","missing"),
          function(x, i) x@models[i])

#' Get Names of Models in a ubmsFitList
#'
#' @param x A \code{ubmsFitList} object
#'
#' @return A character vector of model names.
#'
#' @export
setMethod("names", "ubmsFitList", function(x){
  names(x@models)
})
