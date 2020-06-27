#' Methods for ubmsFitList objects

#' Basic methods for ubmsFitList objects. More complex methods have
#' their own documentation pages.
#'
#' @name ubmsFitList-methods
#'
#' @param x A list of \code{ubmsFit} models of class \code{ubmsFitList}
#' @param name,i A single \code{ubmsFit} model
#'
NULL

#' @rdname ubmsFitList-methods
#' @export
setMethod("$", "ubmsFitList", function(x, name){
  x@models[[name]]
})

#' @rdname ubmsFitList-methods
#' @export
setMethod("names", "ubmsFitList", function(x){
  names(x@models)
})

#' @rdname ubmsFitList-methods
#' @export
setMethod("[[", c("ubmsFitList", "numeric","missing"),
          function(x, i) x@models[[i]])

#' @rdname ubmsFitList-methods
#' @export
setMethod("[", c("ubmsFitList", "numeric","missing","missing"),
          function(x, i) x@models[i])
