#' ubms
#'
#' Unmarked Bayesian Models using Stan
#'
#' @docType package
#'
#' @author Ken Kellner
#'
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom methods new
#' @importFrom unmarked getY siteCovs obsCovs
#' @importFrom stats as.formula model.frame predict rbinom
#' @useDynLib "ubms", .registration = TRUE
#' @name ubms
NULL
