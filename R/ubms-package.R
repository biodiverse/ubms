#' ubms
#'
#' Unmarked Bayesian Models using Stan
#'
#' @author Ken Kellner
#'
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom methods new
#' @importFrom unmarked getY siteCovs obsCovs
#' @importFrom stats as.formula model.frame rbinom dbinom
#' @importFrom rlang .data
#' @useDynLib "ubms", .registration = TRUE
#' @name ubms
#' @aliases ubms-package
NULL
