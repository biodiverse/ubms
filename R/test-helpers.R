expect_RMSE <- function(object, compare, maxRMSE){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
  stopifnot(length(object)==length(compare))
  act$RMSE <- sqrt(mean((compare-object)^2))
  testthat::expect(act$RMSE <= maxRMSE,
         sprintf("%s has RMSE %.3f, greater than max %.3f", act$lab,
                 act$RMSE, maxRMSE))
  invisible(act$val)
}

expect_between <- function(object, lower, upper){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")
  stopifnot(length(object)==1)
  testthat::expect(act$val <= upper & act$val >= lower,
         sprintf("%s is %.2f, outside the allowed range %.2f - %.2f",
                 act$lab, act$val, lower, upper))
  invisible(act$val)
}
