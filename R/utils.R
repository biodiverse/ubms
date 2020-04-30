get_samples <- function(fit, draws){
  nsamp <- nsamples(fit)
  out <- 1:nsamp
  if(!is.null(draws) && draws < nsamp){
    out <- sample(out, draws)
  }
  out
}
