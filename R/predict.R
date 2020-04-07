#' @export
setMethod("predict", "ubmsFit", function(object, type, random=TRUE, 
                                           summary=TRUE, ...){  
  
  #Use submodels for this
  inp <- build_inputs(object@psiformula, object@pformula, object@data) 

  switch(type,
    state = {
      has_random <- as.logical(inp$stan_data$occ_has_random)
      fix_param <- 'beta_occ'
      rand_param <- 'b_occ'
      X <- inp$stan_data$X_occ
      Z <- inp$stan_data$Z_occ
    },
    det = {
      has_random <- as.logical(inp$stan_data$det_has_random)
      fix_param <- 'beta_det'
      rand_param <- 'b_det'
      X <- inp$stan_data$X_det
      Z <- inp$stan_data$Z_det
    })
  
  beta <- rstan::extract(object@stanfit, fix_param)[[1]]
  lp <- X %*% t(beta)
  
  if(has_random & random){
    b <- rstan::extract(object@stanfit, rand_param)[[1]]
    lp <- lp + Z %*% t(b)
  }

  lp <- plogis(lp)

  if(!summary) return(lp)

  stats <- apply(lp, 1, function(x){
      c(Predicted = mean(x),
        SD = stats::sd(x),
        `2.5%` = as.numeric(stats::quantile(x, 0.025)),
        `97.5%` = as.numeric(stats::quantile(x, 0.975))
        )
      })

  as.data.frame(t(stats))

})


