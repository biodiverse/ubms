#' @include fit.R 
#' @export
setMethod("predict", "ubmsFit", 
  function(object, type, fun=mean, random=TRUE, summary=TRUE, ...){  
  
  sm <- object[type]
  include_random <- !is.na(sm@sigma_names) & random

  beta <- extract(object, paste0('beta_',type))[[1]]
  lp <- model.matrix(sm) %*% t(beta)
  
  if(include_random){
    b <- extract(object, paste0('b_', type))[[1]]
    lp <- lp + sm@Z %*% t(b)
  }

  lp <- do.call(sm@link , list(lp))

  if(!summary) return(lp)

  stats <- apply(lp, 1, function(x){
      c(Predicted = fun(x),
        SD = stats::sd(x),
        `2.5%` = as.numeric(stats::quantile(x, 0.025)),
        `97.5%` = as.numeric(stats::quantile(x, 0.975))
        )
      })

  as.data.frame(t(stats))

})

