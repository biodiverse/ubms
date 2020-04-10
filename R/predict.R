#' @include fit.R
#' @importFrom unmarked predict
#' @export
setMethod("predict", "ubmsFit", 
  function(object, type, fun=mean, random=TRUE, summary=TRUE, ...){  
  
  sm <- object[type]

  beta <- extract(object, beta_par(sm))[[1]]
  lp <- model.matrix(sm) %*% t(beta)
  
  if(has_random(sm) & random){
    b <- extract(object, b_par(sm))[[1]]
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
