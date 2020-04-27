#' @include fit.R
#' @importFrom unmarked predict
#' @export
setMethod("predict", "ubmsFit", 
  function(object, type, random=TRUE, summary=TRUE, ...){  
 
  sm <- object[type]

  beta <- extract(object, beta_par(sm))[[1]]
  lp <- model.matrix(sm) %*% t(beta)
 
  if(has_random(sm) & random){
    b <- extract(object, b_par(sm))[[1]]
    lp <- lp + Z_matrix(sm) %*% t(b)
  }

  lp <- do.call(sm@link , list(lp))
  if(!summary) return(lp)

  stats <- apply(lp, 1, function(x){
        quant <- as.numeric(stats::quantile(x, c(0.025,0.975), na.rm=TRUE))
        c(Predicted = mean(x), SD = stats::sd(x), 
          `2.5%` = quant[1], `97.5%` = quant[2])
      })

  as.data.frame(t(stats))
})
