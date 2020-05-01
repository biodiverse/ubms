#' @include fit.R
#' @importFrom lme4 ranef
#' @export
setMethod("ranef","ubmsFit", function(object, submodel, summary=FALSE, ...){ 

  sm <- object[submodel]
  if(!has_random(sm)){
    stop("No random effects terms in this submodel", call.=FALSE)
  }
  re <- get_reTrms(sm@formula, sm@data)
  fl <- re$cnms
  bn <- b_names(sm)

  ran <- lapply(1:length(fl), function(i){
  
    fac <- names(fl)[i]
    trm <- fl[[i]]
    if(length(trm)>1){
      stop("There should only be one term here")
    }

    beta_ind <- which(beta_names(sm) == trm)
    mn_samples <- extract(object, paste0("beta_",submodel))[[1]]
    mn_samples <- mn_samples[,beta_ind]

    b_ind <- re$Gp[i:(i+1)] + c(1,0)
    b_samples <- extract(object, paste0("b_",submodel))[[1]]
    b_samples <- b_samples[,b_ind[1]:b_ind[2]]
    re_samples <- b_samples + mn_samples
    
    fac_lvls <- levels(re$flist[[fac]])
    if(summary){
      qu <- function(x, q) as.numeric(stats::quantile(x, q))
      out <- data.frame(
        Estimate=colMeans(re_samples),
        SD=apply(re_samples, 2, stats::sd),
        `2.5%`=apply(re_samples, 2, qu, q=0.025),
        `97.5%`=apply(re_samples, 2, qu, q=0.975),
        check.names=FALSE
      )
      rownames(out) <- fac_lvls
    } else {
      out <- colMeans(re_samples)
      names(out) <-  fac_lvls
    }
    out <- list(out)
    names(out) <- trm
    out

  })
  names(ran) <- names(fl)
  ran
})
