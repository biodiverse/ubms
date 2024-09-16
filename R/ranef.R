#' Extract Random Effects
#'
#' Extract random effects from a \code{ubmsFit} model. Note that this function
#' works like \code{ranef} for \code{merMod} objects from \code{lme4}, not like
#' \code{ranef} for \code{unmarkedFit} objects. To get functionality similar
#' to that of \code{unmarkedFit}, use \code{posterior_predict}.
#'
#' Note: by default this function adds the overall intercept or slope
#' to the (mean-0) random effect to get the complete random intercept or slope.
#' In this way the output is more like the output of \code{lme4::coef}
#' and not \code{lme4::ranef}. You can turn this off and return just the
#' mean-0 random effect by setting argument \code{add_mean = FALSE}.
#'
#' If you run \code{ranef} on a submodel with a spatial random effect,
#' the function will return estimates of parameter \code{eta}.
#'
#' @param object A fitted model of class \code{ubmsFit}
#' @param submodel The name of the submodel, as a character string, for
#'  which to generate the random effects
#' @param summary If \code{TRUE}, calculate mean, SD, and 95% uncertainty interval
#'  for each random effect term
#' @param add_mean If \code{TRUE} (the default) add the overall intercept or
#'  slope mean and return the complete random intercept or slope.
#' @param ... Currently ignored
#'
#' @return If \code{summary=FALSE}, a list of random effect values; if
#'  \code{TRUE}, a data frame with columns for random effect mean, SD, and
#'  95% uncertainty interval lower and upper bounds.
#'
#' @aliases ranef
#' @seealso \code{\link[lme4]{ranef}}, \code{\link{posterior_predict}}
#' @include fit.R
#' @importFrom unmarked ranef
#' @export
setMethod("ranef","ubmsFit", function(object, submodel, summary=FALSE,
          add_mean = TRUE, ...){

  sm <- object[submodel]

  qu <- function(x, q) as.numeric(stats::quantile(x, q))

  # If this is a spatial submodel return eta
  if(has_spatial(sm)){
    b <- extract(object, "b_state")$b_state
    Kmat <- spatial_matrices(sm)$Kmat
    
    if(summary){
      eta_post <- Kmat %*% t(b)
      out <- data.frame(
        Estimate=rowMeans(eta_post),
        SD=apply(eta_post, 1, stats::sd),
        `2.5%`=apply(eta_post, 1, qu, q=0.025),
        `97.5%`=apply(eta_post, 1, qu, q=0.975),
        check.names=FALSE
      )
    } else {
      out <- drop(Kmat %*% colMeans(b))
    }
    return(list(eta = out))
  }

  if(!has_random(sm)){
    stop("No random effects terms in this submodel", call.=FALSE)
  }
  re <- get_reTrms(sm@formula, sm@data)
  fl <- re$cnms
  bn <- b_names(sm)

  ran <- lapply(1:length(fl), function(i){

    fac <- names(fl)[i]
    trm <- fl[[i]]
    if(length(trm)>1) stop("There should only be one term here")

    b_ind <- re$Gp[i:(i+1)] + c(1,0)
    re_samples <- extract(object, paste0("b_",submodel))[[1]]
    re_samples <- re_samples[,b_ind[1]:b_ind[2]]

    #Add mean value if this is an effects parameterization
    if(trm %in% beta_names(sm) & add_mean){
      message("Adding the mean to get the complete random slope/intercept")
      beta_ind <- which(beta_names(sm) == trm)
      mn_samples <- extract(object, paste0("beta_",submodel))[[1]]
      mn_samples <- mn_samples[,beta_ind]
      re_samples <- re_samples + mn_samples
    }

    fac_lvls <- levels(re$flist[[fac]])
    if(summary){
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
  combine_same_name(ran)
})

combine_same_name <- function(inp_list){
  lnames <- unique(names(inp_list))
  sapply(lnames, function(x){
           list_sub <- inp_list[names(inp_list)==x]
           comb <- unlist(list_sub, recursive=FALSE, use.names=FALSE)
           names(comb) <- sapply(list_sub, names)
           comb
  }, simplify=FALSE)
}
