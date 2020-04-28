setClass("ubmsGOF", slots=c(bpvals="data.frame", samples="data.frame"))

#Constructor
ubmsGOF <- function(df_list){
  bpvals <- lapply(df_list, calc_bpval) 
  new("ubmsGOF",
      bpvals = do.call("rbind", bpvals),
      samples = do.call("rbind", df_list))
}

calc_bpval <- function(df){
  data.frame(Statistic=df$stat[1], `Bayesian P-val`=mean(df$sim > df$obs),
             check.names=FALSE)
}

setMethod("show", "ubmsGOF", function(object) print(object@bpvals))

#' @importFrom ggplot2 ggplot aes geom_abline geom_point theme_bw labs
#' @importFrom ggplot2 facet_wrap theme element_blank element_text element_rect
#' @importFrom ggplot2 geom_label unit
setMethod("plot", "ubmsGOF", function(x, ...){

  bp_dat <- x@bpvals
  bp_dat$lab <- paste("P =", round(bp_dat[,2], 2))
  names(bp_dat)[1] <- 'stat'

  ggplot(x@samples, aes_string(x="obs", y="sim")) +
    geom_abline(aes(intercept=0, slope=1),size=1.2, col='red') +
    geom_point(alpha=0.4) +
    theme_bw() +
    labs(y="Simulated data", x="Observed data") +
    facet_wrap(~stat, scales="free") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=14),
          strip.background=element_rect(fill="transparent"),
          strip.text=element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    geom_label(data=bp_dat, aes_string(x=-Inf, y=Inf, label="lab"),
              hjust=-0.2, vjust=1.4, size=5,
              fill='white', label.size=0, 
              label.padding=unit(0.1, "lines"))
})


#' @export
setGeneric("gof", function(object, nsim=NULL, ...){
             standardGeneric("gof")})

#' @include occu.R
setMethod("gof", "ubmsFitOccu", function(object, draws=NULL, ...){

  y <- getY(object@data)
  ylong <- as.vector(t(y))

  M <- nrow(y)
  J <- ncol(y)

  nsamp <- nsamples(object)
  samples <- 1:nsamp
  if(!is.null(draws) && draws < nsamp){
    samples <- sample(1:nsamp, draws)
  }

  z <- sim_z(object, samples=samples, re.form=NULL)
  ysim_long <- t(sim_y(object, samples=samples, re.form=NULL, z=z))
  z <- t(z)
  p <- t(sim_lp(object, "det", transform=TRUE, newdata=NULL, 
                samples=samples, re.form=NULL))
  
  zp <- z[rep(1:nrow(z), each=J),] * p

  #Deviance
  dev_sim <- dbinom(as.vector(ysim_long), 1, as.vector(zp), log=TRUE)   
  dev_sim <- matrix(dev_sim, nrow=M*J)
  dev_sim <- apply(dev_sim, 2, function(x) -2 * sum(x, na.rm=TRUE)) 
  dev_obs <- apply(zp, 2, function(x){ 
                   -2*sum(dbinom(ylong, 1, x, log=TRUE), na.rm=TRUE)})
  dev_df <- data.frame(stat="Deviance", sim=dev_sim, obs=dev_obs)
  
  #Freeman-Tukey
  ft_sim <- apply((sqrt(ysim_long) - sqrt(zp))^2, 2, sum, na.rm=TRUE)
  ft_obs <- apply(zp, 2, function(x) sum((sqrt(ylong) - sqrt(x))^2, na.rm=TRUE))
  ft_df <- data.frame(stat="Freeman-Tukey", sim=ft_sim, obs=ft_obs)
  
  ubmsGOF(list(dev_df, ft_df))
})


#' @include pcount.R
setMethod("gof", "ubmsFitPcount", function(object, draws=NULL, ...){

  y <- getY(object@data)
  ylong <- as.vector(t(y))

  M <- nrow(y)
  J <- ncol(y)

  nsamp <- nsamples(object)
  samples <- 1:nsamp
  if(!is.null(draws) && draws < nsamp){
    samples <- sample(1:nsamp, draws)
  }  
  draws <- length(samples)

  z <- sim_z(object, samples=samples, re.form=NULL)
  ysim_long <- t(sim_y(object, samples=samples, re.form=NULL, z=z))
  z <- t(z)
  p <- t(sim_lp(object, "det", transform=TRUE, newdata=NULL, 
                samples=samples, re.form=NULL))
  N <- z[rep(1:nrow(z), each=J),]

  #Deviance
  dev_sim <- dbinom(as.vector(ysim_long), as.vector(N), as.vector(p), log=TRUE)   
  dev_sim <- matrix(dev_sim, nrow=M*J)
  dev_sim <- apply(dev_sim, 2, function(x) -2 * sum(x, na.rm=TRUE))  
  dev_obs <- rep(NA, draws)
  for (i in 1:draws){
    dev_obs[i] <- -2*sum(dbinom(ylong, N[,i], p[,i], log=TRUE), na.rm=TRUE) 
  }
  dev_df <- data.frame(stat="Deviance", sim=dev_sim, obs=dev_obs)
  
  #Freeman-Tukey
  Np <- N * p
  ft_sim <- apply((sqrt(ysim_long) - sqrt(Np))^2, 2, sum, na.rm=TRUE)
  ft_obs <- apply(Np, 2, function(x) sum((sqrt(ylong) - sqrt(x))^2, na.rm=TRUE))
  ft_df <- data.frame(stat="Freeman-Tukey", sim=ft_sim, obs=ft_obs)
  
  ubmsGOF(list(dev_df, ft_df))
})
