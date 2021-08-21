check_prior_input <- function(input){
  stopifnot(length(input) %in% c(1,2))
  if(length(input) == 1){
    check_names <- "coef"
  } else {
    check_names <- c("intercept","coef")
  }
  if(!all(check_names %in% names(input))){
    stop(paste0("List of priors for each submodel must have named elements: ",
                paste(check_names, collapse=", ")),
         call.=FALSE)
  }
}

get_prior_info <- function(input, Xmat){

  check_prior_input(input)

  locations <- rep(input$coef$location, ncol(Xmat))
  scales <- rep(input$coef$scale, ncol(Xmat))
  prior_type <- c(1,1) # only normal supported at moment

  if(input$coef$autoscale){
    for (i in 1:ncol(Xmat)){

      if(all(unique(Xmat[,i]) %in% c(0,1))){
        next
      }
      scales[i] <- scales[i] * 1/sd(Xmat[,i])
    }
  }

  if(!is.null(input$intercept) && colnames(Xmat)[1] == "(Intercept)"){
    locations[1] <- input$intercept$location
    scales[1] <- input$intercept$scale
  }

  list(locations=locations, scales=scales, prior_type=prior_type)

}


normal <- function(location=0, scale=10, autoscale=TRUE){
  stopifnot(scale > 0)
  list(dist="normal", location=location, scale=scale, autoscale=autoscale)
}

default_priors <- function(){
  out <- lapply(1:2, function(x){
           list(intercept=normal(0, 10, autoscale=TRUE),
           coef=normal(0, 2.5, autoscale=TRUE))
         })
  names(out) <- c("state","det")
  out
}
