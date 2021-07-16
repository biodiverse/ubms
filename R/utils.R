get_samples <- function(fit, draws){
  nsamp <- nsamples(fit)
  out <- 1:nsamp
  if(!is.null(draws) && draws < nsamp){
    out <- sample(out, draws)
  }
  out
}

submodel_types <- function(object){
  names(object@submodels@submodels)
}

check_type <- function(type, possible_types){
  if(! type %in% possible_types){
    stop(paste("Submodel must be one of:", paste(possible_types, collapse=", ")),
         call.=FALSE)
  }
}

plot_theme <- function(){
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.background=element_blank(),
          strip.text=element_blank()
    )
}

# Remove offset term from a formula
remove_offset <- function(form){
  char <- paste(deparse(form), collapse="")
  char <- gsub(" ", "", char)
  out <- gsub("\\+?offset\\((.*?)\\)", "", char)
  if(out == "~") out <- "~1"
  as.formula(out)
}
