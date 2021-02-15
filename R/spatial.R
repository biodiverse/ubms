#' @importFrom ggplot2 scale_color_manual ggtitle
#' @export
RSR <- function(x, y=NULL, threshold, moran_cut=NULL, plot_site=NULL){
  coords <- cbind(x,y)
  distmat <- as.matrix(dist(coords, upper=TRUE))
  A = matrix(0, nrow(distmat), ncol(distmat))
  A[distmat <= threshold] <- 1
  diag(A) <- 0

  if(!is.null(plot_site)){
    return(plot_RSR(coords, A, threshold, plot_site))
  }

  n_eig <- moran_cut
  if(is.null(n_eig)){
    n_eig <- round(0.1*nrow(A))
  }
  stopifnot(n_eig < nrow(A))

  Q <- A
  Q[Q==1] <- -1
  diag(Q) <- -apply(Q,1,sum)

  list(A=A, Q=Q, n_eig=n_eig)
}

plot_RSR <- function(coords, A, threshold, focal_site){
  plot_dat <- as.data.frame(coords)
  neighbors <- plot_dat[A[focal_site,]==1,,drop=FALSE]
  focal <- plot_dat[focal_site,,drop=FALSE]
  focal$type <- "Focal site"
  if(nrow(neighbors)>0){
    neighbors$type <- "Neighbors"
    focal <- rbind(focal, neighbors)
  }

  ggplot(focal, aes_string(x="x",y="y")) +
    geom_point(data=plot_dat, alpha=0.3) +
    geom_point(aes_string(col="type")) +
    scale_color_manual(values=c("red","blue")) +
    plot_theme() +
    theme(legend.text=element_text(size=14),
          legend.title=element_blank()) +
    ggtitle(paste0("Threshold = ",threshold))
}

remove_RSR <- function(form){
  fc <- as.character(form)[2]
  fc <- gsub(" ", "", fc)
  rem <- gsub("\\+RSR\\((.*)\\)", "", fc)
  out <- as.formula(paste0("~",rem))
  if(length(lme4::findbars(out)) > 0){
    stop("Can't have both regular and spatial random effect in model")
  }
  out
}

setGeneric("has_spatial", function(object, ...) standardGeneric("has_spatial"))

setMethod("has_spatial", "list", function(object, support=TRUE, ...){
  if(!support){
    stop("This model type does not support spatial components", call.=FALSE)
  }
  hs <- sapply(object, function(x){
    any(grepl("RSR(", as.character(x), fixed=TRUE))
  })
  hs <- hs[hs]
  if(length(hs) == 0) return(FALSE)
  if(length(hs) != 1){
    stop("Can only have spatial component in one formula", call.=FALSE)
  }
  possible_mods <- c("state","psi")
  if(!names(hs) %in% possible_mods){
    stop("Can only put spatial components on state model", call.=FALSE)
  }
  return(TRUE)
})

setMethod("has_spatial", "ubmsSubmodel", function(object, ...){
  .hasSlot(object, "spatial")
})

setClass("ubmsSubmodelSpatial", contains = "ubmsSubmodel",
          slots=c(data_aug="data.frame", sites_aug="logical",
                  spatial="formula"))

ubmsSubmodelSpatial <- function(name, type, data, formula, link,
                                sites_aug, data_aug){
  form_no_rsr <- remove_RSR(formula)
  out <- new("ubmsSubmodelSpatial", name=name, type=type, data=data,
             formula=form_no_rsr, link=link, data_aug=data_aug, sites_aug=sites_aug,
             spatial=formula)
  out@missing <- apply(model.matrix(out), 1, function(x) any(is.na(x)))
  out
}

extract_missing_sites <- function(umf){
  sites_augment <- apply(getY(umf), 1, function(x) all(is.na(x)))
  site_cov_noaug <- siteCovs(umf)[!sites_augment,,drop=FALSE]
  site_cov_aug <- siteCovs(umf)[sites_augment,,drop=FALSE]
  if(any(is.na(site_cov_aug))){
    stop("Missing values not allowed in site covariates for unobserved sites", call.=FALSE)
  }

  obs_augment <- rep(sites_augment, each=obsNum(umf))
  obs_cov_noaug <- obsCovs(umf)[!obs_augment,,drop=FALSE]
  y_noaug <- getY(umf)[!sites_augment,,drop=FALSE]
  umf@y <- y_noaug
  siteCovs(umf) <- site_cov_noaug
  obsCovs(umf) <- obs_cov_noaug
  list(umf=umf, data_aug=site_cov_aug, sites_augment=sites_augment)
}

setGeneric("spatial_matrices", function(object, ...) standardGeneric("spatial_matrices"))

#' @importFrom RSpectra eigs
setMethod("spatial_matrices", "ubmsSubmodelSpatial", function(object, ...){
  message("Building RSR matrices")
  form <- object@spatial

  final_rows <- nrow(object@data) + nrow(object@data_aug)
  data <- rbind(object@data, object@data_aug)

  fc <- as.character(form)[2]

  parts <- attr(terms(form), "term.labels")

  has_RSR <- sapply(parts, function(x) grepl("RSR(", x, fixed=TRUE))
  stopifnot(sum(has_RSR) == 1)

  func <- parts[has_RSR]
  rsr_info <- with(data, eval(parse(text=func)))

  A <- rsr_info$A
  nr <- remove_RSR(form)
  X <- model.matrix(nr, data)
  P <- diag(nrow(X)) - X %*% solve(crossprod(X), t(X))
  Op <- (nrow(A)/sum(A))*(P %*% (A %*% P))
  K <- RSpectra::eigs(Op, as.integer(rsr_info$n_eig))$vectors
  Q <- rsr_info$Q
  Qalpha <- as.matrix(t(K) %*% Q %*% K)
  stopifnot(ncol(K) == nrow(Qalpha))
  list(Kmat=K, Qalpha=Qalpha, n_eigen=nrow(Qalpha))
})

setMethod("get_pars", "ubmsSubmodelSpatial", function(object, ...){
  c(paste0(c("beta","b"),"_",object@type), "tau")
})

setMethod("get_stan_data", "ubmsSubmodelSpatial", function(object, ...){
  out <- callNextMethod(object, ...)
  sm <- spatial_matrices(object)
  out$n_random_state <- as.array(sm$n_eig)
  X_aug <- model.matrix(object, newdata=object@data_aug)
  out <- c(out, sm, list(X_aug=X_aug, n_aug_sites=nrow(X_aug)))
  out
})

setMethod("stanfit_names", "ubmsSubmodelSpatial", function(object, ...){
  out <- paste0("beta_",object@type,'[',beta_names(object),']')
  form <- object@spatial
  final_rows <- nrow(object@data) + nrow(object@data_aug)
  data <- rbind(object@data, object@data_aug)
  fc <- as.character(form)[2]
  parts <- attr(terms(form), "term.labels")
  has_RSR <- sapply(parts, function(x) grepl("RSR(", x, fixed=TRUE))
  stopifnot(sum(has_RSR) == 1)
  func <- parts[has_RSR]
  rsr_info <- with(data, eval(parse(text=func)))
  tn <- paste0("b_",object@type,"[theta[",1:rsr_info$n_eig,"]]")
  c(out, tn, "tau")
})
