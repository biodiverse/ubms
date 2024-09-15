#' Get Information for a Restricted Spatial Regression Model
#'
#' A call to \code{RSR} in the formula for a state process adds an autocorrelated
#' spatial random effect to the model in the form of a Restricted Spatial
#' Regression (RSR). For examples of RSRs applied to ecological data, see
#' Johnson et al. (2013) and Broms et al. (2014).
#' The function can also be used outside a formula to visualize the spatial
#' relationships between sites in your data and choose an appropriate
#' distance threshold below which two sites will be considered neighbors, and
#' thus potentially correlated, for the RSR model. For more details, see the example vignette:
#' \code{vignette("spatial-models", package = "ubms")}
#'
#' @param x A vector of coordinates (should be projected)
#' @param y An (optional) second vector of coordinates
#' @param threshold The distance cutoff below which two sites will be
#'  considered neighbors. Should be the same units as the coordinates.
#' @param moran_cut The number of eigenvectors to use in the RSR. The possible
#'  range of values is between 1 and the number of sites. Smaller numbers will
#'  result in faster runtime and smoother map output, and vice-versa. If
#'  not provided (the default), the number of eigenvectors will be set as
#'  10\% of the number of sites which is usually appropriate.
#' @param plot_site If a site number (as an integer) is supplied, the function
#'  returns a plot showing that site and its neighbors under the current settings.
#'  Useful for deciding what to set your threshold at.
#'
#' @return Either a list of spatial matrices used for the RSR (only useful
#'  internally to ubms), or if \code{plot_site} is an integer, a \code{ggplot} object.
#'
#' @examples
#' # Generate some coordinates
#' x <- runif(100,0,10)
#' y <- runif(100,0,10)
#' # Show neighbors of site 10 when threshold is 3 units
#' RSR(x, y, threshold=3, plot_site=10)
#'
#' @references Broms KM, Johnson DS, Altwegg R, Conquest LL. 2014. Spatial
#'  occupancy models applied to atlas data show Southern Ground Hornbills strongly
#'  depend on protected areas. Ecological Applications 24: 363-374.
#'
#' Johnson DS, Conn PB, Hooten MB, Ray JC, Pond BA. 2013. Spatial occupancy
#'  models for large data sets. Ecology 94: 801-808.
#'
#' @importFrom ggplot2 scale_color_manual ggtitle
#' @export
RSR <- function(x, y=NULL, threshold, moran_cut=NULL, plot_site=NULL){
  coords <- cbind(x,y)
  distmat <- as.matrix(stats::dist(coords, upper=TRUE))
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

  list(A=A, Q=Q, n_eig=n_eig, coords=coords)
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

  ggplot(focal, aes(x=.data[["x"]],y=.data[["y"]])) +
    geom_point(data=plot_dat, alpha=0.3) +
    geom_point(aes(col=.data[["type"]])) +
    scale_color_manual(values=c("red","blue")) +
    plot_theme() +
    theme(legend.text=element_text(size=14),
          legend.title=element_blank()) +
    ggtitle(paste0("Threshold = ",threshold))
}

get_rsr_info <- function(object){
  form <- object@spatial

  final_rows <- nrow(object@data) + nrow(object@data_aug)
  data <- rbind(object@data, object@data_aug)

  fc <- as.character(form)[2]

  parts <- attr(stats::terms(form), "term.labels")

  has_RSR <- sapply(parts, function(x) grepl("RSR(", x, fixed=TRUE))
  stopifnot(sum(has_RSR) == 1)

  func <- parts[has_RSR]
  with(data, eval(parse(text=func)))
}

remove_RSR <- function(form){
  fc <- as.character(form)[2]
  fc <- gsub(" ", "", fc)
  if(grepl("|",fc,fixed=TRUE)){
    stop("Can't have both regular and spatial random effect in model")
  }
  rem <- gsub("\\+?RSR\\((.*)\\)", "", fc)
  if(rem == "") rem <- "1"
  out <- as.formula(paste0("~",rem))
  out
}

setGeneric("has_spatial", function(object, ...) standardGeneric("has_spatial"))

setMethod("has_spatial", "list", function(object, support=TRUE, ...){
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
  if(!support){
    stop("This model type does not support spatial components", call.=FALSE)
  }
  return(TRUE)
})

setMethod("has_spatial", "ubmsSubmodel", function(object, ...){
  methods::.hasSlot(object, "spatial")
})

setMethod("has_spatial", "ubmsFit", function(object, ...){
  any(sapply(object@submodels@submodels, has_spatial))
})

setClass("ubmsSubmodelSpatial", contains = "ubmsSubmodel",
          slots=c(data_aug="data.frame", sites_aug="logical",
                  spatial="formula"))

ubmsSubmodelSpatial <- function(name, type, data, formula, link, prior_intercept,
                                prior_coef, prior_sigma, sites_aug, data_aug){
  form_no_rsr <- remove_RSR(formula)
  out <- new("ubmsSubmodelSpatial", name=name, type=type, data=data,
             formula=form_no_rsr, link=link,
             prior_intercept=prior_intercept, prior_coef=prior_coef,
             prior_sigma=prior_sigma,
             data_aug=data_aug, sites_aug=sites_aug, spatial=formula)
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

  obs_augment <- rep(sites_augment, each=unmarked::obsNum(umf))
  obs_cov_noaug <- obsCovs(umf)[!obs_augment,,drop=FALSE]
  y_noaug <- getY(umf)[!sites_augment,,drop=FALSE]
  umf@y <- y_noaug
  siteCovs(umf) <- site_cov_noaug
  if(!is.null(obsCovs(umf))) obsCovs(umf) <- obs_cov_noaug
  list(umf=umf, data_aug=site_cov_aug, sites_augment=sites_augment)
}

setGeneric("spatial_matrices", function(object, ...) standardGeneric("spatial_matrices"))

#' @importFrom RSpectra eigs
setMethod("spatial_matrices", "ubmsSubmodelSpatial", function(object, ...){
  message("Building RSR matrices")
  form <- object@spatial

  final_rows <- nrow(object@data) + nrow(object@data_aug)
  data <- rbind(object@data, object@data_aug)

  rsr_info <- get_rsr_info(object)

  A <- rsr_info$A
  nr <- remove_RSR(form)
  X <- model.matrix(nr, data)
  P <- diag(nrow(X)) - X %*% solve(crossprod(X), t(X))
  Op <- (nrow(A)/sum(A))*(P %*% (A %*% P))
  K <- RSpectra::eigs(Op, as.integer(rsr_info$n_eig))$vectors
  if (all(Im(K <- zapsmall(K))==0)) K <- matrix(as.numeric(K), nrow=nrow(K))
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
  offset_aug <- model_offset(object, newdata=object@data_aug)
  out <- c(out, sm, list(X_aug=X_aug, n_aug_sites=nrow(X_aug), offset_aug=offset_aug))
  out
})

setMethod("stanfit_names", "ubmsSubmodelSpatial", function(object, ...){
  out <- paste0("beta_",object@type,'[',beta_names(object),']')
  form <- object@spatial
  final_rows <- nrow(object@data) + nrow(object@data_aug)
  data <- rbind(object@data, object@data_aug)
  fc <- as.character(form)[2]
  parts <- attr(stats::terms(form), "term.labels")
  has_RSR <- sapply(parts, function(x) grepl("RSR(", x, fixed=TRUE))
  stopifnot(sum(has_RSR) == 1)
  func <- parts[has_RSR]
  rsr_info <- with(data, eval(parse(text=func)))
  tn <- paste0("b_",object@type,"[theta[",1:rsr_info$n_eig,"]]")
  c(out, tn, "tau")
})

#' Plot A Map of the State Parameter Based on a Spatial ubms Model
#'
#' @param object A \code{ubmsFit} model with a spatial random effect
#' @param param The parameter to plot, either \code{"state"} for, e.g.,
#'  mean occupancy or abundance, or \code{"eta"} for the random effect itself
#' @param sites If \code{TRUE}, also plot the locations of sites that
#'  were sampled on the map and if had a detection of the species
#' @param cell_size The size of each side of the (square) cells drawn in the map,
#'  in the same units as the coordinates. If \code{NULL}, \code{ubms} will try
#'  to pick a reasonable cell size for you
#'
#' @importFrom ggplot2 geom_tile scale_fill_gradientn scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_color_manual
#' @importFrom grDevices terrain.colors
#' @export
plot_spatial <- function(object, param=c('state','eta'), sites=TRUE, cell_size=NULL){
  if(!inherits(object, "ubmsFit")){
    stop("Requires ubmsFit object", call.=FALSE)
  }
  param <- match.arg(param)
  sm <- object["state"]
  if(!has_spatial(sm)){
    stop("No spatial random effect in model", call.=FALSE)
  }
  coords <- get_rsr_info(sm)$coords
  nms <- colnames(coords)
  sampled <- !sm@sites_aug

  if(is.null(cell_size)){
    pd <- as.matrix(stats::dist(coords))
    diag(pd) <- NA
    min_dist <- apply(pd, 1, min, na.rm=TRUE)
    cell_size <- mean(min_dist, na.rm=TRUE)
  }

  if(param == "state"){
    est_raw <- predict(object, "state")$Predicted
    est <- c(est_raw[sampled], est_raw[!sampled])
    param <- "lambda"
    if(inherits(object, c("ubmsFitOccu","ubmsFitOccuTTD"))) param <- "psi"
    if(inherits(object, "ubmsFitOccuRN")) param <- "lambda"
  } else if(param == "eta"){
    b <- extract(object, "b_state")$b_state
    Kmat <- spatial_matrices(sm)$Kmat
    est <- Kmat %*% colMeans(b)
  }
  plot_data <- cbind(as.data.frame(coords), est=est)

  out <- ggplot(data=plot_data, aes(x=.data[[nms[1]]], y=.data[[nms[2]]])) +
    geom_tile(aes(fill=.data[["est"]]), width=cell_size, height=cell_size) +
    scale_fill_gradientn(colors=terrain.colors(10)) +
    labs(fill=param) +
    plot_theme() +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

  if(sites){
    coords_samp <- as.data.frame(coords)[1:nrow(sm@data),]
    coords_samp$obs <- get_Kmin(object@response)
    if(inherits(object, "ubmsFitOccuTTD")){
      coords_samp$obs <- as.numeric(coords_samp$obs < object@response@surveyLength)
    }

    if(inherits(object, c("ubmsFitOccu","ubmsFitOccuTTD","ubmsFitOccuRN"))){
      coords_samp$obs <- factor(coords_samp$obs)
      out <- out + geom_point(data=coords_samp, aes(col=.data[["obs"]]),
                            size=1, pch=19) +
        scale_color_manual(values=c("gray","black")) +
        labs(color="Detected")
    } else {
      out <- out + geom_point(data=coords_samp, aes(size=.data[["obs"]]), pch=19) +
        #scale_color_manual(values=c("gray","black")) +
        labs(size="Minimum\ncount")
    }
  }
  out
}
