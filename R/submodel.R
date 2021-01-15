setClass("ubmsSubmodel",
  slots = c(
    name = "character",
    type = "character",
    data = "data.frame",
    formula = "formula",
    link = "character",
    missing = "logical"
  ),
  prototype = list(
    name = NA_character_,
    type = NA_character_,
    data = data.frame(),
    formula = ~1,
    link = NA_character_,
    missing = logical(0)
  )
)

ubmsSubmodel <- function(name, type, data, formula, link){
  out <- new("ubmsSubmodel", name=name, type=type, data=data,
             formula=formula, link=link)
  out@missing <- apply(model.matrix(out), 1, function(x) any(is.na(x)))
  out
}

setClass("ubmsSubmodelTransition", contains = "ubmsSubmodel")

#' @importFrom methods as
ubmsSubmodelTransition <- function(name, type, data, formula, link, T){
  data <- drop_final_year(data, T)
  out <- ubmsSubmodel(name, type, data, formula, link)
  out <- as(out, "ubmsSubmodelTransition")
  if(any(out@missing)){
    stop("Missing values are not allowed in yearlySiteCovs", call.=FALSE)
  }
  out
}

drop_final_year <- function(yr_df, nprimary){
  to_drop <- seq(nprimary, nrow(yr_df), by=nprimary)
  yr_df <- yr_df[-to_drop,,drop=FALSE]
  droplevels(yr_df)
}

setClass("ubmsSubmodelScalar", contains = "ubmsSubmodel")

ubmsSubmodelScalar <- function(name, type, link){
  out <- ubmsSubmodel(name, type, data.frame(1), ~1, link)
  as(out, "ubmsSubmodelScalar")
}

placeholderSubmodel <- function(type){
  ubmsSubmodel("Placeholder", type, data.frame(), ~1, "identity")
}

is_placeholder <- function(submodel){
  length(model.matrix(submodel)) == 0
}

setMethod("model.matrix", "ubmsSubmodel",
          function(object, newdata=NULL, na.rm=FALSE, ...){

  data <- object@data
  formula <- lme4::nobars(object@formula)
  mf <- model.frame(formula, data, na.action=stats::na.pass)

  if(is.null(newdata)){
    out <- model.matrix(formula, mf)
    if(na.rm) out <- out[!object@missing,,drop=FALSE]
    return(out)
  }

  check_newdata(newdata, formula)
  new_mf <- model.frame(stats::terms(mf), newdata, na.action=stats::na.pass,
                        xlev=get_xlev(data, mf))
  model.matrix(formula, new_mf)
})

#Check if all required variables are in newdata
check_newdata <- function(newdata, formula){
  inp_vars <- names(newdata)
  term_vars <- all.vars(formula)
  not_found <- ! term_vars %in% inp_vars
  if(any(not_found)){
    stop(paste0("Required variables not found in newdata: ",
               paste(term_vars[not_found], collapse=", ")), call.=FALSE)
  }
}

get_xlev <- function(data, model_frame){
  fac_col <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_col, levels)
  xlevs[names(xlevs) %in% names(model_frame)]
}

get_reTrms <- function(formula, data, newdata=NULL){
  fb <- lme4::findbars(formula)
  mf <- model.frame(lme4::subbars(formula), data, na.action=stats::na.pass)
  if(is.null(newdata)) return(lme4::mkReTrms(fb, mf))
  new_mf <- model.frame(stats::terms(mf), newdata, na.action=stats::na.pass,
                        xlev=get_xlev(data, mf))
  lme4::mkReTrms(fb, new_mf, drop.unused.levels=FALSE)
}

Z_matrix <- function(object, newdata=NULL, na.rm=FALSE, ...){
  data <- object@data
  formula <- object@formula
  check_formula(formula, data)

  if(is.null(lme4::findbars(formula))) return(matrix(0,0,0))

  Zt <- get_reTrms(formula, data, newdata)$Zt
  Z <- t(as.matrix(Zt))
  if(is.null(newdata) & na.rm){
    Z <- Z[!object@missing,,drop=FALSE]
  }
  Z
}

check_formula <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(invisible())

  char <- paste(deparse(formula))
  if(grepl(":|/", char)){
    stop("Nested random effects (using / and :) are not supported",
         call.=FALSE)
  }
  theta <- get_reTrms(formula, data)$theta
  if(0 %in% theta){
    stop("Correlated slopes and intercepts are not supported. Use || instead of |.",
         call.=FALSE)
  }
}


beta_names <- function(submodel){
  colnames(model.matrix(submodel))
}

b_names <- function(submodel){
  if(!has_random(submodel)) return(NA_character_)
  group <- get_reTrms(submodel@formula, submodel@data)
  group_nms <- names(group$cnms)
  z_nms <- character()
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      z_nms <- c(z_nms, paste0(nms_i, ":", levels(group$flist[[nm]])))
    } else {
      #When there is a random slope for a factor variable
      z_nms <- c(z_nms, c(t(sapply(paste0(nms_i), paste0, ":",
                                   levels(group$flist[[nm]])))))
    }
  }
  z_nms
}

sigma_names <- function(submodel){
  if(!has_random(submodel)) return(NA_character_)
  nms <- get_reTrms(submodel@formula, submodel@data)$cnms
  nms <- paste0(unlist(nms), "|", names(nms))
  nms <- gsub("(Intercept)", "1", nms, fixed=TRUE)
  paste0("sigma [", nms, "]")
}

setGeneric("has_random", function(object){
  standardGeneric("has_random")
})

setMethod("has_random", "ubmsSubmodel", function(object){
  !is.null(lme4::findbars(object@formula))
})

#Check if submodel has intercept term
has_intercept <- function(submodel){
  mm <- model.matrix(submodel)
  "(Intercept)" %in% colnames(mm)
}

#Quickly generate parameter names from ubmsSubmodel
b_par <- function(object){
  paste0("b_", object@type)
}

beta_par <- function(object){
  paste0("beta_", object@type)
}

sig_par <- function(object){
  paste0("sigma_", object@type)
}

setClass("ubmsSubmodelList", slots=c(submodels="list"),
         prototype=list(submodels=list()))

ubmsSubmodelList <- function(...){
  submodels <- list(...)
  names(submodels) <- sapply(submodels, function(x) x@type)
  new("ubmsSubmodelList", submodels=submodels)
}

#' Extract a ubmsSubmodel From a ubmsSubmodelList Object
#'
#' @param x Object of class \code{ubmsSubmodelList}
#' @param i The name of a submodel
#'
#' @return An object of class \code{ubmsSubmodel}.
#'
#' @export
setMethod("[", c("ubmsSubmodelList", "character", "missing", "missing"),
  function(x, i){

  types <- names(x@submodels)
  if(! i %in% types){
    stop(paste("Possible types are:", paste(types, collapse=", ")),
         call. = FALSE)
  }
  x@submodels[[i]]
})
