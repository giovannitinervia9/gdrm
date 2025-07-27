#' Predicted values on link scale
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list of predicted values on link scale.
#' @export
gdrm_predict <- function(mod_comp) {
  lapply(mod_comp, function(mod_par) {
    comps <- lapply(mod_par, function(comp) comp$fitted(comp$par))
    Reduce(`+`, comps)
  })
}


#' Predicted values on parameter scale
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list of predicted values on parameter scale.
#' @export
gdrm_fitted <- function(mod_comp, distrib) {
  eta <- gdrm_predict(mod_comp)
  link <- distrib$link_list
  Map(function(eta, link) link$linkinv(eta), eta = eta, link = link)
}


#' Coefficient of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the coefficients for each component of the model.
#' @export
gdrm_coef <- function(mod_comp) {
  par_list <- vector("list", length(mod_comp))
  names(par_list) <- names(mod_comp)
  
  for (j in 1:length(par_list)) {
    comp <- mod_comp[[j]]
    for (i in 1:length(comp)) {
      par_list[[j]][[i]] <- comp[[i]]$par
      names(par_list[[j]])[i] <- names(comp)[i]
    }
  }
  
  par_list
}


#' Extract coefficient of gdrm model in a numerical vector
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A numeric vector containing the coefficients of the gdrm model.
#'
#' @export
gdrm_coef_vector <- function(mod_comp) {
  coef_vector <- c()
  
  for (i in 1:length(mod_comp)) {
    comp <- mod_comp[[i]]
    
    for (j in 1:length(comp)) {
      build <- comp[[j]]
      # Each build always has a 'par' element
      coef_vector <- c(coef_vector, build$par)
    }
  }
  
  return(coef_vector)
}


#' Hyperparameters of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the hyperparameters for each component of the model.
#' @export
gdrm_hyperpar <- function(mod_comp) {
  par_list <- vector("list", length(mod_comp))
  names(par_list) <- names(mod_comp)

  for (j in 1:length(par_list)) {
    comp <- mod_comp[[j]]
    for (i in 1:length(comp)) {
      par_list[[j]][[i]] <- comp[[i]]$hyperpar
      names(par_list[[j]])[i] <- names(comp)[i]
    }
  }

  par_list
}


#' Extract hyperparameters of gdrm model in a numerical vector
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A numeric vector containing the hyperparameters of the gdrm model.
#'
#' @export
gdrm_hyperpar_vector <- function(mod_comp) {
  coef_vector <- c()
  
  for (i in 1:length(mod_comp)) {
    comp <- mod_comp[[i]]
    
    for (j in 1:length(comp)) {
      build <- comp[[j]]
      # Each build always has a 'par' element
      coef_vector <- c(coef_vector, build$hyperpar)
    }
  }
  
  unlist(coef_vector)
}


#' Build penalty matrix for gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param hyperpar Logical. If `TRUE` (default), the elements of the penalty matrix are multiplied by the respective hyperparameter.
#'
#' @returns A list of penalty matrices for each model parameter.
#' @export
gdrm_penalty <- function(mod_comp, hyperpar = TRUE) {
  npar <- length(mod_comp)
  Plist <- vector("list", npar)
  names(Plist) <- names(mod_comp)

  for (i in 1:npar) {
    ncomp <- length(mod_comp[[i]])
    for (j in 1:ncomp) {
      comp <- mod_comp[[i]][[j]]

      if (inherits(comp, c("linear", "ridge", "nl"))) {
        if (hyperpar) {
          Plist[[i]][[j]] <- comp$hyperpar * comp$P
        } else {
          Plist[[i]][[j]] <- comp$P
        }
      } else if (inherits(comp, c("smooth"))) {
        if (hyperpar) {
          P <- Map(
            function(Psub, hsub) {
              # multiply each corresponding element of the sublists
              Map(function(Pmat, hval) Pmat * hval, Psub, hsub)
            },
            comp$P,
            comp$hyperpar
          )

          P <- lapply(P, function(Pj) Reduce(`+`, Pj))

          Plist[[i]][[j]] <- as.matrix(Matrix::bdiag(P))
        } else {
          P <- lapply(comp$P, function(Pj) Reduce(`+`, Pj))
          Plist[[i]][[j]] <- as.matrix(Matrix::bdiag(P))
        }
      }

      names(Plist[[i]])[[j]] <- names(mod_comp[[i]])[j]
    }
  }

  lapply(Plist, function(P) as.matrix(Matrix::bdiag(P)))
}