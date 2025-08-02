#' Predicted values on link scale
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list of predicted values on link scale.
#' @export
gdrm_predict <- function(mod_comp) {
  out <- vector("list", length(mod_comp))
  names(out) <- names(mod_comp)
  for (i in 1:length(mod_comp)) {
    r <- 0
    for (j in 1:length(mod_comp[[i]])) {
      r <- r + mod_comp[[i]][[j]]$fitted(mod_comp[[i]][[j]]$par)
    }
    out[[i]] <- r
  }
  out
}


#' Fitted values on parameter scale
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param predict Optional. An object returned by `[gdrm_predict()]`.
#' @return A list of fitted values on parameter scale.
#' @export
gdrm_fitted <- function(mod_comp, distrib, predict = NULL) {
  if (is.null(predict)) {
    predict <- gdrm_predict(mod_comp)
  }
  out <- vector("list", length(mod_comp))
  names(out) <- names(mod_comp)
  link <- distrib$link_list
  for (i in seq_len(length(predict))) {
    out[[i]] <- link[[i]]$linkinv(predict[[i]])
  }
  out
}


#' Coefficient names of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A character vector containing names of the parameters.
#'
#' @export
gdrm_coef_names <- function(mod_comp) {
  unlist(
    lapply(names(mod_comp), function(comp_name) {
      comp <- mod_comp[[comp_name]]
      lapply(comp, function(build) {
        paste0(comp_name, ".", build$parameters)
      })
    }),
    use.names = FALSE
  )
}


#' Coefficient of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the coefficients for each component of the model.
#' @export
gdrm_coef <- function(mod_comp) {
  lapply(mod_comp, function(comp) {
    lapply(comp, `[[`, "par")
  })
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
    for (j in 1:length(mod_comp[[i]])) {
      coef_vector <- c(coef_vector, mod_comp[[i]][[j]]$par)
    }
  }

  coef_vector
}


#' Hyperparameters of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the hyperparameters for each component of the model.
#' @export
gdrm_hyperpar <- function(mod_comp) {
  lapply(mod_comp, function(comp) {
    lapply(comp, `[[`, "hyperpar")
  })
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
    for (j in 1:length(mod_comp[[i]])) {
      coef_vector <- c(coef_vector, mod_comp[[i]][[j]]$hyperpar)
    }
  }

  coef_vector
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

  for (i in seq_along(mod_comp)) {
    comps <- mod_comp[[i]]
    Pparam <- vector("list", length(comps))

    for (j in seq_along(comps)) {
      comp <- comps[[j]]

      ## --- linear / ridge / nl -------------------------------------------
      if (inherits(comp, c("linear", "ridge", "nl"))) {
        Pmat <- comp$P
        if (hyperpar) {
          Pmat <- comp$hyperpar * Pmat
        }
        Pparam[[j]] <- Pmat

        ## --- smooth ---------------------------------------------------------
      } else {
        Psub <- comp$P
        if (hyperpar) {
          hsub <- comp$hyperpar[[1L]]
          for (k in seq_along(Psub)) {
            for (m in seq_along(Psub[[k]])) {
              Psub[[k]][[m]] <- Psub[[k]][[m]] * hsub[[m]]
            }
          }
        }

        for (k in seq_along(Psub)) {
          Psub[[k]] <- Reduce(`+`, Psub[[k]])
        }

        Pparam[[j]] <- Matrix::bdiag(Psub)
      }
      names(Pparam)[j] <- names(comps)[j]
    }

    Plist[[i]] <- as.matrix(Matrix::bdiag(Pparam))
  }

  Plist
}