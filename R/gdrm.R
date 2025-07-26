#' Control options for gdrm()
#' @param force_intercept Logical. If `TRUE` (default) a global intercept is forced to the model if none linear part is present.
#'
#' @returns A list of class `gdrm_control`.
#'
#' @export
gdrm_control <- function(force_intercept = TRUE) {
  r <- list(force_intercept = force_intercept)
  class(r) <- "gdrm_control"
  r
}


#' Predicted values on link scale
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list of predicted values on link scale.
#' @export
gdrm_predict <- function(mod_comp) {
  lapply(mod_comp, function(mod_par) {
    comps <- lapply(mod_par, function(comp) comp$fitted)
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
  pred <- gdrm_predict(mod_comp)
  link <- distrib$link_list
  Map(function(eta, link) link$theta.eta(eta), eta = pred, link = link)
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


#' Derivatives of the loglikelihood function wrt theta of gdrm model
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list containing the the derivatives of the loglikelihood wrt eta for each model parameter.
#' @export
gdrm_l_theta <- function(response, mod_comp, distrib) {
      distrib$grad(response, gdrm_fitted(mod_comp, distrib), sum = FALSE)
    }


#' Derivatives of theta wrt eta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list containing the the derivatives of theta wrt eta for each model parameter.
#' @export
gdrm_theta_eta <- function(mod_comp, distrib) {
  Map(
    function(link, theta) {
      link$theta.eta(theta)
    },
    link = distrib$link_list,
    theta = gdrm_fitted(mod_comp, distrib)
  )
}


#' Derivatives of eta wrt beta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the the derivatives of eta wrt beta for each component of the model.
#' @export
gdrm_eta_beta <- function(mod_comp) {
  eta_beta <- vector("list", length(mod_comp))
  names(eta_beta) <- names(mod_comp)

  for (j in 1:length(eta_beta)) {
    comp <- mod_comp[[j]]
    for (i in 1:length(comp)) {
      if (inherits(comp[[i]], c("linear", "smooth", "ridge"))) {
        eta_beta[[j]][[i]] <- comp[[i]]$X
      } else if (inherits(comp[[i]], c("nl"))) {
        eta_beta[[j]][[i]] <- comp[[i]]$jac(comp[[i]]$par)
      }

      names(eta_beta[[j]])[i] <- names(comp)[i]
    }
  }
  eta_beta
}


#' Gradient of loglikelihood of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param sum Logical. If `TRUE` (default) the gradient is return, else if `FALSE` the individual contributions to gradient are returned.
#'
#' @returns
#'
#' @export
#' @examples
gdrm_grad <- function(response, distrib, mod_comp, sum = TRUE) {
  
    l_theta <- gdrm_l_theta(response, mod_comp, distrib)
    
    theta_eta <- gdrm_theta_eta(mod_comp, distrib)

    eta_beta <- gdrm_eta_beta(mod_comp)

    eta_beta <- lapply(eta_beta, function(x) do.call(cbind, x))

    grad_list <- Map(function(l_theta, theta_eta, eta_beta) {
      l_theta*theta_eta*eta_beta
    }, l_theta = l_theta, theta_eta = theta_eta, eta_beta = eta_beta)
    
    if (sum) {
      grad_list <- lapply(grad_list, colSums)
    }

    grad_list
}


#' Generalized distributional regression models
#' Fits generalized distributional regression models, which are models in which any parameters of a distribution can be modelled as a function (linear, non-linear or both) of covariates.
#'
#' @param formulae A `formulae` object, which consists of formulas for different parameters of the model separated by `&` operator.
#' @param distrib A `[distrib]` object specifying the distribution assumed for the response variable.
#' @param data A `data.frame`.
#' @param gdrm_control_list A list created by `[gdrm_control()]` function.
#'
#' @returns An object of class `gdrm`.
#'
#' @export
gdrm <- function(
  formulae,
  distrib = normal1(),
  data,
  gdrm_control_list = gdrm_control()
) {

  force_intercept <- gdrm_control_list$force_intercept

  # get response
  response <- get_response(formulae, data)

  # get model components
  mod_comp <- interpret_formulae(
    formulae,
    distrib,
    data,
    force_intercept = force_intercept
  )

  # get npar for each model parameter
  npars <- lapply(mod_comp, function(f) {
    sum(sapply(f, function(ff) length(ff$par)))
  })



  # function to build hessian
}