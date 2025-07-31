#' Gradient of the loglikelihood function wrt theta of gdrm model
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param fitted Optional. An object returned by `[gdrm_fitted()]`.
#' @return A list containing the gradient of the loglikelihood wrt eta for each model parameter.
#' @export
gdrm_l_theta <- function(response, mod_comp, distrib, fitted = NULL) {
  if (is.null(fitted)) {
    fitted <- gdrm_fitted(mod_comp, distrib)
  }
  distrib$grad(response, fitted, sum = FALSE)
}


#' Hessian of the loglikelihood function wrt theta of gdrm model
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#' @param fitted Optional. An object returned by `[gdrm_fitted()]`.
#' @return A list containing the hessian of the loglikelihood wrt eta for each model parameter.
#' @export
gdrm_l2_theta2 <- function(response, mod_comp, distrib, expected = TRUE, fitted = NULL) {
  if (is.null(fitted)) {
    fitted <- gdrm_fitted(mod_comp, distrib)
  }
  distrib$hess(response, fitted, sum = FALSE, expected = expected)
}


#' First derivative of theta wrt eta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param predict Optional. An object returned by `[gdrm_predict()]`.
#' @return A list containing the first derivative of theta wrt eta for each model parameter.
#' @export
gdrm_theta_eta <- function(mod_comp, distrib, predict = NULL) {
  if (is.null(predict)) {
    predict <- gdrm_predict(mod_comp) 
  }
  for(i in 1:length(mod_comp)) {
    predict[[i]] <- distrib$link_list[[i]]$theta.eta(predict[[i]])
  }
  predict
}


#' Second derivative of theta wrt eta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param predict Optional. An object returned by `[gdrm_predict()]`.
#' @return A list containing the second derivative of theta wrt eta for each model parameter.
#' @export
gdrm_theta2_eta2 <- function(mod_comp, distrib, predict = NULL) {
  if (is.null(predict)) {
    predict <- gdrm_predict(mod_comp) 
  }
  for(i in 1:length(mod_comp)) {
    predict[[i]] <- distrib$link_list[[i]]$theta2.eta2(predict[[i]])
  }
  predict
}


#' First derivative of eta wrt beta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the first derivatives of eta wrt beta for each component of the model.
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


#' Second derivative of eta wrt beta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @return A list containing the second derivatives of eta wrt beta for each component of the model.
#' @export
gdrm_eta2_beta2 <- function(mod_comp) {
  eta2_beta2 <- vector("list", length(mod_comp))
  names(eta2_beta2) <- names(mod_comp)

  for (j in 1:length(eta2_beta2)) {
    comp <- mod_comp[[j]]
    for (i in 1:length(comp)) {
      if (inherits(comp[[i]], c("linear", "smooth", "ridge"))) {
        eta2_beta2[[j]][[i]] <- matrix(0, ncol(comp[[i]]$X), ncol(comp[[i]]$X))
      } else if (inherits(comp[[i]], c("nl"))) {
        eta2_beta2[[j]][[i]] <- comp[[i]]$hes(comp[[i]]$par)
      }

      names(eta2_beta2[[j]])[i] <- names(comp)[i]
    }
  }
  eta2_beta2
}


#' Gradient of loglikelihood of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param sum Logical. If `TRUE` (default) the gradient is return, else if `FALSE` the individual contributions to gradient are returned.
#' @param penalty Logical. If `TRUE` (default) the penalized gradient is returned.
#' @param P Optional. An object returned by `[gdrm_penalty()]`.
#' @returns Gradient of loglikelihood function.
#'
#' @export
gdrm_grad <- function(response, distrib, mod_comp, sum = TRUE, penalty = TRUE, P = NULL) {
  predict <- gdrm_predict(mod_comp)
  fitted <- gdrm_fitted(mod_comp, distrib, predict)
  l_theta <- gdrm_l_theta(response, mod_comp, distrib, fitted)

  theta_eta <- gdrm_theta_eta(mod_comp, distrib, predict)

  eta_beta <- gdrm_eta_beta(mod_comp)

  eta_beta <- lapply(eta_beta, function(x) do.call(cbind, x))

  grad_list <- Map(
    function(l_theta, theta_eta, eta_beta) {
      l_theta * theta_eta * eta_beta
    },
    l_theta = l_theta,
    theta_eta = theta_eta,
    eta_beta = eta_beta
  )

  if (penalty) {
    if (is.null(P)) {
      P <- gdrm_penalty(mod_comp)
    }
  }

  if (sum) {
    grad_list <- lapply(grad_list, colSums)
    if (penalty) {
      par <- lapply(gdrm_coef(mod_comp), unlist)
      P <- Map(function(P, par) drop(2*P%*%par), P = P, par = par)
      grad_list <- Map(function(g, P) g - P, g = grad_list, P = P)
    }
  } else {

    if (penalty) {
      n <- length(l_theta$mu)
      par <- lapply(gdrm_coef(mod_comp), unlist)
      P <- Map(function(P, par) drop(2*P%*%par)/n, P = P, par = par)
      grad_list <- Map(function(g_matrix, P_vector) sweep(g_matrix, 2, P_vector, "-"), grad_list, P)
    }
  }

  grad_list
}


#' Hessian of loglikelihood of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param sum Logical. If `TRUE` (default) the hessian is returned, else if `FALSE` the individual contributions to the hessian are returned.
#' @param penalty Logical. If `TRUE` (default) the penalized hessian is returned.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#'
#' @returns Hessian of loglikelihood function.
#'
#' @export
gdrm_hessian <- function(response, distrib, mod_comp, sum = TRUE, penalty = TRUE, expected = TRUE) {
  l_theta <- gdrm_l_theta(response, mod_comp, distrib)
  l2_theta2 <- gdrm_l2_theta2(response, mod_comp, distrib, expected = expected)

  theta_eta <- gdrm_theta_eta(mod_comp, distrib)
  theta2_eta2 <- gdrm_theta2_eta2(mod_comp, distrib)

  eta_beta <- gdrm_eta_beta(mod_comp)
  eta2_beta2 <- gdrm_eta2_beta2(mod_comp)

  n <- length(response)
  p <- length(l_theta)
  k <- lapply(eta_beta, length)

  l <- lapply(eta_beta, function(comp) {
    sapply(comp, ncol)
  })
  nl <- sapply(l, length)

  npar <- sum(unlist(l))
  pars <- unlist(lapply(eta_beta, function(comp) {
    sapply(comp, colnames)
  }))

  hi <- matrix(0, nrow = npar, ncol = npar, dimnames = list(pars, pars))
  h_list <- rep(list(hi), n)

  # Helper function to safely extract eta2_beta2
  get_eta2_beta2 <- function(eta2_list, j, l, i) {
    eta2_array <- eta2_list[[j]][[l]]
    if (length(dim(eta2_array)) == 3) {
      return(eta2_array[,, i])
    } else {
      return(eta2_array)
    }
  }

  # Helper function for parameter indexing
  get_param_indices <- function(j, l, l_structure) {
    if (j == 1 && l == 1) {
      start_idx <- 1
    } else if (j == 1) {
      start_idx <- sum(l_structure[[j]][1:(l - 1)]) + 1
    } else if (l == 1) {
      start_idx <- sum(unlist(l_structure[1:(j - 1)])) + 1
    } else {
      start_idx <- sum(unlist(l_structure[1:(j - 1)])) +
        sum(l_structure[[j]][1:(l - 1)]) +
        1
    }

    end_idx <- sum(unlist(l_structure[1:(j - 1)])) + sum(l_structure[[j]][1:l])
    if (j == 1) {
      end_idx <- sum(l_structure[[j]][1:l])
    }

    return(c(start_idx, end_idx))
  }

  for (i in 1:n) {
    for (j1 in 1:p) {
      for (j2 in j1:p) {
        if (j1 == j2) {
          # Same distribution parameter
          for (l1 in 1:nl[j1]) {
            for (l2 in l1:nl[j1]) {
              if (l1 == l2) {
                # Case 1: same component, same predictor
                eta2_term <- get_eta2_beta2(eta2_beta2, j1, l1, i)

                h_block <- eta_beta[[j1]][[l1]][i, ] %*%
                  t(eta_beta[[j1]][[l1]][i, ]) *
                  (l2_theta2[j1, j1, i] *
                    (theta_eta[[j1]][i])^2 +
                    l_theta[[j1]][i] * theta2_eta2[[j1]][i]) +
                  l_theta[[j1]][i] * theta_eta[[j1]][i] * eta2_term
              } else {
                # Case 2: different components, same predictor
                h_block <- eta_beta[[j1]][[l1]][i, ] %*%
                  t(eta_beta[[j1]][[l2]][i, ]) *
                  (l2_theta2[j1, j1, i] *
                    (theta_eta[[j1]][i])^2 +
                    l_theta[[j1]][i] * theta2_eta2[[j1]][i])
              }

              # Get parameter indices
              indices_l1 <- get_param_indices(j1, l1, l)
              indices_l2 <- get_param_indices(j1, l2, l)

              # Fill the Hessian block
              h_list[[i]][
                indices_l1[1]:indices_l1[2],
                indices_l2[1]:indices_l2[2]
              ] <- h_block
              if (l1 != l2) {
                h_list[[i]][
                  indices_l2[1]:indices_l2[2],
                  indices_l1[1]:indices_l1[2]
                ] <- t(h_block)
              }
            }
          }
        } else {
          # Different distribution parameters (j1 != j2)
          for (l1 in 1:nl[j1]) {
            for (l2 in 1:nl[j2]) {
              # Cases 3 & 4: cross-parameter terms
              h_block <- eta_beta[[j1]][[l1]][i, ] %*%
                t(eta_beta[[j2]][[l2]][i, ]) *
                theta_eta[[j1]][i] *
                theta_eta[[j2]][i] *
                l2_theta2[j1, j2, i]

              # Get parameter indices for cross-parameter blocks
              indices_l1 <- get_param_indices(j1, l1, l)
              indices_l2 <- get_param_indices(j2, l2, l)

              # Fill the Hessian block
              h_list[[i]][
                indices_l1[1]:indices_l1[2],
                indices_l2[1]:indices_l2[2]
              ] <- h_block
              h_list[[i]][
                indices_l2[1]:indices_l2[2],
                indices_l1[1]:indices_l1[2]
              ] <- t(h_block)
            }
          }
        }
      }
    }
  }

  # Sum over observations if requested
  if (sum) {
    h <- Reduce("+", h_list)

    if (penalty) {
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))
      h <- h - 2*P
    }


  } else {

    if (penalty) {
      n <- length(l_theta[[1]])
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))/n
      h_list <- lapply(h_list, function(h) h - 2*P)
    }

    h <- simplify2array(h_list)
  }

  return(h)
}


#' Mapped Gradient of loglikelihood of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#' @param sum Logical. If `TRUE` (default) the gradient is return, else if `FALSE` the individual contributions to gradient are returned.
#' @param penalty Logical. If `TRUE` (default) the penalized gradient is returned.
#'
#' @returns Mapped Gradient of loglikelihood function.
#'
#' @export
gdrm_grad_map <- function(response, distrib, mod_comp, map_functions, sum = TRUE, penalty = TRUE) {

  par <- lapply(gdrm_coef(mod_comp), unlist)

  invert_jacobian <- lapply(map_functions, function(map) map$invert_jacobian)

  j <- Map(function(j, par) j(par), j = invert_jacobian, par = par)
  g <- gdrm_grad(response, distrib, mod_comp, sum, penalty = FALSE)

  if (sum) {
    
    g <- Map(function(g, j) g*j, g = g, j = j)

    if (penalty) {
      P <- Map(function(P, par) drop(2*P%*%par), P = gdrm_penalty(mod_comp), par = par)
      g <- Map(function(g, P) g - P, g = g, P = P)
    }

  } else {
    g <- Map(function(g, j) t(t(g) * j), g, j)

    if (penalty) {
      n <- NROW(g[[1]])
      P <- Map(function(P, par) drop(2*P%*%par)/n, P = gdrm_penalty(mod_comp), par = par)
      g <- Map(function(g_matrix, P_vector) sweep(g_matrix, 2, P_vector, "-"), g, P)
    }
  }

  g
}


#' Mapped Hessian of loglikelihood of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#' @param sum Logical. If `TRUE` (default) the hessian is returned, else if `FALSE` the individual contributions to the hessian are returned.
#' @param penalty Logical. If `TRUE` (default) the penalized hessian is returned.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#'
#' @returns Mapped Hessian of loglikelihood function.
#'
#' @export
gdrm_hessian_map <- function(response, distrib, mod_comp, map_functions, sum = TRUE, penalty = TRUE, expected = TRUE) {

  par <- lapply(gdrm_coef(mod_comp), unlist)

  invert_jacobian <- lapply(map_functions, function(map) map$invert_jacobian)
  invert_hessian <- lapply(map_functions, function(map) map$invert_hessian)

  j <- unlist(Map(function(j, par) j(par), j = invert_jacobian, par = par))
  dj <- diag(j)
  h <- unlist(Map(function(h, par) h(par), h = invert_hessian, par = par))
  
  hess <- gdrm_hessian(response, distrib, mod_comp, sum, penalty = FALSE)

  if (sum) {
    g <- unlist(gdrm_grad(response, distrib, mod_comp, sum))
    hess <- dj%*%hess%*%dj + diag(g*h)

    if (penalty) {
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))
      hess <- hess - 2*P
    }
  } else {

    g <- gdrm_grad(response, distrib, mod_comp, sum)
    g <- do.call(cbind, g)
    k <- ncol(g)
    n <- nrow(g)
    result <- array(dim = c(k, k, n))
    for(i in 1:n) {
      result[,,i] <- dj %*% hess[,,i] %*% dj + diag(g[i,] * h)
    }
    hess <- result

    if (penalty) {
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))/n
      hess - array(2*P, dim = c(k, k, n))
    }
  }
  hess
}


#' Loglikelihood function of gdrm model
#'
#' @param response Response variable.
#' @param distrib A `[distrib]` object.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param penalty Logical. If `TRUE` (default) the penalty term is added to the loglikelihood.
#' @param fitted Optional. An object returned by `[gdrm_fitted()]`.
#' @param predict Optional. An object returned by `[gdrm_predict()]`.
#' @param P Optional. An object returned by `[gdrm_penalty()]`.
#' 
#' @returns A numeric value which is the loglikelihood of the gdrm model.
#'
#' @export
gdrm_loglik <- function(response, distrib, mod_comp, penalty = TRUE, fitted = NULL, predict = NULL, P = NULL) {
  ll <- distrib$loglik

  if (is.null(fitted)) {
    if (is.null(predict)) {
      predict <- gdrm_predict(mod_comp)
    }
    fitted <- gdrm_fitted(mod_comp, distrib, predict)
  } 
    
  if (penalty) {

    if (is.null(P)) {
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp))) 
    } else {
      P <- as.matrix(Matrix::bdiag(P))
    }
    
    par <- gdrm_coef_vector(mod_comp)
    pen <- drop(t(par)%*%P%*%par)
  } else {
    pen <- 0
  }

  ll(response, fitted) - pen
}