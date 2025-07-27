#' Control options for gdrm()
#' @param force_intercept Logical. If `TRUE` (default) a global intercept is forced to the model if none linear part is present.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#' @returns A list of class `gdrm_control`.
#'
#' @export
gdrm_control <- function(force_intercept = TRUE, expected = TRUE) {
  r <- list(
    force_intercept = force_intercept,
    expected = expected)
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

#' Update parameters of a gdrm model
#'
#' @param par A numeric vector of updated coefficients of a gdrm model.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A new `mod_comp` object with updated parameters.
#'
#' @export
gdrm_update_coef <- function(par, mod_comp) {
  # Create a deep copy of mod_comp to avoid modifying the original
  updated_mod <- mod_comp
  
  # Initialize parameter index counter
  par_idx <- 1
  
  # Function to get the total number of parameters needed
  count_parameters <- function(mod_comp) {
    total_params <- 0
    
    for (i in 1:length(mod_comp)) {
      comp <- mod_comp[[i]]
      
      for (j in 1:length(comp)) {
        build <- comp[[j]]
        # Each build always has a 'par' element
        total_params <- total_params + length(build$par)
      }
    }
    
    return(total_params)
  }
  
  # Check if the parameter vector length matches expected
  expected_params <- count_parameters(mod_comp)
  if (length(par) != expected_params) {
    stop(paste("Parameter vector length (", length(par), 
               ") does not match expected number of parameters (", 
               expected_params, ")", sep = ""))
  }
  
  # Update parameters - simplified since each component always has 'par'
  for (i in 1:length(updated_mod)) {
    comp <- updated_mod[[i]]
    
    for (j in 1:length(comp)) {
      build <- comp[[j]]
      
      # Each build always has a 'par' element
      n_params <- length(build$par)
      if (n_params > 0) {
        # Extract the corresponding slice from par
        new_params <- par[par_idx:(par_idx + n_params - 1)]
        
        # Update the parameters (keep names if they exist)
        if (!is.null(names(build$par))) {
          names(new_params) <- names(build$par)
        }
        
        updated_mod[[i]][[j]]$par <- new_params
        par_idx <- par_idx + n_params
      }
    }
  }
  
  return(updated_mod)
}


#' Map coefficients of the model to unconstrained space
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#'
#' @returns A vector containing mapped coefficients of the gdrm model.
#'
#' @export
gdrm_map_coef <- function(mod_comp, map_functions) {
  par <- lapply(gdrm_coef(mod_comp), unlist)
  map <- lapply(map_functions, function(map) map$map)
  unlist(Map(function(map, par) map(par), map, par))
}


#' Adds mapped coefficients to a gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#'
#' @returns The `mod_comp` passed as input with added the `par_map` object to each component of the model parameters.
#'
#' @export
gdrm_add_map_coef <- function(mod_comp, map_functions) {
  # Create a deep copy of mod_comp to avoid modifying the original
  updated_mod <- mod_comp
  par <- gdrm_map_coef(updated_mod, map_functions)
  
  # Initialize parameter index counter
  par_idx <- 1
  
  # Function to get the total number of parameters needed
  count_parameters <- function(mod_comp) {
    total_params <- 0
    
    for (i in 1:length(mod_comp)) {
      comp <- mod_comp[[i]]
      
      for (j in 1:length(comp)) {
        build <- comp[[j]]
        # Each build always has a 'par' element
        total_params <- total_params + length(build$par)
      }
    }
    
    return(total_params)
  }
  
  # Check if the parameter vector length matches expected
  expected_params <- count_parameters(mod_comp)
  if (length(par) != expected_params) {
    stop(paste("Parameter vector length (", length(par), 
               ") does not match expected number of parameters (", 
               expected_params, ")", sep = ""))
  }
  
  # Update parameters - simplified since each component always has 'par'
  for (i in 1:length(updated_mod)) {
    comp <- updated_mod[[i]]
    
    for (j in 1:length(comp)) {
      build <- comp[[j]]
      
      # Each build always has a 'par' element
      n_params <- length(build$par)
      if (n_params > 0) {
        # Extract the corresponding slice from par
        new_params <- par[par_idx:(par_idx + n_params - 1)]
        
        # Update the parameters (keep names if they exist)
        if (!is.null(names(build$par))) {
          names(new_params) <- names(build$par)
        }
        
        updated_mod[[i]][[j]]$par_map <- new_params
        par_idx <- par_idx + n_params
      }
    }
  }
  
  return(updated_mod)
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


#' Update hyperparameters of a gdrm model
#'
#' @param par A numeric vector of updated hyperparameters of a gdrm model.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A new `mod_comp` object with updated hyperparameters.
#'
#' @export
gdrm_update_hyperpar <- function(par, mod_comp) {
  # Create a deep copy of mod_comp to avoid modifying the original
  updated_mod <- mod_comp
  
  # Initialize parameter index counter
  par_idx <- 1
  
  # Helper function to recursively update hyperparameters
  update_hyperpar_recursive <- function(hyperpar_structure, par_vector, start_idx) {
    if (is.numeric(hyperpar_structure) && length(hyperpar_structure) == 1) {
      # Simple scalar hyperparameter
      return(list(value = par_vector[start_idx], next_idx = start_idx + 1))
    } else if (is.list(hyperpar_structure)) {
      # Nested list structure
      updated_structure <- hyperpar_structure
      current_idx <- start_idx
      
      for (i in 1:length(hyperpar_structure)) {
        result <- update_hyperpar_recursive(hyperpar_structure[[i]], par_vector, current_idx)
        updated_structure[[i]] <- result$value
        current_idx <- result$next_idx
      }
      
      return(list(value = updated_structure, next_idx = current_idx))
    } else {
      # Handle other cases (shouldn't occur with proper hyperpar structure)
      return(list(value = hyperpar_structure, next_idx = start_idx))
    }
  }
  
  # Function to get the total number of hyperparameters needed
  count_hyperparameters <- function(mod_comp) {
    total_params <- 0
    
    for (i in 1:length(mod_comp)) {
      comp <- mod_comp[[i]]
      
      for (j in 1:length(comp)) {
        build <- comp[[j]]
        # Each build always has a 'hyperpar' element
        total_params <- total_params + length(unlist(build$hyperpar))
      }
    }
    
    return(total_params)
  }
  
  # Check if the parameter vector length matches expected
  expected_params <- count_hyperparameters(mod_comp)
  if (length(par) != expected_params) {
    stop(paste("Parameter vector length (", length(par), 
               ") does not match expected number of hyperparameters (", 
               expected_params, ")", sep = ""))
  }
  
  # Update hyperparameters following the same order as gdrm_hyperpar_vector
  for (i in 1:length(updated_mod)) {
    comp <- updated_mod[[i]]
    
    for (j in 1:length(comp)) {
      build <- comp[[j]]
      
      # Update hyperparameters using recursive function
      result <- update_hyperpar_recursive(build$hyperpar, par, par_idx)
      updated_mod[[i]][[j]]$hyperpar <- result$value
      par_idx <- result$next_idx
    }
  }
  
  return(updated_mod)
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


#' Gradient of the loglikelihood function wrt theta of gdrm model
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list containing the gradient of the loglikelihood wrt eta for each model parameter.
#' @export
gdrm_l_theta <- function(response, mod_comp, distrib) {
  distrib$grad(response, gdrm_fitted(mod_comp, distrib), sum = FALSE)
}


#' Hessian of the loglikelihood function wrt theta of gdrm model
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#' @return A list containing the hessian of the loglikelihood wrt eta for each model parameter.
#' @export
gdrm_l2_theta2 <- function(response, mod_comp, distrib, expected = TRUE) {
  distrib$hess(response, gdrm_fitted(mod_comp, distrib), sum = FALSE, expected = expected)
}


#' First derivative of theta wrt eta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list containing the first derivative of theta wrt eta for each model parameter.
#' @export
gdrm_theta_eta <- function(mod_comp, distrib) {
  Map(
    function(link, eta) {
      link$theta.eta(eta)
    },
    link = distrib$link_list,
    eta = gdrm_predict(mod_comp)
  )
}


#' Second derivative of theta wrt eta of gdrm model
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object.
#' @return A list containing the second derivative of theta wrt eta for each model parameter.
#' @export
gdrm_theta2_eta2 <- function(mod_comp, distrib) {
  Map(
    function(link, eta) {
      link$theta2.eta2(eta)
    },
    link = distrib$link_list,
    eta = gdrm_predict(mod_comp)
  )
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
#'
#' @returns Gradient of loglikelihood function.
#'
#' @export
gdrm_grad <- function(response, distrib, mod_comp, sum = TRUE, penalty = TRUE) {
  l_theta <- gdrm_l_theta(response, mod_comp, distrib)

  theta_eta <- gdrm_theta_eta(mod_comp, distrib)

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

  if (sum) {
    grad_list <- lapply(grad_list, colSums)
    if (penalty) {
      par <- lapply(gdrm_coef(mod_comp), unlist)
      P <- Map(function(P, par) drop(2*P%*%par), P = gdrm_penalty(mod_comp), par = par)
      grad_list <- Map(function(g, P) g + P, g = grad_list, P = P)
    }
  } else {

    if (penalty) {
      n <- length(l_theta$mu)
      par <- lapply(gdrm_coef(mod_comp), unlist)
      P <- Map(function(P, par) drop(2*P%*%par)/n, P = gdrm_penalty(mod_comp), par = par)
      grad_list <- Map(function(g_matrix, P_vector) sweep(g_matrix, 2, P_vector, "+"), grad_list, P)
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
      h <- h + 2*P
    }


  } else {

    if (penalty) {
      n <- length(l_theta[[1]])
      P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))/n
      h_list <- lapply(h_list, function(h) h + 2*P)
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
      g <- Map(function(g, P) g + P, g = g, P = P)
    }

  } else {
    g <- Map(function(g, j) t(t(g) * j), g, j)

    if (penalty) {
      n <- NROW(g[[1]])
      P <- Map(function(P, par) drop(2*P%*%par)/n, P = gdrm_penalty(mod_comp), par = par)
      g <- Map(function(g_matrix, P_vector) sweep(g_matrix, 2, P_vector, "+"), g, P)
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
      hess <- hess + 2*P
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
      hess + array(2*P, dim = c(k, k, n))
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
#'
#' @returns A numeric value which is the loglikelihood of the gdrm model.
#'
#' @export
gdrm_loglik <- function(response, distrib, mod_comp, penalty = TRUE) {
  ll <- distrib$loglik
  
  theta <- gdrm_fitted(mod_comp, distrib)
  
  if (penalty) {
    P <- as.matrix(Matrix::bdiag(gdrm_penalty(mod_comp)))
    par <- gdrm_coef_vector(mod_comp)
    pen <- drop(t(par)%*%P%*%par)
  } else {
    pen <- 0
  }

  ll(response, theta) + pen
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
  expected <- gdrm_control_list$expected

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
