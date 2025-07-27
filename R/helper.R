#' Extract coefficient names
#'
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#'
#' @returns A list containing the coefficient names for each model parameters.
#'
#' @export
get_coef_names <- function(mod_comp) {
  r <- lapply(mod_comp, function(comp){
    lapply(comp, function(build) {
      names(build$par)
    })
  })
  lapply(r, unlist)
}


#' Create list of coefficient bounds
#'
#' @param coef_bounds A list of coefficient bounds. The list should contain an element for each model parameters, and any of those element should be another list containing named elements being numerical vector of lower and upper bounds for parameters.
#' @param coef_names A list of coefficient names created by [`get_coef_names()`] function.
#'
#' @returns A list of coefficient bounds to pass to [`make_map_function()`] function.
#'
#' @export
make_coef_bounds <- function(coef_bounds, coef_names) {
  
  npar <- length(coef_names)
  ncoefs <- sapply(coef_names, length)

  if (is.null(coef_bounds)) {
  coef_bounds <- vector("list", npar)
  names(coef_bounds) <- names(coef_names)
  
  for(i in 1:npar) {
    for(j in 1:ncoefs[[i]]) {
      coef_bounds[[i]][[j]] <- c(-Inf, Inf)
      names(coef_bounds[[i]])[j] <- coef_names[[i]][j]
    }
  }
    coef_bounds
  } else {
    # Initialize the result list
  result <- list()
  component_names <- names(coef_names)
  
  # Loop through each component in coef_names using numerical index
  for (i in 1:length(coef_names)) {
    component <- component_names[i]
    
    # Get coefficient names for this component
    coef_list <- coef_names[[i]]
    
    # Initialize bounds list for this component
    result[[component]] <- list()
    
    # Loop through each coefficient in this component using numerical index
    for (j in 1:length(coef_list)) {
      coef <- coef_list[j]
      
      # Check if bounds exist for this coefficient in this component
      if (component %in% names(coef_bounds) && coef %in% names(coef_bounds[[component]])) {
        # Use existing bounds
        result[[component]][[coef]] <- coef_bounds[[component]][[coef]]
      } else {
        # Use default bounds
        result[[component]][[coef]] <- c(-Inf, Inf)
      }
    }
  }
  
  return(result)
  }
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


#' Starting values of coefficients for gdrm models
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A `[distrib]` object specifying the distribution assumed for the response variable.
#' @param coef_bounds A list of coefficient bounds. The list should contain an element for each model parameters, and any of those element should be another list containing named elements being numerical vector of lower and upper bounds for parameters.
#'
#' @returns Update `mod_comp` setting starting values for the coefficients.
#'
#' @export
gdrm_start <- function(response, mod_comp, distrib, coef_bounds) {
  start <- distrib$starting_values(response)
  link_list <- distrib$link_list
  start <- Map(function(start, link) link$linkfun(start), start = start, link = link_list)
  
  for(i in 1:length(mod_comp)) {
    comp <- mod_comp[[i]]
    k <- 1
    for(j in 1:length(comp)) {
      build <- comp[[j]]
      for(par_idx in 1:length(build$par)) {
        bounds <- coef_bounds[[i]][[k]]
        p <- build$par[par_idx]
        
        # Only adjust if bounds are finite and parameter is out of bounds
        if (!all(is.infinite(bounds))) {
          if (p <= bounds[1]) {
            p <- bounds[1] + 0.1
          }
          if (p >= bounds[2]) {
            p <- bounds[2] - 0.01
          }
        }
        
        build$par[par_idx] <- p
        k <- k + 1
      }
      
      # update intercept with starting value
      if(inherits(build, "linear") && "(Intercept)" %in% names(build$par)) {
        build$par["(Intercept)"] <- start[[i]]
      }
      
      mod_comp[[i]][[j]] <- build
    }
  }
  
  mod_comp
}
