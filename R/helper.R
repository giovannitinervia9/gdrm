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


