#' Control options for gdrm()
#'
#' @param force_intercept Logical. If `TRUE` (default) a global intercept is 
#'   forced to the model if none linear part is present.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information 
#'   if available.
#' @param optimizer Character string specifying the optimization method. 
#'   Options are "adam" (default) and "sa" (simulated annealing).
#' @param optimizer_control List of control parameters for the chosen optimizer.
#'   Use `adam_control()` for ADAM or `sa_control()` for simulated annealing.
#'   If NULL, default controls are used.
#'
#' @returns A list of class `gdrm_control`.
#'

#' @export
gdrm_control <- function(
  force_intercept = TRUE, 
  expected = TRUE,
  optimizer = c("adam", "sa"),
  optimizer_control = NULL
) {
  
  # Match optimizer argument
  optimizer <- match.arg(optimizer)
  
  # Set default optimizer control if not provided
  if (is.null(optimizer_control)) {
    optimizer_control <- switch(
      optimizer,
      "adam" = adam_control(),
      "sa" = sa_control()
    )
  }
  
  r <- list(
    force_intercept = force_intercept,
    expected = expected,
    optimizer = optimizer,
    optimizer_control = optimizer_control
  )
  
  class(r) <- "gdrm_control"
  r
}

#' Generalized distributional regression models
#'
#' Fits generalized distributional regression models, which are models in which 
#' any parameters of a distribution can be modelled as a function (linear, 
#' non-linear or both) of covariates.
#'
#' @param formulae A `formulae` object, which consists of formulas for different 
#'   parameters of the model separated by `&` operator.
#' @param distrib A `[distrib]` object specifying the distribution assumed for 
#'   the response variable.
#' @param data A `data.frame`.
#' @param gdrm_control_list A list created by `[gdrm_control()]` function.
#' @param coef_bounds A list of coefficient bounds. The list should contain an element for each model parameters, and any of those element should be another list containing named elements being numerical vector of lower and upper bounds for parameters.
#' @returns An object of class `gdrm`.
#'
#' @details
#' The function supports multiple optimization algorithms:
#' \describe{
#'   \item{ADAM}{Adaptive Moment Estimation - efficient for large datasets}
#'   \item{SA}{Simulated Annealing - robust for complex surfaces}
#' }
#' 
#' The optimizer and its control parameters can be specified through the 
#' `gdrm_control_list` parameter.
#'
#' @export
gdrm <- function(
  formulae,
  distrib = normal1(),
  data,
  gdrm_control_list = gdrm_control(),
  coef_bounds = NULL
) {
  
  force_intercept <- gdrm_control_list$force_intercept
  expected <- gdrm_control_list$expected
  optimizer <- gdrm_control_list$optimizer
  optimizer_control <- gdrm_control_list$optimizer_control
  
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
  ncoefs <- sapply(mod_comp, function(f) {
    sum(sapply(f, function(ff) length(ff$par)))
  })
  
  # npar
  npar <- length(ncoefs)
  
  # id of coefficients relative to its model parameter
  coef_id <- rep(1:npar, ncoefs)
  
  # get coef names
  coef_names <- get_coef_names(mod_comp)
  
  # create bound list
  coef_bounds <- make_coef_bounds(coef_bounds, coef_names)
  
  # extract lower bounds
  bound_lower <- lapply(coef_bounds, function(ll) {
    unlist(lapply(ll, function(x) x[1]))
  })
  
  # extract upper bounds
  bound_upper <- lapply(coef_bounds, function(ll) {
    unlist(lapply(ll, function(x) x[2]))
  })
  
  # create map_functions
  map_functions <- Map(
    function(lower, upper) {
      make_map_function(lower, upper)
    },
    lower = bound_lower,
    upper = bound_upper
  )
  
  # set starting values
  mod_comp <- gdrm_start(response, mod_comp, distrib, coef_bounds)
  
  # create map functions for hyperparameters
  nhyper <- lapply(mod_comp, function(comp) {
    sum(sapply(comp, function(build) length(build$hyperpar)))
  })
  
  hyper_bounds <- lapply(nhyper, function(x) {
    rep(list(c(0, Inf)), x)
  })
  
  hyper_bound_lower <- lapply(hyper_bounds, function(ll) {
    unlist(lapply(ll, function(x) x[1]))
  })
  
  hyper_bound_upper <- lapply(hyper_bounds, function(ll) {
    unlist(lapply(ll, function(x) x[2]))
  })
  
  hyper_map_functions <- Map(
    function(lower, upper) {
      make_map_function(lower, upper)
    },
    lower = hyper_bound_lower,
    upper = hyper_bound_upper
  )
  
  # OPTIMIZATION - Track elapsed time
  start_time <- Sys.time()
  
  optimization_result <- switch(
    optimizer,
    "adam" = gdrm_adam(
      response = response,
      distrib = distrib,
      mod_comp = mod_comp,
      map_functions = map_functions,
      adam_control_list = optimizer_control
    ),
    "sa" = gdrm_sa(
      response = response,
      distrib = distrib,
      mod_comp = mod_comp,
      map_functions = map_functions,
      sa_control_list = optimizer_control
    ),
    stop("Unknown optimizer: ", optimizer)
  )
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  # Extract optimized model components
  mod_comp_optimized <- optimization_result$mod_comp
  
  # Build return object
  result <- list(
    mod_comp = mod_comp_optimized,
    formulae = formulae,
    distrib = distrib,
    data = data,
    optimization = list(
      optimizer = optimizer,
      optimizer_control = optimizer_control,
      iterations = optimization_result$it,
      elapsed_time = elapsed_time
    ),
    call = match.call()
  )
  
  class(result) <- "gdrm"
  result
}

#' Print method for gdrm objects
#' 
#' @param x A gdrm object
#' @param ... Additional arguments (ignored)
#' @export
print.gdrm <- function(x, ...) {
  cat("\nGeneralized Distributional Regression Model\n")
  cat("===========================================\n\n")
  
  cat("Distribution:", x$distrib$distrib, "\n")
  cat("Optimizer:", x$optimization$optimizer, "\n")
  cat("Iterations:", x$optimization$iterations, "\n")
  cat("Elapsed time:", format(x$optimization$elapsed_time, digits = 3), "\n")
  
  invisible(x)
}