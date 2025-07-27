#' Control options for Simulated Annealing algorithm
#'
#' @param iter Integer. Number of iterations.
#' @param t_start Numeric. Temperature at the start of algorithm.
#' @param t_end Numeric. Temperature at the end of algorithm.
#' @param sd Numeric. Standard deviation used for sampling new parameters value.
#' @param save_history Logical. If `TRUE` (default) loglikelihood values at each iterations are stored and returned.
#' @param verbose Logical. If `TRUE` (default) progress is printed.
#'
#' @returns A list containing Simulated annealing options.
#' @export
sa_control <- function(iter = 30000, t_start = 100, t_end = 1, sd = 1, save_history = TRUE, verbose = TRUE) {
  list(
    iter = iter,
    t_start = t_start,
    t_end = t_end,
    sd = sd,
    save_history = save_history,
    verbose = verbose
  )
}


#' Simulated annealing optimization for gdrm
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A [`distrib`] object.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#' @param sa_control_list A list of options returned by [`sa_control()`] function.
#'
#' @returns A list containing:
#' \item{`mod_comp`}{`mod_comp` object updated with new parameters value.}
#' \item{`history`}{Numeric vector containing the loglikelihood values at each iterations.}
#'
#' @export
#' @importFrom stats runif
gdrm_sa <- function(
  response,
  mod_comp,
  distrib,
  map_functions,
  sa_control_list = sa_control()) {
  
  # get npar for each model parameter
  ncoefs <- sapply(mod_comp, function(f) {
    sum(sapply(f, function(ff) length(ff$par)))
  })
  
  # number of model parameters
  n_param_groups <- length(ncoefs)
  
  # id of coefficients relative to its model parameter
  coef_id <- rep(1:n_param_groups, ncoefs)
  
  # id of coefficients relative to its model parameter
  coef_id <- rep(1:npar, ncoefs)

  # coef names
  coef_names <- get_coef_names(mod_comp)
  
  # options
  iter <- sa_control_list$iter
  t_start <- sa_control_list$t_start
  t_end <- sa_control_list$t_end
  sd <- sa_control_list$sd
  save_history <- sa_control_list$save_history
  verbose <- sa_control_list$verbose

  # set temperature
  const_temp <- log(t_start/t_end)/iter
  temp <- t_start*exp(-const_temp*0:(iter - 1))

  # set start point
  l0 <- l_best <- gdrm_loglik(response, distrib, mod_comp, penalty = FALSE)
  map_par0 <- map_par_best <- gdrm_map_coef(mod_comp, map_functions)
  npar <- length(map_par0)

  if (save_history) {
    history <- numeric(iter + 1)
    history[1] <- l0
  } else {
    history <- NULL
  }
  
  for (i in 1:iter) {
    
    # sample new par
    map_par1 <- map_par0 + rnorm(npar, 0, sd)

    # invert the par to the original space
    par1 <- gdrm_invert_coef(map_par1, coef_id, map_functions)

    # create a temporal mod_comp1
    mod_comp1 <- gdrm_update_coef(par1, mod_comp)

    # compute loglikelihood with new par
    l1 <- gdrm_loglik(response, distrib, mod_comp1, penalty = FALSE)

    # compute delta (we want it to be greater than 0)
    delta <- l1 - l0

    if (delta > 0 | runif(1) < exp(delta/temp[i])) {
      mod_comp <- mod_comp1
      map_par0 <- map_par1
      l0 <- l1

      if (l0 > l_best) {
        l_best <- l0
        map_par_best <- map_par0
      }
    }

    if (save_history) {
      history[i + 1] <- l0
    }

    if (verbose && i %% 100 == 0) {
      cat(sprintf("Iter %d, LogLik = %.4f\n", i, l0))
    }
  }
  
  # invert the best par to the original space
  par_best <- gdrm_invert_coef(map_par_best, coef_id, map_functions)
  # create a temporal mod_comp1
  mod_comp <- gdrm_update_coef(par_best, mod_comp)

  o <- list(mod_comp = mod_comp, history = history)
  invisible(o)
}