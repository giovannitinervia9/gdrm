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
sa_control <- function(iter = 30000, t_start = 100, t_end = 1, sd = .1, save_history = TRUE, verbose = TRUE) {
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


#' Control Parameters for ADAM Optimization
#'
#' @param maxit Integer. Maximum number of iterations. Default is 100000.
#' @param grad_tol Numeric. Gradient tolerance for convergence. Default is 1e-06.
#' @param param_tol Numeric. Parameter tolerance for convergence. Default is 1e-10.
#' @param require_both Logical. If TRUE, both gradient and parameter tolerances 
#'   must be satisfied for convergence. If FALSE, either condition is sufficient.
#'   Default is FALSE.
#' @param alpha Numeric. Learning rate (step size) for ADAM optimization. 
#'   Default is 0.1.
#' @param beta1 Numeric. Exponential decay rate for the first moment estimates.
#'   Should be in [0, 1). Default is 0.9.
#' @param beta2 Numeric. Exponential decay rate for the second moment estimates.
#'   Should be in [0, 1). Default is 0.999.
#' @param epsilon Numeric. Small constant added to denominator for numerical 
#'   stability. Default is 1e-8.
#' @param verbose Logical. If TRUE, prints iteration progress. Default is TRUE.
#' @param expected Logical. Currently unused parameter for potential future 
#'   extensions. Default is TRUE.
#' @param lr_decay Logical. If TRUE, applies learning rate decay. Default is TRUE.
#' @param lr_decay_it Integer. Iteration interval for learning rate decay.
#'   Default is floor(0.5 * maxit).
#' @param lr_decay_rate Numeric. Rate of learning rate decay. Default is 0.01.
#'
#' @return A named list containing all control parameters for ADAM optimization.
#'
#' @details
#' 
#' The learning rate decay is applied logarithmically when \code{lr_decay = TRUE}:
#' \deqn{\alpha_{new} = \frac{\alpha}{1 + \text{lr\_decay\_rate} \times \log(\text{it})}}

#' @seealso \code{\link{gdrm_adam}} for the main optimization function.
#'
#' @export
adam_control <- function(
  maxit = 100000,
  grad_tol = 1e-06,
  param_tol = 1e-10,
  require_both = FALSE,
  alpha = 0.1,
  beta1 = 0.9,
  beta2 = 0.999,
  epsilon = 1e-8,  
  verbose = TRUE,
  expected = TRUE,
  lr_decay = TRUE,
  lr_decay_it = floor(.5*maxit),
  lr_decay_rate = .01
) {
  list(
    alpha = alpha,
    beta1 = beta1,
    beta2 = beta2,
    epsilon = epsilon,
    grad_tol = grad_tol,
    maxit = maxit,
    verbose = verbose,
    lr_decay = lr_decay,
    lr_decay_it = lr_decay_it,
    lr_decay_rate = lr_decay_rate,
    param_tol = param_tol,
    require_both = require_both
  )
}


#' ADAM Optimization for gdrm models
#'
#' Fits generalized distributional regression models using the ADAM optimization
#' algorithm. This function optimizes model parameters by maximizing the 
#' log-likelihood using adaptive moment estimation.
#'
#' @param response Response variable.
#' @param mod_comp A list of model components from `[interpret_formulae()]`.
#' @param distrib A [`distrib`] object.
#' @param map_functions A list of functions to map parameters from constrained to real line created by [`make_map_function()`].
#' @param adam_control_list A list of options returned by [`adam_control()`] function.
#'
#' @return A list with components:
#'   \item{mod_comp}{Updated model components with optimized parameter values.}
#'   \item{it}{Integer. Number of iterations until convergence.}
#'
#' @details
#' The function implements the ADAM optimization algorithm for maximum likelihood
#' estimation in generalized distributional regression models. ADAM maintains
#' separate adaptive learning rates for each parameter by computing exponential
#' moving averages of gradients and their squared values.
#' 
#' The algorithm updates parameters using bias-corrected moment estimates:
#' \deqn{m_t = \beta_1 m_{t-1} + (1-\beta_1) g_t}
#' \deqn{v_t = \beta_2 v_{t-1} + (1-\beta_2) g_t^2}
#' \deqn{\hat{m}_t = \frac{m_t}{1-\beta_1^t}}
#' \deqn{\hat{v}_t = \frac{v_t}{1-\beta_2^t}}
#' \deqn{\theta_{t+1} = \theta_t + \frac{\alpha \hat{m}_t}{\sqrt{\hat{v}_t} + \epsilon}}
#' 
#' Convergence is determined by either gradient tolerance, parameter change
#' tolerance, or both (depending on \code{require_both} setting).
#'
#' @seealso \code{\link{adam_control}} for control parameter specification.
#'
#' @export
gdrm_adam <- function(
  response,
  distrib,
  mod_comp,
  map_functions,
  adam_control_list = adam_control()
) {
  alpha <- adam_control_list$alpha
  beta1 <- adam_control_list$beta1
  beta2 <- adam_control_list$beta2
  epsilon <- adam_control_list$epsilon
  grad_tol <- adam_control_list$grad_tol
  maxit <- adam_control_list$maxit
  verbose <- adam_control_list$verbose
  lr_decay <- adam_control_list$lr_decay
  lr_decay_it <- adam_control_list$lr_decay_it
  lr_decay_rate <- adam_control_list$lr_decay_rate
  param_tol <- adam_control_list$param_tol
  require_both <- adam_control_list$require_both

  # get npar for each model parameter
  ncoefs <- sapply(mod_comp, function(f) {
    sum(sapply(f, function(ff) length(ff$par)))
  })

  # npar
  npar <- length(ncoefs)
  # id of coefficients relative to its model parameter
  coef_id <- rep(1:npar, ncoefs)

  ncoef <- length(gdrm_coef_vector(mod_comp))

  # Initialize Adam parameters
  m <- rep(0, ncoef) # first moment vector
  v <- rep(0, ncoef) # second moment vector

  it <- 0
  eps <- 2
  param_change <- Inf

  # Penalty matrix
  P <- gdrm_penalty(mod_comp)
  predict <- gdrm_predict(mod_comp)
  fitted <- gdrm_fitted(mod_comp, distrib, predict)
  g_map <- unlist(gdrm_grad_map(response, distrib, mod_comp, map_functions, P = P,
    predict = predict, fitted = fitted))

  while (it < maxit) {
    it <- it + 1
    
    par0 <- gdrm_map_coef(mod_comp, map_functions)

    # Update biased first moment estimate
    m <- beta1 * m + (1 - beta1) * g_map

    # Update biased second raw moment estimate
    v <- beta2 * v + (1 - beta2) * (g_map^2)

    # Compute bias-corrected first moment estimate
    m_hat <- m / max(1 - beta1^it, 1e-12)

    # Compute bias-corrected second raw moment estimate
    v_hat <- v / max(1 - beta2^it, 1e-12)

    # Update parameters (+ for maximization, - for minimization)

    if (lr_decay & it %% lr_decay_it == 0) {
      alpha <- alpha / (1 + lr_decay_rate * log(it))
    }
    
    
    delta <- alpha * m_hat / (sqrt(v_hat) + epsilon)

    par1 <- par0 + delta

    mod_comp <- gdrm_update_coef(
      gdrm_invert_coef(par1, coef_id, map_functions),
      mod_comp
    )
    predict <- gdrm_predict(mod_comp)
    fitted <- gdrm_fitted(mod_comp, distrib, predict)
    g_map <- unlist(gdrm_grad_map(response, distrib, mod_comp, map_functions,
      P = P, predict = predict, fitted = fitted))
  
    # Convergence check
    eps <- max(abs(unlist(g_map)))
    param_change <- max(abs(par1 - par0))
    
    # Combined convergence condition
    if (require_both) {
      converged <- (eps < grad_tol) && (param_change < param_tol)
    } else {
      converged <- (eps < grad_tol) || (param_change < param_tol)
    }
    
    if (converged) break
    
    if (verbose) {

      cat("it: ", it,
      ", max(|grad|) = ", eps,
      ", max(|changepar|) = ", param_change,
      "\n", sep = "")
    }
  }

  invisible(list(mod_comp = mod_comp, it = it))
}
