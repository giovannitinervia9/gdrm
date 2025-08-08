#' `distrib` objects
#'
#' Objects of class **`distrib`** store all analytic ingredients of a
#' probability distribution (mean, variance, density, random generator,
#' links, &c.) in a single list.
#'
#' @section Components:
#' \describe{
#'   \item{`distrib`}{`character(1)`. Name of the distribution.}
#'   \item{`link_list`}{`list` of `length(parameters)` element containing `link-gdrm` objects created with [gdrm::make_link()].}
#'   \item{`parameters`}{`character` of `length(parameters)` containing the names of model parameters.}
#'   \item{`parameters_description`}{`character` of `length(parameters)` containing the description of model parameters.}
#'   \item{`parameters_bounds`}{`list` of `length(parameters)` element. Each element contains a `numeric(2)` vector specifying lower
#' and upper bounds for the relative parameter.}
#'   \item{`E`, `V`, `A`, `K`, `Me`, `Mo`}{Functions returning mean, variance, skewness,
#'         kurtosis, median and mode. Those functions take as argument a list of `length(parameters)` element `par`, in which each
#' element is a numeric `n`-dimensional vector containing the values of the relative parameter.}
#'   \item{`pdf`}{`function` for computing the probability density function, which takes three arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the pdf is to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements;
#'   \item `log`: logical (default = `FALSE`) specifying if the log of the pdf is to be returned.
#' }}
#'   \item{`loglik`, `lli`}{Functions to compute the loglikelihood of the model and the individual contributions to loglikelihood.
#' Those functions take two arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the likelihood is to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements.
#' }}
#'   \item{`grad`}{Functions to compute the gradient of the loglikelihood with respect to the model parameters, with three arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the gradient is to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those must be single element vectors;
#'   \item `sum`: logical indicating whether the sum over all observation should be returned (default, `TRUE`) or if the individual contributes to the gradient should be returned (`FALSE`).
#' }}
#'   \item{`hess`}{Functions to compute the hessian of the loglikelihood with respect to the model parameters, with three arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the hessian is to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those must be single element vectors;
#'   \item `sum`: logical indicating whether the sum over all observation should be returned (default, `TRUE`) or if the individual contributes to the hessian should be returned (`FALSE`).
#' }}
#'   \item{`cdf`}{`function` for computing the cumulative distribution function, which takes four arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the cdf is to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements;
#'   \item `lower.tail`: logical; if `TRUE` (default) probabilities are \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)};
#'   \item `log.p`: logical (default = `FALSE`) specifying if the log of the cdf is to be returned.
#' }}
#'   \item{`qf`}{`function` for computing the quantile function, which takes four arguments:
#' \itemize{
#'   \item `p`: a numeric vector of `n` elements containing the probability values for which the quantiles are to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements;
#'   \item `lower.tail`: logical; if `TRUE` (default) probabilities are \eqn{P(X \leq x)}, otherwise \eqn{P(X > x)};
#'   \item `log.p`: logical (default = `FALSE`) specifying if probabilities `p` are given as `log(p)`.
#' }}
#'   \item{`rng`}{`function` for random number generations which take two arguments:
#' \itemize{
#'   \item `n`: a numeric vector specifying how many random numbers are to be generated;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements.
#' }}
#'   \item{`mgf`, `cml`, `cf`}{Functions to compute the moment generating function, the cumulant generating function and the characteristic function.
#' Those functions take two arguments:
#' \itemize{
#'   \item `x`: a numeric vector of `n` elements containing the values at which the functions are to be computed;
#'   \item `par`: a list of `length(parameters)` element, one for each parameter of the model; those could be either single element vectors of numeric vectors of `n` elements.
#' }}
#' #' \item{`starting_values`}{List containing starting values for model parameters.}
#' \item{`distrib_args`}{List containing distribution specific arguments.}
#' }
#'
#' @name distrib
#' @docType class
NULL


#' Print method for distrib objects
#'
#' @param x A `distrib` object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.distrib <- function(x, ...) {
  
  # Extract link names
  link_names <- sapply(x$link_list, function(link) link$name)
  
  # Format parameter info
  param_info <- paste0(x$parameters, "(", link_names, ")", collapse = ", ")
  
  cat("\nFamily:", x$distrib, "\n")
  cat("Link functions:", param_info, "\n\n")
  
  invisible(x)
}


#' Compare Analytical and Numerical Derivatives of the Log-Likelihood
#'
#' Evaluates and compares the analytical (explicit) and numerical gradients 
#' and Hessians of the log-likelihood function for a given distribution.
#'
#' @param distrib An object of class [`distrib`].
#' @param x A numeric vector of observations at which the derivatives are evaluated.
#' @param par A list containing the values of the model parameters.
#' @param round_digits An integer specifying the number of decimal places 
#'   to which the difference between analytical and numerical derivatives is rounded.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{\code{gradient}}{A data frame comparing the analytical and numerical gradients.}
#'   \item{\code{hessian}}{A data frame comparing the analytical and numerical Hessians 
#'         (diagonal and upper off-diagonal elements only).}
#' }
#' @examples 
#' check_derivatives(normal1(), rnorm(100), par = list(0, 1))
#' check_derivatives(normal2(), rnorm(100), par = list(0, 1))
#' check_derivatives(normal3(), rnorm(100), par = list(0, 1))
#' @export
#' @importFrom stats setNames
check_derivatives <- function(distrib, x, par, round_digits = 5) {
  par_vec <- unlist(par)

  ll <- function(x, par_vec) {
    par <- as.list(par_vec)
    distrib$loglik(x, par)
  }

  grad_dist <- unlist(distrib$grad(x, par))
  grad_num <- numDeriv::grad(function(par_vec) ll(x, par_vec), par_vec)

  hess_dist <- distrib$hess(x, par, expected = FALSE)
  hess_num <- numDeriv::hessian(function(par_vec) ll(x, par_vec), par_vec)

  # Parameter names
  param_names <- distrib$parameters

  # Gradient comparison
  grad_df <- data.frame(
    param = param_names,
    grad_dist = as.numeric(grad_dist),
    grad_num = as.numeric(grad_num)
  )

  # Function to extract diagonal + upper off-diagonal with names
  get_upper_named <- function(mat, names) {
    idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
    data.frame(
      row = names[idx[, 1]],
      col = names[idx[, 2]],
      value = mat[idx]
    )
  }

  # Hessian comparison
  hess_df <- merge(
    get_upper_named(hess_dist, param_names) |>
      setNames(c("row", "col", "hess_dist")),
    get_upper_named(hess_num, param_names) |>
      setNames(c("row", "col", "hess_num")),
    by = c("row", "col")
  )

  grad_df$diff <- round(grad_df$grad_dist - grad_df$grad_num, round_digits)
  hess_df$diff <- round(hess_df$hess_dist - hess_df$hess_num, round_digits)

  colnames(grad_df) <- c("par", "explicit", "numerical", "difference")
  colnames(hess_df) <- c("rowpar", "colpar", "explicit", "numerical", "difference")
  list(
    gradient = grad_df,
    hessian = hess_df
  )
}


#' Normal distribution (\eqn{\mu, \sigma^2} parametrization)
#'
#' @description
#' Builds a `distrib` object representing a Normal distribution with mean
#' \eqn{\mu} and variance \eqn{\sigma^2}.
#'
#' @param link_mu    Character string naming the link for the mean parameter
#'                   (default `"identity"`).
#' @param link_sigma2 Character string naming the link for the variance parameter
#'                   (default `"log"`).
#'
#' @return An object of class **`distrib`**.
#'   See [distrib] for the structure common to all distribution objects..
#'
#' @seealso [distrib] for the structure common to all distribution objects.
#' @details The `normal1()` function generates a `distrib` object representing a normal distribution with parameters:
#' \deqn{\bm{\theta} = (\mu, \sigma^2)}
#' and pdf:
#' \deqn{f(y) = \dfrac{1}{\sqrt{2 \pi \sigma^2}} \exp \left\{-\dfrac{(y - \mu)^2}{2\sigma^2}\right\}}
#' @export
#' @importFrom stats rnorm var
normal1 <- function(link_mu = "identity", link_sigma2 = "log") {
  distrib <- "normal1"

  parameters <- c("mu", "sigma2")

  parameters_description <- c(mu = "mean", sigma2 = "variance")

  parameters_bounds <- list(
    mu = c(-Inf, Inf),
    sigma2 = c(0, Inf)
  )

  link_ok <- list(
    mu = c("identity", "log", "sqrt"),
    sigma2 = c("log", "sqrt")
  )

  link_list <- list(mu = link_mu, sigma2 = link_sigma2)

  are_link_ok <- mapply(
    function(link, link_ok) link %in% link_ok,
    link = link_list,
    link_ok = link_ok
  )

  if (any(!are_link_ok)) {
    stop(paste0(
      "admissible links are:\n",
      paste(
        paste0(
          parameters,
          ": ",
          sapply(link_ok, function(link) paste0(link, collapse = ", "))
        ),
        sep = "\n",
        collapse = "\n"
      )
    ))
  }

  link_list <- lapply(link_list, gdrm::make_link)

  E <- function(par) {
    par[[1]]
  }

  V <- function(par) {
    par[[2]]
  }

  A <- function(par) {
    0
  }

  K <- function(par) {
    3
  }

  Mo <- function(par) {
    E(par)
  }

  Me <- function(par) {
    E(par)
  }

  pdf <- function(x, par, log = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    dnorm(x, mu, sqrt(sigma2), log = log)
  }

  lli <- function(x, par) {
    pdf(x, par, log = TRUE)
  }

  loglik <- function(x, par) {
    sum(pdf(x, par, log = TRUE))
  }

  grad <- function(x, par, sum = TRUE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]

    grad_mu <- (x - mu) / sigma2

    grad_sigma2 <- .5 * ((x - mu)^2 - sigma2) / sigma2^2

    if (sum) {
      list(mu = sum(grad_mu), sigma2 = sum(grad_sigma2))
    } else {
      list(mu = grad_mu, sigma2 = grad_sigma2)
    }
  }

  hess <- function(x, par, sum = TRUE, expected = TRUE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]

    n <- length(x)
    hess <- array(0, dim = c(2, 2, n))

    if (expected) {
      hess[1, 1, ] <- -1 / sigma2
      hess[2, 2, ] <- - 1 / (2 * sigma2^2)
      hess[1, 2, ] <- hess[2, 1, ] <- 0
    } else {
      hess[1, 1, ] <- -1 / sigma2
      hess[2, 2, ] <- (sigma2 - 2 * (x - mu)^2) / (2 * sigma2^3)
      hess[1, 2, ] <- hess[2, 1, ] <- -(x - mu) / sigma2^2      
    }

    

    if (sum) {
      apply(hess, c(1, 2), sum)
    } else {
      hess
    }
  }

  cdf <- function(x, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    pnorm(x, mu, sqrt(sigma2), lower.tail = lower.tail, log.p = log.p)
  }

  qf <- function(p, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    qnorm(p, mu, sqrt(sigma2), lower.tail = lower.tail, log.p = log.p)
  }

  rng <- function(n = 1, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    rnorm(n, mu, sqrt(sigma2))
  }

  cml <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    mu * x + sigma2 * x^2 / 2
  }

  mgf <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    exp(mu * x + sigma2 * x^2 / 2)
  }

  cf <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    exp((1i) * mu * x - sigma2 * x^2 / 2)
  }

  starting_values <- function(x){
    list(mu = mean(x), sigma2 = var(x))
  }

  r <- list(
    distrib = distrib,
    link_list = link_list,
    parameters = parameters,
    parameters_description = parameters_description,
    parameters_bounds = parameters_bounds,
    E = E,
    V = V,
    A = A,
    K = K,
    Me = Me,
    Mo = Mo,
    pdf = pdf,
    loglik = loglik,
    lli = lli,
    grad = grad,
    hess = hess,
    cdf = cdf,
    qf = qf,
    rng = rng,
    mgf = mgf,
    cml = cml,
    cf = cf,
    starting_values = starting_values,
    distrib_args = NULL
  )

  class(r) <- "distrib"
  r
}


#' Normal distribution (\eqn{\mu, \sigma} parametrization)
#'
#' @description
#' Builds a `distrib` object representing a Normal distribution with mean
#' \eqn{\mu} and standard deviation \eqn{\sigma}.
#'
#' @param link_mu    Character string naming the link for the mean parameter
#'                   (default `"identity"`).
#' @param link_sigma Character string naming the link for the standard deviation parameter
#'                   (default `"log"`).
#'
#' @return An object of class **`distrib`**.
#'   See [distrib] for the structure common to all distribution objects..
#'
#' @seealso [distrib] for the structure common to all distribution objects.
#' @details The `normal2()` function generates a `distrib` object representing a normal distribution with parameters:
#' \deqn{\bm{\theta} = (\mu, \sigma)}
#' and pdf:
#' \deqn{f(y) = \dfrac{1}{\sqrt{2 \pi} \sigma} \exp \left\{-\dfrac{1}{2}\left(\dfrac{y - \mu}{\sigma}\right)^2\right\}}
#' @export
#' @importFrom stats rnorm sd
normal2 <- function(link_mu = "identity", link_sigma = "log") {
  distrib <- "normal2"

  parameters <- c("mu", "sigma")

  parameters_description <- c(mu = "mean", sigma = "std.dev")

  parameters_bounds <- list(
    mu = c(-Inf, Inf),
    sigma = c(0, Inf)
  )

  link_ok <- list(
    mu = c("identity", "log", "sqrt"),
    sigma = c("log", "sqrt")
  )

  link_list <- list(mu = link_mu, sigma = link_sigma)

  are_link_ok <- mapply(
    function(link, link_ok) link %in% link_ok,
    link = link_list,
    link_ok = link_ok
  )

  if (any(!are_link_ok)) {
    stop(paste0(
      "admissible links are:\n",
      paste(
        paste0(
          parameters,
          ": ",
          sapply(link_ok, function(link) paste0(link, collapse = ", "))
        ),
        sep = "\n",
        collapse = "\n"
      )
    ))
  }

  link_list <- lapply(link_list, gdrm::make_link)

  E <- function(par) {
    par[[1]]
  }

  V <- function(par) {
    par[[2]]^2
  }

  A <- function(par) {
    0
  }

  K <- function(par) {
    3
  }

  Mo <- function(par) {
    E(par)
  }

  Me <- function(par) {
    E(par)
  }

  pdf <- function(x, par, log = FALSE) {
    mu <- par[[1]]
    sigma <- par[[2]]
    dnorm(x, mu, sigma, log = log)
  }

  lli <- function(x, par) {
    pdf(x, par, log = TRUE)
  }

  loglik <- function(x, par) {
    sum(pdf(x, par, log = TRUE))
  }

  grad <- function(x, par, sum = TRUE) {
    mu <- par[[1]]
    sigma <- par[[2]]

    grad_mu <- (x - mu) / sigma^2

    grad_sigma <- ((x - mu)^2 - sigma^2) / sigma^3

    if (sum) {
      list(mu = sum(grad_mu), sigma = sum(grad_sigma))
    } else {
      list(mu = grad_mu, sigma = grad_sigma)
    }
  }

  hess <- function(x, par, sum = TRUE, expected = TRUE) {
    mu <- par[[1]]
    sigma <- par[[2]]

    n <- length(x)
    hess <- array(0, dim = c(2, 2, n))

    if (expected) {
      hess[1, 1, ] <- -1 / sigma^2
      hess[2, 2, ] <- - 2 / (sigma^2)
      hess[1, 2, ] <- hess[2, 1, ] <- 0
    } else {
      hess[1, 1, ] <- -1 / sigma^2
      hess[2, 2, ] <- (sigma^2 - 3 * (x - mu)^2) / (sigma^4)
      hess[1, 2, ] <- hess[2, 1, ] <- -2 * (x - mu) / sigma^3

    }

    

    if (sum) {
      apply(hess, c(1, 2), sum)
    } else {
      hess
    }
  }

  cdf <- function(x, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma <- par[[2]]
    pnorm(x, mu, sigma, lower.tail = lower.tail, log.p = log.p)
  }

  qf <- function(p, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma <- par[[2]]
    qnorm(p, mu, sigma, lower.tail = lower.tail, log.p = log.p)
  }

  rng <- function(n = 1, par) {
    mu <- par[[1]]
    sigma <- par[[2]]
    rnorm(n, mu, sigma)
  }

  cml <- function(x, par) {
    mu <- par[[1]]
    sigma <- par[[2]]
    mu * x + sigma^2 * x^2 / 2
  }

  mgf <- function(x, par) {
    mu <- par[[1]]
    sigma <- par[[2]]
    exp(mu * x + sigma^2 * x^2 / 2)
  }

  cf <- function(x, par) {
    mu <- par[[1]]
    sigma <- par[[2]]
    exp((1i) * mu * x - sigma^2 * x^2 / 2)
  }

  starting_values <- function(x){
    list(mu = mean(x), sigma = sd(x))
  }

  r <- list(
    distrib = distrib,
    link_list = link_list,
    parameters = parameters,
    parameters_description = parameters_description,
    parameters_bounds = parameters_bounds,
    E = E,
    V = V,
    A = A,
    K = K,
    Me = Me,
    Mo = Mo,
    pdf = pdf,
    loglik = loglik,
    lli = lli,
    grad = grad,
    hess = hess,
    cdf = cdf,
    qf = qf,
    rng = rng,
    mgf = mgf,
    cml = cml,
    cf = cf,
    starting_values = starting_values,
    distrib_args = NULL
  )

  class(r) <- "distrib"
  r
}


#' Normal distribution (\eqn{\mu, \tau^2} parametrization)
#'
#' @description
#' Builds a `distrib` object representing a Normal distribution with mean
#' \eqn{\mu} and precision \eqn{\tau^2}.
#'
#' @param link_mu    Character string naming the link for the mean parameter
#'                   (default `"identity"`).
#' @param link_tau2 Character string naming the link for the precision parameter
#'                   (default `"log"`).
#'
#' @return An object of class **`distrib`**.
#'   See [distrib] for the structure common to all distribution objects..
#'
#' @seealso [distrib] for the structure common to all distribution objects.
#' @details The `normal3()` function generates a `distrib` object representing a normal distribution with parameters:
#' \deqn{\bm{\theta} = (\mu, \tau^2)}
#' and pdf:
#' \deqn{f(y) = \sqrt{\dfrac{\tau^2}{2 \pi}} \exp \left\{-\dfrac{\tau^2 (y - \mu)^2}{2}\right\}}
#' @export
#' @importFrom stats rnorm var
normal3 <- function(link_mu = "identity", link_tau2 = "log") {
  distrib <- "normal3"

  parameters <- c("mu", "tau2")

  parameters_description <- c(mu = "mean", tau2 = "precision")

  parameters_bounds <- list(
    mu = c(-Inf, Inf),
    tau2 = c(0, Inf)
  )

  link_ok <- list(
    mu = c("identity", "log", "sqrt"),
    tau2 = c("log", "sqrt")
  )

  link_list <- list(mu = link_mu, tau2 = link_tau2)

  are_link_ok <- mapply(
    function(link, link_ok) link %in% link_ok,
    link = link_list,
    link_ok = link_ok
  )

  if (any(!are_link_ok)) {
    stop(paste0(
      "admissible links are:\n",
      paste(
        paste0(
          parameters,
          ": ",
          sapply(link_ok, function(link) paste0(link, collapse = ", "))
        ),
        sep = "\n",
        collapse = "\n"
      )
    ))
  }

  link_list <- lapply(link_list, gdrm::make_link)

  E <- function(par) {
    par[[1]]
  }

  V <- function(par) {
    1 / par[[2]]
  }

  A <- function(par) {
    0
  }

  K <- function(par) {
    3
  }

  Mo <- function(par) {
    E(par)
  }

  Me <- function(par) {
    E(par)
  }

  pdf <- function(x, par, log = FALSE) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    dnorm(x, mu, sqrt(1 / tau2), log = log)
  }

  lli <- function(x, par) {
    pdf(x, par, log = TRUE)
  }

  loglik <- function(x, par) {
    sum(pdf(x, par, log = TRUE))
  }

  grad <- function(x, par, sum = TRUE) {
    mu <- par[[1]]
    tau2 <- par[[2]]

    grad_mu <- tau2 * (x - mu)

    grad_tau2 <- .5 * (1 / tau2 - (x - mu)^2)

    if (sum) {
      list(mu = sum(grad_mu), tau2 = sum(grad_tau2))
    } else {
      list(mu = grad_mu, tau2 = grad_tau2)
    }
  }

  hess <- function(x, par, sum = TRUE, expected = TRUE) {
    mu <- par[[1]]
    tau2 <- par[[2]]

    n <- length(x)
    hess <- array(0, dim = c(2, 2, n))

    if (expected) {
      hess[1, 1, ] <- -tau2
      hess[2, 2, ] <- -1 / (2 * tau2^2)
      hess[1, 2, ] <- hess[2, 1, ] <- 0

    } else {
      hess[1, 1, ] <- -tau2
      hess[2, 2, ] <- -1 / (2 * tau2^2)
      hess[1, 2, ] <- hess[2, 1, ] <- x - mu
    }

    

    if (sum) {
      apply(hess, c(1, 2), sum)
    } else {
      hess
    }
  }

  cdf <- function(x, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    pnorm(x, mu, sqrt(1 / tau2), lower.tail = lower.tail, log.p = log.p)
  }

  qf <- function(p, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    qnorm(p, mu, sqrt(1 / tau2), lower.tail = lower.tail, log.p = log.p)
  }

  rng <- function(n = 1, par) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    rnorm(n, mu, sqrt(1 / tau2))
  }

  cml <- function(x, par) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    mu * x + (1 / tau2) * x^2 / 2
  }

  mgf <- function(x, par) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    exp(mu * x + (1 / tau2) * x^2 / 2)
  }

  cf <- function(x, par) {
    mu <- par[[1]]
    tau2 <- par[[2]]
    exp((1i) * mu * x - (1 / tau2) * x^2 / 2)
  }

  starting_values <- function(x){
    list(mu = mean(x), tau2 = 1/var(x))
  }

  r <- list(
    distrib = distrib,
    link_list = link_list,
    parameters = parameters,
    parameters_description = parameters_description,
    parameters_bounds = parameters_bounds,
    E = E,
    V = V,
    A = A,
    K = K,
    Me = Me,
    Mo = Mo,
    pdf = pdf,
    loglik = loglik,
    lli = lli,
    grad = grad,
    hess = hess,
    cdf = cdf,
    qf = qf,
    rng = rng,
    mgf = mgf,
    cml = cml,
    cf = cf,
    starting_values = starting_values,
    distrib_args = NULL
  )

  class(r) <- "distrib"
  r
}


#' Gamma distribution (\eqn{\mu, \sigma^2} parametrization)
#'
#' @description
#' Builds a `distrib` object representing a Gamma distribution with mean
#' \eqn{\mu} and variance \eqn{\sigma^2}.
#'
#' @param link_mu    Character string naming the link for the mean parameter
#'                   (default `"log"`).
#' @param link_sigma2 Character string naming the link for the variance parameter
#'                   (default `"log"`).
#'
#' @return An object of class **`distrib`**.
#'   See [distrib] for the structure common to all distribution objects..
#'
#' @seealso [distrib] for the structure common to all distribution objects.
#' @details The `gamma1()` function generates a `distrib` object representing a Gamma distribution with parameters:
#' \deqn{\bm{\theta} = (\mu, \sigma^2)}
#' and pdf:
#' \deqn{f(y) = \dfrac{\left(\dfrac{\mu}{\sigma^2}\right)^{\frac{\mu^2}{\sigma^2}} y ^ {\frac{\mu^2}{\sigma^2}-1} \exp\left\{ - \dfrac{\mu}{\sigma^2}y\right\}}{\Gamma\left(\dfrac{\mu^2}{\sigma^2}\right)}}
#' @export
#' @importFrom stats qgamma pgamma dgamma rgamma var
gamma1 <- function(link_mu = "log", link_sigma2 = "log") {
  distrib <- "gamma1"

  parameters <- c("mu", "sigma2")

  parameters_description <- c(mu = "mean", sigma2 = "variance")

  parameters_bounds <- list(
    mu = c(0, Inf),
    sigma2 = c(0, Inf)
  )

  link_ok <- list(
    mu = c("identity", "log", "sqrt"),
    sigma2 = c("log", "sqrt")
  )

  link_list <- list(mu = link_mu, sigma2 = link_sigma2)

  are_link_ok <- mapply(
    function(link, link_ok) link %in% link_ok,
    link = link_list,
    link_ok = link_ok
  )

  if (any(!are_link_ok)) {
    stop(paste0(
      "admissible links are:\n",
      paste(
        paste0(
          parameters,
          ": ",
          sapply(link_ok, function(link) paste0(link, collapse = ", "))
        ),
        sep = "\n",
        collapse = "\n"
      )
    ))
  }

  link_list <- lapply(link_list, gdrm::make_link)

  E <- function(par) {
    par[[1]]
  }

  V <- function(par) {
    par[[2]]
  }

  A <- function(par) {
    2 * sqrt(par[[2]]) / par[[1]]
  }

  K <- function(par) {
    6 * par[[2]] / par[[1]]^2
  }

  Mo <- function(par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    ifelse(
      mu^2 > sigma2,
      (mu^2 - sigma2) / mu,
      0
    )
  }

  Me <- function(par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    qgamma(.5, shape = mu^2 / sigma2, rate = mu / sigma2)
  }

  pdf <- function(x, par, log = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    dgamma(x, shape = mu^2 / sigma2, rate = mu / sigma2, log = log)
  }

  lli <- function(x, par) {
    pdf(x, par, log = TRUE)
  }

  loglik <- function(x, par) {
    sum(pdf(x, par, log = TRUE))
  }

  grad <- function(x, par, sum = TRUE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]

    grad_mu <- (-2 *
      mu *
      digamma(mu^2 / sigma2) +
      2 * mu * log(mu / sigma2) +
      mu +
      2 * mu * log(x) -
      x) /
      sigma2

    grad_sigma2 <- -(mu *
      (-mu *
        digamma(mu^2 / sigma2) +
        mu * log(mu / sigma2) +
        mu +
        mu * log(x) -
        x)) /
      sigma2^2
    if (sum) {
      list(mu = sum(grad_mu), sigma2 = sum(grad_sigma2))
    } else {
      list(mu = grad_mu, sigma2 = grad_sigma2)
    }
  }

  hess <- function(x, par, sum = TRUE, expected = TRUE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]

    n <- length(x)
    hess <- array(0, dim = c(2, 2, n))

    # Common quantities
    k <- mu^2 / sigma2
    psi <- digamma(k) # ψ(k)
    psi1 <- trigamma(k) # ψ₁(k)
    lmu <- log(mu / sigma2)

    if (expected) {
      hess[1, 1, ] <- 3 / sigma2 - 4 * mu^2 * psi1 / sigma2^2
      hess[1, 2, ] <- 2 * mu^3 * psi1 / sigma2^3 - 2 * mu / sigma2^2
      hess[2, 1, ] <- hess[2, 2, ] <- mu^2 / sigma2^3 - mu^4 * psi1 / sigma2^4
    } else {
      logx <- log(x)

      hess[1, 1, ] <- (-(4 * mu^2 * psi1) /
        sigma2 -
        2 * psi +
        2 * lmu +
        2 * logx +
        3) /
        sigma2

      hess[2, 2, ] <- -(mu *
        (2 *
          mu *
          sigma2 *
          psi +
          mu^3 * psi1 +
          sigma2 * (-2 * mu * lmu - 3 * mu - 2 * mu * logx + 2 * x))) /
        sigma2^4

      hess[1, 2, ] <- hess[2, 1, ] <-
        (2 *
          mu *
          sigma2 *
          psi +
          2 * mu^3 * psi1 +
          sigma2 * (-2 * mu * lmu - 3 * mu - 2 * mu * logx + x)) /
        sigma2^3
    }

    if (sum) {
      apply(hess, c(1, 2), sum)
    } else {
      hess
    }
  }

  cdf <- function(x, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    pgamma(
      x,
      shape = mu^2 / sigma2,
      rate = mu / sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  qf <- function(p, par, lower.tail = TRUE, log.p = FALSE) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    qgamma(
      p,
      shape = mu^2 / sigma2,
      rate = mu / sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  rng <- function(n = 1, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    rgamma(n, shape = mu^2 / sigma2, rate = mu / sigma2)
  }

  cml <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    ifelse(x < mu / sigma2, (-mu^2 / sigma2) * log(1 - x * sigma2 / mu), NA)
  }

  mgf <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    ifelse(x < mu / sigma2, (1 - x * sigma2 / mu)^(-mu^2 / sigma2), NA)
  }

  cf <- function(x, par) {
    mu <- par[[1]]
    sigma2 <- par[[2]]
    (1 - (1i) * x * sigma2 / mu)^(-mu^2 / sigma2)
  }

  starting_values <- function(x) {
    list(mu = mean(x), sigma2 = var(x))
  }

  r <- list(
    distrib = distrib,
    link_list = link_list,
    parameters = parameters,
    parameters_description = parameters_description,
    parameters_bounds = parameters_bounds,
    E = E,
    V = V,
    A = A,
    K = K,
    Me = Me,
    Mo = Mo,
    pdf = pdf,
    loglik = loglik,
    lli = lli,
    grad = grad,
    hess = hess,
    cdf = cdf,
    qf = qf,
    rng = rng,
    mgf = mgf,
    cml = cml,
    cf = cf,
    starting_values = starting_values,
    distrib_args = NULL
  )

  class(r) <- "distrib"
  r
}