#' Create mapping functions between open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' an open interval `(lower, upper)` and the real line, together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, lower, upper)}{Maps from (lower, upper) to the real line: \deqn{f(x) = \log\left(\dfrac{x - \texttt{lower}}{\texttt{upper} - x}\right)}}
#'   \item{invert(x, lower, upper)}{Inverse mapping: \deqn{f^{-1}(z) = \dfrac{\texttt{lower} + \texttt{upper} \cdot e^z}{1 + e^z}}}
#'   \item{map_jacobian(x, lower, upper)}{First derivative of the forward map: \deqn{f'(x) = \dfrac{\texttt{upper} - \texttt{lower}}{(x - \texttt{lower})(\texttt{upper} - x)}}}
#'   \item{map_hessian(x, lower, upper)}{Second derivative of the forward map: \deqn{f''(x) = -\dfrac{(\texttt{lower} - \texttt{upper})(\texttt{lower} + \texttt{upper} - 2x)}{(x - \texttt{lower})^2 (\texttt{upper} - x)^2}}}
#'   \item{invert_jacobian(x, lower, upper)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = \dfrac{e^z(\texttt{upper} - \texttt{lower})}{\left(1 + e^z\right)^2}}}
#'   \item{invert_hessian(x, lower, upper)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = \dfrac{e^z(e^z - 1)(\texttt{lower} - \texttt{upper})}{\left(1 + e^z\right)^3}}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_interval()
#' x <- seq(0.1, 0.9, length.out = 10)
#' f$map(x, lower = 0, upper = 1)
#' f$invert(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
#' f$map_jacobian(x, lower = 0, upper = 1)
#' f$invert_jacobian(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
#' f$map_hessian(x, lower = 0, upper = 1)
#' f$invert_hessian(f$map(x, lower = 0, upper = 1), lower = 0, upper = 1)
map_interval <- function() {

  map <- function(x, lower, upper) {
    log((x - lower) / (upper - x))
  }

  invert <- function(x, lower, upper) {
    (lower + upper * exp(x)) / (1 + exp(x))
  }

  map_jacobian <- function(x, lower, upper) {
    (upper - lower) / ((x - lower) * (upper - x))
  }

  map_hessian <- function(x, lower, upper) {
    ((lower - upper)*(lower + upper - 2*x))/(((x - lower)^2)*((upper - x)^2))
  }

  invert_jacobian <- function(x, lower, upper) {
    (exp(x)*(upper - lower))/(1 + exp(x))^2
  }

  invert_hessian <- function(x, lower, upper) {
    (exp(x)*(exp(x) - 1)*(lower - upper))/(1 + exp(x))^3
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)
}


#' Create mapping functions between half-open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' a half-open interval `(lower, Inf)` and the real line using a logarithmic transformation,
#' together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, lower)}{Maps from (lower, Inf) to the real line: \deqn{f(x) = \log(x - \texttt{lower})}}
#'   \item{invert(x, lower)}{Inverse mapping: \deqn{f^{-1}(z) = \texttt{lower} + \exp(z)}}
#'   \item{map_jacobian(x, lower)}{First derivative of the forward map: \deqn{f'(x) = 1/(x - \texttt{lower})}}
#'   \item{map_hessian(x, lower)}{Second derivative of the forward map: \deqn{f''(x) = -1/(x - \texttt{lower})^2}}
#'   \item{invert_jacobian(x, lower)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = \exp(z)}}
#'   \item{invert_hessian(x, lower)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = \exp(z)}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_positive()
#' x <- 1:5
#' f$map(x, lower = 0)
#' f$invert(f$map(x, lower = 0), lower = 0)
#' f$map_jacobian(x, lower = 0)
#' f$invert_jacobian(f$map(x, lower = 0), lower = 0)
#' f$map_hessian(x, lower = 0)
#' f$invert_hessian(f$map(x, lower = 0), lower = 0)
map_positive <- function() {

  map <- function(x, lower) {
    log(x - lower)
  }

  invert <- function(x, lower) {
    lower + exp(x)
  }

  map_jacobian <- function(x, lower) {
    1/(x - lower)
  }

  map_hessian <- function(x, lower) {
    -1 / (x - lower)^2
  }

  invert_jacobian <- function(x, lower) {
    exp(x)
  }

  invert_hessian <- function(x, lower) {
    exp(x)
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)

}


#' Create mapping functions between half-open intervals and the real line
#'
#' Returns a list of functions implementing a bijective transformation between
#' a half-open interval `(-Inf, upper)` and the real line using a logarithmic transformation,
#' together with their derivatives and inverses.
#'
#' @details
#' The returned list contains:
#' \describe{
#'   \item{map(x, upper)}{Maps from (-Inf, upper) to the real line: \deqn{f(x) = \log(\texttt{upper} - x)}}
#'   \item{invert(x, upper)}{Inverse mapping: \deqn{f^{-1}(z) = \texttt{upper} - \exp(z)}}
#'   \item{map_jacobian(x, upper)}{First derivative of the forward map: \deqn{f'(x) = -1/(\texttt{upper} - x)}}
#'   \item{map_hessian(x, upper)}{Second derivative of the forward map: \deqn{f''(x) = -1/(\texttt{upper} - x)^2}}
#'   \item{invert_jacobian(x, upper)}{First derivative of the inverse map: \deqn{(f^{-1})'(z) = -\exp(z)}}
#'   \item{invert_hessian(x, upper)}{Second derivative of the inverse map: \deqn{(f^{-1})''(z) = -\exp(z)}}
#' }
#'
#' @return A list of six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' f <- map_negative()
#' x <- -(5:1)
#' f$map(x, upper = 0)
#' f$invert(f$map(x, upper = 0), upper = 0)
#' f$map_jacobian(x, upper = 0)
#' f$invert_jacobian(f$map(x, upper = 0), upper = 0)
#' f$map_hessian(x, upper = 0)
#' f$invert_hessian(f$map(x, upper = 0), upper = 0)
map_negative <- function() {

  map <- function(x, upper) {
    log(upper - x)
  }

  invert <- function(x, upper) {
    upper - exp(x)
  }

  map_jacobian <- function(x, upper) {
    -1 / (upper - x)
  }

  map_hessian <- function(x, upper) {
    -1 / (upper - x)^2
  }

  invert_jacobian <- function(x, upper) {
    -exp(x)
  }

  invert_hessian <- function(x, upper) {
    -exp(x)
  }

  list(map = map, invert = invert,
       map_jacobian = map_jacobian, map_hessian = map_hessian,
       invert_jacobian = invert_jacobian, invert_hessian = invert_hessian)

}


#' Create mapping functions between constrained parameter spaces and the real line
#'
#' Returns a set of transformation functions that map parameters from various constrained spaces
#' to the real line and vice versa, along with their first and second derivatives.
#'
#' @param lower A numeric vector containing the lower bounds for each parameter.
#'   Use \code{-Inf} for unbounded below parameters.
#' @param upper A numeric vector containing the upper bounds for each parameter.
#'   Use \code{Inf} for unbounded above parameters.
#'
#' @details
#' The function automatically selects the appropriate transformation for each parameter based on its bounds:
#' \itemize{
#'   \item For parameters bounded both below and above \code{(lower[i], upper[i])}, [map_interval()] is used.
#'   \item For parameters bounded only below \code{(lower[i], Inf)}, [map_positive()] is used.
#'   \item For parameters bounded only above \code{(-Inf, upper[i])}, [map_negative()] is used.
#'   \item For completely unbounded parameters \code{(-Inf, Inf)}, the identity function is used.
#' }
#'
#' The returned list contains six functions:
#' \describe{
#'   \item{map(par)}{Maps parameters from the constrained space to the real line.}
#'   \item{invert(par)}{Maps parameters from the real line back to the constrained space.}
#'   \item{map_jacobian(par)}{Returns the first derivative (Jacobian) of the forward transformation.}
#'   \item{map_hessian(par)}{Returns the second derivative (Hessian) of the forward transformation.}
#'   \item{invert_jacobian(par)}{Returns the first derivative (Jacobian) of the inverse transformation.}
#'   \item{invert_hessian(par)}{Returns the second derivative (Hessian) of the inverse transformation.}
#' }
#' Each function expects a numeric vector \code{par} of the same length as \code{lower} and \code{upper}.
#'
#' @return A list containing six functions: \code{map}, \code{invert}, \code{map_jacobian},
#'   \code{map_hessian}, \code{invert_jacobian}, \code{invert_hessian}.
#'
#' @export
#'
#' @examples
#' # Define constraints for parameters
#' lower <- c(-Inf, 0, -1)
#' upper <- c(Inf, Inf, 1)
#'
#' # Create mapping functions
#' map_functions <- make_map_function(lower, upper)
#'
#' # Extract the individual functions
#' map <- map_functions$map
#' invert <- map_functions$invert
#' map_jacobian <- map_functions$map_jacobian
#' map_hessian <- map_functions$map_hessian
#' invert_jacobian <- map_functions$invert_jacobian
#' invert_hessian <- map_functions$invert_hessian
#'
#' # Define parameter values in the constrained space
#' x <- c(0, 3, 0.2)
#'
#' # Map to unconstrained space
#' y <- map(x)
#'
#' # Invert mapping back to constrained space
#' x_recovered <- invert(y)
#'
#' # Verify that the mapping is bijective
#' all.equal(x, x_recovered)
#'
#' # Compute derivatives
#' j_map <- map_jacobian(x)
#' h_map <- map_hessian(x)
#' j_invert <- invert_jacobian(y)
#' h_invert <- invert_hessian(y)
make_map_function <- function(lower, upper) {
  # Input validation
  if (length(lower) != length(upper)) {
    stop("lower and upper must have the same length")
  }
  
  npar <- length(lower)
  
  # Pre-allocate result vectors for better performance
  result_template <- numeric(npar)
  
  # Classify transformation types once
  lower_inf <- is.infinite(lower)
  upper_inf <- is.infinite(upper)
  
  # Create transformation type indicators (avoids repeated is.infinite calls)
  # 1 = unbounded, 2 = lower bounded, 3 = upper bounded, 4 = interval bounded
  transform_type <- ifelse(lower_inf & upper_inf, 1L,
                    ifelse(!lower_inf & upper_inf, 2L,
                    ifelse(lower_inf & !upper_inf, 3L, 4L)))
  
  # Pre-compute transformation objects to avoid repeated function calls
  tf_positive <- NULL
  tf_negative <- NULL
  tf_interval <- NULL
  
  # Only create transformation objects if needed
  if (any(transform_type == 2L)) tf_positive <- map_positive()
  if (any(transform_type == 3L)) tf_negative <- map_negative()
  if (any(transform_type == 4L)) tf_interval <- map_interval()
  
  # Create optimized transformation functions using vectorized operations where possible
  list(
    map = function(par) {
      result <- result_template
      
      # Handle each transformation type in batch
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- par[unbounded_idx]
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$map(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$map(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$map(par[interval_bounded_idx], 
                                                       lower[interval_bounded_idx], 
                                                       upper[interval_bounded_idx])
      }
      
      result
    },
    
    invert = function(par) {
      result <- result_template
      
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- par[unbounded_idx]
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$invert(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$invert(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$invert(par[interval_bounded_idx], 
                                                          lower[interval_bounded_idx], 
                                                          upper[interval_bounded_idx])
      }
      
      result
    },
    
    map_jacobian = function(par) {
      result <- result_template
      
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- 1.0
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$map_jacobian(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$map_jacobian(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$map_jacobian(par[interval_bounded_idx], 
                                                               lower[interval_bounded_idx], 
                                                               upper[interval_bounded_idx])
      }
      
      result
    },
    
    map_hessian = function(par) {
      result <- result_template
      
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- 0.0
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$map_hessian(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$map_hessian(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$map_hessian(par[interval_bounded_idx], 
                                                              lower[interval_bounded_idx], 
                                                              upper[interval_bounded_idx])
      }
      
      result
    },
    
    invert_jacobian = function(par) {
      result <- result_template
      
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- 1.0
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$invert_jacobian(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$invert_jacobian(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$invert_jacobian(par[interval_bounded_idx], 
                                                                  lower[interval_bounded_idx], 
                                                                  upper[interval_bounded_idx])
      }
      
      result
    },
    
    invert_hessian = function(par) {
      result <- result_template
      
      unbounded_idx <- transform_type == 1L
      if (any(unbounded_idx)) {
        result[unbounded_idx] <- 0.0
      }
      
      lower_bounded_idx <- transform_type == 2L
      if (any(lower_bounded_idx)) {
        result[lower_bounded_idx] <- tf_positive$invert_hessian(par[lower_bounded_idx], lower[lower_bounded_idx])
      }
      
      upper_bounded_idx <- transform_type == 3L
      if (any(upper_bounded_idx)) {
        result[upper_bounded_idx] <- tf_negative$invert_hessian(par[upper_bounded_idx], upper[upper_bounded_idx])
      }
      
      interval_bounded_idx <- transform_type == 4L
      if (any(interval_bounded_idx)) {
        result[interval_bounded_idx] <- tf_interval$invert_hessian(par[interval_bounded_idx], 
                                                                 lower[interval_bounded_idx], 
                                                                 upper[interval_bounded_idx])
      }
      
      result
    }
  )
}