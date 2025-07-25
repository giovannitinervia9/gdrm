#' Specification of non linear component
#'
#' @param formula A right-sided formula specifying a non linear expression.
#' @param parameters A list with elements of the same name of the parameters specifying a vector of lower and upper bounds.
#' @param data A data.frame containing variables used to compute the non linear function.
#'
#' @returns A list with tree elements:
#' \describe{
#' \item{`f`}{A function with argument the parameters `par` to compute the non linear expression.}
#' \item{`jac`}{A function with argument the parameters `par` to compute the jacobian of `f`.}
#' \item{`hes`}{A function with argument the parameters `par` to compute the hessian of `f`.}
#' }
#'
#' @export
#' 
#' @importFrom numDeriv jacobian hessian
#' @importFrom Deriv Deriv
#' @examples
#' formula <- ~ theta0 * (x1^theta1) * (x2^theta2)
#' parameters <- list(theta0 = c(0, Inf), theta1 = c(0, Inf), theta2 = c(0, Inf))
#' data <- abs(data.frame(x1 = rnorm(5), x2 = rnorm(5)))
#' r <- nl(formula, parameters, data)
#' par <- list(1, 2, 3)
#' r$f(par)
#' r$jac(par)
#' r$hes(par)
nl <- function(formula, parameters, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  ff <- deparse1(formula)
  ff <- gsub("~", "", ff)

  vars <- colnames(data)
  pars <- names(parameters)
  npars <- length(pars)

  if (any(vars %in% pars) | any(pars %in% vars)) {
    stop(
      "Some variables and parameters share same name! Make sure variables and parameters are unique"
    )
  }

  assign_str <- character(npars)
  for (i in 1:npars) {
    assign_str[i] <- paste0(pars[i], " <- par[[", i, "]]")
  }

  expr <- paste0(
    paste0(paste0(assign_str, collapse = ";"), ";"),
    ff,
    collapse = ";"
  )

  create_function <- function(expr_string, data_env) {
    expr_parsed <- parse(text = paste("with(data_env, {", expr_string, "})"))
    f <- function(par) {}
    body(f) <- expr_parsed[[1]]
    environment(f)$data_env <- data_env
    f
  }

  f <- create_function(expr, data)

  jac <- tryCatch(
    expr = {
      expr_j <- Deriv::Deriv(f = ff, x = pars, nderiv = 1, combine = "cbind")

      if (!grepl("cbind", expr_j)) {
        expr_j <- paste0("cbind(", expr_j, ")")
      }

      expr_j <- sub("cbind", "jj <- cbind", expr_j)

      expr_j <- paste0(
        paste0(paste0(assign_str, collapse = ";"), ";"),
        expr_j,
        collapse = ";"
      )

      expr_j <- paste0(expr_j, "; colnames(jj) <- pars; jj")

      create_function(expr_j, data)
    },
    error = function(e) {
      jac <- function(par) {
        par <- unlist(par)
        numDeriv::jacobian(func = f, x = par)
      }
      environment(jac) <- environment()
      jac
    }
  )

  hes <- tryCatch(
    expr = {
      expr_h <- Deriv::Deriv(f = ff, x = pars, nderiv = 2, combine = "cbind")

      if (!grepl("cbind", expr_h)) {
        expr_h <- paste0("cbind(", expr_h, ")")
      }
      
      expr_h <- sub("cbind", "hh <- cbind", expr_h)

      expr_h <- paste0(
        paste0(paste0(assign_str, collapse = ";"), ";"),
        expr_h,
        collapse = ";"
      )
      expr_h <- paste0(
        expr_h,
        "; hh <- simplify2array(apply(hh, 1, FUN = function(x) matrix(x, npars, npars, byrow = TRUE), simplify = FALSE));",
        "if (is.null(dim(hh))) hh <- array(hh, dim = c(1, 1, length(hh)));",
        "dimnames(hh) <- list(pars, pars, NULL); hh"
      )
      create_function(expr_h, data)
    },
    error = function(e) {
      hes <- function(par) {
        par <- unlist(par)
        hlist <- vector("list", nrow(data))
        for (i in 1:nrow(data)) {
          xj <- data[i, , drop = FALSE]
          f_temp <- create_function(expr, xj)
          hlist[[i]] <- numDeriv::hessian(func = f_temp, x = par)
        }
        simplify2array(hlist)
      }
      environment(hes) <- environment()
      hes
    }
  )

  P <- matrix(0, npars, npars)
  
  # initialize parameters value
  par <- sapply(parameters, function(x) x[1]) + 1e-10
  names(par) <- pars
  fitted = function(par) {drop(f(par))}
  
  r <- list(f = f,
    jac = jac,
    hes = hes, 
    parameters = parameters,
    par = par,
    P = P,
    fitted = fitted)
  class(r) <- "nl"
  r
}


#' Build non linear components
#'
#' @param object A character vector representing a call to [`nl()`] function.
#' @param data A data.frame containign the variables in the formula for [`nl()`].
#' @returns An object created by [`nl()`] function.
#'
#' @export
build_nl <- function(object, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  if(!grepl("data", object)) {
    object <- sub("\\)$", ", data = data)", object)
  }

  eval(parse(text = object))
}