#' Specification of a ridge component
#'
#' @param formula A formula specifying a linear predictor to be penalized with ridge penalty.
#' @param scale_numerical Logical. If `TRUE` (default), numerical predictors are standardized.
#' @param data A data.frame containign the variables in the formula.
#'
#' @returns A list containing the `model.matrix` `X`, the diagonal penalty matrix `P` and center and scale values of numerical variables saved in the list `scale_info`.
#'
#' @export
#' @importFrom stats model.matrix as.formula terms
ridge <- function(formula, scale_numerical = TRUE, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  formula <- as.formula(formula)
  
  # get variables names
  vars <- attr(terms(formula), "variables")
  vars <- as.list(vars)[-1]
  vars <- vapply(vars, deparse, character(1))
  num_vars <- vars[sapply(data[vars], class) %in% c("numeric", "integer")]
  
  if (scale_numerical) {
    num_scaled <- scale(data[num_vars])
    center <- attr(num_scaled, "scaled:center")
    scale <- attr(num_scaled, "scaled:scale")
    data[num_vars] <- num_scaled
  } else {
    center <- rep(0, length(num_vars))
    scale <- rep(1, length(num_vars))
    names(center) <- names(scale) <- num_vars
  }

  X <- model.matrix(formula, data)[,-1]

  par <- numeric(ncol(X))
  names(par) <- colnames(X)
  hyperpar <- .001

  P <- diag(1, ncol(X), ncol(X))
  fitted <- function(par) {drop(X%*%par)}

  r <- list(X = X,
    P = P,
    par = par,
    parameters = colnames(X),
    hyperpar = hyperpar,
    fitted = fitted,
    scale_info = list(center = center, scale = scale))
  class(r) <- "ridge"
  r
}


#' Build ridge components
#'
#' @param object A character vector representing a call to [`ridge()`] function.
#' @param data A data.frame containign the variables in the formula for [`ridge()`].
#' @returns An object created by [`ridge()`] function.
#'
#' @export
build_ridge <- function(object, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  if(!grepl("data", object)) {
    object <- sub("\\)$", ", data = data)", object)
  }

  eval(parse(text = object))
}