#' Build linear components
#'
#' @param formula A formula specifying a linear predictor.
#' @param data A data.frame containign the variables in the formula.
#'
#' @returns A list containing the `model.matrix` `X` and the zero Penalty
#' matrix `P`.
#'
#' @export
#' @importFrom stats model.matrix as.formula terms
build_linear <- function(formula, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  formula <- as.formula(formula)
  X <- model.matrix(formula, data)
  par <- numeric(ncol(X))
  names(par) <- colnames(X)
  P <- matrix(0, ncol(X), ncol(X))
  fitted <- drop(X%*%par)
  r <- list(X = X,
    P = P,
    par = par,
    fitted = fitted)
  class(r) <- "linear"
  r
}