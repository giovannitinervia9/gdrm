#' Build smooth components
#'
#' @param smooth_term A smooth term specified by [mgcv::s()].
#' @param data A data.frame containign the variables to build the smooth components.
#'
#' @returns A list containing the `model.matrix` `X` and the penalty matrix `P`.
#'
#' @export
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
build_smooth <- function(smooth_term, data = NULL) {

  if (is.null(data)) {
    data <- get("data", envir = parent.frame())
  }

  if (is.character(smooth_term)) {
    smooth_term <- eval(parse(text = smooth_term))
  }
  r <- mgcv::smoothCon(
    smooth_term,
    data = data,
    absorb.cons = TRUE,
    scale.penalty = FALSE
  )
  nbasis <- length(r)
  Xlist <- vector("list", nbasis)
  Plist <- vector("list", nbasis)
  for (i in 1:nbasis) {
    Xj <- r[[i]]$X
    Pj <- r[[i]]$S[[1]]
    colnames(Xj) <- paste0(r[[i]]$label, "_", 1:ncol(Xj))
    Xlist[[i]] <- Xj
    Plist[[i]] <- Pj
  }
  X <- do.call(cbind, Xlist)
  par <- numeric(ncol(X))
  names(par) <- colnames(X)
  P <- as.matrix(Matrix::bdiag(Plist))
  id <- rep(1:nbasis, sapply(Xlist, ncol))
  colnames(P) <- rownames(P) <- colnames(X)
  fitted <- drop(X%*%par)
  r <- list(X = X,
    P = P,
    par = par,
    fitted = fitted,
    id = id)
  class(r) <- "smooth"
  r
}