#' Classify terms in formula
#'
#' @param x A character vector indicating a formula term.
#'
#' @returns The type of the formula term indicated.
#' @export
classify_term <- function(x) {
  if (grepl("^\\s*s\\s*\\(", x)) {
    "smooth"
  } else if (grepl("^\\s*nl\\s*\\(", x)) {
    "nl"
  } else if (grepl("^\\s*ar\\s*\\(", x)) {
    "ar"
  } else if (grepl("^\\s*ridge\\s*\\(", x)) {
    "ridge"
  } else if (grepl("^\\s*lasso\\s*\\(", x)) {
    "lasso"
  } else if (grepl("^\\s*re\\s*\\(", x)) {
    "re"
  } else if (grepl("^\\s*enet\\s*\\(", x)) {
    "enet"
  } else if (grepl("^\\s*penp\\s*\\(", x)) {
    "penp"
  } else {
    "linear"
  }
}


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
build_linear <- function(formula, data) {
  formula <- as.formula(formula)
  X <- model.matrix(formula, data)
  P <- matrix(0, ncol(X), ncol(X))
  list(X = X, P = P)
}


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
build_smooth <- function(smooth_term, data) {
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
  P <- as.matrix(Matrix::bdiag(Plist))
  id <- rep(1:nbasis, sapply(Xlist, ncol))
  colnames(P) <- rownames(P) <- colnames(X)
  list(X = X, P = P, id = id)
}


#' Interpret a formula
#'
#' @param formula An object of class `formula`.
#'
#' @returns A list of two elements:
#' \describe{
#' \item{`parts`}{Single terms appearing in the formula.}
#' \item{`types`}{Type of each term.}
#' }
#'
#' @export
interpret_formula <- function(formula) {
  if (gsub(" ", "", trimws(deparse(formula))) == "~1") {
    return(list(
      parts = "~1",
      types = "linear"
    ))
  }

  parts <- attr(terms(formula), "term.labels")

  types <- sapply(parts, classify_term, USE.NAMES = FALSE)

  if (any(types == "linear")) {
    linpart <- paste0("~", paste(parts[types == "linear"], collapse = " + "))
    parts <- parts[types != "linear"]
    types <- types[types != "linear"]
    parts <- c(linpart, parts)
    types <- c("linear", types)
  }

  list(
    parts = parts,
    types = types
  )
}


#' Interpret a formulae
#'
#' @param formulae A formula with different parts related to different parameters separated by `&` operator.
#' @param distrib An object of class [distrib].
#' @param data A data.frame containign the variables in the formula.
#'
#' @returns A list for each parameter of the distribution containing model matrices and penalty matrices
#' for each term of the formula
#'
#' @export
interpret_formulae <- function(formulae, distrib, data) {
  npar <- length(distrib$parameters)
  fchr <- deparse(formulae)
  f_split <- as.list(trimws(strsplit(fchr, "&")[[1]]))
  nform <- length(f_split)

  if (nform < npar) {
    f_split[(nform + 1):npar] <- "~1"
  }

  names(f_split) <- distrib$parameters

  f_split <- lapply(f_split, function(f) {
    if (!grepl("~", f)) {
      paste0("~ ", f)
    } else {
      f
    }
  })

  f_list <- lapply(f_split, as.formula)

  response <- deparse(f_list[[1]][[2]])

  f_int <- lapply(f_list, interpret_formula)

  f_comp <- vector("list", nform)
  names(f_comp) <- distrib$parameters

  for (i in 1:nform) {
    parts <- f_int[[i]]$parts
    types <- f_int[[i]]$types
    funs <- paste0("build_", types)
    f_comp[[i]] <- Map(
      function(funs, parts) match.fun(funs)(parts, data),
      funs = funs,
      parts = parts
    )
  }

  f_comp
}
