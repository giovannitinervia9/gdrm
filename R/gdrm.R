#' Control options for gdrm()
#' @param force_intercept Logical. If `TRUE` (default) a global intercept is forced to the model if none linear part is present.
#' @param expected Logical. If `TRUE` (default) use expected Fisher Information if available.
#' @returns A list of class `gdrm_control`.
#'
#' @export
gdrm_control <- function(force_intercept = TRUE, expected = TRUE) {
  r <- list(
    force_intercept = force_intercept,
    expected = expected)
  class(r) <- "gdrm_control"
  r
}


#' Generalized distributional regression models
#' Fits generalized distributional regression models, which are models in which any parameters of a distribution can be modelled as a function (linear, non-linear or both) of covariates.
#'
#' @param formulae A `formulae` object, which consists of formulas for different parameters of the model separated by `&` operator.
#' @param distrib A `[distrib]` object specifying the distribution assumed for the response variable.
#' @param data A `data.frame`.
#' @param gdrm_control_list A list created by `[gdrm_control()]` function.
#'
#' @returns An object of class `gdrm`.
#'
#' @export
gdrm <- function(
  formulae,
  distrib = normal1(),
  data,
  gdrm_control_list = gdrm_control()
) {
  force_intercept <- gdrm_control_list$force_intercept
  expected <- gdrm_control_list$expected

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
  npars <- lapply(mod_comp, function(f) {
    sum(sapply(f, function(ff) length(ff$par)))
  })

  # function to build hessian
}
