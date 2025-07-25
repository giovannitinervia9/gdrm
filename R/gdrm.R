
#' Generalized distributional regression models
#' Fits generalized distributional regression models, which are models in which any parameters of a distribution can be modelled as a function (linear, non-linear or both) of covariates.
#'
#' @param formulae A `formulae` object, which consists of formulas for different parameters of the model separated by `&` operator. 
#' @param distrib A `[distrib]` object specifying the distribution assumed for the response variable.
#' @param data A `data.frame`.
#'
#' @returns An object of class `gdrm`.
#'
#' @export
gdrm <- function(formulae, distrib = normal1(), data) {

  

}