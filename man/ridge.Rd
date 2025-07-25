% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridge.R
\name{ridge}
\alias{ridge}
\title{Specification of a ridge component}
\usage{
ridge(formula, scale_numerical = TRUE, data = NULL)
}
\arguments{
\item{formula}{A formula specifying a linear predictor to be penalized with ridge penalty.}

\item{scale_numerical}{Logical. If \code{TRUE} (default), numerical predictors are standardized.}

\item{data}{A data.frame containign the variables in the formula.}
}
\value{
A list containing the \code{model.matrix} \code{X}, the diagonal penalty matrix \code{P} and center and scale values of numerical variables saved in the list \code{scale_info}.
}
\description{
Specification of a ridge component
}
