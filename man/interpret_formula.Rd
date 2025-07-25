% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulae.R
\name{interpret_formula}
\alias{interpret_formula}
\title{Interpret a formula}
\usage{
interpret_formula(formula, force_intercept = TRUE)
}
\arguments{
\item{formula}{An object of class \code{formula}.}

\item{force_intercept}{Logical. If \code{TRUE} (default) a global intercept is forced to the model if none linear part is present.}
}
\value{
A list of two elements:
\describe{
\item{\code{parts}}{Single terms appearing in the formula.}
\item{\code{types}}{Type of each term.}
}
}
\description{
Interpret a formula
}
