% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdrm.R
\name{gdrm}
\alias{gdrm}
\title{Generalized distributional regression models
Fits generalized distributional regression models, which are models in which any parameters of a distribution can be modelled as a function (linear, non-linear or both) of covariates.}
\usage{
gdrm(formulae, distrib = normal1(), data, gdrm_control_list = gdrm_control())
}
\arguments{
\item{formulae}{A \code{formulae} object, which consists of formulas for different parameters of the model separated by \code{&} operator.}

\item{distrib}{A \verb{[distrib]} object specifying the distribution assumed for the response variable.}

\item{data}{A \code{data.frame}.}

\item{gdrm_control_list}{A list created by \verb{[gdrm_control()]} function.}
}
\value{
An object of class \code{gdrm}.
}
\description{
Generalized distributional regression models
Fits generalized distributional regression models, which are models in which any parameters of a distribution can be modelled as a function (linear, non-linear or both) of covariates.
}
