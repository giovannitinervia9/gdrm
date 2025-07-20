#-------------------------------------------------------------------------------

#' Create a link function object for gdrm
#'
#' This function creates a link function object similar to `stats::make.link()` but extends it by
#' including the second derivative of theta with respect to eta (`theta2.eta2`).
#'
#' @param link Character string specifying the link function. Options include:
#'   \itemize{
#'     \item \code{"logit"}: log(theta/(1-theta))
#'     \item \code{"probit"}: inverse of the standard normal CDF
#'     \item \code{"cauchit"}: inverse of the Cauchy CDF
#'     \item \code{"cloglog"}: log(-log(1-theta))
#'     \item \code{"identity"}: no transformation
#'     \item \code{"log"}: natural logarithm
#'     \item \code{"sqrt"}: square root
#'     \item \code{"1/theta^2"}: inverse square
#'     \item \code{"inverse"}: reciprocal
#'   }
#'
#' @return A list with S3 class \code{"link-gnlm"} containing the following components:
#'   \itemize{
#'     \item \code{linkfun}: Function that transforms the theta to eta
#'     \item \code{linkinv}: Function that transforms eta to theta
#'     \item \code{theta.eta}: Function returning the derivative of theta with respect to eta (\eqn{d\theta/d\eta})
#'     \item \code{theta2.eta2}: Function returning the second derivative of theta with respect to eta (\eqn{d^2\theta/d^2\eta^2})
#'     \item \code{valideta}: Function that checks if the eta values are valid
#'     \item \code{name}: Character string naming the link function
#'   }
#'
#' @note This function extends \code{stats::make.link()} by including the \code{theta2.eta2} component,
#' which calculates the second derivative of theta with respect to eta.
#'
#'
#' @export
#' @importFrom stats dcauchy dnorm pcauchy pnorm qcauchy qnorm
make_link <- function(link) {
  switch(
    link,

    logit = {
      linkfun <- function(theta) .Call(C_logit_link, theta)
      linkinv <- function(eta) .Call(C_logit_linkinv, eta)
      theta.eta <- function(eta) .Call(C_logit_theta_eta, eta)
      theta2.eta2 <- function(eta) {
        2 * exp(2 * eta) / (1 + exp(eta))^3 - exp(eta) / (1 + exp(eta))^2
      }
      valideta <- function(eta) TRUE
    },

    probit = {
      linkfun <- function(theta) qnorm(theta)
      linkinv <- function(eta) {
        thresh <- -qnorm(.Machine$double.eps)
        eta <- pmin(pmax(eta, -thresh), thresh)
        pnorm(eta)
      }
      theta.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
      theta2.eta2 <- function(eta) pmax(eta * dnorm(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    },

    cauchit = {
      linkfun <- function(theta) qcauchy(theta)
      linkinv <- function(eta) {
        thresh <- -qcauchy(.Machine$double.eps)
        eta <- pmin(pmax(eta, -thresh), thresh)
        pcauchy(eta)
      }
      theta.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
      theta2.eta2 <- function(eta) {
        pmax(-2 * eta / (pi * (1 + eta^2)^2), .Machine$double.eps)
      }
      valideta <- function(eta) TRUE
    },

    cloglog = {
      linkfun <- function(theta) log(-log(1 - theta))
      linkinv <- function(eta) {
        pmax(
          pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps),
          .Machine$double.eps
        )
      }
      theta.eta <- function(eta) {
        eta <- pmin(eta, 700)
        pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
      }
      theta2.eta2 <- function(eta) {
        eta <- pmin(eta, 700)
        pmax(-exp(eta - exp(eta)) * (exp(eta) - 1), .Machine$double.eps)
      }
      valideta <- function(eta) TRUE
    },

    identity = {
      linkfun <- function(theta) theta
      linkinv <- function(eta) eta
      theta.eta <- function(eta) rep.int(1, length(eta))
      theta2.eta2 <- function(eta) rep.int(0, length(eta))
      valideta <- function(eta) TRUE
    },

    log = {
      linkfun <- function(theta) log(theta)
      linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
      theta.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
      theta2.eta2 <- function(eta) pmax(exp(eta), .Machine$double.eps)
      valideta <- function(eta) TRUE
    },

    sqrt = {
      linkfun <- function(theta) sqrt(theta)
      linkinv <- function(eta) eta^2
      theta.eta <- function(eta) 2 * eta
      theta2.eta2 <- function(eta) 2
      valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
    },

    `1/theta^2` = {
      linkfun <- function(theta) 1 / theta^2
      linkinv <- function(eta) 1 / sqrt(eta)
      theta.eta <- function(eta) -1 / (2 * eta^1.5)
      theta2.eta2 <- function(eta) 3 / (4 * eta^2.5)
      valideta <- function(eta) all(is.finite(eta)) && all(eta > 0)
    },

    inverse = {
      linkfun <- function(theta) 1 / theta
      linkinv <- function(eta) 1 / eta
      theta.eta <- function(eta) -1 / (eta^2)
      theta2.eta2 <- function(eta) 2 / eta^3
      valideta <- function(eta) all(is.finite(eta)) && all(eta != 0)
    },
    stop(gettextf("%s link not recognised", sQuote(link)), domain = NA)
  )
  environment(linkfun) <- environment(linkinv) <- environment(
    theta.eta
  ) <- environment(valideta) <- asNamespace("stats")
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      theta.eta = theta.eta,
      theta2.eta2 = theta2.eta2,
      valideta = valideta,
      name = link
    ),
    class = "link-gdrm"
  )
}
