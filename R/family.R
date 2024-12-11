#' Functions for specifying a sampling distribution and link function
#'
#' These functions are intended for use in the \code{family} argument of \code{\link{create_sampler}}.
#' In future versions these functions may gain additional arguments, but currently the corresponding
#' functions \code{gaussian} and \code{binomial} can be used as well.
#'
#' @param link the name of a link function. Currently the only allowed link functions are:
#'  \code{"identity"} for (log-)Gaussian sampling distributions, \code{"logit"} (default) and \code{"probit"}
#'  for binomial distributions and \code{"log"} for negative binomial sampling distributions.
#' @param size size or dispersion parameter of the negative binomial distribution used
#'  internally to approximate the Poisson distribution. This should be set to a relatively
#'  large value (default is 100), corresponding to negligible overdispersion, to obtain a
#'  good approximation. However, too large values may cause slow MCMC exploration of the
#'  posterior distribution.
#' @param K number of categories for multinomial model; this must be specified for prior predictive sampling.
#' @param shape.vec optional formula specification of unequal shape parameter for gamma family
#' @param shape.prior prior for gamma shape parameter. Supported prior distributions:
#'  \code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_exp}} and
#'  \code{\link{pr_gamma}}. The current default is \code{pr_gamma(shape=0.1, rate=0.1)}.
#' @param control options for the Metropolis-Hastings algorithm employed
#'  in case the shape parameter is to be inferred. Function \code{\link{set_MH}}
#'  can be used to change the default options. The two choices of proposal
#'  distribution type supported are "RWLN" for a random walk proposal on the
#'  log-shape scale, and "gamma" for an approximating gamma proposal, found using
#'  an iterative algorithm. In the latter case, a Metropolis-Hastings accept-reject
#'  step is currently omitted, so the sampling algorithm is an approximate one,
#'  though often quite accurate and efficient.
#' @param var.data the (variance) data for the gamma part of family \code{gaussian_gamma}.
#' @param ... further arguments passed to \code{f_gamma}.
#' @returns A family object.
#' @name mcmcsae-family
#' @references
#'  J.W. Miller (2019).
#'    Fast and Accurate Approximation of the Full Conditional for Gamma Shape Parameters.
#'    Journal of Computational and Graphical Statistics 28(2), 476-480.
NULL

#' @export
#' @rdname mcmcsae-family
f_gaussian <- function(link="identity") {
  link <- match.arg(link)
  list(family="gaussian", link=link, linkinv=identity)
}

#' @export
#' @rdname mcmcsae-family
f_binomial <- function(link=c("logit", "probit")) {
  link <- match.arg(link)
  if (link == "logit")
    linkinv <- make.link(link)$linkinv
  else
    linkinv <- pnorm
  list(family="binomial", link=link, linkinv=linkinv)
}

#' @export
#' @rdname mcmcsae-family
f_negbinomial <- function(link="logit") {
  link <- match.arg(link)
  list(family="negbinomial", link=link, linkinv=make.link(link)$linkinv)
}

#' @export
#' @rdname mcmcsae-family
f_multinomial <- function(link="logit", K=NULL) {
  link <- match.arg(link)
  if (!is.null(K)) {
    K <- as.integer(K)
    if (length(K) != 1L) stop("number of categories 'K' must be a scalar integer")
    if (K < 2L) stop("number of categories 'K' must be at least 2")
  }
  list(family="multinomial", link=link, linkinv=make.link(link)$linkinv, K=K)
}
