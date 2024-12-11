#' Draw from a multivariate normal distribution, given the Cholesky decomposition of its precision matrix
#'
#' @noRd
#' @param ch a cholesky closure for a precision matrix Q.
#' @param xy numeric vector, defaults to 0. The mean of the multivariate
#'           normal distribution is computed as Q^{-1} xy.
#' @param sd scalar standard deviation factor passed to \code{rnorm}.
#' @returns A draw from the multivariate normal distribution N(Q^{-1} xy, sd^2 Q^{-1}).
drawMVN_cholQ <- function(ch, xy=NULL, sd=1) {
  z <- Crnorm(ch$size, sd=sd)  # NB sd must be scalar in Crnorm
  if (is.null(xy))
    ch$solve(z, system="Lt")
  else  # draw from multivariate normal N(Q^-1 xy, sd^2 Q^-1)
    ch$solve(ch$solve(xy, system="L") + z, system="Lt")
}

drawMVN_Q <- function(M, sd=1, perm=FALSE)
  drawMVN_cholQ(build_chol(M, control=chol_control(perm=perm)), sd=sd)

# q: dimension of M
# M: (symmetric) matrix / dsCMatrix
# sd: standard deviation
drawMVN_cholV <- function(q, M, sd=1) {
  z <- Crnorm(q, sd=sd)  # NB sd must be scalar in Crnorm
  crossprod_mv(M, z)
}

#' Draw from a scaled chi-square distribution with degrees of freedom and inverse scale parameters or their product
#' 
#' @noRd
#' @param n number of draws.
#' @param nu degrees of freedom (vector) parameter.
#' @param s2 inverse scale (vector) parameter.
#' @param psi as an alternative to \code{s2}, the product of \code{nu} and \code{s2} can be passed as \code{psi}.
#' @param A (vector of) scaled chi-square draws.
rchisq_scaled <- function(n, nu, s2, psi) {
  if (missing(s2))
    rchisq(n, nu) * (1 / psi)
  else
    rchisq(n, nu) * (1 / (nu * s2))
}

# draw from (scaled) beta prime distribution
# n number of draws
# a, b shape parameters
# q scale parameter
draw_betaprime <- function(n, a, b, q=1) {
  z <- rbeta(n, a, b)
  q * z / (1-z)
}

#' Compute binomial coefficients
#' 
#' @noRd
#' @param n vector of numbers of trials.
#' @param y vector of numbers of successes.
#' @param log whether to compute the log binomial coefficient (default).
#' @returns Vector of binomial coefficients 'n over y'.
binomial_coef <- function(n, y, log=TRUE) {
  bc <- lgamma(n + 1) - lgamma(y + 1) - lgamma(n - y + 1)
  if (log) bc else exp(bc)
}

#' Compute binomial coefficients for negative binomial distribution
#' 
#' @noRd
#' @param r vector of dispersion parameters.
#' @param y vector of numbers of successes.
#' @param log whether to compute the log binomial coefficient (default).
#' @returns Vector of negative binomial coefficients 'y+r-1 over y'.
negbinomial_coef <- function(r, y, log=TRUE) {
  bc <- lgamma(y + r) - lgamma(y + 1) - lgamma(r)
  if (log) bc else exp(bc)
}

#' Draw a sample from a Multivariate-Log-inverse-Gamma distribution
#'
#' @noRd
#' @param m dimension.
#' @param alpa shape (vector) parameter.
#' @param kappa another (vector) parameter.
#' @param log.kappa logarithm of \code{kappa}. This can be specified instead of
#'  \code{kappa} if numerical computation of the latter would overflow.
#' @returns A draw from the MLiG distribution.
rMLiG <- function(m, alpha, kappa, log.kappa) {
  if (missing(kappa)) {
    # version where log.kappa is specified; useful in case kappa overflows
    if (any(alpha < 0.1)) {
      # prevent underflow, https://stats.stackexchange.com/questions/7969/how-to-quickly-sample-x-if-expx-gamma
      # TODO only for those components with alpha < 0.1
      # NB kappa can still underflow to zero causing unbounded results
      log.kappa - log(rgamma(m, shape = alpha + 1, rate = 1)) - log(runif(m))/alpha
    } else
      log.kappa - log(rgamma(m, shape = alpha, rate = 1))
  } else {
    if (any(alpha < 0.1)) {
      # prevent underflow, https://stats.stackexchange.com/questions/7969/how-to-quickly-sample-x-if-expx-gamma
      # TODO only for those components with alpha < 0.1
      # NB kappa can still underflow to zero causing unbounded results
      -(log(rgamma(m, shape = alpha + 1, rate = kappa)) + log(runif(m))/alpha)
    } else
      -log(rgamma(m, shape = alpha, rate = kappa))
  }
}
