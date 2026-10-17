
#' Set options for random number generation
#'
#' @export
#' @param seed random seed.
#' @param use.randompack if \code{TRUE} use R package randompack
#'  for faster (vectorised) random number generation for a subset
#'  of the distributions used.
#' @param randompack.engine the rng engine used if \code{use.randompack=TRUE}.
#'  See \code{\link[randompack]{randompack_engines}}.
#' @return an environment with rng options and functions.
set_rng <- function(seed=NULL, use.randompack=FALSE, randompack.engine="x256++simd") {

  engine <- NULL

  if (!is.null(seed)) set.seed(seed)

  if (use.randompack) {
    draw_unif   <- function(n, min=0, max=1) engine$unif(n, min, max)
    draw_normal <- function(n, mean=0, sd=1) engine$normal(n, mean, sd)
    draw_t      <- function(n, df, ncp) engine$t(n, df)
    draw_chisq  <- function(n, df, ncp = 0) engine$chi2(n, df)
    draw_gamma  <- function(n, shape, rate=NULL, scale=1) {
      mult <- if (is.null(rate)) scale else 1/rate
      mult * engine$gamma(n, shape)
    }
    draw_exp    <- function(n, rate = 1) if (rate == 1) engine$exp(n) else (1/rate) * engine$exp(n)
    draw_MLiG <- function(n, alpha, kappa, log.kappa) {
      if (missing(kappa)) {
        if (any(alpha < 0.1))
          log.kappa - logC(draw_gamma(n, shape = alpha + 1, rate = 1)) - logC(draw_unif(n))/alpha
        else
          log.kappa - logC(draw_gamma(n, shape = alpha, rate = 1))
      } else {
        if (any(alpha < 0.1))
          -(logC(draw_gamma(n, shape = alpha + 1, rate = kappa)) + logC(draw_unif(n))/alpha)
        else
          -logC(draw_gamma(n, shape = alpha, rate = kappa))
      }
    }
  } else {
    draw_unif   <- runif
    draw_normal <- Crnorm
    draw_t      <- rt
    draw_chisq  <- rchisq
    draw_gamma  <- rgamma
    draw_exp    <- rexp
    draw_MLiG   <- rMLiG
  }

  env <- environment()

  if (use.randompack) {
    init_engine(env, worker_id = 1)
  }

  env
}

#' Initialize or revive a randompack engine safely
#'
#' @noRd
#' @param rng_env The environment returned by set_rng()
#' @param worker_id Integer stream ID for parallel workers
init_engine <- function(rng_env, worker_id = 1) {

  if (!rng_env[["use.randompack"]]) return(invisible(NULL))

  if (!requireNamespace("randompack", quietly=TRUE)) stop("package randompack required")

  eng_type <- rng_env[["randompack.engine"]]
  m_seed <- rng_env[["seed"]]

  new_engine <- randompack::randompack_rng(engine = eng_type)

  # set seed and make sure different workers have independent streams in case
  if (!is.null(m_seed)) {
    if (grepl("x256|xoro|x128", eng_type, ignore.case = TRUE)) {
      new_engine$seed(m_seed)
      if (worker_id > 1) {
        for (j in seq_len(worker_id - 1)) new_engine$jump()
      }
    } else if (grepl("pcg|philox|chacha|squares", eng_type, ignore.case = TRUE)) {
      new_engine$seed(m_seed, stream = worker_id)
    } else {
      derived_seed <- (m_seed + worker_id * 1000007L) %% .Machine$integer.max
      new_engine$seed(derived_seed)
    }
  }

  # inject the created engine into the closure's environment
  rng_env$engine <- new_engine

  invisible(NULL)
}


#' Draw from a multivariate normal distribution, given the Cholesky decomposition of its precision matrix
#'
#' @noRd
#' @param ch a cholesky closure for a precision matrix Q.
#' @param xy numeric vector, defaults to 0. The mean of the multivariate
#'           normal distribution is computed as Q^{-1} xy.
#' @param sd scalar standard deviation factor passed to \code{rnorm}.
#' @returns A draw from the multivariate normal distribution N(Q^{-1} xy, sd^2 Q^{-1}).
drawMVN_cholQ <- function(ch, xy=NULL, sd=1) {
  z <- Crnorm(ch[["size"]], sd=sd)  # NB sd must be scalar in Crnorm
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
drawMVN_cholV <- function(q, M, sd=1) tm_op(M, Crnorm(q, sd=sd))

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
  if (log) bc else expC(bc)
}

#' Compute binomial coefficients for negative binomial distribution
#' 
#' @noRd
#' @param r vector of dispersion shape parameters.
#' @param y vector of numbers of successes.
#' @param log whether to compute the log binomial coefficient (default).
#' @returns Vector of negative binomial coefficients 'y+r-1 over y'.
negbinomial_coef <- function(r, y, log=TRUE) {
  bc <- lgamma(y + r) - lgamma(y + 1) - lgamma(r)
  if (log) bc else expC(bc)
}

#' Draw a sample from a Multivariate-Log-inverse-Gamma distribution
#'
#' @noRd
#' @param n dimension.
#' @param alpa shape (vector) parameter.
#' @param kappa another (vector) parameter.
#' @param log.kappa logarithm of \code{kappa}. This can be specified instead of
#'  \code{kappa} if numerical computation of the latter would overflow.
#' @returns A draw from the MLiG distribution.
rMLiG <- function(n, alpha, kappa, log.kappa) {
  # prevent underflow for small alpha: https://stats.stackexchange.com/questions/7969/how-to-quickly-sample-x-if-expx-gamma
  # TODO only for those components with alpha < 0.1
  # NB kappa can still underflow to zero causing unbounded results
  if (missing(kappa)) {
    # version where log.kappa is specified; useful in case kappa overflows
    if (any(alpha < 0.1))
      log.kappa - logC(rgamma(n, shape = alpha + 1, rate = 1)) - logC(runif(n))/alpha
    else
      log.kappa - logC(rgamma(n, shape = alpha, rate = 1))
  } else {
    if (any(alpha < 0.1))
      -(logC(rgamma(n, shape = alpha + 1, rate = kappa)) + logC(runif(n))/alpha)
    else
      -logC(rgamma(n, shape = alpha, rate = kappa))
  }
}

get_PG_sampler <- function(n, PG.approx, PG.approx.m) {
  if (PG.approx) {
    mPG <- as.integer(PG.approx.m)
    if (all(length(mPG) != c(1L, n))) stop("invalid value for option 'PG.approx.m'")
    function(b, c) CrPGapprox(n, b, c, mPG)
  } else {
    if (!requireNamespace("BayesLogit", quietly=TRUE)) stop("please install package 'BayesLogit' and try again")
    rpg <- BayesLogit::rpg
    function(b, c) rpg(n, b, c)
  }
}
