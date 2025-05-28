
#' Create an object representing a possibly multivariate normal prior distribution
#'
#' @export
#' @param mean scalar or vector mean parameter.
#' @param precision scalar, vector or matrix precision parameter.
#' @param labels optional character vector with coefficient labels. If specified,
#'  it should have the same length as at least one of \code{mean} and \code{precision},
#'  and in that case the normal prior with these parameters is assigned to these coefficients,
#'  while any coefficients not present in labels will be assigned a non-informative
#'  prior with mean 0 and precision 0.
#' @returns An environment representing the specified prior, for internal use.
# TODO allow to not implicitly include a factor of 1/sigma_^2 in precision for gaussian model
#      including this factor only has a computational advantage in case of gaussian variance sigma^2 Sigma0, where Sigma0 is fixed
pr_normal <- function(mean=0, precision=0, labels=NULL) {
  mean <- as.numeric(mean)
  if (!is.null(labels) && all(length(mean) != c(1L, length(labels))))
    stop("if 'labels' are specified, parameter 'mean' must have the same length, or be a scalar")
  if (is.vector(precision)) {
    if (!is.null(labels) && all(length(precision) != c(1L, length(labels)))) stop("if 'labels' are specified, the size of 'precision' must be compatible if it is not a scalar")
    if (any(precision < 0)) stop("parameter 'precision' must be nonnegative")
    informative <- any(precision != 0)
  } else if (is_a_matrix(precision)) {
    if (!is.null(labels) && all(nrow(precision) != c(1L, length(labels)))) stop("if 'labels' are specified, the size of 'precision' must be compatible if it is not a scalar")
    precision <- economizeMatrix(precision, symmetric=TRUE, check=TRUE)
    informative <- !is_zero_matrix(precision)
    # TODO check for positive semi-definiteness
  } else stop("unexpected precision parameter")
  n <- NULL
  rprior <- function(p) stop("please call method 'init' first")
  sigma <- FALSE  # does the prior variance contain (implicitly) a scalar factor sigma_^2
  init <- function(n=1L, coefnames=NULL, sparse=NULL, sigma=FALSE) {
    n <<- as.integer(n)
    if (is.null(labels)) {
      if (all(length(mean) != c(1L, n))) stop("parameter 'mean' has wrong length")
      if (all(NROW(precision) != c(1L, n))) stop("parameter 'precision' has wrong size")
      if (is.vector(precision)) {
        if (length(precision) == 1L)
          precision <<- Cdiag(rep.int(precision, n))
        else
          precision <<- Cdiag(precision)
      }
    } else {
      if (length(coefnames) != n) stop("'coefnames' must have length 'n'")
      m <- fmatch(labels, coefnames)
      if (anyNA(m)) stop("non-matching labels: ", paste0(labels[is.na(m)], collapse=", "))
      temp <- numeric(n)
      temp[m] <- mean
      mean <<- temp
      if (is.vector(precision)) {
        temp <- numeric(n)
        temp[m] <- precision
        precision <- Cdiag(temp)
      } else {
        temp <- new("dsCMatrix", p=rep(0L, n + 1L), uplo="U", Dim=c(n, n))
        temp[m, m] <- precision
        precision <- temp
      }
    }
    # sparse: block sampler expects sparse matrix with x-slot
    precision <<- economizeMatrix(precision, sparse=sparse, symmetric=TRUE)
    if (isTRUE(sparse) && is_unit_ddi(precision)) precision <<- expand_unit_ddi(precision)
    sigma <<- sigma
    if (sigma) {
      # combined conjugate prior for ordinary gaussian regression --> sigma^2 factor in coefficients' prior variance
      rprior <<- function(p) mean + drawMVN_Q(precision, sd=p[["sigma_"]])
    } else
      rprior <<- function(p) mean + drawMVN_Q(precision)
  }
  log_prior_0 <- NULL
  log_prior <- function(p) NULL
  setup_logprior <- function(name) {
    log_prior_0 <<- -0.5 * n * log(2*pi)
    if (informative) {
      d <- diag(suppressWarnings(chol(precision, pivot=TRUE)))
      log_prior_0 <<- log_prior_0 + sum(log(d[d > 0]))
    }
    log_prior <<- function(p) {}
    log_prior <<- add(bquote(delta.beta <- p[[.(name)]] - mean))
    if (sigma) {
      log_prior <<- log_prior |>
        add(quote(sigma <- p[["sigma_"]])) |>
        add(quote(- n * log(sigma) - 0.5 * dotprodC(delta.beta, precision %m*v% delta.beta) / sigma^2))
    } else {
      log_prior <<- add(log_prior, quote(- 0.5 * dotprodC(delta.beta, precision %m*v% delta.beta)))
    }
  }
  type <- "normal"
  environment()
}

#' Create an object representing a Multivariate Log inverse Gamma (MLiG) prior distribution
#'
#' @export
#' @param mean scalar or vector parameter for the mean in the large
#'  \code{a} limit, when the distribution approaches a normal distribution.
#' @param precision scalar or vector parameter for the precision in the
#'  large \code{a} limit, when the distribution approaches a normal
#'  distribution.
#' @param labels optional character vector with coefficient labels. If specified,
#'  it should have the same length as at least one of \code{mean} and \code{precision},
#'  and in that case the MLiG prior with these parameters is assigned to these coefficients,
#'  while any coefficients not present in labels will be assigned a non-informative
#'  prior with mean 0 and precision 0.
#' @param a scalar parameter that controls how close the prior is to independent
#'  normal priors with \code{mean} and \code{precision} parameters. The larger
#'  this value (default is 1000), the closer.
#' @returns An environment representing the specified prior, for internal use.
#' @references
#'  J.R. Bradley, S.H. Holan and C.K. Wikle (2018).
#'    Computationally efficient multivariate spatio-temporal models for
#'    high-dimensional count-valued data (with discussion).
#'    Bayesian Analysis 13(1), 253-310.
pr_MLiG <- function(mean=0, precision=0, labels=NULL, a=1000) {
  if (length(a) != 1L || a <= 0) stop("parameter 'a' must be a positive scalar")
  if (any(precision < 0)) stop("parameter 'precision' must be nonnegative")
  informative <- any(precision > 0)
  n <- NULL
  rprior <- function(p) stop("please call method 'init' first")
  init <- function(n=1L, coefnames) {
    n <<- as.integer(n)
    if (is.null(labels)) {
      if (all(length(mean) != c(1L, n))) stop("parameter 'mean' has wrong length")
      if (all(length(precision) != c(1L, n))) stop("parameter 'precision' has wrong length")
    } else {
      if (all(length(mean) != c(1L, length(labels)))) stop("parameter 'mean' has wrong length")
      if (all(length(precision) != c(1L, length(labels)))) stop("parameter 'precision' has wrong length")
      m <- fmatch(labels, coefnames)
      if (anyNA(m)) stop("non-matching coefficient names: ", paste0(labels[is.na(m)], collapse=", "))
      temp <- numeric(n)
      temp[m] <- mean
      mean <<- temp
      temp <- numeric(n)
      temp[m] <- precision
      precision <<- temp
    }
    rprior <<- function(p) mean + sqrt(a/precision) * rMLiG(n, a, a)
  }
  type <- "MLiG"
  environment()
}

#' Create an object representing a degenerate prior fixing a
#' parameter (vector) to a fixed value
#'
#' @export
#' @param value scalar or vector value parameter.
#' @returns An environment representing the specified prior, for internal use.
pr_fixed <- function(value=1) {
  value <- as.numeric(value)
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(value) != c(1L, n))) stop("value parameter has wrong length")
    if (n > length(value)) value <- rep.int(value, n)
    rprior <<- function() value
  }
  type <- "fixed"
  environment()
}

#' Create an object representing uniform prior distributions
#'
#' @export
#' @param min lower limit.
#' @param max upper limit.
#' @returns An environment representing the specified prior, for internal use.
pr_unif <- function(min=0, max=1) {
  if (!all(min <= max)) stop("'min' should not exceed 'max'")
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  # ldr: log-density-ratio, i.e. log(dbeta(x1)/dbeta(x2)), used in MH log-acceptance-rate
  ldr <- NULL
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(min) != c(1L, n))) stop("parameter 'min' has wrong length")
    if (all(length(max) != c(1L, n))) stop("parameter 'max' has wrong length")
    rprior <<- function() runif(n, min, max)
    ldr <<- function(x1, x2) 0
  }
  type <- "unif"
  environment()
}

#' Create an object representing beta prior distributions
#'
#' @export
#' @param a positive shape parameter.
#' @param b positive shape parameter.
#' @returns An environment representing the specified prior, for internal use.
pr_beta <- function(a=1, b=1) {
  if (!all(a > 0 & b > 0)) stop("'a' and 'b' must be positive")
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  # ldr: log-density-ratio, i.e. log(dbeta(x1)/dbeta(x2)), used in MH log-acceptance-rate
  ldr <- NULL
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(a) != c(1L, n))) stop("parameter 'a' has wrong length")
    if (all(length(b) != c(1L, n))) stop("parameter 'b' has wrong length")
    rprior <<- function() rbeta(n, a, b)
    ldr <<- function(x1, x2) {
      out <- if (a == 1) 0 else (a - 1) * log(x1/x2)
      if (b != 1) out <- out + (b - 1) * log((1 - x1)/(1 - x2))
      out
    }
  }
  type <- "beta"
  environment()
}

#' Create an object representing truncated normal prior distributions
#'
#' @export
#' @param mean scalar or vector mean parameter.
#' @param precision scalar, vector or matrix precision parameter.
#' @param lower lower limit of the truncated interval.
#' @param upper lower limit of the truncated interval.
#' @returns An environment representing the specified prior, for internal use.
pr_truncnormal <- function(mean=0, precision=1, lower=0, upper=Inf) {
  if (!all(lower <= upper)) stop("'lower' should not exceed 'upper'")
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(mean) != c(1L, n))) stop("mean parameter has wrong length")
    if (all(length(precision) != c(1L, n))) stop("precision parameter has wrong length")
    if (all(length(lower) != c(1L, n))) stop("lower parameter has wrong length")
    if (all(length(upper) != c(1L, n))) stop("upper parameter has wrong length")
    stdev <- sqrt(1/precision)
    rprior <<- function() mean + stdev * Crtuvn((lower - mean)/stdev, (upper - mean)/stdev)
  }
  type <- "truncnormal"
  environment()
}

#' Create an object representing exponential prior distributions
#'
#' @export
#' @param scale scalar or vector scale parameter.
#' @returns An environment representing the specified prior, for internal use.
# TODO sample precision instead of variance parameters; the default value
#      for p then leads to inverse Gaussian posterior, which may be faster
#      to sample from
pr_exp <- function(scale=1) {
  if (!all(scale > 0)) stop("scale parameter must be positive")
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(scale) != c(1L, n))) stop("scale parameter has wrong length")
    rprior <<- function() scale * rexp(n)
  }
  type <- "exp"
  environment()
}

#' Create an object representing gamma prior distributions
#'
#' @export
#' @param shape scalar or vector shape parameter.
#' @param rate scalar or vector rate, i.e. inverse scale, parameter.
#' @returns An environment representing the specified prior, for internal use.
pr_gamma <- function(shape=1, rate=1) {
  if (!all(shape > 0 & rate > 0)) stop("shape and rate parameters must be positive")
  n <- NULL
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(shape) != c(1L, n))) stop("scale parameter has wrong length")
    if (all(length(rate) != c(1L, n))) stop("scale parameter has wrong length")
  }
  rprior <- function() rgamma(n, shape, rate)
  type <- "gamma"
  environment()
}

#' Create an object representing Generalised Inverse Gaussian (GIG) prior distributions
#'
#' @export
#' @param a scalar or vector parameter.
#' @param b scalar or vector parameter.
#' @param p scalar or vector parameter.
#' @returns An environment representing the specified prior, for internal use.
pr_gig <- function(a, b, p) {
  if (any(a < 0)) stop("parameter 'a' must be nonnegative")
  if (any(b < 0)) stop("parameter 'b' must be nonnegative")
  if (any(a[p >= 0] == 0)) stop("parameter 'a' should not be 0 when p >= 0")
  if (any(b[p <= 0] == 0)) stop("parameter 'b' should not be 0 when p <= 0")
  n <- NULL
  rprior <- function() stop("please call method 'init' first")
  init <- function(n=1L) {
    n <<- as.integer(n)
    if (all(length(a) != c(1L, n))) stop("parameter 'a' has wrong length")
    if (all(length(b) != c(1L, n))) stop("parameter 'b' has wrong length")
    if (all(length(p) != c(1L, n))) stop("parameter 'p' has wrong length")
    rprior <<- function() Crgig(n, p, a, b)
  }
  type <- "gig"
  environment()
}

#' Create an object representing inverse chi-squared priors
#' with possibly modelled degrees of freedom and scale parameters
#'
#' @export
#' @param df degrees of freedom parameter. This can be a numeric scalar or
#'  vector of length \code{n}, the dimension of the parameter vector.
#'  Alternatively, for a scalar degrees of freedom parameter,
#'  \code{df="modeled"} or \code{df="modelled"} assign a default (gamma) prior
#'  to the degrees of freedom parameter. For more control of this gamma prior a
#'  list can be passed with some of the following components:
#'  \describe{
#'    \item{alpha0}{shape parameter of the gamma distribution}
#'    \item{beta0}{rate parameter of the gamma distribution}
# TODO MH parameters below should not be part of the prior...
#'    \item{proposal}{"RW" for random walk Metropolis-Hastings
#'      or "mala" for Metropolis-adjusted Langevin}
#'    \item{tau}{(starting) scale of Metropolis-Hastings update}
#'    \item{adapt}{whether to adapt the scale of the proposal distribution
#'      during burnin to achieve better acceptance rates.}
#'   }
#' @param scale scalar or vector scale parameter. Alternatively,
#'  \code{scale="modeled"} or \code{scale="modelled"} puts a default
#'  chi-squared prior on the scale parameter. For more control on this
#'  chi-squared prior a list can be passed with some of the following components:
#'  \describe{
#'    \item{df}{degrees of freedom (scalar or vector)}
#'    \item{scale}{scale (scalar or vector)}
#'    \item{common}{whether the modelled scale parameter of the inverse chi-squared
#'      distribution is (a scalar parameter) common to all \code{n} parameters.}
#'  }
#' @returns An environment representing the specified prior, for internal use.
# TODO sample precisions instead of variances
pr_invchisq <- function(df=1, scale=1) {
  if (is.character(df) && any(df == c("modeled", "modelled"))) df <- list()
  if (is.character(scale) && any(scale == c("modeled", "modelled"))) scale <- list()
  n <- rprior <- draw <- NULL
  if (is.list(df)) {
    if (is.null(df[["alpha0"]])) df$alpha0 <- 2  # default in prior Gamma(alpha0, beta0)
    if (is.null(df[["beta0"]])) df$beta0 <- 0.1  # default in prior Gamma(alpha0, beta0)
    # TODO tau, proposal, adapt: move to another object
    if (is.null(df[["tau"]])) df$tau <- 1  # (starting) scale of MH update
    if (is.null(df[["proposal"]])) df$proposal <- "RW"
    if (is.null(df[["adapt"]])) df$adapt <- TRUE
    rprior_df <- draw_df <- NULL
  }
  if (is.list(scale)) {
    scale <- local({
      defaults <- list(df=1, scale=1, common=FALSE)
      if (!all(names(scale) %in% names(defaults))) stop("invalid 'scale' options list")
      modifyList(defaults, scale)
    })
  }
  if (is.list(scale) || (!is.list(scale) && !is.list(df))) psi0 <- NULL
  init <- function(n=1L) {
    n <<- as.integer(n)
    rprior <<- function() {}
    if (is.list(df)) {
      rprior_df <<- function() rgamma(1L, df[["alpha0"]], df[["beta0"]])
      formals(rprior) <<- c(alist(df=), formals(rprior))  # different signature without check NOTE
    } else {
      if (all(length(df) != c(1L, n))) stop("degrees of freedom parameter has wrong length")
      # do not enforce df > 0 to allow improper prior for data scale variance
    }
    if (is.list(scale)) {
      if (all(length(scale[["df"]]) != c(1L, n))) stop("degrees of freedom parameter has wrong length")
      if (all(length(scale[["scale"]]) != c(1L, n))) stop("scale parameter has wrong length")
      if (scale[["common"]] && n == 1L) scale$common <<- FALSE
      psi0 <<- scale[["df"]] / scale[["scale"]]
      if (scale[["common"]]) {
        if (length(scale[["df"]]) != 1L || length(scale[["scale"]]) != 1L) stop("scalar 'df' and 'scale' expected in common scale model")
        rprior <<- rprior |>
          add(quote(scale <- rchisq_scaled(1L, scale[["df"]], psi=psi0))) |>
          add(bquote(1 / rchisq_scaled(.(n), df, scale)))
      } else {
        rprior <<- add(rprior, bquote(draw_betaprime(.(n), 0.5*scale[["df"]], 0.5*df, df/psi0)))
      }
    } else {
      if (all(length(scale) != c(1L, n))) stop("scale parameter has wrong length")
      if (is.list(df)) {
        rprior <<- add(rprior, bquote(1 / rchisq_scaled(.(n), df, scale)))
      } else {
        psi0 <<- df * scale
        rprior <<- add(rprior, bquote(1 / rchisq_scaled(.(n), df, psi=psi0)))
      }
    }
  }
  make_draw <- function() {
    # function draw to sample from full conditional posterior
    # assumes the invchisq prior is for variance parameters of a gaussian distribution
    draw <<- function(df.data, SSR) {}
    if (is.list(df)) {
      draw_df <<- function(df.current, Q.current) {}  # Q.current is current precision parameter (vector)
      switch(df[["proposal"]],
        RW = draw_df <<- add(draw_df, bquote(draw_df_MH_RW(.(as.numeric(n)), df.current, Q.current, df))),
        mala = draw_df <<- add(draw_df, bquote(draw_df_MH_mala(.(as.numeric(n)), df.current, Q.current, df)))
      )
      formals(draw) <<- c(alist(df=), formals(draw))
    }
    if (is.list(scale)) {
      formals(draw) <<- c(formals(draw), alist(Q=))
      # draw 1/kappa from its full conditional posterior
      if (scale[["common"]]) {
        if (is.list(df) || length(df) == 1L)
          draw <<- add(draw, bquote(kappa_inv <- rchisq_scaled(1L, .(n) * df + scale$df, psi = psi0 + sum(df * Q))))
        else
          draw <<- add(draw, quote(kappa_inv <- rchisq_scaled(1L, sum(df) + scale$df, psi = psi0 + sum(df * Q))))
      } else {
        draw <<- add(draw, bquote(kappa_inv <- rchisq_scaled(.(n), df + scale$df, psi = psi0 + df * Q)))
      }
      draw <<- add(draw, quote(psi <- df * kappa_inv))
    } else {
      draw <<- add(draw, quote(psi <- df * scale))
    }
    draw <<- add(draw, bquote(1 / rchisq_scaled(.(n), df + df.data, psi=psi + SSR)))
  }
  type <- "invchisq"
  environment()
}

#' Create an object representing an inverse Wishart prior,
#' possibly with modelled scale matrix
#'
#' @export
#' @param df Degrees of freedom parameter. This should be a scalar numeric value.
#'  Default value is the dimension plus one.
#' @param scale Either a (known) scale matrix, or
#'  \code{scale="modeled"} or \code{scale="modelled"}, which puts default
#'  chi-squared priors on the diagonal elements of the inverse Wishart scale matrix.
#'  For more control on these chi-squared priors a list can be passed with some of the
#'  following components:
#'  \describe{
#'    \item{df}{degrees of freedom (scalar or vector) of the chi-squared distribution(s)}
#'    \item{scale}{scale parameter(s) of the chi-squared distribution(s)}
#'    \item{common}{whether the modelled scale parameter of the inverse chi-squared
#'      distribution is (a scalar parameter) common to all \code{n} diagonal elements.}
#'  }
#' @returns An environment representing the specified prior, for internal use.
#' @references
#'  A. Huang and M.P. Wand (2013).
#'    Simple marginally noninformative prior
#'    distributions for covariance matrices.
#'    Bayesian Analysis 8, 439-452.
pr_invwishart <- function(df=NULL, scale=NULL) {
  if (is_a_matrix(scale)) {
    n <- dim(scale)[1L]
  } else {
    n <- NULL
  }
  if (is.null(df)) {
    if (!is.null(n)) df <- n + 1  # default number of degrees of freedom
  } else {
    if (length(df) != 1L) stop("degrees of freedom parameter for inverse Wishart must be scalar")
    # to allow improper priors, do not enforce df > n - 1
  }
  if (!is.null(scale)) {
    if (is.character(scale) && any(scale == c("modeled", "modelled"))) scale <- list()
    if (is.list(scale)) {  # Huang-Wand prior
      if (is.null(scale$df)) scale$df <- 1
      if (is.null(scale$scale)) scale$scale <- 1
      if (is.null(scale$common)) scale$common <- FALSE  # by default (HW prior)
    }
  }
  init <- function(n=2L) {
    n <<- as.integer(n)
    if (n <= 1L) stop("invalid argument; n must be larger than 1")
    if (is.null(df)) {
      df <<- n + 1  # default number of degrees of freedom
    }
    if (is.null(scale)) {
      scale <<- diag(n)
    } else {
      if (is.list(scale)) {  # Huang-Wand prior
        if (all(length(scale$df) != c(1L, n))) stop("degrees of freedom parameter has wrong length")
        if (all(length(scale$scale) != c(1L, n))) stop("scale parameter has wrong length")
        if (scale$common) {
          if (length(scale$df) != 1L || length(scale$scale) != 1L) stop("scalar 'df' and 'scale' expected in common scale model")
        }
      } else {
        if (!identical(dim(scale), c(n, n))) stop("incompatible scale matrix")
      }
    }
  }
  type <- "invwishart"
  environment()
}
