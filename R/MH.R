#' Set options for Metropolis-Hastings sampling
#'
#' @export
#' @param type a character string defining the proposal distribution.
#'  Among supported types are random walk proposals "RWTN",
#'  "RWN" and "RWLN" with truncated normal, normal and log-normal
#'  proposal distributions. Other choices correspond to independence
#'  proposals: "TN" for a truncated normal proposal, "unif" for a
#'  uniform proposal, and "beta" and "gamma" for specific beta and
#'  gamma proposal distributions. Not all types are supported for
#'  a particular parameter; see the specific help of the function
#'  defining the model component of interest to see which proposal
#'  distribution types are supported.
#' @param adaptive in case of the random walk "RWTN" or "RWN"
#'  proposals, whether the scale parameter is adapted based on
#'  acceptance rates during the burnin phase of the MCMC simulation.
#'  The default is \code{TRUE} in these cases.
#' @param scale in case of the "RWTN" proposal, the (initial) scale of the
#'  distribution.
#' @param ... additional parameters depending on the proposal type.
#'  Supported arguments are 'l' and 'u' to pass the lower and upper
#'  limits of uniform or random walk truncated normal proposals
#'  (defaults l=0 and u=1), and 'a' and 'b' to pass the shape parameters
#'  of a beta proposal distribution (defaults a = b = 0.5).
#' @returns An environment with variables and methods for Metropolis-Hastings
#'  sampling, for use by other package functions.
set_MH <- function(type = "RWTN", scale = 0.025, adaptive = NULL, ...) {
  type <- match.arg(type, c("RWTN", "RWN", "RWLN", "TN", "unif", "beta", "gamma"))
  scale <- as.numeric(scale)
  if (length(scale) != 1L || scale <= 0) stop("'scale' must be a single positive value")
  if (is.null(adaptive)) {
    adaptive <- any(type == c("RWTN", "RWN", "RWLN"))
  } else {
    if (!is.logical(adaptive) || length(adaptive) != 1L)
      stop("argument 'adaptive' must be a single logical value, or 'NULL' for a default choice")
    if (adaptive && all(type != c("RWTN", "RWN", "RWLN")))
      stop("adaptive=TRUE only supported for 'RWTN', 'RWN' and 'RWLN' proposals")
  }
  extra.args <- list(...)
  if (!all(names(extra.args) %in% c("l", "u", "a", "b")))
    stop("invalid argument(s) for set_MH: ", paste(setdiff(names(extra.args), c("l", "u", "a", "b")), collapse=", "))
  if (any(type == c("RWTN", "unif"))) {
    if (any("l" == names(extra.args))) {
      l <- as.numeric(extra.args[["l"]])
      if (length(l) != 1L) stop("lower bound 'l' must be a single numeric value")
    } else {
      l <- NULL  # default can vary e.g. -1 for AR1, 0 for bym2
    }
    if (any("u" == names(extra.args))) {
      u <- as.numeric(extra.args[["u"]])
      if (length(u) != 1L) stop("upper bound 'u' must be a single numeric value")
    } else {
      u <- NULL
    }
  }
  switch(type,
    RWTN = {
      # random walk truncated normal proposal
      propose <- function(x) x + scale * Crtuvn((l - x)/scale, (u - x)/scale)
      log.ar_prop <- function(x.star, x)
        log(
          (pnorm((u - x)/scale) - pnorm((l - x)/scale)) /
          (pnorm((u - x.star)/scale) - pnorm((l - x.star)/scale))
        )
    },
    RWN = {
      # random walk normal proposal
      propose <- function(x) x + scale * rnorm(1L)
      log.ar_prop <- function(x.star, x) 0
    },
    RWLN = {
      propose <- function(x) x * exp(scale * rnorm(1L))
      log.ar_prop <- function(x.star, x) log(x.star/x)
    },
    TN = {
      # independence truncated normal proposal
      # this is currently specific to AR1 parameter,
      # so methods propose and log.ar_prop are not defined here
    },
    unif = {
      # uniform proposal
      propose <- function(x) runif(1L, l, u)
      log.ar_prop <- function(x.star, x) 0
    },
    beta = {
      # beta proposal
      if (any("a" == names(extra.args))) {
        a <- as.numeric(extra.args[["a"]])
        if (length(a) != 1L) stop("argument 'a' of beta proposal must be a single numeric value")
        if (a <= 0) stop("argument 'a' of beta proposal must be positive")
      } else {
        a <- 0.5
      }
      if (any("b" == names(extra.args))) {
        b <- as.numeric(extra.args[["b"]])
        if (length(b) != 1L) stop("argument 'b' of beta proposal must be a single numeric value")
        if (b <= 0) stop("argument 'b' of beta proposal must be positive")
      } else {
        b <- 0.5
      }
      propose <- function(x) rbeta(1L, a, b)
      log.ar_prop <- function(x.star, x) {
        out <- if (a == 1) 0 else (a - 1) * log(x/x.star)
        if (b != 1) out <- out + (b - 1) * log((1 - x)/(1 - x.star))
        out
      }
    },
    gamma = {
      # independence gamma proposal
      # this is currently specific to gamma family shape parameter,
      # so methods propose and log.ar_prop are not defined here
    },
    stop("MH proposal type should be one of 'RWTN', 'unif' or 'beta'")
  )
  rm(extra.args)
  MH_accept <- function(x.star, x, log.ar.post)
    log(runif(1L)) < log.ar.post + log.ar_prop(x.star, x)
  if (adaptive) {
    # adaptation based on acceptance rates, called every 50 iterations (cf. Roberts, Rosenthal)
    adapt <- function(ar) {
      if (ar < .2)
        scale <<- scale * runif(1L, 0.6, 0.9)
      else if (ar > .7)
        scale <<- scale * runif(1L, 1.1, 1.5)
    }
  }
  environment()
}


###############################
# degrees of freedom parameter

# RW proposal on log(df) with scale tau
draw_df_MH_RW <- function(r, df.current, Q.current, df.mod) {
  # draw new degrees of freedom value from proposal
  df.star <- exp(rnorm(1L, sd=df.mod$tau)) * df.current
  # compute log-acceptance-ratio
  log.ar <- r * (  0.5 * df.star * log(0.5 * df.star) - 0.5 * df.current * log(0.5 * df.current)
                   + lgamma(0.5 * df.current) - lgamma(0.5 * df.star) ) +
    df.mod$alpha0 * log(df.star / df.current) +
    (df.current - df.star) * (df.mod$beta0 + 0.5 * sum(Q.current - log(Q.current)))
  if (log(runif(1L)) < log.ar) df.star else df.current
}

# random walk with drift: Metropolis Adjusted Langevin Algorithm
draw_df_MH_mala <- function(r, df.current, Q.current, df.mod) {
  # compute mala drift
  temp <- df.mod$beta0 + 0.5 * sum(Q.current - log(Q.current))
  drift <- df.mod$alpha0 + 0.5 * r * df.current * (1 - log(0.5 * df.current) - digamma(0.5 * df.current)) - df.current * temp
  # draw new degrees of freedom value from mala proposal
  df.star <- exp(df.mod$tau * rnorm(1L, mean=0.5 * drift)) * df.current
  # compute log-acceptance-ratio
  log.ar <- r * ( 0.5 * df.star * log(0.5 * df.star) - 0.5 * df.current * log(0.5 * df.current)
                  + lgamma(0.5 * df.current) - lgamma(0.5 * df.star) ) +
    + (df.current - df.star) * temp + df.mod$alpha0 * log(df.star/df.current)
  if (log(runif(1L)) < log.ar) df.star else df.current
}
