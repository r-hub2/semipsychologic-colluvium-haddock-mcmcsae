
#' Specify a Poisson sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Poisson sampling distribution.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed
#'  link function for the Poisson distribution is \code{"log"}.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{poisson_control}}.
#' @returns A family object.
f_poisson <- function(link="log", control=poisson_control()) {
  family <- "poisson"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  control <- check_poisson_control(control)
  # NB size is only used as an internal negbinomial parameter, and only for model fitting
  log.shape <- log(control$nb.shape)
  init <- function(data, y=NULL) {
    if (!is.null(y)) {
      if (!is.numeric(y)) stop("non-numeric target value not allowed in case of Poisson sampling distribution")
      if (any(y < 0)) warn("negative response value(s)")
      if (any(abs(round(y) - y) > sqrt(.Machine$double.eps))) warn("non-integral values modelled by Poisson sampling distribution")
    }
    y
  }
  make_llh <- function(y) {
    llh_0 <- -sum(lgamma(y + 1))
    # Poisson log-likelihood, first remove internal offset from linear predictor p[["e_"]]
    function(p) {
      eta <- p[["e_"]] + log.shape
      llh_0 + sum(y * eta - exp(eta))
    }
  }
  make_llh_i <- function(y) {
    llh_0_i <- -lgamma(y + 1)
    # NB e_i must be the linear predictor excluding the internal offset
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      rep_each(llh_0_i[i], nr) + rep_each(y[i], nr) * e_i - exp(e_i)
    }
  }
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
    } else {
      nn <- nrow(newdata)
    }
    if (is.null(weights))
      function(p, lp) rpois(nn, lambda=exp(lp))
    else
      function(p, lp) rpois(nn, lambda=weights * exp(lp))
  }
  environment()
}

#' Set computational options for the sampling algorithms
#'
#' @export
#' @param nb.shape shape parameter of the negative binomial distribution used
#'  internally to approximate the Poisson distribution. This should be set to a relatively
#'  large value (default is 100), corresponding to negligible overdispersion, to obtain a
#'  good approximation to the Poisson sampling distribution. However, note that very large
#'  values may cause slow MCMC exploration of the posterior distribution.
#' @returns A list with computational options for the sampling algorithm.
poisson_control <- function(nb.shape=100) {
  list(nb.shape=nb.shape)
}

check_poisson_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- poisson_control()
  w <- which(!(names(control) %in% names(defaults)))
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$nb.shape <- as.numeric(control$nb.shape)
  if (length(control$nb.shape) != 1L || is.na(control$nb.shape) || control$nb.shape <= 0)
    stop("'nb.shape' must be a positive numerical scalar")
  control
}
