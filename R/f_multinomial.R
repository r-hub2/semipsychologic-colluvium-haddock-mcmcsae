
#' Specify a multinomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a multinomial sampling distribution. This
#' includes the special case of categorical (multinoulli) data.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for the multinomial distribution is \code{"logit"}.
#' @param n.trial the number of multinomial trials. This can be specified either as
#'  a formula for a variable number of trials, or as a scalar value for a common
#'  number of trials for all units.
#' @param K number of categories for multinomial model; only used for prior predictive sampling.
#' @returns A family object.
f_multinomial <- function(link="logit", n.trial=NULL, K=NULL) {
  family <- "multinomial"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  if (!is.null(K)) {
    K <- as.integer(K)
    if (length(K) != 1L) stop("number of categories 'K' must be a scalar integer")
    if (K < 2L) stop("number of categories 'K' must be at least 2")
  }
  cats <- NULL
  ny0 <- ny <- NULL
  init <- function(data, y=NULL) {
    # response can be either
    # 1. a single vector x (factor, character, boolean or integer, but not numeric!) --> as.factor(x)
    #    --> single trial
    #    "factor"
    # 2. n x (K-1) numeric matrix x, interpreted as proportion of successes; requires n.trial
    #    "fraction"
    # 3. a K-column integer matrix; 1st column is number of successes, 2nd column number of failures
    #    "matrix"
    if (is.null(n.trial)) {
      ny0 <<- 1L
    } else {
      if (is_numeric_scalar(n.trial)) {
        ny0 <<- n.trial
      } else {
        if (!inherits(n.trial, "formula")) stop("'n.trial' must be either a single numeric value (applying to all observations), or a formula")
        ny0 <<- get_var_from_formula(n.trial, data)
        if (all(length(ny0) != c(1L, n_row(data)))) stop("wrong length for number of binomial trials")
      }
      if (anyNA(ny0)) stop("missing(s) in binomial number of trials")
      if (!is.numeric(ny0)) stop("non-numeric number of binomial trials")
      if (any(ny0 < 0)) stop("negative number of binomial trials")
    }
    if (!is.null(y)) {
      if (is.vector(y) || is.factor(y)) {
        # case 1, categorical/multinoulli
        if (is.numeric(y) && any(round(y) != y)) stop("response variable does not seem to be multinomial")
        y <- qF(y)
        cats <<- levels(y)
        y <- model_matrix(~ 0 + y, sparse=FALSE)
      } else {
        if (!is.matrix(y)) stop("unexpected response variable for multinomial family")
        if (any(y < 0)) stop("negative response value(s)")
        if (all(y == 0 | y >= 1)) {
          # case 3
          ny0 <<- dapply(y, sum, MARGIN=1L)  # ny0 will be int if y is (rowSums would convert to double)
          cats <<- dimnames(y)[[2L]]
          if (is.null(cats)) cats <<- as.character(seq_len(ncol(y)))
        } else {
          if (any(y > 1)) stop("ambiguous multinomial response data: contains both fractions and values greater than 1")
          # case 2, assume that only the first K-1 columns are specified; here we require n.trial argument
          cats <<- dimnames(y)[[2L]]
          if (is.null(cats)) cats <<- as.character(seq_len(ncol(y)))
          cats <<- c(cats, "_last_")
          y <- ny0 * y
          y <- cbind(y, ny0 - rowSums(y))  # rowSums converts int to double, but that may be required in samplers.R anyway
        }
      }
      if (!is.null(K) && K != ncol(y)) warn("argument 'K' of f_multinomial differs from number of categories inferred from response vector; 'K' will be ignored")
      K <<- ncol(y)
      # construct ny variable according to stick-breaking representation
      ny <<- if (length(ny0) == 1L) rep.int(ny0, n_row(data)) else ny0
      ny <<- get_sequential_trials(ny, y, K - 1L)
      # return y as a vector minus the part corresponding to the last category
      y <- as.vector(y[, -length(cats)])
    } else {
      if (is.null(K)) stop("for prior multinomial sampling the number of categories must be specified through argument 'K' of f_multinomial")
      cats <<- as.character(seq_len(K))
    }
    if (any(abs(ny0 - round(ny0)) > sqrt(.Machine$double.eps))) {
      warn("one or more non-integral number of trials")
    } else {
      ny0 <<- as.integer(ny0)  # ensures that generated multinomial data is integer!
    }
    if (!is.null(y)) {
      if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) warn("one or more non-integral number of successes")
      if (any(y > ny)) stop("number of successes must not exceed number of trials")  # NB algorithm may still run
    }
    y
  }
  make_llh <- function(y) {
    llh_0 <- sum(binomial_coef(ny, y))  # zero in case of binary data
    function(p) {
      neg_fitted <- -p[["e_"]]
      llh_0 + sum((ny - y) * neg_fitted - ny * log1pexpC(neg_fitted))
    }
  }
  make_llh_i <- function(y) {
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      neg_fitted_i <- -e_i
      rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(ny[i] - y[i], nr) * neg_fitted_i - rep_each(ny[i], nr) * log1pexpC(neg_fitted_i)
    }
  }
  # weights: passed from predict.mcdraws
  #   can be either a numeric scalar, or a vector of length n or nrow(newdata) if the latter is provided
  # long vector format
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata %/% (K - 1L)
      size <- ny0
    } else {
      nn <- nrow(newdata)
      if (is_numeric_scalar(n.trial)) {
        size <- n.trial
      } else if (inherits(n.trial, "formula")) {
        size <- get_var_from_formula(n.trial, newdata)
      } else {
        if (length(ny0) != 1L) stop("number of multinomial trials for prediction cannot be derived")
        size <- ny0
      }
    }
    if (!is.null(weights)) size <- weights * size
    if (any(abs(size - round(size)) > sqrt(.Machine$double.eps))) {
      warn("non-integral values for number of trials are rounded")
      size <- round(size)
    }
    size <- as.integer(size)
    n.out <- nn * (K - 1L)
    loop.range <- 2:(K - 1L)
    function(p, lp) {
      ptilde <- linkinv(lp)
      out <- integer(n.out)
      ind <- seq_len(nn)
      # assume size has length nn or 1
      temp <- rbinom(nn, size=size, prob=ptilde[ind])
      out[ind] <- temp
      for (k in loop.range) {
        size <- size - temp
        ind <- ind + nn
        temp <- rbinom(nn, size=size, prob=ptilde[ind])
        out[ind] <- temp
      }
      out
    }
  }
  # short vector format, only possible for categorica/multinoulli data
  make_rpredictive_cat <- function(newdata, weights=NULL) {
    # assume all ny 1 (or 0) and any weights are ignored (TODO check/warn)
    nn <- n_row(newdata)
    loop.range <- 2:(K - 1L)
    function(p, lp) {
      pSB <- linkinv(lp)
      out <- rep.int(K, nn)  # baseline category
      out[ny0 == 0L] <- NA_integer_
      ind <- seq_len(nn)
      temp <- rbinom(nn, size=ny0, prob=pSB[ind])
      out[temp == 1L] <- 1L
      for (k in loop.range) {
        ny0 <- ny0 - temp
        ind <- ind + nn
        temp <- rbinom(nn, size=ny0, prob=pSB[ind])
        out[temp == 1L] <- k
      }
      out
    }
  }
  environment()
}

#' Get the number of sequential trials for the stick-breaking representation of the multinomial distribution
#'
#' @noRd
#' @param ny vector of total number of trials.
#' @param y matrix of multinomial data.
#' @param ncols the number of columns of y to use.
#' @returns A vector containing the number of sequential trials for the
#'  stick-breaking representation of the multinomial distribution.
get_sequential_trials <- function(ny, y, ncols) {
  out <- matrix(NA_integer_, length(ny), ncols)
  out[, 1L] <- ny
  for (j in seq_len(ncols - 1L)) out[, j + 1L] <- out[, j] - y[, j]
  as.vector(out)
}
