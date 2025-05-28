
#' Specify a binomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a binomial sampling distribution. This
#' includes the special case of binary (Bernoulli) data.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link functions
#'  for the binomial distribution are \code{"logit"} (default) and \code{"probit"}.
#' @param n.trial the number of binomial trials. This can be specified either as
#'  a formula for a variable number of trials, or as a scalar value for a common
#'  number of trials for all units.
#' @returns A family object.
f_binomial <- function(link=c("logit", "probit"), n.trial=NULL) {
  family <- "binomial"
  link <- match.arg(link)
  if (link == "logit")
    linkinv <- make.link(link)$linkinv
  else
    linkinv <- pnorm
  ny <- NA_integer_
  init <- function(data, y=NULL) {
    # y=NULL in case of prior sampling/prediction
    # response can be either
    # 1. a single vector x (factor, character, boolean or integer, but not numeric!) where successes are interpreted as all cases not having the first level of as.factor(x)
    #    --> only binary data, transformed to 0/1 integer vector
    #    "binary"
    # 2. numeric vector x, interpreted as proportion of successes; requires n.trial
    #    "fraction"
    # 3. a two-column integer matrix; 1st column is number of successes, 2nd column number of failures
    #    "matrix"
    if (is.null(n.trial)) {
      ny <<- 1
    } else {
      if (is_numeric_scalar(n.trial)) {
        ny <<- n.trial
      } else {
        if (!inherits(n.trial, "formula")) stop("'n.trial' must be either a single numeric value (applying to all observations), or a formula")
        ny <<- get_var_from_formula(n.trial, data)
        if (all(length(ny) != c(1L, n_row(data)))) stop("wrong length for number of binomial trials")
        if (anyNA(ny)) stop("missing(s) in binomial number of trials")
      }
      if (!is.numeric(ny)) stop("non-numeric number of binomial trials")
      if (any(ny < 0)) stop("negative number of binomial trials")
    }
    if (!is.null(y)) {
      if (is.vector(y) || is.factor(y)) {
        if (is.logical(y)) {
          y <- as.integer(y)
        } else {
          if (is.character(y)) y <- qF(y)
          if (is.factor(y)) {
            if (nlevels(y) > 2L) stop("response factor variable with more than two levels; please use family='multinomial'")
            y <- as.integer(y) - 1L
          } else {
            if (!is.numeric(y)) stop("unexpected response input")
            if (any(y < 0)) stop("negative response value(s)")
            if (!all(y == 0 | y >= 1)) {
              if (any(y > 1)) stop("ambiguous binomial response data: contains both fractions and values greater than 1")
              # case 2, here we require n.trial argument
              y <- ny * y
            }
          }
        }
      } else if (is.matrix(y)) {
        # case 3
        # even though n.trial is ignored in this case for model fitting, it may still be useful for prediction
        if (ncol(y) != 2L) stop("matrix response for binomial family must have 2 columns")
        ny <<- rowSums(y)
        y <- y[, 1L]
      } else {
        stop("unexpected type of response")
      }
      if (link == "probit") {
        if (!allv(ny, 1L)) stop("only binary data (number of trials = 1) supported for binomial probit model")
      } else {
        ny <<- as.numeric(ny)  # both CrPGapprox and BayesLogit::rpg require doubles
      }
    }
    if (any(abs(ny - round(ny)) > sqrt(.Machine$double.eps))) warn("one or more non-integral number of trials")
    if (!is.null(y)) {
      if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) warn("one or more non-integral number of successes")
      if (any(y > ny)) stop("number of successes must not exceed number of trials")  # NB algorithm may still run
    }
    y
  }
  make_llh <- function(y) {
    llh_0 <- sum(binomial_coef(ny, y))  # zero in case of binary data
    if (link == "probit") {
      function(p) {
        sum(pnorm((2*y - 1) * p[["e_"]], log.p=TRUE))
      }
    } else {
      function(p) {
        neg_fitted <- -p[["e_"]]
        llh_0 + sum((ny - y) * neg_fitted - ny * log1pexpC(neg_fitted))
      }
    }
    # faster than dbinom since normalization constant computed only once  
    #sum(dbinom(y, ny, 1 / (1 + exp(neg_fitted)), log=TRUE))
  }
  make_llh_i <- function(y) {
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      neg_fitted_i <- -e_i
      if (link == "logit") {
        if (length(ny) == 1L)  # typically binary regression, ny=1
          rep_each(binomial_coef(ny, y[i]), nr) + rep_each(ny - y[i], nr) * neg_fitted_i - ny * log1pexpC(neg_fitted_i)
        else  # general (negative) binomial regression, ny has length n
          rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(ny[i] - y[i], nr) * neg_fitted_i - rep_each(ny[i], nr) * log1pexpC(neg_fitted_i)
      } else {  # probit
        pnorm(rep_each(1 - 2*y[i], nr) * neg_fitted_i, log.p=TRUE)
      }
    }
  }
  # weights: passed from predict.mcdraws
  #   can be either a numeric scalar, or a vector of length n or nrow(newdata) if the latter is provided
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
      size <- ny
    } else {
      nn <- nrow(newdata)
      if (is.null(n.trial))
        size <- 1
      else if (is_numeric_scalar(n.trial))
        size <- n.trial
      else
        size <- get_var_from_formula(n.trial, newdata)
    }
    if (!is.null(weights)) size <- weights * size
    if (any(abs(size - round(size)) > sqrt(.Machine$double.eps))) {
      warn("non-integral values for number of trials are rounded")
      size <- round(size)
    }
    size <- as.integer(size)
    function(p, lp) rbinom(nn, size, prob=linkinv(lp))
  }
  environment()
}
