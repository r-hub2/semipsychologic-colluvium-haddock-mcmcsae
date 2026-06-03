
#' Specify a multinomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a multinomial sampling distribution. This
#' includes the special case of categorical (multinoulli) data.
#' 
#' For the multinomial family, the left hand side of the \code{formula} argument of
#' \code{\link{create_sampler}} can be specified in one of the following ways:
#' \describe{
#'   \item{1}{as a single factor, character, boolean or integer variable, say \code{y}.
#'     The categories then correspond to the levels of \code{as.factor(y)}. This option
#'     can only be used for categorical data.}
#'   \item{2}{as a n x (K-1) numeric matrix with values between 0 and 1, where n is the
#'     number of observations and K the number of categories. The values are interpreted
#'     as proportions of observations in each category. This requires specifying the number
#'     of multinomial trials through argument \code{n.trial}.}
#'   \item{3}{a K-column integer matrix, where K is the number of categories, each
#'     column contaning the number of 'successes' for the corresponding category.}
#' }
#'
#' @examples
#' y <- factor(sample(c("a", "b", "c"), 800, prob=c(0.3, 0.5, 0.2), replace=TRUE))
#' sampler <- create_sampler(y ~ 0 + cat_, family=f_multinomial())
#' sim <- MCMCsim(sampler, n.chain=2, burnin=200, n.iter=300, verbose=FALSE)
#' summary(sim)
#' summary(predict(sim, newdata=data.frame(id=1:5)))
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for the multinomial distribution is \code{"logit"}.
#' @param n.trial the number of multinomial trials. This can be specified either as
#'  a formula for a variable number of trials, or as a scalar value for a common
#'  number of trials for all units.
#' @param K number of categories for multinomial model; only used for prior predictive sampling.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{multinomial_control}}.
#' @returns A family object.
f_multinomial <- function(link="logit", n.trial=NULL, K=NULL,
                          control=multinomial_control()) {
  link <- match.arg(link)
  control <- check_multinomial_control(control)
  list(family="multinomial", link=link, n.trial=n.trial, K=K,
       control=control, `_raw_`=TRUE)
}

ff_multinomial <- function(link="logit", n.trial=NULL, K=NULL,
                           control=multinomial_control(),
                           sm, data, y=NULL) {
  family <- "multinomial"
  linkinv <- make.link(link)$linkinv
  e.is.res <- FALSE
  sigma.fixed <- TRUE
  modeled.Q <- TRUE
  prior.only <- is.null(y)
  n <- n_row(data)
  Q0.type <- "unit"
  if (!is.null(K)) {
    K <- as.integer(K)
    if (length(K) != 1L) stop("number of categories 'K' must be a scalar integer")
    if (K < 2L) stop("number of categories 'K' must be at least 2")
  }
  scale.e <- 2.5
  scale.sigma <- scale.e
  sc <- sm[["control"]]
  store_default <- function() NULL
  if (is.null(n.trial)) {
    ny0 <- 1L
  } else {
    if (is_numeric_scalar(n.trial)) {
      ny0 <- n.trial
    } else {
      if (!inherits(n.trial, "formula")) stop("'n.trial' must be either a single numeric value (applying to all observations), or a formula")
      ny0 <- get_var_from_formula(n.trial, data)
      if (all(length(ny0) != c(1L, n))) stop("wrong length for number of multinomial trials")
    }
    if (anyNA(ny0)) stop("missing(s) in multinomial number of trials")
    if (!is.numeric(ny0)) stop("non-numeric number of multinomial trials")
    if (any(ny0 < 0)) stop("negative number of multinomial trials")
  }
  if (prior.only) {
    if (is.null(K)) stop("for prior multinomial sampling the number of categories must be specified through argument 'K' of f_multinomial")
    cats <- as.character(seq_len(K))
  } else {
    if (is.vector(y) || is.factor(y)) {
      # case 1, categorical/multinoulli
      if (is.numeric(y) && any(round(y) != y)) stop("response variable does not seem to be multinomial")
      y <- qF(y)
      cats <- levels(y)
      y <- model_matrix(~ 0 + y, sparse=FALSE)
    } else {
      if (!is.matrix(y)) stop("unexpected response variable for multinomial family")
      if (any(y < 0)) stop("negative response value(s)")
      if (all(y == 0 | y >= 1)) {
        # case 3
        ny0 <- dapply(y, sum, MARGIN=1L)  # ny0 will be int if y is (rowSums would convert to double)
        cats <- dimnames(y)[[2L]]
        if (is.null(cats)) cats <- as.character(seq_len(ncol(y)))
      } else {
        if (any(y > 1)) stop("ambiguous multinomial response data: contains both fractions and values greater than 1")
        # case 2, assume that only the first K-1 columns are specified; here we require n.trial argument
        cats <- dimnames(y)[[2L]]
        if (is.null(cats)) cats <- as.character(seq_len(ncol(y)))
        cats <- c(cats, "_last_")
        y <- ny0 * y
        y <- cbind(y, ny0 - rowSums(y))  # rowSums converts int to double, but that may be required in samplers.R anyway
      }
    }
    if (!is.null(K) && K != ncol(y)) warn("argument 'K' of f_multinomial differs from number of categories inferred from response vector; 'K' will be ignored")
    K <- ncol(y)
    # construct ny variable according to stick-breaking representation
    ny <- if (length(ny0) == 1L) rep.int(ny0, n_row(data)) else ny0
    ny <- get_sequential_trials(ny, y, K - 1L)
    # remove from y the part corresponding to the last category
    y <- as.vector(y[, -length(cats)])
  }
  if (any(abs(ny0 - round(ny0)) > .tol)) {
    warn("one or more non-integral number of trials")
  } else {
    ny0 <- as.integer(ny0)  # ensures that generated multinomial data is integer!
  }
  Km1 <- K - 1L
  n0 <- n
  n <- n0 * Km1
  Q0 <- CdiagU(n)
  if (!prior.only) {
    y_shifted <- y - 0.5 * ny
    if (sc[["single.block"]])
      Q_e <- function(p) y_shifted
    else
      Q_e <- function(p) y_shifted - p[["Q_"]] * p[["e_"]]
    rPolyaGamma <- get_PG_sampler(n, control[["PG.approx"]], control[["PG.approx.m"]])
    draw <- function(p) {
      p$llh_ <- llh(p)
      p$Q_ <- rPolyaGamma(ny, p[["e_"]])
      p
    }
    start <- function(p) {
      p$Q_ <- check_and_get(p, "Q_", n, \() rPolyaGamma(ny, p[["e_"]]), pos=TRUE)
      p
    }
    if (!is.null(sc[["CG"]]) || sc[["cMVN.sampler"]]) {
      # set up a function that multiplies by L Chol factor of Q, for sampling from N(., Q)
      cholQ <- build_chol(runif(n, 0.9, 1.1))
      # draw from MVN with variance(!) Q
      drawMVNvarQ <- function(p) {
        cholQ$update(p[["Q_"]])
        cholQ$Ltimes(Crnorm(n), transpose=FALSE)
      }
    }
    if (any(abs(y - round(y)) > .tol)) warn("one or more non-integral number of successes")
    if (any(y > ny)) stop("number of successes must not exceed number of trials")  # NB algorithm may still run

    llh_0 <- sum(binomial_coef(ny, y))  # zero in case of binary data
    if (all(ny == 1)) {
      llh <- function(p) llh_0 + sum(y * p[["e_"]] - log1pexpC(p[["e_"]]))
      llh_i <- function(draws, i, e_i) {
        nr <- dim(e_i)[1L]
        rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(y[i], nr) * e_i - log1pexpC(e_i)
      }
    } else {
      llh <- function(p) llh_0 + sum(y * p[["e_"]] - ny * log1pexpC(p[["e_"]]))
      llh_i <- function(draws, i, e_i) {
        nr <- dim(e_i)[1L]
        rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(y[i], nr) * e_i - rep_each(ny[i], nr) * log1pexpC(e_i)
      }
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
    if (any(abs(size - round(size)) > .tol)) {
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
  rm(data, sm)
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

#' Set computational options for the sampling algorithms
#'
#' @export
#' @param PG.approx whether Polya-Gamma draws for logistic binomial models are
#'  approximated by a hybrid gamma convolution approach. If not, \code{BayesLogit::rpg}
#'  is used, which is exact for some values of the shape parameter.
#' @param PG.approx.m if \code{PG.approx=TRUE}, the number of explicit gamma draws in the
#'  sum-of-gammas representation of the Polya-Gamma distribution. The remainder (infinite)
#'  convolution is approximated by a single moment-matching gamma draw. Special values are:
#'  \code{-2L} for a default choice depending on the value of the shape parameter
#'  balancing performance and accuracy, \code{-1L} for a moment-matching normal approximation,
#'  and \code{0L} for a moment-matching gamma approximation.
#' @returns A list with computational options.
multinomial_control <- function(PG.approx=TRUE, PG.approx.m=-2L) {
  list(PG.approx=PG.approx, PG.approx.m=PG.approx.m)
}

check_multinomial_control <- function(control) {
  if (is.null(control)) return(multinomial_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- binomial_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  if (!is_logical_scalar(control[["PG.approx"]])) stop("'PG.approx' must be TRUE or FALSE")
  control$PG.approx.m <- as.integer(control[["PG.approx.m"]])
  if (!length(control[["PG.approx.m"]])) stop("unexpected input for 'PG.approx.m'")
  if (any(is.na(control[["PG.approx.m"]]) | control[["PG.approx.m"]] < -2L | is.infinite(control[["PG.approx.m"]])))
    stop("'PG.approx.m' value(s) out of range or missing")
  control
}
