#' Specify a binomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a binomial sampling distribution.
#'
#' For the binomial family, the left hand side of the \code{formula} argument of
#' \code{\link{create_sampler}} can be specified in several ways (similar to the
#' options allowed by \code{\link[stats]{glm}} for binomial family):
#' \describe{
#'   \item{1}{as a factor, character, or boolean variable, say \code{y}.
#'     In this case the first level of \code{as.factor(y)} is interpreted as 'failures' and all
#'     other values as 'successes'. This option can only be used for binary data.}
#'   \item{2a}{as a numeric or integer vector with no values strictly between 0 and 1.
#'     The values are then interpreted as numbers of successes. This requires specifying
#'     the number of trials through argument \code{n.trial}, unless the data is
#'     binary (values 0 and 1 only). For non-integer values a warning is issued.}
#'   \item{2b}{as a numeric vector with values between 0 and 1. The values are then interpreted
#'     as the proportion of the number of successes. This requires specifying the number of trials
#'     through argument \code{n.trial}.}
#'   \item{3}{as a two-column integer matrix. The first column is interpreted as the number of
#'     successes, the second column as the number of failures.}
#' }
#'
#' @examples
#' y <- c(TRUE, FALSE, TRUE)
#' sampler <- create_sampler(y ~ 1, family=f_binomial())
#' # logistic binary regression is the default, and may also be
#' # specified as family="binomial"
#' sim <- MCMCsim(sampler, n.chain=2, n.iter=100, burnin=50, verbose=FALSE)
#' summary(predict(sim, newdata=data.frame(id=1)))
#'
#' y <- c(0, 0, 1, 1, 0, 0, 0, 0, 1)
#' sampler <- create_sampler(y ~ 1, family=f_binomial(link="probit"))
#' sim <- MCMCsim(sampler, n.chain=2, n.iter=100, burnin=50, verbose=FALSE)
#' summary(predict(sim, newdata=data.frame(id=1)))
#'
#' sampler <- create_sampler(c(0.2, 0.3, 0.5, 0.0) ~ 1, family=f_binomial(n.trial=10))
#' # is interpreted as
#' sampler <- create_sampler(c(2, 3, 5, 0) ~ 1, family=f_binomial(n.trial=10))
#' sim <- MCMCsim(sampler, n.chain=2, n.iter=100, burnin=50, verbose=FALSE)
#' summary(predict(sim, newdata=data.frame(id=1)))
#'
#' n <- 1000
#' dat <- data.frame(
#'   x = runif(n),
#'   g = factor(sample(1:10, n, replace=TRUE)),
#'   trials = sample(1:5, n, replace=TRUE)
#' )
#' v <- rnorm(10)
#' dat$y <- rbinom(n, size=dat$trials, prob=1 / (1 + exp(-(1 - dat$x + v[dat$g]))))
#' sampler <- create_sampler(y ~ x + (1|g), data=dat, family=f_binomial(n.trial = ~ trials))
#' # alternatively:
#' sampler <- create_sampler(cbind(y, trials - y) ~ x + (1|g), data=dat, family="binomial")
#' sim <- MCMCsim(sampler, store.all=TRUE, burnin=200, n.iter=300, n.chain=2, verbose=FALSE)
#' summ <- summary(sim)
#' plot(v, summ$gen2[, "Mean"]); abline(0, 1)
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link functions
#'  for the binomial distribution are \code{"logit"} (default) and \code{"probit"}.
#' @param n.trial the number of binomial trials. This can be specified either as
#'  a formula for a variable number of trials, or as a scalar value for a common
#'  number of trials for all units.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{binomial_control}}.
#' @returns A family object.
f_binomial <- function(link=c("logit", "probit"), n.trial=NULL,
                       control=binomial_control()) {
  link <- match.arg(link)
  control <- check_binomial_control(control)
  list(family="binomial", link=link, n.trial=n.trial, control=control,
       `_raw_`=TRUE)
}

ff_binomial <- function(link=c("logit", "probit"), n.trial=NULL,
                        control=binomial_control(),
                        sm, data, y=NULL, famid=NULL, sub=NULL) {
  family <- "binomial"
  if (link == "logit") {
    linkinv <- make.link(link)$linkinv
    f_mean <- function(eta) ny / (1 + exp(-eta))
  } else {
    linkinv <- pnorm
    f_mean <- function(eta) ny * pnorm(eta)
  }
  if (link == "probit" && !is.null(sub)) stop("probit link as part of multi-response family not supported")
  e.is.res <- FALSE
  sigma.fixed <- TRUE
  modeled.Q <- link != "probit"
  prior.only <- is.null(y)
  multifam <- !is.null(sub)
  sc <- sm[["control"]]
  store_default <- function() NULL
  n <- n_row(data)
  Q0.type <- "unit"
  Q0 <- CdiagU(n)
  scale.e <- 2.5
  scale.sigma <- scale.e
  if (is.null(n.trial)) {
    ny <- 1
  } else {
    if (is_numeric_scalar(n.trial)) {
      ny <- n.trial
    } else {
      if (!inherits(n.trial, "formula")) stop("'n.trial' must be either a single numeric value (applying to all observations), or a formula")
      ny <- get_var_from_formula(n.trial, data)
      if (all(length(ny) != c(1L, n_row(data)))) stop("wrong length for number of binomial trials")
      if (anyNA(ny)) stop("missing(s) in binomial number of trials")
    }
    if (!is.numeric(ny)) stop("non-numeric number of binomial trials")
    if (any(ny < 0)) stop("negative number of binomial trials")
  }
  if (!prior.only) {
    if (link == "probit") {
      if (sc[["single.block"]])
        Q_e <- function(p) p[["z_"]]
      else
        Q_e <- function(p) p[["z_"]] - p[["e_"]]
      draw <- function(p) {
        p$llh_ <- llh(p)
        p$z_ <- CrTNprobit(p[["e_"]], y)
        p
      }
      start <- function(p) {
        p$z_ <- check_and_get(p, "z_", n, \() CrTNprobit(p[["e_"]], y))
        p
      }
    } else {  # logit
      if (sc[["single.block"]]) {
        Q_e <- function(p) y_shifted
      } else {
        if (multifam)
          Q_e <- function(p) y_shifted - p[["Q_"]][sub] * p[["e_"]][sub]
        else
          Q_e <- function(p) y_shifted - p[["Q_"]] * p[["e_"]]
      }
      if (!control[["PG.approx"]]) {
        ny[ny == 0] <- .Machine$double.eps  # prevent generation of NAs by BayesLogit::rpg
      }
      rPolyaGamma <- get_PG_sampler(n, control[["PG.approx"]], control[["PG.approx.m"]])
      if (multifam) {
        draw <- function(p) {
          p$llh_ <- p$llh_ + llh(p)
          p$Q_[sub] <- rPolyaGamma(ny, p[["e_"]][sub])
          p
        }
      } else {
        draw <- function(p) {
          p$llh_ <- llh(p)
          p[["Q_"]] <- rPolyaGamma(ny, p[["e_"]])
          p
        }
        start <- function(p) {
          p$Q_ <- check_and_get(p, "Q_", n, \() rPolyaGamma(ny, p[["e_"]]), pos=TRUE)
          p
        }
      }
    }
    if (is.vector(y) || is.factor(y)) {
      # case 1: logical, factor or character data
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
            # case 2b, here we require n.trial argument
            y <- ny * y
          }
          # case 2a: integer data
        }
      }
    } else if (is.matrix(y)) {
      # case 3
      # even though n.trial is ignored in this case for model fitting, it may still be useful for prediction
      if (ncol(y) != 2L) stop("matrix response for binomial family must have 2 columns")
      ny <- rowSums(y)
      y <- y[, 1L]
    } else {
      stop("unexpected type of response")
    }
    if (link == "probit") {
      if (!allv(ny, 1L)) stop("only binary data (number of trials = 1) supported for binomial probit model")
    } else {
      ny <- as.numeric(ny)  # both CrPGapprox and BayesLogit::rpg require doubles
      y_shifted <- y - 0.5 * ny
    }
  }
  if (any(abs(ny - round(ny)) > .tol)) warn("one or more non-integral number of trials")
  if (!prior.only) {
    if (any(abs(y - round(y)) > .tol)) warn("one or more non-integral number of successes")
    if (any(y > ny)) stop("number of successes must not exceed number of trials")  # NB algorithm may still run
    if (!is.null(sc[["CG"]]) || sc[["cMVN.sampler"]]) {
      # set up a function that multiplies by L Chol factor of Q, for sampling from N(., Q)
      if (link == "probit") {
        cholQ <- build_chol(CdiagU(n))
        drawMVNvarQ <- function(p) cholQ$Ltimes(Crnorm(n), transpose=FALSE)
      } else {
        cholQ <- build_chol(runif(n, 0.9, 1.1))
        drawMVNvarQ <- function(p) {
          cholQ$update(if (multifam) p[["Q_"]][sub] else p[["Q_"]])
          cholQ$Ltimes(Crnorm(n), transpose=FALSE)
        }
      }
    }
    if (link == "probit") {
      llh <- function(p) sum(pnorm((2*y - 1) * p[["e_"]], log.p=TRUE))
      llh_i <- function(draws, i, e_i) {
        nr <- dim(e_i)[1L]
        pnorm(rep_each(2*y[i] - 1, nr) * e_i, log.p=TRUE)
      }
    } else {
      llh_0 <- sum(binomial_coef(ny, y))  # zero in case of binary data
      if (all(ny == 1))
        llh <- function(p) {
          e_ <- if (multifam) p[["e_"]][sub] else p[["e_"]]
          llh_0 + sum(y * e_ - log1pexpC(e_))
        }
      else
        llh <- function(p) {
          e_ <- if (multifam) p[["e_"]][sub] else p[["e_"]]
          llh_0 + sum(y * e_ - ny * log1pexpC(e_))
        }
      llh_i <- function(draws, i, e_i) {
        nr <- dim(e_i)[1L]
        if (length(ny) == 1L)  # typically binary regression, ny=1
          rep_each(binomial_coef(ny, y[i]), nr) + rep_each(y[i], nr) * e_i - ny * log1pexpC(e_i)
        else  # general (negative) binomial regression, ny has length n
          rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(y[i], nr) * e_i - rep_each(ny[i], nr) * log1pexpC(e_i)
      }
    }
    # faster than dbinom since normalisation constant computed only once
    #sum(dbinom(y, ny, 1 / (1 + exp(neg_fitted)), log=TRUE))
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
    if (any(abs(size - round(size)) > .tol)) {
      warn("non-integral values for number of trials are rounded")
      size <- round(size)
    }
    size <- as.integer(size)
    function(p, lp) rbinom(nn, size, prob=linkinv(lp))
  }
  rm(data, sm)
  environment()
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
#' @param probit.HaarPXDA only used for binomial models with probit link. If supported,
#'  \code{probit.HaarPXDA=TRUE} and the Haar PX-DA sandwich step is added to the Albert-Chib
#'  data augmentation scheme for probit multilevel models. This will usually result in
#'  a faster mixing MCMC algorithm. Currently supported when all coefficients are sampled
#'  in a single Gibbs block and all coefficients' prior means equal zero.
#' @returns A list with computational options.
binomial_control <- function(PG.approx=TRUE, PG.approx.m=-2L, probit.HaarPXDA=TRUE) {
  list(PG.approx=PG.approx, PG.approx.m=PG.approx.m, probit.HaarPXDA=probit.HaarPXDA)
}

check_binomial_control <- function(control) {
  if (is.null(control)) return(binomial_control())
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
  if (!is_logical_scalar(control[["probit.HaarPXDA"]])) stop("'probit.HaarPXDA' must be TRUE or FALSE")
  control
}
