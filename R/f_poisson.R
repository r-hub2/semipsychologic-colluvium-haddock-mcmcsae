#' Specify a Poisson sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Poisson sampling distribution.
#'
#' @examples
#' \dontrun{
#' n <- 3000
#' m <- 25
#' dat <- data.frame(
#'   x = rnorm(n),
#'   g = factor(sample(1:m, n, replace=TRUE), levels=1:m)
#' )
#' v <- rnorm(m, sd=0.6)
#' dat$y <- rpois(n, lambda=exp(with(dat, 1 - 0.5*x + v[g])))
#'
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") +       # fixed effects
#'       gen(factor = ~ g, name="v"),  # random intercepts
#'   data=dat, family="poisson"
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' waic(sim)
#' summary(sim)
#' 
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), c(1, -0.5))
#' bayesplot::mcmc_recover_intervals(as.array(sim$v_sigma), 0.6)
#' yrep <- predict(sim, iters=sample(1:1000, 10))
#' bayesplot::pp_check(dat$y, as.matrix(yrep), bayesplot::ppc_dens_overlay)
#' }
#'
#' @export
#' @param link the name of a link function. Currently the only allowed
#'  link function for the Poisson distribution is \code{"log"}.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{poisson_control}}.
#' @returns A family object.
f_poisson <- function(link="log", control=poisson_control()) {
  link <- match.arg(link)
  control <- check_poisson_control(control)
  lnb <- f_negbinomial(link="log", shape.vec = ~ 1,
    inv.shape.prior=pr_fixed(value=1/control[["nb.shape"]]),
    control=negbinomial_control(PG.approx=control[["PG.approx"]],
      PG.approx.m=control[["PG.approx.m"]]
    )
  )
  lnb <- lnb[-which(names(lnb) %in% c("family", "link"))]
  names(lnb)[names(lnb) == "control"] <- "nb.control"
  c(list(family="poisson", link=link, control=control), lnb)
}

ff_poisson <- function(link="log", control=poisson_control(),
                       shape.vec, inv.shape.prior, nb.control,
                       sm, data, y=NULL, famid=NULL, sub=NULL) {
  family <- "poisson"
  linkinv <- make.link(link)$linkinv
  f_mean <- function(eta) exp(eta)
  # NB shape0 is only used as an internal negbinomial parameter, and only for model fitting
  shape0 <- control[["nb.shape"]]
  log.shape <- log(shape0)
  f <- ff_negbinomial(link="log", shape.vec=shape.vec,
    inv.shape.prior=inv.shape.prior, control=nb.control,
    sm=sm, data=data, y=y, famid=famid, sub=sub
  )
  scale.e <- f[["scale.e"]]
  scale.sigma <- f[["scale.sigma"]]
  e.is.res <- f[["e.is.res"]]
  sigma.fixed <- f[["sigma.fixed"]]
  modeled.Q <- f[["modeled.Q"]]
  multifam <- f[["multifam"]]
  prior.only <- is.null(y)
  n <- f[["n"]]
  Q0.type <- f[["Q0.type"]]
  Q0 <- f[["Q0"]]
  sc <- sm[["control"]]
  store_default <- function() NULL
  if (!prior.only) {
    y_shifted <- f[["y_shifted"]]
    Q_e <- f[["Q_e"]]
    rPolyaGamma <- f[["rPolyaGamma"]]
    draw <- f[["draw"]]
    start <- f[["start"]]
    if (is.function(f[["drawMVNvarQ"]])) {
      cholQ <- f[["cholQ"]]
      drawMVNvarQ <- f[["drawMVNvarQ"]]
    }
    llh0_i <- -lgamma(y + 1)
    llh0 <- sum(llh0_i)
    # Poisson log-likelihood, first remove internal offset from linear predictor p[["e_"]]
    llh <- function(p) {
      eta <- if (multifam) p[["e_"]][sub] + log.shape else p[["e_"]] + log.shape
      llh0 + sum(y * eta - exp(eta))
    }
    assign("llh", llh, environment(draw))  # draw was defined by ff_negbinomial so we should replace its llh function
    # NB e_i must be the linear predictor excluding the internal offset
    llh_i <- function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      rep_each(llh0_i[i], nr) + rep_each(y[i], nr) * e_i - exp(e_i)
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
  rm(data)
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
#' @param PG.approx whether Polya-Gamma draws for logistic binomial models are
#'  approximated by a hybrid gamma convolution approach. If not, \code{BayesLogit::rpg}
#'  is used, which is exact for some values of the shape parameter.
#' @param PG.approx.m if \code{PG.approx=TRUE}, the number of explicit gamma draws in the
#'  sum-of-gammas representation of the Polya-Gamma distribution. The remainder (infinite)
#'  convolution is approximated by a single moment-matching gamma draw. Special values are:
#'  \code{-2L} for a default choice depending on the value of the shape parameter
#'  balancing performance and accuracy, \code{-1L} for a moment-matching normal approximation,
#'  and \code{0L} for a moment-matching gamma approximation.
#' @returns A list with computational options for the sampling algorithm.
poisson_control <- function(nb.shape=100, PG.approx=TRUE, PG.approx.m=-2L) {
  list(nb.shape=nb.shape, PG.approx=PG.approx, PG.approx.m=PG.approx.m)
}

check_poisson_control <- function(control) {
  if (is.null(control)) return(poisson_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- poisson_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$nb.shape <- as.numeric(control[["nb.shape"]])
  if (length(control[["nb.shape"]]) != 1L || is.na(control[["nb.shape"]]) || control[["nb.shape"]] <= 0)
    stop("'nb.shape' must be a positive numerical scalar")
  if (!is_logical_scalar(control[["PG.approx"]])) stop("'PG.approx' must be TRUE or FALSE")
  control$PG.approx.m <- as.integer(control[["PG.approx.m"]])
  if (!length(control[["PG.approx.m"]])) stop("unexpected input for 'PG.approx.m'")
  if (any(is.na(control[["PG.approx.m"]]) | control[["PG.approx.m"]] < -2L | is.infinite(control[["PG.approx.m"]])))
    stop("'PG.approx.m' value(s) out of range or missing")
  control
}
