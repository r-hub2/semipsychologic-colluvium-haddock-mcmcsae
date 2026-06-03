#' Specify a negative binomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a negative binomial sampling distribution.
#'
#' The negative binomial distribution with shape r and probability p
#' has density
#' \deqn{p(y|r, p) = {\Gamma(y + r)\over y!\Gamma(r)}(1-p)^r p^y}
#' with mean \eqn{\mu = E(y|r,p) = {rp\over 1-p}} and variance
#' \eqn{V(y|r,p) = \mu(1 + \mu/r)}. The second term of the variance
#' can be interpreted as overdispersion with respect to a Poisson
#' distribution, which would correspond to the limit \eqn{r \rightarrow \infty}.
#' So the reciprocal shape \eqn{1/r} is an overdispersion parameter,
#' which typically is inferred. It is assigned a default prior, which
#' may be changed through argument \code{inv.shape.prior}.
#'
#' The only supported link function is \code{"log"}. Strictly speaking
#' the relation between mean \eqn{\mu} and linear predictor \eqn{\eta} is
#' \deqn{\log\mu = \log r + \log{p\over 1-p} = \log r + \eta}
#' This way the likelihood function has the same form as that of logistic
#' binomial regression, so that a Polya-Gamma data augmentation sampling
#' algorithm can be employed. Note that the fact that the linear predictor
#' \eqn{\eta} does not include \eqn{\log r} effectively changes the
#' interpretation of its intercept.
# NB negative binomial regression with unknown shape does NOT belong to the
#    class of generalised linear models
#'
#' @examples
#' \dontrun{
#' n <- 1000
#' nT <- 40
#' dat <- data.frame(
#'   x = rnorm(n),
#'   t = factor(sample(1:nT, n, replace=TRUE), levels=1:nT)
#' )
#' model <- ~ reg(~ x, prior=pr_normal(precision=1), name="beta") +
#'            gen(factor = ~ RW1(t), name="v")
#' gd <- generate_data(model, dat, family=f_negbinomial())
#' str(gd)
#' 
#' dat$y <- gd$y
#' sampler <- create_sampler(
#'   model <- y ~ reg(~ x, name="beta") +
#'              gen(factor = ~ RW1(t), name="v"),
#'   data=dat, family=f_negbinomial()
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' summ <- summary(sim)
#' loo(sim)
#' plot.ts(gd$pars$v)
#' lines(summ$v[, "Mean"], col=2)
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), gd$pars$beta)
#' bayesplot::mcmc_recover_hist(as.array(sim$negbin_shape_), gd$pars$negbin_shape_)
#' }
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for the negative binomial sampling distribution is \code{"log"}.
#' @param shape.vec optional formula specification of unequal shape values. The
#'  negative binomial (vector) shape parameter is then equal to this vector of
#'  shape values, multiplied by the scalar shape parameter, whose prior is
#'  specified through \code{inv.shape.prior}.
#' @param inv.shape.prior Prior on the (scalar) \emph{reciprocal} shape parameter,
#'  i.e. the overdispersion parameter. Supported prior distributions are
#'  \code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_invchisq}} and
#'  \code{\link{pr_gig}}. The current default is \code{pr_invchisq(df=1, scale=1)}.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{negbinomial_control}}.
#' @returns A family object.
#' @references
#'  N. Polson, J.G. Scott and J. Windle (2013).
#'    Bayesian Inference for Logistic Models Using Polya-Gamma Latent Variables.
#'    Journal of the American Statistical Association 108(504), 1339-1349.
#'
#'  M. Zhou and L. Carin (2015).
#'    Negative Binomial Process Count and Mixture Modeling.
#'    IEEE Transactions on Pattern Analysis and Machine Intelligence 37(2), 307-320.
f_negbinomial <- function(link="log", shape.vec = ~ 1, inv.shape.prior = pr_invchisq(df=1, scale=1),
                          control = negbinomial_control()) {
  link <- match.arg(link)
  control <- check_negbinomial_control(control)
  if (!inherits(shape.vec, "formula"))
    stop("'shape.vec' must be a formula")
  if (is_numeric_scalar(inv.shape.prior))
    inv.shape.prior <- pr_fixed(value = inv.shape.prior)
  else
    if (!is.environment(inv.shape.prior)) stop("'inv.shape.prior' must either be a numeric scalar or a prior specification")
  switch(inv.shape.prior[["type"]],
    fixed = {
      if (inv.shape.prior[["value"]] <= 0)
        stop("negative binomial inverse shape (dispersion) parameter must be positive")
    },
    invchisq = {
      if (is.list(inv.shape.prior[["df"]]))
        stop("modelled degrees of freedom parameter not supported for negative binomial (inverse) shape parameter")
    },
    gig = {},
    stop("unsupported prior")
  )
  inv.shape.prior$init(n=1L)
  list(family="negbinomial", link=link, shape.vec=shape.vec,
       inv.shape.prior=inv.shape.prior, control=control, `_raw_`=TRUE)
}

ff_negbinomial <- function(link="log", shape.vec = ~ 1, inv.shape.prior = pr_invchisq(df=1, scale=1),
                           control = negbinomial_control(),
                           sm, data, y=NULL, famid=NULL, sub=NULL) {
  family <- "negbinomial"
  linkinv <- make.link(link)$linkinv
  e.is.res <- FALSE
  sigma.fixed <- TRUE
  modeled.Q <- TRUE
  prior.only <- is.null(y)
  n <- n_row(data)
  Q0.type <- "unit"
  Q0 <- CdiagU(n)
  shape.fixed <- inv.shape.prior[["type"]] == "fixed"
  shape.scalar <- intercept_only(shape.vec)
  if (shape.fixed)
    f_mean <- function(eta) shape0 * exp(eta)
  else
    f_mean <- function(eta) stop("TODO fitted values or residuals for negative binomial family with inferred shape parameter")
  scale.e <- 2.5
  scale.sigma <- scale.e
  multifam <- !is.null(famid)
  if (shape.fixed) {
    nb.shape.name <- NULL
  } else {
    if (multifam)
      nb.shape.name <- paste0(famid, "_negbin_shape_")
    else
      nb.shape.name <- "negbin_shape_"
  }
  sc <- sm[["control"]]
  store_default <- function() nb.shape.name
  shape0 <- NULL
  make_get_shape <- function(data) {
    if (shape.scalar) {
      shape0 <<- 1
    } else {
      shape0 <<- get_var_from_formula(shape.vec, data)
      if (all(length(shape0) != c(1L, n_row(data))))
        stop("wrong length for shape vector")
      if (anyNA(shape0)) stop("missings in shape vector not allowed")
      if (any(shape0 <= 0)) stop("shape vector must be positive")
    }
    if (shape.fixed) {
      shape0 <<- shape0 / inv.shape.prior[["value"]]
      function(p) shape0
    } else {
      if (shape.scalar)
        function(p) p[[nb.shape.name]]
      else
        function(p) shape0 * p[[nb.shape.name]]
    }
  }
  get_shape <- make_get_shape(data)
  if (!shape.fixed) {
    rprior <- function(p) {
      p[[nb.shape.name]] <- 1/inv.shape.prior$rprior()
      p
    }
  }
  if (!prior.only) {
    if (multifam) {
      draw <- function(p) {
        p$llh_ <- p[["llh_"]] + llh(p)
        e_ <- p[["e_"]][sub]
      }
    } else {
      draw <- function(p) {
        p$llh_ <- llh(p)
        e_ <- p[["e_"]]
      }
    }
    start <- function(p) {}
    if (shape.fixed) {
      y_shifted <- 0.5 * (y - shape0)
      if (sc[["single.block"]]) {
        Q_e <- function(p) y_shifted
      } else {
        if (multifam)
          Q_e <- function(p) y_shifted - p[["Q_"]][sub] * p[["e_"]][sub]
        else
          Q_e <- function(p) y_shifted - p[["Q_"]] * p[["e_"]]
      }
    } else {  # inferred overdispersion
      if (sc[["single.block"]]) {
        Q_e <- function(p) 0.5 * (y - get_shape(p))
      } else {
        if (multifam)
          Q_e <- function(p) 0.5 * (y - get_shape(p)) - p[["Q_"]][sub] * p[["e_"]][sub]
        else
          Q_e <- function(p) 0.5 * (y - get_shape(p)) - p[["Q_"]] * p[["e_"]]
      }
      make_draw_shape <- function(y) {
        # draw latent L_i (i=1:n) from its CRT f.c.
        mCRT <- control[["CRT.approx.m"]]
        draw <- function(p, e_) {L <- CrCRT(y, get_shape(p), mCRT)}
        switch(inv.shape.prior[["type"]],
          # TODO check shape0 appearance in modeled scale invchisq and gig cases
          invchisq = {
            inv.shape.prior$make_draw()
            if (is.list(inv.shape.prior[["scale"]])) {
              if (shape.scalar)
                draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / inv.shape.prior$draw(2*sum(L), 2*shape0*sum(log1pexpC(e_)), p[[.(nb.shape.name)]])))
              else
                draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / inv.shape.prior$draw(2*sum(L), 2*sum(shape0*log1pexpC(e_)), p[[.(nb.shape.name)]])))
            } else {
              if (shape.scalar)
                draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / inv.shape.prior$draw(2*sum(L), 2*shape0*sum(log1pexpC(e_)))))
              else
                draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / inv.shape.prior$draw(2*sum(L), 2*sum(shape0*log1pexpC(e_)))))
            }
          },
          gig = {
            prior.p <- inv.shape.prior[["p"]]
            prior.a <- inv.shape.prior[["a"]]
            prior.b <- inv.shape.prior[["b"]]
            rgig <- GIGrvg::rgig
            if (shape.scalar)
              draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / rgig(1L, prior.p - sum(L), prior.b + 2*shape0*sum(log1pexpC(e_)), prior.a)))
            else
              draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- 1 / rgig(1L, prior.p - sum(L), prior.b + 2*sum(shape0*log1pexpC(e_)), prior.a)))
          }
        )
        draw
      }
      draw_shape <- make_draw_shape(y)
      if (inv.shape.prior[["type"]] == "fixed") {
        start <- add(start, bquote(
          p[[.(nb.shape.name)]] <- check_and_get(p, .(nb.shape.name), 1L,
            \() 1/inv.shape.prior$rprior(), pos=TRUE)
        ))
      } else {
        start <- add(start, bquote(
          p[[.(nb.shape.name)]] <- check_and_get(p, .(nb.shape.name), 1L,
            \() runif(1L, 0.1, 10), pos=TRUE)
        ))
      }
      draw <- add(draw, bquote(p[[.(nb.shape.name)]] <- draw_shape(p, e_)))
    }
    rPolyaGamma <- get_PG_sampler(n, control[["PG.approx"]], control[["PG.approx.m"]])
    if (shape.fixed) {
      ny <- y + shape0
    } else {
      y <- as.numeric(y)  # for C++ rCRT function; alternatively force to be integer
      draw <- add(draw, quote(ny <- y + get_shape(p)))
    }
    start <- add(start, quote(ny <- y + get_shape(p)))
    if (multifam) {
      draw <- add(draw, quote(p$Q_[sub] <- rPolyaGamma(ny, e_)))
    } else {
      draw <- add(draw, quote(p$Q_ <- rPolyaGamma(ny, e_)))
      start <- add(start, quote(p$Q_ <- check_and_get(p, "Q_", n, \() rPolyaGamma(ny, p[["e_"]]), pos=TRUE)))
    }
    draw <- add(draw, quote(p))
    start <- add(start, quote(p))
    if (!is.null(sc[["CG"]]) || sc[["cMVN.sampler"]]) {
      # set up a function that multiplies by L Chol factor of Q, for sampling from N(., Q)
      cholQ <- build_chol(runif(n, 0.9, 1.1))
      # draw from MVN with variance(!) Q
      drawMVNvarQ <- function(p) {
        cholQ$update(if (multifam) p[["Q_"]][sub] else p[["Q_"]])
        cholQ$Ltimes(Crnorm(n), transpose=FALSE)
      }
    }
    if (!is.numeric(y)) stop("non-numeric target value not allowed in case of negative binomial sampling distribution")
    # NB the algorithm still runs with negative responses, but it probably makes not much sense
    if (any(y < 0)) warn("negative response value(s)")
    if (any(abs(round(y) - y) > .tol)) warn("non-integral values modelled by negative binomial sampling distribution")
  }

  if (shape.fixed) {
    llh_0 <- sum(negbinomial_coef(shape0, y))
    llh <- function(p) {
      ny <- y + shape0  # maybe precompute?
      e_ <- if (multifam) p[["e_"]][sub] else p[["e_"]]
      llh_0 + sum(y * e_ - ny * log1pexpC(e_))
    }
    llh_i <- function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      if (shape.scalar)
        rep_each(negbinomial_coef(shape0, y[i]), nr) - shape0 * e_i - rep_each(y[i] + shape0, nr) * log1pexpC(-e_i)
      else
        rep_each(negbinomial_coef(shape0[i], y[i]), nr) - rep_each(shape0[i], nr) * e_i - rep_each(y[i] + shape0[i], nr) * log1pexpC(-e_i)
    }
  } else {
    llh_0 <- -sum(lgamma(y + 1))
    llh <- function(p) {
      r <- get_shape(p)
      ny <- y + r
      e_ <- if (multifam) p[["e_"]][sub] else p[["e_"]]
      llh_0 + sum(lgamma(ny) - lgamma(r)) + sum(y * e_ - ny * log1pexpC(e_))
    }
    llh_i <- function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      r <- as.numeric(as.matrix.dc(draws[[nb.shape.name]], colnames=FALSE))
      if (shape.scalar)
        r <- shape0 * r
      else
        r <- r * rep_each(shape0[i], nr)
      yi <- rep_each(y[i], nr)
      nyi <- yi + r
      negbinomial_coef(r, yi) + yi * e_i - nyi * log1pexpC(e_i)
    }
  }
  make_rpredictive <- function(newdata, weights=NULL) {
    # NB definition of rnbinom has p <-> 1-p
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
    } else {
      nn <- nrow(newdata)
      get_shape <- make_get_shape(newdata)
    }
    if (is.null(weights)) {
      function(p, lp) {
        size <- get_shape(p)
        rnbinom(nn, size, mu=size*exp(lp))
      }
      #rnbinom(length(lp), size=get_shape(p), prob=1/(1 + exp(lp)))
    } else {
      function(p, lp) {
        size <- weights * get_shape(p)
        rnbinom(nn, size, mu=size*exp(lp))
      }
    }
  }
  rm(data, sm)
  if (!multifam) rm(famid, sub)
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
#' @param CRT.approx.m scalar integer specifying the degree of approximation to sampling
#'  from a Chinese Restaurant Table distribution. The approximation is based on Le Cam's theorem.
#'  Larger values yield a slower but more accurate sampler.
#' @returns A list with computational options for the sampling algorithm.
negbinomial_control <- function(PG.approx=TRUE, PG.approx.m=-2L, CRT.approx.m=20L) {
  list(PG.approx=PG.approx, PG.approx.m=PG.approx.m, CRT.approx.m=CRT.approx.m)
}

check_negbinomial_control <- function(control) {
  if (is.null(control)) return(negbinomial_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- negbinomial_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  if (!is_logical_scalar(control[["PG.approx"]])) stop("'PG.approx' must be TRUE or FALSE")
  control$PG.approx.m <- as.integer(control[["PG.approx.m"]])
  if (!length(control[["PG.approx.m"]])) stop("unexpected input for 'PG.approx.m'")
  if (any(is.na(control[["PG.approx.m"]]) | control[["PG.approx.m"]] < -2L | is.infinite(control[["PG.approx.m"]])))
    stop("'PG.approx.m' value(s) out of range or missing")
  control$CRT.approx.m <- as.integer(control[["CRT.approx.m"]])
  if (length(control[["CRT.approx.m"]]) != 1L || control[["CRT.approx.m"]] < 1L)
    stop("'CRT.approx.m' must be a positive scalar integer")
  control
}
