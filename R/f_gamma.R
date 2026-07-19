#' Specify a Gamma sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Gamma sampling distribution.
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
#' alpha <- 1
#' mu <- exp(with(dat, 1 - 0.5*x + v[g]))
#' dat$y <- rgamma(n, shape=alpha, rate=alpha/mu)
#'
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") +       # fixed effects
#'       gen(factor = ~ g, name="v"),  # random intercepts
#'   data=dat, family="gamma"
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' waic(sim)
#' summary(sim)
#'
#' bayesplot::mcmc_recover_intervals(as.array(sim$gamma_shape_), alpha)
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), c(1, -0.5))
#' bayesplot::mcmc_recover_intervals(as.array(sim$v_sigma), 0.6)
#' yrep <- predict(sim, iters=sample(1:1000, 10))
#' bayesplot::pp_check(dat$y, as.matrix(yrep), bayesplot::ppc_dens_overlay) + ggplot2::xlim(0, 25)
#' }
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for the gamma distribution is \code{"log"}.
#' @param shape.vec optional formula specification of unequal shape parameter.
#' @param shape.prior prior for gamma shape parameter. Supported prior distributions:
#'  \code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_exp}} and
#'  \code{\link{pr_gamma}}. The current default is \code{pr_gamma(shape=0.1, rate=0.1)}.
#' @param control options for the Metropolis-Hastings algorithm employed
#'  in case the shape parameter is to be inferred. Function \code{\link{set_MH}}
#'  can be used to change the default options. The two choices of proposal
#'  distribution type supported are "RWLN" for a random walk proposal on the
#'  log-shape scale, and "gamma" for an approximating gamma proposal, found using
#'  an iterative algorithm. In the latter case, a Metropolis-Hastings accept-reject
#'  step is currently omitted, so the sampling algorithm is an approximate one,
#'  though often quite accurate and efficient.
#' @returns A family object.
#' @references
#'  J.W. Miller (2019).
#'    Fast and Accurate Approximation of the Full Conditional for Gamma Shape Parameters.
#'    Journal of Computational and Graphical Statistics 28(2), 476-480.
f_gamma <- function(link="log", shape.vec = ~ 1, shape.prior = pr_gamma(0.1, 0.1),
                    control = set_MH(type="RWLN", scale=0.2, adaptive=TRUE)) {
  link <- match.arg(link)
  if (!inherits(shape.vec, "formula")) stop("'shape.vec' must be a formula")
  if (is_numeric_scalar(shape.prior))
    shape.prior <- pr_fixed(value = shape.prior)
  else
    if (!is.environment(shape.prior)) stop("'shape.prior' must either be a numeric scalar or a prior specification")
  switch(shape.prior[["type"]],
    fixed = {
      if (shape.prior[["value"]] <= 0)
        stop("gamma shape parameter must be positive")
    },
    exp = {  # special case of gamma
      shape.prior <- pr_gamma(shape=1, rate=1/shape.prior[["scale"]])
    },
    gamma = {},
    stop("unsupported prior for gamma shape parameter")
  )
  shape.prior$init(1L)  # scalar parameter
  list(family="gamma", link=link, shape.vec=shape.vec,
       shape.prior=shape.prior, control=control, `_raw_`=TRUE)
}

ff_gamma <- function(link="log", shape.vec = ~ 1, shape.prior = pr_gamma(0.1, 0.1),
                     control = set_MH(type="RWLN", scale=0.2, adaptive=TRUE),
                     sm, data, y=NULL) {
  family <- "gamma"
  linkinv <- make.link(link)$linkinv
  f_mean <- function(eta) exp(eta)
  e.is.res <- FALSE
  sigma.fixed <- TRUE
  modeled.Q <- FALSE
  n <- n_row(data)
  Q0.type <- "unit"
  Q0 <- CdiagU(n)
  scale.e <- 1
  scale.sigma <- scale.e
  prior.only <- is.null(y)
  sc <- sm[["control"]]
  store_default <- function() if (alpha.fixed) NULL else "gamma_shape_"
  if (!prior.only) {
    if (!is.numeric(y)) stop("non-numeric target value not allowed in case of Gamma sampling distribution")
    if (any(y <= 0)) stop("response variable modelled by Gamma distribution must be strictly positive")
    draw <- if (sc[["compute.llh"]]) function(p) {p$llh_ <- llh(p)} else function(p) {}
  }
  alpha.fixed <- shape.prior[["type"]] == "fixed"
  alpha.scalar <- intercept_only(shape.vec)
  g <- function(y, p) y * exp(-p[["e_"]]) + p[["e_"]]  # TODO p[["Q_"]] case
  if (!alpha.fixed) {
    pred_pars <- function() "gamma_shape_"
    if (!is.environment(control)) stop("f_gamma: 'control' argument must be an environment created with function set_MH")
    control$type <- match.arg(control[["type"]], c("RWLN", "gamma"))
    rprior <- function(p) {
      p[["gamma_shape_"]] <- shape.prior$rprior()
      p
    }
    if (!prior.only) {
      make_draw_shape <- function() {
        # set up sampler for full conditional posterior for alpha, given linear predictor
        # MH within Gibbs
        switch(control[["type"]],
          RWLN = {
            f <- function(p) {
              alpha <- p[["gamma_shape_"]]
              alpha.star <- control$propose(alpha)
            }
            if (alpha.scalar) {
              sumlogy <- sum(log(y))
              f <- add(f, quote(
                log.ar.post <-
                  (shape.prior[["shape"]] - 1) * log(alpha.star/alpha) - shape.prior[["rate"]] * (alpha.star - alpha) +
                  n * (lgamma(alpha) - lgamma(alpha.star) + alpha.star * log(alpha.star) - alpha * log(alpha)) +
                  (alpha.star - alpha) * (sumlogy - sum(g(y, p)))
              ))
            } else {
              f <- f |>
                add(quote(alpha.vec <- alpha0 * alpha)) |>
                add(quote(alpha.star.vec <- alpha0 * alpha.star)) |>
                add(quote(
                  log.ar.post <-
                    (shape.prior[["shape"]] - 1) * log(alpha.star/alpha) - shape.prior[["rate"]] * (alpha.star - alpha) +
                    sum(lgamma(alpha.vec) - lgamma(alpha.star.vec)) +
                    sum(alpha.star.vec * log(alpha.star.vec * y)) -
                    sum(alpha.vec * log(alpha.vec * y)) +
                    sum((alpha.vec - alpha.star.vec) * g(y, p))
                ))
            }
            add(f, quote(if (control$MH_accept(alpha.star, alpha, log.ar.post)) alpha.star else alpha))
          },
          gamma = {
            # Miller's gamma approximation of the shape's full conditional
            #   iteration starts with approximate gamma density derived using Stirling's formula
            if (alpha.scalar) {
              A0 <- shape.prior[["shape"]] + 0.5 * n
              B00 <- shape.prior[["rate"]] - sum(log(y)) - n
              function(p) {
                A <- A0
                B0 <- B00 + sum(g(y, p))
                B <- B0
                for (i in 1:10) {
                  a <- A/B
                  A <- shape.prior[["shape"]] + n*a*(a * trigamma(a) - 1)
                  B <- B0 + (A - shape.prior[["shape"]])/a + n*(digamma(a) - log(a))
                  if (abs(a/(A/B) - 1) < 1e-8) break
                }
                rgamma(1L, A, B)
              }
              # To add an MH correction step:
              # (does not seem necessary, as the approximation is often excellent)
              #alpha.star <- rgamma(1L, A, B)
              #alpha <- p[["gamma_shape_"]]
              #log.ar <- n * (lgamma(alpha) - lgamma(alpha.star) + alpha.star * log(alpha.star) - alpha * log(alpha)) +
              #  (alpha.star - alpha) * (sumlogy - sum(p[["e_"]]) - sum(y * exp(-p[["e_"]]))) +
              #  (A - shape.prior$shape) * log(alpha/alpha.star) - (B - shape.prior$rate) * (alpha - alpha.star)
              #if (log(runif(1L)) < log.ar) alpha.star else alpha
            } else {
              ala0 <- sum(alpha0 * log(alpha0))
              A0 <- shape.prior[["shape"]] + 0.5 * n
              B00 <- shape.prior[["rate"]] - sum(alpha0 * (log(y) + 1))
              function(p) {
                A <- A0
                B0 <- B00 + sum(alpha0 * g(y, p))
                B <- B0
                for (i in 1:10) {
                  a <- A/B
                  a.vec <- a * alpha0
                  A <- shape.prior[["shape"]] - sum(a.vec) + sum(a.vec * a.vec * trigamma(a.vec))
                  B <- B0 + (A - shape.prior[["shape"]])/a - log(a) * sum(alpha0) +
                  sum(alpha0 * digamma(a.vec)) - ala0
                  if (abs(a/(A/B) - 1) < 1e-8) break
                }
                rgamma(1L, A, B)
              }
              # possibly add MH correction step
            }
          }
        )
      }
      draw_shape <- make_draw_shape()
      MHpars <- "gamma_shape_"
      if (control[["type"]] == "RWLN" && control[["adaptive"]])
        adapt <- function(ar) control$adapt(ar[["gamma_shape_"]])
      draw <- add(draw, quote(p$gamma_shape_ <- draw_shape(p)))
      # shape start value should not be too small
      start <- function(p) {
        p$gamma_shape_ <- check_and_get(p, "gamma_shape_", 1L,
          \() exp(runif(1L, -3, 4)), pos=TRUE)
        p
      }
    }
  }
  if (!prior.only) {
    draw <- add(draw, quote(p))
  }
  alpha0 <- NULL
  make_get_shape <- function(data) {
    if (alpha.scalar)
      alpha0 <<- 1
    else
      alpha0 <<- get_var_from_formula(shape.vec, data)
    rm(data)
    if (alpha.fixed) {
      alpha0 <<- alpha0 * shape.prior[["value"]]
      function(p) alpha0
    } else {
      if (alpha.scalar)
        function(p) p[["gamma_shape_"]]
      else
        function(p) alpha0 * p[["gamma_shape_"]]
    }
  }
  get_shape <- make_get_shape(data)

  if (!prior.only) {
    if (alpha.fixed) {
      alpha <- get_shape()
      pllh0 <- alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * log(y)
      if (alpha.scalar) {
        llh0 <- (alpha - 1) * sum(log(y))
        llh0 <- llh0 + n * (alpha * log(alpha) - lgamma(alpha))
        llh <- function(p) llh0 - alpha * sum(g(y, p))
        llh_i <- function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh0[i], nr) - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      } else {
        llh0 <- sum((alpha - 1) * log(y))
        llh0 <- llh0 + sum(alpha * log(alpha) - lgamma(alpha))
        llh <- function(p) llh0 - sum(alpha * g(y, p))
        llh_i <- function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh0[i], nr) - rep_each(alpha[i], nr) * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      }
    } else {
      llh0 <- -sum(log(y))
      if (alpha.scalar) {
        llh <- function(p) {
          alpha <- get_shape(p)
          (1 - alpha) * llh0 + n * (alpha * log(alpha) - lgamma(alpha)) - alpha * sum(g(y, p))
        }
        llh_i <- function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha <- as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE))
          alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * rep_each(log(y[i]), nr) +
            - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      } else {
        llh <- function(p) {
          alpha <- get_shape(p)
          llh0 + sum(alpha * log(alpha * y) - lgamma(alpha)) - sum(alpha * g(y, p))
        }
        llh_i <- function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha_i <- outer(as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE)), alpha0[i])
          alpha_i * log(alpha_i) - lgamma(alpha_i) + (alpha_i - 1) * rep_each(log(y[i]), nr) +
            - alpha_i * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      }
    }
  }
  make_rpredictive <- function(newdata, weights=NULL) {
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
        alpha <- get_shape(p)
        rgamma(nn, shape=alpha, rate=alpha*exp(-lp))
      }
    } else {
      function(p, lp) {
        alpha <- get_shape(p)
        rgamma(nn, shape=weights*alpha, rate=alpha*exp(-lp))
      }
    }
  }
  rm(data, sm)
  environment()
}

#' Specify a Gaussian-Gamma sampling distribution
#'
#' This function can be used in the \code{family} argument of
#' \code{\link{create_sampler}} or \code{\link{generate_data}} to specify a
#' Gaussian-Gamma sampling distribution, i.e., a Gaussian sampling distribution
#' whose variances are observed subject to error according to a Gamma
#' distribution.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for this distribution family is \code{"identity"}.
#' @param var.model a formula specifying the terms of a variance model.
#'  The left-hand side of the formula should specify the observed variances,
#'  unless the family object is used for data generation only.
#'  Several types of model terms on the right-hand side of the formula are supported:
#'  a regression term for the log-variance specified with \code{\link{vreg}(...)},
#'  and a term \code{\link{vfac}(...)} for multiplicative modelled factors
#'  at a certain level specified by a factor variable. In addition, \code{\link{reg}} and \code{\link{gen}}
#'  can be used to specify regression or random effect terms. In that case the prior distribution
#'  of the coefficients is not normal, but Multivariate Log inverse Gamma (MLiG),
#'  see also \code{\link{pr_MLiG}}.
#' @param gaussian.control a list with computational options passed to \code{\link{f_gaussian}}.
#' @param ... further arguments passed to \code{\link{f_gamma}}.
#' @returns A family object.
f_gaussian_gamma <- function(link="identity", var.model,
                             gaussian.control=gaussian_control(), ...) {
  lgau <- f_gaussian(link=link, var.prior=1, var.model=var.model,
             control=gaussian.control)
  lgau <- lgau[-which(names(lgau) %in% c("family", "_raw_"))]
  names(lgau)[names(lgau) == "control"] <- "gaussian.control"
  lgam <- f_gamma(...)
  lgam <- lgam[-which(names(lgam) %in% c("family", "link"))]
  c(list(family="gaussian_gamma"), lgau, lgam)
}

ff_gaussian_gamma <- function(link="identity", var.prior, var.vec,
                              prec.mat, var.model, logJacobian,
                              gaussian.control=gaussian_control(),
                              shape.vec, shape.prior, control,
                              sm, data, y=NULL) {
  family <- "gaussian_gamma"
  linkinv <- make.link(link)$linkinv
  prior.only <- is.null(y)
  y.family <- ff_gaussian(link=link, var.prior=var.prior, var.vec=var.vec,
    prec.mat=prec.mat, var.model=var.model, logJacobian=logJacobian,
    control=gaussian.control, sm=sm, data=data, y=y
  )
  if (prior.only) {
    sigmasq <- NULL
  } else {
    sigmasq <- get_response(var.model, data)
    if (is.null(sigmasq)) stop("no variance data vector specified")
  }
  var.family <- ff_gamma(shape.vec=shape.vec, shape.prior=shape.prior,
                         control=control, sm=sm, data=data, y=sigmasq)
  sc <- sm[["control"]]
  store_default <- function() {
    c(y.family$store_default(), var.family$store_default())
  }
  if (!prior.only) {
    pred_pars <- function() c(y.family$pred_vars(), var.family$pred_vars())
    var.family$g <- function(y, p) y * p[["Q_"]] - log(p[["Q_"]])
    environment(var.family[["g"]]) <- var.family
  }
  if (is.function(y.family$rprior) || is.function(var.family$rprior)) {
    rprior <- function(p) {
      if (is.function(y.family$rprior)) p <- y.family$rprior(p)
      if (is.function(var.family$rprior)) p <- var.family$rprior(p)
      p
    }
  }
  if (is.function(y.family[["draw"]])) {
    if (is.function(var.family[["draw"]])) {
      draw <- function(p) {
        p <- y.family$draw(p)
        if (sc[["compute.llh"]]) gaussian.llh <- p[["llh_"]]
        p <- var.family$draw(p)
        if (sc[["compute.llh"]]) p[["llh_"]] <- gaussian.llh + p[["llh_"]]
        p
      }
    } else {
      draw <- y.family[["draw"]]
    }
  } else if (is.function(var.family[["draw"]])) {
    draw <- var.family[["draw"]]
  }
  MHpars <- c(y.family[["MHpars"]], var.family[["MHpars"]])
  if (is.function(y.family[["adapt"]]) || is.function(var.family[["adapt"]])) {
    adapt <- function(ar) {
      if (is.function(y.family[["adapt"]])) y.family$adapt(ar)
      if (is.function(var.family[["adapt"]])) var.family$adapt(ar)
    }
  }
  if (is.function(y.family[["start"]]) || is.function(var.family[["start"]])) {
    start <- function(p) {
      if (is.function(y.family[["start"]])) p <- y.family$start(p)
      if (is.function(var.family[["start"]])) p <- var.family$start(p)
      p
    }
  }
  if (is.function(y.family[["drawMVNvarQ"]])) {
    cholQ <- y.family[["cholQ"]]
    drawMVNvarQ <- y.family[["drawMVNvarQ"]]
  }
  self <- environment()
  copy_refs(y.family, self,
    c("n", "sigma.fixed", "modeled.Q", "Q0", "Q0.type",
      "e.is.res", "Q_e", "Vmod", "scale.e", "scale.sigma",
      "SSR.name", "compute_SSR")
  )
  copy_refs(var.family, self,
    c("alpha.fixed", "control", "get_shape",
      if (var.family[["alpha.fixed"]]) NULL else "shape.prior")
  )
  compute_Q <- function(p, Qfactor=NULL) y.family$compute_Q(p, Qfactor)
  llh <- function(p) y.family$llh(p) + var.family$llh(p)
  llh_i <- function(draws, i, e_i) y.family$llh_i(draws, i, e_i) + var.family$llh_i(draws, i, e_i)
  rm(data, sm)
  self
}
