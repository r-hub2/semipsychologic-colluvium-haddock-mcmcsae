
#' Specify a Gamma sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Gamma sampling distribution.
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
  family <- "gamma"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
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
  alpha.fixed <- shape.prior[["type"]] == "fixed"
  shape.prior$init(1L)  # scalar parameter
  alpha.scalar <- intercept_only(shape.vec)
  alpha0 <- NULL
  get_shape <- function() stop("please call method 'init' first")
  init <- function(data, y=NULL) {
    get_shape <<- make_get_shape(data)
    if (!is.null(y)) {
      if (!is.numeric(y)) stop("non-numeric target value not allowed in case of Gamma sampling distribution")
      if (any(y <= 0)) stop("response variable modelled by Gamma distribution must be strictly positive")
    }
    y
  }
  make_get_shape <- function(data) {
    if (alpha.scalar) {
      alpha0 <<- 1
    } else {
      alpha0 <<- get_var_from_formula(shape.vec, data)
    }
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
  g <- function(y, p) y * exp(-p[["e_"]]) + p[["e_"]]
  if (!alpha.fixed) {
    if (!is.environment(control)) stop("f_gamma: 'control' argument must be an environment created with function set_MH")
    control$type <- match.arg(control[["type"]], c("RWLN", "gamma"))
    make_draw_shape <- function(y) {
      # set up sampler for full conditional posterior for alpha, given linear predictor
      n <- length(y)
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
  }
  make_llh <- function(y) {
    n <- length(y)
    if (alpha.fixed) {
      alpha <- get_shape()
      if (alpha.scalar) {
        llh_0 <- (alpha - 1) * sum(log(y))
        llh_0 <- llh_0 + n * (alpha * log(alpha) - lgamma(alpha))
        function(p) llh_0 - alpha * sum(g(y, p))
      } else {
        llh_0 <- sum((alpha - 1) * log(y))
        llh_0 <- llh_0 + sum(alpha * log(alpha) - lgamma(alpha))
        function(p) llh_0 - sum(alpha * g(y, p))
      }
    } else {
      llh_0 <- -sum(log(y))
      if (alpha.scalar) {
        function(p) {
          alpha <- get_shape(p)
          (1 - alpha) * llh_0 + n * (alpha * log(alpha) - lgamma(alpha)) - alpha * sum(g(y, p))
        }
      } else
        function(p) {
          alpha <- get_shape(p)
          llh_0 + sum(alpha * log(alpha * y) - lgamma(alpha)) - sum(alpha * g(y, p))
        }
    }
  }
  make_llh_i <- function(y) {
    n <- length(y)
    if (alpha.fixed) {
      alpha <- get_shape()
      pllh_0 <- alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * log(y)
      if (alpha.scalar)
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh_0[i], nr) - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      else
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh_0[i], nr) - rep_each(alpha[i], nr) * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
    } else {
      if (alpha.scalar)
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha <- as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE))
          alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * rep_each(log(y[i]), nr) +
            - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      else
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha_i <- outer(as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE)), alpha0[i])
          alpha_i * log(alpha_i) - lgamma(alpha_i) + (alpha_i - 1) * rep_each(log(y[i]), nr) +
            - alpha_i * (e_i + rep_each(y[i], nr) * exp(-e_i))
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
#'  of the coefficients is not exactly normal, but instead Multivariate Log inverse Gamma (MLiG),
#'  see also \code{\link{pr_MLiG}}.
#' @param ... further arguments passed to \code{\link{f_gamma}}.
#' @returns A family object.
f_gaussian_gamma <- function(link="identity", var.model, ...) {
  family <- "gaussian_gamma"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  y.family <- f_gaussian(link=link, var.prior=1, var.model=var.model)
  var.family <- f_gamma(...)
  self <- environment()
  sigmasq <- NULL  # placeholder for variance data vector
  init <- function(data, y=NULL) {
    if (is.null(y)) {
      y.family$init(data)
      var.family$init(data)
    } else {
      y <- y.family$init(data, y)
      sigmasq <<- get_response(var.model, data)
      if (is.null(sigmasq)) stop("no variance data vector specified")
      if (!all(is.finite(sigmasq))) stop(sum(!is.finite(sigmasq)), " missing or infinite value(s) in variance data vector")
      sigmasq <<- var.family$init(data, sigmasq)
    }
    copy_objects(y.family, self,
      c("sigma.fixed", "modeled.Q", "Q0", "Q0.type", "Vmod", "types")
    )
    copy_objects(var.family, self,
      c("alpha.fixed", "control", "get_shape",
        if (var.family[["alpha.fixed"]]) NULL else "shape.prior")
    )
    y
  }
  compute_Q <- NULL
  set_Vmod <- function(Vmod) {
    y.family$set_Vmod(Vmod)
    compute_Q <<- function(p, Qfactor=NULL) y.family$compute_Q(p, Qfactor)
  }
  g <- function(y, p) y * p[["Q_"]] - log(p[["Q_"]])
  make_draw_shape <- function(y) {
    draw_shape <- var.family$make_draw_shape(y)
    assign("g", g, envir=environment(draw_shape))
    draw_shape
  }
  make_llh <- function(y) {
    llh_gaussian <- y.family$make_llh(y)
    llh_gamma <- var.family$make_llh(sigmasq)
    assign("g", g, envir=environment(llh_gamma))
    function(p, SSR) llh_gaussian(p, SSR) + llh_gamma(p)
  }
  make_llh_i <- function(y) {
    llh_i_gaussian <- y.family$make_llh_i(y)
    llh_i_gamma <- var.family$make_llh_i(sigmasq)
    function(draws, i, e_i) llh_i_gaussian(draws, i, e_i) + llh_i_gamma(draws, i, e_i)
  }
  self
}
