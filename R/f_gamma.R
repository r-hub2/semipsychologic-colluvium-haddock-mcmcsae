
#' @export
#' @rdname mcmcsae-family
f_gamma <- function(link="log", shape.vec = ~ 1, shape.prior = pr_gamma(0.1, 0.1),
                    control = set_MH(type="RWLN", scale=0.2, adaptive=TRUE)) {
  family <- "gamma"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  if (!inherits(shape.vec, "formula")) stop("'shape.vec' must be a formula")
  if (shape.prior$type == "fixed") {
    if (shape.prior$value <= 0) stop("gamma shape parameter must be positive")
  } else if (shape.prior$type == "exp") {  # special case of gamma
    shape.prior <- pr_gamma(shape=1, rate=1/shape.prior$scale)
  } else if (shape.prior$type != "gamma") {
    stop("unsupported prior for gamma shape parameter")
  }
  alpha.fixed <- shape.prior$type == "fixed"
  shape.prior$init(1L)  # scalar parameter
  alpha.scalar <- intercept_only(shape.vec)
  alpha0 <- NULL
  get_shape <- function() stop("please call method 'init' first")
  init <- function(data) get_shape <<- make_get_shape(data)
  make_get_shape <- function(data) {
    if (alpha.scalar) {
      alpha0 <<- 1
    } else {
      alpha0 <<- model_matrix(update.formula(shape.vec, ~ . - 1), data)
      if (ncol(alpha0) != 1L) stop("'shape.vec' must contain a single numeric variable")
      alpha0 <<- as.numeric(alpha0)
    }
    if (alpha.fixed) {
      alpha0 <<- alpha0 * shape.prior$value
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
    control$type <- match.arg(control$type, c("RWLN", "gamma"))
    make_draw_shape <- function(y) {
      # set up sampler for full conditional posterior for alpha, given linear predictor
      n <- length(y)
      # MH within Gibbs
      switch(control$type,
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
                #(alpha.star - alpha) * (sumlogy - sum(p[["e_"]]) - sum(y * exp(-p[["e_"]])))
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
                  #sum((alpha.vec - alpha.star.vec) * (y * exp(-p[["e_"]]) + p[["e_"]]))
              ))
          }
          add(f, quote(if (control$MH_accept(alpha.star, alpha, log.ar.post)) alpha.star else alpha))
          #add(f, quote(if (log(runif(1L)) < log.ar) alpha.star else alpha))
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
              #B0 <- B00 + sum(y * exp(-p[["e_"]])) + sum(p[["e_"]])
              B <- B0
              #A <- shape.prior$shape + n
              #B0 <- shape.prior$rate + sum(y * exp(-p[["e_"]])) - sumlogy + sum(p[["e_"]]) - n
              #B <- B0 + n
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
              #B0 <- B00 + sum(alpha0 * (y * exp(-p[["e_"]]) + p[["e_"]]))
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
  make_rpredictive <- function(newdata) {
    if (is.null(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      rpredictive <- function(p, lp) {
        alpha <- get_shape(p)
        rgamma(length(lp), shape=alpha, rate=alpha * exp(-lp))
      }
    } else {
      newn <- nrow(newdata)
      get_newshape <- make_get_shape(newdata)
      rpredictive <- function(p, lp) {
        alpha <- get_newshape(p)
        rgamma(newn, shape=alpha, rate=alpha * exp(-lp))
      }
    }
  }
  environment()
}

#' @export
#' @rdname mcmcsae-family
f_gaussian_gamma <- function(link="identity", var.data, ...) {
  family <- "gaussian_gamma"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  if (!inherits(var.data, "formula")) stop("'var.data' must be a formula")
  var.family <- f_gamma(...)
  self <- environment()
  sigmasq <- NULL
  init <- function(data) {
    var.family$init(data)
    copy_objects(var.family, self,
      c("alpha.fixed", "control", "get_shape",
        if (var.family$alpha.fixed) NULL else "shape.prior"))
    sigmasq <<- model_matrix(update.formula(var.data, ~ . - 1), data)
    if (ncol(sigmasq) != 1L) stop("'var.data' must contain a single numeric variable")
    sigmasq <<- as.numeric(sigmasq)
  }
  g <- function(y, p) y * p[["Q_"]] - log(p[["Q_"]])
  make_draw_shape <- function(y) {
    draw_shape <- var.family$make_draw_shape(y)
    assign("g", g, envir=environment(draw_shape))
    draw_shape
  }
  make_llh_gamma <- function(y) {
    llh_gamma <- var.family$make_llh(sigmasq)
    assign("g", g, envir=environment(llh_gamma))
    # NB gaussian contribution is added in samplers.R (until all llh computations have been moved to family objects)
    llh_gamma
  }
  make_llh_i <- function(y) {
    stop("TBI: pointwise log-likelihood for gaussian-gamma family")
    # TODO gaussian contribution
    var.family$make_llh_i(y)
  }
  self
}
