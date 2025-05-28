
#' Specify a Gaussian sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Gaussian sampling distribution.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link functions
#'  for the Gaussian distribution is \code{"identity"}.
#' @param var.prior prior for the variance parameter of a Gaussian sampling distribution.
#'  This can be specified by a call to one of the prior specification functions
#'  \code{\link{pr_invchisq}}, \code{\link{pr_exp}}, \code{\link{pr_gig}} or \code{\link{pr_fixed}} for
#'  inverse chi-squared, exponential, generalised inverse gaussian or degenerate prior distribution,
#'  respectively. The default is an improper prior \code{pr_invchisq(df=0, scale=1)}. A half-t prior on the
#'  standard deviation can be specified using \code{\link{pr_invchisq}} with a chi-squared distributed scale
#'  parameter.
#' @param var.vec a formula to specify unequal variances, i.e. heteroscedasticity.
#'  The default corresponds to equal variances.
#' @param prec.mat a possibly non-diagonal positive-definite symmetric matrix
#'  interpreted as the precision matrix, i.e. inverse of the covariance matrix.
#'  If this argument is specified \code{var.vec} is ignored.
#' @param var.model a formula specifying the terms of a variance model in the case of a Gaussian likelihood.
#'  Several types of terms are supported: a regression term for the log-variance
#'  specified with \code{\link{vreg}(...)}, and a term \code{\link{vfac}(...)} for multiplicative modelled factors
#'  at a certain level specified by a factor variable. By using unit-level inverse-chi-squared factors the marginal
#'  sampling distribution becomes a Student-t distribution, and by using unit-level exponential factors it becomes
#'  a Laplace or double exponential distribution. In addition, \code{\link{reg}} and \code{\link{gen}}
#'  can be used to specify regression or random effect terms. In that case the prior distribution
#'  of the coefficients is not exactly normal, but instead Multivariate Log inverse Gamma (MLiG),
#'  see also \code{\link{pr_MLiG}}.
#' @param logJacobian if the data are transformed the logarithm of the Jacobian can be supplied so that it
#'  is incorporated in all log-likelihood computations. This can be useful for comparing information criteria
#'  for different transformations. It should be supplied as a vector of the same size as the response variable,
#'  and is currently only supported if \code{family="gaussian"}.
#'  For example, when a log-transformation is used on response vector \code{y}, the vector \code{-log(y)}
#'  should be supplied.
#' @returns A family object.
f_gaussian <- function(link="identity", var.prior = pr_invchisq(df=0, scale=1),
                       var.vec = ~ 1, prec.mat=NULL, var.model=NULL, logJacobian=NULL) {
  family <- "gaussian"
  link <- match.arg(link)
  linkinv <- identity
  if (is_numeric_scalar(var.prior))
    var.prior <- pr_fixed(value = var.prior)
  else
    if (!is.environment(var.prior)) stop("'var.prior' must either be a numeric scalar or a prior specification")
  switch(var.prior[["type"]],
    fixed = {
      if (var.prior[["value"]] <= 0)
        stop("gaussian variance parameter must be positive")
      },
    invchisq = {
      if (is.list(var.prior[["df"]])) stop("modelled degrees of freedom parameter not supported in 'var.prior'")
    },
    exp=, gig = {},
    stop("unsupported prior")
  )
  var.prior$init(n=1L)
  sigma.fixed <- var.prior[["type"]] == "fixed" && var.prior[["value"]] == 1
  if (!is.null(prec.mat)) {
    if (!is_a_matrix(prec.mat)) stop("'prec.mat' must be a matrix object")
    var.vec <- ~ 1  # ignored in this case
  }
  Q0 <- NULL  # precision matrix
  Q0.type <- NULL
  modeled.Q <- !is.null(var.model)
  if (modeled.Q) {
    if (!inherits(var.model, "formula")) stop("'var.model' must be a formula")
    types <- NULL
    compute_Qfactor <- NULL
    compute_Q <- NULL
  }
  Vmod <- NULL
  n <- NULL
  if (!sigma.fixed || modeled.Q) {
    if (sigma.fixed) {
      rprior <- function(p) {}
    } else {
      rprior <- function(p) {
        p$sigma_ <- sqrt(var.prior$rprior())
      }
    }
  }
  init <- function(data, y=NULL) {
    # y=NULL in case of prior sampling/prediction
    if (!is.null(y)) {
      if (is.logical(y)) y <- as.integer(y)
      if (!is.numeric(y)) stop("non-numeric response variable")
      if (var.prior[["type"]] == "invchisq") var.prior$make_draw()
    }
    n <<- n_row(data)
    if (!is.null(prec.mat)) {
      Q0 <<- economizeMatrix(prec.mat, symmetric=TRUE, check=TRUE)
      if (!identical(dim(Q0), c(n, n))) stop("incompatible precision matrix")
      prec.mat <<- "prec.mat.used"
    } else if (is.vector(var.vec)) {
      # allow vector input for backward compatibility
      if (!is.numeric(var.vec) || anyNA(var.vec) || any(var.vec <= 0)) stop("invalid input for 'var.vec'")
      Q0 <<- Cdiag(if (length(var.vec) == 1L) rep.int(1/var.vec, n) else 1/var.vec)
    } else if (intercept_only(var.vec)) {
      Q0 <<- CdiagU(n)
    } else {
      temp <- get_var_from_formula(var.vec, data)
      if (any(temp <= 0)) stop("non-positive variance(s) in 'var.vec'")
      if (length(temp) == 1L) {
        Q0 <<- Cdiag(rep.int(1/temp, n))
      } else {
        if (length(temp) != n) stop("variance vector has wrong length")
        Q0 <<- Cdiag(1/temp)
      }
    }
    if (isDiagonal(Q0)) {
      if (is_unit_diag(Q0))
        Q0.type <<- "unit"
      else
        Q0.type <<- "diag"
    } else {
      Q0.type <<- "symm"  # non-diagonal precision matrix
    }
    if (modeled.Q) {
      var.model <<- standardize_formula(var.model, "reg", data=data)  # pass data to interpret '.'
      Vmod <<- to_mclist(var.model, prefix="v")
      if (!length(Vmod)) stop("empty 'var.model'")
      types <<- get_types(Vmod)
      if (any(types %in% c("mec", "brt"))) stop("'mec' and 'brt' can only be used in mean model specification")
    }
    if (!sigma.fixed && !modeled.Q) rprior <<- add(rprior, quote(p))
    y
  }
  set_Vmod <- function(Vmod) {
    Vmod <<- Vmod
    for (k in seq_along(Vmod)) {
      mc <- Vmod[[k]]
      switch(types[k],
        vreg=, reg = {
          rprior <<- add(rprior, bquote(p[[.(mc[["name"]])]] <- Vmod[[.(k)]]$rprior(p)))
        },
        vfac=, gen = {
          rprior <<- add(rprior, bquote(p <- Vmod[[.(k)]]$rprior(p)))
        }
      )
    }
    rprior <<- add(rprior, quote(p))
    # compute product of precision factors
    compute_Qfactor <<- function(p) {
      out <- Vmod[[1L]]$compute_Qfactor(p)
      for (mc in Vmod[-1L]) out <- out * mc$compute_Qfactor(p)
      out
    }
    # compute data-level precision matrix from Q0 and scale factor computed by compute_Qfactor
    switch(Q0.type,
      unit =
        compute_Q <<- function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_ <- compute_Qfactor(p)
          p
        },
      diag =
        compute_Q <<- function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_ <- Q0@x * compute_Qfactor(p)
          p
        },
      symm =  # store full precision matrix
        compute_Q <<- function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) Qfactor <- compute_Qfactor(p)
          p$QM_ <- block_scale_dsCMatrix(Q0, Qfactor)
          p$Q_ <- Qfactor  # by default in store.mean for use in compute_DIC
          p
        }
    )
  }
  make_llh <- function(y) {
    llh_0 <- -0.5 * n * log(2*pi)  # constant term of log-likelihood
    if (!is.null(logJacobian)) {
      logJacobian <- as.numeric(logJacobian)
      if (length(logJacobian) != n) stop("'logJacobian' should be a vector of the same size as the response vector")
      if (anyNA(logJacobian)) stop("missing values in 'logJacobian'")
      llh_0 <- llh_0 + sum(logJacobian)
    }
    if (!modeled.Q || Q0.type == "symm")
      llh_0 <- llh_0 + 0.5 * as.numeric(determinant(1 * Q0, logarithm=TRUE)$modulus)  # 1 * Q0 to ensure no chol object is stored with Q0 and possibly other refs to Q0
    # log-likelihood function; SSR need not be recomputed if it is provided
    function(p, SSR) {
      out <- llh_0
      if (modeled.Q) out <- out + 0.5 * sum(log(p[["Q_"]]))
      if (sigma.fixed)
        out - 0.5 * SSR
      else
        out - n * log(p[["sigma_"]]) - 0.5 * SSR / p[["sigma_"]]^2
    }
  }
  make_llh_i <- function(y) {
    # for WAIC computation: compute log-likelihood for each observation/batch of observations, vectorized over parameter draws and observations
    llh_0_i <- -0.5 * log(2*pi)
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      all.units <- length(i) == n
      if (sigma.fixed) {
        q <- 1
      } else {
        sigma <- as.numeric(as.matrix.dc(draws[["sigma_"]], colnames=FALSE))
        q <- 1/sigma^2
      }
      if (modeled.Q) {
        q <- matrix(q, nr, length(i))
        for (mc in Vmod) {
          if (all.units) Xi <- mc[["X"]] else Xi <- mc[["X"]][i, , drop=FALSE]
          if (mc[["type"]] == "vfac")
            q <- q * tcrossprod(1 / as.matrix.dc(draws[[mc[["name"]]]], colnames=FALSE), Xi)
          else
            q <- q * exp(-tcrossprod(as.matrix.dc(draws[[mc[["name"]]]], colnames=FALSE), Xi))
        }
      }
      switch(Q0.type,
        diag = q <- q * rep_each(Q0@x[i], each=nr),
        symm = {
          dQ0 <- rep_each(diag(Q0)[i], each=nr)
          q <- q * dQ0
          if (all.units) {
            # pointwise loo-llh for non-factorizable models: p(y_i|y_{-i},theta)
            # if length(i) < n we currently ignore this correction; should warn in the calling function
            e_i <- (1 / dQ0) * (e_i %m*m% Q0)
          }
        }
      )
      if (is.null(logJacobian))
        llh_0_i + 0.5 * ( log(q) - q * e_i^2 )
      else
        llh_0_i + matrix(rep_each(logJacobian[i], nr), nr, length(i)) + 0.5 * ( log(q) - q * e_i^2 )
    }
  }
  make_get_sds <- function(newdata) {
    if (identical(prec.mat, "prec.mat.used")) stop("out-of-sample prediction not supported if 'prec.mat' is used")
    if (intercept_only(var.vec)) {
      sds0 <- 1
    } else {
      sds0 <- get_var_from_formula(var.vec, newdata)
      if (any(sds0 <= 0)) stop("non-positive variance(s) in 'var.vec' for prediction")
      if (all(length(sds0) != c(1L, nrow(newdata)))) stop("wrong length for prediction variance vector")
      sds0 <- invsqrt(sds0)
    }
    if (modeled.Q) {
      V <- list()
      for (Vmc in Vmod) {
        if (any(Vmc[["type"]] == c("vfac", "vreg"))) {
          V[[Vmc[["name"]]]] <- Vmc$make_predict_Vfactor(newdata)
        } else {
          V[[Vmc[["name"]]]] <- Vmc$make_predict(newdata)
        }
      }
      function(p) {
        var <- sds0^2
        for (Vmc in Vmod) {
          if (any(Vmc[["type"]] == c("vfac", "vreg")))
            var <- var * V[[Vmc[["name"]]]](p)
          else
            var <- var * exp(V[[Vmc[["name"]]]](p))
        }
        sqrt(var)
      }
    } else {
      function(p) sds0
    }
  }
  # weights: passed from predict.mcdraws
  #   can be either a numeric scalar, or a vector of length n or nrow(newdata) if the latter is provided
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
      # need cholQ
      if (modeled.Q) {
        # can only be in-sample prediction
        cholQ <- build_chol(
          switch(Q0.type,
            unit = rep.int(1, n),
            diag = Q0@x,
            symm = Q0
          )
        )
      } else {
        # either in-sample or out-of-sample prediction
        if (nn != n) {
          # in this case Q0.type == "unit" (checked in predict)
          cholQ <- build_chol(CdiagU(nn))
        } else {
          # in-sample prediction
          cholQ <- build_chol(Q0)
        }
      }
      function(p, lp) {
        sigma <- if (sigma.fixed) 1 else p[["sigma_"]]
        if (modeled.Q) {
          if (is.null(p[["Q_"]])) p <- compute_Q(p)
          if (Q0.type == "symm")
            cholQ$update(p[["QM_"]])
          else
            cholQ$update(p[["Q_"]])
        }
        if (is.null(weights))
          lp + drawMVN_cholQ(cholQ, sd=sigma)
        else
          weights * lp + sqrt(weights) * drawMVN_cholQ(cholQ, sd=sigma)
      }
    } else {
      nn <- nrow(newdata)
      get_sds <- make_get_sds(newdata)
      # here we assume that Q0 is diagonal
      if (is.null(weights))
        function(p, lp) {
          sigma <- if (sigma.fixed) 1 else p[["sigma_"]]
          lp + sigma * get_sds(p) * Crnorm(nn)
        }
      else
        function(p, lp) {
          sigma <- if (sigma.fixed) 1 else p[["sigma_"]]
          weights * lp + sigma * sqrt(weights) * get_sds(p) * Crnorm(nn)
        }
    }
  }
  environment()
}
