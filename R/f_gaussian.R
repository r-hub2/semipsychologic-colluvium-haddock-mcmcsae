#' Specify a Gaussian sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Gaussian sampling distribution.
#'
#' @examples
#' \dontrun{
#' n <- 4000
#' m <- 25
#' dat <- data.frame(
#'   x = rnorm(n),
#'   g = factor(sample(1:m, n, replace=TRUE), levels=1:m)
#' )
#' v <- rnorm(m, sd=0.6)
#' dat$y <- rnorm(n, mean = with(dat, 1 - 0.5*x + v[g]), sd=0.4)
#'
#' sampler <- create_sampler(
#'   y ~ x + (1|g), data=dat
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' summary(sim)
#'
#' # more explicit specification, allowing non-default names, priors etc.
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") + gen(~1, factor = ~ g, name="v"),
#'   data=dat,
#'   family = f_gaussian(var.prior = pr_fixed(value = 0.4^2))
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' summary(sim)
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), c(1, -0.5))
#' bayesplot::mcmc_recover_hist(as.array(sim$v_sigma), 0.6)
#' bayesplot::mcmc_recover_scatter(as.array(sim$v), v)
#'
#' # heteroscedastic data
#' dat$y <- rnorm(n,
#'   mean = with(dat, 1 - 0.5*x + v[g]),
#'   sd = with(dat, exp(0.5*(0.2 + 0.5*x)))
#' )
#'
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") + gen(~1, factor = ~ g, name="v"),
#'   data=dat,
#'   family = f_gaussian(
#'     var.prior = pr_fixed(value = 1),
#'     var.model = ~ 1 + x
#'   )
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' compute_WAIC(sim)
#' summary(sim)
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), c(1, -0.5))
#' bayesplot::mcmc_recover_intervals(as.array(sim$vreg1), c(0.2, 0.5))
#' bayesplot::mcmc_recover_hist(as.array(sim$v_sigma), 0.6)
#' bayesplot::mcmc_recover_scatter(as.array(sim$v), v)
#' }
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for a Gaussian distribution is \code{"identity"}.
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
#'  for different transformations. It should be supplied as a vector of the same size as the response variable.
#'  For example, when a log-transformation is used on response vector \code{y}, the vector \code{-log(y)}
#'  should be supplied.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{gaussian_control}}.
#' @returns A family object.
f_gaussian <- function(link="identity", var.prior = pr_invchisq(df=0, scale=1),
                       var.vec = ~ 1, prec.mat=NULL, var.model=NULL, logJacobian=NULL,
                       control = gaussian_control()) {
  link <- match.arg(link)
  control <- check_gaussian_control(control)
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
  if (!is.null(prec.mat)) {
    if (!is_a_matrix(prec.mat)) stop("'prec.mat' must be a matrix object")
    var.vec <- ~ 1  # ignored in this case
  }
  if (!is.null(var.model) && !inherits(var.model, "formula"))
    stop("'var.model' must be a formula")
  # or, make it an environment, and check it, cf. set_MH
  list(family="gaussian", link=link, var.prior=var.prior, var.vec=var.vec,
       prec.mat=prec.mat, var.model=var.model,
       logJacobian=logJacobian, control=control, `_raw_`=TRUE)
}

# additional arguments sc (sampler.control), data, y
ff_gaussian <- function(link="identity", var.prior = pr_invchisq(df=0, scale=1),
                        var.vec = ~ 1, prec.mat=NULL, var.model=NULL, logJacobian=NULL,
                        control = gaussian_control(),
                        sm, data, y=NULL, famid=NULL, sub=NULL) {
  family <- "gaussian"
  linkinv <- identity
  f_mean <- identity  # mean function acting on linear predictor
  e.is.res <- TRUE
  sigma.fixed <- var.prior[["type"]] == "fixed" && var.prior[["value"]] == 1
  modeled.Q <- !is.null(var.model)
  n <- n_row(data)
  prior.only <- is.null(y)
  multifam <- !is.null(famid)
  if (multifam) {
    sd.name <- if (sigma.fixed) NULL else paste0(famid, "_sigma_")
    SSR.name <- paste0(famid, "_SSR_")
  } else {
    sd.name <- if (sigma.fixed) NULL else "sigma_"
    SSR.name <- "SSR_"
  }
  sc <- sm[["control"]]
  store_default <- function() {
    #out <- if (var.prior[["type"]] == "fixed") NULL else sd.name
    out <- sd.name
    if (modeled.Q) for (mc in Vmod) out <- c(out, mc[["store.default"]])
    out
  }
  if (!prior.only) {
    if (is.logical(y) || is.integer(y))
      y <- as.numeric(y)
    else if (!is.numeric(y)) stop("non-numeric response variable")
    if (var.prior[["type"]] == "invchisq") var.prior$make_draw()
  }
  if (!is.null(prec.mat)) {
    Q0 <- economizeMatrix(prec.mat, symmetric=TRUE, check=TRUE)
    if (!identical(dim(Q0), c(n, n))) stop("incompatible precision matrix")
    prec.mat <- "prec.mat.used"
  } else if (is.vector(var.vec)) {
    # allow vector input for backward compatibility
    if (!is.numeric(var.vec) || anyNA(var.vec) || any(var.vec <= 0)) stop("invalid input for 'var.vec'")
    if (!is.null(names(var.vec))) names(var.vec) <- NULL
    Q0 <- Cdiag(if (length(var.vec) == 1L) rep.int(1/var.vec, n) else 1/var.vec)
  } else if (intercept_only(var.vec)) {
    Q0 <- CdiagU(n)
  } else {
    temp <- get_var_from_formula(var.vec, data)
    if (!is.null(names(temp))) names(temp) <- NULL
    if (any(temp <= 0)) stop("non-positive variance(s) in 'var.vec'")
    if (length(temp) == 1L) {
      Q0 <- Cdiag(rep.int(1/temp, n))
    } else {
      if (length(temp) != n) stop("variance vector has wrong length")
      Q0 <- Cdiag(1/temp)
    }
  }
  if (isDiagonal(Q0)) {
    if (is_unit_diag(Q0))
      Q0.type <- "unit"
    else
      Q0.type <- "diag"
  } else {
    Q0.type <- "symm"  # non-diagonal precision matrix
  }
  if (prior.only) {
    scale.e <- 1
  } else {
    scale.e <- 0.5 * sd(y)
    if (scale.e == 0 || !is.finite(scale.e)) scale.e <- 1
    # Q_e computes Q times the partial residual for a model component
    if (modeled.Q) {
      if (multifam) {
        # factor 1/sigma^2 included in p$Q_ or p$QM_
        if (sc[["single.block"]])
          Q_e <- switch(Q0.type,
            unit=, diag = function(p) p[["Q_"]][sub] * y,
            symm = function(p) p[["QM_"]][sub, sub] %m*v% y
          )
        else
          Q_e <- switch(Q0.type,
            unit=, diag = function(p) p[["Q_"]][sub] * (y - p[["e_"]][sub]),
            # in case of non-diagonal Q0, full precision matrix stored as p[["QM_"]]
            symm = function(p) p[["QM_"]][sub, sub] %m*v% (y - p[["e_"]][sub])
          )
        compute_SSR <- switch(Q0.type,
          unit=, diag = function(p) {
            res <- y - p[["e_"]][sub]
            dotprodC(res, p[["Q_"]][sub] * res)
          },
          symm = function(p) {
            res <- y - p[["e_"]][sub]
            dotprodC(res, p[["QM_"]][sub, sub] %m*v% res)
          }
        )
      } else {
        # factor 1/sigma^2 not included in p$Q_ or p$QM_
        Q_e <- switch(Q0.type,
          unit=, diag = function(p) p[["Q_"]] * p[["e_"]],
          # in case of non-diagonal Q0, full precision matrix stored as p[["QM_"]]
          symm = function(p) p[["QM_"]] %m*v% p[["e_"]]
        )
        compute_SSR <- function(p) dotprodC(p[["e_"]], Q_e(p))
      }
    } else {
      if (multifam) {
        if (sc[["single.block"]]) {
          if (sigma.fixed) {
            Q_e <- switch(Q0.type,
              unit = function(p) y,
              diag = function(p) Q0@x * y,
              symm = function(p) Q0 %m*v% y
            )
          } else {
            Q_e <- switch(Q0.type,
              unit = function(p) (1/p[[sd.name]]^2) * y,
              diag = function(p) (1/p[[sd.name]]^2) * Q0@x * y,
              symm = function(p) (1/p[[sd.name]]^2) * (Q0 %m*v% y)
            )
          }
        } else {
          if (sigma.fixed) {
            Q_e <- switch(Q0.type,
              unit = function(p) y - p[["e_"]][sub],
              diag = function(p) Q0@x * (y - p[["e_"]][sub]),
              symm = function(p) Q0 %m*v% (y - p[["e_"]][sub])
            )
          } else {
            Q_e <- switch(Q0.type,
              unit = function(p) (1/p[[sd.name]]^2) * (y - p[["e_"]][sub]),
              diag = function(p) (1/p[[sd.name]]^2) * Q0@x * (y - p[["e_"]][sub]),
              symm = function(p) (1/p[[sd.name]]^2) * (Q0 %m*v% (y - p[["e_"]][sub]))
            )
          }
        }
        # SSR for multifam; note that SSR uses full residuals whereas Q_e uses partial residuals
        # SSR is always computed without 1/sigma^2 factor
        compute_SSR <- switch(Q0.type,
          unit = function(p) {
            e <- y - p[["e_"]][sub]
            dotprodC(e, e)
          },
          diag = function(p) {
            e <- y - p[["e_"]][sub]
            dotprodC(e, Q0@x * e)
          },
          symm = function(p) {
            res <- y - p[["e_"]][sub]
            dotprodC(res, Q0 %m*v% res)
          }
        )
      } else {
        Q_e <- switch(Q0.type,
          unit = function(p) p[["e_"]],
          diag = function(p) Q0@x * p[["e_"]],
          symm = function(p) Q0 %m*v% p[["e_"]]
        )
        compute_SSR <- function(p) dotprodC(p[["e_"]], Q_e(p))
      }
    }
  }
  # scale.e used to generate default starting values for residuals (or fitted values for non-gaussian models)
  scale.sigma <- scale.e * fmean.default(sqrt(diag(Q0)), na.rm=FALSE)

  if (!prior.only) {
    draw <- function(p) {}
    draw <- add(draw, bquote(p[[.(SSR.name)]] <- compute_SSR(p)))
    if (sc[["compute.llh"]]) {
      if (multifam)
        draw <- add(draw, quote(p$llh_ <- p[["llh_"]] + llh(p)))
      else
        draw <- add(draw, quote(p$llh_ <- llh(p)))
    }
    if (!sigma.fixed) {
      df.sigma <- n  # likelihood contribution; any prior contributions are added by reg, gen etc components
      # any prior contributions to SSR_ are added by the model components' draw functions
      # so that the correct total SSR_ is used for drawing sigma
      draw <- add(draw, bquote(
        for (mc in sm[["mod"]]) if (is.function(mc[["SSR_sigma"]])) p[[.(SSR.name)]] <- p[[.(SSR.name)]] + mc$SSR_sigma(p)
      ))
    }
  }
  if (!sigma.fixed || modeled.Q) {
    rprior <- if (sigma.fixed)
      function(p) {}
    else
      function(p) {p[[sd.name]] <- sqrt(var.prior$rprior())}
    if (!prior.only) {
      start <- function(p) {}
      # sufficient statistics SSR and df.sigma used for sigma draws, contain contributions from
      #   likelihood as well as coefficient priors; only SSR is updated in each draw
      if (!sigma.fixed) {
        start <- add(start, bquote(p[[.(SSR.name)]] <- n * scale.sigma^2))
        switch(var.prior[["type"]],
          fixed = {
            draw <- add(draw, bquote(p[[.(sd.name)]] <- sqrt(var.prior[["value"]])))
          },
          invchisq = {
            if (is.list(var.prior[["scale"]]))
              draw <- add(draw, bquote(p[[.(sd.name)]] <- sqrt(var.prior$draw(df.sigma, p[[.(SSR.name)]], 1 / p[[.(sd.name)]]^2))))
            else
              draw <- add(draw, bquote(p[[.(sd.name)]] <- sqrt(var.prior$draw(df.sigma, p[[.(SSR.name)]]))))
          },
          exp = {
            rgig <- GIGrvg::rgig
            draw <- add(draw, bquote(p[[.(sd.name)]] <- sqrt(rgig(1L, 1 - 0.5*df.sigma, p[[.(SSR.name)]], 2/var.prior[["scale"]]))))
          },
          gig = {
            rgig <- GIGrvg::rgig
            draw <- add(draw, bquote(p[[.(sd.name)]] <- sqrt(rgig(1L, var.prior[["p"]] - 0.5*df.sigma, var.prior[["b"]] + p[[.(SSR.name)]], var.prior[["a"]]))))
          }
        )
        if (var.prior[["type"]] == "fixed") {
          start <- add(start, bquote(
            p[[.(sd.name)]] <- check_and_get(p, .(sd.name), 1L,
              \() sqrt(var.prior[["value"]]), pos=TRUE)
          ))
        } else {
          start <- add(start, bquote(
            p[[.(sd.name)]] <- check_and_get(p, .(sd.name), 1L,
              \() runif(1L, 0.1 * scale.sigma, scale.sigma), pos=TRUE)
          ))
        }
      }
    }
  }

  self <- environment()

  if (modeled.Q) {
    var.model <- standardise_formula(var.model, "reg", data=data)  # pass data to interpret '.'
    Vmod <- to_mclist(var.model, prefix="v")
    if (!length(Vmod)) stop("empty 'var.model'")
    types <- get_types(Vmod)
    if (any(types %in% c("mec", "brt"))) stop("'mec' and 'brt' can only be used in mean model specification")
    if (is.logical(control[["block.V"]])) {
      if (control[["block.V"]]) {
        # all components of type reg and gen in a single block
        control$block.V <- list(names(Vmod)[types %in% c("reg", "gen")])
        # a single component by default not handled as a block
        if (length(control[["block.V"]][[1L]]) <= 1L) control$block.V <- NULL
      } else {
        control$block.V <- NULL
      }
    } else {
      for (bl in control[["block.V"]]) {
        if (!all(bl %in% names(Vmod)))
          stop("invalid name(s) '", paste0(setdiff(bl, names(Vmod)), collapse="', '"), "' in 'block.V'")
        if (any(types[bl] %in% c("vreg", "vfac"))) stop("'vreg' and 'vfac' components cannot be part of a Gibbs block")
      }
      if (any_duplicated(unlst(control[["block.V"]]))) stop("duplicate model component name in 'block.V'")
    }
    single.V.block <- any(length(Vmod) == c(1L, length(unlst(control[["block.V"]]))))

    for (k in seq_along(Vmod)) {
      mc <- Vmod[[k]]
      mc$name <- names(Vmod)[k]
      mc$sc <- sc
      mc$fam <- self
      mc$in.block <- any(mc[["name"]] == unlst(control[["block.V"]]))
      mc$prior.only <- prior.only
      mc$data <- data
      mc <- as.list(mc)[-1L]
      Vmod[[mc[["name"]]]] <- do.call(
        getFromNamespace(paste0("mc_", types[k]), "mcmcsae"),
        mc, envir=environment(var.model)
      )
    }
    for (k in seq_along(Vmod))
      rprior <- add(rprior, bquote(p <- Vmod[[.(k)]]$rprior(p)))
    # compute product of precision factors
    compute_Qfactor <- function(p) {
      out <- Vmod[[1L]]$compute_Qfactor(p)
      for (mc in Vmod[-1L]) out <- out * mc$compute_Qfactor(p)
      out
    }
    # compute data-level precision matrix from Q0 and scale factor computed by compute_Qfactor
    # p$Q_ by default in store.mean for use in compute_DIC
    if (multifam) {
      compute_Q <- switch(Q0.type,
        unit = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_[sub] <- compute_Qfactor(p)
          p
        },
        diag = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_[sub] <- Q0@x * compute_Qfactor(p)
          p
        },
        symm = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) Qfactor <- compute_Qfactor(p)
          p$QM_[sub, sub] <- block_scale_dsCMatrix(Q0, Qfactor)  # store full precision matrix
          p$Q_[sub] <- Qfactor
          p
        }
      )
    } else {
      compute_Q <- switch(Q0.type,
        unit = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_ <- compute_Qfactor(p)
          p
        },
        diag = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) p$Q_ <- Q0@x * compute_Qfactor(p)
          p
        },
        symm = function(p, Qfactor=NULL) {
          if (is.null(Qfactor)) Qfactor <- compute_Qfactor(p)
          p$QM_ <- block_scale_dsCMatrix(Q0, Qfactor)  # store full precision matrix
          p$Q_ <- Qfactor
          p
        }
      )
    }

    if (!prior.only) {
      MHpars <- NULL
      adapt <- function(ar) {}
      for (k in seq_along(Vmod)) {
        mc <- Vmod[[k]]
        switch(mc[["type"]],
          vreg = MHpars <- c(MHpars, mc[["name"]]),
          vfac = {
            if (mc$prior[["type"]] == "invchisq" && is.list(mc$prior[["df"]])) {
              MHpars <- c(MHpars, mc[["name_df"]])
              if (mc$prior$df[["adapt"]])
                adapt <- add(adapt, bquote(Vmod[[.(k)]]$adapt(ar)))
            }
          },
          gen = {
            MHpars <- c(MHpars, if (mc[["usePX"]]) mc[["name_sigma_raw"]] else mc[["name_sigma"]])
            if (mc$control[["MHprop"]] == "LNRW")
              adapt <- add(adapt, bquote(Vmod[[.(k)]]$adapt(ar)))
          }
        )
        if (!(mc[["type"]] == "reg" && mc[["in.block"]])) {
          draw <- add(draw, bquote(p <- Vmod[[.(k)]]$draw(p)))
          start <- add(start, bquote(p <- Vmod[[.(k)]]$start(p)))
        }
      }  # END for (k in seq_along(Vmod))
      if (length(control[["block.V"]])) {
        mbs.V <- list()
        for (k in seq_along(control[["block.V"]])) {
          mbs.V[[k]] <- create_mc_block(Vmod[control[["block.V"]][[k]]], self, sc, prior.only=FALSE)
          draw <- add(draw, bquote(p <- mbs.V[[.(k)]]$draw(p)))
          start <- add(start, bquote(p <- mbs.V[[.(k)]]$start(p)))
        }
      }
      if (!single.V.block) {
        # recompute precision Q for numerical stability
        draw <- add(draw, quote(p <- compute_Q(p)))
      }
      start <- add(start, quote(p <- compute_Q(p)))
      if (length(body(adapt)) <= 1L) adapt <- NULL
    }
  }
  if (!sigma.fixed || modeled.Q) {
    rprior <- add(rprior, quote(p))
    if (!prior.only) {
      start <- add(start, quote(p))
    }
  }
  if (!prior.only) {
    if (multifam && !sigma.fixed) {
      if (modeled.Q)
        draw <- add(draw, bquote(p$Q_[sub] <- p[["Q_"]][sub] * (1 / p[[.(sd.name)]]^2)))
      else
        draw <- add(draw, bquote(p$Q_[sub] <- 1 / p[[.(sd.name)]]^2))
    }
    draw <- add(draw, quote(p))
  }
  if (!prior.only && sc[["cMVN.or.CG"]]) {
    # set up a function that multiplies by L Chol factor of Q, for sampling from N(., Q)
    if (modeled.Q) {
      cholQ <- switch(Q0.type,
        unit=, diag = build_chol(runif(n, 0.9, 1.1)),
        symm = build_chol(crossprod_sym(Cdiag(runif(0.9, 1.1)), Q0), control=sc[["chol.control"]])
      )
    } else {
      cholQ <- switch(Q0.type,
        unit = build_chol(CdiagU(n)),
        diag = build_chol(Cdiag(Q0@x)),
        symm = build_chol(Q0, control=sc[["chol.control"]])
      )
    }
    # draw from MVN with variance(!) Q
    if (modeled.Q) {
      drawMVNvarQ <- function(p) {
        if (Q0.type == "symm")
          cholQ$update(p[["QM_"]])
        else
          cholQ$update(p[["Q_"]])
        cholQ$Ltimes(Crnorm(n), transpose=FALSE)
      }
    } else {
      drawMVNvarQ <- function(p) cholQ$Ltimes(Crnorm(n), transpose=FALSE)
    }
  }
  pred_pars <- function() {
    out <- if (sigma.fixed) NULL else sd.name
    if (modeled.Q) for (mc in Vmod) out <- c(out, mc[["name"]])
    out
  }
  if (!prior.only) {
    llh0 <- -0.5 * n * log(2*pi)  # constant term of log-likelihood
    if (!is.null(logJacobian)) {
      logJacobian <- as.numeric(logJacobian)
      if (length(logJacobian) != n) stop("'logJacobian' should be a vector of the same size as the response vector")
      if (anyNA(logJacobian)) stop("missing values in 'logJacobian'")
      llh0 <- llh0 + sum(logJacobian)
    }
    if (!modeled.Q || Q0.type == "symm")
      llh0 <- llh0 + 0.5 * as.numeric(determinant(1 * Q0, logarithm=TRUE)$modulus)  # 1 * Q0 to ensure no chol object is stored with Q0 and possibly other refs to Q0
    if (modeled.Q) {
      if (multifam) {
        if (sigma.fixed)
          llh <- function(p) llh0 + 0.5 * sum(log(p[["Q_"]][sub])) - 0.5 * p[[SSR.name]]
        else
          llh <- function(p) llh0 + 0.5 * sum(log(p[["Q_"]][sub])) - n * log(p[[sd.name]]) - 0.5 * p[[SSR.name]] / p[[sd.name]]^2
      } else {
        if (sigma.fixed)
          llh <- function(p) llh0 + 0.5 * sum(log(p[["Q_"]])) - 0.5 * p[[SSR.name]]
        else
          llh <- function(p) llh0 + 0.5 * sum(log(p[["Q_"]])) - n * log(p[[sd.name]]) - 0.5 * p[[SSR.name]] / p[[sd.name]]^2
      }
    } else {
      if (sigma.fixed)
        llh <- function(p) llh0 - 0.5 * p[[SSR.name]]
      else
        llh <- function(p) llh0 - n * log(p[[sd.name]]) - 0.5 * p[[SSR.name]] / p[[sd.name]]^2
    }
    # for WAIC computation: compute log-likelihood for each observation/batch of observations, vectorized over parameter draws and observations
    llh0_i <- -0.5 * log(2*pi)
    llh_i <- function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      all.units <- length(i) == n
      if (sigma.fixed) {
        q <- 1
      } else {
        sigma <- as.numeric(as.matrix.dc(draws[[sd.name]], colnames=FALSE))
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
        llh0_i + 0.5 * ( log(q) - q * e_i^2 )
      else
        llh0_i + matrix(rep_each(logJacobian[i], nr), nr, length(i)) + 0.5 * ( log(q) - q * e_i^2 )
    }
  }
  # make_get_sds used for prediction; allow 0 variances
  make_get_sds <- function(newdata) {
    if (identical(prec.mat, "prec.mat.used")) stop("out-of-sample prediction not supported if 'prec.mat' is used")
    if (intercept_only(var.vec)) {
      sds0 <- 1
    } else {
      sds0 <- get_var_from_formula(var.vec, newdata)
      if (any(sds0 < 0)) stop("negative variance(s) in 'var.vec' for prediction")
      if (all(length(sds0) != c(1L, nrow(newdata)))) stop("wrong length for prediction variance vector")
      sds0 <- sqrt(sds0)
    }
    if (modeled.Q) {
      V <- list()
      for (Vmc in Vmod) {
        if (any(Vmc[["type"]] == c("vfac", "vreg")))
          V[[Vmc[["name"]]]] <- Vmc$make_predict_Vfactor(newdata)
        else
          V[[Vmc[["name"]]]] <- Vmc$make_predict(newdata)
      }
      function(p) {
        var <- 1
        for (Vmc in Vmod) {
          if (any(Vmc[["type"]] == c("vfac", "vreg")))
            var <- var * V[[Vmc[["name"]]]](p)
          else
            var <- var * exp(V[[Vmc[["name"]]]](p))
        }
        sds0 * sqrt(var)
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
        sigma <- if (sigma.fixed) 1 else p[[sd.name]]
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
      if (is.null(weights)) {
        if (sigma.fixed)
          function(p, lp) lp + get_sds(p) * Crnorm(nn)
        else
          function(p, lp) lp + p[[sd.name]] * get_sds(p) * Crnorm(nn)
      } else {
        if (sigma.fixed)
          function(p, lp) weights * lp + sqrt(weights) * get_sds(p) * Crnorm(nn)
        else
          function(p, lp) weights * lp + p[[sd.name]] * sqrt(weights) * get_sds(p) * Crnorm(nn)
      }
    }
  }
  rm(data)
  self
}

#' Set computational options for the sampling algorithms
#'
#' @export
#' @param block.V if \code{TRUE}, the default, all coefficients of \code{reg}
#'  and \code{gen} components in a variance model formula are sampled in a
#'  single block. Alternatively, a list of character vectors with names of
#'  model components whose coefficients should be sampled together in blocks.
#' @returns A list with computational options.
gaussian_control <- function(block.V = TRUE) {
  list(block.V = block.V)
}

check_gaussian_control <- function(control) {
  if (is.null(control)) return(gaussian_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- gaussian_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  if (is.logical(control[["block.V"]])) {
    if (length(control[["block.V"]]) != 1L) stop("unexpected input for 'block.V'")
  } else {
    if (!is.list(control[["block.V"]])) stop("'block.V' should be either a scalar logical or a list of variance model component name vectors")
    if (!length(control[["block.V"]])) stop("'block.V' must contain at least one character vector")
  }
  control
}
