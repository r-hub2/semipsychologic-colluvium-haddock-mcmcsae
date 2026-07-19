#' Specify a Student-t sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Student-t sampling distribution.
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
#' dat$y <- with(dat, 1 - 0.5*x + v[g]) + 0.4 * rt(n, df=2.5)
#' 
#' sampler <- create_sampler(
#'   y ~ x + (1|g), data=dat, family="student_t"
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' compute_DIC(sim)
#' waic(sim)
#' summary(sim)
#' 
#' #' # more explicit specification, allowing non-default names, priors etc.
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") + gen(~1, factor = ~ g, name="v"),
#'   data=dat,
#'   family = f_student_t(df.prior=pr_gamma(2, 0.2))
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' loo(sim)
#' summary(sim)
#' bayesplot::mcmc_recover_hist(as.array(sim$student_t_df), 2.5)
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), c(1, -0.5))
#' bayesplot::mcmc_recover_hist(as.array(sim$v_sigma), 0.6)
#' bayesplot::mcmc_recover_scatter(as.array(sim$v), v)
#' }
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for a Student-t distribution is \code{"identity"}.
#' @param df.prior prior specification for the (scalar) number of degrees of
#'  freedom of the Studen-t sampling distribution. Currently allowed priors are
#'  gamma priors specified using \code{\link{pr_gamma}}, or a degenerate prior
#'  fixing the degrees of freedom parameter at a fixed value.
#' @param control sampling options for the degrees of freedom parameter
#'  (if not fixed).
#' @param ... further arguments passed to \code{\link{f_gaussian}}. Note that
#'  \code{var.prior}, \code{var.vec} and \code{var.model} now refer to the
#'  Student-t scale parameters.
#' @returns A family object.
f_student_t <- function(link="identity", df.prior=pr_gamma(2, 0.1),
                        control=student_t_control(), ...) {
  link <- match.arg(link)
  df.prior$init(n=1L)
  control <- check_student_t_control(control)
  dotargs <- list(...)
  # TODO
  # - allow PX
  # - allow (fixed) vector df
  # - pass on gaussian_control options
  #c(list(family="student_t", df.prior=df.prior, var.model=var.model), lgau)
  c(list(family="student_t", link=link, df.prior=df.prior, control=control, `_raw_`=TRUE), dotargs)
}

ff_student_t <- function(link="identity", df.prior=pr_gamma(2, 0.1),
                         control=student_t_control(),
                         sm, data, y=NULL, famid=NULL, sub=NULL, ...) {
  family <- "student_t"
  linkinv <- make.link(link)$linkinv
  prior.only <- is.null(y)
  sc <- sm[["control"]]
  multifam <- !is.null(famid)
  name.t.vfac <- if (multifam) paste0(famid, "_student_t") else "student_t"
  dotargs <- list(...)
  if (df.prior[["type"]] == "fixed") {
    var.model <- eval(bquote(~ vfac(prior=pr_invchisq(df=df.prior[["value"]]), name=.(name.t.vfac))))
  } else if (df.prior[["type"]] == "gamma") {
    name.df <- paste0(name.t.vfac, "_df")
    var.model <- eval(bquote(~ vfac(factor="local_", prior=pr_invchisq(df=list(
        alpha=df.prior[["shape"]], beta0=df.prior[["rate"]],
        proposal=control[["proposal"]], tau=control[["MHscale"]], adapt=control[["adapt"]]
      )), name=.(name.t.vfac))
    ))
  } else {
    stop("unsupported prior for Student-t degrees of freedom parameter")
  }
  lgau <- do.call(f_gaussian, dotargs)
  if (!is.null(lgau[["prec.mat"]])) {
    stop("argument 'prec.mat' for specifying a non-diagonal sampling ",
         "covariance matrix, is not supported for Student-t family")
  }
  if (!is.null(lgau[["var.model"]])) {
    var.model <- as.formula(
      call("~", call("+", var.model[[2L]], lgau[["var.model"]][[2L]])),
      env=environment(lgau[["var.model"]])
    )
  }
  y.family <- ff_gaussian(link=link, var.prior=lgau[["var.prior"]],
    var.vec=lgau[["var.vec"]], prec.mat=NULL, var.model=var.model,
    logJacobian=if (prior.only) NULL else lgau[["logJacobian"]],
    #control=dotargs[["gaussian.control"]],
    sm=sm, data=data, y=y, famid=famid, sub=sub
  )
  rm(lgau)
  store_default <- function() y.family$store_default()
  if (is.function(y.family[["adapt"]]))
    adapt <- function(ar) y.family$adapt(ar)
  if (is.function(y.family[["drawMVNvarQ"]])) {
    cholQ <- y.family[["cholQ"]]
    drawMVNvarQ <- y.family[["drawMVNvarQ"]]
  }
  n <- y.family[["n"]]
  Q0 <- y.family[["Q0"]]
  Q0.type <- y.family[["Q0.type"]]
  sd.name <- y.family[["sd.name"]]
  sigma.fixed <- y.family[["sigma.fixed"]]
  Vmod <- y.family[["Vmod"]]
  modeled.Q <- y.family[["modeled.Q"]]
  logJacobian <- y.family[["logJacobian"]]
  self <- environment()
  copy_refs(y.family, self, c(
    "y", "e.is.res", "rprior", "SSR.name", "compute_Q"
  ))
  t.modeled.Q <- modeled.Q && length(Vmod) > 1L
  if (t.modeled.Q) {
    Vmod.t <- Vmod[-which(names(Vmod) == name.t.vfac)]
    compute_Qfactor <- function(p) {
      out <- Vmod.t[[1L]]$compute_Qfactor(p)
      for (mc in Vmod.t[-1L]) out <- out * mc$compute_Qfactor(p)
      out
    }
    environment(compute_Q) <- self
  }
  vfac.t <- Vmod[[which(names(Vmod) == name.t.vfac)]]
  df.fixed <- !is.list(vfac.t$prior[["df"]])
  if (df.fixed) df.value <- vfac.t$prior[["df"]]
  pred_pars <- function() {
    out <- if (sigma.fixed) NULL else sd.name
    if (!df.fixed) out <- c(out, vfac.t[["name_df"]])
    if (t.modeled.Q) for (mc in Vmod.t) out <- c(out, mc[["name"]])
    out
  }
  if (!prior.only) {
    draw <- y.family[["draw"]]  # avoid check NOTE
    copy_refs(y.family, self, c(
      "Q_e", "scale.e", "scale.sigma", "start",
      "MHpars", "compute_SSR"
    ))
    if (df.fixed)
      llh0 <- n*(lgamma(0.5*(df.value + 1)) - lgamma(0.5*df.value) - 0.5*log(pi*df.value)) 
    else
      llh0 <- 0
    if (!is.null(logJacobian)) llh0 <- llh0 + sum(logJacobian)
    llh <- function(p) {
      q <- if (sigma.fixed) 1 else 1/p[["sigma_"]]^2
      if (t.modeled.Q) for (mc in Vmod.t) q <- q * mc$compute_Qfactor(p)
      if (Q0.type == "diag") q <- q * Q0@x
      sumlogq <- if (length(q) == 1) n*log(q) else sum(log(q))
      e2 <- if (multifam) p[["e_"]][sub]^2 else p[["e_"]]^2
      if (df.fixed) {
        llh0 + 0.5 * sumlogq - 0.5 * (df.value + 1) * sum(log1p((q * e2) / df.value))
      } else {
        df <- p[[name.df]]
        llh0 + n*(lgamma(0.5*(df + 1)) - lgamma(0.5*df) - 0.5*log(pi*df)) + 0.5 * sumlogq - 0.5 * (df + 1) * sum(log1p((q/df) * e2))
      }
    }
    assign("llh", llh, environment(draw))  # draw was defined by ff_gaussian so we should replace its llh function
    # for WAIC computation: compute log-likelihood for each observation/batch of observations, vectorized over parameter draws and observations
    # use Rao-Blackwellised llh0_i instead of gaussian family's llh_i
    if (df.fixed) {
      llh0_i <- lgamma(0.5*(df.value + 1)) - lgamma(0.5*df.value) - 0.5*log(pi*df.value) 
    }
    llh_i <- function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      all.units <- length(i) == n
      if (df.fixed)
        df <- df.value
      else
        df <- as.numeric(as.matrix.dc(draws[[name.df]], colnames=FALSE))
      if (sigma.fixed) {
        q <- 1
      } else {
        sigma <- as.numeric(as.matrix.dc(draws[[sd.name]], colnames=FALSE))
        q <- 1/sigma^2
      }
      if (t.modeled.Q) {
        q <- matrix(q, nr, length(i))
        for (mc in Vmod.t) {
          if (all.units) Xi <- mc[["X"]] else Xi <- mc[["X"]][i, , drop=FALSE]
          if (mc[["type"]] == "vfac")
            q <- q * tcrossprod(1 / as.matrix.dc(draws[[mc[["name"]]]], colnames=FALSE), Xi)
          else
            q <- q * exp(-tcrossprod(as.matrix.dc(draws[[mc[["name"]]]], colnames=FALSE), Xi))
        }
      }
      if (Q0.type == "diag") q <- q * rep_each(Q0@x[i], each=nr)
      if (df.fixed) {
        out <- llh0_i + 0.5 * log(q) - 0.5 * (df + 1) * log1p((q * e_i^2) / df)
      } else {
        out <- lgamma(0.5*(df + 1)) - lgamma(0.5*df) - 0.5*log(pi*df) + 0.5 * log(q) - 0.5 * (df + 1) * log1p((q * e_i^2) / df)
      }
      if (is.null(logJacobian))
        out
      else
        matrix(rep_each(logJacobian[i], nr), nr, length(i)) + out
    }
  }
  # make_get_sds used for prediction; allow 0 variances
  make_get_sds <- function(newdata) {
    if (intercept_only(y.family[["var.vec"]])) {
      sds0 <- 1
    } else {
      sds0 <- get_var_from_formula(y.family[["var.vec"]], newdata)
      if (any(sds0 < 0)) stop("negative variance(s) in 'var.vec' for prediction")
      if (all(length(sds0) != c(1L, nrow(newdata)))) stop("wrong length for prediction variance vector")
      sds0 <- sqrt(sds0)
    }
    if (t.modeled.Q) {
      V <- list()
      for (Vmc in Vmod.t) {
        if (any(Vmc[["type"]] == c("vfac", "vreg")))
          V[[Vmc[["name"]]]] <- Vmc$make_predict_Vfactor(newdata)
        else
          V[[Vmc[["name"]]]] <- Vmc$make_predict(newdata)
      }
      function(p) {
        var <- 1
        for (Vmc in Vmod.t) {
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
  # Q0.type "symm" currently not supported (needs whitening transformation)
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
      function(p, lp) {
        Qfac <- if (Q0.type == "diag") Q0@x else 1
        if (t.modeled.Q) {
          if (is.null(p[["Q_"]])) p <- compute_Q(p)
          Qfac <- Qfac * p[["Q_"]]
        }
        df <- if (df.fixed) df.value else p[[name.df]]
        if (is.null(weights)) {
          if (sigma.fixed)
            lp + invsqrt(Qfac) * rt(nn, df)
          else
            lp + p[[sd.name]] * invsqrt(Qfac) * rt(nn, df)
        } else {
          if (sigma.fixed)
            weights * lp + sqrt(weights) * invsqrt(Qfac) * rt(nn, df)
          else
            weights * lp + p[[sd.name]] * sqrt(weights) * invsqrt(Qfac) * rt(nn, df)
        }
      }
    } else {
      nn <- nrow(newdata)
      get_sds <- make_get_sds(newdata)
      get_df <- if (df.fixed) function(p) df.value else function(p) p[[name.df]]
      if (is.null(weights)) {
        if (sigma.fixed)
          function(p, lp) lp + get_sds(p) * rt(nn, get_df(p))
        else
          function(p, lp) lp + p[[sd.name]] * get_sds(p) * rt(nn, get_df(p))
      } else {
        if (sigma.fixed)
          function(p, lp) weights * lp + sqrt(weights) * get_sds(p) * rt(nn, get_df(p))
        else
          function(p, lp) weights * lp + p[[sd.name]] * sqrt(weights) * get_sds(p) * rt(nn, get_df(p))
      }
    }
  }
  rm(data, sm)
  self
}

#' Set (MCMC) computational options for a Student-t likelihood model
#'
#' @export
#' @param proposal "RW" for random walk Metropolis-Hastings or "mala"
#'  for Metropolis-adjusted Langevin.
#' @param MHscale (starting) scale of Metropolis-Hastings update.
#' @param adapt whether to adapt the scale of the proposal distribution
#'  during burnin to achieve better acceptance rates.
#' @returns A list with computational options.
student_t_control <- function(proposal="RW", MHscale=1, adapt=TRUE) {
  list(proposal=proposal, MHscale=MHscale, adapt=adapt)
}

check_student_t_control <- function(control) {
  if (is.null(control)) return(student_t_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- student_t_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$proposal <- match.arg(control[["proposal"]], c("RW", "mala"))
  if (!is_numeric_scalar(control[["MHscale"]]) || control[["MHscale"]] <= 0)
    stop("'MHscale' must be a positive numeric value")
  if (!is_logical_scalar(control[["adapt"]])) stop("'adapt' must be TRUE or FALSE")
  control
}
