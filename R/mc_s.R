#' Specify a smooth term component of the linear predictor
#'
#' This function can be used inside the formula specification
#' of the linear predictor in \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. The smooth term is set up by
#' the smooth term specification function \code{\link[mgcv]{s}}
#' of package \pkg{mgcv}. The smooth term is usually composed
#' of random (penalised) effects as well as a few fixed
#' (unpenalised) effects, not including an intercept.
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' set.seed(0)
#' dat <- gamSim(5, n=200, scale=2)
#' b <- gam(y ~ x0 + s(x1) + s(x2) + s(x3), data=dat)
#'
#' sampler <- create_sampler(
#'   y ~ x0 + s(x1) + s(x2) + s(x3), data=dat
#' )
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' (summ <- summary(sim))
#' plot(
#'   coef(b),
#'   c(summ$reg1[, "Mean"],
#'     summ$s2_r[, "Mean"], summ$s2_f[, "Mean"],
#'     summ$s3_r[, "Mean"], summ$s3_f[, "Mean"],
#'     summ$s4_r[, "Mean"], summ$s4_f[, "Mean"]
#'    )
#' ); abline(0, 1)
#' predb <- predict(b, newdata=dat[1:5, ])
#' pred <- predict(sim, newdata=dat[1:5, ], type="response")
#' (summpred <- summary(pred))
#' plot(predb, summpred[, "Mean"]); abline(0, 1)
#' }
#'
#' @export
#' @param ... parameters passed to \code{mgcv::s}. Note that
#'  variables appearing in \code{...} must be present in the data frame
#'  passed to the data argument of \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}.
#' @param unit.precision to be implemented.
#' @param name the name of the model component. By default  the name will be 's'
#'  with the number of the model term attached. This name is used in the output
#'  of the MCMC simulation function \code{\link{MCMCsim}}. The name is appended
#'  by '_f' for any unpenalised (fixed) effects, if any, and by '_r' for the
#'  penalised (random) effects.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component,
#'  intended for internal use by other package functions.
#' @references
#'  S.N. Wood (2017).
#'    Generalized additive models: an introduction with R.
#'    Chapman and Hall/CRC.
s <- function(..., unit.precision=FALSE, name="", debug=FALSE) {
  stop("function 'mcmcsae::s' should only be used inside a formula")
}

# additional argument e to pass sampler environment, and in.block
mc_s <- function(..., unit.precision=FALSE, name="",
                 debug=FALSE, e, in.block) {
  type <- "s"
  if (name == "") stop("missing model component name")

  if (!requireNamespace("mgcv", quietly=TRUE)) stop("package mgcv required for a model including smooth terms defined through s()")
  sm <- mgcv::smoothCon(mgcv::s(...), data=e[["data"]], absorb.cons=TRUE)[[1L]]
  if (length(sm[["S"]]) > 1L) stop("smooth terms with multiple penalties are not supported")
  if (sm[["by"]] != "NA") stop("'by' argument of s() not yet supported")
  X <- sm[["X"]]
  Q0 <- economizeMatrix(sm[["S"]][[1L]], sparse=TRUE, symmetric=TRUE, drop.zeros=TRUE)

  if (unit.precision) {
    re <- mgcv::smooth2random(sm, "", type=2L)
    stop("unit.precision=TRUE transformation not currently supported")
  } else {
    if (nrow(X) != e[["n"]]) stop("design matrix with incompatible number of rows")
    e$coef.names[[name]] <- dimnames(X)[[2L]]
    # fixed effects model component
    qf <- sm[["null.space.dim"]]
    qr <- ncol(X) - qf
    q <- qf + qr
    name.r <- paste0(name, "_r")
    i.r <- seq_len(qr)
    if (qf > 0L) {
      name.f <- paste0(name, "_f")
      i.f <- (qr + 1L):q
      mf <- mc_reg(X = X[, i.f, drop=FALSE], name=name.f, e=e, in.block=in.block)
    }
    mr <- mc_gen(X = X[, i.r, drop=FALSE],
      factor = ~ custom(Q = Q0[i.r, i.r]),
      name=name.r, e=e, in.block=in.block
    )
    store.default <- mr[["store.default"]]
    if (qf > 0L) {
      store.default <- c(store.default, mf[["store.default"]])
      Q0 <- bdiag_ddidsC(list(mr[["Q"]], mf[["Q0"]]))
    }

    if (in.block && !e[["prior.only"]])
      get_Q <- function(p) {
        if (qf > 0L)
          c(mr$get_Q(p), mf$get_Q(p))
        else
          mr$get_Q(p)
      }

    lp_update <- function(x, plus=TRUE, p) {
      if (qf > 0L) mf$lp_update(x, plus, p)
      mr$lp_update(x, plus, p)
    }

    rprior <- function(p) {
      if (qf > 0L) p <- mf$rprior(p)
      mr$rprior(p)
    }

    # draws_linpred method acts on (subset of) mcdraws object, used in fitted() and pointwise log-likelihood llh_i functions
    draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
      out <- mr$draws_linpred(obj, units, chains, draws, matrix)
      if (qf > 0L) {
        outr <- mr$draws_linpred(obj, units, chains, draws, TRUE)
        if (matrix)
          out <- out + outr
        else
          for (ch in seq_along(chains))
            out[[ch]] <- out[[ch]] + outr[[ch]]
      }
      out
    }

    make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
      if (is.null(newdata) && is.null(Xnew)) {
        # in-sample prediction
        pred.r <- mr$make_predict(verbose=verbose)
        if (q > 0L) pred.f <- mf$make_predict(verbose=verbose)
      } else if (is.null(Xnew)) {
        Xnew <- mgcv::PredictMat(sm, data=newdata)
        pred.r <- mr$make_predict(Xnew=Xnew[, i.r, drop=FALSE], verbose=verbose)
        if (q > 0L) pred.f <- mf$make_predict(Xnew=Xnew[, i.f, drop=FALSE], verbose=verbose)
      }
      linpred <- function(p) {
        out <- pred.r$linpred(p)
        if (qf > 0L)
          pred.f$linpred_update(out, TRUE, p)
        else
          out
      }
      linpred_update <- function(x, plus=TRUE, p) {
        pred.r$linpred_update(x, plus, p)
        if (qf > 0L) pred.f$linpred_update(x, plus, p)
        x
      }
      environment()
    }

    if (!e[["prior.only"]]) {

      if (!e[["sigma.fixed"]]) {
        df.add <- mr[["df.add"]]
        if (qf > 0L) df.add <- df.add + mf[["df.add"]]
        SSR_add <- function(p)
          if (qf == 0L || is.null(mf$SSR_add))
            mr$SSR_add(p)
          else
            mf$SSR_add(p) + mr$SSR_add(p)
      }

      draw <- function(p) {
        if (qf > 0L && !in.block) p <- mf$draw(p)
        mr$draw(p)
      }
      start <- function(p) {
        if (qf > 0L && !in.block) p <- mf$start(p)
        p <- mr$start(p)
        p
      }
    }

  }

  environment()
}
