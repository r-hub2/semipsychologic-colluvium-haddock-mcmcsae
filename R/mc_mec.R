#' Create a model component object for a regression (fixed effects) component
#' in the linear predictor with measurement errors in quantitative covariates
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates an additive regression term in the
#' model's linear predictor. Covariates are assumed to be measured subject
#' to normally distributed errors with zero mean and variance specified using
#' the \code{formula} or \code{V} arguments. Note that this means that \code{formula}
#' should only contain quantitative variables, and no intercept.
#' By default, the prior for the regression
#' coefficients is improper uniform. A proper normal prior can be set up
#' using function \code{\link{pr_normal}}, and passed to argument \code{prior}.
#' It should be noted that \code{\link{pr_normal}} expects a precision matrix
#' as input for its second argument, and that the prior variance (matrix) is
#' taken to be the inverse of this precision matrix, where in case the
#' model's family is \code{"gaussian"} this matrix is additionally
#' multiplied by the residual scalar variance parameter \code{sigma_^2}.
#'
#' @examples
#' \donttest{
#' # example of Ybarra and Lohr (2008)
#' m <- 50
#' X <- rnorm(m, mean=5, sd=3)  # true covariate values
#' v <- rnorm(m, sd=2)
#' theta <- 1 + 3*X + v  # true values
#' psi <- rgamma(m, shape=4.5, scale=2)
#' e <- rnorm(m, sd=sqrt(psi))  # sampling error
#' y <- theta + e  # direct estimates
#' C <- c(rep(3, 10), rep(0, 40))  # measurement error for first 10 values
#' W <- X + rnorm(m, sd=sqrt(C))  # covariate subject to measurement error
#'
#' # fit Ybarra-Lohr model
#' sampler <- create_sampler(
#'   y ~ 1 + mec(~ 0 + W, V=C) + gen(factor=~local_),
#'   Q0=1/psi, sigma.fixed=TRUE, linpred="fitted"
#' )
#' sim <- MCMCsim(sampler, n.iter=800, n.chain=2, store.all=TRUE, verbose=FALSE)
#' (summ <- summary(sim))
#' plot(X, W, xlab="true X", ylab="inferred X")
#' points(X, summ$mec2_X[, "Mean"], col="green")
#' abline(0, 1, col="red")
#' legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
#' }
#'
#' @export
#' @param formula a formula specifying the predictors subject to measurement error
#'  and possibly their variances as well. In the latter case the formula syntax
#'  \code{~ (x1 | V.x1) + (x2 | V.x2) + ...} should be used where \code{x1, x2, ...}
#'  are the names of (quantitative) predictors and \code{V.x1, V.x2, ...} are the names
#'  of the variables holding the corresponding measurement error variances.
#'  If only the predictors are specified
#'  the formula has the usual form \code{~ x1 + x2 + ...}. In that case variances
#'  should be specified using argument \code{V}.
#'  All variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param sparse whether the model matrix associated with \code{formula} should
#'  be sparse. The default is to base this on a simple heuristic.
#' @param X a (possibly sparse) design matrix can be specified directly, as an
#'  alternative to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
#' @param V measurement error variance; can contain zeros
# TODO in general case this can be a list of q x q precision matrices where q
#      is the number of covariates
#' @param prior prior specification for the regression coefficients. Currently only
#'  normal priors are supported, specified using function \code{\link{pr_normal}}.
#' @param Q0 prior precision matrix for the regression effects. The default is a
#'  zero matrix corresponding to a noninformative improper prior.
#'  It can be specified as a scalar value, as a numeric vector of appropriate
#'  length, or as a matrix object. DEPRECATED, please use argument \code{prior}
#'  instead, i.e. \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param b0 prior mean for the regression effect. Defaults to a zero vector.
#'  It can be specified as a scalar value or as a numeric vector of
#'  appropriate length. DEPRECATED, please use argument \code{prior}
#'  instead, i.e. \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param R optional constraint matrix for equality restrictions R'x = r where
#'  \code{x} is the vector of regression effects.
#' @param r right hand side for the equality constraints.
#' @param S optional constraint matrix for inequality constraints S'x >= s where
#'  x is the vector of regression effects.
#' @param s right hand side for the inequality constraints.
#' @param lower as an alternative to \code{s}, \code{lower} and \code{upper} may be specified
#'  for two-sided constraints lower <= S'x <= upper.
#' @param upper as an alternative to \code{s}, \code{lower} and \code{upper} may be specified
#'  for two-sided constraints lower <= S'x <= upper.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'reg'
#'  with the number of the model term attached.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
#' @references
#'  L.M. Ybarra and S.L. Lohr (2008).
#'    Small area estimation when auxiliary information is measured with error.
#'    Biometrika 95(4), 919-931.
#'
#'  S. Arima, G.S. Datta and B. Liseo (2015).
#'    Bayesian estimators for small area models when auxiliary information is measured with error.
#'    Scandinavian Journal of Statistics 42(2), 518-529.
mec <- function(formula = ~ 1, sparse=NULL, X=NULL, V=NULL,
                prior=NULL, Q0=NULL, b0=NULL,
                R=NULL, r=NULL, S=NULL, s=NULL, lower=NULL, upper=NULL,
                name="", debug=FALSE) {

  e <- sys.frame(-2L)
  type <- "mec"
  if (name == "") stop("missing model component name")
  remove.redundant <- FALSE

  if (is.null(X)) {
    # TODO: warn if formula contains explicit intercept or categorical terms
    vs <- as.list(attr(terms(formula), "variables")[-1L])
    if (all(lengths(vs) == 3L)) {
      if (!all(b_apply(vs, function(x) identical(x[[1L]], as.name("|"))))) stop("invalid 'formula' argument")
      V.in.formula <- TRUE
      formula.X <- as.formula(paste0("~ 0 + ", paste0(s_apply(vs, function(x) deparse(x[[2L]])), collapse=" + ")), env=environment(formula))
      if (grepl("|", as.character(formula.X)[2L], fixed=TRUE)) stop("invalid 'formula'")  # maybe '()' forgotten?
      formula.V <- as.formula(paste0("~ 0 + ", paste0(s_apply(vs, function(x) deparse(x[[3L]])), collapse=" + ")), env=environment(formula))
      X <- model_matrix(formula.X, e[["data"]], sparse=sparse)
      V <- model_matrix(formula.V, e[["data"]], sparse=sparse)
    } else {
      if (all(lengths(vs) <= 2L)) {
        V.in.formula <- FALSE
        formula.X <- as.formula(paste0("~ 0 + ", paste0(s_apply(vs, deparse), collapse=" + ")), env=environment(formula))
        X <- model_matrix(formula.X, e[["data"]], sparse=sparse)
      } else {
        stop("invalid 'formula' argument")
      }
    }
  } else {
    if (is.null(colnames(X))) colnames(X) <- seq_len(ncol(X))
    V.in.formula <- FALSE
  }
  if (nrow(X) != e[["n"]]) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- colnames(X)
  X <- economizeMatrix(X, sparse=sparse, strip.names=TRUE, check=TRUE)
  q <- ncol(X)
  in_block <- any(name == unlst(e$block))

  if (!V.in.formula) {
    if (is.null(V)) stop("no measurement error variances supplied in mec component")
    if (is.vector(V)) {
      if (q != 1L) stop("'V' should be a matrix")
      V <- matrix(V, ncol=1L)
    }
    if (nrow(V) != e[["n"]] || ncol(V) != q) stop("wrong size of variance matrix 'V'")
    # TODO allow n x q matrix (uncorrelated measurement errors), or list of n q x q symmetric matrices
  }
  # determine units with measurement error
  V <- unname(V)
  if (any(V < 0)) stop("negative measurement error variance(s)")
  i.me <- which(apply(V, 1L, function(x) all(x > 0)))  # TODO add tolerance + allow list of matrix
  nme <- length(i.me)
  if (nme == e[["n"]]) i.me <- 1:e[["n"]]  # use ALTREP (?)
  if (q == 1L) {
    QME <- 1 / V[i.me, ]
    rm(V)
  } else {
    V <- V[i.me, ]
  }

  if (e$Q0.type == "unit") {
    if (q == 1L) v0 <- 1 else v0 <- rep.int(1, nme)
  } else {
    if (e$Q0.type == "symm") stop("not supported: measurement error component and non-diagonal sampling variance")
    v0 <- 1/diag(e[["Q0"]])[i.me]
  }

  if (!is.null(b0) || !is.null(Q0)) {
    warn("arguments 'b0' and 'Q0' are deprecated; please use argument 'prior' instead")
    if (is.null(prior)) prior <- pr_normal(mean = if (is.null(b0)) 0 else b0, precision = if (is.null(Q0)) 0 else Q0)
  } else {
    if (is.null(prior)) prior <- pr_normal(mean=0, precision=0)
  }
  if (prior$type != "normal") stop("only a normal prior is currently supported for 'mec' model component effects")
  prior$init(q, e$coef.names[[name]], sparse=if (in_block) TRUE else NULL, sigma=!e$sigma.fixed)  # mc_block uses x-slot of precision matrix
  informative.prior <- prior$informative
  Q0 <- prior$precision
  zero.mean <- !informative.prior || all(prior$mean == 0)
  if (zero.mean) {
    prior$mean <- 0
    Q0b0 <- 0
  } else {
    if (length(prior$mean) == 1L)
      Q0b0 <- Q0 %m*v% rep.int(prior$mean, q)
    else
      Q0b0 <- Q0 %m*v% prior$mean
  }
  # TODO TMVN prior

  if (!zero.mean && e$compute.weights) stop("weights cannot be computed if some coefficients have non-zero prior means")

  name_X <- paste0(name, "_X")
  lp <- function(p) p[[name_X]] %m*v% p[[name]]
  lp_update <- function(x, plus=TRUE, p) mv_update(x, plus, p[[name_X]], p[[name]])
  # draws_linpred method acts on (subset of) mcdraws object, used in fitted() and pointwise log-likelihood llh_i functions
  draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
    Xi <- get_from(obj[[name_X]], vars=units, chains=chains, draws=draws)
    if (matrix) {
      out <- tcrossprod(as.matrix.dc(get_from(obj[[name]], chains=chains, draws=draws), colnames=FALSE), Xi)
    } else {
      out <- list()
      for (ch in seq_along(chains))
        out[[ch]] <- tcrossprod(get_from(obj[[name]], chains=chains[ch], draws=draws)[[1L]], Xi)
    }
    out
  }

  # function that creates an oos prediction function (closure) based on new data
  # this computes the contribution of this component to the linear predictor
  make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
    if (is.null(newdata)) {
      if (is.null(Xnew)) {
        # in-sample prediction
        Xnew <- X
        linpred <- function(p) p[[name_X]] %m*v% p[[name]]
        linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, p[[name_X]], p[[name]])
      } else {
        stop("for out-of-sample prediction with 'mec' components a data frame ",
             "must be supplied")
      }
    } else {
      if (!V.in.formula) stop("for out-of-sample prediction with 'mec' components make sure to specify ",
                              "the measurement error variances using the formula argument")
      Xnew <- model_matrix(formula.X, newdata, sparse=sparse)
      dimnames(Xnew) <- list(NULL, NULL)
      SEnew <- model_matrix(formula.V, newdata, sparse=sparse)
      dimnames(SEnew) <- list(NULL, NULL)
      if (any(SEnew < 0)) stop("negative measurement error variance(s) in 'newdata'")
      SEnew <- sqrt(SEnew)
      if (nrow(SEnew) < 0.5 * nrow(Xnew)) {
        linpred <- function(p) (Xnew + SEnew * Crnorm(length(SEnew))) %m*v% p[[name]]
        linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, Xnew + SEnew * Crnorm(length(SEnew)), p[[name]])
      } else {
        i.me.new <- which(apply(SEnew, 1L, function(x) all(x > 0)))  # TODO add tolerance + allow list of matrix
        if (length(i.me.new) == nrow(newdata)) i.me.new <- seq_len(nrow(newdata))
        SEnew <- SEnew[i.me.new, ]
        linpred <- function(p) {
          out <- Xnew
          out[i.me.new, ] <- out[i.me.new, ] + SEnew * Crnorm(length(SEnew))
          out %m*v% p[[name]]
        }
        linpred_update <- function(x, plus=TRUE, p) {
          out <- Xnew
          out[i.me.new, ] <- out[i.me.new, ] + SEnew * Crnorm(length(SEnew))
          mv_update(x, plus, out, p[[name]])
        }
      }
    }
    rm(newdata, verbose)
    environment()
  }

  rprior <- function(p) {}
  # X, QME are assumed to be part of the prior information
  if (nme == e[["n"]])
    rprior <- add(rprior, bquote(p[[.(name_X)]] <- X + sqrt(.(if (q == 1L) quote(1/QME) else quote(V))) * Crnorm(.(length(QME)))))
  else {
    rprior <- rprior |>
      add(bquote(p[[.(name_X)]] <- X)) |>
      add(bquote(p[[.(name_X)]][i.me, ] <- sqrt(.(if (q == 1L) quote(1/QME) else quote(V))) * Crnorm(.(length(QME)))))
  }
  rprior <- rprior |>
    add(bquote(p[[.(name)]] <- prior$rprior(p))) |>
    add(quote(p))


  if (!e$prior.only) {

    if (q > 1L) base_tcrossprod <- base::tcrossprod

    draw <- if (debug) function(p) {browser()} else function(p) {}
    if (e$single.block && length(e$mod) == 1L) {
      draw <- add(draw, quote(p$e_ <- e$y_eff()))
    } else {
      if (e$e.is.res)
        draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, p[[.(name_X)]], p[[.(name)]])))
      else
        draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, p[[.(name_X)]], p[[.(name)]])))
    }

    # 1. draw covariates p[[name_X]]
    # for now we assume that X is dense
    draw <- add(draw, bquote(scale <- .(if (e$sigma.fixed) quote(v0) else quote(v0 * p[["sigma_"]]^2))))
    if (q == 1L) {
      draw <- add(draw, bquote(Vscaled <- 1 / (QME * scale + p[[.(name)]]^2)))
      if (nme == e[["n"]]) {
        draw <- add(draw, bquote(p[[.(name_X)]] <- X + Vscaled * p[[.(name)]] * (p[["e_"]] - p[[.(name)]] * X) + sqrt(scale * Vscaled) * Crnorm(.(nme))))
      } else {
        draw <- add(draw, bquote(p[[.(name_X)]][i.me] <- X[i.me] + Vscaled * p[[.(name)]] * (p[["e_"]][i.me] - p[[.(name)]] * X[i.me]) + sqrt(scale * Vscaled) * Crnorm(.(nme))))
      }
    } else {
      # V-form (independent m.e.)
      draw <- add(draw, bquote(
        for (i in i.me) {
          Vi <- V[i, ]
          Cb <- Vi * p[[.(name)]]
          f <- 1 / (scale[i] + dotprodC(p[[.(name)]], Cb))
          Vx <- add_diagC(- f * base_tcrossprod(Cb), Vi)
          Xi <- X[i, ]
          p[[.(name_X)]][i, ] <- Xi + Cb * f * (p[["e_"]][i] - dotprodC(p[[.(name)]], Xi)) + crossprod_mv(chol.default(Vx), rnorm(q))
        }
      ))
      # TODO use Rcpp, and use Q-form at least for cases where all QME[i, ] > 0, and handle list-of-matrix (correlated m.e.)
    }

    if (!in_block) {
      # 2. draw coefficients p[[name]]
      if (all(Q0b0 == 0))
        draw <- add(draw, bquote(Xy <- crossprod_mv(p[[.(name_X)]], e$Q_e(p))))
      else
        draw <- add(draw, bquote(Xy <- crossprod_mv(p[[.(name_X)]], e$Q_e(p)) + Q0b0))

      # TODO modeled.Q case: use p[["Q_"]]
      #if (e$modeled.Q) {  # in this case both XX and Xy must be updated in each iteration
      if (e$Q0.type == "symm") {
        # TODO this case not (yet) supported!
        draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], p[["QM_"]])))
      } else {
        #draw <- add(draw, quote(XX <- crossprod_sym(p[[name_X]], p[["Q_"]])))
        if (informative.prior) {
          # TODO matsum template? or can we always assume dense matrices?
          draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], e[["Q0"]]) + Q0))
        } else
          draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], e[["Q0"]])))
      }
      if (informative.prior) {
        # TODO instead of runif(n, ...) account for the structure in Qmod; lambda may not vary per unit!
        # TODO here we need M1 to be dense without zeros
        mat_sum <- make_mat_sum(M0=Q0, M1=crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]])))
        MVNsampler <- create_TMVN_sampler(
          Q=mat_sum(crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]]))),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=mat_sum(XX), Xy=Xy)[[.(name)]]))
      } else {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]])),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=XX, Xy=Xy)[[.(name)]]))
      }
    }  # END if (!in_block)

    # compute full residuals/fitted values
    if (e$e.is.res)
      draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, p[[.(name_X)]], p[[.(name)]])))
    else
      draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, p[[.(name_X)]], p[[.(name)]])))
    draw <- add(draw, quote(p))

    start <- function(p) {
      if (!in_block) p <- MVNsampler$start(p)
      p[[name_X]] <- X
      p
    }
  }  # END if (!e$prior.only)

  rm(R, r, S, s, lower, upper)

  environment()
}
