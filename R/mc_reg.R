#' Create a model component object for a regression (fixed effects) component
#' in the linear predictor
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates an additive regression term in the
#' model's linear predictor. By default, the prior for the regression
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
#' data(iris)
#' # default: flat priors on regression coefficients
#' sampler <- create_sampler(Sepal.Length ~
#'     reg(~ Petal.Length + Species, name="beta"),
#'   data=iris
#' )
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' # (weakly) informative normal priors on regression coefficients
#' sampler <- create_sampler(Sepal.Length ~
#'     reg(~ Petal.Length + Species, prior=pr_normal(precision=1e-2), name="beta"),
#'   data=iris
#' )
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' # binary regression
#' sampler <- create_sampler(Species == "setosa" ~
#'     reg(~ Sepal.Length, prior=pr_normal(precision=0.1), name="beta"),
#'   family="binomial", data=iris)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' pred <- predict(sim)
#' str(pred)
#' # example with equality constrained regression effects
#' n <- 500
#' df <- data.frame(x=runif(n))
#' df$y <- rnorm(n, 1 + 2*df$x)
#' R <- matrix(1, 2, 1)
#' r <- 3
#' C <- set_constraints(R=R, r=r)
#' sampler <- create_sampler(y ~ reg(~ 1 + x, constraints=C, name="beta"), data=df)
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' plot(sim, "beta")
#' summary(transform_dc(sim$beta, fun=function(x) crossprod_mv(R, x) - r))
#' }
#'
#' @export
#' @param formula a formula specifying the predictors to be used in the model,
#'  in the same way as the right hand side of the \code{formula} argument of
#'  R's \code{lm} function. Variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param remove.redundant whether redundant columns should be removed from the
#'  design matrix. Default is \code{FALSE}. But note that treatment contrasts are
#'  automatically applied to all factor variables in \code{formula}.
# TODO allow to pass contrasts.arg to model_matrix (e.g. "contr.none", and later "contr.sparse")
#' @param sparse whether the model matrix associated with \code{formula} should
#'  be sparse. The default is to base this on a simple heuristic.
#' @param X a (possibly sparse) design matrix can be specified directly, as an
#'  alternative to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
#' @param prior prior specification for the regression coefficients. Supported
#'  priors can be specified using functions \code{\link{pr_normal}},
#'  \code{\link{pr_fixed}}, or \code{\link{pr_MLiG}}. The latter prior is only
#'  available in conjunction with a gamma family sampling distribution.
#' @param Q0 prior precision matrix for the regression effects. The default is a
#'  zero matrix corresponding to a noninformative improper prior.
#'  It can be specified as a scalar value, as a numeric vector of appropriate
#'  length, or as a matrix object. DEPRECATED, please use argument \code{prior}
#'  instead, i.e. \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param b0 prior mean for the regression effect. Defaults to a zero vector.
#'  It can be specified as a scalar value or as a numeric vector of
#'  appropriate length. DEPRECATED, please use argument \code{prior}
#'  instead, i.e. \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param constraints optional linear equality and/or inequality constraints
#'  imposed on the vector of regression coefficients. Use function
#'  \code{\link{set_constraints}} to specify the constraint matrices and
#'  right-hand sides.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'reg'
#'  with the number of the model term attached.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
reg <- function(formula = ~ 1, remove.redundant=FALSE, sparse=NULL, X=NULL,
                prior=NULL, Q0=NULL, b0=NULL,
                constraints=NULL,
                name="", debug=FALSE) {

  e <- sys.frame(-2L)
  if (is.null(e[["family"]]) || is.null(e$family[["family"]])) stop("'reg()' should only be used within a formula")
  type <- "reg"
  if (name == "") stop("missing model component name")

  if (!inherits(formula, "formula")) stop("formula argument of reg() must be a formula")

  if (e$family[["family"]] == "gamma") {
    modus <- "gamma"       # model for log(mean) of gamma
  } else if (any(name == names(e[["Vmod"]]))) {
    if (e$family[["family"]] == "gaussian_gamma")
      modus <- "vargamma"  # model for log(var) of gaussian and log(mean) of gamma
    else
      modus <- "var"       # model for log(var) of gaussian
  } else {
    modus <- "regular"     # model for mean of gaussian/binomial/...
  }

  if (e$family[["family"]] == "multinomial") {
    edat <- new.env(parent = environment(formula))
    environment(formula) <- edat
  }

  if (is.null(X)) {
    if (e$family[["family"]] == "multinomial") {
      for (k in seq_len(e[["Km1"]])) {
        edat$cat_ <- factor(rep.int(e$cats[k], e[["n0"]]), levels=e$cats[-length(e$cats)])
        X <- rbind(X, model_matrix(formula, data=e[["data"]], sparse=sparse))
      }
      edat$cat_ <- NULL
    } else {
      X <- model_matrix(formula, e[["data"]], sparse=sparse)
    }
    if (remove.redundant) X <- remove_redundancy(X)
  } else {
    if (is.null(dimnames(X)[[2L]])) colnames(X) <- seq_len(ncol(X))
  }
  if (nrow(X) != e[["n"]]) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- dimnames(X)[[2L]]
  X <- economizeMatrix(X, sparse=sparse, strip.names=TRUE, check=TRUE)
  q <- ncol(X)
  in_block <- any(name == unlst(e[["block"]])) || any(name == unlst(e[["block.V"]]))

  if (!is.null(b0) || !is.null(Q0)) {
    warn("arguments 'b0' and 'Q0' are deprecated; please use argument 'prior' instead")
    if (is.null(prior)) {
      if (modus == "regular")
        prior <- pr_normal(mean = if (is.null(b0)) 0 else b0, precision = if (is.null(Q0)) 0 else Q0)
      else
        prior <- pr_MLiG(mean = if (is.null(b0)) 0 else b0, precision = if (is.null(Q0)) 0 else Q0)
    }
  } else if (is.null(prior)) {
    if (modus == "regular")
      prior <- pr_normal(mean=0, precision=0)
    else
      prior <- pr_MLiG(mean=0, precision=0)
  } else if (is.numeric(prior)) {
    prior <- pr_fixed(value=prior)
  } else {
    if (!is.environment(prior) || !is_character_scalar(prior[["type"]]))
      stop("unexpected prior input")
  }

  switch(prior[["type"]],
    fixed = {
      prior$init(q)
      informative.prior <- TRUE
      zero.mean <- allv(prior[["value"]], 0)
      rprior <- function(p) prior$rprior()
    },
    normal = {
      prior$init(q, e$coef.names[[name]], sparse=if (in_block) TRUE else NULL, sigma=!e[["sigma.fixed"]])
      informative.prior <- prior[["informative"]]
      if (modus != "regular" && informative.prior) stop("please use 'pr_MLiG' to specify a (conjugate) prior for gamma family or variance modelling")
      Q0 <- prior[["precision"]]
      zero.mean <- !informative.prior || allv(prior[["mean"]], 0)
      if (zero.mean) {
        prior$mean <- 0
        Q0b0 <- 0
      } else {
        if (length(prior[["mean"]]) == 1L)
          Q0b0 <- Q0 %m*v% rep.int(prior[["mean"]], q)
        else
          Q0b0 <- Q0 %m*v% prior[["mean"]]
      }
      rprior <- function(p) prior$rprior(p)
    },
    MLiG = {
      if (modus == "regular") stop("MLiG prior only available in combination with gamma family or variance modelling")
      if (!is.null(constraints)) stop("constraints not supported in combination with MLiG prior")
      prior$init(q, e$coef.names[[name]])
      informative.prior <- prior[["informative"]]
      zero.mean <- allv(prior[["mean"]], 0)
      rprior <- function(p) prior$rprior()
      if (in_block) {  # block sampler sparse Q0 for template
        if (length(prior[["precision"]]) == 1L)
          Q0 <- Cdiag(rep.int(prior[["precision"]]/prior[["a"]], q))
        else
          Q0 <- Cdiag(prior[["precision"]]/prior[["a"]])
      }
    },
    stop("'", prior[["type"]], "' is not a supported prior for regression coefficients '", name, "'")
  )

  is.proper <- informative.prior

  if (e[["compute.weights"]]) {
    # TODO see if following restriction can be relaxed
    if (prior[["type"]] != "normal") stop("weights computation not supported for priors of type ", prior$type)
    if (!zero.mean) stop("weights cannot be computed if some coefficients have non-zero prior means")
  }

  if (modus == "var" || modus == "vargamma") {
    if (is_ind_matrix(X) && q < e[["n"]])
      compute_Qfactor <- function(p) X %m*v% exp(-p[[name]])
    else
      compute_Qfactor <- function(p) exp(X %m*v% (-p[[name]]))
  }
  lp <- function(p) X %m*v% p[[name]]
  lp_update <- function(x, plus=TRUE, p) mv_update(x, plus, X, p[[name]])
  # draws_linpred method acts on (subset of) mcdraws object, used in fitted() and pointwise log-likelihood llh_i functions
  draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
    if (is.null(units)) Xi <- X else Xi <- X[units, , drop=FALSE]
    if (matrix) {
      out <- tcrossprod(as.matrix.dc(get_from(obj[[name]], chains=chains, draws=draws), colnames=FALSE), Xi)
    } else {
      out <- list()
      for (ch in seq_along(chains))
        out[[ch]] <- tcrossprod(get_from(obj[[name]], chains=chains[ch], draws=draws)[[1L]], Xi)
    }
    out
  }

  # for regression components assume that categorical variables in newdata
  #   carry the exact same levels as those in the training data
  make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
    if (modus == "vargamma") stop("prediction for 'gaussian_gamma' family not supported")
    if (is.null(newdata) && is.null(Xnew)) {
      # in-sample prediction
      Xnew <- X
    } else {
      if (is.null(newdata)) {
        if (!is_a_matrix(Xnew)) stop("not a matrix object")
        if (is.null(dimnames(Xnew)[[2L]]) && ncol(Xnew) != q)
          stop("wrong number of columns for Xnew matrix of component ", name)
      } else {
        nnew <- nrow(newdata)
        if (e$family[["family"]] == "multinomial") {
          Xnew <- NULL
          for (k in seq_len(e[["Km1"]])) {
            edat$cat_ <- factor(rep.int(e$cats[k], nnew), levels=e$cats[-length(e$cats)])
            Xnew <- rbind(Xnew, model_matrix(formula, data=newdata, sparse=sparse))
          }
          edat$cat_ <- NULL
        } else {
          Xnew <- model_matrix(formula, newdata, sparse=sparse)
        }
      }
      if (is.null(dimnames(Xnew)[[2L]])) {
        if (ncol(Xnew) != q) stop("'newdata' yields ", ncol(Xnew), " predictor column(s) for model term '", name, "' versus ", q, " originally")
      } else if (!identical(dimnames(Xnew)[[2L]], e$coef.names[[name]])) {
        if (!all(e$coef.names[[name]] %in% dimnames(Xnew)[[2L]])) {
          stop("the following prediction matrix columns are missing: ",
            paste0(setdiff(e$coef.names[[name]], colnames(Xnew)), collapse=", "))
        }
        if (!remove.redundant) warn("prediction matrix has additional columns")
        # TODO check that any columns removed are redundant in newdata as well
        Xnew <- Xnew[, e$coef.names[[name]], drop=FALSE]
      }
      Xnew <- economizeMatrix(Xnew, sparse=sparse, strip.names=TRUE, check=TRUE)
    }
    rm(newdata, verbose)
    linpred <- function(p) Xnew %m*v% p[[name]]
    linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, Xnew, p[[name]])
    environment()
  }

  if (!in_block && !e[["prior.only"]] && prior[["type"]] != "fixed") {

    draw <- if (debug) function(p) {browser()} else function(p) {}
    if (modus == "var" || modus == "vargamma") {
      if (!e[["single.V.block"]]) {
        if (is_ind_matrix(X) && q < e[["n"]])
          draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * (X %m*v% exp(p[[.(name)]]))))
        else
          draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(X %m*v% p[[.(name)]])))
      }
    } else {
      if (e[["single.block"]]) {
        if (e[["e.is.res"]])
          draw <- add(draw, quote(p$e_ <- e$y_eff()))
      } else {
        if (e[["e.is.res"]])
          draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, X, p[[.(name)]])))
        else
          draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, X, p[[.(name)]])))
      }
    }
    if (e[["single.block"]] && e$family[["link"]] != "probit" && modus == "regular" &&
        (!e[["modeled.Q"]] || (all(e$family[["family"]] != c("gaussian", "gaussian_gamma")) && !(e$family$family == "negbinomial" && !e$family$shape.fixed)))) {
      # single regression component, no variance modelling, Xy fixed
      if (e$family[["family"]] == "gaussian")
        Xy <- crossprod_mv(X, e[["Q0"]] %m*v% e$y_eff()) + Q0b0
      else
        Xy <- crossprod_mv(X, e$Q_e(e$y_eff())) + Q0b0
    } else if (modus == "regular") {
      # Xy updated in each iteration
      if (allv(Q0b0, 0))
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
      else
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p)) + Q0b0))
    } else {  # modus gamma, var, vargamma
      if (modus == "var" || modus == "vargamma") {  # variance modelling
        if (e[["single.V.block"]])
          draw <- add(draw, bquote(vkappa <- .(if (e[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2))
        else
          draw <- add(draw, bquote(vkappa <- .(if (e[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2 * p[["Q_"]]))
        draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), 0.5, vkappa))))
      }
      if (modus == "gamma") {
        if (e$family[["alpha.fixed"]]) {
          alpha <- e$family$get_shape()
          if (e[["single.block"]])
            kappa <- alpha * e[["y"]]
          else {
            kappa0 <- alpha * e[["y"]]
            draw <- add(draw, quote(kappa <- kappa0 * exp(-p[["e_"]])))
          }
        } else {
          draw <- add(draw, quote(alpha <- e$family$get_shape(p)))
          if (e[["single.block"]]) {
            draw <- add(draw, quote(kappa <- alpha * e[["y"]]))
          } else {
            draw <- add(draw, quote(kappa <- alpha * e[["y"]] * exp(-p[["e_"]])))
          }
        }
        draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
      } else if (modus == "vargamma") {
        if (e$family[["alpha.fixed"]]) {
          alpha <- e$family$get_shape()
          if (e[["single.V.block"]]) {
            kappa <- alpha * e$family[["sigmasq"]]
          } else {
            kappa0 <- alpha * e$family[["sigmasq"]]
            draw <- add(draw, quote(kappa <- kappa0 * p[["Q_"]]))
          }
        } else {
          draw <- add(draw, quote(alpha <- e$family$get_shape(p)))
          if (e[["single.V.block"]]) {
            draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]]))
          } else {
            draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]] * p[["Q_"]]))
          }
        }
        draw <- add(draw, bquote(Hz <- Hz + crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
      }
      if (informative.prior) {
        if (zero.mean)
          draw <- add(draw, bquote(Hz <- Hz + sqrt(prior$precision/prior$a) * rMLiG(.(q), prior$a, prior$a)))
        else
          draw <- add(draw, bquote(Hz <- Hz + sqrt(prior$precision/prior$a) * rMLiG(.(q), prior$a, prior$a * exp(sqrt(prior$precision/prior$a) * prior$mean))))
      }
    }

    if (e[["modeled.Q"]] && modus == "regular") {  # in this case both XX and Xy must be updated in each iteration
      if (e[["Q0.type"]] == "symm")
        draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
      else
        draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
      if (informative.prior) {
        # TODO instead of runif(n, ...) account for the structure in Qmod; lambda may not vary per unit!
        mat_sum <- make_mat_sum(M0=Q0, M1=crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]])))
        MVNsampler <- create_TMVN_sampler(
          Q=mat_sum(crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]]))),
          update.Q=TRUE, name=name, constraints=constraints,
          chol.control=e$control[["chol.control"]]
        )
        draw <- add(draw, bquote(p <- MVNsampler$draw(p, .(if (e[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])), Q=mat_sum(XX), Xy=Xy)))
      } else {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]])),
          update.Q=TRUE, name=name, constraints=constraints,
          chol.control=e$control[["chol.control"]]
        )
        draw <- add(draw, bquote(p <- MVNsampler$draw(p, .(if (e[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])), Q=XX, Xy=Xy)))
      }
    } else {  # precision matrix XX + Q0 not updated
      if (modus == "regular") {
        if (e[["single.block"]] && e$family[["link"]] != "probit") {
          MVNsampler <- create_TMVN_sampler(
            Q=crossprod_sym(X, e[["Q0"]]) + Q0, Xy=Xy, constraints=constraints,
            name=name, chol.control=e$control[["chol.control"]]
          )
          draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])))[[.(name)]]))
          rm(Xy)
        } else {
          MVNsampler <- create_TMVN_sampler(
            Q=crossprod_sym(X, e[["Q0"]]) + Q0,
            update.mu=TRUE, constraints=constraints,
            name=name, chol.control=e$control[["chol.control"]]
          )
          draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])), Xy=Xy)[[.(name)]]))
        }
      } else {
        if (modus == "vargamma")
          cholHH <- build_chol(economizeMatrix(2 * crossprod_sym(X, CdiagU(e[["n"]])) + (prior$precision/prior$a) * CdiagU(q), symmetric=TRUE))
        else
          cholHH <- build_chol(economizeMatrix(crossprod_sym(X, CdiagU(e[["n"]])) + (prior$precision/prior$a) * CdiagU(q), symmetric=TRUE))
        draw <- add(draw, bquote(p[[.(name)]] <- cholHH$solve(Hz)))
      }
    }
    if (modus == "var" || modus == "vargamma") {
      if (e[["single.V.block"]]) {
        if (is_ind_matrix(X) && q < e[["n"]])
          draw <- add(draw, bquote(p[["Q_"]] <- X %m*v% exp(-p[[.(name)]])))
        else
          draw <- add(draw, bquote(p[["Q_"]] <- exp(X %m*v% (-p[[.(name)]]))))
      } else {
        if (is_ind_matrix(X) && q < e[["n"]])
          draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * (X %m*v% exp(-p[[.(name)]]))))
        else
          draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(X %m*v% (-p[[.(name)]]))))
      }
    } else {
      if (e[["e.is.res"]]) {
        draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, X, p[[.(name)]])))
      } else {
        if (e[["single.block"]])
          draw <- add(draw, bquote(p$e_ <- X %m*v% p[[.(name)]]))
        else
          draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, X, p[[.(name)]])))
      }
    }
    draw <- add(draw, quote(p))

    if (modus == "regular") {
      start <- function(p) MVNsampler$start(p)
    } else {
      start <- function(p) {
        if (is.null(p[[name]])) {
          # account for the scale of covariates to avoid over/underflow of exp(linpred)
          p[[name]] <- rnorm(q) / colwise_maxabs(X)
        } else {
          p[[name]] <- as.numeric(p[[name]])
          if (length(p[[name]]) != q) stop("start value for '", name, "' has wrong length")
        }
        p
      }
    }
  } else if (!e[["prior.only"]] && prior[["type"]] == "fixed") {
    if (in_block) stop("'pr_fixed' not supported for components in a Gibbs block")
    draw <- function(p) {
      if (debug) browser()
      p[[name]] <- prior$rprior()
      p
    }
    start <- function(p) {
      if (is.null(p[[name]])) {
        p[[name]] <- prior$rprior()
      } else {
        p[[name]] <- as.numeric(p[[name]])
        if (length(p[[name]]) != q) stop("start value for '", name, "' has wrong length")
      }
      if (e[["single.block"]]) {
        # compute residuals here once and for all
        p[["e_"]] <- e$compute_e(p)
      }
      p
    }
  }

  if (in_block && (!is.null(e$control[["CG"]]) || e$control[["cMVN.sampler"]])) {
    if (informative.prior) {
      cholQV <- build_chol(Q0, control=e$control[["chol.control"]])
      drawMVNvarQ <- function(p) cholQV$Ltimes(Crnorm(q), transpose=FALSE)
    }
  }

  environment()
}

# # experimental: display the reg component's prior
# show_reg <- function(mc) {
#   mod_str <- paste(mc$name, "~")
#   mod_str <- paste(mod_str,
#     if (mc$informative.prior) {
#       if (mc$e$sigma.fixed) {
#         "N(b0, Q0^-1)"
#       } else {
#         "N(b0, sigma_^2 Q0^-1)"
#       }
#     } else {
#       "1"
#     }
#   )
#   if (!is.null(mc$MVNsampler)) {
#     # TODO MVNsampler currently only used for posterior sampling
#     if (mc$MVNsampler$eq) mod_str <- paste0(mod_str, ",  R'", mc$name, " = r")
#     if (mc$MVNsampler$ineq) mod_str <- paste0(mod_str, ",  S'", mc$name, " >= s")
#   }
#   cat(mod_str, "\n")
# }
