#' Create a model component object for a regression component in the variance function of a
#' gaussian sampling distribution
#'
#' This function is intended to be used on the right hand side of the \code{formula.V} argument to
#' \code{\link{create_sampler}} or \code{\link{generate_data}}.
#'
#' @export
#' @param formula a formula for the regression effects explaining the log-variance.
#'  Variable names are looked up in the data frame passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param remove.redundant whether redundant columns should be removed from the design matrix.
#'  Default is \code{FALSE}.
#' @param sparse whether the model matrix associated with \code{formula} should be sparse.
#'  The default is determined by a simple heuristic based on storage size.
#' @param X a (possibly sparse) design matrix can be specified directly, as an alternative to the
#'  creation of one based on \code{formula}. If \code{X} is specified \code{formula} is ignored.
#' @param prior prior specification for the coefficients. A normal prior can be
#'  specified using function \code{\link{pr_normal}}. Alternatively, fixed values for
#'  the coefficients can be specified using \code{\link{pr_fixed}}, e.g. to generate
#'  data with given coefficients.
#' @param Q0 prior precision matrix for the regression effects. The default is a
#'  zero matrix corresponding to a noninformative improper prior.
#'  DEPRECATED, please use argument \code{prior} instead, i.e.
#'  \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param b0 prior mean for the regression effect. Defaults to a zero vector.
#'  DEPRECATED, please use argument \code{prior} instead, i.e.
#'  \code{prior = pr_normal(mean = b0.value, precision = Q0.value)}.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'vreg'
#'  with the number of the variance model term attached.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component,
#'  intended for internal use by other package functions.
#' @references
#'  E. Cepeda and D. Gamerman (2000).
#'    Bayesian modeling of variance heterogeneity in normal regression models.
#'    Brazilian Journal of Probability and Statistics, 207-221.
#'
#'  T.I. Lin and W.L. Wang (2011).
#'    Bayesian inference in joint modelling of location and scale parameters
#'    of the t distribution for longitudinal data.
#'    Journal of Statistical Planning and Inference 141(4), 1543-1553.
vreg <- function(formula=NULL, remove.redundant=FALSE, sparse=NULL, X=NULL,
                 prior=NULL, Q0=NULL, b0=NULL, name="") {
  stop("function 'vreg' should only be used inside a formula")
}

mc_vreg <- function(formula=NULL, remove.redundant=FALSE, sparse=NULL, X=NULL,
                    prior=NULL, Q0=NULL, b0=NULL, name="",
                    sc,  # unused
                    fam,
                    in.block,  # unused
                    prior.only, data) {
  type <- "vreg"
  if (name == "") stop("missing model component name")
  store.default <- name

  if (fam[["Q0.type"]] == "symm") stop("TBI: vreg component with (compatible) non-diagonal sampling variance matrix")

  if (is.null(X)) {
    X <- model_matrix(formula, data, sparse=sparse)
    if (remove.redundant) X <- remove_redundancy(X)
  } else {
    if (is.null(dimnames(X)[[2L]])) colnames(X) <- seq_len(ncol(X))
  }
  if (nrow(X) != fam[["n"]]) stop("design matrix with incompatible number of rows")
  coef.names <- dimnames(X)[[2L]]
  label_funs <- list()
  label_funs[[name]] <- function() coef.names
  X <- economizeMatrix(X, sparse=sparse, strip.names=FALSE, check=TRUE)
  q <- ncol(X)

  if (!is.null(b0) || !is.null(Q0)) {
    warn("arguments 'b0' and 'Q0' are deprecated; please use argument 'prior' instead")
    if (is.null(prior))
      prior <- pr_normal(mean = if (is.null(b0)) 0 else b0, precision = if (is.null(Q0)) 0 else Q0)
  } else {
    if (is.null(prior))
      prior <- pr_normal(mean=0, precision=0)
  }
  switch(prior[["type"]],
    fixed = prior$init(q),
    normal = {
      prior$init(q, coef.names)
      informative.prior <- prior[["informative"]]
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
    },
    stop("'vreg' priors must be specified using one of pr_normal and pr_fixed functions")
  )
  rprior <- function(p) {
    p[[name]] <- prior$rprior()
    p
  }

  if (is_ind_matrix(X) && q < fam[["n"]])
    compute_Qfactor <- function(p) X %m*v% exp(-p[[name]])
  else
    compute_Qfactor <- function(p) exp(X %m*v% (-p[[name]]))

  # function that creates an oos prediction function (closure) based on new data
  make_predict_Vfactor <- function(newdata) {
    Xnew <- model_matrix(formula, newdata)
    # TODO remove columns not corresponding to columns in X
    Xnew <- unname(Xnew)
    rm(newdata)
    if (is_ind_matrix(Xnew) && q < nrow(Xnew))
      function(p) Xnew %m*v% exp(p[[name]])
    else
      function(p) exp(Xnew %m*v% p[[name]])
  }

  if (prior.only) {
    rm(sc, fam)
    return(environment())
  }

  if (prior[["type"]] == "fixed") {
    start <- draw <- prior[["rprior"]]
    rm(sc, fam)
    return(environment())
  }

  sumX <- 0.5 * colSums(X) - Q0b0
  MVNsampler <- create_TMVN_sampler(Q=0.5*crossprod(X) + Q0, name=name)

  proposal_scale <- 2.4/sqrt(q)
  sigma.fixed <- fam[["sigma.fixed"]]
  if (!sigma.fixed) sd.name <- fam[["sd.name"]]
  draw <- function(p) {
    # TODO consider case Q0.type="symm"
    # random walk MH proposal; see Lin and Wang (2011)
    gamma <- p[[name]]
    gamma.star <- gamma + MVNsampler$draw(p, proposal_scale)[[name]]
    if (sigma.fixed)
      Qres <- p[["Q_"]] * p[["e_"]]^2
    else
      Qres <- p[["Q_"]] * (1 / p[[sd.name]]^2) * p[["e_"]]^2
    Qratio <- exp(X %m*v% (gamma - gamma.star))
    log.ar <- 0.5 * (dotprodC(gamma, Q0 %m*v% gamma) - dotprodC(gamma.star, Q0 %m*v% gamma.star)) +
              dotprodC(gamma - gamma.star, sumX) + 0.5 * sum((1 - Qratio) * Qres)
    if (log(runif(1L)) < log.ar) {
      p[[name]] <- gamma.star
      p$Q_ <- p[["Q_"]] * Qratio
    }
    p
  }

  start <- function(p) {
    p[[name]] <- check_and_get(p, name, q, \() runif(q, -1e-3, 1e-3))
    p
  }

  rm(sc)
  environment()
}
