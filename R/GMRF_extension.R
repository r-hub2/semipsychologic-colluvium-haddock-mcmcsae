#' Set a GMRF structure for a generic model component
#'
#' This function is used to specify a (non-default) GMRF structure
#' to pass to argument \code{strucA} of function \code{\link{gen}}.
#'
#' @export
#' @param type one of "default", "bym2" or "leroux". The default choice
#'  corresponds to the precision matrix \eqn{Q_A} as specified by argument \code{factor}
#'  of \code{gen}. Type "bym2" modifies the default structure to one with
#'  covariance matrix \eqn{\phi \tilde{Q}_{A}^- + (1 - \phi) I} where
#'  \eqn{\tilde{Q}_{A*}^-} is the generalized inverse of \eqn{Q_A}, by default
#'  scaled such that the geometric mean of the marginal variances equals 1.
#'  Type "leroux" modifies the default structure to one with precision matrix
#'  \eqn{\phi Q_A + (1 - \phi) I}.
#' @param scale.precision whether to scale the structured precision matrix. By default
#'  set to \code{TRUE} only for type "bym2".
#' @param prior prior for the parameter \eqn{phi} in the "bym2" or "leroux" extension.
#'  Supported priors can be set using functions \code{\link{pr_fixed}} or \code{\link{pr_unif}}.
#' @param control options for the Metropolis-Hastings sampler used to sample
#'  from the full conditional distribution of parameter \eqn{phi} in case of
#'  "bym2" or "leroux" extensions. If \code{NULL} a reasonable default configuration
#'  is used. A user can change these settings using function \code{\link{set_MH}}.
#'  Supported proposal distribution types are "RWTN", "RWN", "unif" and "beta".
#' @returns An environment defining the desired GMRF structure, for use by other
#'  package functions.
#' @references
#'  B. Leroux, X. Lei and N. Breslow (1999).
#'    Estimation of Disease Rates in Small Areas: A New Mixed Model for Spatial Dependence.
#'    In M. Halloran and D. Berry (Eds.), Statistical Models in Epidemiology,
#'    the Environment and Clinical Trials, 135-178.
#'
#'  A. Riebler, S.H. Sorbye, D. Simpson and H. Rue (2016).
#'    An intuitive Bayesian spatial model for disease mapping that accounts for scaling.
#'    Statistical methods in medical research, 25(4), 1145-1165.
GMRF_structure <- function(type=c("default", "bym2", "leroux"),
                           scale.precision = (type == "bym2"),
                           prior=NULL, control=NULL) {
  type <- match.arg(type)
  if (!is.logical(scale.precision) || length(scale.precision) != 1L) stop("wrong input for 'scale.precision'")
  if (type == "default") {
    prior <- NULL
    control <- NULL
    return(environment())
  }
  if (is.null(prior)) {
    prior <- pr_unif(0, 1)
  } else if (is_numeric_scalar(prior)) {
    prior <- pr_fixed(value=prior)
  } else {
    if (!is.environment(prior) || is.null(prior$type)) stop("invalid prior specification")
    # TODO allow truncnormal prior
    if (!any(prior$type == c("fixed", "beta", "unif"))) stop("unsupported prior")
  }
  if (prior$type == "fixed") {
    if (prior$value <= 0 || prior$value >= 1) stop("fixed prior value for 'bym2' or 'leroux' extension must be strictly between 0 and 1")
  }
  if (prior$type == "fixed") {
    control <- NULL
  } else if (is.null(control)) {
    control <- set_MH(type="RWTN", scale=0.025, adaptive=TRUE, l=0, u=1)
  } else {
    if (!is.environment(control))
      stop("GMRF_structure: 'control' must be an environment created with function set_MH")
    control$type <- match.arg(control$type, c("RWTN", "RWN", "unif", "beta"))
    if (any(control$type == c("RWTN", "RWN", "unif"))) {
      if (is.null(control$l)) control$l <- 0
      if (is.null(control$u)) control$u <- 1
      if (control$l > control$u) stop("lower bound of MH truncated proposal exceeds upper bound")
    }
  }
  environment()
}


create_GMRF_structure <- function(settings, mc, prior.only=FALSE) {

  type <- settings$type
  prior <- settings$prior
  control <- settings$control
  scale.precision <- settings$scale.precision
  rm(settings)
  if (type == "default") {
    update.Q <- FALSE
    return(environment())
  }

  update.Q <- prior$type != "fixed"
  if (is_unit_diag(mc$QA)) warn("possibly inefficient model structure: GMRF extension for iid random effects")

  prior$init(1L)
  rprior <- prior$rprior
  if (update.Q) name_ext <- paste0(mc[["name"]], "_GMRFext")
  if (type == "leroux") {
    if (update.Q) {
      idL <- CdiagU(mc[["l"]])
      mat_sum <- make_mat_sum(M1=mc$QA, M2=idL)
      update_Q <- function(QA, phi) mat_sum(QA, idL, phi, 1 - phi)
      if (!prior.only) {
        # base the determinant template on 0.5*(QA + I) to ensure it is non-singular
        detQ <- make_det(update_Q(mc$QA, 0.5))
        name_detQ <- paste0(mc[["name"]], "_detQ_")
        start <- function(p) {
          if (is.null(p[[name_ext]])) p[[name_ext]] <- runif(1L)
          if (is.null(p[[name_detQ]])) p[[name_detQ]] <- detQ(2*p[[name_ext]], 1 - 2*p[[name_ext]])
          p
        }
        draw <- function(p, coef.raw) {
          phi <- p[[name_ext]]
          # draw a candidate value (proposal)
          #   and initialize log acceptance probability
          phi.star <- control$propose(phi)
          phi.diff <- phi.star - phi
          QA.diff <- mat_sum(mc[["QA"]], idL, phi.diff, -phi.diff)
          tr.diff <- switch(mc[["var"]],
            unstructured = sum(crossprod_sym(coef.raw, QA.diff) * p[[mc$name_Qraw]]),
            diagonal = sum(.colSums(coef.raw * (QA.diff %m*m% coef.raw), mc[["l"]], mc[["q0"]]) * p[[mc$name_Qraw]]),
            scalar =
              if (mc[["q0"]] == 1L)
                dotprodC(coef.raw, QA.diff %m*v% coef.raw) * p[[mc$name_Qraw]]
              else
                sum(coef.raw * (QA.diff %m*m% coef.raw)) * p[[mc$name_Qraw]]
          )
          detQ.star <- detQ(2*phi.star, 1 - 2*phi.star)
          log.ar.post <- 0.5 * (mc[["q0"]] * (detQ.star - p[[name_detQ]]) - tr.diff)
          if (control$MH_accept(phi.star, phi, log.ar.post)) {
            p[[name_ext]] <- phi.star
            p[[name_detQ]] <- detQ.star
          }
          p
        }
      }
    } else {
      mc$QA <- economizeMatrix(
        prior$value * mc$QA + (1 - prior$value) * CdiagU(mc[["l"]]),
        symmetric=TRUE, sparse=if (mc$in_block) TRUE else NULL
      )
    }
  } else {  # bym2
    l <- mc[["l"]] %/% 2L  # initial value of 'l'
    if (!is.null(mc$RA)) {
      mc$RA <- economizeMatrix(
        rbind(zeroMatrix(nrow(mc$RA), ncol(mc$RA)), mc$RA),
        sparse = if (mc$in_block) TRUE else NULL, allow.tabMatrix = FALSE
      )
    }
    if (update.Q) {
      idL <- CdiagU(l)
      QA.template <- economizeMatrix(
        rbind(
          cbind((1/(1-0.5))*idL, -(sqrt(0.5)/(1-0.5))*idL),
          cbind(-(sqrt(0.5)/(1-0.5))*idL, mc$QA + (0.5/(1-0.5))*idL)
        ),
        symmetric=TRUE, sparse=if (mc$in_block) TRUE else NULL
      )
      ind1 <- seq_len(l)  # indices for 1/(1-phi) I in upper-left quadrant
      j <- get_col_ind(QA.template)
      ind2 <- which(QA.template@i == j - l)  # indices for -sqrt(phi)/(1-phi) I in upper-right quadrant
      ind3 <- which(QA.template@i == j & j >= l)  # indices for term phi/(1-phi) I in lower-right quadrant
      rm(j)
      update_Q <- function(QA, phi) {
        out <- QA.template
        out@x[ind1] <- 1/(1-phi)
        out@x[ind2] <- -sqrt(phi)/(1-phi)
        out@x[ind3] <- out@x[ind3] + (phi/(1 - phi) - 1)
        out
      }
      # template for computing QA(phi.star) - QA(phi)
      QA.diff.template <- economizeMatrix(
        rbind(
          cbind(idL, idL),
          cbind(idL, idL)
        ),
        symmetric=TRUE, sparse=if (mc$in_block) TRUE else NULL
      )
      ind.diff2 <- seq.int(l + 1L, by=2L, length.out=l)  # upper-right quadrant
      ind.diff3 <- seq.int(l + 2L, by=2L, length.out=l)  # lower-right quadrant
      QA_diff <- function(phi.star, phi) {
        out <- QA.diff.template
        out@x[ind1] <- 1/(1 - phi.star) - 1/(1 - phi)
        out@x[ind.diff2] <- - sqrt(phi.star)/(1 - phi.star) + sqrt(phi)/(1 - phi)
        out@x[ind.diff3] <- phi.star/(1 - phi.star) - phi/(1 - phi)
        out
      }
      if (!prior.only) {
        start <- function(p) {
          if (is.null(p[[name_ext]])) p[[name_ext]] <- runif(1L)
          p
        }
        draw <- function(p, coef.raw) {
          phi <- p[[name_ext]]
          phi.star <- control$propose(phi)
          QA.diff <- QA_diff(phi.star, phi)
          tr.diff <- switch(mc[["var"]],
            unstructured = sum(crossprod_sym(coef.raw, QA.diff) * p[[mc$name_Qraw]]),
            diagonal = sum(.colSums(coef.raw * (QA.diff %m*m% coef.raw), mc[["l"]], mc[["q0"]]) * p[[mc$name_Qraw]]),
            scalar =
              if (mc[["q0"]] == 1L)
                dotprodC(coef.raw, QA.diff %m*v% coef.raw) * p[[mc$name_Qraw]]
              else
                sum(coef.raw * (QA.diff %m*m% coef.raw)) * p[[mc$name_Qraw]]
          )
          log.ar.post <- -0.5 * (mc[["q0"]] * l * log((1 - phi.star)/(1 - phi)) + tr.diff)
          if (control$MH_accept(phi.star, phi, log.ar.post))
            p[[name_ext]] <- phi.star
          p
        }
      }
    } else {
      mc$QA <- economizeMatrix(
        rbind(
          cbind((1/(1-prior$value))*CdiagU(l), -(sqrt(prior$value)/(1-prior$value))*CdiagU(l)),
          cbind(-(sqrt(prior$value)/(1-prior$value))*CdiagU(l), mc$QA + (prior$value/(1-prior$value))*CdiagU(l))
        ),
        symmetric=TRUE, sparse=if (mc$in_block) TRUE else NULL
      )
    }
  }
  environment()
}
