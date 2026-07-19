
#' Set up a function for direct sampling from a constrained multivariate normal distribution
#'
#'
#' @export
#' @param D factor of precision matrix Q such that Q=D'D.
#' @param Q precision matrix.
#' @param update.Q whether to update (D and) Q for each draw.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param eps1 scalar parameter to control numerical robustness against singularity of Q.
#'  Currently ignored in case of no constraints.
#' @param eps2 scalar parameter associated with the constraint part to control numerical robustness.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @returns An environment with precomputed quantities and a method 'draw' for sampling
#'   from a multivariate normal distribution subject to equality constraints.
create_cMVN_sampler <- function(D=NULL, Q=NULL, update.Q=FALSE, R=NULL, r=NULL,
                                eps1=sqrt(.Machine$double.eps), eps2=sqrt(.Machine$double.eps),
                                chol.control=chol_control(perm=TRUE)) {
  if (is.null(Q)) {
    if (is.null(D)) stop("at least one of 'D' and 'Q' must be provided")
    Q <- crossprod(D)
  }
  q <- nrow(Q)
  if (is.null(R)) {
    # do we need LDL=TRUE here too in case of singular Q?
    cholQ <- build_chol(Q, control=chol.control)
  } else {
    # TODO create a function for handling R, as this code is duplicated from TMVN_sampler
    R <- economizeMatrix(R, vec.as.diag=FALSE, check=TRUE)
    if (nrow(R) != q || ncol(R) > q) stop("incompatible constraint matrix 'R'")
    if (is.null(r)) {
      r <- numeric(ncol(R))
    } else {
      r <- as.numeric(r)
      if (length(r) != ncol(R)) stop("length of 'r' should equal the number of columns of 'R'")
    }
    # TODO check that R has full column rank
    # parameters eps1, eps2 modify Q.expanded for more numerical stability
    # seems to improve the accuracy with which the constraints are satisfied, especially if update.Q=TRUE
    Q.expanded <- economizeMatrix(
      rbind(
        cbind(if (eps1 == 0) Q else Q + eps1 * CdiagU(q), R),
        cbind(t(R), - eps2 * CdiagU(ncol(R)))
      ), sparse=TRUE, symmetric=TRUE  # for now only support dsCMatrix
    )
    # need LDL=TRUE as the system is indefinite
    cholQ <- build_chol(Q.expanded, control=chol.control, LDL=TRUE)
    if (update.Q) {
      switch(class(Q)[1L],
        matrix = {
          if (eps1 > 0) I.diag <- seq.int(1L, by=q+1L, length.out=q)
          I.nz <- whichv(abs(as.numeric(if (eps1 == 0) Q else Q + eps1 * CdiagU(q))) > .Machine$double.eps, TRUE)
          upper.ind <- upper.tri(Q, diag=TRUE)
          update_Q.expanded <- function(Q) {
            x <- Q[upper.ind]
            if (eps1 > 0) x[I.diag] <- x[I.diag] + eps1
            attr(Q.expanded, "x")[seq_along(I.nz)] <<- x[I.nz]
          }
        },
        ddiMatrix = {
          if (is_unit_ddi(Q)) stop("unit Diagonal matrix not allowed if update.Q = TRUE")
          if (eps1 == 0) I.nz <- whichv(abs(Q@x) > .Machine$double.eps, TRUE)
          update_Q.expanded <- function(Q) {
            if (eps1 == 0)
              attr(Q.expanded, "x")[seq_along(I.nz)] <<- Q@x[I.nz]
            else
              attr(Q.expanded, "x")[seq_len(q)] <<- Q@x + eps1
          }
        },
        dsCMatrix = {
          if (eps1 > 0)
            matsum <- make_mat_sum(M0 = eps1 * CdiagU(q), M1 = Q, force.sparse=TRUE)
          update_Q.expanded <- function(Q) {
            if (eps1 > 0) {
              Qplus <- matsum(Q)
              attr(Q.expanded, "x")[seq_along(Qplus@x)] <<- Qplus@x
            } else
              attr(Q.expanded, "x")[seq_along(Q@x)] <<- Q@x
          }
        },
        stop("'Q' has unsupported type")
      )
    }
  }
  rhs <- c(numeric(q), r)
  Iq <- seq_len(q)
  if (update.Q) {
    if (is.null(R)) {
      update <- function(Q, Imult=0) cholQ$update(Q, Imult)
    } else {
      update <- function(Q, Imult=0) {
        update_Q.expanded(Q)
        cholQ$update(Q.expanded, Imult)
      }
    }
    draw <- function(u) {
      rhs[Iq] <- u
      cholQ$solve(rhs)[Iq]
    }
  } else {
    if (is.null(D)) stop("'D' must be provided when update.Q=FALSE")
    qstar <- nrow(D)
    rm(Q.expanded)
    draw <- function(u) {
      rhs[Iq] <- crossprod_mv(D, Crnorm(qstar))
      cholQ$solve(rhs)[Iq]
    }
    formals(draw) <- alist()
  }
  environment()
}


#' Set up a a function for direct sampling from a constrained multivariate normal distribution
#'
#' @keywords internal
#' @param mbs block component containing several model components.
#' @param X design matrix.
#' @param Q precision matrix.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param fam family object.
#' @param name name of the cMVN vector parameter.
#' @param control list with numerical control options, which can be set using function \code{\link{cMVN_control}}.
#' @returns An environment with precomputed quantities and functions for sampling
#'   from a multivariate normal distribution subject to equality constraints.
create_block_cMVN_sampler <- function(mbs, X, Q, R=NULL, r=NULL, fam, name="x", control) {

  if (name == "") stop("empty name")

  chol.control <- control[["chol.control"]]
  smplr <- create_cMVN_sampler(D=NULL, Q=Q, update.Q=TRUE, R=R, r=r,
    eps1=control[["eps1"]], eps2=control[["eps2"]], chol.control=chol.control
  )
  update <- function(Q, Imult=0) smplr$update(Q, Imult)

  # Xy is rhs, i.e. X' Qn ytilde (+ possibly prior reg term)
  # for multi-response family there is no global sigma_ parameter
  draw <- function(p, Xy, X) {
    if (is.null(p[["sigma_"]])) sigma <- 1 else sigma <- p[["sigma_"]]
    u <- Xy + sigma * crossprod_mv(X, fam$drawMVNvarQ(p))
    for (mc in mbs) {
      if (mc[["type"]] == "gen")
        u[mc$block.i] <- u[mc$block.i] + sigma^2 * mc$drawMVNvarQ(p)
      else if (mc[["informative.prior"]])
        u[mc$block.i] <- u[mc$block.i] + sigma * mc$drawMVNvarQ(p)
    }
    p[[name]] <- smplr$draw(u)
    p
  }

  environment()
}

#' Compute a Monte Carlo estimate of the marginal variances of a (I)GMRF
#'
#' Estimate marginal variances of a (I)GMRF prior defined in terms
#' of a sparse precision matrix and possibly a set of equality constraints.
#' The marginal variances might be used to rescale the precision matrix
#' such that a default prior for a corresponding variance component is
#' more appropriate.
#'
# @examples
# n <- 100
# mv <- sim_marg_var(D=D_RW2(n), R=R_RW2(n))
# plot(seq_len(n), mv)
# mv <- sim_marg_var(D=D_season(n, 12), R=R_season(n, 12))
# plot(seq_len(n), mv)
#
#' @export
#' @param D factor of precision matrix Q such that Q=D'D.
#' @param Q precision matrix.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param eps1 passed to \code{\link{create_cMVN_sampler}}.
#' @param eps2 passed to \code{\link{create_cMVN_sampler}}.
#' @param nSim number of Monte Carlo samples used to estimate the marginal variances.
#' @returns A vector of Monte Carlo estimates of the marginal variances.
#' @references
#'  S.H. Sorbye and H. Rue (2014).
#'    Scaling intrinsic Gaussian Markov random field priors in spatial modelling.
#'    Spatial Statistics, 8, 39-51.
sim_marg_var <- function(D, Q=NULL, R=NULL, r=NULL, eps1=1e-9, eps2=1e-9, nSim=100L) {
  smplr <- create_cMVN_sampler(D, Q, R=R, r=r, eps1=eps1, eps2=eps2)
  res <- matrix(NA_real_, ncol(D), nSim)
  for (i in seq_len(nSim)) res[, i] <- smplr$draw()
  rowVarsC(res)
}

#' Set options for fast sampling from a degenerate multivariate normal distribution
#'
#' @export
#' @param eps1 scalar parameter to control numerical robustness against singularity of Q.
#' @param eps2 scalar parameter associated with the constraint part to control numerical robustness.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @returns A list of options used by the constrained MVN sampler.
cMVN_control <- function(eps1=sqrt(.Machine$double.eps), eps2=sqrt(.Machine$double.eps),
                         chol.control=chol_control()) {
  list(eps1=eps1, eps2=eps2, chol.control=chol.control)
}

# a manually specified list of control options might not be complete or consistent --> check it
check_cMVN_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- cMVN_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$chol.control <- check_chol_control(control[["chol.control"]])
  control
}
