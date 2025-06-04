#' Solve Ax=b by preconditioned conjugate gradients
#'
#' @export
#' @keywords internal
#' @param b right hand side vector.
#' @param env environment containing at least a function \code{A_times} that computes
#'  the matrix-vector product Ax for some input vector x, and a function \code{M_solve}
#'  that computes M^-1 x for some preconditioner matrix M.
#' @param x start value for the conjugate gradient algorithm.
#' @param max.it maximum number of iterations.
#' @param e total squared error stop criterion.
#' @param verbose whether progress information is shown.
#' @param ... any parameters passed to \code{A_times} and \code{M_solve}.
#' @returns The (approximated) solution to Ax=b.
#' @references
#'  M.R. Hestenes and E. Stiefel (1952).
#'    Methods of conjugate gradients for solving linear systems.
#'    Journal of Research of the National Bureau of Standards 49(6), 409-436.
# TODO
#   - matrix rhs
#   - maybe in Rcpp
CG <- function(b, env, x=0*b, max.it=length(b), e = 1e6 * length(b), verbose=FALSE, ...) {
  r <- b - env$A_times(x, ...)
  z <- env$M_solve(r, ...)
  p <- z
  k <- 0L
  repeat {
    Ap <- env$A_times(p, ...)
    rz <- dotprodC(r, z)
    alpha <-  rz / dotprodC(p, Ap)
    x <- x + alpha * p
    # if alpha*p sufficiently small exit; TODO better criterion
    #if (max(abs(alpha * p)) < 1e-6 || k > max.it) break; this may prevent alpha becoming NA (for small e)?
    r <- r - alpha * Ap
    # if r sufficiently small exit; THIS criterion does not seem to work for non-psd preconditioner
    #if (mean(abs(r)) < sqrt(.Machine$double.eps)) break
    z <- env$M_solve(r, ...)
    discr <- dotprodC(z, z)
    if (discr < e || k > max.it) {
      # criterion used in Nishimura and Suchard
      if (discr > e) warn("max.it iterations reached with discrepancy ", discr)
      break
    }
    if (verbose) cat("iteration ", k, "   |r| = ", sqrt(sum(r^2)), "   |z| = ", sqrt(sum(z^2)), "\n")
    beta <- dotprodC(r, z) / rz
    p <- z + beta * p
    k <- k + 1L
  }
  x
}

#' Set up conjugate gradient sampler
#'
# for N(Q^-1 X' Q0 y, Q^-1)
# where Q = X' Q0 X + oplus_k QA_k x QV_k is q x q and (for mc_gen) QA_k = DA_k' DA_k may be singular
# NB experimental feature for now, needs more testing
#' @keywords internal
#' @param mbs block component containing several model components.
#' @param X design matrix.
#' @param sampler sampler object as created by \code{\link{create_sampler}}.
#' @param control a list of options for the conjugate gradient algorithm that can be passed
#'   using function \code{\link{CG_control}}.
#' @returns An environment with precomputed quantities and functions for multiplication
#'   by A and by an (inverse) preconditioning matrix.
#' @references
#'  A. Nishimura and M.A. Suchard (2022).
#'    Prior-preconditioned conjugate gradient method for accelerated Gibbs sampling in
#'    "large n, large p" Bayesian sparse regression.
#'    Journal of the American Statistical Association, 1-14.
# TODO
# - allow updates of QA (priorA or GMRF extension), i.p. QA = DA' diag(w) DA where w is updated
# - allow updates of X (model component mec)
# - use X's of underlying model components
#   and for each component use X0 and XA instead of X, using mixed product relations
#   for Khatri-Rao etc.; X0 will typically be dense and XA tabMatrix
# - p >> n case
setup_CG_sampler <- function(mbs, X, sampler, control=CG_control()) {
  q <- ncol(X)
  control <- check_CG_control(control)
  if (any(b_apply(mbs, function(mc) isTRUE(mc[["strucA"]][["type"]] == "bym2")))) {
    if (control[["preconditioner"]] != "identity")
      stop("conjugate gradient sampler for models with 'bym2' component currently only allows 'identity' preconditioner")
  }
  if (is.null(control[["max.it"]])) {
    control$max.it <- q
  } else {
    if (!is_numeric_scalar(control[["max.it"]]) || control[["max.it"]] <= 0) stop("max.it must be a single positive integer")
    control$max.it <- as.integer(control[["max.it"]])
    # warn if max.it > q ?
  }
  if (is.null(control[["stop.criterion"]])) {
    control$stop.criterion <- 1e-6 * q
  } else {
    if (!is_numeric_scalar(control[["stop.criterion"]]) || control[["stop.criterion"]] <= 0) stop("stop.criterion must be a single positive integer")
  }

  if (any(sampler$family[["family"]] == c("gaussian", "gaussian_gamma"))) {
    if (sampler[["modeled.Q"]]) {
      Q_x <- switch(sampler[["Q0.type"]],
        unit=, diag = function(p, x) p[["Q_"]] * x,
        symm = function(p, x) p[["QM_"]] %m*v% x
      )
    } else {
      Q_x <- switch(sampler[["Q0.type"]],
        unit = function(p, x) x,
        diag = function(p, x) sampler[["Q0"]]@x * x,
        symm = function(p, x) sampler[["Q0"]] %m*v% x
      )
    }
  } else {
    if (sampler$family[["link"]] == "probit") {
      Q_x <- function(p, x) x
    } else {
      Q_x <- function(p, x) p[["Q_"]] * x
    }
  }

  # set up / precompute Cholesky factors
  # for gen these factors are updated in each MCMC iteration
  cholQV <- list()
  for (mc in mbs) {
    if (mc[["type"]] == "gen") {
      if (mc[["var"]] == "unstructured")
        cholQV[[mc$name]] <- build_chol(rWishart(1L, mc[["q0"]], diag(mc[["q0"]]))[,,1L])
      else
        cholQV[[mc$name]] <- build_chol(runif(mc[["q0"]], 0.5, 1.5))
    } else {
      if (mc[["type"]] != "reg") stop("TBI")
    }
  }

  # QT = oplus_k QA x QV is created in parent's draw function
  if (sampler[["sigma.fixed"]])
    A_times <- function(x, X, QT, cholQV, p)
      crossprod_mv(X, Q_x(p, X %m*v% x)) + QT %m*v% x
  else
    A_times <- function(x, X, QT, cholQV, p)
      crossprod_mv(X, Q_x(p, X %m*v% x)) + (QT %m*v% x) * p[["sigma_"]]^2

  for (mc in mbs) {
    if (mc[["type"]] == "gen") {
      # TODO
      # - remove duplicate code, see setup_priorGMRFsampler in mc_gen
      # - update cholDD in case of local shrinkage priorA
      # cholDD is only used in the GMRF preconditioners
      if (mc[["lD"]] <= mc[["l"]]) {
        if (is.null(mc[["AR1.inferred"]]))
          mc$cholDD <- build_chol(tcrossprod(mc[["DA"]]), control=chol_control(perm=FALSE))  # sometimes perm=TRUE may be more efficient
        else
          mc$cholDD <- build_chol(tcrossprod(mc[["DA.template"]][["DA0.5"]]), control=chol_control(perm=FALSE))  # sometimes perm=TRUE may be more efficient
      } else {
        # more edges than vertices, as usual in e.g. spatial models --> adjust for non-singularity
        if (is.null(mc[["AR1.inferred"]]))
          mc$cholDD <- build_chol(tcrossprod(mc[["DA"]]) + 0.1 * CdiagU(mc[["lD"]]), control=chol_control(perm=FALSE))
        else
          mc$cholDD <- build_chol(tcrossprod(mc[["DA.template"]][["DA0.5"]]) + 0.1 * CdiagU(mc[["lD"]]), control=chol_control(perm=FALSE))
      }
    } else {
      # regression usually uses non- or weakly-informative prior -->
      # use simple regression posterior variances in preconditioner, as suggested in NS paper
      # TODO more efficient computation of diagonal values only
      mc$gamma <- 2*diag(solve(crossprod_sym(mc[["X"]], sampler[["Q0"]])))
    }
  }

  switch(control[["preconditioner"]],
    GMRF =
      # multiply by D+ DA x V, the (pseudo-)inverse of preconditioner M
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- numeric(q)
        for (mc in mbs) {
          if (mc[["type"]] == "gen") {
            DA <- if (is.null(mc[["AR1.inferred"]])) mc[["DA"]] else mc[["DA.template"]]$update(p[[mc$name_AR1]])
            temp <- DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(x[mc$block.i], nrow=mc[["q0"]])))
            out[mc$block.i] <- as.numeric(t.default(crossprod_mm(DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$block.i] <- mc[["gamma"]] * x[mc$block.i]
          }
        }
        out
      },
    GMRF2 =
      # multiply by D+ D+' x QV^-1, cf. Nishimura and Suchard
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- numeric(q)
        for (mc in mbs) {
          if (mc[["type"]] == "gen") {
            DA <- if (is.null(mc[["AR1.inferred"]])) mc[["DA"]] else mc[["DA.template"]]$update(p[[mc$name_AR1]])
            temp <- DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(x[mc$block.i], nrow=mc[["q0"]])))
            temp <- mc$cholDD$solve(temp)
            out[mc$block.i] <- as.numeric(t.default(crossprod_mm(DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$block.i] <- mc[["gamma"]] * x[mc$block.i]
          }
        }
        out
      },
    GMRF3 = {
      # multiply by D+ D+' x QV^-1, cf. Nishimura and Suchard + scaling
      scale <- control[["scale"]]
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- numeric(q)
        for (mc in mbs) {
          if (mc[["type"]] == "gen") {
            DA <- if (is.null(mc[["AR1.inferred"]])) mc[["DA"]] else mc[["DA.template"]]$update(p[[mc$name_AR1]])
            temp <- DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(scale * x[mc$block.i], nrow=mc[["q0"]])))
            temp <- mc$cholDD$solve(temp)
            out[mc$block.i] <- as.numeric(t.default(crossprod_mm(DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$block.i] <- mc[["gamma"]] * x[mc$block.i]
          }
        }
        out
      }
    },
    identity =
      M_solve <- function(x, X, QT, cholQV, p) x
  )
  self <- environment()
  # X, QT passed from block's draw function
  draw <- function(p, Xy, X, QT, sampler, start=NULL) {
    # Xy is rhs, i.e. X' Qn ytilde (+ possibly prior reg term)
    if (is.null(p[["sigma_"]])) sigma <- 1 else sigma <- p[["sigma_"]]
    u <- Xy + sigma * crossprod_mv(X, sampler$drawMVNvarQ(p))
    for (mc in mbs) {
      if (mc[["type"]] == "gen")
        u[mc$block.i] <- u[mc$block.i] + sigma^2 * mc$drawMVNvarQ(p)
      else if (mc[["informative.prior"]])
        u[mc$block.i] <- u[mc$block.i] + sigma * mc$drawMVNvarQ(p)
    }
    # solve Qtot v = u using preconditioned conjugate gradients, where Qtot = X'QX + sigma^2 QT
    out <- CG(u, self, start, max.it=control[["max.it"]], e=control[["stop.criterion"]], verbose=control[["verbose"]],
              X=X, QT=QT, cholQV=cholQV, p=p)
    if (control[["verbose"]]) cat("discrepancies: ", summary(A_times(out, X, QT, cholQV, p) - u), "\n")
    out
  }  # END function draw

  self
}

#' Set options for the conjugate gradient (CG) sampler
#'
#' @export
#' @param max.it maximum number of CG iterations.
#' @param stop.criterion total squared error stop criterion for the CG algorithm.
#' @param preconditioner one of  "GMRF", "GMRF2", "GMRF3" and "identity".
#' @param scale scale parameter; only used by the "GMRF3" preconditioner.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @param verbose whether diagnostic information about the CG sampler is shown.
#' @returns A list of options used by the conjugate gradients algorithm.
CG_control <- function(max.it=NULL, stop.criterion=NULL,
                       preconditioner=c("GMRF", "GMRF2", "GMRF3", "identity"),
                       scale=1, chol.control=chol_control(),
                       verbose=FALSE) {
  preconditioner <- match.arg(preconditioner)
  list(max.it=max.it, stop.criterion=stop.criterion, preconditioner=preconditioner, scale=scale,
       chol.control=chol.control, verbose=verbose)
}

# a manually specified list of control options might not be complete or consistent --> check it
check_CG_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- CG_control()
  w <- which(!(names(control) %in% names(defaults)))
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$chol.control <- check_chol_control(control[["chol.control"]])
  control
}
