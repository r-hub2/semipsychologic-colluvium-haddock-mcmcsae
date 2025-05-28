

# function that creates a projection function, based on conditioning by Kriging
# project function can take M/matrix argument, but then silently assumes NULL rhs
# TODO more crossprod methods --> generalize crossprod_mm to other combinations without matrix factor
create_projector <- function(cholQ=NULL, V=NULL, update.Q=FALSE, R) {
  # assumptions: V and R are constant
  if (is.null(cholQ)) {
    if (is.null(V)) stop("either 'cholQ' or 'V' must be provided")
    if (update.Q) stop("'update.Q=TRUE' currently only supported when 'cholQ' is provided")
  } else {
    if (!is.null(V)) warn("'V' is ignored if 'cholQ' is provided")
  }
  if (update.Q || prod(dim(R)) > 1e8) {
    # TODO allow chol control options to be passed
    #cholRVR <- build_chol(crossprod_sym2(R, cholQ$solve(R)))
    cholRVR <- build_chol(crossprod_sym2(cholQ$solve(R, system="L")))
    cholQ.changed <- FALSE
    if (update.Q)
      signal_cholQ_change <- function() cholQ.changed <<- TRUE
    if (prod(dim(R)) > 1e5 && ncol(R) > 5L) {
      # this seems faster for large R
      project <- function(x, cholQ, r=NULL) {
        if (cholQ.changed) {
          cholRVR$update(crossprod_sym2(cholQ$solve(R, system="L")))
          cholQ.changed <<- FALSE
        }
        if (is.vector(x)) {
          if (is.null(r))
            x - cholQ$solve(R %m*v% cholRVR$solve(crossprod_mv(R, x)))
          else
            x + cholQ$solve(R %m*v% cholRVR$solve(r - crossprod_mv(R, x)))
        } else {
          x - cholQ$solve(R %*% cholRVR$solve(crossprod(R, x)))
        }
      }
      # alternative (not faster, unless R is dense, as cholV.R is usually less sparse than R):
      # x + cholQ$solve(cholV.R %m*v% cholRVR$solve(r - crossprod_mv(R, x)), system="Lt")
    } else {
      VR <- cholQ$solve(R)
      project <- function(x, cholQ, r=NULL) {
        if (cholQ.changed) {
          VR <<- cholQ$solve(R)
          cholRVR$update(crossprod_sym2(R, VR))
          cholQ.changed <<- FALSE
        }
        if (is.vector(x)) {
          if (is.null(r))
            x - VR %m*v% cholRVR$solve(crossprod_mv(R, x))
          else
            x + VR %m*v% cholRVR$solve(r - crossprod_mv(R, x))
        } else {
          x - VR %*% cholRVR$solve(crossprod(R, x))
        }
      }
    }
  } else {
    # NB for large problems VR.RVRinv (often dense) might become too large
    VR.RVRinv <- local({
      VR <- if (is.null(V)) cholQ$solve(R) else economizeMatrix(V %*% R, allow.tabMatrix=FALSE)
      economizeMatrix(VR %*% solve(crossprod_sym2(R, VR)), drop.zeros=TRUE)
    })
    project <- function(x, cholQ, r=NULL) {
      if (is.vector(x)) {
        if (is.null(r))
          x - VR.RVRinv %m*v% crossprod_mv(R, x)
        else
          x + VR.RVRinv %m*v% (r - crossprod_mv(R, x))
      } else {
        x - VR.RVRinv %*% crossprod(R, x)
      }
    }
  }
  environment()
}


eliminate_equalities <- function(R, r=NULL, Q, mu=NULL, keep.Q=FALSE, chol.control=chol_control()) {

  QRofR <- qr(economizeMatrix(R, sparse=FALSE))  # check the rank
  # transformation z = Q'x where Q is the orthogonal Q matrix of QRofR
  z1 <- if (is.null(r)) NULL else solve(t(qr.R(QRofR, complete=FALSE)), r)  # first, fixed part of z
  QofR <- economizeMatrix(qr.Q(QRofR, complete=TRUE))
  n1 <- ncol(R)  # nr of equality constraints, n1 >= 1
  n2 <- nrow(R) - n1
  if (n2 < 1L) stop("degenerate case: as many equality restrictions as variables")
  Q1 <- QofR[, seq_len(n1), drop=FALSE]
  Q2 <- QofR[, n1 + seq_len(n2), drop=FALSE]
  zero.x0 <- is.null(mu) || allv(mu, 0)
  if (!zero.x0) {
    mu_z1 <- crossprod_mv(Q1, mu)
    mu_z2 <- crossprod_mv(Q2, mu)
  }
  rm(R, r, QofR, mu)
  # NB cholQ will now refer to z2 subspace
  cholQ <- build_chol(crossprod_sym(Q2, Q), control=chol.control)  # chol factor L
  if (!zero.x0) {
    if (is.null(z1)) {
      x0 <- Q2 %m*v% (mu_z2 + cholQ$solve(crossprod_mv(Q2, (Q %m*v% (Q1 %m*v% mu_z1)))))
    } else {
      # fixed part of x: backtransform mu_z2_given_z1
      Q1 %m*v% z1 + Q2 %m*v% (mu_z2 - cholQ$solve(crossprod_mv(Q2, (Q %m*v% (Q1 %m*v% (z1 - mu_z1))))))
    }
    rm(mu_z1, mu_z2)
    zero.x0 <- allv(x0, 0)
    if (zero.x0) rm(x0)
  }
  if (keep.Q)
    Q <- crossprod_sym(Q2, Q)
  else
    rm(Q)
  n <- n2
  rm(Q1, z1, n1, n2)

  # now we have
  #   x0 as offset
  #   Q2 to backtransform
  #   chol_Q, e.g. to draw MVN variates, and optionally Q
  #   size n of the reduced frame

  if (zero.x0) {
    transform <- function(x) crossprod_mv(Q2, x)
    untransform <- function(z) Q2 %m*v% z
  } else {
    transform <- function(x) crossprod_mv(Q2, x - x0)
    untransform <- function(z) x0 + Q2 %m*v% z
  }

  # transform restriction matrix S and rhs s to the reduced frame
  if (zero.x0) {
    transform_restriction_rhs <- function(S, s) s
  } else {
    # NB here S is the untransformed restriction matrix
    transform_restriction_rhs <- function(S, s) if (is.null(s)) -crossprod_mv(S, x0) else s - crossprod_mv(S, x0)
  }
  transform_restriction_matrix <- function(S) economizeMatrix(crossprod(Q2, S), drop.zeros=TRUE, allow.tabMatrix=FALSE)

  environment()
}
