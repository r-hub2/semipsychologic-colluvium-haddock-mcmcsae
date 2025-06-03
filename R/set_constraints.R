#' Set up a system of linear equality and/or inequality constraints
#'
#' Two-sided inequalities specified by \code{S2, l2, u2} are currently
#' transformed into the one-sided form \eqn{S'x >= s}, combined with
#' any directly specified constraints of this form. Some basic consistency
#' checks are carried out, notably regarding the dimensions of the inputs.
#'
#' @export
#' @param R equality constraint matrix each column of which corresponds to a constraint.
#' @param r right-hand side vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param S inequality constraint matrix each column of which corresponds to an inequality constraint.
#' @param s rhs vector for inequality constraints \eqn{S'x >= s}, where \eqn{S'} denotes the transpose of S.
#' @param S2 inequality constraint matrix each column of which corresponds to a two-sided inequality constraint.
#' @param l2 vector of lower bounds for two-sided inequality constraints \eqn{l_2 <= S_2'x <= u_2}.
#' @param u2 vector of upper bounds for two-sided inequality constraints \eqn{l_2 <= S_2'x <= u_2}.
#' @param scale whether to scale the columns of all constraint matrices to unit Euclidean norm.
#' @returns An environment with constraint matrices and vectors and a method to check whether a numeric vector
#'  satisfies all constraints. Returns \code{NULL} in case of no constraints.
set_constraints <- function(R=NULL, r=NULL,
                            S=NULL, s=NULL, S2=NULL, l2=NULL, u2=NULL,
                            scale=FALSE) {

  if (length(R) == 0L) {
    rm(R, r)
    n <- NULL
    eq <- FALSE
  } else {
    Rnames <- dimnames(R)[[2L]]
    R <- economizeMatrix(R, vec.as.diag=FALSE, check=TRUE)
    n <- nrow(R)
    ncR <- ncol(R)
    if (ncR > n) stop("number of columns of equality constraint matrix must not exceed its number of rows")
    if (!is.null(r)) {
      if (length(r) != ncR) stop("length of 'r' should equal the number of columns of 'R'")
      if (allv(r, 0)) r <- NULL  # NULL r indicates zero-vector right-hand side
    }
    # TODO check that R has full column rank
    if (scale) {
      f <- 1 / sqrt(colSums(R * R))
      # TODO fast scaling for supported matrix types
      R <- economizeMatrix(R %*% Cdiag(f))
      if (!is.null(r)) r <- f * r
      rm(f)
    }
    check_equalities <- function(x, tol=sqrt(.Machine$double.eps)) {
      if (!is_numeric_scalar(tol) || tol < 0) stop("'tol' must be a single nonnegative number")
      if (length(x) != n) stop("wrong length of input vector")
      dif <- if (is.null(r)) crossprod_mv(R, x) else crossprod_mv(R, x) - r
      viol <- which(abs(dif) > tol)
      if (length(viol)) {
        out <- data.frame(
          constraint_nr = viol, violation = dif[viol]
        )
        rownames(out) <- Rnames[viol]
      } else {
        out <- NULL
      }
      out
    }
    eq <- TRUE
  }

  if (length(S) == 0L) {
    rm(S, s)
    ineq1 <- FALSE
  } else {
    Snames <- dimnames(S)[[2L]]
    S <- economizeMatrix(S, vec.as.diag=FALSE, check=TRUE)
    if (is.null(n)) {
      n <- nrow(S)
    } else {
      if (nrow(S) != n) stop("incompatible constraint matrix 'S'")
    }
    ncS <- ncol(S)
    if (!is.null(s)) {
      s <- as.numeric(s)
      if (length(s) != ncS) stop("length of 's' should equal the number of columns of 'S'")
      if (allv(s, 0)) s <- NULL  # NULL s indicates zero-vector right-hand side
    }
    if (scale) {
      f <- 1 / sqrt(colSums(S * S))
      S <- economizeMatrix(S %*% Cdiag(f))
      if (!is.null(s)) s <- f * s
      rm(f)
    }
    # rename to S1 etc for use in check_inequalities; S etc may be updated below if two-sided inequalities are explicitly defined
    S1 <- S; S1names <- Snames; s1 <- s
    check_inequalities <- function(x, tol=sqrt(.Machine$double.eps)) {
      if (!is_numeric_scalar(tol) || tol < 0) stop("'tol' must be a single nonnegative number")
      if (length(x) != n) stop("wrong length of input vector")
      dif <- if (is.null(s1)) crossprod_mv(S1, x) else crossprod_mv(S1, x) - s1
      viol <- which(dif < -tol)
      if (length(viol)) {
        out <- data.frame(
          constraint_nr = viol, violation = dif[viol]
        )
        rownames(out) <- S1names[viol]
      } else {
        out <- NULL
      }
      out
    }
    ineq1 <- TRUE
  }

  if (is.null(S2) || length(S2) == 0L) {
    rm(S2, l2, u2)
    ineq2 <- FALSE
  } else {
    S2names <- dimnames(S2)[[2L]]
    S2 <- economizeMatrix(S2, vec.as.diag=FALSE, check=TRUE)
    if (is.null(n)) {
      n <- nrow(S2)
    } else {
      if (nrow(S2) != n) stop("incompatible constraint matrix 'S2'")
    }
    ncS2 <- ncol(S2)
    if (is.null(l2) || is.null(u2)) stop("l2 and u2 must be provided for inequalities corresponding to S2")
    if (length(l2) != length(u2) || any(l2 >= u2)) stop("'l2' and 'u2' incompatible")
    if (length(l2) != ncS2) stop("'S2' incompatible with 'l2' and 'u2'")
    if (scale) {
      f <- 1 / sqrt(colSums(S2 * S2))
      S2 <- economizeMatrix(S2 %*% Cdiag(f))
      l2 <- f * l2
      u2 <- f * u2
      rm(f)
    }
    check_inequalities2 <- function(x, tol=sqrt(.Machine$double.eps)) {
      if (!is_numeric_scalar(tol) || tol < 0) stop("'tol' must be a single nonnegative number")
      if (length(x) != n) stop("wrong length of input vector")
      temp <- crossprod_mv(S2, x)
      dif.l <- temp - l2
      dif.u <- temp - u2
      viol.l <- which(dif.l < -tol)
      viol.u <- which(dif.u > tol)
      out <- NULL
      if (length(viol.l) || length(viol.u)) {
        if (length(viol.l)) {
          out <- data.frame(
            constraint_nr = viol.l, violation = dif.l[viol.l]
          )
          rownames(out) <- S2names[viol.l]
        }
        if (length(viol.u)) {
          out <- rbind(out, data.frame(
            constraint_nr = viol.u, violation = dif.u[viol.u]
          ))
        }
      }
      out
    }
    ineq2 <- TRUE
    # for now use S'x >= s form
    if (ineq1) {
      S <- cbind(S, S2, -S2)
      s <- c(if (is.null(s)) numeric(ncS) else s, l2, -u2)
      ncS <- ncS + ncS2
    } else {
      S <- cbind(S2, -S2)
      s <- c(l2, -u2)
      ncS <- ncS2
    }
  }

  ineq <- ineq1 || ineq2

  check <- function(x) {
    out <- if (eq) {
      temp <- check_equalities(x)
      if (!is.null(temp)) out <- cbind(type="eq", temp)
    } else {
      out <- NULL
    }
    if (ineq1) {
      temp <- check_inequalities(x)
      if (!is.null(temp)) out <- rbind(out, cbind(type="ineq", temp))
    }
    if (ineq2) {
      temp <- check_inequalities2(x)
      if (!is.null(temp)) out <- rbind(out, cbind(type="ineq2", temp))
    }
    out
  }

  if (eq || ineq) environment() else NULL
}
