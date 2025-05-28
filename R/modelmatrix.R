
# memory-efficient version of model.matrix
# does not first compute a full model frame
# simplifications:
# - any missing value raises an error (cf. na.action=fail in model.frame)
# - only treatment constrasts supported (i.p. no polynomial contrasts for ordered factors)
# very efficient in sparse case
# optionally uses consistent naming, e.g., var1$cat1:var2$cat2
# forbid use of category separator and : in variable names to guarantee unique column names
# option to aggregate by a factor variable (term by term, i.e. in a memory-efficient way)
# another difference compared to model_matrix is that model_matrix is tolerant
# to factors with a single level: in case contrasts are applied a warning is issued
# and the corresponding term is removed; otherwise the term is kept without warning

# TODO
# - composite matrix as list of tabMatrices
# - contrasts: "contr.sparse" for removing level with most observations


#' Compute possibly sparse model matrix
#'
#' @export
#' @param formula model formula.
#' @param data data frame containing all variables used in \code{formula}.
#'  These variables should not contain missing values. An error is raised in case any of them does.
#' @param contrasts.arg specification of contrasts for factor variables. Currently supported are
#'  "contr.none" (no contrasts applied), "contr.treatment" (first level removed) and "contr.SAS" (last level removed).
#'  Alternatively, a named list specifying a single level per factor variable can be passed.
#' @param drop.unused.levels whether empty levels of individual factor variables should be removed.
#' @param sparse if \code{TRUE} a sparse matrix of class \code{dgCMatrix} is returned. This can be efficient
#'  for large datasets and a model containing categorical variables with many categories. If \code{sparse=NULL}, the default,
#'  whether a sparse or dense model matrix is returned is based on a simple heuristic.
#' @param drop0 whether to drop any remaining explicit zeros in resulting sparse matrix.
#' @param catsep separator for concatenating factor variable names and level names.
#'  By default it is the empty string, reproducing the labels of \code{model.matrix}.
#' @param by a vector by which to aggregate the result.
#' @param tabM if \code{TRUE} return a list of tabMatrix objects.
#' @returns Design matrix X, either an ordinary matrix or a sparse \code{dgCMatrix}.
model_matrix <- function(formula, data=NULL, contrasts.arg=NULL,
                         drop.unused.levels=FALSE, sparse=NULL, drop0=TRUE, catsep="",
                         by=NULL, tabM=FALSE) {
  enclos <- environment(formula)
  if (is.null(data)) {
    trms <- terms(formula)
    if (!NROW(attr(trms, "factors"))) stop("cannot determine the sample size")
    n <- NROW(eval_in(dimnames(attr(trms, "factors"))[[1L]][1L], data, enclos))
  } else if (is.integer(data) && length(data) == 1L) {
    # data is interpreted as sample size
    trms <- terms(formula)
    n <- data
  } else {
    trms <- terms(formula, data=data)
    n <- nrow(data)
  }
  if (!is.null(by)) {
    if (length(by) != n) stop("incompatible vector 'by'")
    Maggr <- aggrMatrix(by)
  }
  has_intercept <- attr(trms, "intercept") == 1L
  tmat <- terms_matrix(trms)
  if (!length(tmat)) {
    if (has_intercept) {
      if (isTRUE(sparse))
        out <- new("dgCMatrix", i=0:(n - 1L), p=c(0L, n), x=rep.int(1, n), Dim=c(n, 1L), Dimnames=list(NULL, "(Intercept)"))
      else
        out <- matrix(1, n, 1L, dimnames=list(NULL, "(Intercept)"))
      if (is.null(by))
        return(out)
      else
        return(crossprod(Maggr, out))
    } else {
      # empty design matrix
      if (isTRUE(sparse))
        return(zeroMatrix(if (is.null(by)) n else ncol(by), 0L))
      else
        return(matrix(0, nrow=if (is.null(by)) n else ncol(by), ncol=0L))
    }
  }
  tnames <- dimnames(tmat)[[2L]]
  vnames <- dimnames(tmat)[[1L]]
  # 1. analyse
  # which variables are quantitative
  qvar <- !catvars(trms, data)
  qvar <- vnames[which(qvar)]
  q <- if (has_intercept) 1L else 0L  # nr of columns
  qd <- q  # equivalent nr of dense columns, for estimation of sparseness
  if (!is.list(contrasts.arg)) {
    contrasts.arg <- match.arg(contrasts.arg, c("contr.treatment", "contr.SAS", "contr.none"))
    if (contrasts.arg == "contr.none") contrasts.arg <- NULL
  }
  contr.list <- NULL
  first <- TRUE  # keep track of first categorical variable if the model has no intercept
  terms.removed <- rep.int(FALSE, length(tnames))  # terms with a single-level factor are removed when contrasts are applied
  # loop over terms to determine nr of columns, and store levels to be removed
  for (k in seq_along(tnames)) {
    nc <- 1L  # number of cells (or columns in design matrix)
    # loop over quantitative variables, including possibly matrices
    for (v in intersect(vnames[tmat[, k] > 0L], qvar)) {
      qv <- eval_in(v, data, enclos)
      if (NROW(qv) != n) {
        if (length(qv) == 1L)
          qv <- rep.int(qv, n)
        else
          stop("variable '", v, "' has different size ", NROW(qv))
      }
      if (inherits(qv, "data.frame")) qv <- as.matrix(qv)
      if (!all(is.finite(qv))) stop("total of ", sum(!is.finite(qv)), " NA/NaN/Inf in variable ", v)
      if (is.matrix(qv)) nc <- nc * ncol(qv)
    }
    qd <- qd + nc
    rmlevs <- NULL
    # loop over the factor vars and compute number of levels accounting for removed cols
    for (v in setdiff(vnames[tmat[, k] > 0L], qvar)) {
      fac <- eval_in(v, data, enclos)
      if (length(fac) != n) stop("variable '", v, "' has different length ", length(fac))
      if (is.factor(fac)) {
        if (drop.unused.levels) fac <- fdroplevels(fac)
      } else {
        fac <- qF(fac)
      }
      if (anyNA(fac)) stop("missing(s) in variable ", v)
      levs <- attr(fac, "levels")
      ncat <- length(levs)
      # NB if tmat[v, k] == 2L a main effect is missing and we should not remove a level!
      # also, if there is no intercept, no level should be removed for the first main factor effect
      # see the documentation of terms.formula
      if (!has_intercept && first && attr(trms, "order")[k] == 1L) {
        tmat[v, k] <- 2L  # next time we can just check for value 2, meaning not to remove a level
        first <- FALSE
      }
      # store levels to remove to be passed to fac2tabM
      if (!is.null(contrasts.arg) && tmat[v, k] == 1L) {
        if (is.list(contrasts.arg)) {
          rml <- contrasts.arg[[v]]
          if (is.null(rml))
            rml <- levs[1L]
          else
            if (all(rml != levs)) stop("invalid contrasts.arg: cannot remove level '", rml, "' from factor ", v)
        } else if (contrasts.arg == "contr.treatment") {
          rml <- levs[1L]  # remove first category
        } else if (contrasts.arg == "contr.SAS") {
          rml <- levs[ncat]  # remove last category
        }
        ncat <- ncat - 1L
        if (ncat == 0L) {
          warn("model term '", tnames[k], "' is removed as contrasts are applied to the single-level factor variable '", v, "'")
          terms.removed[k] <- TRUE  # contrasts applied to single-level factor
        }
        rmlevs <- c(rmlevs, setNames(rml, v))
      }
      nc <- nc * ncat
    }
    if (!is.null(rmlevs))
      contr.list <- c(contr.list, setNames(list(rmlevs), tnames[k]))
    q <- q + nc
    # TODO also count nonzeros, i.e. the zeros in quantitative variables!
  }
  if (catsep == ":") stop("':' is not allowed as category separator in column labels")
  if (is.null(sparse)) sparse <- qd < 0.5 * q
  # 2. construct
  if (!is.null(by)) {
    n <- ncol(Maggr)
  }
  if (sparse) {
    # allocate memory for i, x slots
    i <- integer(n * (ncol(tmat) + has_intercept))
    x <- numeric(length(i))
    p <- integer(q + 1L)
  } else {
    out <- matrix(NA_real_, n, q)
  }
  lab <- character(q)
  # loop over terms and fill in i,x slots of sparse model matrix
  if (has_intercept) {
    lab[1L] <- "(Intercept)"
    xk <- if (is.null(by)) 1 else colSums(Maggr)
    if (sparse) {
      i[seq_len(n)] <- 0:(n - 1L)
      x[seq_len(n)] <- xk
      p[2L] <- n
      at <- n + 1L
    } else {
      out[, 1L] <- xk
    }
    col <- 2L
  } else {
    if (sparse) at <- 1L
    col <- 1L
  }
  for (k in seq_along(tnames)) {
    if (terms.removed[k]) next
    countvars <- intersect(vnames[tmat[, k] > 0L], qvar)
    if (length(countvars)) {
      xk <- eval_in(countvars[1L], data, enclos)
      if (inherits(xk, "data.frame")) xk <- as.matrix(xk)
      if (is.matrix(xk))
        labk <- paste(countvars[1L], if (is.null(dimnames(xk)[[2L]])) seq_len(ncol(xk)) else dimnames(xk)[[2L]], sep=catsep)
      else
        labk <- countvars[1L]
      for (v in countvars[-1L]) {
        temp <- eval_in(v, data, enclos)
        if (is.matrix(temp)) {
          xk <- t(KhatriRao(t(xk), t(temp)))
          labk <- outer(labk, paste(v, dimnames(xk)[[2L]], sep=catsep), paste, sep=":")
        } else {
          xk <- xk * temp
          labk <- paste(labk, v, sep=":")
        }
      }
    } else {
      xk <- NULL
      labk <- NULL
    }
    facvars <- setdiff(vnames[tmat[, k] > 0L], qvar)
    if (length(facvars)) {
      if (length(countvars) && !is.matrix(xk)) {
        # TODO allow matrix; --> generalize x-slot in tabMatrix to matrix (and even dgCMatrix)
        fk <- fac2tabM(facvars, data, enclos, x=xk, drop.unused.levels=drop.unused.levels, contrasts=contr.list[[tnames[k]]], catsep=catsep)
      } else {
        fk <- fac2tabM(facvars, data, enclos, drop.unused.levels=drop.unused.levels, contrasts=contr.list[[tnames[k]]], catsep=catsep)
      }
      if (is.matrix(xk)) {
        lab[col:(col + ncol(fk)*ncol(xk) - 1L)] <- outer(dimnames(fk)[[2L]], labk, paste, sep=":")
        fk <- t(KhatriRao(t(xk), t(fk)))  # col-index of fk runs fastest
      } else {
        lab[col:(col + ncol(fk) - 1L)] <- paste0(labk, if (!is.null(labk)) ":" else "", dimnames(fk)[[2L]])
      }
      if (!is.null(by)) {
        fk <- crossprod(Maggr, fk)
      }
      if (sparse) {
        if (class(fk)[1L] != "dgCMatrix") {
          if (class(fk)[1L] == "tabMatrix")
            fk <- Ctab2dgC(fk)
          else
            fk <- as(fk, "CsparseMatrix")
        }
        if (l <- length(fk@i)) {
          i[at:(at + l - 1L)] <- fk@i
          x[at:(at + l - 1L)] <- fk@x
          at <- at + l
        }
        p[(col + 1L):(col + ncol(fk))] <- p[col] + fk@p[-1L]
      } else {
        if (class(fk)[1L] == "tabMatrix")
          out[, col:(col + ncol(fk) - 1L)] <- Ctab2mat(fk)
        else
          out[, col:(col + ncol(fk) - 1L)] <- as(fk, "matrix")
      }
      col <- col + ncol(fk)
    } else {
      if (is.matrix(xk)) {
        lab[col:(col + ncol(xk) - 1L)] <- labk
        if (!is.null(by)) {
          xk <- crossprod(Maggr, xk)
        }
        if (sparse) {
          i[at:(at + length(xk) - 1L)] <- rep.int(0:(n - 1L), ncol(xk))
          x[at:(at + length(xk) - 1L)] <- xk
          p[(col + 1L):(col + ncol(xk))] <- seq.int(at + n - 1L, by=n, length.out=ncol(xk))
          at <- at + length(xk)
        } else {
          out[, col:(col + ncol(xk) - 1L)] <- xk
        }
        col <- col + ncol(xk)
      } else {
        lab[col] <- labk
        if (!is.null(by)) {
          xk <- crossprod_mv(Maggr, xk)
        }
        if (sparse) {
          i[at:(at + length(xk) - 1L)] <- 0:(n - 1L)
          x[at:(at + length(xk) - 1L)] <- xk
          at <- at + length(xk)
          p[col + 1L] <- at - 1L
        } else {
          out[, col] <- xk
        }
        col <- col + 1L
      }
    }
  }
  if (sparse) {
    i <- i[seq_len(at - 1L)]
    x <- x[seq_len(at - 1L)]
    out <- new("dgCMatrix", i=i, p=p, x=x, Dim=c(n, q), Dimnames=list(NULL, lab))
    if (drop0) out <- drop0(out, is.Csparse=TRUE)
  } else {
    dimnames(out) <- list(NULL, lab)
  }
  out
}


# return a named logical vector indicating for each variable in
#   the model terms object whether it is categorical
catvars <- function(trms, data) {
  enclos <- environment(trms)
  vnames <- dimnames(terms_matrix(trms))[[1L]]
  vapply(vnames, function(x) {
      temp <- eval_in(x, data, enclos)
      is.factor(temp) || is.character(temp)
    }, FALSE
  )
}

# extract the independent variables' terms matrix from a terms object
terms_matrix <- function(trms) {
  tmat <- attr(trms, "factor")
  if (attr(trms, "response")) {
    w <- whichv(rowSums(tmat), 0)
    if (length(w)) tmat <- tmat[-w, , drop=FALSE]
  }
  tmat
}
