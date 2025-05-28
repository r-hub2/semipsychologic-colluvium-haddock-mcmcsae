# tabMatrix class: generalization of indMatrix from Matrix package
# it has at most one non-zero value in each row
setClass("tabMatrix",
  slots = c(
    perm = "integer",  # as in indMatrix, but 0-based
    reduced = "logical",  # whether columns are removed (contrasts); a -1 in perm represents such columns
    # if reduced there can be all-0 rows
    num = "logical",  # whether the matrix is numeric or 0/1
    # TODO generalise x to matrix/dgCMatrix
    x = "numeric",  # the numeric values in case num=TRUE
    xlab = "character"  # name associated with numeric x component (scalar, later to be generalized to colnames)
  ),
  contains = c("sparseMatrix", "generalMatrix"),
  validity = function(object) {
    n <- object@Dim[1L]
    d <- object@Dim[2L]
    perm <- object@perm
    num <- object@num
    reduced <- object@reduced
    xlab <- object@xlab
    if (length(perm) != n) return(paste0("length of 'perm' slot must be ", n))
    if (num && length(object@x) != n) return(paste0("length of 'x' slot must be ", n))
    if (length(xlab) > 1L) return("length of 'xlab' must not be greater than 1")
    if (n > 0L && (any(perm > d - 1L) || (!reduced && any(perm < 0L)))) return("'perm' slot is not a valid index")
    # NB 'num' can represent zero rows, so it is more general than 'reduced'
    #    so we disallow the case that both are TRUE
    #    'reduced' is used for matrices having only 0/1 entries and at least one row without 1
    if (num && reduced) return("'num' and 'reduced' cannot both be TRUE")
    TRUE
  }
)

# general constructor
tabMatrix <- function(Dim, Dimnames=NULL, reduced=FALSE, perm=integer(), num=FALSE, x=numeric(), xlab=NULL) {
  # TODO check that elements in perm are between 0 and Dim[2] - 1
  out <- Ctab(Dim=Dim, reduced=reduced, perm=perm, num=num, x=x)
  if (!is.null(Dimnames)) attr(out, "Dimnames") <- Dimnames
  if (!is.null(xlab)) attr(out, "xlab") <- xlab
  out
}

# constructor from factor
setAs("factor", "tabMatrix",
  function(from) {
    levs <- attr(from, "levels")
    tabMatrix(
      Dim = c(length(from), length(levs)),
      Dimnames = list(NULL, levs),
      reduced = FALSE,
      perm = as.integer(from) - 1L,
      num = FALSE,
      x = numeric(0L),
      xlab = character(0L)
    )
  }
)

# index matrix has exactly one 1 in each row and otherwise 0s
is_ind_matrix <- function(X) class(X)[1L] == "tabMatrix" && !X@num && !X@reduced

setMethod("colSums", "tabMatrix", \(x, na.rm=FALSE, dims=1, ...)
  if (x@num)
    fast_aggrC(x@x, x@perm + 1L, x@Dim[2L])
  else
    tabulate(x@perm + 1L, nbins=x@Dim[2L])
)

setMethod("rowSums", "tabMatrix", \(x, na.rm=FALSE, dims=1, ...) {
  if (x@num) return(x@x)
  out <- rep.int(1, x@Dim[1L])
  if (x@reduced) out[x@perm == -1L] <- 0
  out
})

setMethod("isDiagonal", "tabMatrix", \(object) {
  d <- object@Dim
  if (d[1L] != d[2L]) return(FALSE)
  if (object@reduced) {
    sub <- whichv(object@perm, -1L, invert=TRUE)
    identical(object@perm[sub], (0:(d[1L] - 1L))[sub])
  } else {
    identical(object@perm, 0:(d[1L] - 1L))
  }
})

setMethod("diag", "tabMatrix", \(x=1, nrow, ncol) {
  d <- x@Dim
  if (d[1L] <= d[2L]) {
    if (x@num)
      out <- ifelse(x@perm == 0:(d[1L] - 1L), x@x, 0)
    else
      out <- ifelse(x@perm == 0:(d[1L] - 1L), 1, 0)
  } else {
    if (x@num)
      out <- ifelse(x@perm[seq_len(d[2L])] == 0:(d[2L] - 1L), x@x[seq_len(d[2L])], 0)
    else
      out <- ifelse(x@perm[seq_len(d[2L])] == 0:(d[2L] - 1L), 1, 0)
  }
  if (x@reduced) out[x@perm == -1L] <- 0
  out
})

setMethod("t", "tabMatrix", \(x) {
  if (x@reduced) {
    sub <- whichv(x@perm, -1L, invert=TRUE)
    sparseMatrix(i=x@perm[sub], j=(0:(x@Dim[1L] - 1L))[sub], x=rep.int(1, length(sub)), dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
  } else if (x@num) {
    sparseMatrix(i=x@perm, j=0:(x@Dim[1L] - 1L), x=x@x, dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
  } else {
    sparseMatrix(i=x@perm, j=0:(x@Dim[1L] - 1L), x=rep.int(1, x@Dim[1L]), dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
  }
})

setMethod("coerce", c("tabMatrix", "ddiMatrix"), \(from, to="ddiMatrix", strict=TRUE) {
  if (!isDiagonal(from)) stop("not a diagonal tabMatrix")
  if (from@reduced) {
    x <- rep.int(1, from@Dim[1L])
    x[from@perm == -1L] <- 0
    out <- Cdiag(x)
  } else if (from@num) {
    out <- Cdiag(from@x)
  } else {
    out <- CdiagU(from@Dim[1L])
  }
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "Dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("tabMatrix", "matrix"), \(from, to="matrix", strict=TRUE) {
  out <- Ctab2mat(from)
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("tabMatrix", "CsparseMatrix"), \(from, to="CsparseMatrix", strict=TRUE) {
  out <- Ctab2dgC(from)
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "Dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("ddiMatrix", "tabMatrix"), \(from, to="tabMatrix", strict=TRUE)
  tabMatrix(Dim=from@Dim, Dimnames=from@Dimnames,
    reduced=FALSE, perm=0:(from@Dim[1L] - 1L), num = from@diag == "N",
    x = if (from@diag == "N") from@x else numeric()
  )
)

setMethod("coerce", c("matrix", "tabMatrix"), \(from, to="tabMatrix", strict=TRUE) {
  perm <- apply(from, 1L, \(x) {
    nz <- whichv(x, 0, invert=TRUE)
    if (length(nz) == 1L) nz else if (length(nz) == 0L) NA_integer_ else stop("not a tabMatrix")
  })
  num <- reduced <- FALSE
  x <- from[cbind(seq_len(nrow(from)), perm)]
  if (!all(x[!is.na(perm)] == 1)) {
    x[is.na(perm)] <- 0
    perm[is.na(perm)] <- 1L
    num <- TRUE
  } else {
    x <- numeric()
    if (anyNA(perm)) {
      perm[is.na(perm)] <- 0L
      reduced <- TRUE
    }
  }
  tabMatrix(Dim=dim(from), Dimnames=dimnames(from),
            reduced=reduced, perm=perm - 1L, num=num, x=x
  )
})

setMethod("coerce", c("dgCMatrix", "tabMatrix"), \(from, to="tabMatrix", strict=TRUE) {
  perm <- rep.int(-1L, nrow(from))
  x <- numeric(nrow(from))
  nnz <- from@p[-1L] - from@p[-length(from@p)]
  ind <- 1L
  for (co in seq_along(nnz)) {
    if (nnz[co]) {
      ix_ind <- ind:(ind + nnz[co] - 1L)
      perm[from@i[ix_ind][from@x[ix_ind] != 0] + 1L] <- co - 1L
      x[from@i[ix_ind][from@x[ix_ind] != 0] + 1L] <- from@x[ix_ind][from@x[ix_ind] != 0]
      ind <- ind + nnz[co]
    }
  }
  num <- reduced <- FALSE
  if (!all(x[perm >= 0L] %in% c(0, 1))) {
    x[perm < 0L] <- 0
    perm[perm < 0L] <- 0L
    num <- TRUE
  } else {
    if (any(perm < 0L)) {
      perm[x == 0] <- -1L
      reduced <- TRUE
    }
    x <- numeric(0L)
  }
  tabMatrix(Dim=from@Dim, Dimnames=from@Dimnames,
    reduced=reduced, perm=perm, num=num, x=x
  )
})

tab_is_zero <- function(x) (x@reduced && all(x@perm == -1L)) || (x@num && all(x@x == 0))

dgC_is_tabMatrix <- function(M) if (length(M@i) > nrow(M)) FALSE else any_duplicated(M@i)
  
# see whether tabMatrix is actually a permutation matrix
tab_isPermutation <- function(M) {
  M@Dim[1L] == M@Dim[2L] && !M@num && identical(sort(M@perm), 0:(M@Dim[1L] - 1L))
}

#' S4 method for row and column subsetting a 'tabMatrix'
#'
#' @keywords internal
#' @param x a tabMatrix object.
#' @param i integer vector indicating the rows to select. Not used for column subsetting.
#' @param j integer vector indicating the columns to select. Not used for row subsetting.
#' @param ... not used.
#' @param drop whether to return a vector in case of a single selected row or column.
#' @returns The selected rows/columns as a tabMatrix or a vector.
#' @name tabMatrix-indexing
NULL

# auxiliary function to get a positive integer index vector
get_ind <- function(index, M, type="row") {
  if (is.character(index)) {
    ind <- fmatch(index, if (type == "row") dimnames(M)[[1L]] else dimnames(M)[[2L]])
  } else if (is.logical(index)) {
    if (length(index) != M@Dim[if (type == "row") 1L else 2L]) stop("incompatible index vector")
    ind <- base::which(index)
  } else {
    ind <- as.integer(index)
    if (anyNA(ind)) stop("index vector with NAs not allowed")
    if (any(ind < 0L)) {
      if (all(ind < 0L)) {
        ind <- if (type == "row") seq_len(nrow(M))[ind] else seq_len(ncol(M))[ind]
      } else stop("mixed negative and nonnegative index vector")
    }
  }
  if (anyNA(ind) | any(ind == 0L | ind > M@Dim[if (type == "row") 1L else 2L])) stop("index out of range")
  ind
}

tab_row_select <- function(M, i, drop) {
  i <- get_ind(i, M, "row")
  if (drop && length(i) == 1L) {
    out <- numeric(M@Dim[2L])
    if (M@perm[i] >= 0L)
      out[M@perm[i] + 1L] <- if (M@num) M@x[i] else 1
    out
  } else {
    tabMatrix(Dim=c(length(i), M@Dim[2L]), Dimnames=list(M@Dimnames[[1L]][i], M@Dimnames[[2L]]),
      reduced=M@reduced, perm=M@perm[i], num=M@num, x=if (M@num) M@x[i] else numeric(0L)
    )
  }
}

# row selection
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="index", j="missing", drop="logical"), \(x, i, j, ..., drop=TRUE)
  tab_row_select(x, i, drop)
)

# row selection, missing drop
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="index", j="missing", drop="missing"), \(x, i, j, ..., drop=TRUE)
  tab_row_select(x, i, drop=TRUE)
)

tab_col_select <- function(M, j, drop) {
  j <- get_ind(j, M, "col")
  out <- Ctab2dgC(M)[, j, drop=drop]
  if (is.null(dim(out))) out else as(out, "tabMatrix")
}

# column selection
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="missing", j="index", drop="logical"), \(x, i, j, ..., drop=TRUE)
  tab_col_select(x, j, drop)
)

# column selection, missing drop
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="missing", j="index", drop="missing"), \(x, i, j, ..., drop=TRUE)
  tab_col_select(x, j, drop=TRUE)
)


# rbind of tabMatrices is tabMatrix
setMethod("rbind2", c("tabMatrix", "tabMatrix"), \(x, y, ...) {
  if (x@Dim[2L] != y@Dim[2L]) stop("number of columns of matrices must match for rbind2")
  perm <- c(x@perm, y@perm)
  num <- x@num || y@num
  if (num) {
    outx <- c(
      (if (x@reduced) x@perm >= 0L else 1) * (if (x@num) x@x else rep.int(1, x@Dim[1L])),
      (if (y@reduced) y@perm >= 0L else 1) * (if (y@num) y@x else rep.int(1, y@Dim[1L]))
    )
    if (x@reduced || y@reduced) perm[perm < 0L] <- 0L
  }
  out <- Ctab(
    Dim = c(x@Dim[1L] + y@Dim[1L], x@Dim[2L]),
    reduced = !num && (x@reduced || y@reduced),
    perm = perm,
    num = num,
    x = if (num) outx else numeric()
  )
  dnx <- dimnames(x)
  dny <- dimnames(y)
  if (identical(c(dnx, dny), list(NULL, NULL, NULL, NULL))) return(out)
  cn <- if (!is.null(dnx[[2L]])) dnx[[2L]] else dny[[2L]]
  rx <- dnx[[1L]]
  ry <- dny[[1L]]
  rn <- if (is.null(rx) && is.null(ry))
    NULL
  else c(if (!is.null(rx)) rx else character(x@Dim[1L]),
         if (!is.null(ry)) ry else character(y@Dim[1L]))
  out@Dimnames <- list(rn, cn)
  out
})

setMethod("rbind2", c("tabMatrix", "Matrix"), \(x, y, ...)
  rbind2(Ctab2dgC(x), y, ...)
)
setMethod("rbind2", c("Matrix", "tabMatrix"), \(x, y, ...)
  rbind2(x, Ctab2dgC(y), ...)
)
setMethod("cbind2", c("tabMatrix", "Matrix"), \(x, y, ...)
  cbind2(Ctab2dgC(x), y, ...)
)
setMethod("cbind2", c("Matrix", "tabMatrix"), \(x, y, ...)
  cbind2(x, Ctab2dgC(y), ...)
)
setMethod("rbind2", c("tabMatrix", "matrix"), \(x, y, ...)
  rbind2(Ctab2dgC(x), y, ...)
)
setMethod("rbind2", c("matrix", "tabMatrix"), \(x, y, ...)
  rbind2(x, Ctab2dgC(y), ...)
)
setMethod("cbind2", c("tabMatrix", "matrix"), \(x, y, ...)
  cbind2(Ctab2dgC(x), y, ...)
)
setMethod("cbind2", c("matrix", "tabMatrix"), \(x, y, ...)
  cbind2(x, Ctab2dgC(y), ...)
)


setMethod("nnzero", "tabMatrix", \(x, na.counted=NA) {
  if (x@reduced)
    sum(x@perm != -1L)
  else if (x@num)
    sum(x@x != 0)
  else
    x@Dim[1L]
})

#' S4 method for generic 'anyNA' and signature 'tabMatrix'
#'
#' @keywords internal
#' @param x a tabMatrix object.
#' @param recursive not used.
#' @returns Whether the tabMatrix object contains missings or not.
setMethod("anyNA", "tabMatrix", \(x, recursive=FALSE)
  if (x@num) anyNA(x@x) else FALSE
)

setMethod("show", "tabMatrix", \(object) str(object))

# expand a vector, using the indicator part of tabM
expand_mv <- function(tabM, v) Ctab_numeric_prod(tabM, v, TRUE)

remove_levels <- function(f, l=1L) {
  lvs <- attr(f, "levels")
  if (length(l) != 1L || l < 1L || l > length(lvs)) stop("incorrect input")
  lvsnew <- lvs[-l]
  structure(fmatch(lvs, lvsnew)[f], levels=lvsnew, class="factor")
}

# factor (interaction) to tabMatrix
# fvars looked up in data, which can either be a data.frame or an environment
# contrasts: either "contr.treatment", "contr.SAS" or an integer vector of the same length as
#   fvars specifying for each factor variable which level (by name) is considered the baseline.
#   If left unspecified no levels are removed.
# drop.unused.levels: drop unused levels of each factor
# drop: drop empty cells in the interaction (which may result in more columns to be dropped)
fac2tabM <- function(fvars, data, enclos=emptyenv(), x=numeric(), xlab=character(), drop.unused.levels=FALSE, drop=FALSE, contrasts=NULL, varsep=":", catsep="$", lex.order=FALSE) {
  if (missing(fvars) || !length(fvars)) stop("unexpected 'fvars' argument")
  for (f in seq_along(fvars)) {
    fac <- eval_in(fvars[f], data, enclos)
    if (is.factor(fac)) {
      if (drop.unused.levels) fac <- fdroplevels(fac)
    } else {
      fac <- qF(fac)
    }
    if (!is.null(contrasts)) {
      # TODO instead of changing levels here, compute the indices of the final codes vector
      #      which should be set to -1; this would be a lot faster for long factors
      if (length(contrasts) == 1L && contrasts == "contr.treatment") {  # R default: first level is baseline
        fac <- remove_levels(fac, 1L)
      } else if (length(contrasts) == 1L && contrasts == "contr.SAS") {  # 'SAS' convention: last level is baseline
        fac <- remove_levels(fac, length(attr(fac, "levels")))
      } else if (any(fvars[f] == names(contrasts))) {  # user-specified base
        m <- fmatch(contrasts[fvars[f]], attr(fac, "levels"))
        if (length(m) != 1L || is.na(m)) stop("invalid contrasts")
        fac <- remove_levels(fac, m)
      }
    }
    if (f == 1L) {
      out <- fac
      attr(out, "levels") <- paste(fvars[1L], attr(out, "levels"), sep=catsep)
    } else {
      out <- interaction(out, fac, drop=drop, sep=paste0(varsep, fvars[f], catsep), lex.order=lex.order)
    }
  }
  if (length(xlab))
    labs <- paste(attr(out, "levels"), xlab, sep=":")
  else
    labs <- attr(out, "levels")
  out <- as.integer(out) - 1L

  num <- reduced <- FALSE
  if (length(x) && !all(x[!is.na(out)] == 1)) {
    x[is.na(out)] <- 0
    out[is.na(out)] <- 0L
    num <- TRUE
  } else {
    x <- numeric()
    if (anyNA(out)) {
      out[is.na(out)] <- -1L
      reduced <- TRUE
    }
  }

  tabMatrix(Dim=c(length(out), length(labs)), Dimnames=list(NULL, labs),
    reduced=reduced, perm=out, num=num, x=as.double(x), xlab=xlab
  )
}

# formula to list of tabMatrices
# ... passed to fac2tabM
tables2tabM <- function(formula, data, ...) {
  n <- nrow(data)
  trms <- terms(formula, data=data)
  tmat <- terms_matrix(trms)
  if (!length(tmat)) {
    if (intercept_only(formula))
      return(list(new("tabMatrix", perm=integer(n), reduced=FALSE, num=TRUE, x=rep.int(1, n), Dim=c(n, 1L))))
    else
      stop("empty formula")
  }
  tnames <- dimnames(tmat)[[2L]]
  vnames <- dimnames(tmat)[[1L]]
  qvar <- !catvars(trms, data)  # quantitative variables
  qvar <- vnames[which(qvar)]
  out <- setNames(vector(mode="list", length(tnames)), tnames)
  for (k in seq_along(tnames)) {
    countvars <- intersect(vnames[tmat[, k] > 0L], qvar)
    if (length(countvars)) {
      xk <- eval_in(countvars[1L], data)
      if (is.matrix(xk)) stop("matrix variables not yet allowed")
      labk <- countvars[1L]
      for (v in countvars[-1L]) {
        temp <- eval_in(v, data)
        xk <- xk * temp
        labk <- paste(labk, v, sep=":")
      }
    } else {
      xk <- NULL
      labk <- NULL
    }
    facvars <- setdiff(vnames[tmat[, k] > 0L], qvar)
    if (length(facvars)) {
      enclos <- environment(formula)
      if (length(countvars)) {
        fk <- fac2tabM(facvars, data, enclos, x=xk, xlab=labk, contrasts=NULL, ...)
      } else {
        fk <- fac2tabM(facvars, data, enclos, contrasts=NULL, ...)
      }
    } else {
      if (length(countvars)) {
        fk <- new("tabMatrix", perm=integer(n), reduced=FALSE, num=TRUE,
                  x=xk, xlab=labk, Dim=c(n, 1L), Dimnames=list(NULL, labk))
      } else {
        # this should not happen, as the intercept only case is handled above
        fk <- new("tabMatrix", perm=integer(n), reduced=FALSE, num=TRUE,
                  x=rep.int(1, n), Dim=c(n, 1L))
      }
    }
    out[[k]] <- fk
  }
  out
}
