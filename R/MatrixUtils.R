#' Detect and remove redundancy in columns of a matrix object
#'
#' @noRd
#' @param X a (possibly sparse) matrix. If \code{method="chol"} the input matrix \code{X}
#'  for function \code{detect_redundancy} must be symmetric.
#' @param method either "chol" (Cholesky with pivoting) or "qr" (QR decomposition).
#'  Defaults to "chol".
#' @param tol a numeric tolerance value used to test for redundancy. If \code{NULL},
#'  a default value depending on \code{method} is used.
# NB a negative tol value for method "chol" -->
#   default RELATIVE tolerance of ncol(X) * .Machine$double.eps * max(diag(X))
#' @returns For function \code{detect_redundancy} an integer vector indicating columns of \code{X}
#'  that are numerically redundant, and \code{NULL} in case no redundancy is detected.
#'  For function \code{remove_redundancy} the matrix object with redundant columns,
#'  if any, removed.
#' @name redundancy
NULL

## @export
## @rdname redundancy
detect_redundancy <- function(X, method="chol", tol=NULL) {
  if (method == "chol") {
    if (is.null(tol)) tol <- -1
    test <- suppressWarnings(chol.default(X, pivot=TRUE, tol=tol))
    rank <- attr(test, "rank")
    pivot <- attr(test, "pivot")
  } else {
    if (is.null(tol)) tol <- 1e-9
    test <- qr.default(X, tol=tol)
    rank <- test$rank
    pivot <- test$pivot
  }
  if (rank < ncol(X))
    rcols <- pivot[(rank + 1L):ncol(X)]
  else
    rcols <- NULL
  rcols
}

## @export
## @rdname redundancy
remove_redundancy <- function(X, method="chol", tol=NULL) {
  # first remove trivial redundancy: zero columns
  zerocols <- which(zero_col(X))
  if (length(zerocols)) X <- X[, -zerocols, drop=FALSE]
  # remove remaining redundancy
  # assuming X is n x p with n > p --> detect redundancy in X'X typically more efficient
  rcols <- detect_redundancy(crossprod(X), method, tol)
  if (is.null(rcols)) X else X[, -rcols, drop=FALSE]
}


is_a_matrix <- function(x) is.matrix(x) || inherits(x, "Matrix")

is_zero_matrix <- function(x) {
  if (is.matrix(x) || is.vector(x)) return(all(x == 0))
  if (inherits(x, "Matrix")) {
    if (class(x)[1L] == "tabMatrix") return(tab_is_zero(x))
    if (length(x@x))
      all(x@x == 0)
    else
      FALSE  # unit-diagonal
  } else {
    stop("not a supported matrix object")
  }
}

# create a dgCMatrix of zeros
zeroMatrix <- function(nr, nc)
  new("dgCMatrix", p=integer(nc + 1L), Dim=c(as.integer(nr), as.integer(nc)))


################################################################################
# matrix algebra: efficient matrix-vector products and some new Matrix methods

# RcppEigen 'knows' dgC, but not dsC --> represent dsC as dgC
#.dgC.class.attr <- class(.m2sparse(matrix(0), "dgC"))
.dgC.class.attr <- class(as(matrix(c(0, 0)), "CsparseMatrix"))
.dsC.class.attr <- class(forceSymmetric(as(matrix(0), "CsparseMatrix")))

#' Fast matrix-vector multiplications
#' 
#' Functions for matrix-vector multiplies like \code{\%*\%} and \code{crossprod},
#' but often faster for the matrix types supported. The return value is always a
#' numeric vector.
#'
#' @examples
#' M <- matrix(rnorm(10*10), 10, 10)
#' x <- rnorm(10)
#' M %m*v% x
#' crossprod_mv(M, x)
#' M <- Matrix::rsparsematrix(100, 100, nnz=100)
#' x <- rnorm(100)
#' M %m*v% x
#' crossprod_mv(M, x)
#'
#' @param M a matrix of class 'matrix', 'dgCMatrix', 'dsCMatrix', 'tabMatrix', or 'ddiMatrix'.
#' @param v a numeric vector.
#' @returns For \code{\%m*v\%} the vector \eqn{Mv} and for \code{crossprod_mv} the vector
#'  \eqn{M'v} where \eqn{M'} denotes the transpose of \eqn{M}.
#' @name matrix-vector
NULL

#' @export
#' @rdname matrix-vector
`%m*v%` <- function(M, v) {
  switch(class(M)[1L],
    matrix = Cdense_numeric_prod(M, v),
    ddiMatrix = if (length(M@x)) M@x * v else copy_vector(v),
    tabMatrix = Ctab_numeric_prod(M, v),
    dgCMatrix = Csparse_numeric_prod(M, v),
    dsCMatrix = {
      class(M) <- .dgC.class.attr
      CsparseS_numeric_prod(M, v)
    },
    numeric = M * v,  # interpreted as diagonal M, including scalar (1 x 1) case
    stop("unsupported class '", class(M)[1L], "'")
  )
}

#' @export
#' @rdname matrix-vector
crossprod_mv <- function(M, v) {
  switch(class(M)[1L],
    matrix = Cdense_numeric_crossprod(M, v),
    ddiMatrix = if (length(M@x)) M@x * v else copy_vector(v),
    tabMatrix = Ctab_numeric_crossprod(M, v),
    dgCMatrix = Csparse_numeric_crossprod(M, v),
    dsCMatrix = {
      class(M) <- .dgC.class.attr
      CsparseS_numeric_prod(M, v)
    },
    dtCMatrix = {  # possibly used by Ltimes method of cholesky object
      class(M) <- .dgC.class.attr  # convert M to dgC
      Csparse_numeric_crossprod(M, v)
    },
    numeric = M * v,  # interpreted as diagonal M, including scalar (1 x 1) case
    stop("unsupported class '", class(M)[1L], "'")
  )
}


# compute D %*% Q %*% D, where D=Diagonal(x=scale)
# if scale has length 1, it is computed as scale^2 * Q
scale_dsCMatrix <- function(Q, scale) {
  class(Q) <- .dgC.class.attr  # convert Q to dgC
  Q <- Cscale_sparse(Q, scale)
  # convert back to dsC
  attr(Q, "uplo") <- "U"
  class(Q) <- .dsC.class.attr
  Q
}

# faster alternative for scale_dsCMatrix for scale factor compatible with symmetric matrix
# compatible diagonal matrix D=diag(d): [D, M] = 0 s.t. D M = sD M sD where sD = diag(sqrt(d))
# NB low-level function: no checks!
# Q a dsCMatrix and scale a numeric vector
block_scale_dsCMatrix <- function(Q, scale) {
  attr(Q, "x") <- scale[Q@i + 1L] * Q@x
  Q
}

# compute D %*% Q %*% D, where D=Diagonal(x=scale)
# if scale has length 1, it is computed as scale^2 * Q
scale_mat <- function(Q, scale) {
  switch(class(Q)[1L],
    matrix = Cscale_dense(Q, scale),
    numeric = scale^2 * Q,
    ddiMatrix =
      if (Q@diag == "U") {
        attr(Q, "x") <- scale^2
        attr(Q, "diag") <- "N"
        Q
      } else {
        attr(Q, "x") <- scale^2 * Q@x
        Q
      },
    dgCMatrix = Cscale_sparse(Q, scale),
    dsCMatrix = scale_dsCMatrix(Q, scale)
    # TODO tabMatrix
  )
}

# product of two matrices, at least one of which is a dense matrix
# assume that either M1 or M2 is matrix
# numeric vector input represents diagonal matrix
# length 1 numeric vector is allowed: scalar multiplication
`%m*m%` <- function(M1, M2) {
  switch(class(M1)[1L],
    matrix =
      switch(class(M2)[1L],
        matrix = Cdense_dense_prod(M1, M2),
        dgCMatrix = Cdense_sparse_prod(M1, M2),
        dsCMatrix = {
          class(M2) <- .dgC.class.attr
          Cdense_sparseS_prod(M1, M2)
        },
        numeric =
          if (length(M2) == 1L)
            if (M2 == 1) M1 else M2 * M1
          else
            Cdense_diag_prod(M1, M2),
        ddiMatrix =
          if (length(M2@x))
            Cdense_diag_prod(M1, M2@x)
          else
            M1,
        stop("unsupported matrix product")
      ),
    dgCMatrix = Csparse_dense_prod(M1, M2),
    dsCMatrix = {
      class(M1) <- .dgC.class.attr
      CsparseS_dense_prod(M1, M2)
    },
    tabMatrix = Ctab_dense_prod(M1, M2),
    ddiMatrix = if (length(M1@x)) M1@x * M2 else M2,
    stop("unsupported matrix product")
  )
}

# crossprod of two matrices, at least one of which is a dense matrix
# numeric vector input represents diagonal matrix
# length 1 numeric vector is allowed: scalar multiplication
crossprod_mm <- function(M1, M2) {
  switch(class(M1)[1L],
    matrix =
      switch(class(M2)[1L],
        matrix = Cdense_dense_crossprod(M1, M2),
        dgCMatrix = Cdense_sparse_crossprod(M1, M2),
        ddiMatrix =
          if (length(M2@x))
            Cdense_diag_crossprod(M1, M2@x)
              else
            t.default(M1),
        stop("unsupported matrix crossproduct")
      ),
    dgCMatrix = Csparse_dense_crossprod(M1, M2),
    tabMatrix = Ctab_dense_crossprod(M1, M2),
    ddiMatrix = if (length(M1@x)) M1@x * M2 else M2,
    stop("unsupported matrix crossproduct")
  )
}

# Compute the symmetric product M'QM where Q is symmetric
# M can be matrix, ddiMatrix, tabMatrix or dgCMatrix
# Q can be matrix, dsCMatrix, ddiMatrix, or a (nonnegative) vector representing a diagonal matrix
# the result is matrix (if and only if M or Q is matrix), ddiMatrix or dsCMatrix
crossprod_sym <- function(M, Q) {
  classQ <- class(Q)[1L]
  switch(class(M)[1L],
    matrix =
      switch(classQ,
        matrix = Cdense_crossprod_sym2(M, Cdense_dense_prod(Q, M)),
        dsCMatrix = Cdense_crossprod_sym2(M, Q %m*m% M),
        ddiMatrix =
          if (Q@diag == "U")
            Cdense_crossprod_sym0(M)
          else
            Cdense_crossprod_sym(M, Q@x),
        numeric = Cdense_crossprod_sym(M, Q)
      ),
    ddiMatrix = {
      if (M@diag == "U")
        if (classQ == "numeric") {
          attr(M, "diag") <- "N"
          attr(M, "x") <- Q
          return(M)
        } else {
          return(Q)
        }
      switch(classQ,
        dsCMatrix = scale_dsCMatrix(Q, M@x),
        matrix = Cscale_dense(Q, M@x),
        ddiMatrix = {
          attr(M, "x") <- if (Q@diag == "U") M@x^2 else Q@x * M@x^2
          M
        },
        numeric = {
          attr(M, "x") <- Q * M@x^2
          M
        }
      )
    },
    tabMatrix =
      switch(classQ,
        dsCMatrix =
          if (M@isPermutation) {
            class(Q) <- .dgC.class.attr
            out <- Csparse_sym_twist(Q, M@perm)
            attr(out, "uplo") <- "U"
            class(out) <- .dsC.class.attr
            out
          } else {
            crossprod_sym(Ctab2dgC(M), Q)
          },
        matrix = crossprod_sym(Ctab2dgC(M), Q),
        ddiMatrix = {
          if (M@num)
            Q <- if (Q@diag == "N") Q@x * M@x else M@x
          else
            Q <- if (Q@diag == "N") Q@x else rep.int(1, Q@Dim[1L])
          Cdiag(Ctab_numeric_crossprod(M, Q))
        },
        numeric = {
          if (M@num) Q <- Q * M@x
          Cdiag(Ctab_numeric_crossprod(M, Q))
        }
      ),
    dgCMatrix =
      if (classQ == "dsCMatrix") {
        class(Q) <- .dgC.class.attr  # convert Q to dgC
        Q <- Csparse_crossprod_sym(M, Q)
        # convert back to dsC
        attr(Q, "uplo") <- "U"
        class(Q) <- .dsC.class.attr
        Q
      } else if (classQ == "matrix") {
        Csparse_dense_crossprod_sym(M, Q)
      } else {
        if (classQ == "ddiMatrix")
          if (Q@diag == "N")
            Q <- Csparse_diag_crossprod_sym(M, Q@x)
          else
            Q <- Csparse_crossprod_sym2(M, M)
        else  # numeric Q
          Q <- Csparse_diag_crossprod_sym(M, Q)
        attr(Q, "uplo") <- "U"
        class(Q) <- .dsC.class.attr
        Q
      },
    stop("unsupported class '", class(M)[1L], "'")
  )
}

# create a template for fast updating of X'QX under changing
# diagonal matrix Q, and fixed sparse (dgC) X
sparse_crossprod_sym_template <- function(X, max.size.cps.template=100) {
  if (!inherits(X, "dgCMatrix")) stop("matrix of class 'dgCMatrix' expected")
  XtQX <- crossprod_sym(X, 0.1 + runif(nrow(X)))
  XtQX_nsT <- as(as(XtQX, "nMatrix"), "TsparseMatrix")
  nnz_per_col <- Cnnz_per_col_scps_template(X, XtQX_nsT@i, XtQX_nsT@j)
  # check whether the total number of nonzeros is above a certain threshold
  # if so, then we revert to direct computation of the sparse symmetric crossproduct
  # approximate size of the sparse template matrix in MB:
  mem <- ((8 + 4) * sum(nnz_per_col) + length(nnz_per_col)) / 1e6
  if (mem > max.size.cps.template) {
    stop("sparse symmetric crossproduct template requires ", mem, " MB of memory")
  }
  spX <- Ccreate_sparse_crossprod_sym_template(X, XtQX_nsT@i, XtQX_nsT@j, nnz_per_col)
  rm(X, XtQX_nsT, nnz_per_col, mem)
  function(Q) {
    out <- XtQX
    attr(out, "x") <- Csparse_numeric_crossprod(spX, Q)
    out
  }
}

# symmetric crossprod: compute M1'M2, known to be symmetric because M2 = Q M1 with Q symmetric
crossprod_sym2 <- function(M1, M2=M1) {
  switch(class(M2)[1L],
    matrix = Cdense_crossprod_sym2(M1, M2),  # assume that M1 is matrix as well! exploits symmetry of the resulting matrix
    dgCMatrix = {
      # assume that M1 is dgC as well!
      out <- Csparse_crossprod_sym2(M1, M2)  # dgCMatrix with only upper elements stored
      attr(out, "uplo") <- "U"
      class(out) <- .dsC.class.attr
      out
    },
    stop("unsupported class '", class(M2)[1L], "'")
  )
}


ddi_diag <- function(M) if (M@diag == "N") M@x else rep.int(1, M@Dim[1L])


commutator <- function(M1, M2) M1 %*% M2 - M2 %*% M1


#' S4 methods for products of matrix objects
#'
#' Several methods for products of matrix objects. Here a matrix object can be
#' an ordinary (dense) \code{matrix} or a (sparse) \code{\link[Matrix]{Matrix}} of class
#' \code{\link[Matrix]{ddiMatrix-class}}, \code{tabMatrix},
#' \code{\link[Matrix]{dgCMatrix-class}}, or \code{\link[Matrix]{dsCMatrix-class}}.
#' The return types are restricted to \code{matrix} or \code{numeric} in case of dense objects, and
#' \code{dgCMatrix}, \code{dsCMatrix}, \code{ddiMatrix} or \code{tabMatrix} in case of sparse objects.
#'
#' @include tabMatrix.R
#' @keywords internal
#' @param x a matrix object.
#' @param y a matrix object.
#' @returns A matrix object. In case one of the arguments is a regular (dense) \code{matrix}
#'  the result is a \code{matrix} as well.
#' @name Matrix-methods
NULL

#' @rdname Matrix-methods
setMethod("%*%", signature("ddiMatrix", "matrix"), function(x, y) {
  if (x@Dim[2L] != dim(y)[1L]) stop("incompatible dimensions")
  if (x@diag == "U") y else x@x * y
})

# NB the RcppEigen version is faster (especially for small x)
#' @rdname Matrix-methods
setMethod("%*%", signature("dgCMatrix", "matrix"),
  function(x, y) Csparse_dense_prod(x, y)
)

#' @rdname Matrix-methods
setMethod("%*%", signature("dsCMatrix", "matrix"), function(x, y) {
  if (x@uplo != "U") stop("uplo must be 'U'")
  class(x) <- .dgC.class.attr
  CsparseS_dense_prod(x, y)
})

#' @rdname Matrix-methods
setMethod("%*%", signature("tabMatrix", "matrix"),
  function(x, y) Ctab_dense_prod(x, y)
)

#' @rdname Matrix-methods
setMethod("%*%", signature("matrix", "tabMatrix"),
  function(x, y) x %*% Ctab2dgC(y)
)

#' @rdname Matrix-methods
setMethod("%*%", signature("tabMatrix", "numLike"),
  function(x, y) Ctab_dense_prod(x, matrix(y, ncol=1L))
)

#' @rdname Matrix-methods
setMethod("%*%", signature("ddiMatrix", "dgCMatrix"),
  function(x, y) if (x@diag == "U") y else Cdiag_sparse_prod(x@x, y)
)

#' @rdname Matrix-methods
setMethod("%*%", signature("ddiMatrix", "tabMatrix"), function(x, y) {
  if (x@diag == "U") return(y)
  if (y@num) {
    attr(y, "x") <- x@x * y@x
  } else {
    yx <- x@x
    if (y@reduced) {
      yx[y@perm < 0L] <- 0
      y@perm[y@perm < 0L] <- 0L
      attr(y, "reduced") <- FALSE
    }
    attr(y, "x") <- yx
    attr(y, "num") <- TRUE
  }
  y
})

#' @rdname Matrix-methods
setMethod("%*%", signature("CsparseMatrix", "tabMatrix"),
  function(x, y) x %*% Ctab2dgC(y)
)

#' @rdname Matrix-methods
setMethod("%*%", signature("tabMatrix", "CsparseMatrix"),
  function(x, y) Ctab2dgC(x) %*% y
)

#' @rdname Matrix-methods
setMethod("tcrossprod", signature("matrix", "dgCMatrix"),
  function(x, y) Cdense_sparse_tcrossprod(x, y)
)

#' @rdname Matrix-methods
setMethod("tcrossprod", signature("matrix", "tabMatrix"),
  function(x, y) Cdense_tab_tcrossprod(x, y)
)

#' @rdname Matrix-methods
setMethod("tcrossprod", signature("matrix", "ddiMatrix"), function(x, y) {
  if (dim(x)[2L] != y@Dim[2L]) stop("incompatible dimensions")
  if (y@diag == "U") x else x * rep_each(y@x, dim(x)[1L])
})

# (unary) crossprod method for tabMatrix
#' @rdname Matrix-methods
setMethod("crossprod", signature=c("tabMatrix", "missing"),
  function(x, y) Cdiag(Ctab_unary_crossprod(x))
)

#' @rdname Matrix-methods
setMethod("crossprod", signature("tabMatrix", "matrix"),
  function(x, y) Ctab_dense_crossprod(x, y)
)

#' @rdname Matrix-methods
setMethod("crossprod", signature("tabMatrix", "dgCMatrix"),
  function(x, y) crossprod(Ctab2dgC(x), y)
)

#' @rdname Matrix-methods
setMethod("crossprod", signature("tabMatrix", "tabMatrix"),
  function(x, y) crossprod(Ctab2dgC(x), Ctab2dgC(y))
)

# used in fast GMRF prior sampler
#' @rdname Matrix-methods
setMethod("crossprod", signature=c("dgCMatrix", "matrix"),
  function(x, y) Csparse_dense_crossprod(x, y)
)


is_unit_ddi <- function(M) class(M)[1L] == "ddiMatrix" && M@diag == "U"

# test for unit diagonal matrix
# assume that matrix is economized, i.e. either ddi (sparse) or matrix (dense)
is_unit_diag <- function(M) {
  if (is.matrix(M))
    all(diag(M) == 1) && isDiagonal(M)
  else
    is_unit_ddi(M)
}

# make sure ddiMatrix has x-slot
expand_unit_ddi <- function(M) {
  attr(M, "x") <- rep.int(1, M@Dim[1L])
  attr(M, "diag") <- "N"
  M
}

rm_Matrix_names <- function(M) {
  if (!all(b_apply(M@Dimnames, is.null))) attr(M, "Dimnames") <- list(NULL, NULL)
  M
}

# 'economize' a diagonal matrix; M is assumed to be a diagonal matrix
# if (vec.diag) returns a numeric vector, in case all elements equal a numeric scalar
# else returns a ddiMatrix, if possible a unit-ddiMatrix
economizeDiagMatrix <- function(M, strip.names=TRUE, vec.diag=FALSE) {
  M <- as(M, "diagonalMatrix")
  if (M@diag == "N" && all(M@x == 1)) {
    attr(M, "x") <- numeric(0L)
    attr(M, "diag") <- "U"
  }
  if (vec.diag) {  # NB in scalar case names are not returned
    if (M@diag == "U")
      1
    else if (all(M@x == M@x[1L]))
      M@x[1L]
    else if (strip.names)
      M@x
    else setNames(M@x, colnames(M))
  } else {
    if (strip.names)
      M <- rm_Matrix_names(M)
    M
  }
}

#' Transform a matrix into an efficient representation, if possible.
#'
#' Supported matrix representations are \code{ddiMatrix} (unit or not),
#' tabMatrix (only if \code{symmetric=FALSE}), \code{dgCMatrix}
#' (\code{dsCMatrix} if \code{symmetric=TRUE}) and matrix.
#' For sparse matrices the order of decreasing efficiency is: unit-ddi,
#' ddi, tab, dgC/dsC.
# TODO add option tol to pass to drop0 to set numerically very small coefficients to zero (or zapsmall for relative criterion)
#'
#' @noRd
#' @param M the matrix to be 'economized'.
#' @param sparse whether the matrix should be sparse. If \code{NULL} (default) a simple heuristic is used to make a choice.
#' @param symmetric whether the matrix is symmetric. This is relevant for dgC matrices, which are then better represented as dsC.
#' @param strip.names whether to remove dimnames.
#' @param allow.tabMatrix whether \code{tabMatrix} is allowed as return type. Default is \code{TRUE}.
#' @param drop.zeros whether explicit zeros should be dropped. Only relevant for dgC/dsCMatrix. \code{FALSE} by default.
#' @param vec.diag if \code{TRUE} a diagonal matrix is represented by a vector. In addition, if all elements
#'  of the vector are equal then it is reduced to a scalar.
#' @param vec.as.diag if \code{TRUE}, the default, vector input is expanded to a diagonal matrix with
#'  this vector at its diagonal. Otherwise, vector \code{M} is interpreted as a single-column matrix.
#' @param check if \code{TRUE}, do some checks on the input matrix, such as a check for NAs and a
#'  symmetry check in case \code{symmetric=TRUE}.
#' @returns The matrix in an efficient (memory and performance-wise) format.
economizeMatrix <- function(M, sparse=NULL, symmetric=FALSE, strip.names=TRUE,
                            allow.tabMatrix=TRUE, drop.zeros=FALSE, vec.diag=FALSE, vec.as.diag=TRUE, check=FALSE) {
  if (is.vector(M)) {
    if (vec.as.diag) M <- Cdiag(M) else M <- matrix(M, ncol=1L)
  }
  if (check && !is_a_matrix(M)) stop("not a matrix or Matrix")
  if (check && anyNA(M)) stop("economizeMatrix: missing values in input matrix")
  if (symmetric && !isFALSE(sparse) && isDiagonal(M))
    return(economizeDiagMatrix(M, strip.names, vec.diag))
  if (is.null(sparse)) sparse <- better_sparse(M, symmetric)
  if (sparse) {
    #if (is.matrix(M)) M <- .m2sparse(M, "dgC")
    if (is.matrix(M)) M <- as(as(as(M, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  } else {
    if (!is.matrix(M)) M <- as.matrix(M)
    # make sure the matrix has mode 'double', otherwise some C++ matrix routines might not work
    if (!is.double(M)) storage.mode(M) <- "double"
    if (check && symmetric && !isSymmetric(M)) stop("economizeMatrix: matrix is not symmetric")
    return(if (strip.names) unname(M) else M)
  }
  if (isDiagonal(M)) return(economizeDiagMatrix(M, strip.names, vec.diag))
  if (allow.tabMatrix && !symmetric) {
    if (class(M)[1L] != "tabMatrix") {
      M <- as(as(as(M, "CsparseMatrix"), "generalMatrix"), "dMatrix")
      if (dgC_is_tabMatrix(M)) M <- as(M, "tabMatrix")
    }
    if (class(M)[1L] == "tabMatrix") {
      attr(M, "isPermutation") <- tab_isPermutation(M)
      return(if (strip.names) rm_Matrix_names(M) else M)
    }
  }
  if (all(class(M)[1L] != c("dgCMatrix", "dsCMatrix")))
    M <- as(as(as(M, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  if (drop.zeros) M <- drop0(M, is.Csparse=TRUE)
  if (symmetric && class(M)[1L] == "dgCMatrix") {
    if (check && !isSymmetric(M)) stop("economizeMatrix: matrix is not symmetric")
    M <- forceSymmetric(M, uplo="U")
  }
  if (class(M)[1L] == "dsCMatrix" && M@uplo != "U") stop("need uplo='U' for dsCMatrix")
  if (strip.names) rm_Matrix_names(M) else M
}

#' Compute fraction of zero matrix elements
#'
#' @noRd
#' @param x a (possibly sparse) matrix object.
#' @returns The fraction of zero matrix elements.
sparsity <- function(x) 1 - nnzero(x) / prod(dim(x))

#' Heuristic for choosing between sparse Matrix and (dense) matrix
#'
#' @noRd
#' @param x a matrix or Matrix.
#' @returns Whether, based on a simple heuristic, a sparse matrix representation would probably be more efficient.
better_sparse <- function(x, symmetric=FALSE) {
  if (symmetric) {
    # the packed dsC format is another advantage compared to matrix
    sparsity(x) > 0.25 && prod(dim(x)) > 500
  } else {
    sparsity(x) > 0.5 && prod(dim(x)) > 500
  }
}

# assume x is a matrix object; if dense, it is assumed to be matrix
large_and_sparse <- function(x) !is.matrix(x) && prod(dim(x)) > 1e6

setMethod("solve", signature("ddiMatrix", "numeric"), function(a, b)
  switch(a@diag,
    U = b,
    N = b / a@x
  )
)

setMethod("solve", signature("ddiMatrix", "matrix"), function(a, b)
  switch(a@diag,
    U = b,
    N = b / a@x
  )
)

#' Combine factors into a single factor variable representing the interaction of the factors
#'
#' @noRd
#' @param fvars vector of factor names.
#' @param data data frame in which to look up the factor names.
#' @param drop whether to drop unobserved levels in the resulting factor.
#' @param sep the levels of the combined factor are obtained by concatenating the levels
#'  of all factors using this separator. Default is ':'.
#' @param lex.order passed to \code{\link[base]{interaction}}.
#' @param enclos enclosure to look for objects not found in \code{data}.
#' @returns The combined (interaction) factor variable.
combine_factors <- function(fvars, data, drop=FALSE, sep=":", lex.order=FALSE, enclos=emptyenv()) {
  if (length(fvars)) {
    if (!is.character(fvars)) stop("'fvars' must be a character vector")
    fac <- eval_in(fvars[1L], data, enclos)
    # keep all levels, including those of empty combinations, even combinations of empty levels
    for (f in fvars[-1L])
      fac <- interaction(fac, eval_in(f, data, enclos), drop=drop, sep=sep, lex.order=lex.order)
  } else {
    fac <- NULL
  }
  fac
}

# take the kronecker product of two sparse precision or incidence matrices
# the result is always sparse (ddi, dsC in case of precision matrix or dgC in case of non-symmetric matrix)
# in kronecker product the indices of the first factor run slowest
#   --> kronecker(Q2, Q1) to match the order of interaction(f1, f2)
cross <- function(Q1, Q2) {
  c1 <- class(Q1)[1L]
  c2 <- class(Q2)[1L]
  if ((c1 == "ddiMatrix") && (c2 == "ddiMatrix")) {
    if (Q1@diag == "U" && Q2@diag == "U")
      return(CdiagU(nrow(Q1)*nrow(Q2)))
    else
      return(Cdiag(as.numeric(base::tcrossprod(ddi_diag(Q1), ddi_diag(Q2)))))
  }
  if (all(c(c1, c2) %in% c("ddiMatrix", "dsCMatrix"))) return(forceSymmetric(as(kronecker(Q2, Q1), "CsparseMatrix"), uplo="U"))
  return(as(kronecker(Q2, Q1), "CsparseMatrix"))  # at least one matrix is dgC
}

#' Utility function to construct a sparse aggregation matrix from a factor
#'
#' @examples
#' n <- 1000
#' f <- sample(1:100, n, replace=TRUE)
#' x <- runif(n)
#' M <- aggrMatrix(f)
#' all.equal(crossprod_mv(M, x), as.vector(tapply(x, f, sum)))
#'
#' @export
#' @param fac factor variable.
#' @param w vector of weights associated with the levels of \code{fac}.
#' @param mean if \code{TRUE}, aggregation will produce (weighted) means instead of sums.
#' @param facnames whether the factor levels should be used as column names for the aggregation matrix.
#' @returns A sparse aggregation matrix of class \code{tabMatrix}.
aggrMatrix <- function(fac, w=1, mean=FALSE, facnames=FALSE) {
  n <- length(fac)
  if (n == 0L) stop("empty factor variable")
  if (!identical(w, 1) && length(w) != n) stop("fac and w must have the same length")
  fac <- as.factor(fac)
  levs <- attr(fac, "levels")
  if (n == length(levs) && all(fac == levs)) {
    Maggr <- if (all(w == 1) || mean) CdiagU(n) else Cdiag(w)
  } else {
    if (mean) {
      if (length(w) == n)
        w <- w / fast_aggrC(w, fac, length(levs))[as.integer(fac)]
      else
        w <- 1 / tabulate(fac)[as.integer(fac)]
    }
    if (all(w == 1)) {
      Maggr <- tabMatrix(Dim=c(n, length(levs)), reduced=FALSE, perm=as.integer(fac) - 1L, num=FALSE)
      attr(Maggr, "isPermutation") <- tab_isPermutation(Maggr)
    } else {
      Maggr <- tabMatrix(Dim=c(n, length(levs)), reduced=FALSE, perm=as.integer(fac) - 1L, num=TRUE, x=w)
    }
  }
  if (facnames) attr(Maggr, "Dimnames") <- list(NULL, levs)
  Maggr
}

# apply a function to the NON-ZERO values by column of a dgCMatrix
# M: dgCMatrix
# fun: vector to scalar function
dgC_colwise <- function(M, fun) {
  out <- numeric(ncol(M))
  colsizes <- diff(M@p)
  ind <- 1L
  for (i in seq_along(colsizes)) {
    if (colsizes[i] > 0L)
      out[i] <- fun(M@x[ind:(ind + colsizes[i] - 1L)])
    ind <- ind + colsizes[i]
  }
  out
}

# returns a logical vector indicating whether a column is zero or not
# TODO use relative instead of absolute tolerance (?), or mean (over nonzero elements) instead of sum?
zero_col <- function(M, tol=sqrt(.Machine$double.eps)) {
  f <- function(x) sum(abs(x))
  switch(class(M)[1L],
    matrix = apply(M, 2L, f) < tol,
    ddiMatrix = if (M@diag == "U") rep(FALSE, M@Dim[2L]) else abs(M@x) < tol,
    tabMatrix =
      if (M@num)
        as.numeric(tapply(M@x, 0:(M@Dim[2L] - 1L), f)) < tol
      else
        tabulate(M@perm + 1L, nbins=M@Dim[2L]) == 0L,
    dgCMatrix = dgC_colwise(M, f) < tol,
    stop("unsupported class '", class(M)[1L], "'")
  )
}

colwise_maxabs <- function(M) {
  f <- function(x) max(abs(x))
  switch(class(M)[1L],
    matrix = apply(M, 2L, f),
    ddiMatrix = if (M@diag == "U") rep.int(1, M@Dim[2L]) else abs(M@x),
    tabMatrix =
      if (M@num)
        as.numeric(tapply(M@x, 0:(M@Dim[2L] - 1L), f))
      else
        as.numeric(tabulate(M@perm + 1L, nbins=M@Dim[2L]) != 0L),
    dgCMatrix = dgC_colwise(M, f),
    stop("unsupported class '", class(M)[1L], "'")
  )
}
