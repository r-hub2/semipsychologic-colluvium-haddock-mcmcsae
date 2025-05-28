
.solve.systems <- setNames(c(1L, 5L, 6L), c("A", "L", "Lt"))

# define R function to set defaults
# ordering=1L (AMD only) seems to give same permutations as Matrix::Cholesky
# ordering=0L is the cholmod default and also tries nested dissection
Cholesky_dsC <- function(A, perm=TRUE, super=NA, Imult=0, ordering=0L, LDL=FALSE)
  cCHM_dsC_Cholesky(A, perm, super, Imult, ordering, LDL)

#' Set options for Cholesky decomposition
#'
#' These options are only effective in case the matrix to be decomposed is sparse, i.p.
#' of class \code{\link[Matrix]{dsCMatrix-class}}.
#'
#' @export
#' @param perm logical scalar, see \code{\link[Matrix]{Cholesky}}. If \code{NULL}, the default,
#'  the choice is left to a simple heuristic.
#' @param super logical scalar, see \code{\link[Matrix]{Cholesky}}.
#' @param ordering an integer scalar passed to CHOLMOD routines determining which reordering
#'  schemes are tried to limit sparse Cholesky fill-in.
#' @param inplace whether sparse Cholesky updates should re-use the same memory location.
#' @returns A list with specified options used for Cholesky decomposition.
chol_control <- function(perm=NULL, super=NA, ordering=0L, inplace=TRUE)
  list(perm=perm, super=super, ordering=ordering, inplace=inplace)

check_chol_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- chol_control()
  w <- which(!(names(control) %in% names(defaults)))
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$ordering <- as.integer(control[["ordering"]])
  control
}

# in case of permutation/pivoting, the decomposition is PtLLtP
# for LDL there is currently no Ltimes method and no solve method for systems other than "A"
build_chol <- function(M, Imult=0, control=chol_control(), LDL=FALSE) {
  type <- class(M)[1L]
  if (all(type != c("numeric", "matrix", "ddiMatrix", "dsCMatrix"))) stop("'", type, "' is not a supported matrix class")
  size <- if (type == "numeric") length(M) else dim(M)[1L]

  if (type == "dsCMatrix") {
    if (LDL) control$super <- FALSE  # not supported by Matrix/CHOLMOD
    if (is.null(control[["perm"]])) {
      # is there a better decision rule for using permutation in cholesky factorisation?
      control$perm <- size > 100L
    }
  } else {
    if (LDL) stop("LDL only supported for dsCMatrix")
    control$perm <- FALSE
  }

  # define update, solve and Ltimes methods with arguments
  # update
  #   parent: new matrix to be decomposed; in case of sparse matrix the same zero pattern is assumed
  #   mult: mult*I is added to parent before decomposition
  # solve
  #   rhs: right hand side, of type numeric, matrix or dgCMatrix
  #   system: the type of system to be solved, see Matrix package
  # Ltimes left-multiplies a vector/matrix by L (in case of permutations P'L which need not be lower triangular!) or its transpose
  # logdet: return logarithm of determinant of the Cholesky factor
  #         NB 1) in case of LDL the logarithm of the absolute value of the determinant is returned
  #         NB 2) for the determinant of M itself, multiply the return value by 2
  # inverse: return inverse of the original matrix
  switch(type,
    numeric=, ddiMatrix = {  # diagonal represented as vector case (including trivial scalar case)
      if (type == "ddiMatrix" && M@diag == "U") {  # trivial case, identity matrix, no updates
        if (Imult != 0) stop("unit ddiMatrix assumed fixed")
        cholM <- M
        update <- function(parent, mult=0) if (mult != 0) stop("unit ddiMatrix assumed fixed")
        Ltimes <- function(x, transpose=TRUE) x
        logdet <- function() 0
        solve <- function(rhs, system="A") rhs
        inverse <- function() cholM
      } else {
        if (type == "ddiMatrix") {
          cholM <- sqrt(M@x + Imult)  # a vector is sufficient as internal representation
          update <- function(parent, mult=0) cholM <<- sqrt(parent@x + mult)
        } else {
          cholM <- sqrt(M + Imult)
          update <- function(parent, mult=0) cholM <<- sqrt(parent + mult)
        }
        Ltimes <- function(x, transpose=TRUE) cholM * x
        logdet <- function() sum(log(cholM))
        solve <- function(rhs, system="A") {
          switch(class(rhs)[1L],
            numeric=, matrix =
              switch(system,
                A = rhs / cholM^2,
                Lt=, L = rhs / cholM
              ),
            dgCMatrix =
              switch(system,
                A = {
                  attr(rhs, "x") <- rhs@x / cholM[rhs@i + 1L]^2
                  rhs
                },
                Lt=, L = {
                  attr(rhs, "x") <- rhs@x / cholM[rhs@i + 1L]
                  rhs
                }
              ),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
        }
        inverse <- function() Cdiag(1 / cholM^2)
      }
    },
    dsCMatrix = {
      cholM <- Cholesky_dsC(M, perm=control$perm, super=control$super, Imult=Imult, ordering=control$ordering, LDL=LDL)
      perm <- is.unsorted(cholM@perm, strictly = TRUE)
      if (perm) {  # permutation matrix unaffected by updates
        P <- cholM@perm + 1L
        iP <- invertPerm(cholM@perm, off=0L, ioff=1L)
      }
      if (control[["inplace"]]) {
        # requires that zero-pattern does not change!
        update <- function(parent, mult=0) cCHM_update_inplace(cholM, parent, mult)
      } else {
        # for forked parallel processing cholM should probably be stored in the state list p
        update <- function(parent, mult=0) cholM <<- .updateCHMfactor(cholM, parent, mult)
      }
      if (perm) {
        # TODO
        # - check whether the permutation is handled correctly in the transpose=FALSE case
        # - use update flag to avoid unnecessary calls to expand1
        Ltimes <- function(x, transpose=TRUE) {
          L <- expand1(cholM, which="L")
          if (transpose)
            if (is.vector(x)) crossprod_mv(L, x[P]) else crossprod_mv(L, x[P, , drop=FALSE])
          else
            if (is.vector(x)) (L %m*v% x)[iP] else (L %m*v% x)[iP, , drop=FALSE]
        }
        solve <- function(rhs, system="A") {
          switch(class(rhs)[1L],
            numeric =
              switch(system,
                A = cCHMf_solve(cholM, rhs, 1L),
                Lt = cCHMf_solve(cholM, rhs, 6L)[iP],
                L = cCHMf_solve(cholM, rhs[P], 5L)
              ),
            matrix =
              switch(system,
                A = cCHMf_solve_matrix(cholM, rhs, 1L),
                Lt = cCHMf_solve_matrix(cholM, rhs, 6L)[iP, , drop=FALSE],
                L = cCHMf_solve_matrix(cholM, rhs[P, , drop=FALSE], 5L)
              ),
            dgCMatrix =
              switch(system,
                A = cCHMf_spsolve(cholM, rhs, 1L),
                Lt = cCHMf_spsolve(cholM, rhs, 6L)[iP, , drop=FALSE],
                L = cCHMf_spsolve(cholM, rhs[P, , drop=FALSE], 5L)
              ),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
        }
      } else {
        Ltimes <- function(x, transpose=TRUE) {
          L <- expand1(cholM, which="L")
          if (transpose)
            if (is.vector(x)) crossprod_mv(L, x) else crossprod_mm(L, x)
          else
            if (is.vector(x)) L %m*v% x else L %m*m% x
        }
        solve <- function(rhs, system="A")
          switch(class(rhs)[1L],
            numeric = cCHMf_solve(cholM, rhs, .solve.systems[system]),
            matrix = cCHMf_solve_matrix(cholM, rhs, .solve.systems[system]),
            dgCMatrix = cCHMf_spsolve(cholM, rhs, .solve.systems[system]),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
      }
      # sqrt=TRUE (without factor 0.5) yields NaN for LDL factors of negative definite matrices
      logdet <- function() 0.5 * as.numeric(determinant(cholM, logarithm=TRUE, sqrt=FALSE)[["modulus"]])
      inverse <- function() Matrix::solve(cholM)
    },
    matrix = {
      # use chol.default to initialize as it signals non-positive-definiteness
      # use faster Ccholesky for update
      # both return upper triangular matrix, i.e. Lt
      if (Imult == 0)
        cholM <- chol.default(M)
      else
        cholM <- chol.default(add_diagC(M, rep.int(Imult, size)))
      update <- function(parent, mult=0) {
        if (mult == 0)
          cholM <<- Ccholesky(parent)
        else
          cholM <<- Ccholesky(add_diagC(parent, rep.int(mult, size)))
      }
      Ltimes <- function(x, transpose=TRUE) {
        if (transpose)
          if (is.vector(x)) cholM %m*v% x else cholM %m*m% x
        else
          if (is.vector(x)) crossprod_mv(cholM, x) else crossprod_mm(cholM, x)
      }
      logdet <- function() as.numeric(determinant.matrix(cholM, logarithm=TRUE)[["modulus"]])
      solve <- function(rhs, system="A") {
        switch(class(rhs)[1L],
          numeric =
            switch(system,
              A = Cbacksolve(cholM, Cforwardsolve(cholM, rhs)),
              Lt = Cbacksolve(cholM, rhs),
              L = Cforwardsolve(cholM, rhs),
              stop("unsupported solve system")
            ),
          matrix =
            switch(system,
              A = CbacksolveM(cholM, CforwardsolveM(cholM, rhs)),
              Lt = CbacksolveM(cholM, rhs),
              L = CforwardsolveM(cholM, rhs),
              stop("unsupported solve system")
            ),
          # disallow dgCMatrix rhs with dense M; in this case it is better to force rhs to be dense as well, as the result is typically dense
          #dgCMatrix = {  # this works, but may be slow; TODO instead of dgC, return the result as a dense matrix
          #  nc <- dim(rhs)[2L]
          #  out <- rhs
          #  for (i in seq_len(nc)) out[, i] <- solve(rhs[, i], system)
          #  out
          #}
          stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
        )
      }
      inverse <- function() chol2inv(cholM)
    }
  )  # END switch(type, ...)

  rm(M, Imult, control)
  environment()
}
