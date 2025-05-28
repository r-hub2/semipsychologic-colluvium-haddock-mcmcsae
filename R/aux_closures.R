
# get column indices of CsparseMatrix (compressed column-oriented)
get_col_ind <- function(M, zero.based=TRUE) {
  if (zero.based)
    rep(seq_len(dim(M)[2L]) - 1L, diff.default(M@p))
  else
    rep(seq_len(dim(M)[2L]), diff.default(M@p))
}

# output type determined by input types and force.sparse, not by economizeMatrix
make_mat_sum <- function(M0 = NULL, M1, M2 = NULL, force.sparse = FALSE) {
  q <- nrow(M1)
  if (ncol(M1) != q) stop("M1 is not symmetric")
  doM0 <- !is.null(M0)
  doM2 <- !is.null(M2)
  if (doM0 && !identical(dim(M0), c(q, q))) stop("dimension of M0 differs from dimension of M1")
  if (doM2 && !identical(dim(M2), c(q, q))) stop("dimension of M2 differs from dimension of M1")
  # check that matrix input is symmetric?
  class0 <- substr(class(M0)[1L], 1L, 3L)
  class1 <- substr(class(M1)[1L], 1L, 3L)
  class2 <- substr(class(M2)[1L], 1L, 3L)
  any.mat <- any(c(class0, class1, class2) == "mat")
  sparse.output <- force.sparse || !any.mat
  need.template <- doM0 || sparse.output
  if (need.template) {
    if (sparse.output) {
      template <- economizeMatrix(
        if (doM0 && doM2)
          M0 + M1 + M2
        else if (doM0)
          M0 + M1
        else
          M1 + M2,
        sparse=sparse.output, symmetric=TRUE
      )
    } else {
      template <- as.matrix(M0)
    }
  }

  if (need.template && class(template)[1L] == "ddiMatrix") {
    # ddi output
    # some input matrices might be (diagonal) dsC, but not (diagonal) matrix
    if (any.mat) stop("unexpected input for make_mat_sum")
    if (is_unit_ddi(template)) stop("unexpected unit diagonal matrix sum")
    UnitdiagM1 <- is_unit_ddi(M1)
    sub1 <- class1 == "dsC" && length(M1@x) < q
    if (!UnitdiagM1 && sub1) {
      # diagonal dsC containing zeros
      ind1 <- fmatch(M1@i, seq_len(q) - 1L)
    }
    if (doM2) {
      UnitdiagM2 <- is_unit_ddi(M2)
      sub2 <- class2 == "dsC" && length(M2@x) < q
      if (!UnitdiagM2 && sub2) {
        # diagonal dsC containing zeros
        ind2 <- fmatch(M2@i, seq_len(q) - 1L)
      }
    }
    if (doM0) {
      if (class0 == "dsC")
        template <- as(M0, "diagonalMatrix")
      else
        template <- M0
      if (is_unit_ddi(template)) template <- expand_unit_ddi(template)
    } else {
      template <- Cdiag(numeric(q))
    }
    update <- function(M1, M2, w1=1, w2=1) {
      out <- template
      if (sub1) {
        x <- numeric(q)
        x[ind1] <- w1 * M1@x
      } else {
        x <- if (UnitdiagM1) w1 else w1 * M1@x
      }
      if (doM2) {
        if (sub2)
          x[ind2] <- x[ind2] + w2 * M2@x
        else
          x <- x + if (UnitdiagM2) w2 else w2 * M2@x
      }
      attr(out, "x") <- if (doM0) template@x + x else x
      out
    }
  } else if (sparse.output) {
    # dsC output
    UnitDiag1 <- is_unit_ddi(M1)
    UnitDiag2 <- doM2 && is_unit_ddi(M2)
    nx <- length(template@i)
    j <- get_col_ind(template)
    switch(class1,
      ddi = {
        ind1 <- which(template@i == j & ddi_diag(M1)[j + 1L] != 0) - 1L
        if (anyv(M1@x, 0)) {
          nz1 <- whichv(M1@x, 0, invert=TRUE)
          getM1x <- function(M) M@x[nz1]
        } else {
          getM1x <- function(M) M@x
        }
      },
      mat = {
        Mind1 <- which(upper.tri(M1, diag=TRUE) & M1 != 0)
        ind1 <- fmatch(Mind1 - 1L, template@i + q * j) - 1L
        getM1x <- function(M) M[Mind1]
      },
      dsC = {
        ind1 <- fmatch(list(M1@i, get_col_ind(M1)), list(template@i, j)) - 1L
        getM1x <- function(M) M@x
      }
    )
    if (doM2) {
      switch(class2,
        ddi = {
          ind2 <- which(template@i == j & ddi_diag(M2)[j + 1L] != 0) - 1L
          if (anyv(M2@x, 0)) {
            nz2 <- whichv(M2@x, 0, invert=TRUE)
            getM2x <- function(M) M@x[nz2]
          } else {
            getM2x <- function(M) M@x
          }
        },
        mat = {
          Mind2 <- which(upper.tri(M2, diag=TRUE) & M2 != 0)
          ind2 <- fmatch(Mind2 - 1L, template@i + q * j) - 1L
          getM2x <- function(M) M[Mind2]
          if (class1 == "mat") {
            # additional check for cancellations in M1 + M2 --> invalid sum template
            if (anyNA(ind2)) stop("invalid matrix sum template")
            # maybe randomize the (nonzeros of the) input matrices to prevent cancellations
          }
        },
        dsC = {
          ind2 <- fmatch(list(M2@i, get_col_ind(M2)), list(template@i, j)) - 1L
          getM2x <- function(M) M@x
        }
      )
    } else {
      ind2 <- integer()
      getM2x <- function(M) numeric()
    }
    if (doM0) {
      attr(template, "x") <- template@x - sparse_sum_x(nx, ind1, ind2, getM1x(M1), getM2x(M2), UnitDiag1, UnitDiag2, 1, 1)
      update <- function(M1, M2, w1=1, w2=1) {
        out <- template
        attr(out, "x") <- out@x + sparse_sum_x(nx, ind1, ind2, getM1x(M1), getM2x(M2), UnitDiag1, UnitDiag2, w1, w2)
        out
      }
    } else {
      attr(template, "x") <- numeric()
      update <- function(M1, M2, w1=1, w2=1) {
        out <- template
        attr(out, "x") <- sparse_sum_x(nx, ind1, ind2, getM1x(M1), getM2x(M2), UnitDiag1, UnitDiag2, w1, w2)
        out
      }
    }
  } else {
    # dense (matrix) output
    update <- function(M1, M2, w1=1, w2=1) {}
    if (doM0)
      update <- add(update, quote(x <- template))
    else if (class1 == "mat")
      update <- add(update, quote(x <- w1 * M1))
    else
      update <- add(update, quote(x <- w2 * M2))
    if (doM0 || class1 != "mat") {
      if (class1 == "mat") {
        update <- add(update, quote(x <- x + w1 * M1))
      } else if (class1 == "ddi") {
        update <- add(update, quote(x <- add_diagC(x, w1 * M1@x)))
      } else {  # dsC
        j1 <- get_col_ind(M1)
        ind1 <- 1L + M1@i + q * j1
        update <- add(update, quote(x[ind1] <- x[ind1] + w1 * M1@x))
        ind1.nd <- which(j1 != M1@i)
        ind1.lower <- 1L + j1[ind1.nd] + q * M1@i[ind1.nd]
        update <- add(update, quote(x[ind1.lower] <- x[ind1.lower] + w1 * M1@x[ind1.nd]))
      }
    }
    if ((doM0 && doM2) || (!doM0 && class1 == "mat")) {
      if (class2 == "mat") {
        update <- add(update, quote(x <- x + w2 * M2))
      } else if (class2 == "ddi") {
        update <- add(update, quote(x <- add_diagC(x, w2 * M2@x)))
      } else {  # dsC
        j2 <- get_col_ind(M2)
        ind2 <- 1L + M2@i + q * j2
        update <- add(update, quote(x[ind2] <- x[ind2] + w2 * M2@x))
        ind2.nd <- which(j2 != M2@i)
        ind2.lower <- 1L + j2[ind2.nd] + q * M2@i[ind2.nd]
        update <- add(update, quote(x[ind2.lower] <- x[ind2.lower] + w2 * M2@x[ind2.nd]))
      }
    }
    update <- add(update, quote(x))
  }
  rm(M0, M1, M2, class0, class1, class2)
  update
}

# creates a function to update the log-determinant
make_det <- function(M, chol.control=chol_control(perm=FALSE)) {
  if (nrow(M) <= 1000L) {
    # for moderate dimensions eigenvalues are the fastest way for simple determinant updates
    ev0 <- eigen(M, only.values=TRUE)$values
    rm(M, chol.control)
    function(w1, w2) sum(log(w1 * ev0 + w2))
  } else {
    d <- nrow(M)
    if (is.matrix(M)) {
      function(w1, w2) d * log(w1) +
        as.numeric(determinant.matrix(Ccholesky(add_diagC(M, rep.int(w2/w1, d))), logarithm=TRUE)[["modulus"]])
    } else if (class(M)[1L] == "dsCMatrix") {
      ch <- build_chol(M, control=chol.control)
      rm(chol.control)
      function(w1, w2) d * log(w1) +
        as.numeric(determinant(.updateCHMfactor(ch[["cholM"]], M, mult=w2/w1), logarithm=TRUE, sqrt=FALSE)[["modulus"]])
    } else stop("unexpected matrix type")
    #function(w1, w2) d * log(w1) + 2 * detchol$logdet(mult=w2/w1)
  }
}
