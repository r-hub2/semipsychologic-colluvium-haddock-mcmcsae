
# allowed M1: ddi, dsC, mat
# allowed M2: num, mat (numeric is interpreted as a diagonal matrix)
# q2 is the dimension of M2 (necessary in case M2 is represented by a numeric scalar)
# returns a closure which yields the kronecker product as a ddi, dsC/dgC or ordinary matrix
# all matrices assumed symmetric; one exception: dgC M1, and num M2; in that case dgC result
# NB zero structure assumed to be constant, even for matrix M1
build_kron <- function(M1, M2, q2, M1.fixed=FALSE) {

  if (is_unit_ddi(M1) && !M1.fixed) stop("unit ddiMatrix assumed to be fixed")
  if (q2 == 1L) {  # scalar multiplication
    if (is.matrix(M1)) {
      update <- function(M1, M2x, values.only=FALSE) M2x * M1
    } else if (is_unit_ddi(M1)) {
      template <- expand_unit_ddi(M1)
      update <- function(M1, M2x, values.only=FALSE) {
        x <- M2x * template@x
        if (values.only)
          x
        else {
          out <- template
          attr(out, "x") <- x
          out
        }
      }
    } else {
      update <- function(M1, M2x, values.only=FALSE) {
        x <- M2x * M1@x
        if (values.only) {
          x
        } else {
          out <- M1
          attr(out, "x") <- x
          out
        }
      }
    }
    rm(M1, M2, M1.fixed)
    return(update)
  }
  if (is.matrix(M2) && nrow(M1) == 1L) {
    if (is.matrix(M1)) {
      update <- function(M1, M2x, values.only=FALSE) M1[1L, 1L] * M2x
    } else {
      if (is_unit_ddi(M1))
        update <- function(M1, M2x, values.only=FALSE) M2x
      else
        update <- function(M1, M2x, values.only=FALSE) M1@x * M2x
    }
    rm(M1, M2, M1.fixed)
    return(update)
  }
  base_tcrossprod <- base::tcrossprod
  type <- paste0(substr(class(M1)[1L], 1L, 3L), substr(class(M2)[1L], 1L, 3L))
  switch(type,
    matmat = update <- function(M1, M2x, values.only=FALSE) Cdense_kron(M1, M2x),
    ddinum = {  # result ddi
      expand <- q2 > length(M2)  # scalar M2
      template <- Cdiag(as.numeric(base_tcrossprod(if (expand) rep.int(M2, q2) else M2, if (is_unit_ddi(M1)) rep.int(1, nrow(M1)) else M1@x)))
      attr(template, "x") <- NULL  # save some space
      if (is_unit_ddi(M1)) {
        ones <- rep.int(1, nrow(M1))
        update <- function(M1, M2x, values.only=FALSE) {
          x <- as.numeric(base_tcrossprod(if (expand) rep.int(M2x, q2) else M2x, ones))
          if (values.only)
            x
          else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      } else {
        update <- function(M1, M2x, values.only=FALSE) {
          x <- as.numeric(base_tcrossprod(if (expand) rep.int(M2x, q2) else M2x, M1@x))
          if (values.only)
            x
          else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      }
    },
    ddimat = {  # block diagonal case --> dsCMatrix
      template <- kron_ddimat(M1, M2)
      attr(template, "x") <- NULL
      upper <- which(row(M2) <= col(M2))
      if (is_unit_ddi(M1)) {
        q1 <- nrow(M1)
        update <- function(M1, M2x, values.only=FALSE) {
          x <- rep.int(M2x[upper], q1)
          if (values.only)
            x
          else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      } else {
        update <- function(M1, M2x, values.only=FALSE) {
          x <- as.numeric(base_tcrossprod(M2x[upper], M1@x))
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      }
    },
    matnum = {
      expand <- q2 > length(M2)  # scalar M2
      template <- forceSymmetric(as(kronecker(M1, Cdiag(if (expand) rep.int(M2, q2) else M2)), "CsparseMatrix"), uplo="U")
      attr(template, "x") <- NULL
      M1dsC <- forceSymmetric(as(M1, "CsparseMatrix"), uplo="U")
      w <- which(as.matrix(M1dsC) != 0 & row(M1) <= col(M1))
      d <- diff.default(M1dsC@p)
      d <- d[d > 0L]
      rm(M1dsC)
      update <- function(M1, M2x, values.only=FALSE) {
        x <- Crepgen(M1[w], d, if (expand) rep.int(M2x, q2) else M2x)
        if (values.only) {
          x
        } else {
          out <- template
          attr(out, "x") <- x
          out
        }
      }
    },
    ddidsC = {
      template <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")
      attr(template, "x") <- NULL
      if (is_unit_ddi(M1)) {
        q1 <- nrow(M1)
        update <- function(M1, M2x, values.only=FALSE) {
          x <- rep.int(M2x, q1)
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      } else {
        update <- function(M1, M2x, values.only=FALSE) {
          x <- as.numeric(base_tcrossprod(M2x, M1@x))
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      }
    },
    dsCnum=, dgCnum = {
      expand <- q2 > length(M2)  # scalar M2
      template <- as(kronecker(M1, Cdiag(if (expand) rep.int(M2, q2) else M2)), "CsparseMatrix")
      if (class(M1)[1L] == "dsCMatrix") template <- forceSymmetric(template, uplo="U")
      attr(template, "x") <- NULL
      d <- diff.default(M1@p)
      d <- d[d > 0L]
      update <- function(M1, M2x, values.only=FALSE) {
        x <- Crepgen(M1@x, d, if (expand) rep.int(M2x, q2) else M2x)
        if (values.only) {
          x
        } else {
          out <- template
          attr(out, "x") <- x
          out
        }
      }
    },
    matdsC = {
      template <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")
      upper <- which(row(M1) <= col(M1))
      prod.table <- base_tcrossprod(M1[upper], M2@x)
      ind <- arrayInd(fmatch(template@x, prod.table), dim(prod.table))
      ind1 <- upper[ind[, 1L]]
      ind2 <- ind[, 2L]
      if (!isTRUE(all.equal(M1[ind1] * M2@x[ind2], template@x))) stop("incorrect sparse kronecker template")
      attr(template, "x") <- NULL
      rm(upper, ind, prod.table)
      if (M1.fixed) {
        attr(M2, "x") <- rep.int(1, length(M2@x))
        x0 <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")@x
        rm(ind1)
        update <- function(M1, M2x, values.only=FALSE) {
          x <- x0 * M2x[ind2]
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      } else {
        update <- function(M1, M2x, values.only=FALSE) {
          x <- M1[ind1] * M2x[ind2]
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      }
    },
    dsCmat = {
      template <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")
      upper <- which(row(M2) <= col(M2))
      x1 <- unique(M1@x)
      x1.ind <- fmatch(x1, M1@x)
      prod.table <- base_tcrossprod(x1, M2[upper])
      ind <- arrayInd(fmatch(template@x, prod.table), dim(prod.table))
      ind1 <- x1.ind[ind[, 1L]]
      ind2 <- upper[ind[, 2L]]
      if (!isTRUE(all.equal(M1@x[ind1] * M2[ind2], template@x))) stop("incorrect sparse kronecker template")
      attr(template, "x") <- NULL
      rm(upper, x1, x1.ind, ind, prod.table)
      if (M1.fixed) {
       x0 <- forceSymmetric(as(kronecker(M1, matrix(1, nrow(M2), ncol(M2))), "CsparseMatrix"), uplo="U")@x
       rm(ind1)
       update <- function(M1, M2x, values.only=FALSE) {
         x <- x0 * M2x[ind2]
         if (values.only)
           x
         else {
           out <- template
           attr(out, "x") <- x
           out
         }
       }
      } else {
       update <- function(M1, M2x, values.only=FALSE) {
         x <- M1@x[ind1] * M2x[ind2]
         if (values.only) {
           x
         } else {
           out <- template
           attr(out, "x") <- x
           out
         }
       }
      }
    },
    dsCdsC = {
      template <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")
      x1 <- unique(M1@x)
      x1.ind <- fmatch(x1, M1@x)
      x2 <- unique(M2@x)
      x2.ind <- fmatch(x2, M2@x)
      prod.table <- base_tcrossprod(x1, x2)
      ind <- arrayInd(fmatch(template@x, prod.table), dim(prod.table))
      ind1 <- x1.ind[ind[, 1L]]
      ind2 <- x2.ind[ind[, 2L]]
      if (!isTRUE(all.equal(M1@x[ind1] * M2@x[ind2], template@x))) stop("incorrect sparse kronecker template")
      attr(template, "x") <- NULL
      rm(ind, x1, x1.ind, x2, x2.ind, prod.table)
      if (M1.fixed) {
        attr(M2, "x") <- rep.int(1, length(M2@x))
        x0 <- forceSymmetric(as(kronecker(M1, M2), "CsparseMatrix"), uplo="U")@x
        rm(ind1)
        update <- function(M1, M2x, values.only=FALSE) {
          x <- x0 * M2x[ind2]
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      } else {
        update <- function(M1, M2x, values.only=FALSE) {
          x <- M1@x[ind1] * M2x[ind2]
          if (values.only) {
            x
          } else {
            out <- template
            attr(out, "x") <- x
            out
          }
        }
      }
    },
    stop("unexpected input")
  )
  rm(M1, M2, type)
  update
}

# Kronecker product of ddiMatrix with matrix
# NB Mmat assumed symmetric
kron_ddimat <- function(Mddi, Mmat) {
  upper <- which(row(Mmat) <= col(Mmat))
  rmat <- as.integer(row(Mmat) - 1L)[upper]
  n.ddi <- nrow(Mddi)
  n.mat <- nrow(Mmat)
  i <- rep.int(rmat, n.ddi) + nrow(Mmat) * rep_each(seq_len(n.ddi) - 1L, length(rmat))
  p <- c(0L, cumsum(rep.int(seq_len(n.mat), n.ddi)))
  if (is_unit_ddi(Mddi))
    x <- rep.int(Mmat[upper], n.ddi)
  else
    x <- as.numeric(tcrossprod(Mmat[upper], Mddi@x))
  size <- n.ddi * n.mat
  new("dsCMatrix", i=i, p=p, x=x, uplo="U", Dim=c(size, size))
}
