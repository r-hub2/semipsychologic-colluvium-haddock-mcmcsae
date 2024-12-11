
context("Matrix algebra")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("checking for redundant columns works", {
  M <- matrix(rnorm(5*2), 5, 2)
  expect_null(detect_redundancy(crossprod(M)))
  expect_null(detect_redundancy(M, method="qr"))
  expect_identical(remove_redundancy(M), M)
  Mext <- M[, c(1, 1, 1, 2)]
  expect_identical(remove_redundancy(Mext), M)
})

test_that("is_zero_matrix works", {
  expect_false(is_zero_matrix(matrix(0.1, 2, 1)))
  expect_true(is_zero_matrix(matrix(0, 2, 1)))
  expect_false(is_zero_matrix(Diagonal(3)))
  expect_true(is_zero_matrix(0*Diagonal(3)))
  expect_false(is_zero_matrix(as(Diagonal(1), "tabMatrix")))
  expect_true(is_zero_matrix(0*as(Diagonal(1), "tabMatrix")))
})

test_that("zero_col works", {
  expect_true(!any(zero_col(Diagonal(3))))
  expect_identical(zero_col(Diagonal(x=c(1, 0, 2))), c(FALSE, TRUE, FALSE))
  M <- as(matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  expect_identical(zero_col(M), c(FALSE, FALSE, TRUE))
  expect_identical(zero_col(Ctab2dgC(M)), zero_col(Ctab2mat(M)))
  M <- as(matrix(c(-2,0,0,0,1,0,0,0,0), 3, 3), "tabMatrix")
  expect_identical(zero_col(M), c(FALSE, FALSE, TRUE))
  expect_identical(zero_col(Ctab2dgC(M)), zero_col(Ctab2mat(M)))
  M <- as(matrix(c(1,0,1,0,0,0,0,0,0), 3, 3), "tabMatrix")
  expect_identical(zero_col(M), c(FALSE, TRUE, TRUE))
  expect_identical(zero_col(Ctab2dgC(M)), zero_col(Ctab2mat(M)))
})

test_that("inverseSPD works", {
  n <- 4
  M <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)
  expect_equal(solve(M), inverseSPD(M))
})

test_that("inverse method of cholesky object works", {
  n <- 20L
  M <- crossprod(matrix(rnorm(n*n), n, n)) + 2*diag(n)
  MdsC <- as(M, "CsparseMatrix")
  cholM <- build_chol(MdsC)
  MdsCinv <- cholM$inverse()
  expect_is(MdsCinv, "dsCMatrix")
  expect_equal(as.matrix(MdsCinv), inverseSPD(M))
})

test_that("dotprodC works", {
  x <- rnorm(10)
  y <- runif(10)
  expect_equal(sum(x*y), dotprodC(x, y))
})

test_that("add_diagC works", {
  n <- 7
  M <- matrix(rnorm(n*n), n, n)
  d <- rnorm(n)
  Md <- add_diagC(M, d)
  expect_equal(Md, M + diag(d))
  expect_equal(diag(M), diag(Md) - d)
})

test_that("matrix-vector products work", {
  n <- 4
  M <- matrix(rnorm(n*n), n)
  x <- rnorm(n)
  expect_equal(as.numeric(M %*% x), M %m*v% x)
  M <- as(M, "CsparseMatrix")
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- crossprod(M)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- Diagonal(n)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- Diagonal(x=rnorm(n))
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- as(matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- as(rnorm(3)*matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
  M <- as(matrix(c(1,0,0), 3, 1), "tabMatrix")
  x <- 2
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
  M <- as(matrix(c(1,0,0,0,2,0), 3, 2), "tabMatrix")
  x <- rnorm(2)
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
})

test_that("matrix-vector crossproducts work", {
  n <- 4
  M <- matrix(rnorm(n*n), n)
  x <- rnorm(n)
  expect_equal(as.numeric(crossprod(M, x)), crossprod_mv(M, x))
  M <- as(M, "CsparseMatrix")
  expect_equal(as.vector(crossprod(M, x)), crossprod_mv(M, x))
  M <- crossprod(M)
  expect_equal(as.vector(crossprod(M, x)), crossprod_mv(M, x))
  M <- Diagonal(n)
  expect_equal(x, crossprod_mv(M, x))
  M <- Diagonal(x=rnorm(n))
  expect_equal(M@x * x, crossprod_mv(M, x))
  M <- as(matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(as.vector(t(M) %*% x), crossprod_mv(M, x))
  M <- as(rnorm(3)*matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
  x <- rnorm(2)
  expect_error(crossprod_mv(M, x))  # incompatible dimensions
  M <- as(matrix(c(1,0,0), 3, 1), "tabMatrix")
  x <- rnorm(3)
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
  M <- as(matrix(c(1,0,0,0,2,0), 3, 2), "tabMatrix")
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
})

test_that("diagonal-dgC product works", {
  nr <- 10; nc <- 25
  Q <- Diagonal(x=rnorm(nr))
  X <- rsparsematrix(nr, nc, 0.01)
  expect_equal(Q %*% X, Cdiag_sparse_prod(Q@x, X))
})

test_that("tab to dgC conversion works", {
  m1 <- matrix(c(1,0,1,0,1,0,0,0,0), 3, 3)
  M1 <- as(m1, "tabMatrix")
  expect_false(M1@reduced)
  expect_false(M1@num)
  expect_equal(Ctab2dgC(M1), as(as(m1, "CsparseMatrix"), "generalMatrix"))
  m1 <- matrix(c(2,0,0), 3, 1)
  M1 <- as(m1, "tabMatrix")
  expect_false(M1@reduced)
  expect_true(M1@num)
  expect_equal(Ctab2dgC(M1), as(as(m1, "CsparseMatrix"), "generalMatrix"))
  M2 <- aggrMatrix(sample(1:7, 11, replace=TRUE))
  expect_equal(Ctab2dgC(M2), as(as(M2, "CsparseMatrix"), "generalMatrix"))
  M2 <- aggrMatrix(sample(1:7, 11, replace=TRUE), w=runif(11))
  expect_equivalent(as(Ctab2dgC(M2), "matrix"), Ctab2mat(M2))
  expect_error(as(matrix(c(1,1,1,0), 2, 2), "tabMatrix"))
})

test_that("tab <-> matrix conversion works", {
  m1 <- matrix(0, 2, 3)
  expect_equal(as(as(m1, "tabMatrix"), "matrix"), m1)
  m1 <- matrix(c(0, 0, 0, 1.1), 2, 2)
  expect_equal(as(as(m1, "tabMatrix"), "matrix"), m1)
})

test_that("tabMatrix row selection works", {
  m <- matrix(c(1,0,1,0,1,0,0,0,0), 3, 3)
  M <- as(m, "tabMatrix")
  expect_identical(M[1, ], m[1, ])
  expect_identical(Ctab2mat(M[1, , drop=FALSE]), m[1, , drop=FALSE])
  expect_identical(Ctab2mat(M[c(1, 3), ]), m[c(1, 3), ])
  expect_identical(Ctab2mat(M[-2, ]), m[-2, ])
  expect_identical(M[c(-2, -3), ], m[c(-2, -3), ])
  expect_identical(Ctab2mat(M[, c(-2, -3), drop=FALSE]), m[, c(-2, -3), drop=FALSE])
  rownames(M) <- paste0("r", seq_len(ncol(M)))
  expect_identical(Ctab2mat(M[c("r2"), , drop=FALSE]), m[2, , drop=FALSE])
})

test_that("tabMatrix column selection works", {
  m <- matrix(c(1,0,1,0,1,0,0,0,0), 3, 3)
  M <- as(m, "tabMatrix")
  expect_identical(M[, 1], m[, 1])
  expect_identical(Ctab2mat(M[, 1, drop=FALSE]), m[, 1, drop=FALSE])
  expect_identical(Ctab2mat(M[, c(1, 3)]), m[, c(1, 3)])
  expect_identical(Ctab2mat(M[, -2]), m[, -2])
  expect_identical(M[, c(-2, -3)], m[, c(-2, -3)])
  expect_identical(Ctab2mat(M[, c(-2, -3), drop=FALSE]), m[, c(-2, -3), drop=FALSE])
  colnames(M) <- paste0("c", seq_len(ncol(M)))
  expect_identical(Ctab2mat(M[, c("c2"), drop=FALSE]), m[, 2, drop=FALSE])
})

test_that("diag works for tabMatrix", {
  Mm <- matrix(c(1,0,1,0,2.3,0,0,0,0), 3, 3)
  M <- as(Mm, "tabMatrix")
  expect_equal(diag(M), diag(Mm))
})

test_that("crossprod_sym works", {
  n <- 25L
  q <- runif(n)
  Q <- Diagonal(x=q)
  Qsym <- crossprod(rsparsematrix(n, n, density=0.1))
  Qmat <- as(Qsym, "matrix")
  M <- matrix(rnorm(n^2), n, n)
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(crossprod(M, Qsym %*% M), crossprod_sym(M, Qsym))
  expect_equal(crossprod(M, Qmat %*% M), crossprod_sym(M, Qmat))
  M <- Diagonal(n)
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(crossprod(M, Qsym %*% M), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- Diagonal(x=rnorm(n))
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- aggrMatrix(sample(1:n, n))
  expect_true(tab_isPermutation(M))
  expect_equal(as(crossprod(M, Q %*% M), "diagonalMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- M[, -c(10, 12, 21)]
  expect_false(tab_isPermutation(M))
  expect_equal(as(crossprod(M, Q %*% M), "diagonalMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  M <- economizeMatrix(M, sparse=TRUE)  # this adds slot isPermutation again
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- aggrMatrix(sample(1:n, n), w=runif(n))
  expect_false(tab_isPermutation(M))
  expect_equal(as(crossprod(M, Q %*% M), "diagonalMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  attr(M, "isPermutation") <- FALSE
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  m <- 12L
  M <- rsparsematrix(n, m, density=0.1)
  expect_equal(as(crossprod(M, Q %*% M), "symmetricMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
})

test_that("crossprod_sym2 works", {
  n <- 11L; m <- 5L
  M1 <- matrix(runif(n*m), n, m)
  expect_equal(crossprod_sym2(M1), crossprod(M1))
  M2 <- diag(runif(n)) %*% M1
  expect_equal(crossprod_sym2(M1, M2), crossprod(M1, M2))
  M1 <- as(as(M1, "CsparseMatrix"), "generalMatrix")
  M2 <- as(as(M2, "CsparseMatrix"), "generalMatrix")
  expect_equal(crossprod_sym2(M1), crossprod(M1))
  expect_equal(crossprod_sym2(M1, M2), as(crossprod(M1, M2), "symmetricMatrix"))
})

test_that("(re)defined S4 methods work", {
  n <- 10L; m <- 5L
  M <- matrix(rnorm(n*m), n, m)
  expect_equal(Diagonal(n) %*% M, M)
  expect_equal(Diagonal(x=1:n) %*% M, (1:n) * M)
  Ms <- rsparsematrix(n, n, density=0.2)
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
  Ms <- crossprod(Ms)
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
  Ms <- aggrMatrix(sample(1:n, n))
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
  M <- Diagonal(n)
  x <- matrix(rnorm(n*m), nrow=n)
  expect_equal(solve(M, x), x)
  expect_equal(solve(M, x[, 1L]), x[, 1L])
  dM <- 0.1 + runif(n)
  M <- Diagonal(x=dM)
  expect_equal(solve(M, x), diag(1/dM) %*% x)
  expect_equal(solve(M, x[, 1L]), x[, 1L] / dM)
})

test_that("sparse symmetric weighted crossprod template works", {
  n <- 1000
  p <- 100
  X <- rsparsematrix(n, p, density = 1e-2)
  W <- 0.1 + runif(n)
  XWX_updater <- sparse_crossprod_sym_template(X)
  expect_equal(crossprod_sym(X, W), XWX_updater(W))
})

test_that("row and column selection of tabMatrix works", {
  n <- 100
  m <- 50
  M <- aggrMatrix(sample(1:m, n, replace=TRUE))
  expect_equal(M[11, ], as.matrix(M)[11, ])
  expect_equal(as.matrix(M[11, , drop=FALSE]), as.matrix(M)[11, , drop=FALSE])
  expect_equal(as.matrix(M[1:10, ]), as.matrix(M)[1:10, ])
  expect_equal(M[, 11], as.matrix(M)[, 11])
  expect_equal(as.matrix(M[, 11, drop=FALSE]), as.matrix(M)[, 11, drop=FALSE])
  expect_equal(as.matrix(M[, c(11, 15, 20)]), as.matrix(M)[, c(11, 15, 20)])
  expect_equal(as.matrix(M[, -2]), as.matrix(M)[, -2])
  M <- aggrMatrix(sample(1:m, n, replace=TRUE), w=runif(n))
  expect_equal(M[11, ], as.matrix(M)[11, ])
  expect_equal(as.matrix(M[11, , drop=FALSE]), as.matrix(M)[11, , drop=FALSE])
  expect_equal(as.matrix(M[1:10, ]), as.matrix(M)[1:10, ])
  expect_equal(M[, 11], as.matrix(M)[, 11])
  expect_equal(as.matrix(M[, 11, drop=FALSE]), as.matrix(M)[, 11, drop=FALSE])
  expect_equal(as.matrix(M[, c(11, 15, 20)]), as.matrix(M)[, c(11, 15, 20)])
  expect_equal(as.matrix(M[, -2]), as.matrix(M)[, -2])
})

test_that("rbinding tabMatrices works", {
  n <- 100
  m <- 50
  M1 <- aggrMatrix(factor(sample(1:m, n, replace=TRUE), levels=1:50))
  M2 <- aggrMatrix(factor(sample(1:m, 2*n, replace=TRUE), levels=1:50))
  expect_identical(rbind(M1, M2), as(rbind(Ctab2dgC(M1), Ctab2dgC(M2)), "tabMatrix"))
  M2 <- aggrMatrix(factor(sample(1:m, 2*n, replace=TRUE), levels=1:50), w=runif(2*n))
  expect_identical(rbind(M1, M2), as(rbind(Ctab2dgC(M1), Ctab2dgC(M2)), "tabMatrix"))
  M1 <- M1[, 1:10]
  M2 <- M2[, 1:10]
  expect_identical(rbind(M1, M2), as(rbind(Ctab2dgC(M1), Ctab2dgC(M2)), "tabMatrix"))
  colnames(M1) <- paste0("c", 1:ncol(M1))
  expect_identical(colnames(rbind(M1, M2)), colnames(M1))
  rownames(M2) <- paste0("r", 1:nrow(M2))
  expect_identical(rownames(rbind(M1, M2)), c(rep.int("", nrow(M1)), rownames(M2)))
})

test_that("economizeMatrix works", {
  M <- matrix(1L, 2, 3)  # integer matrix
  expect_is(economizeMatrix(M, sparse=TRUE), "dgCMatrix")
  n <- 25
  M <- Diagonal(x=rnorm(n))
  colnames(M) <- paste0("c", 1:n)
  expect_equal(colnames(economizeMatrix(M)), NULL)
  expect_equal(colnames(economizeMatrix(M, strip.names=FALSE)), colnames(M))
  expect_is(economizeMatrix(M, vec.diag=TRUE, strip.names=FALSE), "numeric")
  expect_equal(names(economizeMatrix(M, vec.diag=TRUE, strip.names=FALSE)), colnames(M))
  M <- diag(10)
  colnames(M) <- paste("c", 1:10)
  expect_is(economizeMatrix(M, sparse=TRUE), "ddiMatrix")
  expect_equal(colnames(economizeMatrix(M, sparse=TRUE, strip.names=FALSE)), colnames(M))
})
