
context("Cholesky solves")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("Cholesky for ddiMatrix works", {
  n <- 100L
  M <- Diagonal(x=runif(n, 0.1, 1))
  y0 <- rnorm(n)
  y1 <- matrix(rnorm(3*n), n, 3)
  y2 <- rsparsematrix(n, 10, 0.1)
  ch <- build_chol(M)
  M <- Diagonal(x = 0.1 + 10*rnorm(n)^2)
  ch$update(M)
  expect_equal(ch$solve(y0), y0/M@x)
  expect_equal(ch$solve(y1), y1/M@x)
  expect_equal(ch$solve(y2), Diagonal(x=1/M@x) %*% y2)
  expect_equal(ch$solve(y0, system="Lt"), y0/sqrt(M@x))
  expect_equal(ch$solve(y0, system="Lt"), y0/sqrt(M@x))
  expect_equal(ch$solve(y0, system="L"), y0/sqrt(M@x))
  expect_equal(ch$solve(y1, system="Lt"), y1/sqrt(M@x))
  expect_equal(ch$solve(y2, system="Lt"), Diagonal(x=1/sqrt(M@x)) %*% y2)
})

test_that("Cholesky for matrix works", {
  n <- 12
  M <- crossprod(matrix(rnorm(n^2), n))
  diag(M) <- diag(M) + runif(n, 0.1, 1)
  ch <- build_chol(M)
  M <- crossprod(matrix(rnorm(n^2), n))
  diag(M) <- diag(M) + runif(n, 0.1, 1)
  ch$update(M)
  y0 <- rnorm(n)
  y1 <- matrix(rnorm(3*n), n, 3)
  expect_equal(ch$solve(y0), solve(M, y0))
  expect_equal(ch$solve(y1), solve(M, y1))
  expect_equal(ch$solve(y0, system="Lt"), solve(chol(M), y0))
  expect_equal(ch$solve(y0, system="L"), solve(t(chol(M)), y0))
  expect_equal(ch$solve(y1, system="Lt"), solve(chol(M), y1))
  x <- rnorm(n)
  expect_equal(ch$Ltimes(x), ch$cholM %m*v% x)
})

test_that("Cholesky for dsCMatrix works", {
  n <- 100
  M <- crossprod(rsparsematrix(n, n, density=0.01)) + Diagonal(n)
  ch <- build_chol(M)
  # in-place update only works for same sparsity pattern!!
  M <- M + Diagonal(n)
  ch$update(M)
  cholM <- Cholesky_dsC(M, perm=FALSE, super=NA)
  y0 <- rnorm(n)
  y1 <- matrix(rnorm(3*n), n, 3)
  y2 <- rsparsematrix(n, 10, 0.1)
  expect_equal(ch$solve(y0), as.vector(solve(cholM, y0)))
  expect_equal(ch$solve(y1), array(solve(cholM, y1)@x, dim(y1)))
  expect_equal(ch$solve(y2), solve(cholM, y2))
  expect_equal(ch$solve(y0, system="Lt"), as.vector(solve(cholM, y0, system="Lt")))
  expect_equal(ch$solve(y1, system="Lt"), array(solve(cholM, y1, system="Lt")@x, dim(y1)))
  expect_equal(ch$solve(y2, system="Lt"), solve(cholM, y2, system="Lt"))
  expect_equal(ch$solve(y0, system="L"), as.vector(solve(cholM, y0, system="L")))
  expect_equal(ch$solve(y1, system="L"), array(solve(cholM, y1, system="L")@x, dim(y1)))
  expect_equal(ch$solve(y2, system="L"), solve(cholM, y2, system="L"))
  # with permutation
  ch <- build_chol(M, control=chol_control(perm=TRUE))
  M <- M + Diagonal(x=runif(n))
  ch$update(M)
  cholM <- Cholesky_dsC(M, perm=TRUE, super=NA)
  expect_equal(ch$solve(y0), as.vector(solve(cholM, y0)))
  expect_equal(ch$solve(y1), array(solve(cholM, y1)@x, dim(y1)))
  expect_equal(ch$solve(y2), solve(cholM, y2))
  expect_equal(ch$solve(y0, system="Lt"), as.vector(solve(cholM, solve(cholM, y0, system="Lt"), system="Pt")))
  expect_equal(ch$solve(y1, system="Lt"), array(solve(cholM, solve(cholM, y1, system="Lt"), system="Pt")@x, dim(y1)))
  expect_equal(ch$solve(y2, system="Lt"), solve(cholM, solve(cholM, y2, system="Lt"), system="Pt"))
  expect_equal(ch$solve(y0, system="L"), as.vector(solve(cholM, solve(cholM, y0, system="P"), system="L")))
  expect_equal(ch$solve(y1, system="L"), array(solve(cholM, solve(cholM, y1, system="P"), system="L")@x, dim(y1)))
  expect_equal(ch$solve(y2, system="L"), solve(cholM, solve(cholM, y2, system="P"), system="L"))
})

test_that("Determinants are computed correctly", {
  n <- 30
  M <- Diagonal(x=runif(n, 0.1, 1))
  ch <- build_chol(M)
  expect_equal(ch$logdet(), 0.5*c(determinant(M, logarithm = TRUE)$modulus))
  M <- as.matrix(crossprod(matrix(rnorm(n^2), n, n)) + 10*CdiagU(n))
  ch <- build_chol(M)
  expect_equal(ch$logdet(), 0.5*c(determinant(M, logarithm = TRUE)$modulus))
  M <- rsparsematrix(n, n, 0.1, symmetric=TRUE) + 10*Diagonal(n)
  ch <- build_chol(M)
  expect_equal(ch$logdet(), 0.5*c(determinant(M, logarithm = TRUE)$modulus))
  # determinant templates
  n <- 999  # for n < 1000 an approach based on eigen decomposition is used
  QA <- Q_RW1(n)
  Q <- 0.5 * QA + 0.5 * CdiagU(n)
  detchol <- make_det(Q)
  L <- runif(1L)
  w1 <- 2*L; w2 <- 1 - 2*L
  expect_equal(detchol(w1, w2), c(determinant(L * Q_RW1(n) + (1 - L) * CdiagU(n), logarithm=TRUE)$modulus))
  n <- 1001  # for n > 1000 (sparse) Cholesky updates are used
  QA <- Q_RW1(n)
  Q <- 0.5 * QA + 0.5 * CdiagU(n)
  detchol <- make_det(Q)
  L <- runif(1L)
  w1 <- 2*L; w2 <- 1 - 2*L
  expect_equal(detchol(w1, w2), c(determinant(L * Q_RW1(n) + (1 - L) * CdiagU(n), logarithm=TRUE)$modulus))
})
