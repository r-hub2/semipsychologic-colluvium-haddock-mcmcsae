
context("Conjugate gradient algorithm")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 100
A <- crossprod(matrix(rnorm(n*n), n, n)) + 10*diag(n)
b <- rnorm(n)

create_env <- function() {
  A_times <- function(x) A %m*v% x
  #M_solve <- identity
  idA <- 1/diag(A)
  M_solve <- function(x) x * idA 
  environment()
}

test_that("conjugate gradient solver for dense symmetric matrix works", {
  sol <- solve.default(A, b)
  approx <- CG(b, env=create_env(), e=1e-6, verbose=FALSE)
  expect_equal(approx, sol, tolerance=1e-2)
  #plot(sol, approx); abline(0, 1, col="red")
})

n <- 1000
A <- rsparsematrix(n, n, 0.02, symmetric=TRUE) + 10*Diagonal(n) + Q_RW2(n)
b <- rnorm(n)

test_that("conjugate gradient solver for sparse symmetric matrix works", {
  sol <- build_chol(A)$solve(b)
  approx <- CG(b, env=create_env(), e=1e-6, verbose=FALSE)
  expect_equal(approx, sol, tolerance=1e-2)
  #plot(sol, approx); abline(0, 1, col="red")
  #microbenchmark(
  #  sol <- build_chol(A)$solve(b),
  #  sol <- solve(A, b),
  #  approx <- CG(b, env=create_env(), e=1e-6)
  #)
})
