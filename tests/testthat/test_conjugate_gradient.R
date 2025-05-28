
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


# generate some data
n <- 1000
df <- data.frame(
  x1 = runif(n),
  x2 = rnorm(n)
)
dat <- generate_data(
  ~ reg(~ x1 + x2, prior=pr_normal(mean=c(1,2,3), precision=100), name="b"), data=df,
  family = f_gaussian(var.prior = pr_fixed(1e-6))
)
df$y <- dat$y

test_that("conjugate gradient sampler for a simple regression model works", {
  sampler <- create_sampler(
    y ~ reg(~ x1 + x2), data=df
  )
  expect_length(sampler$block, 0L)  # by default mc_block not used in case of a single block
  sampler <- create_sampler(
    y ~ reg(~ x1 + x2, name="b"), data=df,
    control=list(CG=TRUE)
  )
  expect_length(sampler$block, 1L)  # for CG need mc_block
  sim <- MCMCsim(sampler, verbose=FALSE, n.chain=2L, n.iter=700L)
  summ <- summary(sim)
  expect_between(summ$b[, "Mean"], 0.5*dat$pars$b, 2*dat$pars$b)
  expect_between(summ$sigma_[, "Mean"], 0.9 * 1e-3, 1.1 * 1e-3)
  expect_gt(compute_DIC(sim)["p_DIC"], 0)
})
