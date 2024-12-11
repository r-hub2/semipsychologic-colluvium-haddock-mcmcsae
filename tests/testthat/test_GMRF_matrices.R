
context("GMRF matrices")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("GMRF matrices have expected types", {
  n <- 10
  expect_is(D_AR1(n, 0), "dgCMatrix")
  expect_is(Q_AR1(n, 0), "dsCMatrix")
  expect_equal(Q_AR1(n, 1), Q_RW1(n))
  expect_is(D_AR1(n, 0.1), "dgCMatrix")
  expect_is(Q_AR1(n, 0.12), "dsCMatrix")
  expect_is(D_AR1(n, 0.1, w=runif(n - 1)), "dgCMatrix")
  expect_is(Q_AR1(n, 0.12, w=runif(n - 1)), "dsCMatrix")
})

test_that("crossprod of incidence matrix yields precision matrix", {
  expect_equal(crossprod(D_AR1(14, 0.1)), Q_AR1(14, 0.1))
  expect_equal(crossprod(D_season(120, 12)), Q_season(120, 12))
})

x <- seq(0, 1, 0.01)
dat <- data.frame(x, y=x*sin(5*x))

test_that("spline model works", {
  sampler <- create_sampler(y ~ gen(factor = ~ spline(x, knots=15)),
                            linpred="fitted", data=dat)
  sim <- MCMCsim(sampler, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_true(max(abs(dat$y - summ$linpred_[, "Mean"])) < 0.02)
})


n <- 400L
dat <- data.frame(
  x = runif(n),
  f = as.factor(sample(1:10, n, replace = TRUE))
)
v <- rnorm(20, sd=0.5)
dat$y <- 1 + 2*dat$x + v[dat$f] + rnorm(n, sd=0.25)
test_that("custom precision matrix specification works", {
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ custom(f, D=diag(10))),
    data=dat
  )
  Q <- diag(10)
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ custom(f, Q=Q)),
    data=dat
  )
  Q2 <- Q
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ custom(f, Q=Q2)),
    data=dat
  )
  sim <- MCMCsim(sampler, verbose=FALSE)
  summ <- summary(sim)
  expect_true(0.1 < summ$sigma_[, "Mean"] && summ$sigma_[, "Mean"] < 0.4)
})


test_that("function compute_GMRF_matrices works", {
  dat <- data.frame(f = 1:2)
  expect_identical(dim(compute_GMRF_matrices(~ f, data=dat)$D), c(2L, 2L))
  Q0 <- matrix(c(5,2,2,4), 2)
  expect_identical(compute_GMRF_matrices(~ custom(f, Q=Q0), D=FALSE, data=dat)$Q, Q0)
})
