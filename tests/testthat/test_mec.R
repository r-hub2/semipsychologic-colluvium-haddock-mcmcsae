
context("Measurement error component")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

m <- 50
X <- rnorm(m, mean=5, sd=3)  # true covariate values
theta <- 1 + 3*X  # true values
psi <- rgamma(m, shape=4.5, scale=2)
e <- rnorm(m, sd=sqrt(psi))  # sampling error
y <- theta + e  # direct estimates
C <- c(rep(3, 10), rep(0, 40))  # measurement error for first 10 values
W <- X + rnorm(m, sd=sqrt(C))  # covariate subject to measurement error

test_that("measurement error component in linear regression works", {
  sampler <- create_sampler(
    y ~ 1 + mec(~ 0 + W, V=C, name="ME"),
    family=f_gaussian(var.prior=pr_fixed(1), var.vec=~psi),
    linpred="fitted"
  )
  sim <- MCMCsim(sampler, burnin=200, n.iter=700, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  #plot(X, W, xlab="true X", ylab="inferred X")
  #points(X, summ$ME_X[, "Mean"], col="green"); abline(0, 1, col="red")
  #legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
  expect_lt(mean(abs(summ$ME_X[, "Mean"] - X)), mean(abs(W - X)))
})

test_that("measurement error component in linear regression with proper priors works", {
  sampler <- create_sampler(
    y ~ 1 + mec(~ 0 + W, V=C, name="ME", prior=pr_normal(precision=1)),
    family=f_gaussian(var.prior=pr_fixed(1), var.vec=~psi),
    linpred="fitted"
  )
  sim <- MCMCsim(sampler, burnin=200, n.iter=700, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_lt(summ$ME[, "Mean"], 3)
  #plot(X, W, xlab="true X", ylab="inferred X")
  #points(X, summ$ME_X[, "Mean"], col="green"); abline(0, 1, col="red")
  #legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
  expect_lt(mean(abs(summ$ME_X[, "Mean"] - X)), mean(abs(W - X)))
})

# modeled variances
m <- 400
X <- rnorm(m, mean=1, sd=3)  # true covariate values
z <- rnorm(m)
eta <- -1 + 3*X  # true values
y <- rnorm(m, mean=eta, sd = exp(0.5*(1 + z)))
C <- c(rep(3, 100), rep(0, 300))  # measurement error for first 100 values
W <- X + rnorm(m, sd=sqrt(C))  # covariate subject to measurement error

test_that("measurement error component in linear regression with modeled variances works", {
  sampler <- create_sampler(
    y ~ 1 + mec(~ 0 + W, V=C, name="ME"),
    family=f_gaussian(var.prior=pr_fixed(1), var.model = ~ reg(~ z)),
    linpred="fitted"
  )
  expect_true(sampler$modeled.Q)
  sim <- MCMCsim(sampler, burnin=150, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 2*(-1), 0.5*(-1))
  expect_between(summ$vreg1[, "Mean"], 0.5*c(1, 1), 2*c(1, 1))
  expect_between(summ$ME[, "Mean"], 0.5*(3), 2*(3))
  #plot(X, W, xlab="true X", ylab="inferred X")
  #points(X, summ$ME_X[, "Mean"], col="green"); abline(0, 1, col="red")
  #legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
  expect_lt(mean(abs(summ$ME_X[, "Mean"] - X)), mean(abs(W - X)))
})

# generate binomial data
m <- 400
X <- rnorm(m, mean=1, sd=3)  # true covariate values
eta <- -1 + 3*X  # true values
y <- rbinom(m, 100, 1/(1 + exp(-eta)))
C <- c(rep(3, 100), rep(0, 300))  # measurement error for first 100 values
W <- X + rnorm(m, sd=sqrt(C))  # covariate subject to measurement error

test_that("measurement error component in non-gaussian models works", {
  expect_error(
    sampler <- create_sampler(
      y ~ 1 + mec(~ 0 + W, V=C, name="ME"),
      family=f_gaussian(var.model = ~ reg(~ 1) + mec(~ 0 + W, V=C)),
      linpred="fitted"
    ), "mean model"
  )
  sampler <- create_sampler(
    y ~ reg(~1, prior=pr_normal(mean=-1, precision=100)) +
      mec(~ 0 + W, V=C, name="ME", debug=FALSE),
    family=f_binomial(n.trial = 100),
    linpred="fitted",
    control=sampler_control(block=FALSE)
  )
  sim <- MCMCsim(sampler, burnin=150, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 2*(-1), 0.5*(-1))
  expect_between(summ$ME[, "Mean"], 0.5*(3), 2*(3))
  #plot(X, W, xlab="true X", ylab="inferred X")
  #points(X, summ$ME_X[, "Mean"], col="green"); abline(0, 1, col="red")
  #legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
  expect_lt(mean(abs(summ$ME_X[, "Mean"] - X)), mean(abs(W - X)))
  
  # fully blocked Gibbs sampler:
  sampler <- create_sampler(
    y ~ reg(~1, prior=pr_normal(mean=-1, precision=100)) +
      mec(~ 0 + W, V=C, name="ME", debug=FALSE),
    family=f_binomial(n.trial = 100),
    linpred="fitted",
    control=sampler_control(block=TRUE)
  )
  sim <- MCMCsim(sampler, burnin=150, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 2*(-1), 0.5*(-1))
  expect_between(summ$ME[, "Mean"], 0.5*(3), 2*(3))
  #plot(X, W, xlab="true X", ylab="inferred X")
  #points(X, summ$ME_X[, "Mean"], col="green"); abline(0, 1, col="red")
  #legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
  expect_lt(mean(abs(summ$ME_X[, "Mean"] - X)), mean(abs(W - X)))
})
