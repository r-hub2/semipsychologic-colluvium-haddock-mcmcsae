
context("Truncated multivariate normal sampling")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

mu <- c(-1, 1)
V <- matrix(c(1, 0.99, 0.99, 1), 2, 2)
Q <- solve(V)

test_that("direct unrestricted MVN sampling method works", {
  sampler <- create_TMVN_sampler(Q=Q, mu=mu)
  expect_true(sampler$method$use.cholV)
  sim <- MCMCsim(sampler, burnin=0, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$x[, "Mean"], mu - 0.1, mu + 0.1)
  expect_between(cov(as.matrix(sim$x)), 0.8 * V, 1.2 * V)
  sampler <- create_TMVN_sampler(Q=Q, mu=mu, method=m_direct(use.cholV=FALSE))
  expect_false(sampler$method$use.cholV)
  sim <- MCMCsim(sampler, burnin=0, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$x[, "Mean"], mu - 0.1, mu + 0.1)
  expect_between(cov(as.matrix(sim$x)), 0.8 * V, 1.2 * V)
})

test_that("equality-restricted MVN sampler works", {
  R <- set_constraints(R = cbind(c(1, 1)))
  # conditioning by kriging
  sampler <- create_TMVN_sampler(Q=Q, mu=mu, constraints=R)
  sim <- MCMCsim(sampler, burnin=0, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_lt(mean(abs(c(as.matrix(sim$x) %*% R$R))), 1e-10)
  R$check(summ$x[, "Mean"])
  # sample from a reduced system
  sampler <- create_TMVN_sampler(Q=Q, mu=mu, constraints=R, reduce=TRUE)
  sim <- MCMCsim(sampler, burnin=0, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_lt(mean(abs(c(as.matrix(sim$x) %*% R$R))), 1e-10)
  R$check(summ$x[, "q0.5"])
})

# 5d example, 1 equality, one inequality
C <- set_constraints(
  R = cbind(c(1,-1,-1,-1,1)),  # equalities R'x = r
  r = 0,
  S = cbind(c(0,0,0,1,-1)),    # inequalities S'x >= s
  s = 1.7
)
mu0 <- c(6.9, 4.2, 1.0, 4.9, 4.3)
sd0 <- c(0.3, 0.2, 0.1, 0.3, 0.3)
Q0 <- diag(1/sd0^2)  # precision

# 'exact' answer
answer <- c(6.974, 4.167, 0.9918, 5.507, 3.693)
answer.sd <- c(0.191, 0.172, 0.097, 0.218, 0.218)


test_that("HMC TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C)
  expect_equal(sampler$method$name, "HMC")
  sim <- MCMCsim(sampler, burnin=1000, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})

test_that("HMC TMVN method works after projection on equality constraint surface", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C, reduce=TRUE)
  expect_equal(sampler$method$name, "HMC")
  sim <- MCMCsim(sampler, burnin=1000, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})

test_that("HMCZigZag TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C,
                                 method=m_HMCZigZag(scale=5/sqrt(diag(Q0))))
  expect_equal(sampler$method$name, "HMCZigZag")
  sim <- MCMCsim(sampler, burnin=40, n.iter=80, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0, tolerance=0.2)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.2)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.2)
  # also test with projection on equality constraint surface
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C,
    reduce=TRUE, method=m_HMCZigZag(scale=2),
  )
  sim <- MCMCsim(sampler, burnin=40, n.iter=80, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0, tolerance=0.2)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.2)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.2)
})

test_that("Gibbs TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C, method="Gibbs")
  expect_equal(sampler$method$name, "Gibbs")
  sim <- MCMCsim(sampler, burnin=500, n.iter=3000, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
  # does it work without equalities?
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0,
    constraints=set_constraints(S=C$S, s=C$s), method="Gibbs"
  )
  expect_false(sampler$eq)
  sim <- MCMCsim(sampler, burnin=500, n.iter=2000, verbose=FALSE)
  summ <- summary(sim)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
})

test_that("slice-Gibbs TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C, method=m_Gibbs(slice=TRUE))
  expect_true(sampler$method$slice)
  sim <- MCMCsim(sampler, burnin=500, n.iter=3000, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.05)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.05)
  # does it work without equalities?
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0,
    constraints=set_constraints(S=C$S, s=C$s), method=m_Gibbs(slice=TRUE)
  )
  expect_false(sampler$eq)
  sim <- MCMCsim(sampler, burnin=500, n.iter=2000, verbose=FALSE)
  summ <- summary(sim)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
})

test_that("Soft TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, constraints=C, method="softTMVN")
  expect_equal(sampler$method$name, "softTMVN")
  sim <- MCMCsim(sampler, burnin=500, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(C$R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
  # does it work without equalities?
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0,
    constraints=set_constraints(S=C$S, s=C$s), method="softTMVN"
  )
  sim <- MCMCsim(sampler, burnin=500, n.iter=2000, verbose=FALSE)
  summ <- summary(sim)
  expect_true(crossprod_mv(C$S, summ$x[, "Mean"]) >= C$s)
})

test_that("inequality-constrained linear regression works", {
  n <- 200
  dat <- data.frame(
    x = rnorm(n),
    z = runif(n)
  )
  dat$y <- 0.1 + 0.05 * dat$x + 0.01 * dat$z + rnorm(n, sd=0.1)
  #lm(y ~ x + z, data=dat)
  sampler <- create_sampler(
    y ~ reg(~ 1 + x + z, constraints=set_constraints(S=diag(3)), name="beta"),
    data=dat
  )
  sim <- MCMCsim(sampler, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all(summ$beta[, "Mean"] >= 0))
})

test_that("inequality-constrained random effects work", {
  # two time series where the first has to be larger than the second
  nT <- 25L  # length of time series
  nx <- 2L   # number of variables
  df <- data.frame(time=rep(1:nT, nx), idx=as.factor(rep(1:nx, each=nT)))
  Q0 <- rep(2, nT*nx)
  dat <- generate_data(
    ~ reg(~ 0 + idx, prior=pr_normal(mean=c(10, 9), precision=c(1, 1))) +
      gen(~ 0 + idx, var="diagonal", factor = ~ RW1(time), prior=pr_invchisq(df=5, scale=0.5)),
    data=df, family = f_gaussian(var.prior=1, var.vec = ~ I(1/Q0))
  )
  df$y <- dat$y
  sampler <- create_sampler(
    y ~ 0 + idx + gen(~ 0 + idx, var="diagonal", factor=~RW1(time)),
    family=f_gaussian(var.prior=1, var.vec= ~ I(1/Q0)), data=df, linpred="fitted"
  )
  sim <- MCMCsim(sampler, n.iter=500, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.9 * dat$pars$reg1, 1.1 * dat$pars$reg1)
  df$Q0 <- Q0
  sampler <- create_sampler(
    y ~ gen(~ 0 + idx, var="diagonal", factor = ~RW1(time),
      constraints0 = set_constraints(S=matrix(c(1, -1), 2, 1)),
      GMRFconstr=FALSE, PX=FALSE
    ),
    family=f_gaussian(var.prior=1, var.vec = ~ I(1/Q0)), data=df, linpred="fitted"
  )
  expect_identical(dim(sampler$mod$gen1$S), c(nx*nT, nT))
  start <- function() {
    # use feasible start values
    list(gen1 = rep(c(9.75, 8.75), nT) + 0.5*runif(2*nT))
  }
  sim <- MCMCsim(sampler, n.iter=600, n.chain=2, verbose=FALSE, start=start)
  summ <- summary(sim)
  expect_true(all(summ$linpred_[1:nT, "Mean"] + 0.1 >= summ$linpred_[(nT+1):(2*nT), "Mean"]))
})
