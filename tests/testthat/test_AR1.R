
context("Models with AR1 component")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 500L
dat <- data.frame(
  x = runif(n),
  t = 1:n
)
gd <- generate_data(
  ~ reg(~ 1 + x, prior=pr_fixed(c(1,2)), name="beta") +
    gen(factor = ~ AR1(t, -0.9), name="v"),
  sigma.mod = pr_fixed(1), data=dat
)
dat$y <- gd$y
#ts.plot(dat$y)

test_that("data generation for AR1 model works", {
  expect_equal(gd$pars$beta, c(`(Intercept)`=1, x=2))
  expect_equal(unname(gd$pars$sigma_), 1)
  expect_between(cor(gd$pars$v[-1], gd$pars$v[-n]), -1, -0.5)
})

test_that("AR1 parameter can be inferred well with blocked Gibbs sampler", {
  expect_error(Q_AR1(10, 2), "between")
  expect_error(D_AR1(7, -1.1), "between")
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ AR1(t), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=TRUE)
  )
  sim <- MCMCsim(sampler, n.iter=500, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5*c(1,2), 2*c(1,2))
  expect_between(summ$v_sigma[, "Mean"], 0.5*gd$pars$v_sigma, 2*gd$pars$v_sigma)
  expect_between(summ$v_AR1[, "Mean"], -1, -0.3)
})

test_that("AR1 parameter can be inferred well without Gibbs blocking", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ AR1(t), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=FALSE)
  )
  sim <- MCMCsim(sampler, n.iter=750, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5*c(1,2), 2*c(1,2))
  expect_between(summ$v_sigma[, "Mean"], 0.5*gd$pars$v_sigma, 2*gd$pars$v_sigma)
  expect_between(summ$v_AR1[, "Mean"], -1, -0.3)
})

test_that("uniform prior for AR1 parameter works", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + 
        gen(factor = ~ AR1(t, pr_unif(-0.7, 1)), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=TRUE)
  )
  expect_false(sampler$mod$v$AR1sampler$tnprior)
  expect_equal(sampler$mod$v$info$extra[[1]]$prior$min, -0.7)
  sim <- MCMCsim(sampler, n.iter=600, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5*c(1,2), 2*c(1,2))
  expect_gte(summ$v_AR1[, "Mean"], -0.7)
})

test_that("truncated normal prior for AR1 parameter works", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + 
      gen(factor = ~ AR1(t, pr_truncnormal(mean=0, precision=1e6, lower=-1, upper=1)), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=TRUE)
  )
  expect_true(sampler$mod$v$AR1sampler$tnprior)
  expect_equal(sampler$mod$v$info$extra[[1]]$prior$lower, -1)
  sim <- MCMCsim(sampler, n.iter=500, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.25*c(1,2), 4*c(1,2))
  expect_between(summ$v_AR1[, "Mean"], -0.2, 0.2)
})

test_that("AR1 parameter can be inferred under non-normal random effects' prior", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ AR1(t), priorA=pr_exp(), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=FALSE)
  )
  sim <- MCMCsim(sampler, n.iter=500, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5*c(1,2), 2*c(1,2))
  expect_between(summ$v_sigma[, "Mean"], 0.5*gd$pars$v_sigma, 2*gd$pars$v_sigma)
  expect_between(summ$v_AR1[, "Mean"], -1, -0.3)
  # next with truncated normal prior
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") +
      gen(factor = ~ AR1(t, pr_truncnormal(mean=0, precision=1e6, lower=-1, upper=1)), priorA=pr_exp(), name="v"),
    data=dat, sigma.fixed=TRUE,
    control=sampler_control(block=FALSE)
  )
  sim <- MCMCsim(sampler, n.iter=500, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.25*c(1,2), 4*c(1,2))
  expect_between(summ$v_AR1[, "Mean"], -0.2, 0.2)
})

n <- 1000
dat <- data.frame(
  f1 = factor(sample(10, n, replace=TRUE)),
  f2 = factor(sample(21, n, replace=TRUE)),
  f3 = factor(sample(114, n, replace=TRUE))
)
info <- get_factor_info(~ f1 * AR1(f2, phi=0.5) * RW1(f3), data=dat)

test_that("DA_AR1_template works", {
  DA0.5 <- compute_GMRF_matrices(info, Q=FALSE, R=FALSE, D=TRUE)$D
  DA.templ <- DA_AR1_template(info, DA0.5 = DA0.5, nr=2)
  expect_equal(
    DA.templ$update(0.23),
    compute_GMRF_matrices(~ f1 * AR1(f2, phi=0.23) * RW1(f3), Q=FALSE, R=FALSE, D=TRUE, data=dat)$D
  )
})

test_that("QA_AR1_template works", {
  QA0.5 <- compute_GMRF_matrices(info, Q=TRUE, R=FALSE, D=FALSE)$Q
  QA.templ <- QA_AR1_template(info, QA0.5 = QA0.5, nr=2)
  expect_equal(
    QA.templ$update(0.66),
    compute_GMRF_matrices(~ f1 * AR1(f2, phi=0.66) * RW1(f3), Q=TRUE, R=FALSE, D=FALSE, data=dat)$Q
  )
})
