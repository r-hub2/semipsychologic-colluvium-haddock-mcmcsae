
context("Output from a simple linear regression example")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("output of create_sampler and MCMCsim is as expected", {
  n <- 1000L
  df <- data.frame(x = rnorm(n))
  df$y <- 1 + df$x + rnorm(n)
  sampler <- create_sampler(y ~ reg(~ x, name="beta"), data=df)
  expect_false(is.null(sampler$draw))
  expect_false(is.null(sampler$start))
  sim <- MCMCsim(sampler, n.iter=2000, n.chain=2, burnin=1000, verbose=FALSE)
  expect_identical(sort(par_names(sim)), c("beta", "llh_", "sigma_"))
  expect_is(n_eff(sim$beta), "numeric")
  expect_is(R_hat(sim$beta), "numeric")
  expect_is(subset(sim$beta, vars=1L, chains=2:1, draws=c(100:n_draws(sim$beta))), "dc")
  expect_is(get_means(sim), "list")
  expect_is(get_sds(sim), "list")
  summ <- summary(sim)
  expect_is(summ, "mcdraws_summary")
  expect_true(all(summ$beta[, "Mean"] - coef(lm(y ~ x, data=df)) < 0.1))
  DIC <- compute_DIC(sim)
  expect_true(DIC["p_DIC"] > 2 && DIC["p_DIC"] < 4)
  WAIC <- compute_WAIC(sim, n)
  expect_true(WAIC["p_WAIC1"] > 2 && WAIC["p_WAIC1"] < 4)
})
