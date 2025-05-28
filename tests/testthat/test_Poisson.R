
context("Poisson regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# generate Poisson data
n <- 1000
x <- rnorm(n)
eta <- 1 + 2*x
y <- rpois(n, lambda = exp(eta))

test_that("negative binomial approximation to Poisson regression works", {
  r <- 100  # should be large for negligible overdispersion, but too large values result in slow MCMC mixing
  o <- rep(-log(r), n)
  sampler <- create_sampler(y ~ 1 + x + offset(o), family=f_negbinomial(inv.shape.prior = pr_fixed(1/r)))
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(summ$reg1[, "Mean"], c(`(Intercept)`=1, x=2), tolerance=0.25)
  predsumm <- summary(predict(sim, show.progress = FALSE))
  expect_between(mean(predsumm[, "Mean"]), 0.75*mean(y), 1.3*mean(y))
  compute_DIC(sim)
  compute_WAIC(sim)
})

test_that("Poisson shortcut works", {
  sampler <- create_sampler(y ~ 1 + x, family="poisson")
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(summ$reg1[, "Mean"], c(`(Intercept)`=1, x=2), tolerance=0.25)
  predsumm <- summary(predict(sim, show.progress = FALSE))
  expect_between(mean(predsumm[, "Mean"]), 0.75*mean(y), 1.3*mean(y))
  compute_DIC(sim)
  compute_WAIC(sim)
})

test_that("f_poisson works", {
  expect_error(sampler <- create_sampler(y ~ 1 + x, ry=100, family="poisson"), "no longer")
  sampler <- create_sampler(y ~ 1 + x, family=f_poisson(control=poisson_control(nb.shape=200)))
  expect_equal(sampler$ry, 200)
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(summ$reg1[, "Mean"], c(`(Intercept)`=1, x=2), tolerance=0.25)
  predsumm <- summary(predict(sim, show.progress = FALSE))
  expect_between(mean(predsumm[, "Mean"]), 0.75*mean(y), 1.3*mean(y))
  compute_DIC(sim)
  compute_WAIC(sim)
})

test_that("in-sample prediction works", {
  sampler <- create_sampler(y ~ 1 + x, linpred = "fitted",
                            family=f_poisson(control=poisson_control(nb.shape=200)))
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=400, verbose=FALSE)
  summ <- summary(sim)
  fitted <- transform_dc(sim$linpred_, fun = exp)
  fittedsumm <- summary(fitted)
  expect_between(mean(fittedsumm[, "Mean"]), 0.75*mean(y), 1.3*mean(y))
  fitted <- fitted(sim, type="response")
  fittedsumm <- summary(fitted)
  expect_between(mean(fittedsumm[, "Mean"]), 0.75*mean(y), 1.3*mean(y))
})
