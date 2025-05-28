
context("Offsets")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000
y <- 1 + rnorm(n)

test_that("offset-only mean model works", {
  sampler <- create_sampler(y ~ offset(1))
  expect_length(sampler$mod, 1L)
  expect_equal(sampler$mod[[1]]$offset, 1)
  sim <- MCMCsim(sampler, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$sigma_[, "Mean"], 0.75, 1.3)
  summfitted <- summary(fitted(sim, units = 1:10))
  expect_equal(as.numeric(summfitted[, "Mean"]), rep(1, 10))
  expect_equal(as.numeric(summfitted[, "SD"]), rep(0, 10))
  # prediction
  summpred <- summary(predict(sim, iters=sample(n, 250), show.progress=FALSE))
  expect_between(summpred[, "Mean"], 0.75, 1.3)
  expect_between(summpred[, "SD"], 0.75, 1.3)
  summpred <- summary(predict(sim, newdata=data.frame(x=runif(10)), show.progress=FALSE))
  expect_between(summpred[, "Mean"], 0.75, 1.3)
  expect_between(summpred[, "SD"], 0.75, 1.3)
  summpred <- summary(predict(sim, newdata=data.frame(x=runif(10)), type="response", show.progress=FALSE))
  expect_equal(as.numeric(summpred[, "Mean"]), rep(1, 10))
  expect_equal(as.numeric(summpred[, "SD"]), rep(0, 10))

  # zero mean model
  sampler <- create_sampler(y ~ offset(0))
  expect_equal(sampler$mod[[1]]$offset, 0)
  sim <- MCMCsim(sampler, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$sigma_[, "Mean"], 0.75*sqrt(mean(y^2)), 1.3*sqrt(mean(y^2)))
  sampler <- create_sampler(y ~ mc_offset(value=0))
  expect_equal(sampler$mod[[1]]$offset, 0)
  sim <- MCMCsim(sampler, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$sigma_[, "Mean"], 0.75*sqrt(mean(y^2)), 1.3*sqrt(mean(y^2)))
  compute_DIC(sim)
  compute_WAIC(sim)
  # prediction
  summpred <- summary(predict(sim, iters=sample(n, 250), show.progress=FALSE))
  expect_between(summpred[, "Mean"], -0.4, 0.4)
  expect_between(summpred[, "SD"], 0.75*sqrt(mean(y^2)), 1.3*sqrt(mean(y^2)))
  summpred <- summary(predict(sim, newdata=data.frame(x=runif(10)), show.progress=FALSE))
  expect_between(summpred[, "Mean"], -0.4, 0.4)
  expect_between(summpred[, "SD"], 0.75*sqrt(mean(y^2)), 1.3*sqrt(mean(y^2)))
})
