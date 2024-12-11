
context("Gamma regression")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n)
)
b <- c(0.8, 2, 1)
alpha <- 1
mu <- exp(b[1] + b[2]*df$x1 + b[3]*df$x2)
df$y <- rgamma(n, shape=alpha, rate=alpha/mu)

test_that("gamma regression works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family="gamma", data=df)
  expect_equal(sampler$family$shape.prior$type, "gamma")
  expect_true("gamma_shape_" %in% names(sampler$rprior()))
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.2 * b, 5 * b)
  expect_length(acceptance_rates(sim)[["gamma_shape_"]], 2L)
  expect_between(summ$gamma_shape_[, "Mean"], 0.25 * 1, 4 * 1)
  #plot(sim, "gamma_shape_")
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

test_that("gamma regression prediction works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family="Gamma", data=df[1:900, ])
  # both "gamma" and "Gamma" (as used in stats) are allowed
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$reg1[, "Mean"]), b, tolerance=1)
  pred <- predict(sim, newdata = df[901:1000, ], show.progress = FALSE)
  summpred <- summary(pred)
  #plot(log(summpred[, "Mean"]), log(df$y[901:1000])); abline(0, 1)
})

alpha <- 10
df$y <- rgamma(n, shape=alpha, rate=alpha/mu)
test_that("gamma regression with fixed shape works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family=f_gamma(shape.prior = pr_fixed(10)), data=df)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.2 * b, 5 * b)
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

test_that("exponential prior on shape works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family=f_gamma(shape.prior = pr_exp(1)), data=df)
  expect_equal(sampler$family$shape.prior$type, "gamma")
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$gamma_shape_[, "Mean"], 0.2 * alpha, 5 * alpha)
  expect_between(summ$reg1[, "Mean"], 0.2 * b, 5 * b)
})

test_that("gamma regression with offset works", {
  alpha <- 1
  df$o <- rexp(n)
  mu <- exp(df$o + b[1] + b[2]*df$x1 + b[3]*df$x2)
  df$y <- rgamma(n, shape=alpha, rate=alpha/mu)
  sampler <- create_sampler(
    y ~ offset(o) + reg(~ x1 + x2, name="beta"),
    family="gamma", data = df
  )
  expect_length(sampler$block, 0L)
  expect_equal(sampler$mod[[length(sampler$mod)]]$offset, df$o)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.25 * b, 4 * b)
  expect_warning(sampler <- create_sampler(
      y ~ offset(o) + reg(~ x1 + x2, name="beta"),
      family="gamma", data = df, control=sampler_control(block=list("beta"))
    ), "only one"
  )
  expect_length(sampler$block, 1L)
  expect_equal(sampler$mod[[length(sampler$mod)]]$offset, df$o)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.25 * b, 4 * b)
})

n <- 1000
dat <- data.frame(
  f = as.factor(sample(1:30, n, replace=TRUE))
)
vf <- rnorm(30, sd=10)
sh <- 0.6
dat$y <- rgamma(n, shape=sh, rate = sh * exp(- vf[dat$f]))  # NB rate = shape / mean = shape * exp(-eta)
test_that("gamma family with single gen component works", {
  sampler <- create_sampler(
    y ~ gen(factor = ~ f, name="v"),
    family="gamma", data=dat
  )
  expect_false(sampler$mod[["v"]]$usePX)
  sim <- MCMCsim(sampler, n.iter=800, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 10, 3 * 10)
  expect_between(summ$gamma_shape_[, "Mean"], 0.25 * 0.6, 4 * 0.6)
})

n <- 1000
dat <- data.frame(
  x = runif(n),
  z = rnorm(n),
  f = as.factor(sample(1:30, n, replace=TRUE))
)
sh <- 0.6
b <- c(1, 1, 1)
test_that("gamma multilevel model works", {
  gd <- generate_data(
    ~ reg(~ 1 + x + z, prior=pr_fixed(c(1, 1, 1)), name="beta") + gen(factor = ~ f, prior=pr_invchisq(1e6, 2), name="v"),
    family = f_gamma(shape.prior=pr_fixed(sh)), data=dat
  )
  expect_equal(unname(gd$pars$beta), b)
  expect_between(gd$pars$v_sigma, 0.99 * sqrt(2), 1.01 * sqrt(2))
  sampler <- create_sampler(
    gd$y ~ reg(~ x + z, name="beta") + gen(factor = ~ f, name="v"),
    family="gamma", data=dat
  )
  expect_length(sampler$block, 1L)
  expect_false(sampler$mod[["v"]]$usePX)
  sim <- MCMCsim(sampler, n.iter=800, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.3 * b, 3 * b)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * sqrt(2), 3 * sqrt(2))
  expect_between(summ$gamma_shape_[, "Mean"], 0.25 * sh, 4 * sh)
})
