
context("Input checks")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
dat <- data.frame(
  x = rnorm(n),
  z = runif(n)
)
dat$y <- rnorm(n)

test_that("missing values in y are flagged", {
  dat$y[sample.int(n, 2L)] <- NA
  expect_error(create_sampler(y ~ x + z, data=dat), "2 missing")
  expect_error(create_sampler(y ~ x + z, family="binomial", data=dat), "2 missing")
})

test_that("range for binomial response is checked", {
  dat$y <- rnorm(n)
  expect_error(sampler <- create_sampler(y ~ x + z, data=dat, family="binomial"), "negative")
  expect_error(create_sampler(y ~ x + z, data=dat, family=f_binomial(n.trial=100)), "negative")
  expect_error(create_sampler(y ~ x + z, data=dat, family=f_binomial(n.trial=-1)), "negative")
})

test_that("logistic binomial model works for non-integral data", {
  dat$y <- runif(n)
  expect_warning(sampler <- create_sampler(y ~ 1, data=dat, family="binomial"), "non-integral")
  sim <- MCMCsim(sampler, n.chain=2, burnin=100, n.iter=300, verbose=FALSE)
  expect_equal(summary(sim$reg1)[, "Mean"], 0, tolerance=0.2)
})

test_that("for negative binomial and Poisson families negative response values are flagged", {
  dat$y <- rnorm(n)
  expect_warning(sampler <- create_sampler(y ~ x + z, data=dat, family="negbinomial"), "negative")
  # this still runs:
  sim <- MCMCsim(sampler, n.iter=200, n.chain=2, verbose=FALSE)
  summary(sim)
  expect_warning(create_sampler(y ~ x + z, data=dat, family="poisson"), "negative")
})

test_that("response variable for multinomial family is checked", {
  dat$y <- rnorm(n)
  expect_error(create_sampler(y ~ x + z, data=dat, family="multinomial"))
  dat$y <- cbind(rbinom(n, 2, prob=0.2), rbinom(n, 2, prob=0.2), rep(1, n))
  sampler <- create_sampler(y ~ x + z, data=dat, family="multinomial")
  expect_identical(sampler$Km1, 2L)
  dat$y[1, 1] <- -1L
  expect_error(create_sampler(y ~ x + z, data=dat, family="multinomial"), "negative")
})

test_that("zeroes are flagged for gamma family", {
  dat$y <- rgamma(n, shape=1e-4)  # this yields many 0s due to numerical underflow
  expect_error(create_sampler(y ~ 1, data=dat, family="gamma"), "strictly positive")
})

test_that("using 'vreg' or 'vfac' model component in formula is flagged", {
  dat$y <- rnorm(n)
  dat$g <- sample(1:10, n, replace=TRUE)
  expect_error(create_sampler(y ~ 1 + gen(factor = ~ g) + vreg(factor="g")), "variance model")
  expect_error(create_sampler(y ~ 1 + x + vfac(factor="g")), "variance model")
})

test_that("empty formulas are handled correctly", {
  expect_equal(standardize_formula(y ~ 0), y ~ 0)
  expect_equal(standardize_formula(~ 0), ~ 0)
  expect_error(sampler <- create_sampler( ~ 0), "empty")
})

test_that("family argument can also be a function", {
  dat <- data.frame(
    y = c(1, 0, 1),
    x = c(0.5, 0.6, 2)
  )
  sampler <- create_sampler(y ~ x, family=binomial, data=dat)
  expect_is(sampler$family, "environment")
  expect_identical(sampler$family$family, "binomial")
  gd <- generate_data(~ reg(~ x, prior=pr_fixed(1), name="beta"), family=binomial, data=dat)
  expect_equal(unname(gd$pars$beta), c(1, 1))
  sampler <- create_sampler(y ~ x, family=f_gaussian, data=dat)
  expect_is(sampler$family, "environment")
  expect_identical(sampler$family$family, "gaussian")
  gd <- generate_data(~ reg(~ x, prior=1, name="beta"), family=f_gaussian, data=dat)
  expect_equal(unname(gd$pars$beta), c(1, 1))
  expect_length(gd$y, nrow(dat))
})
