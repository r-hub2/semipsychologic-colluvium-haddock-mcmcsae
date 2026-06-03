
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
  expect_error(
    sim <- MCMCsim(sampler, n.iter=200, n.chain=2,
      start = list(list(negbin_shape_ = -1), list(negbin_shape_ = 1))
    ), "positive"
  )
  expect_error(
    sim <- MCMCsim(sampler, n.iter=200, n.chain=2,
      start = list(list(negbin_shape_ = NA), list(negbin_shape_ = 1))
    ), "NA"
  )
  expect_error(
    sim <- MCMCsim(sampler, n.iter=200, n.chain=2,
      start = list(list(negbin_shape_ = 1:2), list(negbin_shape_ = 3:4))
    ), "length"
  )
  expect_warning(create_sampler(y ~ x + z, data=dat, family="poisson"), "negative")
})

test_that("response variable for multinomial family is checked", {
  dat$y <- rnorm(n)
  expect_error(create_sampler(y ~ x + z, data=dat, family="multinomial"))
  dat$y <- cbind(rbinom(n, 2, prob=0.2), rbinom(n, 2, prob=0.2), rep(1, n))
  sampler <- create_sampler(y ~ x + z, data=dat, family="multinomial")
  expect_identical(sampler$family$Km1, 2L)
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
  expect_equal(standardise_formula(y ~ 0), y ~ 0)
  expect_equal(standardise_formula(~ 0), ~ 0)
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

test_that("use of same name for multiple components is not allowed", {
  n <- 100
  dat <- data.frame(
    y = rnorm(n),
    x = runif(n),
    g = factor(sample(letters, n, replace=TRUE))
  )
  expect_error(
    sampler <- create_sampler(
      ~ reg(~ x, name="beta") + gen(~ x, factor = ~ g, name="beta"),
      data = dat
    ), "duplicate"
  )
  expect_error(
    sampler <- create_sampler(
      y ~ reg(~ x, name="beta") + gen(~ x, factor = ~ g, name="beta"),
      data = dat
    ), "duplicate"
  )
  expect_error(
    sampler <- create_sampler(
      ~ reg(~ x, name="beta_sigma") + gen(~ x, factor = ~ g, name="beta"),
      data = dat
    ), "not allowed"
  )
  expect_error(
    sampler <- create_sampler(
      y ~ reg(~ x, name="beta") + gen(~ x, factor = ~ g, name="v"),
       data = dat,
       family=f_gaussian(var.model = ~ vreg(~ 0 + x, name="v"))
    ), "must be distinct"
  )
})

test_that("wrong input for 'family' triggers an appropriate error", {
  expect_error(sampler <- create_sampler(y ~ x, family=33), "list")
  expect_error(sampler <- create_sampler(y ~ x, family=list(family="gaussian", link="identity")), "list")
})

test_that("gen components' factor term is checked for invalid factor components", {
  expect_error(
    sampler <- create_sampler(
      y ~ x + gen(factor = ~ seq_len(n)), data=dat
    ), "unsupported"
  )
  expect_error(
    sampler <- create_sampler(
      y ~ x + gen(factor = ~ x * RWq(z)), data=dat
    ), "unsupported"
  )
  expect_error(
    sampler <- create_sampler(
      y ~ x + gen(factor = ~ 1:n, name="v"), data=dat
    ), "length"
  )
  # this works:
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ local_, name="v"), data=dat
  )
  expect_equal(sampler$mod[["v"]]$q, n)
  # or this:
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ iid(seq_len(n)), name="v"), data=dat
  )
  expect_equal(sampler$mod[["v"]]$q, n)
})

test_that("starting values are checked", {
  sampler <- create_sampler(y ~ x + z, data=dat)
  expect_error(
    sim <- MCMCsim(sampler, n.chain=1, start=list(list(sigma_ = -0.1))),
    "not positive"
  )
  sampler <- create_sampler(y ~ x + z + gen(factor = ~ iid(local_), name="v"), data=dat)
  v0 <- rnorm(nrow(dat))
  v0[2] <- -Inf
  expect_error(
    sim <- MCMCsim(sampler, n.chain=1,
      start=list(list(sigma_ = 0.1, v = v0))
    ),
    "Inf"
  )
})
