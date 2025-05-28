
context("Group-level covariates")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 100L
df <- data.frame(x=runif(n), t=1:n)
dat <- generate_data(
  family = f_gaussian(var.prior = pr_invchisq(df=1e6, scale=1)),
  formula = ~ reg(~ 1 + x, prior=pr_fixed(1:2), name="beta") + gen(factor = ~ RW1(t), PX=FALSE, prior=pr_invchisq(df=1e6, scale=0.1), name="v"),
  data = df
)
df$y <- dat$y

test_that("generated parameters are as expected", {
  expect_equivalent(dat$pars$sigma_, 1, tol=0.2)
  expect_equivalent(dat$pars$beta, 1:2)
  expect_equivalent(dat$pars$v_sigma, sqrt(0.1), tol=0.2)
  gd <- generate_data(
    formula = ~ reg(~ 1 + x, prior=pr_normal(mean=1:2, precision=1e6), name="beta"),
    data = df
  )
  expect_equivalent(gd$pars$beta, 1:2, tol=0.2)
})

test_that("non-centered sampler runs", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ RW1(t), name="v"),
    data=df
  )
  sim <- MCMCsim(sampler, n.chain=2L, n.iter=250L, verbose=FALSE)
  expect_length(acceptance_rates(sim), 0L)
  summ <- summary(sim)
  expect_equal(rownames(summ$beta), names(dat$pars$beta))
  expect_equivalent(summ$sigma_[, "Mean"], dat$pars$sigma_, tol=1)
  expect_equivalent(summ$beta[, "Mean"], dat$pars$beta, tol=1)
})

test_that("partially centered sampler (intercept) runs", {
  sampler <- create_sampler(y ~ reg(~ 0 + x, name="beta") +
      gen(factor = ~ RW1(t), formula.gl = ~ glreg(~ 1), name="v"),
    data=df
  )
  sim <- MCMCsim(sampler, n.chain=2L, n.iter=250L, verbose=FALSE)
  expect_length(acceptance_rates(sim), 1L)
  expect_gte(acceptance_rates(sim)[[1L]][[1L]], 0)
  expect_lte(acceptance_rates(sim)[[1L]][[1L]], 1)
  summ <- summary(sim)
  expect_equal(rownames(summ$beta), "x")
  expect_equal(rownames(summ$v_gl), "(Intercept)")
  expect_equivalent(summ$sigma_[, "Mean"], dat$pars$sigma_, tol=1)
  expect_equivalent(summ$beta[, "Mean"], dat$pars$beta["x"], tol=1)
  expect_equivalent(summ$v_gl[, "Mean"], dat$pars$beta["(Intercept)"], tol=1)
})

test_that("partially centered sampler (x) runs", {
  sampler <- create_sampler(y ~ reg(~ 1, name="beta") +
      gen(factor = ~ RW1(t), formula.gl = ~ glreg(~ 0 + x), name="v"),
    data=df
  )
  sim <- MCMCsim(sampler, n.chain=2L, n.iter=250L, verbose=FALSE)
  expect_length(acceptance_rates(sim), 1L)
  summ <- summary(sim)
  expect_equal(rownames(summ$beta), "(Intercept)")
  expect_equal(rownames(summ$v_gl), "x")
  expect_equivalent(summ$sigma_[, "Mean"], dat$pars$sigma_, tol=1)
  expect_equivalent(summ$beta[, "Mean"], dat$pars$beta["(Intercept)"], tol=1)
  expect_equivalent(summ$v_gl[, "Mean"], dat$pars$beta["x"], tol=1)
})

test_that("centered sampler runs", {
  sampler <- create_sampler(
    y ~ gen(factor = ~ RW1(t), formula.gl = ~ glreg(~ x), name="v"),
    data=df
  )
  sim <- MCMCsim(sampler, n.chain=2L, n.iter=250L, verbose=FALSE)
  expect_length(acceptance_rates(sim), 1L)
  summ <- summary(sim)
  expect_equal(rownames(summ$v_gl), names(dat$pars$beta))
  expect_equivalent(summ$sigma_[, "Mean"], dat$pars$sigma_, tol=1)
  expect_equivalent(summ$v_gl[, "Mean"], dat$pars$beta, tol=1)
})

test_that("assigning a normal prior precision for group-level effects works", {
   sampler <- create_sampler(
     y ~ gen(factor = ~ RW1(t), formula.gl = ~ glreg(~ x, prior=pr_normal(precision = 1, labels="x")), name="v"),
     data=df
   )
   expect_equal(sampler$mod[[1]]$glp$Q0[2, 2], 1)
   sim <- MCMCsim(sampler, n.chain=1L, n.iter=250L, verbose=FALSE)
   expect_length(acceptance_rates(sim), 1L)
   summ <- summary(sim)
   expect_equal(rownames(summ$v_gl), names(dat$pars$beta))
   expect_equivalent(summ$sigma_[, "Mean"], dat$pars$sigma_, tol=1)
   expect_equivalent(summ$v_gl[, "Mean"], dat$pars$beta, tol=1)
})

n <- 1000
m <- 25
df <- data.frame(
  area = sample(m, n, replace=TRUE),
  x1 = 8*runif(n),
  x2 = 8*runif(n)
)
betaTrue <- c(0.6, 0.3, 0.2)
dat <- generate_data(
  formula = ~ reg(~ x1 + x2, prior=pr_fixed(betaTrue), name="beta") +
              gen(~ x1 + x2, var="diagonal", factor = ~ area, PX=FALSE,
                  prior=pr_invchisq(df=1e6, scale=c(0.5, 0.03, 0.03)), name="v"),
  data = df, family="binomial"
)
df$y <- dat$y

test_that("more complex example works", {
  # uncentered
  sampler <- create_sampler(
    y ~ reg(name="beta", formula=~x1+x2) + gen(name="v", formula=~x1+x2, factor=~area),
    family="binomial", data=df,
    control=sampler_control(PG.approx = TRUE, PG.approx.m = 0L)
  )
  expect_true(sampler$control$PG.approx)
  expect_equal(sampler$control$PG.approx.m, 0L)
  sim <- MCMCsim(sampler, n.iter=700, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(nrow(summ$v_rho) == 3L)
  expect_true(all(abs(summ$v_rho[, "Mean"]) < 0.7))
  expect_between(summ$beta[, "Mean"], 0 * betaTrue, 2 * betaTrue)
  # centered
  sampler <- create_sampler(
    y ~ gen(name="v", formula=~x1+x2, factor=~iid(area), formula.gl=~glreg(~1), debug=FALSE),
    family="binomial", data=df,
    control=sampler_control(PG.approx.m = 0L)
  )
  sim <- MCMCsim(sampler, n.iter=700, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all(abs(summ$v_rho[, "Mean"]) < 0.7))
  expect_between(summ$v_gl[, "Mean"], 0 * betaTrue, 2 * betaTrue)
  sampler <- create_sampler(
    y ~ gen(name="v", factor=~iid(area), formula.gl=~glreg(~x1+x2), debug=FALSE),
    family="binomial", data=df,
    control=sampler_control(PG.approx.m = 0L)
  )
  sim <- MCMCsim(sampler, n.iter=700, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
})
