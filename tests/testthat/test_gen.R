
context("Model component 'gen'")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("generating data for iid random effects at the data level works", {
  dat <- data.frame(t=1:100)
  gd <- generate_data(~ gen(factor=~t, name="v"), data=dat)
  expect_length(gd$pars$v, 100L)
})

n <- 1000L
m <- 25L
df <- data.frame(
  f = factor(sample(m, n, replace=TRUE))
)
vf <- rnorm(m, sd=0.4)

df$y <- with(df, vf[f] + rnorm(n))
test_that("gaussian model with a single gen component works", {
  sampler <- create_sampler(y ~ 1 + gen(factor = ~ f, name="v"), data=df)
  expect_true(sampler$mod$v$PX$data.scale)
  expect_identical(names(sampler$mod), c("reg1", "v"))  # explicit intercept
  sampler <- create_sampler(y ~ gen(factor = ~ f, name="v"), data=df)
  expect_identical(names(sampler$mod), "v")
  sim <- MCMCsim(sampler, n.chain = 2, n.iter=500, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 0.4, 3 * 0.4)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
})

df$y <- with(df, rbinom(n, 1, prob = 1 / (1 + exp(-vf[f]))))
test_that("binomial model with a single gen component works", {
  sampler <- create_sampler(y ~ gen(factor = ~ f, name="v"), data=df, family="binomial")
  sim <- MCMCsim(sampler, n.chain = 2, n.iter=500, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 0.4, 3 * 0.4)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
})

df$x <- runif(n)
df$y <- with(df, rbinom(n, 1, prob = 1 / (1 + exp(-(1 + 0.5*df$x + vf[f])))))
test_that("binomial multilevel model with non-zero prior mean in regression component works", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, prior=pr_normal(mean=c(0, 0.5), precision=c(0, 1e6))) +
        gen(factor = ~ f, name="v"),
    data=df, family="binomial"
  )
  expect_false(sampler$mod$v$PX$data.scale)
  expect_length(sampler$control$block[[1L]], 2L)
  sim <- MCMCsim(sampler, n.chain = 1, n.iter=500, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x", "Mean"], 0.48, 0.52)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 0.4, 3 * 0.4)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
  # no PX:
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, prior=pr_normal(mean=c(0, 0.5), precision=c(0, 1e6))) +
      gen(factor = ~ f, name="v", PX=FALSE),
    data=df, family="binomial"
  )
  expect_false(sampler$mod$v$PX)
  expect_length(sampler$control$block[[1L]], 2L)
  #sampler$mod[[2]]$draw  # no XX and Xy computed
  sim <- MCMCsim(sampler, n.chain = 1, n.iter=500, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x", "Mean"], 0.48, 0.52)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 0.4, 3 * 0.4)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
  expect_warning(
    sampler <- create_sampler(
      y ~ reg(~ 1 + x, prior=pr_normal(mean=c(0, 0.5), precision=c(0, 1e6))) +
        gen(factor = ~ f, name="v", PX=list(data.scale=TRUE)),
      data=df, family="binomial"
    ), "data.scale"
  )
})

test_that("a (binomial) multilevel model with constraints works", {
  RA <- cbind(rep(1, 25), c(1, rep(0, 24)))
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, prior=pr_normal(mean=c(0, 0.5), precision=c(0, 1e6)),
            constraints = set_constraints(R=cbind(c(0, 1)), r=0.45)) +
      gen(factor = ~ f, name="v", constraintsA=set_constraints(R=RA)),
    data=df, family="binomial"
  )
  expect_length(sampler$control$block[[1L]], 2L)
  sim <- MCMCsim(sampler, n.chain=2, burnin=50, n.iter=100, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x", "Mean"], 0.445, 0.455)
  expect_equal(crossprod_mv(RA, summ$v[, "Mean"]), c(0, 0), tolerance=1e-3)
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, prior=pr_normal(mean=c(0, 0.5), precision=c(0, 1e6)),
            constraints = set_constraints(R=cbind(c(0, 1)), r=0.45)) +
      gen(factor = ~ f, name="v", constraintsA=set_constraints(R=RA)),
    data=df, family="binomial", control=sampler_control(block=FALSE)
  )
  expect_length(sampler$block, 0L)
  sim <- MCMCsim(sampler, n.chain=2, burnin=50, n.iter=100, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x", "Mean"], 0.445, 0.455)
  expect_equal(crossprod_mv(RA, summ$v[, "Mean"]), c(0, 0), tolerance=1e-3)
})

n <- 2000L
m <- 100L
df <- data.frame(
  f = factor(sample(m, n, replace=TRUE))
)
vf <- 0.4 * rt(m, df=1)
df$y <- with(df, 2 + vf[f] + rnorm(n))
test_that("model with t-distributed random effects works", {
  sampler <- create_sampler(
    y ~ 1 + gen(factor = ~ f, priorA=pr_invchisq(df=2), name="v"),
    data=df, family="gaussian"
  )
  expect_equal(sampler$mod$v$priorA$type, "invchisq")
  expect_equal(sampler$family$df.sigma, n + sampler$mod[["v"]]$q0)
  sim <- MCMCsim(sampler, store.all=TRUE, n.chain=2, n.iter=700, verbose=FALSE)
  summ <- summary(sim)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
  expect_equal(unname(summ$v[, "Mean"]), vf, tolerance=0.5)
})

n <- 60
df <- data.frame(
  x=runif(n),
  z=rnorm(n),
  u=rgamma(n, shape = 2)
)
df$y <- 1 - df$x + df$z - df$u + rnorm(n, sd=0.4)
test_that("model with gen component with empty factor argument works", {
  sampler <- create_sampler(
    y ~ gen(~ x + z + u, var="diagonal", name="beta"),
    data=df
  )
  expect_equal(sampler$mod[[1]]$l, 1L)
  expect_true(sampler$mod$beta$PX$data.scale)
  expect_equal(sampler$family$df.sigma, n + sampler$mod$beta$q0)
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(c(1, -1, 1, -1), summ$beta[, "Mean"] - 5 * summ$beta[, "SD"], summ$beta[, "Mean"] + 5 * summ$beta[, "SD"])
})

n <- 1000L
m1 <- 10L
m2 <- 21L
df <- data.frame(
  x = rnorm(n),
  f1 = factor(sample(m1, n, replace=TRUE)),
  f2 = factor(sample(m2, n, replace=TRUE))
)
vf1 <- 0.4 * rt(m1, df=1)
vf2 <- 0.7 * rnorm(m2)
df$y <- df$x + vf1[df$f1] + vf2[df$f2]
test_that("specification of random effects using lme4 syntax works", {
  sampler <- create_sampler(
    y ~ x + (1|f1) + (1 | f2), data=df
  )
  expect_length(sampler$mod, 3)
  expect_identical(unname(sapply(sampler$mod, `[[`, "q")), c(2L, m1, m2))
  sim <- MCMCsim(sampler, burnin=10, n.iter=20, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  sampler <- create_sampler(
    y ~ x + (1|f1/f2), data=df
  )
  expect_length(sampler$mod, 3)
  expect_identical(unname(sapply(sampler$mod, `[[`, "q")), c(2L, m1, sum(table(df$f1, df$f2) > 0)))
  sim <- MCMCsim(sampler, burnin=10, n.iter=20, n.chain=2, store.all=TRUE, verbose=FALSE)
  sampler <- create_sampler(
    y ~ x + gen(factor = ~ f1) + (x || f2), data=df
  )
  expect_length(sampler$mod, 3)
  sim <- MCMCsim(sampler, burnin=10, n.iter=20, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_length(summ$gen3_sigma[, "Mean"], 2L)
  expect_null(summ$gen3_rho)
  expect_error(
    create_sampler(y ~ x + us(x | f1), data=df), "unsupported"
  )
})
