
context("Model component 'gen'")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
m <- 25L
df <- data.frame(
  f = factor(sample(m, n, replace=TRUE))
)
vf <- rnorm(m, sd=0.4)

df$y <- with(df, vf[f] + rnorm(n))
test_that("gaussian model with a single gen component works", {
  sampler <- create_sampler(y ~ 1 + gen(factor = ~ f, name="v"), data=df)
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
  expect_length(sampler$block[[1L]], 2L)
  sim <- MCMCsim(sampler, n.chain = 2, n.iter=500, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x", "Mean"], 0.49, 0.51)
  expect_between(summ$v_sigma[, "Mean"], 0.3 * 0.4, 3 * 0.4)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
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
  sim <- MCMCsim(sampler, store.all=TRUE, n.chain=2, n.iter=700, verbose=FALSE)
  summ <- summary(sim)
  #plot(vf, summ$v[, "Mean"]); abline(0, 1)
  expect_equal(unname(summ$v[, "Mean"]), vf, tolerance=0.5)
})
