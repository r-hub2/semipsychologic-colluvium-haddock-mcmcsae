
context("Binomial regression")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n)
)
b <- c(0.8, 2, 1)
df$y <- rbinom(n, 1, prob=1/(1 + exp(-(b[1] + b[2]*df$x1 + b[3]*df$x2))))

test_that("logistic regression works", {
  expect_warning(sampler <- create_sampler(y ~ reg(~ x1+x2, Q0=0.1), family="binomial", data=df), "deprecated")
  expect_equal(sampler$mod[[1L]]$prior$precision, Diagonal(x=rep.int(0.1, 3)))
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.3*b, 3*b)
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

test_that("matrix specification of response variable is possible", {
  sampler0 <- create_sampler(y ~ reg(~ x1+x2, prior=pr_normal(precision=0.1)), family="binomial", data=df)
  sampler <- create_sampler(cbind(y, 1-y) ~ reg(~ x1+x2, prior=pr_normal(precision=0.1)), family="binomial", data=df)
  expect_equal(sampler0$y, sampler$y)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.3*b, 3*b)
})
  
test_that("logistic regression with non-zero prior mean", {
  sampler <- create_sampler(y ~ reg(~ x1+x2, prior=pr_normal(mean=c(0, 0, b[3]), precision=c(0, 0, 1e6))),
    family="binomial", data=df)
  expect_length(sampler$block, 0L)
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1["x2", "Mean"], 0.9*b[3], 1.1*b[3])
})

df$y <- rbinom(n, 1, prob=pnorm(b[1] + b[2]*df$x1 + b[3]*df$x2))
test_that("probit regression works", {
  sampler <- create_sampler(y ~ reg(~ x1 + x2), family=f_binomial(link="probit"), data=df)
  expect_equal(sampler$family$ny, 1)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.3*b, 3*b)
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

y <- rbinom(1000L, 1L, prob=0.3)
test_that("minimal probit binary regression prediction example works", {
  sampler <- create_sampler(y ~ 1, family=f_binomial(link="probit"))
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) pnorm(x)))
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, type="response", iters=sample(1:500, 250), show.progress=FALSE))  # in-sample probability
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), show.progress=FALSE))  # out-of-sample
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), type="response", show.progress=FALSE))  # out-of-sample probability
  expect_between(summ[, "Mean"], 0.2, 0.4)
})

y <- 30; ny <- 100
test_that("minimal logistic binomial regression example works", {
  sampler <- create_sampler(y ~ 1, family=f_binomial(n.trial=ny))
  expect_equal(sampler$family$ny, ny)
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) 1/(1+exp(-x))))
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, show.progress=FALSE))  # in-sample
  expect_between(summ[, "Mean"], 20, 40)
  summ <- summary(predict(sim, type="response", show.progress=FALSE))  # in-sample probability
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), weights=2, show.progress=FALSE))  # out-of-sample
  expect_between(summ[, "Mean"], 40, 80)
  expect_warning(summ <- summary(predict(sim, newdata=data.frame(1), weights=2, type="response", show.progress=FALSE)), "ignored")
  expect_between(summ[, "Mean"], 0.2, 0.4)
})

test_that("minimal logistic binomial regression example works if response is specified as fractions", {
  y <- 0.3; ny <- 100
  sampler <- create_sampler(y ~ 1, family=f_binomial(n.trial=ny))
  expect_equal(sampler$family$ny, ny)
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) 1/(1+exp(-x))))
  expect_between(summ[, "Mean"], 0.2, 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), weights=2, show.progress=FALSE))  # out-of-sample
  # now test the case where ny (and y) are in a data frame, and n.trial a formula
  dat <- data.frame(
    y = c(0.1, 0.1),
    ny = c(100, 200)
  )
  sampler <- create_sampler(y ~ 1, family=f_binomial(n.trial = ~ ny), data=dat)
  expect_equal(sampler$family$ny, dat$ny)
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) 1/(1+exp(-x))))
  expect_between(summ[, "Mean"], 0.05, 0.2)
  summ <- summary(predict(sim, newdata=data.frame(ny=1000), show.progress=FALSE))  # out-of-sample
  expect_between(summ[, "Mean"], 75, 125)
  summ <- summary(predict(sim, newdata=data.frame(ny=1000, COUNT=2), weights = ~ COUNT, show.progress=FALSE))  # out-of-sample
  expect_between(summ[, "Mean"], 150, 250)
})

test_that("binomial multilevel model runs", {
  ex <- mcmcsae_example(1000, family="binomial")
  sampler <- create_sampler(ex$model, data=ex$dat, family="binomial",
    control=sampler_control(max.size.cps.template = 0))
  expect_false(is.function(sampler$mbs[[1]]$cps_template))
  sampler <- create_sampler(ex$model, data=ex$dat, family="binomial",
    control=sampler_control(max.size.cps.template=100))
  expect_length(sampler$mbs, 1L)
  expect_true(is.function(sampler$mbs[[1]]$cps_template))
  sim <- MCMCsim(sampler, burnin=50, n.iter=100, n.chain=2, verbose=FALSE)
  expect_in("v_rho", names(summary(sim)))
  sampler <- create_sampler(ex$model, data=ex$dat, family="binomial",
    control=sampler_control(block=list(c("v", "u"))))
  expect_identical(sort(names(sampler$mbs[[1]]$mcs)), c("u", "v"))
})
