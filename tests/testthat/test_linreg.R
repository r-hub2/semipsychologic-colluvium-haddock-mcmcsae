
context("Linear regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n),
  x3 = rbinom(n, 1, 0.1),
  x4 = rgamma(n, 1, 1)
)
df$y <- with(df, 1 + x1 + 2*x2 + 3*x3 + 4*x4 + rnorm(n))

test_that("linear regression works, as well as fitted, residuals and predict methods", {
  sampler <- create_sampler(y ~ reg(~1+x1+x2+x3+x4, name="beta"), data=df)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5 * c(1,1,2,3,4), 2 * c(1,1,2,3,4))
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
  expect_equal(fitted(sim, mean.only=TRUE), as.vector(summary(fitted(sim))[, "Mean"]))
  expect_equal(residuals(sim, mean.only=TRUE), as.vector(summary(residuals(sim))[, "Mean"]))
  expect_equal(nrow(summary(predict(sim, iters=sample(1:500, 100), show.progress=FALSE))), n)
  expect_equal(n_vars(predict(sim, X.=list(beta=model_matrix(~1+x1+x2+x3+x4, df[1:10, ])), show.progress=FALSE)), 10)
  expect_equal(n_vars(predict(sim, newdata=df[21, ], show.progress=FALSE)), 1)
})

test_that("single-site Gibbs sampler works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
        reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"), data=df,
    control=sampler_control(block=FALSE)
  )
  expect_null(sampler$mbs)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true((0.5 < summ$mu[, "Mean"]) && (summ$mu[, "Mean"] < 2))
  expect_true((0.5 < summ$beta1[, "Mean"]) && (summ$beta1[, "Mean"] < 2))
  expect_true((2*0.5 < summ$beta2[, "Mean"]) && (summ$beta2[, "Mean"] < 2*2))
  expect_true((3*0.5 < summ$beta3[, "Mean"]) && (summ$beta3[, "Mean"] < 3*2))
  expect_true((4*0.5 < summ$beta4[, "Mean"]) && (summ$beta4[, "Mean"] < 4*2))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("fully blocked Gibbs sampler works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
        reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"),
    data=df
  )
  expect_length(sampler$mbs, 1L)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$mu[, "Mean"], 0.5, 2)
  expect_between(summ$beta1[, "Mean"], 0.5, 2)
  expect_between(summ$beta2[, "Mean"], 2*0.5, 2*2)
  expect_between(summ$beta3[, "Mean"], 3*0.5, 3*2)
  expect_between(summ$beta4[, "Mean"], 4*0.5, 4*2)
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
})

test_that("custom blocking works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
      reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"),
    data=df, control=sampler_control(block=list(c("beta4", "beta2"), c("mu", "beta3")))
  )
  expect_length(sampler$mbs, 2L)
  expect_identical(sort(names(sampler$mbs[[1]]$mcs)), c("beta2", "beta4"))
  expect_identical(sort(names(sampler$mbs[[2]]$mcs)), c("beta3", "mu"))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$mu[, "Mean"], 0.5, 2)
  expect_between(summ$beta1[, "Mean"], 0.5, 2)
  expect_between(summ$beta2[, "Mean"], 2*0.5, 2*2)
  expect_between(summ$beta3[, "Mean"], 3*0.5, 3*2)
  expect_between(summ$beta4[, "Mean"], 4*0.5, 4*2)
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
})

test_that("using an offset for linear regression works", {
  sampler <- create_sampler(
    formula = y ~ reg(~ x1 + x2 + x3, name="beta") + offset(4*x4), data=df
  )
  expect_length(sampler$mod, 2L)
  expect_equal(sampler$mod[[length(sampler$mod)]]$offset, 4*df$x4)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5 * c(1,1,2,3), 2 * c(1,1,2,3))
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
})

test_that("simplified formula + offset works for linear regression", {
  sampler <- create_sampler(
    formula = y ~ x1 + x2 + x3 + offset(4*x4), data=df
  )
  expect_length(sampler$mod, 2L)
  expect_equal(sampler$mod[[length(sampler$mod)]]$offset, 4*df$x4)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], 0.5 * c(1,1,2,3), 2 * c(1,1,2,3))
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
})

test_that("model specification in two steps works", {
  mod <- y ~ x1 + x2 + x3 + x4
  ml_mod <- y ~ reg(mod, prior=pr_normal(precision=1e-6), name="beta")
  sampler <- create_sampler(ml_mod, data=df)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5 * c(1,1,2,3,4), 2 * c(1,1,2,3,4))
  expect_between(summ$sigma_[, "Mean"], 0.5, 2)
})

test_that("non-zero prior mean works", {
  mod <- y ~ x1 + x2 + x3 + x4
  ml_mod <- y ~ reg(mod, prior=pr_normal(mean=c(1,1,2,3,4), precision=1e3), name="beta")
  sampler <- create_sampler(ml_mod, data=df)
  expect_length(sampler$mod[[1L]]$prior[["mean"]], 5L)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.75 * c(1,1,2,3,4), 1.5 * c(1,1,2,3,4))
  expect_between(summ$sigma_[, "Mean"], 0.75, 1.5)
})

test_that("fixed prior mean works", {
  mod <- y ~ x1 + x2 + x3 + x4
  f <- c(1,1,2,3,4)
  ml_mod <- y ~ reg(mod, prior=pr_fixed(value=f), name="beta")
  sampler <- create_sampler(ml_mod, data=df)
  sim <- MCMCsim(sampler, n.iter=10, burnin=0, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$beta[, "q0.05"]), f)
  expect_equal(unname(generate_data(ml_mod, data=df)$pars$beta), f)
})

n <- 1000
sd0 <- 0.41
y <- 1 + rnorm(n, sd=sd0)
test_that("different priors for residual variance work", {
  expect_warning(
    sampler <- create_sampler(y ~ 1, sigma.mod=pr_fixed(sd0^2)),
    "deprecated"
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_between(pm_sigma, 0.75*sd0, 1.25*sd0)
  expect_warning(
    sampler <- create_sampler(y ~ 1, sigma.mod=pr_invchisq(df=1, scale="modeled")),
    "deprecated"
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_between(pm_sigma, 0.75*sd0, 1.25*sd0)
  sampler <- create_sampler(y ~ 1, family=f_gaussian(var.prior = pr_exp(scale=1)))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_between(pm_sigma, 0.75*sd0, 1.25*sd0)
  sampler <- create_sampler(y ~ 1, family=f_gaussian(var.prior=pr_gig(a=2, b=1, p=1.4)))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_between(pm_sigma, 0.75*sd0, 1.25*sd0)
})

test_that("redundancy in design matrix is handled well in estimation and prediction", {
  n <- 500L
  dat <- data.frame(
    age5 = factor(sample(5, n, replace=TRUE)),
    sex = factor(sample(c("M", "F"), n, replace=TRUE))
  )
  dat$age3 <- dat$age5
  levels(dat$age3) <- c("1", "1", "2", "3", "3")
  gd <- generate_data(
    ~ reg(~ age5 + sex*age3, name="beta", prior=pr_normal(prec=1)),
    data=dat
  )
  expect_error(
    sampler <- create_sampler(
      gd$y ~ reg(~ age5 + sex*age3, name="beta"),
      data=dat
    ), "Non-positive-definite matrix"
  )
  sampler <- create_sampler(
    gd$y ~ reg(~ age5 + sex*age3, name="beta", prior=pr_normal(precision=1)),
    data=dat
  )
  expect_identical(sampler$mod[["beta"]]$q, 10L)
  sampler <- create_sampler(
    gd$y ~ reg(~ age5 + sex*age3, name="beta", remove.redundant=TRUE),
    data=dat
  )
  expect_lte(sampler$mod[["beta"]]$q, 8L)
  sim <- MCMCsim(sampler, verbose = FALSE)
  summ <- summary(sim)
  pred0 <- predict(sim, show.progress=FALSE)
  summpred0 <- summary(pred0)
  pred <- predict(sim, newdata = dat, show.progress=FALSE)
  summpred <- summary(pred)
  all.equal(summpred0[, "Mean"], summpred[, "Mean"], tol=0.4)
})
