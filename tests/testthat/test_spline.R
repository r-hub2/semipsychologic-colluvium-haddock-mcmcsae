
context("Models with splines component")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 200
x <- seq(0, 1, length.out=n)
sde <- 0.1
dat <- data.frame(x = x, y = 1 - x + x*sin(22*x) + rnorm(n, sd=sde))
#plot(dat$x, dat$y)

test_that("data generation for splines model works", {
  knots <- 25L
  degree <- 2L
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") +
        gen(factor = ~ splines(x, knots=knots, degree=degree), name="v"),
    data=dat
  )
  expect_equal(sampler$mod$v$info$factors[[1]]$degree, degree)
  expect_equal(sampler$mod$v$q, knots - degree - 1L)
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ splines(x, knots=30, degree=2), name="v"),
    data=dat
  )
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE, n.chain=2, n.iter=700)
  summ <- summary(sim)
  pred <- predict(sim, type="response", show.progress=FALSE)
  summpred <- summary(pred)
  #plot(dat$x, dat$y)
  #points(dat$x, summpred[, "Mean"], pch=20, col="red")
  expect_gt(cor(dat$y, summpred[, "Mean"]), 0.8)
  expect_lt(mean(abs(dat$y - summpred[, "Mean"])), 2 * sde)
  newdat <- data.frame(x=runif(100, 0.5, 1.5))
  predn <- predict(sim, newdata=newdat, type="response", show.progress=FALSE)
  summpredn <- summary(predn)
  #plot(newdat$x, summpredn[, "Mean"], pch=20, col="red")
})

test_that("spline modelling using mgcv smooths seems to work", {
  sampler <- create_sampler(
    y ~ reg(~ 1, name="beta") + s(x),
    data=dat
  )
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE, n.chain=2, n.iter=700)
  summ <- summary(sim)
  expect_identical(nrow(summ$s2_f), 1L)
  pred <- predict(sim, type="response", show.progress=FALSE)
  summpred <- summary(pred)
  #plot(dat$x, dat$y)
  #points(dat$x, summpred[, "Mean"], pch=20, col="red")
  expect_gt(cor(dat$y, summpred[, "Mean"]), 0.8)
  expect_lt(mean(abs(dat$y - summpred[, "Mean"])), 2 * sde)
  start.fun <- function() list(s2_f = 1)
  sim <- MCMCsim(sampler, start=start.fun, burnin=0, n.iter=1, verbose=FALSE)
  start.fun <- function() list(s2_f = c(NA))
  expect_error(
    sim <- MCMCsim(sampler, start=start.fun, burnin=0, n.iter=1, verbose=FALSE),
    "is NA"
  )
  start.fun <- function() list(s2_f = c(1, 2))
  expect_error(
    sim <- MCMCsim(sampler, start=start.fun, burnin=0, n.iter=1, verbose=FALSE),
    "starting value"
  )
})

n <- 2000
p <- 10
dat <- as.data.frame(matrix(runif(n * p), n, p))
dat$y <- with(dat, 10 * sin(pi * V1 * V2) +
  20 * (V3 - 0.5)^2 + 10 * V4 + 5 * V5) +
  rnorm(n, sd = 1)
test_that("using multiple mgcv smooths works", {
  sampler <- create_sampler(
    y ~  s(V1) + s(V2) + s(V3) + s(V4) + s(V5) + s(V6) +
      poly(I(V1*V2),3),
    data=dat
  )
  expect_length(sampler$control$block[[1]], 7)
  sim <- MCMCsim(sampler, verbose=FALSE, n.chain=2, burnin=40, n.iter=100, store.all=TRUE)
  summary(sim)
  sampler <- create_sampler(
    y ~ s(V1) + s(V2) + s(V3) + s(V4) + s(V5) + s(V6) +
      poly(I(V1*V2),3),
    data=dat,
    control = sampler_control(block=FALSE)
  )
  expect_length(sampler$control$block, 0)
  sim <- MCMCsim(sampler, verbose=FALSE, n.chain=2, burnin=0, n.iter=10, store.all=TRUE)
})
