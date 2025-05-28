
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
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
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
