
context("Multinomial regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 2000
dat <- data.frame(
  x = rnorm(n)
)
K <- 5  # number of categories

test_that("categorical data generation, model fitting and prediction work", {
  gd <- generate_data(
    ~ reg(~ 1 + x, prior=pr_fixed(c(-1, 1))),
    family=f_multinomial(K=K),
    data=dat
  )
  expect_equal(gd$pars$reg1, setNames(c(-1, 1), c("(Intercept)", "x")))
  expect_true(is.integer(gd$y))
  expect_true(all(gd$y %in% 1:5))
  dat$y <- gd$y
  s <- 1:1000  # subset to use for model fitting
  sampler <- create_sampler(y ~ x, family="multinomial", data=dat[s, , drop=FALSE])
  n.iter <- 400
  sim <- MCMCsim(sampler, burnin=150, n.iter=n.iter, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], c(-2, 0.5), c(-0.5, 2))
  pred <- predict(sim, show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(dim(pred[[1L]]), c(n.iter, (K-1)*length(s)))
  sn <- 1001:2000
  pred <- predict(sim, newdata=dat[sn, , drop=FALSE], show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(nrow(summ), c((K-1)*length(sn)))
  # predict at response scale without multinomial noise:
  pred0 <- predict(sim, newdata=dat[sn, , drop=FALSE], type="response", show.progress=FALSE)
  summ0 <- summary(pred0)
  expect_between(summ0[, "Mean"], 0.5*summ[, "Mean"] - 0.01, 2*summ[, "Mean"] + 0.01)
  expect_true(all(summ0[, "SD"] < summ[, "SD"] + 0.01))
  expect_true(all(rowSums(matrix(summ0[, "Mean"], ncol=K-1)) <= 1))
})

test_that("multinomial data generation, model fitting and prediction work, for different specifications of the response variable", {
  ny <- 10
  K <- 5
  gd <- generate_data(
    ~ reg(~ 1 + x, prior=pr_fixed(c(-1, 1))),
    family=f_multinomial(n.trial = ny, K = K),
    data=dat
  )
  expect_true(is.integer(gd$y))
  expect_equal(dim(gd$y), c(nrow(dat), K))
  expect_true(all(rowSums(gd$y) == ny))
  dat$y <- gd$y
  sampler <- create_sampler(y ~ x, data=dat, family="multinomial")
  expect_equal(sampler$Km1, K - 1L)
  expect_equal(sampler$family$ny0, rep.int(ny, nrow(dat)))
  sim <- MCMCsim(sampler, n.chain=2, burnin=50, n.iter=200, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], c(-2, 0.5), c(-0.5, 2))
  pred <- predict(sim, iters=1, show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(dim(pred[[1]]), c(1L, nrow(dat)*(K - 1L)))
  expect_true(all(rowSums(matrix(pred[[1]][1, ], ncol = K - 1L)) <= ny))
  pred <- predict(sim, iters=2, weights = 2, show.progress=FALSE)
  expect_equal(dim(pred[[1]]), c(1L, nrow(dat)*(K - 1L)))
  expect_true(all(rowSums(matrix(pred[[1]][1, ], ncol = K - 1L)) <= 2*ny))
  # fraction data; remove last column
  dat$y <- dat$y[, -ncol(dat$y)]/rowSums(dat$y)
  sampler <- create_sampler(y ~ x, data=dat, family=f_multinomial(n.trial = ny))
  expect_equal(sampler$Km1, K - 1L)
  expect_equal(sampler$family$ny0, ny)
  sim <- MCMCsim(sampler, n.chain=2, burnin=50, n.iter=200, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[, "Mean"], c(-2, 0.5), c(-0.5, 2))
  pred <- predict(sim, iters=1, show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(dim(pred[[1]]), c(1L, nrow(dat)*(K - 1L)))
  expect_true(all(rowSums(matrix(pred[[1]][1, ], ncol = K - 1L)) <= ny))
  pred <- predict(sim, iters=2, weights = 2, show.progress=FALSE)
  expect_equal(dim(pred[[1]]), c(1L, nrow(dat)*(K - 1L)))
  expect_true(all(rowSums(matrix(pred[[1]][1, ], ncol = K - 1L)) <= 2*ny))
})

n <- 200
X <- data.frame(
  int = rep.int(1, n),
  x = runif(n),
  z = rnorm(n)
)
y <- factor(sample(letters[1:3], n, replace=TRUE))
test_that("multinomial (categorical) model with factor response runs", {
  sampler <- create_sampler(y ~ reg(~ 0 + cat_ / X[, -1]), family="multinomial")
  expect_identical(sampler$family$ny0, 1L)
  sim <- MCMCsim(sampler, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_identical(nrow(summ$reg1), 6L)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_equal(DIC[["DIC"]], WAIC[["WAIC1"]], tol=0.3)
  expect_equal(DIC[["DIC"]], WAIC[["WAIC2"]], tol=0.3)
})

n <- 500
dat <- data.frame(
  x = rnorm(n),
  f = factor(sample(1:10, n, replace=TRUE))
)
K <- 3
test_that("multinomial (categorical) random effects model works", {
  gd <- generate_data(
    ~ reg(~ 0 + cat_ / x, prior=pr_normal(precision=1)) +
      gen(~ 0 + cat_, var="diagonal", factor = ~ f),
    family=f_multinomial(n.trial = 1, K = K),
    data=dat
  )
  dat$y <- gd$y
  sampler <- create_sampler(
    y ~ 0 + cat_ / x + gen(~ 0 + cat_, var="diagonal", factor = ~ f),
    family="multinomial", data=dat
  )
  sim <- MCMCsim(sampler, n.chain=2, burnin=100, n.iter=600, verbose=FALSE, store.all=TRUE)
  summ <- summary(sim)
  expect_between(gd$pars$reg1, summ$reg1[, "Mean"] - 3 * summ$reg1[, "SD"], summ$reg1[, "Mean"] + 3 * summ$reg1[, "SD"])
  pred <- predict(sim, newdata=dat[1:10, ], show.progress=FALSE, verbose=FALSE)
  summpred <- summary(pred)
  expect_equal(nrow(summpred), (K - 1)*10)
  expect_true(
    all(summpred[1:10, "Mean"] + summpred[11:20, "Mean"] >= 0) &&
    all(summpred[1:10, "Mean"] + summpred[11:20, "Mean"] <= 1)
  )
  # use of cat_ in factor
  sampler <- create_sampler(
    y ~ 0 + cat_ / x + gen(factor = ~ f * cat_),
    family="multinomial", data=dat
  )
  sim <- MCMCsim(sampler, n.chain=2, burnin=100, n.iter=100, verbose=FALSE, store.all=TRUE)
  summ <- summary(sim)
  expect_equal(nrow(summ$gen2), (K - 1)*10)
  pred <- predict(sim, newdata=dat[1:10, ], show.progress=FALSE, verbose=FALSE)
  summpred <- summary(pred)
  expect_equal(nrow(summpred), (K - 1)*10)
  expect_true(
    all(summpred[1:10, "Mean"] + summpred[11:20, "Mean"] >= 0) &&
    all(summpred[1:10, "Mean"] + summpred[11:20, "Mean"] <= 1)
  )
})
