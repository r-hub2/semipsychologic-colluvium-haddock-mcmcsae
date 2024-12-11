
context("Multinomial regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 2000
dat <- data.frame(
  x = rnorm(n)
)
K <- 5  # number of categories

test_that("categorical data generation, fitting and prediction work", {
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
  expect_equal(summ$reg1["(Intercept)", "Mean"], -1, tolerance=0.25)
  pred <- predict(sim, show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(dim(pred[[1L]]), c(n.iter, (K-1)*length(s)))
  sn <- 1001:2000
  pred <- predict(sim, newdata=dat[sn, , drop=FALSE], show.progress=FALSE)
  summ <- summary(pred)
  expect_equal(nrow(summ), c((K-1)*length(sn)))
})

test_that("multinomial data generation works", {
  ny <- 10
  gd <- generate_data(
    ~ reg(~ 1 + x, prior=pr_fixed(c(-1, 1))),
    family=f_multinomial(K=K), ny=ny,
    data=dat
  )
  expect_true(is.integer(gd$y))
  expect_equal(dim(gd$y), c(nrow(dat), K))
  expect_true(all(rowSums(gd$y) == ny))
})

n <- 200
X <- data.frame(
  int = rep.int(1, n),
  x = runif(n),
  z = rnorm(n)
)
y <- factor(sample(letters[1:3], n, replace=TRUE))
test_that("multinomial model with unconventional input data runs", {
  sampler <- create_sampler(y ~ reg(~ 0 + cat_ / X[, -1]), family="multinomial")
  sim <- MCMCsim(sampler, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_identical(nrow(summ$reg1), 6L)
})
