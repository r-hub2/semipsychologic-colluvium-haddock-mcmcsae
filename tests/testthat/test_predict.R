
context("Prediction")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# generate a small example dataset
n <- 100
dat <- data.frame(
  x = rnorm(n),
  f1 = factor(sample(1:5, n, replace=TRUE)),
  f2 = sample(letters[1:6], n, replace=TRUE)
)
v1 <- 0.1*(1:5); v1 <- v1 - mean(v1)
v2 <- setNames(5*runif(6), letters[1:6]); v2 <- v2 - mean(v2)
dat$y <- 1 + 2*dat$x + v1[dat$f1] + v2[dat$f2] + 0.1*rnorm(n)

# and fit a model
sampler <- create_sampler(
  y ~ x + gen(factor = ~f1) + gen(factor = ~f2),
  data=dat
)
sim <- MCMCsim(sampler, store.all=TRUE, n.iter=500, verbose=FALSE)


test_that("prediction is reproducible", {
  pred_data <- predict(sim, seed=1, show.progress=FALSE, verbose=FALSE)
  summ <- summary(pred_data)
  expect_equal(summ, summary(predict(sim, seed=1, show.progress=FALSE, verbose=FALSE)))
  # not specifying newdata implies that the training data is used, i.e. same as newdata=dat
  pred_data <- predict(sim, newdata=dat, seed=1, show.progress=FALSE, verbose=FALSE)
  expect_equal(summ, summary(pred_data))
})

test_that("prediction types 'link' and 'response' work", {
  pred_link <- predict(sim, type="link", show.progress=FALSE, verbose=FALSE)
  pred_response <- predict(sim, type="response", show.progress=FALSE, verbose=FALSE)
  # in this example the link function is identity so they should be equal
  expect_equal(summary(pred_link), summary(pred_response))
})

test_that("prediction for data with fewer model factor levels works", {
  pred <- predict(sim, newdata=dat[1:2, ], type="link", show.progress=FALSE, verbose=FALSE)
  summ <- summary(pred)
  pred <- predict(sim, newdata=dat[1:25, ], type="link", show.progress=FALSE, verbose=FALSE)
  expect_equal(summ[1:2, ], summary(pred)[1:2, ])
})

test_that("predict with on-the-fly aggregation works", {
  fun <- function(x) c(sum(x > 2), x[1] < 0)
  pred <- predict(sim, newdata=dat, fun.=fun, show.progress=FALSE, verbose=FALSE)
  expect_identical(ncol(pred[[1]]), 2L)
  expect_identical(storage.mode(pred[[1]]), "integer")
  fun <- function(x) sum(x)
  pred <- predict(sim, newdata=dat, type="link", fun.=fun, show.progress=FALSE, verbose=FALSE)
  expect_identical(ncol(pred[[1]]), 1L)
  expect_identical(storage.mode(pred[[1]]), "double")
  # posterior predictive checks
  fun <- function(x) c(sum = sum(x > 2), neg_x1 = x[1] < 0, tapply(x, dat$f1, mean))
  pred <- predict(sim, newdata=dat, fun.=fun, ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  ppp <- attr(pred, "ppp")
  expect_length(ppp, 7L)
  expect_between(ppp, 0,1)
  # test statistic that also depends on parameters
  fun <- function(x, p) c(p[["gen2_sigma"]] * sum(x > 1), x[1], p[["reg1"]])
  pred <- predict(sim, newdata=dat, fun.=fun, ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  ppp <- attr(pred, "ppp")
  expect_length(ppp, 4L)
  expect_between(ppp, 0, 1)
})

test_that("predict generates out-of-sample random effects", {
  # all levels are different
  newdat <- dat
  newdat$f2 <- sample(LETTERS, nrow(newdat), replace=TRUE)
  pred <- predict(sim, newdata=newdat, show.progress=FALSE, verbose=FALSE)
  ypred <- summary(pred)[, "Mean"]
  ypred2 <- ypred + v2[dat$f2]
  expect_true(cor(dat$y, ypred) < cor(dat$y, ypred2) && cor(dat$y, ypred2) > 0.5)
  # some levels are different
  newdat <- dat
  ind_new <- 10:50
  newdat$f2[ind_new] <- sample(LETTERS[1:5], length(ind_new), replace=TRUE) 
  pred <- predict(sim, newdata=newdat, show.progress=FALSE, verbose=FALSE)
  ypred <- summary(pred)[, "Mean"]
  ypred2 <- ypred; ypred2[ind_new] <- ypred2[ind_new] + v2[dat$f2[ind_new]]
  expect_true(cor(dat$y, ypred) < cor(dat$y, ypred2) && cor(dat$y, ypred2) > 0.5)
})

# fit a model including categorical fixed effects
sampler <- create_sampler(
  y ~ x + f1 + gen(factor = ~f2),
  data=dat
)
sim <- MCMCsim(sampler, store.all=TRUE, n.iter=300, verbose=FALSE)

test_that("predict indicates when not all model matrix columns for fixed effects are present", {
  newdat <- droplevels(dat[1:2, ])
  expect_error(predict(sim, newdata=newdat), "columns are missing")
})
