
context("Bayesian Additive Regression Trees")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 500L
dat <- data.frame(x=runif(n), z=rnorm(n))
fx <- function(x) 1 + 2*x - exp(-5*x)
#plot(dat$x, fx(dat$x))
dat$y <- fx(dat$x) + 0.5*dat$z + rnorm(n, sd=0.01)

testdat <- data.frame(x=runif(10L, -0.4, 1.2), z=0.5)

test_that("BART component 'brt' works", {
  n.tr <- 40L
  sampler <- create_sampler(
    y ~ reg(~ 0 + z, name="beta") +
        brt(~ 0 + x, n.trees=n.tr, debug=FALSE, keepTrees=TRUE, name="bart"),
    data=dat
  )
  n.ch <- 2L
  sim <- MCMCsim(sampler, burnin=400, n.iter=800, n.chain=n.ch, thin=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.5*0.5, 2*0.5)
  expect_identical(nrow(summ$bart), n)
  #plot(dat$x, summ$bart[, "Mean"])
  expect_length(sim$bart_trees_, n.ch)
  expect_length(sim$bart_trees_[[1L]], n_draws(sim))
  expect_equal(max(sim$bart_trees_[[1L]][[1L]]$tree), n.tr)
  compute_DIC(sim)
  waic1 <- compute_WAIC(sim, show.progress=FALSE)
  suppressWarnings(waic2 <- waic(sim))
  expect_equal(waic1[["WAIC2"]], waic2$estimates["waic", "Estimate"])
  expect_equal(waic1[["p_WAIC2"]], waic2$estimates["p_waic", "Estimate"])
  suppressWarnings(loo(sim))

  pred <- predict(sim, newdata=testdat, iters=sample(n_draws(sim), 30), show.progress=FALSE)
  summpred <- summary(pred)
  #plot(testdat$x, summpred[, "Mean"])
  #points(testdat$x, fx(testdat$x) + 0.5*0.5, col="red")

  spl <- split_iters(sim, parts=3)
  expect_length(spl, 3L)
  expect_length(spl[[1]]$bart_trees_, n.ch)
  expect_lt(length(spl[[1]]$bart_trees_[[1L]]), n_draws(sim)/2)
  expect_identical(max(spl[[1]]$bart_trees_[[1L]][[1L]]$tree), n.tr)

  test <- combine_chains(sim, sim)
  expect_length(test$bart_trees_, 2*n.ch)
  expect_length(test$bart_trees_[[1L]], n_draws(sim))
  expect_identical(max(test$bart_trees_[[1L]][[1L]]$tree), n.tr)
})
