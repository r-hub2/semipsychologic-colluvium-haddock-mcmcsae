
context("Computation of weights")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 100L  # sample size
l <- 10L   # number of levels of factor variable
df <- data.frame(
  x = rnorm(n),
  f = sample(seq_len(l), n, replace=TRUE)
)
df$y <- with(df, 1 + x + rnorm(l)[f] + rnorm(n))

test_that("Weights computation is correct", {
  # first linear regression without f effects
  sampler <- create_sampler(y ~ reg(~ 1 + x, name="beta"),
    data=df, linpred=list(beta=matrix(c(1000, 10), nrow=1)),
    compute.weights=TRUE
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  w <- weights(sim)
  expect_equal(w, get_means(sim, "weights_")[[1]])
  # to be 'sure' use equality up to 5 sigmas
  expect_between(sum(w * df$y),
    summ$linpred_[, "Mean"] - 5 * summ$linpred_[, "MCSE"],
    summ$linpred_[, "Mean"] + 5 * summ$linpred_[, "MCSE"]
  )
  # random effects model
  sampler <- create_sampler(y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ iid(f), name="v"),
    data=df, linpred=list(beta=matrix(c(1000, 10), nrow=1), v=matrix(0, nrow=1, ncol=l)),
    compute.weights=TRUE
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  w <- weights(sim)
  expect_equal(w, get_means(sim, "weights_")[[1]])
  # to be 'sure' use equality up to 5 sigmas
  expect_between(sum(w * df$y),
    summ$linpred_[, "Mean"] - 5 * summ$linpred_[, "MCSE"],
    summ$linpred_[, "Mean"] + 5 * summ$linpred_[, "MCSE"]
  )
})
