context("Constrained MVN sampler")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("sampling from IGMRF seems fine", {
  n <- 20L
  R <- R_RW1(n)
  sampler <- create_cMVN_sampler(D = D_RW1(n), R=R)
  expect_length(sampler$draw(), n)
  expect_equal(crossprod_mv(R, sampler$draw()), 0)
  R <- R_RW2(n)
  sampler <- create_cMVN_sampler(D = D_RW2(n), R=R, eps1=1e-10)
  expect_length(sampler$draw(), n)
  expect_equal(crossprod_mv(R, sampler$draw()), c(0, 0))
})

test_that("estimation of marginal variances of IGMRF priors works", {
  n <- 19L
  D <- D_RW2(n)
  R <- R_RW2(n)
  mv <- sim_marg_var(D, R=R)
  expect_length(mv, n)
})

test_that("constrained MVN sampler can be used within block Gibbs sampler", {
  ex <- mcmcsae_example()
  sampler <- create_sampler(
    y ~ reg(~x, name="beta")
      + gen(~x, factor=~fA, name="v")
      + gen(factor = ~RW2(fT), name = "u")
    , data = ex$dat, control = sampler_control(cMVN.sampler = TRUE)
  )
  sim <- MCMCsim(sampler, n.chain=2, burnin=200, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$sigma_[, "Mean"], 0.5*ex$pars$sigma_, 2*ex$pars$sigma_)
})
