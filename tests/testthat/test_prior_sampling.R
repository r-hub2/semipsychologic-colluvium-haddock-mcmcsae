
context("Generate data")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("prior predictive sampling works", {
  df <- data.frame(int=rep(1, 1000))
  dat <- generate_data(~ reg(~ 1, prior=pr_fixed(2)), sigma.mod=pr_fixed(4), data=df)
  expect_between(mean(dat$y), 1.5, 2.5)
  expect_between(sd(dat$y), 1.5, 2.5)
})

test_that("prior specification of regression component works", {
  library(survey)
  data(api)
  mod <- api00 ~
    reg(~ (api99 + stype) * sch.wide,
      prior=pr_normal(mean=1:8, precision=c(1e3,0.01,1,1e6,100,1,1,1e6))) +
    gen(~ api99, factor= ~ cname)
  samplr <- create_sampler(mod, sigma.mod=pr_invchisq(df=10, scale=1),
    data=apisrs, control=sampler_control(block=FALSE))
  sim <- MCMCsim(samplr, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[c(4,8), "Mean"], 0.99*c(4,8), 1.01*c(4,8))
  # prior sampling
  sim <- MCMCsim(samplr, from.prior=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$reg1[c(4,8), "Mean"], 0.99*c(4,8), 1.01*c(4,8))
})
