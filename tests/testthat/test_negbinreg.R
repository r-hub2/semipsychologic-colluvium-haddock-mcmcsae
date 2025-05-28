
context("Negative binomial regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# generate negative binomial data
n <- 1000L
r <- 0.17  # shape or inverse dispersion parameter
df <- data.frame(x1=rnorm(n), x2=runif(n))
mod <- ~ reg(~ 1 + x1 + x2, prior=pr_fixed(c(3, 1, 0.5)), name="beta")
dat <- generate_data(mod, family=f_negbinomial(inv.shape.prior = pr_fixed(value=1/r)), data=df)

test_that("fitting a negative binomial model works", {
  expect_error(create_sampler(dat$y ~ reg(~ 1 + x1 + x2, name="beta"), data=df,
    family = "negbinomial", ry=r), "no longer supported")
  expect_error(create_sampler(dat$y ~ reg(~ 1 + x1 + x2, name="beta"), data=df,
    family = "negbinomial", r.mod=pr_fixed(1/r)), "no longer supported")
  sampler <- create_sampler(dat$y ~ reg(~ 1 + x1 + x2, name="beta"), data=df,
                            family = f_negbinomial(inv.shape.prior = pr_fixed(1/r)))
  sim <- MCMCsim(sampler, n.iter=400L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.4*dat$pars$beta, 2.5*dat$pars$beta)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
  pred <- as.matrix(predict(sim, type="response", iters=1:100, show.progress=FALSE))
  expect_equivalent(r*mean(pred), mean(dat$y), tol=2)
})

test_that("fitting a negative binomial model with scaled beta prime prior on dispersion parameter works", {
  sampler <- create_sampler(dat$y ~ x1 + x2, data=df,
                            family = f_negbinomial(inv.shape.prior = pr_invchisq(df=1, scale="modeled")))
  #summary(replicate(100, sampler$rprior()$negbin_shape_))
  sim <- MCMCsim(sampler, n.iter=400L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_true(abs(summ$negbin_shape_[, "Mean"] - r) < 0.05)
  expect_between(summ$beta[, "Mean"], 0.4*dat$pars$beta, 2.5*dat$pars$beta)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
})

test_that("fitting a negative binomial model with chi-squared prior on dispersion parameter works", {
  sampler <- create_sampler(dat$y ~ x1 + x2, data=df,
                            family = f_negbinomial(inv.shape.prior = pr_invchisq(df=1, scale=1)))
  #summary(replicate(100, sampler$rprior()$negbin_shape_))
  sim <- MCMCsim(sampler, n.iter=400L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_true(abs(summ$negbin_shape_[, "Mean"] - r) < 0.05)
  expect_between(summ$beta[, "Mean"], 0.4*dat$pars$beta, 2.5*dat$pars$beta)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
})

test_that("fitting a negative binomial model with GIG prior on dispersion parameter works", {
  sampler <- create_sampler(dat$y ~ x1 + x2, data=df,
                            family = f_negbinomial(inv.shape.prior=pr_gig(a=0, b=1, p=-1/2)))
  #summary(replicate(100, sampler$rprior()$negbin_shape_))
  sim <- MCMCsim(sampler, n.iter=400L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.4*dat$pars$beta, 2.5*dat$pars$beta)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
})

test_that("fitting a negative binomial model with dispersion parameter fixed by prior works", {
  sampler <- create_sampler(dat$y ~ x1 + x2, data=df,
                            family = f_negbinomial(inv.shape.prior = pr_fixed(value=1/r)))
  expect_equal(sampler$family$shape0, r)
  expect_true(all(replicate(10, sampler$family$inv.shape.prior$rprior()) == 1/r))
  sim <- MCMCsim(sampler, n.iter=400L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "Mean"], 0.4*dat$pars$beta, 2.5*dat$pars$beta)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
})


# generate negative binomial data
n <- 1000L
r <- 1  # dispersion parameter
df <- data.frame(x1=rnorm(n), x2=runif(n), f=factor(sample(1:25, n, replace=TRUE)))
mod <- ~ reg(~ 1 + x1 + x2, prior=pr_fixed(c(3, 1, 0.5)), name="beta") +
  gen(factor=~f, prior=pr_invchisq(df=1000, scale=0.1), name="v")
ry <- c(rep.int(1, 500), rep.int(5, 500))
dat <- generate_data(mod, data=df,
                     family = f_negbinomial(shape.vec=~ry, inv.shape.prior=pr_fixed(1/r)))
test_that("prediction works well for negative binomial data", {
  s <- 1:600
  sampler <- create_sampler(dat$y[s] ~ x1 + x2 + gen(factor=~f), data=df[s, ],
                            family = f_negbinomial(shape.vec = ~ ry[s]))
  sim <- MCMCsim(sampler, n.iter=500L, burnin=200L, n.chain=2L, store.all=TRUE,
                 start=list(list(negbin_shape_=2), list(negbin_shape_=3)), verbose=FALSE)
  #plot(sim, c("negbin_shape_", "reg1[1]", "gen2_sigma"))
  summ <- summary(sim)
  expect_true(abs(summ$negbin_shape_[, "Mean"] - r) < 0.25)
  s <- 601:1000
  pred <- predict(sim, newdata=df[s, ], fun.=function(x) sum(x), show.progress=FALSE)
  summ <- summary(pred)
  expect_true(abs(summ[, "Mean"] - sum(dat$y[s])) < 30000)
  # or ry could be in data:
  df$y <- dat$y
  df$shapevec <- ry
  s <- 1:600
  sampler <- create_sampler(y ~ x1 + x2 + gen(factor=~f), data=df[s, ],
                            family = f_negbinomial(shape.vec = ~ shapevec))
  sim <- MCMCsim(sampler, n.iter=500L, burnin=200L, n.chain=2L, store.all=TRUE,
                 start=list(list(negbin_shape_=2), list(negbin_shape_=3)), verbose=FALSE)
  summ <- summary(sim)
  expect_true(abs(summ$negbin_shape_[, "Mean"] - r) < 0.25)
  s <- 601:1000
  pred <- predict(sim, newdata=df[s, ], fun.=function(x) sum(x), show.progress=FALSE)
  summ <- summary(pred)
  expect_true(abs(summ[, "Mean"] - sum(dat$y[s])) < 30000)
})
