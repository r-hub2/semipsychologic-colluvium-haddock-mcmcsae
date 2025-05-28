
context("Variance modelling")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 250L
df <- data.frame(x=runif(n), f=factor(sample(1:8, n, replace=TRUE)))
Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=5))
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
                     family=f_gaussian(var.prior=1, var.model = Vmodel), data=df)

test_that("generated data based on vfac model is OK", {
  expect_length(dat$y, n)
  expect_false(anyNA(dat$y))
})

test_that("modelling vfac variance structure works", {
  df$y <- dat$y
  Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=2), name="varf")
  expect_warning(
    sampler <- create_sampler(y ~ x + f, sigma.fixed=TRUE, formula.V=Vmodel, data=df),
    "deprecated"
  )
  sampler <- create_sampler(y ~ x + f, family=f_gaussian(var.prior=1, var.model=Vmodel), data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  compute_DIC(sim)
  compute_WAIC(sim)
})

test_that("pr_fixed prior in vfac works", {
  Vmodel <- ~ vfac(factor="f", prior=pr_fixed(value=1:8))
  dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
                       family=f_gaussian(var.prior=1, var.model=Vmodel), data=df)
  expect_equal(unname(dat$pars$vfac1), 1:8)
  expect_length(dat$y, n)
  expect_false(anyNA(dat$y))
  df$y <- dat$y
  sampler <- create_sampler(
    y ~ reg(~ x + f),
    family = f_gaussian(var.prior = 1, var.model =  ~ vfac(factor="f")),
    data=df
  )
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$vfac1[, "Mean"], 0.5 * (1:8) - 0.5, 2*(1:8) + 1)
  sampler <- create_sampler(
    y ~ reg(~ x + f),
    family=f_gaussian(var.prior=1, var.model = ~ vfac(factor="f", prior=pr_fixed(1:8))),
    data=df
  )
  sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$vfac1[, "q0.05"]), 1:8)
})

# non-diagonal sampling covariance matrix (NB very contrived example)
subdiv <- c(110,10,70,45,8,7)
df$f <- as.factor(rep(1:6, subdiv))
Q0 <- bdiag(mapply(Q_RW2, subdiv)) + Diagonal(n)
Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=5))
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
         family=f_gaussian(var.prior=1, prec.mat=Q0, var.model=Vmodel), data=df)

test_that("vfac variance structure works for compatible non-diagonal sampling variance matrix", {
  df$y <- dat$y
  Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=2), name="varf")
  expect_warning(
    sampler <- create_sampler(y ~ x + f, sigma.fixed=TRUE, Q0=Q0, formula.V=Vmodel, data=df),
    "deprecated"
  )
  sampler <- create_sampler(y ~ x + f,
    family = f_gaussian(var.prior=1, prec.mat=Q0, var.model=Vmodel), data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  compute_DIC(sim)
  compute_WAIC(sim)
})


test_that("generated data based on vreg model is OK", {
  Vmodel <- ~ vreg(~ x, prior=pr_normal(precision=1))
  dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
    family=f_gaussian(var.prior=pr_fixed(1), var.model=Vmodel), data=df)
  expect_length(dat$y, n)
  expect_false(anyNA(dat$y))
})

Vmodel <- ~ vreg(~ x, prior=pr_fixed(1))
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
                     family=f_gaussian(var.prior=1, var.model=Vmodel), data=df)
df$y <- dat$y
test_that("vreg variance model works", {
  expect_equal(dat$pars$vreg1, c(`(Intercept)` = 1, x = 1))
  Vmodel <- ~ vreg(formula = ~ x, name="varf")
  sampler <- create_sampler(y ~ x + f, family=f_gaussian(var.prior=1, var.model=Vmodel), data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  expect_between(summ$varf[, "Mean"], 0.4, 2.5)
  compute_DIC(sim)
  compute_WAIC(sim)
})

n <- 900
dat <- data.frame(
  x = runif(n),
  z = rnorm(n)
)
dat$y <- rnorm(n, 1 + dat$x + dat$z, sd=exp(0.5*(0.5 - 0.5*dat$x)))
test_that("reg variance model works", {
  sampler <- create_sampler(y ~
    reg(~ x + z, name="beta"),
    family=f_gaussian(var.prior=1, var.model = ~ reg(~ x, name="vbeta")),
    data=dat
  )
  sim <- MCMCsim(sampler, burnin=150, n.iter=500, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$vbeta, "dc")
  summ <- summary(sim)
  expect_between(summ$beta[, "q0.5"], 0.3 * c(1, 1, 1), 3 * c(1, 1, 1))
  expect_between(summ$vbeta[, "q0.5"], c(0.5 * 0.3, -0.5 * 3), c(0.5 * 3, -0.5 * 0.3))
  compute_DIC(sim)
  compute_WAIC(sim)
})

n <- 1200
dat <- data.frame(
  x = runif(n),
  z = rnorm(n),
  f = factor(sample(22, n, replace=TRUE))
)
test_that("mixed effects for both mean and variance work", {
  gd <- generate_data(
    ~ reg(~ 1 + x + z, prior=pr_fixed(c(1,0.5,-0.3)), name="beta") +
      gen(factor = ~ f, PX=FALSE, prior=pr_invchisq(1e6, 0.5), name="v"),
    family = f_gaussian(
      var.model = ~ reg(~ 1 + x, prior=pr_fixed(c(0.4,1)), name="vbeta") +
                    gen(factor = ~ f, prior=pr_invchisq(1e6, 1), name="vv"),
      var.prior = pr_fixed(1)
    ), data=dat
  )
  expect_equal(unname(gd$pars$beta), c(1,0.5,-0.3))
  expect_equal(unname(gd$pars$vbeta), c(0.4,1))
  expect_between(gd$pars$v_sigma, 0.9*sqrt(0.5), 1.1*sqrt(0.5))
  expect_between(gd$pars$vv_sigma, 0.9*sqrt(1), 1.1*sqrt(1))
  sampler <- create_sampler(
    gd$y ~ reg(~ 1 + x + z, name="beta") +
        gen(factor = ~ f, name="v"),
    family = f_gaussian(
      var.model = ~ reg(~ 1 + x, name="vbeta") +
                    gen(factor = ~ f, name="vv"),
      var.prior = pr_fixed(1)
    ), data=dat
  )
  expect_false(sampler$Vmod$vv$usePX)
  sim <- MCMCsim(sampler, burnin=200, n.iter=600, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "q0.5"], c(1*0.3, 0.5*0.3, -0.3*3), c(1*3, 0.5*3, -0.3*0.3))
  expect_between(summ$vbeta[, "q0.5"], 0.2*c(0.4, 1), 5*c(0.4, 1))
  #plot(gd$pars$v, summ$v[, "Mean"]); abline(0, 1)
  #plot(gd$pars$vv, summ$vv[, "Mean"]); abline(0, 1)
  #mean(summ$vv[, "Mean"])
  compute_DIC(sim)
  compute_WAIC(sim)
  # test with cMVN sampler
  sampler <- create_sampler(
    gd$y ~ reg(~ 1 + x + z, name="beta") +
      gen(factor = ~ f, name="v"),
    family = f_gaussian(
      var.model = ~ reg(~ 1 + x, name="vbeta") +
        gen(factor = ~ f, name="vv"),
      var.prior = pr_fixed(1)
    ), data=dat,
    control = sampler_control(expanded.cMVN.sampler = TRUE)
  )
  expect_true(sampler$control$expanded.cMVN.sampler)
  sim <- MCMCsim(sampler, burnin=200, n.iter=600, n.chain=2, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$beta[, "q0.5"], c(1*0.3, 0.5*0.3, -0.3*3), c(1*3, 0.5*3, -0.3*0.3))
  expect_between(summ$vbeta[, "q0.5"], 0.2*c(0.4, 1), 5*c(0.4, 1))
})
