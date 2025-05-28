
context("GMRF extensions")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 400
dat <- data.frame(
  x=rnorm(n),
  t=1:n
)
gd <- generate_data(
  ~ reg(~ x, name="beta", prior=pr_normal(0, 1)) +
    gen(factor = ~ RW1(t), prior=pr_invchisq(10, 2), name="v",
        strucA=GMRF_structure(type="leroux", prior=0.9)),
  family=f_gaussian(var.prior = pr_fixed(1)), data=dat
)
dat$y <- gd$y
#plot.ts(dat$y)

test_that("Fitting a model with a Leroux-extended random walk works", {
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(type="Leroux")),
    family = f_gaussian(var.prior = pr_fixed(1)), data=dat
  )
  # both "Leroux" and "leroux" should work
  expect_identical(sampler$mod$v$strucA$type, "leroux")
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
        gen(factor = ~ RW1(t), name="v",
            strucA=GMRF_structure(type="leroux")),
    family = f_gaussian(var.prior = pr_fixed(1)), data=dat
  )
  expect_identical(sampler$mod$v$strucA$prior$type, "unif")
  expect_identical(sampler$mod$v$strucA$control$type, "RWTN")
  expect_true(sampler$mod$v$strucA$control$adaptive)
  expect_identical(sampler$MHpars, "v_GMRFext")
  sim <- MCMCsim(sampler, store.all=TRUE, n.chain=2, burnin=300, n.iter=800, verbose=FALSE)
  summ <- summary(sim)
  expect_between(summ$v_GMRFext[, "Mean"], 0.4, 1)
})

test_that("Scaling of precision matrix works", {
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(scale.precision = FALSE)),
    family = f_gaussian(var.prior = 1), data=dat
  )
  expect_equal(sampler$family$var.prior$value, 1)
  expect_equal(sampler$mod$v$QA@x[1], 1)
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
        gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(scale.precision = TRUE)),
    family = f_gaussian(var.prior = 1), data=dat
  )
  expect_between(sampler$mod$v$QA@x[1], 40, 75)
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(type="leroux", prior=0.5, scale.precision = TRUE)),
    family = f_gaussian(var.prior = 1), data=dat
  )
  expect_between(sampler$mod$v$QA@x[1], 20, 40)
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(type="leroux", scale.precision = TRUE)),
    family = f_gaussian(var.prior = 1), data=dat
  )
  expect_between(sampler$mod$v$QA@x[1], 40, 75)
})

test_that("bym2 extension with fixed parameter works", {
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(type="bym2", prior=0.5)),
    family=f_gaussian(var.prior=1), data=dat
  )
  expect_identical(sampler$mod$v$strucA$type, "bym2")
  expect_identical(nrow(sampler$mod$v$QA), 800L)
  expect_identical(dim(sampler$mbs[[1]]$R), c(802L, 1L))
  sim <- MCMCsim(sampler, n.chain=1, burnin=40, n.iter=60, store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_identical(nrow(summ$v), 800L)
  #ts.plot(summ$v[1:400, "Mean"])    # total effect
  #ts.plot(summ$v[401:800, "Mean"])  # structured effect
})

test_that("bym2 extension with inferred parameter works", {
  gd <- generate_data(
    ~ reg(~ x, name="beta", prior=pr_normal(0, 1)) +
      gen(factor = ~ RW1(t), prior=pr_invchisq(10, 2), name="v",
          strucA=GMRF_structure(type="bym2", prior=0.5)),
    family=f_gaussian(var.prior=1), data=dat
  )
  dat$y <- gd$y
  #plot.ts(dat$y)
  sampler <- create_sampler(
    y ~ reg(~ x, name="beta") +
      gen(factor = ~ RW1(t), name="v",
          strucA=GMRF_structure(type="bym2")),
    family=f_gaussian(var.prior=1), data=dat
  )
  expect_identical(sampler$mod$v$strucA$type, "bym2")
  expect_identical(dim(sampler$mbs[[1]]$R), c(802L, 1L))
  expect_identical(sampler$mod$v$strucA$prior$type, "unif")
  expect_identical(sampler$mod$v$strucA$control$type, "RWTN")
  sim <- MCMCsim(sampler, n.chain=1, burnin=40, n.iter=60, store.all=TRUE, verbose=FALSE)
  compute_DIC(sim)
  summ <- summary(sim)
  expect_true(!is.null(summ$v_GMRFext))
  expect_is(acceptance_rates(sim)$v_GMRFext, "list")
  expect_identical(nrow(summ$v), 800L)
  #ts.plot(summ$v[1:400, "Mean"])    # total effect
  #ts.plot(summ$v[401:800, "Mean"])  # structured effect
})
