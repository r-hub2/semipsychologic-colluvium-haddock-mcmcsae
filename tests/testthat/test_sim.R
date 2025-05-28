
context("MCMC simulation")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

ex <- mcmcsae_example()

test_that("MCMCsim arguments work", {
  expect_error(
    sampler <- create_sampler(
      y ~ reg(~x, name="beta") + gen(~x, factor = ~ iid(fA), name="v"),
      linpred = list(beta=matrix(0, 5, 2), v=Matrix(0, 4, 2*nlevels(ex$dat$fA))),
      data=ex$dat), "same number of rows"
  )
  sampler <- create_sampler(
    y ~ reg(~x, name="beta") + gen(~x, factor = ~ iid(fA), name="v"),
    linpred = list(beta=matrix(0, 5, 2), v=Matrix(0, 5, 2*nlevels(ex$dat$fA))),
    data=ex$dat
  )
  # or linear predictor composed of any subset of model components
  sampler <- create_sampler(
    y ~ reg(~x, name="beta") + gen(~x, factor = ~ iid(fA), name="v"),
    linpred = list(v=Diagonal(2*nlevels(ex$dat$fA))),
    data=ex$dat
  )
  sim <- MCMCsim(
    sampler,
    burnin=200, n.iter=500, n.chain=4, thin=2,
    pred = list(v_var="v_sigma^2"),
    trace.convergence=c("v_sigma", "beta", "v[1:3]"),
    verbose=FALSE
  )
  expect_is(sim$linpred_, "dc")
  expect_identical(dim(as.matrix(sim$linpred_)), c(1000L, 2L*nlevels(ex$dat$fA)))
  expect_equal(sim$linpred_, sim$v, check.attributes=FALSE)
  expect_is(sim$v_var, "dc")
  expect_length(sim$v_var, 4)
  expect_identical(dim(sim$v_var[[1]]), c(250L, 2L))
  expect_identical(dim(as.matrix(sim$v_var)), c(1000L, 2L))
  da <- to_draws_array(sim)
  expect_is(da, "draws_array")
  expect_equal(nrow(summary(da)), sum(sapply(sim[par_names(sim)], n_vars)))
  fv <- fitted(sim)
  expect_is(fv, "dc")
  res <- residuals(sim)
  expect_is(res, "dc")
  res <- residuals(sim, matrix=TRUE)
  expect_is(res, "matrix")
  expect_identical(dim(res), c(1000L, nrow(ex$dat)))
})

test_that("combine_chains works", {
  sampler <- create_sampler(
    y ~ reg(~x, name="beta") + gen(~x, factor = ~ iid(fA), name="v"),
    linpred = list(beta=matrix(0, 5, 2), v=Matrix(0, 5, 2*nlevels(ex$dat$fA))),
    data=ex$dat
  )
  sim1 <- MCMCsim(sampler, burnin=10, n.iter=100, store.all=TRUE, verbose=FALSE)
  sim2 <- MCMCsim(sampler, burnin=10, n.iter=100, store.all=TRUE, verbose=FALSE)
  sim <- combine_chains(sim1, sim2)
  expect_identical(n_chains(sim), 6L)
  expect_identical(n_chains(combine_chains_dc(list(sim1$llh_, sim2$llh_))), 6L)
  pred1 <- predict(sim1, ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred2 <- predict(sim2, ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred <- combine_chains_dc(list(pred1, pred2))
  expect_identical(n_chains(pred), 6L)
  expect_length(attr(pred, "ppp"), 100L)
  expect_between(attr(pred, "ppp"), 0, 1)
  f <- function(x) c(mean(x), sd(x))
  pred1 <- predict(sim1, fun.=f, labels=c("mean", "sd"), ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred2 <- predict(sim2, fun.=f, labels=c("mean", "sd"), ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred <- combine_chains_dc(list(pred1, pred2))
  expect_identical(attr(pred, "labels"), c("mean", "sd"))
  f <- function(x) median(x)
  pred1 <- predict(sim1, fun.=f, labels="median", ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred2 <- predict(sim2, fun.=f, labels="median", ppcheck=TRUE, show.progress=FALSE, verbose=FALSE)
  pred <- combine_chains_dc(list(pred1, pred2))
  expect_identical(attr(pred, "labels"), "median")
  pred <- combine_iters_dc(list(pred1, pred2))
})

test_that("Parallel functionality works", {
  skip_if(parallel::detectCores() < 2L, "Skipping test: fewer than 2 cores available")
  sampler <- create_sampler(
    y ~ reg(~x, name="beta") + gen(~x, factor = ~ iid(fA), name="v"),
    linpred = list(beta=matrix(0, 5, 2), v=Matrix(0, 5, 2*nlevels(ex$dat$fA))),
    data=ex$dat
  )
  cl <- setup_cluster(n.cores=2)
  sim <- MCMCsim(sampler, burnin=10, n.iter=100, n.chain=4, store.all=TRUE, cl=cl, verbose=FALSE)
  expect_identical(sim[["_info"]]$n.chain, 4L)
  expect_length(sim[["beta"]], 4L)
  pred0 <- predict(sim, ppcheck=TRUE, n.cores=1L, show.progress=FALSE, verbose=FALSE)
  pred <- predict(sim, ppcheck=TRUE, cl=cl, show.progress=FALSE, verbose=FALSE)
  expect_length(pred, 4L)
  expect_length(attr(pred, "ppp"), 100L)
  expect_between(attr(pred, "ppp"), 0, 1)
  expect_between(mean(attr(pred, "ppp")), 0.8*mean(attr(pred0, "ppp")), 1.2*mean(attr(pred0, "ppp")))
  f <- function(x) sd(x)
  pred0 <- predict(sim, fun.=f, labels="median", ppcheck=TRUE, n.cores=1L, show.progress=FALSE, verbose=FALSE)
  pred <- predict(sim, fun.=f, labels="median", ppcheck=TRUE, cl=cl, show.progress=FALSE, verbose=FALSE)
  expect_length(attr(pred, "ppp"), 1L)
  expect_between(attr(pred, "ppp"), 0.8*attr(pred0, "ppp"), 1.2*attr(pred0, "ppp"))
  stop_cluster(sim[["_cluster"]])
})
