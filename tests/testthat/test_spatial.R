
context("Spatial components")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

nc <- sf::st_read(system.file("shape/nc.shp", package="sf"), quiet=TRUE)
nc.nb <- spdep::poly2nb(nc)

test_that("creating spatial GMRF matrices works", {
  test <- Q_spatial(nc)
  expect_identical(dim(test), c(nrow(nc), nrow(nc)))
  expect_equal(crossprod(D_spatial(nc)), test)
  expect_equal(R_spatial(nc), matrix(1, nrow(nc), 1L))
  expect_identical(Q_spatial(nc.nb), test)
  expect_equal(crossprod(D_spatial(nc.nb)), test)
  expect_equal(R_spatial(nc.nb), matrix(1, nrow(nc), 1L))
  expect_lt(length(Q_spatial(nc, queen=FALSE)@x), length(test@x))
})

test_that("spatial model works", {
  expect_warning(
    sampler <- create_sampler(
      BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, poly.df=nc), name="vs"),
      data=nc
    ), "deprecated"
  )
  expect_warning(
    sampler <- create_sampler(
      BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc, derive.constraints=TRUE), name="vs"),
      data=nc
    ), "deprecated"
  )
  sampler <- create_sampler(
    BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc), constr=FALSE, name="vs"),
    data=nc
  )
  expect_null(sampler$mod[["vs"]][["R"]])
  sampler <- create_sampler(
    BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc), name="vs"),
    data=nc
  )
  sim <- MCMCsim(sampler, burnin=100, n.iter=200, n.chain=2,
                 store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_identical(nrow(summ$vs), nrow(nc))
  expect_lt(sum(summ$vs[, "Mean"]), sqrt(.Machine$double.eps))
  sampler <- create_sampler(
    BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID), name="vs"),
    data=nc
  )
  # if no spatial structure is provided, it is assumed that data contains it
  sim <- MCMCsim(sampler, burnin=100, n.iter=200, n.chain=2,
                 store.all=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_identical(nrow(summ$vs), nrow(nc))
  expect_lt(sum(summ$vs[, "Mean"]), sqrt(.Machine$double.eps))
})

test_that("arguments of spatial() are looked up in the right environment", {
  snap <- 0
  queen <- FALSE
  sampler <- create_sampler(
    BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc, snap=snap, queen=queen), name="vs"),
    data=nc
  )
  expect_identical(sampler$mod[["vs"]]$QA, Q_spatial(nc, snap=snap, queen=queen))
  f <- function() {
    sn <- 1
    qu <- TRUE
    sampler <- create_sampler(
      BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc, snap=sn, queen=qu), name="vs"),
      data=nc
    )
  }
  sampler <- f()
  expect_identical(sampler$mod[["vs"]]$QA, Q_spatial(nc, snap=1, queen=TRUE))
  expect_error(
    sampler <- create_sampler(
      BIR74 ~ SID74 + gen(factor = ~ spatial(CNTY_ID, graph=nc, snap=snap, king=queen), name="vs"),
      data=nc
    ), "unused argument"
  )
})
